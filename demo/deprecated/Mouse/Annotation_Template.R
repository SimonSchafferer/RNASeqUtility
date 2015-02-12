###################################
#   PATH definitions - theese need to be changed!
###################################
rootDir = file.path( "/home","simon","PHDStudies","RNA-Seq","IonProton","CavKO_Striessnig","MouseBrain_wo_Ref_extended")#schaffrr
infernalDB = file.path("/home","simon","dbsOfflineUse","InfernalDB","Rfam.cm.1_1")

# New Version -> one may load the workspace of the mapping assembly and clustering approach MappingAssemblyClustering.rda
load(file.path( rootDir,"MappingAssemblyClustering.rda"))#
#extracting a sample fastq subset! sampleFastq.fastq (100.000 entries)
options(stringsAsFactors = FALSE)
setwd(rootDir)

annotationDir = file.path(rootDir,"annotation")
dir.create(annotationDir)

require(sncRNAannotation)
# require(reshape2)
require(rtracklayer)
# require(ggplot2)
library(biomaRt)
#Range based annotation 
#Fetch sequences from genome
library("BSgenome")
# available.genomes()
library(BSgenome.Mmusculus.UCSC.mm10)
bsgenome = BSgenome.Mmusculus.UCSC.mm10

##########################################################################################################################
#                               Annotation of small ncRNAs that were mapped to ENSEMBL
##########################################################################################################################

###############################################################################################
#         BiomaRt annotation
###############################################################################################
ensReadCountL = lapply( bamToBedAndMergencRNAL, function(x){
  x = x[[2]]
  return( import( file.path( getOutFilePath(getCLIApplication(x)),getOutResultName(getOutResultReference(x))), format="BED",asRangedData=FALSE ) )
})

ensembl=useMart("ensembl")
# human = useDataset("hsapiens_gene_ensembl",mart=ensembl)
mouse = useDataset("mmusculus_gene_ensembl",mart=ensembl)
# listAttributes(mouse)[grep("MGI",listAttributes(mouse)$description),]
# listFilters(mouse)[grep("NM_",listAttributes(mouse)$description),]

#First extract all ENSEMBLTranscriptIDs from all experiments
ensTranscrNames = unique( unlist( lapply(ensReadCountL, function(x){
  return(as.character(seqnames(x)))
}) ) )

#Fetch the annotation to the ensTranscrNames
ensemblAnnotDF = getBM(attributes=c("ensembl_transcript_id","ensembl_gene_id","gene_biotype","external_gene_name","rfam",
                                    "chromosome_name","start_position","end_position","strand",
                                    "mgi_symbol"), filters="ensembl_transcript_id", mart=mouse, 
                       values=ensTranscrNames)
ensemblAnnotDF$chromosome_name = paste0("chr",ensemblAnnotDF$chromosome_name)
#Create a reference dataframe containing all entries (Ensembl ids plus strand as ID)
overallTable = unique(do.call( rbind, lapply( ensReadCountL, function(x){
  xdf = as.data.frame(x)
  colnames(xdf) = c("ensembl_transcript_id", "mapStart","mapEnd","mapWidth","mapStrand","readCount","mapQ")
  xdf$UID = paste0(xdf$ensembl_transcript_id, ifelse(xdf$mapStrand == "-", "_minus", "_plus") )
  return(xdf[,c("ensembl_transcript_id","UID")])
}) ) )

#Now merge with counted reads
ensReadsCountDF = do.call( cbind, lapply( ensReadCountL, function(x){
  xdf = as.data.frame(x)
  colnames(xdf) = c("ensembl_transcript_id", "mapStart","mapEnd","mapWidth","mapStrand","readCount","mapQ")
  xdf$UID = paste0(xdf$ensembl_transcript_id, ifelse(xdf$mapStrand == "-", "_minus", "_plus") )
  #Merge The duplicated UIDs
  xdf = do.call(rbind, lapply( split(xdf, xdf$UID), function(x){
    xnew = x[1,]
    xsize = dim(x)[1]
    if( xsize != 1 ){
      xnew$mapStart = min(as.numeric(x$mapStart))
      xnew$mapEnd = min(as.numeric(x$mapEnd))
      #       xnew$mapWidth = as.numeric(xnew$mapEnd)-as.numeric(xnew$mapStart)+1
      xnew$readCount = sum(as.numeric(x$readCount),na.rm = TRUE)
      xnew$mapQ = mean(as.numeric(x$mapQ),na.rm = TRUE)
    }
    return(xnew)
  }))
  xdf_merged = merge(overallTable, xdf[,c("UID", "readCount","mapStart","mapEnd")], by="UID", all.x=TRUE)
  return(xdf_merged)
} ) )

#Change for mapStart take minimum and map End take maximum? -> else rowMeans should be ok too, but round!
mapStart = ensReadsCountDF[,grep("mapStart", colnames(ensReadsCountDF)) ]
mapStart = apply(mapStart, 1, min, na.rm=TRUE)
mapEnd = ensReadsCountDF[,grep("mapEnd", colnames(ensReadsCountDF)) ]
mapEnd = apply(mapEnd, 1, min, na.rm=TRUE)


ensReadsCountDF = ensReadsCountDF[, c(1, grep("readCount", colnames(ensReadsCountDF))) ]
colnames(ensReadsCountDF)[1] = "UID"
colnames(ensReadsCountDF)[2:length(colnames(ensReadsCountDF))] = samplesInfo$sampleName
#Changing all NA values to zero
ensReadsCountDF[,2:length(colnames(ensReadsCountDF))] = apply(ensReadsCountDF[,2:length(colnames(ensReadsCountDF))], 2, function(x){
  ifelse( is.na(x), 0, x)
})
ensReadsCountDF$mapStart = mapStart
ensReadsCountDF$mapEnd = mapEnd

#can nicely be clustered by the rfam ID!
##########################
# Clustering
##########################

ensReadsCountDF$ensembl_transcript_id = sub("_.*$","",ensReadsCountDF$UID)
ensReadsCountDF_clustered = merge( ensReadsCountDF, ensemblAnnotDF, all.x=TRUE, by="ensembl_transcript_id" )
ensReadsCountDF_clustered$rfam = ifelse(ensReadsCountDF_clustered$rfam == "",ensReadsCountDF_clustered$UID, ensReadsCountDF_clustered$rfam)
ensReadsCountDF_clusteredL = split(ensReadsCountDF_clustered, ensReadsCountDF_clustered$rfam)
ensReadsCountDF_clustered = do.call(rbind, lapply(ensReadsCountDF_clusteredL, function(x){
  readCountColumns = 3:(2+length(samplesInfo$condition))
  xnew = x[1,]
  xsize = dim(x)[1]
  if( xsize != 1 ){
    xnew[,readCountColumns] = colMeans(apply( x[,readCountColumns],2, as.numeric), na.rm = TRUE)   
    xnew[,readCountColumns] = apply( xnew[,readCountColumns],2, round, digits=0)
  }
  xnew$clustered = xsize != 1
  return(xnew)
}) )
ensReadsCountDF_clustered$rfam[grep("ENS",ensReadsCountDF_clustered$rfam)] = ""
###################################
#     Clustered Read Count Table
###################################
head(ensReadsCountDF_clustered)
table(ensReadsCountDF_clustered$gene_biotype)
write.table(ensReadsCountDF_clustered, file.path(annotationDir,"ensReadsCountDF_clustered.csv"),sep="\t",row.names=FALSE, col.names=TRUE )
write.table(ensReadsCountDF, file.path(annotationDir,"ensReadsCountDF.csv"),sep="\t",row.names=FALSE, col.names=TRUE )

##########################################################################################################################
#                               Annotation of small ncRNA contigs assembled
##########################################################################################################################
contigForCountingGR_clustered = contigForCountingGR_clustered[ which(width(contigForCountingGR_clustered) >= 18) ]

#ENSEMBL Range Based annotation
ensembl_annot = EnsemblAnnotation("ensemblAnnot", system.file("resources/ensembl/", package="sncRNAannotation"),"ensembl_gtf_v75_mm10.rda", 
                                  contigForCountingGR_clustered)
ensembl_annotL = annotationSummary(ensembl_annot)
ensembl_protCodAnnot = ensembl_annotL$protCodingDF
ensembl_featureAnnot = ensembl_annotL$featureDF

tmpRelevant = which(ensembl_protCodAnnot$gene_type %in% c("exon","intron","exon/intron"))
ensembl_protCodAnnot$gene_type_ext = ensembl_protCodAnnot$gene_type
ensembl_protCodAnnot$gene_type_ext[tmpRelevant] = mapply( function(x,y,z){
  if( x == y ){
    z 
  } else{
    paste0(z ,"_as")    
  }
}, ensembl_protCodAnnot$gene_strand[tmpRelevant], 
as.character(strand(contigForCountingGR_clustered[tmpRelevant])), 
ensembl_protCodAnnot$gene_type[tmpRelevant]  
)

featureAnnot = ifelse( is.na(ensembl_featureAnnot$feature_biotype), "unknown", ensembl_featureAnnot$feature_biotype )
unknownEntries = which( is.na(ensembl_featureAnnot$feature_biotype))
featureAnnot[unknownEntries] = paste0(ensembl_protCodAnnot$gene_type_ext[unknownEntries] ,"_unknown")
featureAnnot_short = ifelse(is.na(ensembl_featureAnnot$feature_biotype),"unknown",ensembl_featureAnnot$feature_biotype)
ensembl_featureAnnot$featureAnnot_short = featureAnnot_short
table(ensembl_featureAnnot$featureAnnot_short)

#Repeat Masker Annotation
rpmskr_annot = GRangesBasedAnnotation("repeatMasker", "/home/simon/dbsOfflineUse/short_ncRNA_annotation_offlineResources","repeat_mskr_ucsc_mm10.rda", 
                                      contigForCountingGR_clustered)
rpmskr_annotDF = annotationSummary(rpmskr_annot)

##########################
#   Combine Tables
##########################
colnames( ensembl_protCodAnnot ) = paste0("g_",colnames( ensembl_protCodAnnot ))
colnames( ensembl_featureAnnot ) = paste0("f_",colnames( ensembl_featureAnnot ))
colnames( rpmskr_annotDF ) = paste0("r_",colnames( rpmskr_annotDF ))

contigAnnotTable = cbind(ensembl_protCodAnnot, ensembl_featureAnnot, rpmskr_annotDF)
contigAnnotTable$contigID = contigForCountingGR_clustered$name
##################################################
#     Annotate with infernal
##################################################

fastqForAnnotation = "contigSequences.fa"
contigForCountingGR_clustered_seq = getSeq(bsgenome, contigForCountingGR_clustered)
names(contigForCountingGR_clustered_seq) = contigForCountingGR_clustered$name
#write fasta file: 
writeLines( as.vector( t(cbind( paste0(">",names(contigForCountingGR_clustered_seq)), as.character(contigForCountingGR_clustered_seq)  )) ), con = file.path(annotationDir,fastqForAnnotation) )


# #PUT RFAM INTO R PACAKGE -> ONLY 
# ##### infernal runs
# #cmscan --tblout tmp -E 20 --noali Rfam.cm.1_1 fastaFiles/testSeq.fa #NICELY PARSEABLE FILE -> PERL SCRIPT!
########### INFERNAL RUN ####################
infernalOutFN = "contigAnnotation_infernal"
setwd(annotationDir)
cmd1 = paste0( "cmscan --tblout ",file.path(annotationDir, paste0(infernalOutFN,"_orig") ), " --noali --cpu 8 -E 20 ",infernalDB," ", file.path(annotationDir,fastqForAnnotation) )
cmd2 = paste0( "awk -v OFS=\"\\t\" \'$1=$1\' ",file.path(annotationDir, paste0(infernalOutFN,"_orig") )," > ",file.path(annotationDir,infernalOutFN ) )

#Running infernal
system(cmd1)
system(cmd2)

tmp = read.csv(file.path(annotationDir,infernalOutFN), sep="\t", header=FALSE, comment.char="#")
cnames = c("targetName","accession","queryName","accession","mdl","mdlFrom","mdlTo","seqFrom","seqTo","strand","trunc","pass","gc","bias","score","Evalue","inc","descriptionOfTarget")
colnames(tmp) = cnames
bestMatch = tmp[which(!duplicated(tmp[,3])),]
# bestMatch = bestMatch[bestMatch[,17]!="?",] may be too restrictive
descriptions = gsub( " {2,}.*$","", apply( bestMatch[,which(colnames(bestMatch) == "descriptionOfTarget"):dim(bestMatch)[2]], 1, function(x){paste0(x,collapse=" ")} ))
bestMatch = bestMatch[,1:which(colnames(bestMatch) == "descriptionOfTarget")]
bestMatch[,dim(bestMatch)[2]] = descriptions
head(bestMatch)

colnames(bestMatch) = paste0("inf_",colnames(bestMatch))
contigAnnotTable_fin = merge( contigAnnotTable, bestMatch, by.x="contigID", by.y="inf_queryName", all.x=TRUE)

# contigAnnotTable_fin[,c("f_feature_biotype","r_RepeatName","inf_targetName")]

write.table(contigAnnotTable_fin, file.path(annotationDir,"contigAnnotTable_fin.csv"),sep="\t",row.names=FALSE, col.names=TRUE )

ensReadsCountDF_clusteredGR = with( ensReadsCountDF_clustered, GRanges( seqnames=chromosome_name, 
                                                                        IRanges((start_position+mapStart-1),
                                                                                (start_position+mapEnd-1)), strand=strand,
                                                                        ensembl_transcript_id, UID, ensembl_gene_id, 
                                                                        gene_biotype,external_gene_name,rfam,mgi_symbol,clustered ))

# findOverlaps(ensReadsCountDF_clusteredGR, contigForCountingGR_clustered)
# ensReadsCountDF_clusteredGR[591]
# contigForCountingGR_clustered[27]
# contigAnnotTable_fin[contigAnnotTable_fin$contigID == "contig27",]

# ###########################
# load Coverage Matrix
# ###########################
count_contigsDF = read.table( file=file.path( getOutFilePath(getCLIApplication(multiBamCov_CLI_cmdRes)), getOutResultName(getOutResultReference(multiBamCov_CLI_cmdRes))), 
                              sep="\t", header=FALSE )
colnames(count_contigsDF) = c("chr","start","end","contigID","uniqueness","strand",samplesInfo$sampleName)
count_contigsGR = with(count_contigsDF, GRanges(chr, IRanges(start, end), uniqueness=uniqueness, strand=strand,contigID=contigID,
                                                count_contigsDF[,7:dim(count_contigsDF)[2]]))

save.image(file.path(rootDir, "Annotation.rda"))



# 
# 
# chr13:86184594-86184653
# chr5:23216809-23216879
# chr5:23358028-23358098
# 
# chr13 86184594  86184653 . 3
# chr5  23216809 23216879 .  4
# chr5  23358028  23358098  .  4
# #liftover
# chr13  86044989  86045048  .	1 #contig 412!
# chr5	23710991	23711061	.	1 # NOT FOUND
# chr5	23852210	23852280	.	1 # NOT FOUND
# 
# snoRNAs = GRanges(c("chr13","chr5","chr5"), IRanges(c(86044989,23852210,23710991),c(86045048,23852280,23711061)), strand="*", names=c("e307","e470_1","e470_2"))
# 
# snoRNAOverL = lapply( bamToBedAndMergeL, function(x){
#   curr = x[[2]]
#   outP = getOutFilePath(getCLIApplication(curr))
#   outFN = getOutResultName(getOutResultReference(curr))
#   currCont = import( file.path(outP, outFN), asRangedData=FALSE  )
#   return( subsetByOverlaps(currCont,snoRNAs) )
# })
# 
# 
# contigsA = import("/home/simon/PHDStudies/RNA-Seq/IonProton/CavKO_Striessnig/MouseBrain_wo_Ref/bowtieMappings/IonXpressRNA_001.R_2014_08_22_11_48_07_user_PRO-36-2014-08-22_MoBrain_Ion_RNA_-_20-200nt_Transcriptom_trimmed_bowtieOut/IonXpressRNA_001_mapped_s_contigs.bed",asRangedData=FALSE)
# 
# 
# 
# 
# #####################
# # Plotting all Reads with short annotation to get an overview
# #####################
# df = countDF[,READCOLCOVMAT:dim(countDF)[2]]
# df$annotation = featureAnnot_short
# df = melt(df)
# df$group = sub("_.*","",df$variable)
# 
# g4 = ggplot(df, aes(x=annotation, y=log10(value), fill=annotation )) + geom_boxplot() + facet_grid(variable ~.) + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="ncRNA class", fill ="ncRNA class")
# 
# g5 = ggplot(df, aes(x=annotation, y=log10(value), fill=annotation )) + geom_boxplot() + facet_grid(group ~.) + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylab("log10(Read Coverage)") + labs(x="ncRNA class", fill ="ncRNA class")
# 
# 
# ###############################################
# #   Adding repeats to the annotation
# ###############################################
# 
# library(sncRNAannotation)
# rpmskr_annot = GRangesBasedAnnotation("repeatMasker", "/home/simon/dbsOfflineUse/short_ncRNA_annotation_offlineResources","repeat_mskr_ucsc_mm9.rda", 
#                                       readCoverageMatrixGR)
# rpmskr_annotDF = annotationSummary(rpmskr_annot)
# table(is.na(rpmskr_annotDF$RepeatName))
# table(rpmskr_annotDF$repeatClass)
# readCoverageMatrixGR_noRepeats = readCoverageMatrixGR[which(is.na(rpmskr_annotDF$RepeatName))]
# ensembl_annot = EnsemblAnnotation("ensemblAnnot", system.file("resources/ensembl/", package="sncRNAannotation"),"ensembl_gtf_v67_mm9.rda", 
#                                   readCoverageMatrixGR_noRepeats)
# ensembl_annotL = annotationSummary(ensembl_annot)
# ensembl_protCodAnnot = ensembl_annotL$protCodingDF
# ensembl_featureAnnot = ensembl_annotL$featureDF
# 
# #Adding sense or antisense to the annotation
# tmpRelevant = which(ensembl_protCodAnnot$gene_type %in% c("exon","intron","exon/intron"))
# ensembl_protCodAnnot$gene_type_ext = ensembl_protCodAnnot$gene_type
# ensembl_protCodAnnot$gene_type_ext[tmpRelevant] = mapply( function(x,y,z){
#   if( x == y ){
#     z 
#   } else{
#     paste0(z ,"_as")    
#   }
# }, ensembl_protCodAnnot$gene_strand[tmpRelevant], 
# as.character(strand(readCoverageMatrixGR_noRepeats[tmpRelevant])), 
# ensembl_protCodAnnot$gene_type[tmpRelevant]  
# )
# 
# table(ensembl_protCodAnnot$gene_type_ext)
# table(ensembl_featureAnnot$feature_biotype)
# table( paste0( ensembl_protCodAnnot$gene_type_ext,"_", ifelse(is.na(ensembl_featureAnnot$feature_biotype),"unknown",ensembl_featureAnnot$feature_biotype) ) )
# 
# 
# 
# 
# grangesAnnotationDFL = function( grL ){
#   require(sncRNAannotation)
#   require(reshape2)
#   require(rtracklayer)
#   
#   annotL = lapply( names(grL), function(grL_name){
#     
#     ensembl_annot = EnsemblAnnotation("ensemblAnnot", system.file("resources/ensembl/", package="sncRNAannotation"),"ensembl_gtf_v67_mm9.rda", 
#                                       grL[[grL_name]])
#     ensembl_annotL = annotationSummary(ensembl_annot)
#     ensembl_protCodAnnot = ensembl_annotL$protCodingDF
#     ensembl_featureAnnot = ensembl_annotL$featureDF
#     
#     tmpRelevant = which(ensembl_protCodAnnot$gene_type %in% c("exon","intron","exon/intron"))
#     ensembl_protCodAnnot$gene_type_ext = ensembl_protCodAnnot$gene_type
#     ensembl_protCodAnnot$gene_type_ext[tmpRelevant] = mapply( function(x,y,z){
#       if( x == y ){
#         z 
#       } else{
#         paste0(z ,"_as")    
#       }
#     }, ensembl_protCodAnnot$gene_strand[tmpRelevant], 
#     as.character(strand(grL[[grL_name]][tmpRelevant])), 
#     ensembl_protCodAnnot$gene_type[tmpRelevant]  
#     )
#     
#     featureAnnot = ifelse( is.na(ensembl_featureAnnot$feature_biotype), "unknown", ensembl_featureAnnot$feature_biotype )
#     unknownEntries = which( is.na(ensembl_featureAnnot$feature_biotype))
#     featureAnnot[unknownEntries] = paste0(ensembl_protCodAnnot$gene_type_ext[unknownEntries] ,"_unknown")
#     featureAnnot = featureAnnot[which(! featureAnnot %in% c("sense_intronic","sense_overlapping"))]
#     
#     featureAnnot_short = ifelse(is.na(ensembl_featureAnnot$feature_biotype),"unknown",ensembl_featureAnnot$feature_biotype)
#     
#     df = melt(table(featureAnnot))
#     colnames(df)[1] = "variable"
#     df$group = grL_name
#     df = df[df$value > 2,]
#     
#     df_short = melt(table(featureAnnot_short))
#     colnames(df_short)[1] = "variable"
#     df_short$group = grL_name
#     df_short = df_short[df_short$value > 2 ,]
#     return( list("AnnotDF"=df, "AnnotDF_short"=df_short))
#     
#   } )
#   
#   dfAnnot = do.call( rbind, lapply(annotL, function(x){
#     return( x[[1]] )
#   }) )
#   
#   dfAnnot_short = do.call( rbind, lapply(annotL, function(x){
#     return( x[[2]] )
#   }) )
#   
#   return( list("AnnotDF"=dfAnnot, "AnnotDF_short"=dfAnnot_short))
# }
