print("Annotation")
###################################
#   PATH definitions - theese need to be changed!
###################################

args <- commandArgs(TRUE)
if( length(args) == 0 ){
  load( file.path(getwd(),"Configuration.rda") )
} else{
  load(file.path(args[1],"Configuration.rda"))
}

options(stringsAsFactors = FALSE)
setwd(rootDir)
load(file.path( rootDir,"MappingAssemblyClustering.rda"))#


library(sncRNAannotation)
# require(reshape2)
require(rtracklayer)
# require(ggplot2)
require(biomaRt)
#Range based annotation 
#Fetch sequences from genome
library(CLIHelperPackage)
library(RNASeqUtility)
# available.genomes()
# source("http://bioconductor.org/biocLite.R")
# biocLite("BSgenome.Hsapiens.NCBI.GRCh38")
# library(BSgenome.Mmusculus.UCSC.mm10)
# bsgenome = BSgenome.Mmusculus.UCSC.mm10
# library(BSgenome.Hsapiens.NCBI.GRCh38)
# library(BSgenome.Hsapiens.NCBI.GRCh38)
# library(TxDb.Hsapiens.UCSC.hg38.knownGene)
# seqlevelsStyle(bs) = "UCSC" ## convert to UCSC style
# seqlevels(BSgenome.Hsapiens.NCBI.GRCh38)
# seqlevels(bs)
# seqlevels(TxDb.Hsapiens.UCSC.hg38.knownGene) 

#check out more with: listAttributes(useMart("ensembl"))[grep("MGI",listAttributes(useMart("ensembl"))$description),]
switch(organismForAnnotation, 
       hg19={
         print('hg19')
         ensembl=useMart("ENSEMBL_MART_ENSEMBL",host = "www.ensembl.org")
         ensemblMart = useDataset("hsapiens_gene_ensembl",mart=ensembl)
         ensemblAttributes = c("ensembl_transcript_id","ensembl_gene_id","gene_biotype","external_gene_name","rfam",
                               "chromosome_name","start_position","end_position","strand",
                               "hgnc_symbol")
         ensemblRangeBasedAnnotationDir = system.file("resources/ensembl/", package="sncRNAannotation")
         ensemblRangeBasedAnnotationFN = "ensembl_gtf_hg19.rda"
         genomeToAnnotate = "BSgenome.Hsapiens.UCSC.hg19"
       },
       hg38={
         print('hg38')
         ensembl=useMart("ENSEMBL_MART_ENSEMBL", host = "www.ensembl.org")
         ensemblMart = useDataset("hsapiens_gene_ensembl",mart=ensembl)
         ensemblAttributes = c("ensembl_transcript_id","ensembl_gene_id","gene_biotype","external_gene_name","rfam",
                               "chromosome_name","start_position","end_position","strand",
                               "hgnc_symbol")
         ensemblRangeBasedAnnotationDir = system.file("resources/ensembl/", package="sncRNAannotation")
         ensemblRangeBasedAnnotationFN = "ensembl_gtf_hg38_78.rda"
         genomeToAnnotate = "BSgenome.Hsapiens.UCSC.hg38"
       },
       mm10={
         print('mm10')
         ensembl=useMart("ENSEMBL_MART_ENSEMBL", host = "www.ensembl.org")
         ensemblMart = useDataset("mmusculus_gene_ensembl",mart=ensembl)
         ensemblAttributes=c("ensembl_transcript_id","ensembl_gene_id","gene_biotype","external_gene_name","rfam",
                             "chromosome_name","start_position","end_position","strand",
                             "mgi_symbol")
         ensemblRangeBasedAnnotationDir = system.file("resources/ensembl/", package="sncRNAannotation")
         ensemblRangeBasedAnnotationFN = "ensembl_gtf_v75_mm10.rda"
         genomeToAnnotate = "BSgenome.Mmusculus.UCSC.mm10"
         print('bar') 
       }, 
       mm9={
         print('mm9')
         ensembl=useMart("ENSEMBL_MART_ENSEMBL", host = "www.ensembl.org")
         ensemblMart = useDataset("mmusculus_gene_ensembl",mart=ensembl)
         ensemblAttributes=c("ensembl_transcript_id","ensembl_gene_id","gene_biotype","external_gene_name","rfam",
                             "chromosome_name","start_position","end_position","strand",
                             "mgi_symbol")
         ensemblRangeBasedAnnotationDir = system.file("resources/ensembl/", package="sncRNAannotation")
         ensemblRangeBasedAnnotationFN = "ensembl_gtf_v67_mm9.rda"
         genomeToAnnotate = "BSgenome.Mmusculus.UCSC.mm9"
         print('bar') 
       },
       {
         stop('Annotation Type not supported!')
       }
)

library( genomeToAnnotate, character.only=TRUE)
assign("bsgenome", eval(parse( text= genomeToAnnotate)) )
seqlevelsStyle(bsgenome) = "UCSC" #This sets the chrodomosome coordinates from ENSEMBL to UCSC (1 to chr1) -> currently only relevant in case of BSgenome.Hsapiens.UCSC.hg38
##########################################################################################################################
#                               Annotation of small ncRNAs that were mapped to ENSEMBL
##########################################################################################################################

###############################################################################################
#         BiomaRt annotation
###############################################################################################
ncRNAReadCountL = lapply( bamToBedAndMergencRNAL, function(x){
  x = x[[2]]
  return( import( file.path( getOutFilePath(getCLIApplication(x)),getOutResultName(getOutResultReference(x))),asRangedData=FALSE ) )
})

#Subset These to Ensembl only
ensReadCountL = lapply(ncRNAReadCountL, function(x){
  x[grep("ENS.*",seqnames(x))]
})

#First extract all ENSEMBLTranscriptIDs from all experiments
ensTranscrNames = unique( unlist( lapply(ensReadCountL, function(x){
  return(as.character(seqnames(x)))
}) ) )

#Fetch the annotation to the ensTranscrNames
ensemblAnnotDF = getBM(attributes=ensemblAttributes, filters="ensembl_transcript_id", mart=ensemblMart, 
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

##########################
# Clustering by RFAM ID (since RFAM is annotating ncRNA classes and not individual transcripts!)
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

####################################
#     Other ncRNAs (miRNAs, snoRNAs)
####################################
otherncRNAReadCountL = lapply(ncRNAReadCountL, function(x){
  x[-grep("ENS.*",seqnames(x))]
})

#Check if non-ensembl genes are present
otherReadsCountDF = data.frame()
if( length(otherncRNAReadCountL[[1]]) != 0 ){
  
  #Same as before
  overallTable_other = unique(do.call( rbind, lapply( otherncRNAReadCountL, function(x){
    xdf = as.data.frame(x)
    colnames(xdf) = c("transcript_id", "mapStart","mapEnd","mapWidth","mapStrand","readCount","mapQ")
    xdf$UID = paste0(xdf$transcript_id, ifelse(xdf$mapStrand == "-", "_minus", "_plus") )
    return(xdf[,c("transcript_id","UID")])
  }) ) )
  
  #Now merge with counted reads
  otherReadsCountDF = do.call( cbind, lapply( otherncRNAReadCountL, function(x){
    xdf = as.data.frame(x)
    colnames(xdf) = c("transcript_id", "mapStart","mapEnd","mapWidth","mapStrand","readCount","mapQ")
    xdf$UID = paste0(xdf$transcript_id, ifelse(xdf$mapStrand == "-", "_minus", "_plus") )
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
    xdf_merged = merge(overallTable_other, xdf[,c("UID", "readCount","mapStart","mapEnd")], by="UID", all.x=TRUE)
    return(xdf_merged)
  } ) )
  
  #Change for mapStart take minimum and map End take maximum? -> else rowMeans should be ok too, but round!
  mapStart = otherReadsCountDF[,grep("mapStart", colnames(otherReadsCountDF)) ]
  mapStart = apply(mapStart, 1, min, na.rm=TRUE)
  mapEnd = otherReadsCountDF[,grep("mapEnd", colnames(otherReadsCountDF)) ]
  mapEnd = apply(mapEnd, 1, min, na.rm=TRUE)
  
  
  otherReadsCountDF = otherReadsCountDF[, c(1, grep("readCount", colnames(otherReadsCountDF))) ]
  colnames(otherReadsCountDF)[1] = "UID"
  colnames(otherReadsCountDF)[2:length(colnames(otherReadsCountDF))] = samplesInfo$sampleName
  #Changing all NA values to zero
  otherReadsCountDF[,2:length(colnames(otherReadsCountDF))] = apply(otherReadsCountDF[,2:length(colnames(otherReadsCountDF))], 2, function(x){
    ifelse( is.na(x), 0, x)
  })
  otherReadsCountDF$mapStart = mapStart
  otherReadsCountDF$mapEnd = mapEnd
  write.table(otherReadsCountDF, file.path(annotationDir,"otherReadsCountDF.csv"),sep="\t",row.names=FALSE, col.names=TRUE )
}

#Now the question is, if one should eliminate miRNAs and snoRNAs, or just leave them 
write.table(ensReadsCountDF_clustered, file.path(annotationDir,"ensReadsCountDF_clustered.csv"),sep="\t",row.names=FALSE, col.names=TRUE )
write.table(ensReadsCountDF, file.path(annotationDir,"ensReadsCountDF.csv"),sep="\t",row.names=FALSE, col.names=TRUE )

##########################################################################################################################
#                               Annotation of small ncRNA contigs assembled
##########################################################################################################################
contigForCountingGR_clustered = contigForCountingGR_clustered[ which(width(contigForCountingGR_clustered) >= 18) ]

#ENSEMBL Range Based annotation
ensembl_annot = EnsemblAnnotation("ensemblAnnot", ensemblRangeBasedAnnotationDir,ensemblRangeBasedAnnotationFN, 
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

rpmskr_annot = GRangesBasedAnnotation("repeatMasker", repeatMaskerDir, repeatMaskerFN,
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

########### INFERNAL RUN ####################
infernalOutFN = "contigAnnotation_infernal"
setwd(annotationDir)
cmd1 = paste0( "cmscan --tblout ",file.path(annotationDir, paste0(infernalOutFN,"_orig") ), " --noali --cpu 8 -E 20 ",infernalDB," ", file.path(annotationDir,fastqForAnnotation) )
#AWK script to make the file parseable (different number of white spaces -> no csv possible from cmscan output)
cmd2 = paste0( "awk -v OFS=\"\\t\" \'$1=$1\' ",file.path(annotationDir, paste0(infernalOutFN,"_orig") )," > ",file.path(annotationDir,infernalOutFN ) )

#Running infernal
system(cmd1)
system(cmd2)

#THESE colnames are hard coded from the infernal output, should remain the same -> check if infernal version changes!
tmp = read.csv(file.path(annotationDir,infernalOutFN), sep="\t", header=FALSE, comment.char="#")
cnames = c("targetName","accession","queryName","accession","mdl","mdlFrom","mdlTo","seqFrom","seqTo","strand","trunc","pass","gc","bias","score","Evalue","inc","descriptionOfTarget")
colnames(tmp) = cnames
bestMatch = tmp[which(!duplicated(tmp[,3])),]
# bestMatch = bestMatch[bestMatch[,17]!="?",] may be too restrictive -> this can be done afterwards. The ! indicates a sound hit, wheareas the ? marks not sure
descriptions = gsub( " {2,}.*$","", apply( bestMatch[,which(colnames(bestMatch) == "descriptionOfTarget"):dim(bestMatch)[2]], 1, function(x){paste0(x,collapse=" ")} ))
bestMatch = bestMatch[,1:which(colnames(bestMatch) == "descriptionOfTarget")]
bestMatch[,dim(bestMatch)[2]] = descriptions
head(bestMatch)

colnames(bestMatch) = paste0("inf_",colnames(bestMatch))
#Merge the infernal annotations with the range based annotations
contigAnnotTable_fin = merge( contigAnnotTable, bestMatch, by.x="contigID", by.y="inf_queryName", all.x=TRUE)

#Writing the contig annotation table
write.table(contigAnnotTable_fin, file.path(annotationDir,"contigAnnotTable_fin.csv"),sep="\t",row.names=FALSE, col.names=TRUE )

############################
# load Coverage Matrix
############################
count_contigsDF = read.table( file=file.path( getOutFilePath(getCLIApplication(multiBamCov_CLI_cmdRes)), getOutResultName(getOutResultReference(multiBamCov_CLI_cmdRes))), 
                              sep="\t", header=FALSE )
colnames(count_contigsDF) = c("chr","start","end","contigID","uniqueness","strand",samplesInfo$sampleName)
count_contigsGR = with(count_contigsDF, GRanges(chr, IRanges(start, end), uniqueness=uniqueness, strand=strand,contigID=contigID,
                                                count_contigsDF[,7:dim(count_contigsDF)[2]]))

save.image(file.path(rootDir, "Annotation.rda"))


#Testing
# ensReadsCountDF_clusteredGR = with( ensReadsCountDF_clustered, GRanges( seqnames=chromosome_name, 
#                                                                         IRanges((start_position+mapStart-1),
#                                                                                 (start_position+mapEnd-1)), strand=strand,
#                                                                         ensembl_transcript_id, UID, ensembl_gene_id, 
#                                                                         gene_biotype,external_gene_name,rfam,hgnc_symbol,clustered ))
# findOverlaps(ensReadsCountDF_clusteredGR, contigForCountingGR_clustered)
# ensReadsCountDF_clusteredGR[591]
# contigForCountingGR_clustered[27]
# contigAnnotTable_fin[contigAnnotTable_fin$contigID == "contig27",]
# count_contigsGR[count_contigsGR$contigID %in% contigAnnotTable_fin[which(contigAnnotTable_fin$inf_targetName == "tmRNA"),]$contigID]

