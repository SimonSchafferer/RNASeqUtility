#############################################
#     Contig Assembly - among samples - MultiIntersectBed
#   Now a PERL script I have written is doing the multi intersection for plus and minus strand bed files -> and merging them. 
#
#   First usecase: 
#     Basically, if no groups are defined then it uses multiInterSect Bed build in cluster algorithm that searches for an overlap in all samples and only reports the 
#     part of the overlap that is shared with most of the samples -> see bedtools doc... However the perl script adds one stringency that it must be present in all samples 
#     Therefore this script outputs the region that is shared by all samples as contig sequence. 
#   Second usecase: 
#     Groups are defined by e.g. 3 samples -> groups: 1 1 2 (first two in group 1 last in group 2): 
#     When groups are defined multiIntersectBed is called without cluster option -> which reports all intervals that overlap different groups. This is restricted by the fact 
#     that an overlap must be present among all group members of a particular group (e.g. WT WT WT TG TG -> when contig is shared among WT but not among TG it gets reported).
#     This file is them merged to obatin one bed file that contains the longest possible contig over all groups when it is present in all members of a subgroup. 
#     This can be used the same way as the opposite of clustering algorithm explained above by specifying one group only -> the output will then be a region that is the maximum,
#     of all combined: 
#     Grouping:                           No grouping: 
#           A   -------------------         A   -------------------
#           B       ----------              B       ----------
#           C ---------------               C ---------------
#   report:   ---------------------                 ---------
#   Both approaches have their advantages and disadvantages. An advantage of non-grouping would be that when counting only reads that overlap 80-90%, the contigs that are 
#   compared by differential expression will not suffer by a length bias....
#############################################

###################################
#   PATH definitions - theese need to be changed!
###################################

#extracting a sample fastq subset! sampleFastq.fastq (100.000 entries)
options(stringsAsFactors = FALSE)
rootDir = file.path( "/media","KINGSTON","RNA-Seq","MiniAnalysis_multiIntersect")#schaffrr
rawDataDir = file.path(rootDir,"rawData")
gtfFileName = "Mus_musculus.GRCm38.78_mod.gtf"
gtfFilePath = "/home/simon/dbsOfflineUse/GTF_repos/Mus_musculus.GRCm38.78_repo/"
transcriptomeIndex = "/home/simon/dbsOfflineUse/GTF_repos/Mus_musculus.GRCm38.78_repo/Mus_musculus.GRCm38.78"
bowtieIndexFilePath = "/home/simon/dbsOfflineUse/MusMusculus/Mouse_mm10_fasta/bowtie/mm10"
bowtieMappingsDir = file.path(rootDir,"bowtieMappings")
dir.create(bowtieMappingsDir)
#Contig Assembly directory
contigAssemblyDir = file.path(rootDir,"contigAssembly")
dir.create(contigAssemblyDir)
#Contig Clustering Dir
contigClusterDir = file.path(rootDir,"contigClustering")
dir.create(contigClusterDir)
#Read Counts
readCountsDir = file.path(rootDir,"readCounts")
dir.create(readCountsDir)
# samplesInfo = read.table( file.path(rootDir, "samplesInfo.csv") ,sep="\t", header=TRUE)
commandLog = file.path(rootDir, "Commands") 
file.create(commandLog)
executionLog = file.path(rootDir, "executionLog") 
file.create(executionLog)
tmpCommandLog = ""

samplesInfo = read.table(file=file.path(rootDir, "samplesInfo.csv"), sep="\t", header=TRUE) #should contain cloumn condition and column sample name

if(!( "sampleName" %in% colnames(samplesInfo) & "condition"%in% colnames(samplesInfo))) {stop("Please provide a sampleInfo file with column condition and column sampleName")}

############
#   Parameter definition
############
read_threshold = 1 # read threshold for contig assembly
readOverlap_contig = -3 #This parameters defines the number of nt a read has to overlap with the next in order to be clustered (from positive to negative numbers -> APART default: 0)


####################################
#     Adding the PATHS to R environment, if executed within R!
####################################
pathvars = readLines(file.path(path.expand("~"),".bashrc"))
pathvars = pathvars[grep("export PATH\\=\\$PATH:",  pathvars )]
pathvars = sub( "export PATH\\=\\$PATH:", "",pathvars)
Sys.setenv(PATH=paste(Sys.getenv("PATH"),paste(pathvars,collapse=":"),sep=":")) 

#####################################
# Loading the required libraries
#####################################
setwd(rootDir)
library(CLIHelperPackage)
library(RNASeqUtility)


####################################################################################################################################################
#                                                             Start Analysis
####################################################################################################################################################

#######################
#   Cutadapt to trim the fastq files
#######################
cutAdaptCLI = Cutadapt_CLI(inFilePath=rawDataDir, cliParams =  c("--minimum-length 18"), outputFlag = "_trimmed", 
                           outFilePath = file.path(rootDir,"rawDataTrimmed") )
cutAdatptCLI_cmdRes = generateCommandResult(object = cutAdaptCLI )
#logging
tmpCommandLog = getCommandLog(cutAdatptCLI_cmdRes)


#################################################
#     Bowtie1 mapping, since it handles repeating ncRNA sequences more accordingly!!
#################################################
#Also create the bowtie base directory (it is already created within R, but for consistency of command line execution it is also piped to the command list here)

#logging
tmpCommandLog = c(tmpCommandLog, paste0("\nmkdir ", bowtieMappingsDir,"\n"))

fastQFiles = getOutResultName(getOutResultReference(cutAdatptCLI_cmdRes))
bowtieCLI_cmdResL = lapply( fastQFiles, function(fqf){
  bowtieCLI = Bowtie_CLI(inFilePath=getOutFilePath(getCLIApplication(cutAdatptCLI_cmdRes)), inFileNames=fqf, 
                         cliParams="-qSya -n 1 -p 3 --chunkmbs 1024 --best -l 28 -m 100 --strata", #--phred64-quals 
                         outputFlag="_mapped", outFilePath=file.path(bowtieMappingsDir,paste0(sub(".fastq","",fqf),"_bowtieOut")),
                         bowtieIndexFilePath=bowtieIndexFilePath, 
                         addNHTag=TRUE, outputPipeString="| samtools view -uhS -F4 - | samtools sort -", outputFormat="bam")
  return( generateCommandResult(bowtieCLI) )
} )
names(bowtieCLI_cmdResL) = fastQFiles

#logging
for( i in 1:length(bowtieCLI_cmdResL)  ){
  tmpCommandLog = c(tmpCommandLog, getCommandLog(bowtieCLI_cmdResL[[i]]) )
}


######################
#   Samtools commands sorting by Name and coordinates _sn and _s
######################

samToolsHTSeqCmdL = lapply( bowtieCLI_cmdResL, function(x){
  
  inFN = getOutResultName(getOutResultReference(x))
  currPath = getOutFilePath(getCLIApplication(x))
  ######################
  # Sam Commands For htseq count
  ######################
  sortSam1 = Samtools_CLI(inFilePath=currPath, 
                          inFileNames = inFN,
                          cliParams =  c("-n"), 
                          outputFlag = "_sn", 
                          outFilePath = currPath, 
                          samtoolsApplication = "sort", 
                          outputFormat = "bam")
  sortBamCmdRes1 = generateCommandResult( object = sortSam1 )
  
  inFN = getOutResultName(getOutResultReference(sortBamCmdRes1))
  
  #DO NOT GENERATE SAM FILE
  #   samView = Samtools_CLI(inFilePath=currPath, 
  #                          inFileNames = inFN,
  #                          cliParams =  c("-h"), 
  #                          outputFlag = "", 
  #                          outFilePath = currPath, 
  #                          samtoolsApplication = "view" ,
  #                          outputFormat = "sam")
  #   samViewCmdRes = generateCommandResult( object = samView )
  
  ######################
  # Sam Commands For multicov bedtools
  ######################
  
  inFN = getOutResultName(getOutResultReference(x))
  
  sortSam2 = Samtools_CLI(inFilePath=currPath, 
                          inFileNames = inFN,
                          cliParams =  c(""), 
                          outputFlag = "_s", 
                          outFilePath = currPath, 
                          samtoolsApplication = "sort", 
                          outputFormat = "bam")
  sortBamCmdRes2 = generateCommandResult( object = sortSam2 )
  
  
  inFN = getOutResultName(getOutResultReference(sortBamCmdRes2))
  
  samIndex = Samtools_CLI(inFilePath=currPath, 
                          inFileNames = inFN,
                          cliParams =  c(""), 
                          outputFlag = "", 
                          outFilePath = currPath, 
                          samtoolsApplication = "index", 
                          outputFormat = "bai")
  samIndexCmdRes = generateCommandResult( object = samIndex )
  
  return( list(sortBamCmdRes1, sortBamCmdRes2, samIndexCmdRes) )#samViewCmdRes
  
} )

#logging
for( i in 1:length(samToolsHTSeqCmdL)  ){
  for( j in 1:length(samToolsHTSeqCmdL[[i]])){
    tmpCommandLog = c(tmpCommandLog, getCommandLog(samToolsHTSeqCmdL[[i]][[j]]) )
  }
}

#############################################
#     Contig Assembly
#############################################
#ITERATING over the Bowtie Results and create a Bed file for each bam file
bamToBedAndMergeL = lapply( samToolsHTSeqCmdL, function(x){
  #converting bam to bed file!
  currCmdGenResult = x[[2]]
  #The -tag NH reports the number of hits in the genome for each read!
  bamToBed_CLI_cmdRes = generateCommandResult(BamToBed_CLI(inFilePath = getInFilePath(getCLIApplication(currCmdGenResult)),
                                                           inFileNames = getOutResultName(getOutResultReference(currCmdGenResult)), 
                                                           cliParams = "-tag NH", outputFlag = "", 
                                                           outFilePath = getInFilePath(getCLIApplication(currCmdGenResult)),outputFormat = "bed"))
  
  #Then merge the resulting bed files into contigs
  #This produces a bed file with 6 columns: chr 1 3 readCount Average(MAPQ) strand
  
  mergeBedFile_CLI_cmdRes = generateCommandResult(MergeBedFile_CLI(inFilePath = getInFilePath(getCLIApplication(currCmdGenResult)), 
                                                                   inFileNames = getOutResultName(getOutResultReference(bamToBed_CLI_cmdRes)), 
                                                                   cliParams = paste0("-s -d ",readOverlap_contig," -c 4,5,6 -o count,mean,distinct | awk '{if($4 > ",read_threshold,") print }'"), 
                                                                   outputFlag = "_contigs", outFilePath = getInFilePath(getCLIApplication(currCmdGenResult))) 
  )
  
  resL = list(bamToBed_CLI_cmdRes,mergeBedFile_CLI_cmdRes)
  return( resL )
})

for( i in 1:length(bamToBedAndMergeL)  ){
  for( j in 1:length(bamToBedAndMergeL[[i]])){
    tmpCommandLog = c(tmpCommandLog, getCommandLog(bamToBedAndMergeL[[i]][[j]]) )
  }
}

#############################################
#     Contig Assembly - among samples - MultiIntersectBed
#############################################
inFP = sapply( bamToBedAndMergeL, function(x){getInFilePath(getCLIApplication(x[[2]]))})
inFN = sapply( bamToBedAndMergeL, function(x){getOutResultName(getOutResultReference(x[[2]]))})

############ deprecated ##########
# multiInter_CLI = MultiIntersectBed_CLI(inFilePath = inFP, inFileNames = inFN, cliParams = "-cluster",outputFlag = "Contigs_multiInter", 
#                                        outFilePath = contigAssemblyDir, sortInFiles = TRUE, requireStrandness=TRUE, strandColumn = 6, filterConditions=samplesInfo$condition)
# multiInter_CLI_cmdRes = generateCommandResult(multiInter_CLI)
# tmpCommandLog = c(tmpCommandLog, getCommandLog(multiInter_CLI_cmdRes) )
############ deprecated ##########

tmpCommandLog = c(tmpCommandLog, paste0("\nmkdir ",contigAssemblyDir,"\n") )
multiIntersectBed_perl_CLI = MultiIntersectBed_perl_CLI(inFilePath = inFP, inFileNames = inFN, outputFlag = "", outFileName = "multiIntersectClust",outFilePath = contigAssemblyDir)
multiIntersectBed_perl_CLI_cmdRes = generateCommandResult(multiIntersectBed_perl_CLI)
tmpCommandLog = c(tmpCommandLog, getCommandLog(multiIntersectBed_perl_CLI_cmdRes) )

##################################################
# Multi intersect bed with -cluster option returns minimal overlap of all contigs! -> I think this is best suited for comparing many samples -> otherwise gets vague and 
# one has to adjust for length of the contigs when counting reads! -> So now the reads have to overlap at least 80% for counting I would say!
# Otherwise it will report all overlaps and can then be merged by MergeBed
###################################################
#Now we have to intersect these contigs with the reads of the bed (from bam) file and perform the same merge with the counting and reporting of unique reads
#Should be something like this -> this keeps the ID column from the multiIntersectBed!
#NOW Intersecting all -> like this command below and parse these files then they are ready for clustering and DE (one may also count again against clusters by removing dups?!)
# intersectBed -s -f 0.8 -wb -a ../bowtieMappings/test1_trimmed_bowtieOut/test1_trimmed_mapped_s.bed -b multiIntersectClust.bed | mergeBed -s -c 4,5,6,10 -o count,mean,distinct,distinct -i stdin > sameAsmultiIntersectClust.bed

#Then Merge them in R: 
intersectBed_CLI_cmdResL = lapply( bamToBedAndMergeL, function(x){
  curr = x[[1]]
  intersectBed_CLI = IntersectBed_CLI(
    inFilePath = c( getOutFilePath(getCLIApplication(curr)), getOutFilePath(getCLIApplication(multiIntersectBed_perl_CLI_cmdRes))),
    inFileNames = c(getOutResultName(getOutResultReference(curr)), getOutResultName(getOutResultReference(multiIntersectBed_perl_CLI_cmdRes))),
    cliParams = "-s -wb", outputFlag = "_contigsMulti",
    outFilePath = getOutFilePath(getCLIApplication(curr)),
    outFileName = sub(".bed","",getOutResultName(getOutResultReference(curr))))
  
  intersectBed_CLI_cmdRes = generateCommandResult(intersectBed_CLI)
  return(intersectBed_CLI_cmdRes)
})

tmpCommandLog = c(tmpCommandLog, sapply(intersectBed_CLI_cmdResL, getCommandLog) )


mergeBedFile_CLI_cmdResL = lapply(intersectBed_CLI_cmdResL, function(x){
  mergeBedFile_CLI = MergeBedFile_CLI(inFilePath = getOutFilePath(getCLIApplication(x)), 
                                      inFileNames = getOutResultName(getOutResultReference(x)), 
                                      cliParams = "-s -c 4,5,6,10 -o count,mean,distinct,distinct", 
                                      outputFlag = "_merged", 
                                      outFilePath = getOutFilePath(getCLIApplication(x)))  
  mergeBedFile_CLI_cmdRes = generateCommandResult(mergeBedFile_CLI)
  return(mergeBedFile_CLI_cmdRes)
})

tmpCommandLog = c(tmpCommandLog, sapply(mergeBedFile_CLI_cmdResL, getCommandLog) )

######################################
#       Write commands to file
######################################
cat(tmpCommandLog,  file = commandLog, append = FALSE) 

######################################
#       Execute commands!
######################################
cmdExecTime = proc.time()
setwd(rootDir)
system(paste0("bash ", commandLog, " > ", executionLog ))
cmdExecTime = proc.time() - cmdExecTime
cmdExecTime
################################################
#       Contig clustering
################################################
#FIRST Cluster function definition

################################
#     Now combined these files -> in R currently
################################
require(rtracklayer)
generateCountUnclusteredTable = function( contigFile, sampleBedFiles ){
  overallContigPath = file.path( getOutFilePath(getCLIApplication(contigFile)), getOutResultName(getOutResultReference(contigFile)) )
  overallContig = read.table(overallContigPath, sep="\t", header=FALSE)
  colnames(overallContig)  = c("chr","start","end","ID","score","strand")
  overallContig = overallContig[order(overallContig$ID),]
  overallContig = overallContig[,c("chr","start","end","ID","strand")]
  
  mergedDF = do.call( cbind, mapply( function(x, y){
    mergedFilePath = file.path( getOutFilePath(getCLIApplication(x)), getOutResultName(getOutResultReference(x)) )
    
    mergedFile = read.table(mergedFilePath, sep="\t", header=FALSE)
    colnames(mergedFile)  = c("chr","start","end","ReadCount","Uniqueness","strand","ID")
    
    mergedTable = merge(x = overallContig, y = mergedFile[,c("ReadCount","Uniqueness","ID")], by = "ID", all.x = TRUE)
    mergedTable = mergedTable[order(mergedTable$ID),]
    
    mergedTable = mergedTable[,c("ReadCount","Uniqueness")]
    colnames( mergedTable ) = paste0(colnames( mergedTable ),"_",y)
    
    return(mergedTable)
  }, sampleBedFiles, samplesInfo$sampleName, SIMPLIFY = FALSE ))
  
  overallContig = cbind(overallContig, mergedDF)
  return(overallContig)
}

unclDF = generateCountUnclusteredTable(contigFile=multiIntersectBed_perl_CLI_cmdRes, sampleBedFiles=mergeBedFile_CLI_cmdResL)
colnames(unclDF) = gsub(".*\\.","",colnames(unclDF) )

#Now This table is reduced to mean ReadCount and Mean Uniqueness and can then be used for clustering!
toClusterDF = unclDF[,1:5]
toClusterDF$readCount = rowMeans(  unclDF[,which(c(rep(FALSE,5),c(6:dim(unclDF)[2]%%2 == 0)))] ,na.rm = TRUE)
toClusterDF$uniqueness = rowMeans(  unclDF[,which(c(rep(FALSE,5),c(6:dim(unclDF)[2]%%2 == 1)))], na.rm = TRUE )

##################################################
#     Short sanity checking ... should not be necessary however
finiteReadCount = which( is.finite(toClusterDF$readCount) )
if( length(finiteReadCount) != dim(toClusterDF)[1] ){warning("Some contigs do not contain ANY overlapping read -> check!"); toClusterDF = toClusterDF[finiteReadCount,]}
finiteUniq = which( is.finite(toClusterDF$uniqueness) )
if( length(finiteUniq) != dim(toClusterDF)[1] ){warning("Some contigs do not contain ANY overlapping read -> check!"); toClusterDF = toClusterDF[finiteUniq,]}
###################################################
contigForCountingGR = with(toClusterDF, GRanges(seqnames=chr, IRanges(start=start, end=end), strand = strand, readCount = readCount, 
                                                uniqueness =uniqueness))

contigForCountingGR_unique = contigForCountingGR[contigForCountingGR$uniqueness <= 1]
contigForCountingGR_unclustered = contigForCountingGR[contigForCountingGR$uniqueness > 1]

allReads = import(file.path(getOutFilePath(getCLIApplication(bamToBedAndMergeL[[1]][[1]])), 
                            getOutResultName(getOutResultReference(bamToBedAndMergeL[[1]][[1]]))), 
                  format="BED", asRangedData=FALSE)#alignedbwt_s.bed
start(allReads) = start(allReads) - 1 #Strange behaviour in import of rtracklayer!


clusteringTime = proc.time()

clusteredContigs = clusterMultiMappingReads_stringent(contigForCountingGR_unclustered = contigForCountingGR_unclustered,allReads = allReads)
contigForCountingGR_clustered = append(contigForCountingGR_unique, clusteredContigs)
contigForCountingGR_clustered = reduce(contigForCountingGR_clustered)
#get information back!
contigForCountingGR_clustered = subsetByOverlaps( contigForCountingGR, contigForCountingGR_clustered)

clusteringTime = proc.time() - clusteringTime
clusteringTime

clusteredContigsFN = "clusteredContigs.bed"

elementMetadata(contigForCountingGR_clustered)$name = paste0("contig",1:length(contigForCountingGR_clustered))
elementMetadata(contigForCountingGR_clustered)$score = elementMetadata(contigForCountingGR_clustered)$uniqueness
#Writing as bed file (when using export function of rtracklayer -> coordinates get -1?!?)
write.table(data.frame( 
  as.character(seqnames(contigForCountingGR_clustered)), 
  start(contigForCountingGR_clustered), 
  end(contigForCountingGR_clustered), 
  contigForCountingGR_clustered$name, 
  contigForCountingGR_clustered$uniqueness, 
  as.character(strand(contigForCountingGR_clustered))
)
, file.path(contigClusterDir,clusteredContigsFN), sep="\t", quote=FALSE,col.names=FALSE, row.names=FALSE)

##################################################
#         Read counting
##################################################
commandLogCounting = c()
#Definition of the output directories from bowtie first, so they can be iterrated

commandLogCounting = c(commandLogCounting, paste0( "\nmkdir ", readCountsDir ,"\n" ) )

inFN = file.path( sapply( bowtieCLI_cmdResL,function(x){getOutFilePath(getCLIApplication(x))} ),
                  sapply( samToolsHTSeqCmdL, function(x){getOutResultName(getOutResultReference(x[[2]]))} ) )

multiBamCov_CLI = MultiBamCov_CLI(inFilePath = "", 
                                  inFileNames = inFN,
                                  cliParams =  c("-s -f 0.05 -D -q 20"), # -D include duplicated reads 0.05 ~1bp at length 200 - 10 bp -> only good for short read data!
                                  outputFlag = "_multibamcov", 
                                  outFilePath = readCountsDir,
                                  annotationFileMB = file.path(contigClusterDir,clusteredContigsFN), 
                                  annotationType = "bed")

multiBamCov_CLI_cmdRes = generateCommandResult( object = multiBamCov_CLI )
commandLogCounting = c(commandLogCounting, getCommandLog(multiBamCov_CLI_cmdRes) )

setwd(rootDir)
sapply( commandLogCounting, system)


rm(allReads)#to save memory!
save.image(file = "MappingAssemblyClustering.rda")

# readCountTable = read.table( file.path(getOutFilePath(getCLIApplication(multiBamCov_CLI_cmdRes)), getOutResultName(getOutResultReference(multiBamCov_CLI_cmdRes)) ), sep="\t", header=FALSE )

##################################################
#     Annotate with R package
##################################################




##################################################
#     Annotate with infernal
##################################################

# #PUT RFAM INTO R PACAKGE -> ONLY 
# ##### infernal runs
# #cmscan --tblout tmp -E 20 --noali Rfam.cm.1_1 fastaFiles/testSeq.fa #NICELY PARSEABLE FILE -> PERL SCRIPT!








# 
# 
# ###################################################
# #     Now get back the read numbers and uniqueness values of the initial clusters
# ###################################################
# 
# #TODO think about that!!!
# 
# ############################### Now combine all bed files! ##########################
# #x = contigWithinGroupL[[1]]
# inFPcontigs = sapply( contigWithinGroupL, function(x){
#   return( getOutFilePath(getCLIApplication(x[[2]])))
# } )
# inFNcontigs = sapply( contigWithinGroupL, function(x){
#   return( getOutResultName(getOutResultReference(x[[2]])))
# } )
# 
# #Now creating a directory for producting merged contigs!
# IntersectBed_CLI_contigs_cmdRes = generateCommandResult(IntersectBed_CLI(inFilePath = inFPcontigs, 
#                                                                          inFileNames = inFNcontigs, 
#                                                                          cliParams = "-s -wb", outputFlag = "", 
#                                                                          outFilePath = file.path(rootDir, "contigsForCounting"), outFileName = "allIntersection"))
# 
# 
# cat(getCommandLog(IntersectBed_CLI_contigs_cmdRes), file = commandLog, append = TRUE  )   
# 
# #Ok now we obtained a file containing 13 columns (column 7 represents default filename index, else one can use -filenames)
# #Therefore we merge this file to create a merged contig file!
# # mergeBed  -i stdin -s -d -3 -c 4,5,6,11,12 -o mean,mean,distinct,mean,mean #meanReadCounts A, mean MAPQ A, strand distinction, mean ReadCount B..., mean MAPQ B...
# #meanReadCounts A, mean MAPQ A, strand distinction, mean ReadCount B..., mean MAPQ B...
# mergeBedFile_CLI_contig_cmdRes = generateCommandResult(MergeBedFile_CLI(inFilePath = getOutFilePath(getCLIApplication(IntersectBed_CLI_contigs_cmdRes)), 
#                                                                         inFileNames = getOutResultName(getOutResultReference(IntersectBed_CLI_contigs_cmdRes)), 
#                                                                         cliParams = paste0("-s -d -3 -c 4,5,6,11,12 -o mean,mean,distinct,mean,mean"), 
#                                                                         outputFlag = "_contigs", 
#                                                                         outFilePath = getOutFilePath(getCLIApplication(IntersectBed_CLI_contigs_cmdRes))) 
# )
# 
# cat(getCommandLog(mergeBedFile_CLI_contig_cmdRes), file = commandLog, append = TRUE  )
# 
# 
# #################################### star MAPPER: ##################################################
# # Hi Praful,
# # 
# # we are routinely using STAR to map "small RNA" (~<200b) data within the ENCODE project - the miRNA (mostly mature) are a major subclass of these small RNA.
# # We are using STAR with the following parameters:
# #   --outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin 16 --outFilterScoreMinOverLread 0  --outFilterMatchNminOverLread 0 --alignIntronMax 1
# # (>=16b matched to the genome, number of mismatches <= 5% of mapped length, i.e. 0MM for 16-19b, 1MM for 20-39b etc, splicing switched off).
# # 
# # You can clip 3' adapter before feeding the reads to STAR, or you can use simple built-in clipper
# # --clip3pAdapterSeq TGGAATTCTC      --clip3pAdapterMMp 0.1
# # (second parameter is the proportion of mismatches in the matched adapter length).
# # 
# # You would also likely want to filter out reads that STAR "genomically" trims at the 5' (see the discussion about "Soft clipping" here).
# # This simple awk script will filter out all alignments that are trimmed by more than 1 base from the 5'.
# # awk '{S=0; split($6,C,/[0-9]*/); n=split($6,L,/[NMSID]/);  if (and($2,0x10)>0 && C[n]=="S") {S=L[n-1]} else if (and($2,0x10)==0 && C[2]=="S") {S=L[1]}; if (S<=1) print }' Aligned.out.sam > Aligned.filtered.sam
# 
# 
# #i Alex, I can confirm that your warnings are correct - the softclipping will allow the untrimmed adapter sequence to align to the genome if there is enough sequence similarity, resulting in more "unique" alignments.
# #If the trimming is done well (i've tried cutadapt, alientrimmer, and reaper from kraken tools - so far the best one) then the results are comparable between STAR and bowtie. Bowtie is quite common in miRNA pipelines.. 
# #If the trimming is not so accurate (say it missed one or two bases of the adapter), STAR seems to be a little more robust in these cases, as when I've used Cutadapt or Alientrimmer.
# 
# 
# 

# #############################################
# #     Initial Read Counting for clustering!
# #############################################
# 
# ######################
# # MultiBamCoverageCount!
# ######################
# tmpCommandLog = c(tmpCommandLog, paste0( "\nmkdir ", contigClusterDir ,"\n" ) )
# 
# multiIntersectionInFN = file.path( getOutFilePath(getCLIApplication(multiInter_CLI_cmdRes)),
#                                    paste0(getOutputFlag(getCLIApplication(multiInter_CLI_cmdRes)),"_MultiBamIn.bed" ))
# 
# # Modification of the multiIntersect output, so that it is representing the BED file standard (strand in column 6!)
# modifyForCoverageCmd = paste0( "\ncat ",file.path( getOutFilePath(getCLIApplication(multiInter_CLI_cmdRes)),getOutResultName(getOutResultReference(multiInter_CLI_cmdRes)) )
#                                ," | awk \'BEGIN {OFS=\"\\t\"} {print $1,$2,$3,$4,$5,$NF}\' > ",  multiIntersectionInFN,"\n" )
# 
# tmpCommandLog = c(tmpCommandLog,modifyForCoverageCmd)
# 
# multiBamCovInFN = sapply(  samToolsHTSeqCmdL, function(x){
#   curr = x[[3]]
#   file.path( getOutFilePath(getCLIApplication(curr)),getOutResultName(getOutResultReference(curr)) )
# })
# 
# multiBamCovClust_CLI = MultiBamCov_CLI(inFilePath = "", 
#                                        inFileNames = multiBamCovInFN,
#                                        cliParams =  c("-s","-D"), 
#                                        outputFlag = "_multibamcov", 
#                                        outFilePath = contigClusterDir,
#                                        annotationFileMB = multiIntersectionInFN, 
#                                        annotationType = "bed")
# multiBamCovClust_CLI_cmdRes = generateCommandResult( object = multiBamCovClust_CLI )
# 
# 
# tmpCommandLog = c(tmpCommandLog, getCommandLog(multiBamCovClust_CLI_cmdRes) )
