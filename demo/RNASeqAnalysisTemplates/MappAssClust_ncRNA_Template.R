
# running the main scripts with arguments: 
# R CMD BATCH --no-save --no-restore '--args /tmp/test' test.R, or Rscript test.R /tmp/test 
# the second command gives only output of R, the first one outputs all code too

options(stringsAsFactors = FALSE)

#####################################
#     Loading configuration file
#####################################

#If no command line arguments are given search for it in the current directory
args <- commandArgs(TRUE)
if( length(args) == 0 ){
  load( file.path(getwd(),"Configuration.rda") )
} else{
  load(file.path(args[1],"Configuration.rda"))
}

#####################################
# Loading the required libraries
#####################################
setwd(rootDir)
library(CLIHelperPackage)
library(RNASeqUtility)
####################################################################################################################################################
#                                                             Start Analysis
####################################################################################################################################################
tmpCommandLog = ""
#######################
#   Cutadapt to trim the fastq files
#######################
cutAdaptCLI = Cutadapt_CLI(inFilePath=rawDataDir, cliParams =  cutadaptOptions, outputFlag = "_trimmed", 
                           outFilePath = file.path(rootDir,"rawDataTrimmed") )
cutAdatptCLI_cmdRes = generateCommandResult(object = cutAdaptCLI )
#logging
tmpCommandLog = getCommandLog(cutAdatptCLI_cmdRes)


#############################################################
#     Mapping ncRNAs to the ENSEMBL ncRNA fasta file
#############################################################
tmpCommandLog = c(tmpCommandLog, paste0("\nmkdir ", ncRNAmappingDir,"\n"))
fastQFiles = getOutResultName(getOutResultReference(cutAdatptCLI_cmdRes))

mappingncRNACLI_cmdResL = mapply( function(fqf, samplePrefix){
  
  outFN = sub(".fastq","",fqf)
  outFP = file.path( ncRNAmappingDir, outFN )
  
  mappingCLI = RNAStar_CLI(inFilePath = getOutFilePath(getCLIApplication(cutAdatptCLI_cmdRes)), 
                           inFileNames = fqf, 
                           cliParams = rnaStarncRNA_params, 
                           outputFlag = paste0(outFN,"_"), #prefix -> must be distinguishable -> best to use the fastQFileNames ,or the sample Short Names!!
                           outFilePath = outFP, 
                           genomeIndexFilePath = genomeIndexFilePath_ncRNA,
                           filterSam = rnaStarncRNA_filterSam, 
                           outputFormat = rnaStarncRNA_outputFormat)
  
  return( generateCommandResult(mappingCLI) )
}, fastQFiles, samplesInfo$sampleName, SIMPLIFY=FALSE )

names(mappingncRNACLI_cmdResL) = fastQFiles

#logging
for( i in 1:length(mappingncRNACLI_cmdResL)  ){
  tmpCommandLog = c(tmpCommandLog, getCommandLog(mappingncRNACLI_cmdResL[[i]]) )
}

#############################################################
#     Read counting by bamToBed and mergeBed
#############################################################
bamToBedAndMergencRNAL = lapply( mappingncRNACLI_cmdResL, function(x){
  #converting bam to bed file!
  currCmdGenResult = x
  #The -tag NH reports the number of hits in the genome for each read!
  bamToBed_CLI_cmdRes = generateCommandResult(BamToBed_CLI(inFilePath = getOutFilePath(getCLIApplication(currCmdGenResult)),
                                                           inFileNames = getOutResultName(getOutResultReference(currCmdGenResult)), 
                                                           cliParams = "", outputFlag = "", 
                                                           outFilePath = getOutFilePath(getCLIApplication(currCmdGenResult)),outputFormat = "bed"))
  
  #Then merge the resulting bed files into contigs
  #This produces a bed file with 6 columns: chr 1 3 readCount Average(MAPQ) strand
  
  mergeBedFile_CLI_cmdRes = generateCommandResult(MergeBedFile_CLI(inFilePath = getOutFilePath(getCLIApplication(currCmdGenResult)), 
                                                                   inFileNames = getOutResultName(getOutResultReference(bamToBed_CLI_cmdRes)), 
                                                                   cliParams = paste0("-s -d ",readOverlap_contig," -c 4,5,6 -o count,mean,distinct"), 
                                                                   outputFlag = "_counted", outFilePath = getOutFilePath(getCLIApplication(currCmdGenResult))) 
  )
  
  resL = list(bamToBed_CLI_cmdRes,mergeBedFile_CLI_cmdRes)
  return( resL )
})


for( i in 1:length(bamToBedAndMergencRNAL)  ){
  for( j in 1:length(bamToBedAndMergencRNAL[[i]])){
    tmpCommandLog = c(tmpCommandLog, getCommandLog(bamToBedAndMergencRNAL[[i]][[j]]) )
  }
}

##########################################################################################################################
#                   Mapping and assembly of non-canonical-small ncRNAs
##########################################################################################################################

#################################################
#     RNA Star Mapping
#################################################
#logging
tmpCommandLog = c(tmpCommandLog, paste0("\nmkdir ", mappingDir,"\n"))

mappingCLI_cmdResL = mapply( function(x, samplePrefix){
  
    outFN =  paste0(sub( ".fastq$","", getInFileNames(getCLIApplication(x))), "_remapped")
    outFP = file.path( mappingDir, outFN )
    inFN = paste0(getOutputFlag(getCLIApplication(x)),"Unmapped.out.mate1")
    
  mappingCLI = RNAStar_CLI(inFilePath = getOutFilePath(getCLIApplication(x)), 
                           inFileNames = inFN, 
                           cliParams = rnaStarGenome_params,
                           outputFlag = outFN, #prefix -> must be distinguishable -> best to use the fastQFileNames ,or the sample Short Names!!
                           outFilePath = outFP, 
                           genomeIndexFilePath = genomeIndexFilePath,
                           filterSam = rnaStarGenome_filterSam, 
                           outputFormat = rnaStarGenome_outputFormat)
  #"--runThreadN 6 --outFilterMismatchNmax 1 --outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin 16 --outFilterScoreMinOverLread 0  --outFilterMatchNminOverLread 0 --alignIntronMax 1 --outFilterMultimapNmax 100"
  #"--runThreadN 8 --outFilterMismatchNoverReadLmax 0.023 --outFilterMatchNmin 18 --outFilterScoreMinOverLread 0  --outFilterMatchNminOverLread 0 --alignIntronMax 1 --outFilterMultimapNmax 100 --alignEndsType EndToEnd --outSAMprimaryFlag AllBestScore"
  
  return( generateCommandResult(mappingCLI) )
}, mappingncRNACLI_cmdResL, samplesInfo$sampleName, SIMPLIFY=FALSE )

#logging
for( i in 1:length(mappingCLI_cmdResL)  ){
  tmpCommandLog = c(tmpCommandLog, getCommandLog(mappingCLI_cmdResL[[i]]) )
}

######################
#   Samtools commands sorting by Name and coordinates _sn and _s
######################

samToolsHTSeqCmdL = lapply( mappingCLI_cmdResL, function(x){
  
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

if(grouping){
  groupVect = 1:length(samplesInfo$condition)
  names(groupVect) = unique(samplesInfo$condition)
  groupVect = groupVect[samplesInfo$condition]
  multiIntersectBed_perl_CLI = MultiIntersectBed_perl_CLI(inFilePath = inFP, inFileNames = inFN, outputFlag = "",withinGroupTH = withinGroupTH, groupVect = groupVect, outFileName = "multiIntersectClust",outFilePath = contigAssemblyDir)
} else{
  multiIntersectBed_perl_CLI = MultiIntersectBed_perl_CLI(inFilePath = inFP, inFileNames = inFN, outputFlag = "",withinGroupTH = withinGroupTH, outFileName = "multiIntersectClust",outFilePath = contigAssemblyDir)
    
}

tmpCommandLog = c(tmpCommandLog, paste0("\nmkdir ",contigAssemblyDir,"\n") )
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
                                      cliParams = paste0("-s -d ",readOverlap_contig," -c 4,5,6,10 -o count,mean,distinct,distinct"), 
                                      outputFlag = "_merged", 
                                      outFilePath = getOutFilePath(getCLIApplication(x)))  
  mergeBedFile_CLI_cmdRes = generateCommandResult(mergeBedFile_CLI)
  return(mergeBedFile_CLI_cmdRes)
})

tmpCommandLog = c(tmpCommandLog, sapply(mergeBedFile_CLI_cmdResL, getCommandLog) )


#######################################
#     Finally intersect the contig with the bed files of the individual experiments
#     This creates bed files containing all the reads only for the contig regions, to minimize computation
#     One file from each group is chosen for clustering
#######################################
groupsTmp = unique(samplesInfo$condition)
subVect = sapply(groupsTmp, function(x){which(samplesInfo$condition == x )[1]})

intersectBed_CLI_clust_cmdResL = lapply( bamToBedAndMergeL[subVect], function(x){
  curr = x[[1]]
  intersectBed_CLI = IntersectBed_CLI(
    inFilePath = c( getOutFilePath(getCLIApplication(curr)), getOutFilePath(getCLIApplication(multiIntersectBed_perl_CLI_cmdRes))),
    inFileNames = c(getOutResultName(getOutResultReference(curr)), getOutResultName(getOutResultReference(multiIntersectBed_perl_CLI_cmdRes))),
    cliParams = "-s", outputFlag = "_forClustering",
    outFilePath = getOutFilePath(getCLIApplication(multiIntersectBed_perl_CLI_cmdRes)),
    outFileName = sub(".bed","",getOutResultName(getOutResultReference(curr))))
  
  intersectBed_CLI_cmdRes = generateCommandResult(intersectBed_CLI)
  return(intersectBed_CLI_cmdRes)
})

tmpCommandLog = c(tmpCommandLog, sapply(intersectBed_CLI_clust_cmdResL, getCommandLog) )


###################################################
# Concatenate both files
###################################################
readsForClustering = "readsForClustering.bed"
readsForClustering_cmdL = lapply(intersectBed_CLI_clust_cmdResL, function(x){
  if( getOutResultName(getOutResultReference(x)) == getOutResultName(getOutResultReference(intersectBed_CLI_clust_cmdResL[[1]])) ){
    return( paste0( "\ncat ", file.path( getOutFilePath(getCLIApplication(x)), getOutResultName(getOutResultReference(x))), " > ", file.path(contigClusterDir,readsForClustering),"\n") )
  } else{
    return( paste0( "\ncat ", file.path( getOutFilePath(getCLIApplication(x)), getOutResultName(getOutResultReference(x))), " >> ", file.path(contigClusterDir,readsForClustering),"\n") )
  }
})

tmpCommandLog = c(tmpCommandLog, unlist(readsForClustering_cmdL) )


######################################
#       Write commands to file
######################################
cat(tmpCommandLog,  file = commandLog, append = FALSE) 
save.image(file.path(rootDir,"CommandsBeforeClustering.rda"))
######################################
#       Execute commands!
######################################
cmdExecTime = proc.time()
setwd(rootDir)
system(paste0("bash ", commandLog, " > ", executionLog ))
cmdExecTime = proc.time() - cmdExecTime
cmdExecTime


########################################################################################################################################################
#                                                                       Contig clustering
########################################################################################################################################################
#FIRST Cluster function definition


################################
#     Now combined these files -> in R currently
################################
require(rtracklayer)
require(RNASeqUtility)

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
#Delete all sequences that are shorter than 18 nt
contigForCountingGR = contigForCountingGR[width(contigForCountingGR) >= 18]


contigForCountingGR_unique = contigForCountingGR[contigForCountingGR$uniqueness <= 1]
contigForCountingGR_unclustered = contigForCountingGR[contigForCountingGR$uniqueness > 1]

allReads = import(file.path( contigClusterDir ,readsForClustering), 
                  format="BED", asRangedData=FALSE)#alignedbwt_s.bed
start(allReads) = start(allReads) - 1 #Strange behaviour in import of rtracklayer!


clusteringTime = proc.time()

clusteredContigs = clusterMultiMappingReads_stringent(contigForCountingGR_unclustered = contigForCountingGR_unclustered,allReads = allReads, readCompositionIdentity = readCompositionIdentity)
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

inFN = file.path( sapply( mappingCLI_cmdResL,function(x){getOutFilePath(getCLIApplication(x))} ),
                  sapply( samToolsHTSeqCmdL, function(x){getOutResultName(getOutResultReference(x[[2]]))} ) )

multiBamCov_CLI = MultiBamCov_CLI(inFilePath = "", 
                                  inFileNames = inFN,
                                  cliParams =  c("-s -f 0.05 -D"), # -q 20 quality values not supported rna-star output! -D include duplicated reads 0.05 ~1bp at length 200 - 10 bp -> only good for short read data! 
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

