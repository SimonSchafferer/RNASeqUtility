
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
#
#   General structure: 
#   application specific objects containing commands for execution are created e.g. cutAdatptCLI_cmdRes
#   All commands are then written to a command file for documentation and this file is then executed within R
#
####################################################################################################################################################
tmpCommandLog = ""
#######################
#   Cutadapt to trim the fastq files cutadaptOptions are defined in the configuration file
#######################
runCutAdapt = TRUE
bedTools225 = TRUE
if(bedTools225){
  bedTools225Sed = " | sed -r 's/(\\s+)?\\S+//4'"
} else{
  bedTools225Sed = ""
}


if(runCutAdapt){
  cutAdaptCLI = Cutadapt_CLI(inFilePath=rawDataDir, cliParams =  cutadaptOptions, outputFlag = "_trimmed", 
                             outFilePath = file.path(rootDir,"rawDataTrimmed") )
  cutAdatptCLI_cmdRes = generateCommandResult(object = cutAdaptCLI )
  #logging
  tmpCommandLog = getCommandLog(cutAdatptCLI_cmdRes)
  
  #Execute Cutadapt
  cutAdatptCLI_execRes = executeCommandResult(cutAdatptCLI_cmdRes, testing=FALSE)
}


#############################################################
#     Mapping ncRNAs to the ENSEMBL ncRNA fasta file
#       ENSEMBL ncRNA file (see configuration file)
#       mapping parameters are defined in configuration file
#       Important parameters are multi mappings allowed (since the ENSEMBL ncRNA file was also enriched by adding miRBase entries, therefore duplicated annotations are induced)
#       Also important: The unmapped fastq file has to written, in order to proceed with the contig assembly and clustering.
#############################################################
tmpCommandLog = c(tmpCommandLog, paste0("\nmkdir ", ncRNAmappingDir,"\n"))
fastQFiles = getOutResultName(getOutResultReference(cutAdatptCLI_cmdRes))

mappingncRNACLI_cmdResL = mapply( function(fqf, samplePrefix){
  
  outFN = sub(".fastq","",fqf)
  outFP = file.path( ncRNAmappingDir, outFN )
  
  if(runCutAdapt){
    inFilePath = getOutFilePath(getCLIApplication(cutAdatptCLI_cmdRes))
  } else{
    inFilePath = rawDataDir
  }
  
  mappingCLI = RNAStar_CLI(inFilePath = inFilePath, 
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
#Execution
mappingncRNACLI_execResL = lapply( mappingncRNACLI_cmdResL, function(x){
  executeCommandResult(mappingncRNACLI_cmdResL[[i]],  testing=FALSE)
})


#############################################################
#     Read counting by bamToBed and mergeBed
#     Writing bam files to bed files by extracting the <NH> tag (number of multiple mappings) from the bam file. 
#     The bed file (containing all reads) is then merged into contigs, whereby the reads are counted 
#     and the mapping uniqueness is reported by calculating the mean of the <NH> value from all reads in a contig
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
                                                                   cliParams = paste0("-s -d ",readOverlap_contig," -c 4,5,6 -o count,mean,distinct", bedTools225Sed), 
                                                                   outputFlag = "_counted", outFilePath = getOutFilePath(getCLIApplication(currCmdGenResult))) 
  )
  
  resL = list(bamToBed_CLI_cmdRes,mergeBedFile_CLI_cmdRes)
  return( resL )
})


bamToBedAndMerge_execL = list()
for( i in 1:length(bamToBedAndMergencRNAL)  ){
  for( j in 1:length(bamToBedAndMergencRNAL[[i]])){
    tmpCommandLog = c(tmpCommandLog, getCommandLog(bamToBedAndMergencRNAL[[i]][[j]]) )
    bamToBedAndMerge_execL = c( bamToBedAndMerge_execL, executeCommandResult(bamToBedAndMergencRNAL[[i]][[j]],  testing=FALSE) )
  }
}

##########################################################################################################################
#                   Mapping and assembly of potentially novel ncRNAs
##########################################################################################################################

#################################################
#     RNA Star Mapping
#     Mapping previosly unmapped reads to an indexed genome
#################################################
#logging
tmpCommandLog = c(tmpCommandLog, paste0("\nmkdir ", mappingDir,"\n"))
dir.create(mappingDir)
# rnaStarGenome_params = paste0(rnaStarGenome_params," -v 30620442512")
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

mappingWhole_execL = lapply(mappingCLI_cmdResL, function(x){
  executeCommandResult(x, testing=FALSE)
})
######################
#   Samtools commands sorting by Name and coordinates _sn and _s
#   Conversions for later use: 
#     sort bam by name
#     index bam by name
#     sort bam by coordinates
#     index bam by coordinates
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
samToolsHTSeq_execL = list()
for( i in 1:length(samToolsHTSeqCmdL)  ){
  for( j in 1:length(samToolsHTSeqCmdL[[i]])){
    tmpCommandLog = c(tmpCommandLog, getCommandLog(samToolsHTSeqCmdL[[i]][[j]]) )
    samToolsHTSeq_execL = c( samToolsHTSeq_execL, executeCommandResult(samToolsHTSeqCmdL[[i]][[j]],  testing=FALSE) )
  }
}

#############################################
#     Contig Assembly in each sample
#
#     Read counting by bamToBed and mergeBed
#     Writing bam files to bed files by extracting the <NH> tag (number of multiple mappings) from the bam file. 
#     The bed file (containing all reads) is then merged into contigs, whereby the reads are counted 
#     and the mapping uniqueness is reported by calculating the mean of the <NH> value from all reads in a contig
#     Only contigs are kept that are superseding a given threshold (as defined in configuaration file)
#############################################################
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
                                                                   cliParams = paste0("-s -d ",readOverlap_contig," -c 4,5,6 -o count,mean,distinct", bedTools225Sed, " | awk '{if($4 > ",read_threshold,") print }'"), 
                                                                   outputFlag = "_contigs", outFilePath = getInFilePath(getCLIApplication(currCmdGenResult))) 
  )
  
  resL = list(bamToBed_CLI_cmdRes,mergeBedFile_CLI_cmdRes)
  return( resL )
})

bamToBedAndMerge_execL = list()
for( i in 1:length(bamToBedAndMergeL)  ){
  for( j in 1:length(bamToBedAndMergeL[[i]])){
    tmpCommandLog = c(tmpCommandLog, getCommandLog(bamToBedAndMergeL[[i]][[j]]) )
    bamToBedAndMerge_execL = c( bamToBedAndMerge_execL, executeCommandResult(bamToBedAndMergeL[[i]][[j]],  testing=FALSE) )
  }
}

#############################################
#     Contig Assembly - between samples - MultiIntersectBed
#
#     MultiIntersectBed_perl_CLI: 
#     
#     When no grouping is defined: 
#       Contigs must be present in n-withinGroupTH samples to get reported (this makes sure that contigs are not reported if present e.g. in only one sample -> then reead counts in all others would be 0)
#     When grouping is defined: 
#       Contigs must be present in n-withinGroupTH samples within a group (e.g. use case: if one contig is only present in WT but not in TG it will be reported when grouping is defined)
#       Special case grouping is defined but the groupVector is set to one group only, then contigs get reported when it is present in n-withinGroupTH samples 
#
#     CAVE: Interval reporting is different when grouping/no-grouping is defined
#
#     no grouping: (calls Multi intersect bed with -cluster option)
#     ----------          contig A
#         -----------     contig B
#         ------          reported contig
#
#     grouping:
#     ----------          contig A
#         -----------     contig B
#     ---------------     reported contig (currently default)
#
#     grouping is defined in configuration file
#     
#############################################
inFP = sapply( bamToBedAndMergeL, function(x){getInFilePath(getCLIApplication(x[[2]]))})
inFN = sapply( bamToBedAndMergeL, function(x){getOutResultName(getOutResultReference(x[[2]]))})

if(grouping){
  groupVect = 1:length(samplesInfo$condition)
  names(groupVect) = unique(samplesInfo$condition)
  groupVect = groupVect[samplesInfo$condition]
  multiIntersectBed_perl_CLI = MultiIntersectBed_perl_CLI(inFilePath = inFP, inFileNames = inFN, outputFlag = "",withinGroupTH = withinGroupTH, groupVect = groupVect, outFileName = "multiIntersectClust",outFilePath = contigAssemblyDir)
} else{
  multiIntersectBed_perl_CLI = MultiIntersectBed_perl_CLI(inFilePath = inFP, inFileNames = inFN, outputFlag = "",withinGroupTH = withinGroupTH, outFileName = "multiIntersectClust",outFilePath = contigAssemblyDir)
}

tmpCommandLog = c(tmpCommandLog, paste0("\nmkdir ",contigAssemblyDir,"\n") )
dir.create(contigAssemblyDir)
multiIntersectBed_perl_CLI_cmdRes = generateCommandResult(multiIntersectBed_perl_CLI)
tmpCommandLog = c(tmpCommandLog, getCommandLog(multiIntersectBed_perl_CLI_cmdRes) )
multiIntersectBed_perl_CLI_exec = executeCommandResult(multiIntersectBed_perl_CLI_cmdRes, testing=FALSE)

#Execture sed command to eliminate column 4!

#if(bedTools225){
#  sedCmd = paste0("sed -i -r 's/(\\s+)?\\S+//4' ", file.path(multiIntersectBed_perl_CLI_cmdRes@CLIApplication@outFilePath, getOutResultName(getOutResultReference(multiIntersectBed_perl_CLI_cmdRes)))
#  )
#  tmpCommandLog = c(tmpCommandLog, sedCmd)
#  system(sedCmd, intern=TRUE)
#}
#############################################
#   Clusering Preparation
#     Contig Assembly - between samples - MultiIntersectBed
#    
#     Intersect the multiIntersected file with each sample to obtain all reads overlapping contigs (with read names)
#     intermediate step for merging and uniqueness calculation below
#############################################
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

intersectBed_CLI_cmdExecL = lapply( intersectBed_CLI_cmdResL, function(x){
  executeCommandResult(x,  testing=FALSE)
})
#############################################
#   Clusering Preparation
#     Contig Assembly - between samples - MultiIntersectBed
#
#     Merging these to count the reads and also calculate the genomic uniqueness in each contig
#     These values (read counts and uniqueness) are used for clustering, since the multi mapping reads are extracted, and ordered by read number
#############################################
mergeBedFile_CLI_cmdResL = lapply(intersectBed_CLI_cmdResL, function(x){
  mergeBedFile_CLI = MergeBedFile_CLI(inFilePath = getOutFilePath(getCLIApplication(x)), 
                                      inFileNames = getOutResultName(getOutResultReference(x)), 
                                      cliParams = paste0("-s -d ",readOverlap_contig," -c 4,5,6,10 -o count,mean,distinct,distinct",bedTools225Sed), 
                                      outputFlag = "_merged", 
                                      outFilePath = getOutFilePath(getCLIApplication(x)))  
  mergeBedFile_CLI_cmdRes = generateCommandResult(mergeBedFile_CLI)
  return(mergeBedFile_CLI_cmdRes)
})

tmpCommandLog = c(tmpCommandLog, sapply(mergeBedFile_CLI_cmdResL, getCommandLog) )
mergeBedFile_CLI_Exec2L = lapply( mergeBedFile_CLI_cmdResL, function(x){
  executeCommandResult(mergeBedFile_CLI_cmdResL[[i]],  testing=FALSE)
})

#######################################
#   Clusering Preparation
#
#     Finally intersect the contig with the bed files of the individual experiments
#     This creates bed files containing all the reads only for the contig regions, to minimize computation
#     One file from each group is chosen for clustering (since the clustering method makes use of read names in a contig!)
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

intersectBedClust_CLI_cmdExecL = lapply( intersectBed_CLI_clust_cmdResL, function(x){
  executeCommandResult(x,  testing=FALSE)
})

###################################################
# Concatenate files to obtain the file containing reads for clustering
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
readsForClustering_cmdExecL = lapply( readsForClustering_cmdL, function(x){
  system(x, intern=TRUE)
})


######################################
#       Write commands to file
######################################
cat(tmpCommandLog,  file = commandLog, append = FALSE) 
save.image(file.path(rootDir,"CommandsBeforeClustering.rda"))
######################################
#       Execute commands! (already executed -> if test = TRUE in the above examples then this could be uncommented)
######################################
# cmdExecTime = proc.time()
# setwd(rootDir)
# system(paste0("bash ", commandLog, " > ", executionLog ))
# cmdExecTime = proc.time() - cmdExecTime
# cmdExecTime


########################################################################################################################################################
#                                                                       Contig clustering
########################################################################################################################################################

################################
#     Now combined these files -> in R currently
################################
require(rtracklayer)
require(RNASeqUtility)

#########################################
#     This method combines the multiIntersect table with the individual sample files that containg the uniqueness and read counts!
#     The resulting file is is split into contigs compososed solely of unique reads and multi mapping reads
#########################################
unclDF = generateCountUnclusteredTable(contigFile=multiIntersectBed_perl_CLI_cmdRes, sampleBedFiles=mergeBedFile_CLI_cmdResL)
colnames(unclDF) = gsub(".*\\.","",colnames(unclDF) )

#Now This table is reduced to mean ReadCount and Mean Uniqueness and can then be used for clustering!
toClusterDF = unclDF[,1:5]
toClusterDF$readCount = rowMeans(  unclDF[,which(c(rep(FALSE,5),c(6:dim(unclDF)[2]%%2 == 0)))] ,na.rm = TRUE)
toClusterDF$uniqueness = rowMeans(  unclDF[,which(c(rep(FALSE,5),c(6:dim(unclDF)[2]%%2 == 1)))], na.rm = TRUE )

##################################################
#     Short sanity checking ...
finiteReadCount = which( is.finite(toClusterDF$readCount) )
if( length(finiteReadCount) != dim(toClusterDF)[1] ){warning("Some contigs do not contain ANY overlapping read -> check!"); toClusterDF = toClusterDF[finiteReadCount,]}
finiteUniq = which( is.finite(toClusterDF$uniqueness) )
if( length(finiteUniq) != dim(toClusterDF)[1] ){warning("Some contigs do not contain ANY overlapping read -> check!"); toClusterDF = toClusterDF[finiteUniq,]}
###################################################
contigForCountingGR = with(toClusterDF, GRanges(seqnames=chr, IRanges(start=start, end=end), strand = strand, readCount = readCount, 
                                                uniqueness =uniqueness))
#Delete all sequences that are shorter than 18 nt
contigForCountingGR = contigForCountingGR[width(contigForCountingGR) >= 18]

#Combine
contigForCountingGR_unique = contigForCountingGR[contigForCountingGR$uniqueness <= 1]
contigForCountingGR_unclustered = contigForCountingGR[contigForCountingGR$uniqueness > 1]

allReads = import(file.path( contigClusterDir ,readsForClustering), 
                  format="BED", asRangedData=FALSE)#alignedbwt_s.bed
start(allReads) = start(allReads) - 1 #Strange behaviour in import of rtracklayer!


##################################################
#         Contig Clustering (see method description of clusterMultiMappingReads_stringent)
##################################################
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
#Writing the clustered contigs as bed file (when using export function of rtracklayer -> coordinates get -1, since bed files are 0 coordinate based -> in this case would be wrong)
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
#         Counting the reads by employing MultiBamCoverage from bed tools. 
#         Current parameters are: -s for strandness, -f 0.05 (minimum 5% of the read length must overlap a contig -> 1 nt in 20 nt read); -D include duplicated reads (not necessary in this case due to clustering)
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


# rownames(count_contigsDF ) = count_contigsDF$contigID
# rownames(ensReadsCountDF_clustered ) = ensReadsCountDF_clustered$UID
# counts = ensReadsCountDF_clustered[,3:(2+length(samplesInfo$condition))]
# counts = rbind(counts, count_contigsDF[,7:dim(count_contigsDF)[2]])
# counts_tmp = as.data.frame(apply(counts, 2, as.numeric))
# rownames(counts_tmp) = rownames(counts)


setwd(rootDir)
sapply( commandLogCounting, system)

rm(allReads)#to save memory!
save.image(file = "MappingAssemblyClustering.rda")

##################################################
#   Creates a NAME_multibamcov file containing the read counts for the clustered contigs!
#   The ensembl read count table will be created in the Annotation script, since a further clustering by annotation is done
##################################################
