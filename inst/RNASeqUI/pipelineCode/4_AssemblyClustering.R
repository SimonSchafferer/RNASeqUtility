print("Assembly and Clustering")

options(stringsAsFactors = FALSE)

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
