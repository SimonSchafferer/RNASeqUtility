
##############################################
# DOWNLOAD THE ncRNA fasta file from ENSEMBL
# SUBSET THIS FILE BY THE MAXIMUM READ LENGTH -> in our case 200
# bring the file into the right format (no line breaks if there are!)
#/home/simon/dbsOfflineUse/MusMusculus/Mouse_ncRNA_fasta/reformat_fasta.pl
#or python script EliminateNewLinesFasta.py
#/home/simon/dbsOfflineUse/MusMusculus/Mouse_ncRNA_fasta/Mus_musculus.GRCm38.ncrna.mod.fa
# awk '!/^>/ { next } { getline seq } length(seq) <= 200 or 400 for size selection { print $0 "\n" seq }' FastaFileInput.fa > FastaFileOutput.fa
# create genome index file 
# STAR --runThreadN 5 --runMode genomeGenerate --genomeDir /home/simon/dbsOfflineUse/MusMusculus/Mouse_ncRNA_fasta/sizeRestricted --genomeFastaFiles /home/simon/dbsOfflineUse/MusMusculus/Mouse_ncRNA_fasta/sizeRestricted/Mus_musculus.GRCm38.ncrna_sizeRestricted.fa
# Map all reads to this artificial genome containing ncRNAs only!
# proceed with unmapped reads as follows
# This pipeline may be run on unmapped reads from primary mapping
# STAR --genomeDir /home/simon/dbsOfflineUse/MusMusculus/Mouse_ncRNA_fasta --readFilesIn /home/simon/PHDStudies/RNA-Seq/IonProton/CavKO_Striessnig/MouseBrain_wo_Ref/rawData/test/IonXpressRNA_001.R_2014_08_22_11_48_07_user_PRO-36-2014-08-22_MoBrain_Ion_RNA_-_20-200nt_Transcriptom.fastq --outFileNamePrefix test1 --runThreadN 8 --outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin 16 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --alignIntronMax 1 --outFilterMultimapNmax 100 --outSAMprimaryFlag AllBestScore --outReadsUnmapped Fastx
# 
# # --outSAMtype BAM SortedByCoordinate
# 
# samtools view -bS test1Aligned.out.sam > test1Aligned.out.bam
# bamToBed -i test1Aligned.out.bam > test1Aligned.out.bed
# sort -k1,1 -k2,2n test1Aligned.out.bed | mergeBed -s -c 4,5,6 -o count,mean,distinct -i stdin | awk '{if($4 > 5) print }'> test1Aligned.out.merged.bed
###############################################
# STAR --genomeDir /home/simon/dbsOfflineUse/MusMusculus/Mouse_ncRNA_fasta/sizeRestricted --readFilesIn /home/simon/PHDStudies/RNA-Seq/IonProton/CavKO_Striessnig/MouseBrain_wo_Ref/rawData/test/IonXpressRNA_001.R_2014_08_22_11_48_07_user_PRO-36-2014-08-22_MoBrain_Ion_RNA_-_20-200nt_Transcriptom_trimmed.fastq --outFileNamePrefix test1 --runThreadN 8 --outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin 16 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --alignIntronMax 1 --outFilterMultimapNmax 100 --outSAMprimaryFlag AllBestScore --outReadsUnmapped Fastx --outSAMtype BAM Unsorted

#samtools view -bS test1Aligned.out.sam > test1Aligned.out.bam
# bamToBed -i test1Aligned.out.bam > test1Aligned.out.bed
# sort -k1,1 -k2,2n test1Aligned.out.bed | mergeBed -s -c 4,5,6 -o count,mean,distinct -i stdin | awk '{if($4 > 5) print }'> test1Aligned.out.merged.bed

# wc -l test1Aligned.out.bed 14029759 reads
# wc -l test1Aligned.out.merged.bed 1368 ncRNAs
# some stats
#The names need to be modified to provide unique IDs:  _se or _as
# ENSMUST00000082825  31	65	6	0	+
# ENSMUST00000082825	96	169	14	183	-
# in the end get the coordinates and names from biomart: EnsemblGeneID, EnsemblTranscriptID, GeneStart, GeneEnd, ChromosomeName, Strand, MGIsymbol, MGItranscriptname
# MGI symbol and transcript name may be employed for annotation clustering!

# with all unmapped reads perform the following script...






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
rootDir = file.path( "/home","simon","PHDStudies","RNA-Seq","IonProton","CavKO_Striessnig","MouseBrain_wo_Ref_extended")#schaffrr
rawDataDir = file.path(rootDir,"rawData")
gtfFileName = "Mus_musculus.GRCm38.78_mod.gtf"
gtfFilePath = "/home/simon/dbsOfflineUse/GTF_repos/Mus_musculus.GRCm38.78_repo/"
transcriptomeIndex = "/home/simon/dbsOfflineUse/GTF_repos/Mus_musculus.GRCm38.78_repo/Mus_musculus.GRCm38.78"
genomeIndexFilePath = "/home/simon/dbsOfflineUse/MusMusculus/Mouse_mm10_fasta/rnaStar"

# genomeIndexFilePath_ncRNA = createSizeSelectedEnsemblncRNADB(genomeDBdir = <FILEPATH>, ncRNAsize = 400, nrThreads = 6, ensemblFTP = "mouse")
genomeIndexFilePath_ncRNA = "/home/simon/dbsOfflineUse/MusMusculus/Mouse_ncRNA_fasta/size400"#please uncomment above if a db dir has to be created
ncRNAmappingDir = file.path(rootDir, "ncRNAmapping")
dir.create(ncRNAmappingDir)

mappingDir = file.path(rootDir,"mapping")
dir.create(mappingDir)
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
grouping = TRUE #This flag defines if the groups should be considered when intersecting all contigs...
withinGroupTH = 0

samplesInfo = read.table(file=file.path(rootDir, "samplesInfo.csv"), sep="\t", header=TRUE) #should contain cloumn condition and column sample name

if(!( "sampleName" %in% colnames(samplesInfo) & "condition"%in% colnames(samplesInfo))) {stop("Please provide a sampleInfo file with column condition and column sampleName")}

############
#   Parameter definition
############
read_threshold = 5 # read threshold for contig assembly
readOverlap_contig = -3 #This parameters defines the number of nt a read has to overlap with the next in order to be clustered (from positive to negative numbers -> APART default: 0)
readCompositionIdentity = 0.95 #This is used for clustering: If 0.95: At least 95% of the reads have to be shared between a contig with the representative contig (highest read count) in order to be clustered, 
# so in this case if there are 5% reads that are different then the contig will be kept next to the representative contig. When this value is set to 1 then contigs will be clustered to the 
# representative contig if they share the same reads. If this is set to 0, then a contig will be deleted if it shares 1 or more reads with an representative contig!

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
cutAdaptCLI = Cutadapt_CLI(inFilePath=rawDataDir, cliParams =  c("-a ATCACCGACTGCCCATAGAGAGGCTGAGAC --minimum-length 18"), outputFlag = "_trimmed", 
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
                           cliParams = "--runThreadN 8 --outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin 16 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --alignIntronMax 1 --outFilterMultimapNmax 100 --outSAMprimaryFlag AllBestScore --outReadsUnmapped Fastx", 
                           outputFlag = paste0(outFN,"_"), #prefix -> must be distinguishable -> best to use the fastQFileNames ,or the sample Short Names!!
                           outFilePath = outFP, 
                           genomeIndexFilePath = genomeIndexFilePath_ncRNA,
                           filterSam = TRUE, 
                           outputFormat = "sam")
  #"--runThreadN 6 --outFilterMismatchNmax 1 --outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin 16 --outFilterScoreMinOverLread 0  --outFilterMatchNminOverLread 0 --alignIntronMax 1 --outFilterMultimapNmax 100"
  #"--runThreadN 8 --outFilterMismatchNoverReadLmax 0.023 --outFilterMatchNmin 18 --outFilterScoreMinOverLread 0  --outFilterMatchNminOverLread 0 --alignIntronMax 1 --outFilterMultimapNmax 100 --alignEndsType EndToEnd --outSAMprimaryFlag AllBestScore"
  
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
                           cliParams = "--runThreadN 8 --outFilterMismatchNoverReadLmax 0.023 --outFilterMatchNmin 18 --outFilterScoreMinOverLread 0  --outFilterMatchNminOverLread 0 --alignIntronMax 1 --outFilterMultimapNmax 100 --alignEndsType EndToEnd --outSAMprimaryFlag AllBestScore --outSAMtype BAM Unsorted", 
                           outputFlag = outFN, #prefix -> must be distinguishable -> best to use the fastQFileNames ,or the sample Short Names!!
                           outFilePath = outFP, 
                           genomeIndexFilePath = genomeIndexFilePath,
                           filterSam = FALSE, 
                           outputFormat = "bam")
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

