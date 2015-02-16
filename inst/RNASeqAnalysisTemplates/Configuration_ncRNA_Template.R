##########################################################################################################################################################################
#                                                                         Configuration file
#                                                           Mapping Assembly and Clustering Script options
##########################################################################################################################################################################
#Execute Scripts with R CMD BATCH --no-save --no-restore '--args /tmp/test' test.R, or Rscript test.R /tmp/test 

######################################################
#   directories and file paths that NEED to be changed!
######################################################

#Base directory of the analysis
rootDir = file.path( "/tmp","test","","","","","")#schaffrr
#Directory containging fastq files
rawDataDir = file.path(rootDir,"rawData")

#Path to rna star mapper Index file for whole genome (download genome fasta file from UCSC -> twoBit format -> convert with twoBitToFa )
genomeIndexFilePath = "/home/simon/dbsOfflineUse/HomoSapiens/hg19/rnaStarIndex"

#Path to rna star mapper Index file for ncRNAs -> Please use RNASeqUtility::createSizeSelectedEnsemblncRNADB method! *see Prepare ncRNAIndex below
genomeIndexFilePath_ncRNA = "/home/simon/dbsOfflineUse/HomoSapiens/hg19/ncRNA_ENSEMBL/restr400nt_extended"#please uncomment above if a db dir has to be created

##############
#     For the annotation Script
##############
#Path to infernal database (compiled)
infernalDB = file.path("/home","simon","dbsOfflineUse","HomoSapiens","hg19","infernal","Rfam.cm.1_1")

#Path to repeat masker file (Please see below for further information)
repeatMaskerDir = "/home/simon/dbsOfflineUse/HomoSapiens/hg19/repeatMskr"
#Name of the repeat masker file (either .rda (see sncRNAannotation package) or bed file!)
repeatMaskerFN = "rmsk_hg19.bed" #or GRanges object

#Assembly of Annotation:
organismForAnnotation = "hg19" #
#Supported types: hg19, hg38, mm9, mm10


######################################################
#              Default parameters
######################################################
grouping = TRUE #This flag defines if the groups should be considered when intersecting all contigs... *see Contig Assembly below
withinGroupTH = 0 #This is the threshold for contig assembly within the groups, it defines how many samples of a group do not need to have reads for this contig *see Contig Assembly below
read_threshold = 1 # read threshold for contig assembly 
readOverlap_contig = -1 #This parameters defines the number of nt a read has to overlap with the next in order to be clustered (from positive to negative numbers -> APART default: 0)
readCompositionIdentity = 0.95 #This is used for clustering: If 0.95: At least 95% of the reads have to be shared between a contig with the representative contig (highest read count) in order to be clustered, 
# so in this case if there are 5% reads that are different then the contig will be kept next to the representative contig. When this value is set to 1 then contigs will be clustered to the 
# representative contig if they share the same reads. If this is set to 0, then a contig will be deleted if it shares 1 or more reads with an representative contig!


####################################
#     The paths to the command line programs should be set HERE!
####################################
#This code greps all export PATH statement from the bashrc file
pathvars = readLines(file.path(path.expand("~"),".bashrc"))
pathvars = pathvars[grep("export PATH\\=\\$PATH:",  pathvars )]
pathvars = sub( "export PATH\\=\\$PATH:", "",pathvars)
Sys.setenv(PATH=paste(Sys.getenv("PATH"),paste(pathvars,collapse=":"),sep=":")) 


######################################
#       Options that may be changed 
######################################

########################
# Cutadapt Params ncRNA
########################
cutadaptOptions = "-a ATCACCGACTGCCCATAGAGAGGCTGAGAC --minimum-length 18"

########################
# RNAStar Params ncRNA
########################
#BEWARE: if the option --outSAMtype BAM Unsorted is set to bam then the variable rnaStarncRNA_outputFormat has to be set to bam and the rnaStarncRNA_filterSam to FALSE!
#DO NOT SET --outSAMtype BAM  to Sorted!!!
rnaStarncRNA_params = "--runThreadN 8 --outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin 16 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --alignIntronMax 1 --outFilterMultimapNmax 100 --outSAMprimaryFlag AllBestScore --outReadsUnmapped Fastx"
rnaStarncRNA_filterSam = TRUE 
rnaStarncRNA_outputFormat = "sam"
#other param suggestions
#"--runThreadN 6 --outFilterMismatchNmax 1 --outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin 16 --outFilterScoreMinOverLread 0  --outFilterMatchNminOverLread 0 --alignIntronMax 1 --outFilterMultimapNmax 100"
#"--runThreadN 8 --outFilterMismatchNoverReadLmax 0.023 --outFilterMatchNmin 18 --outFilterScoreMinOverLread 0  --outFilterMatchNminOverLread 0 --alignIntronMax 1 --outFilterMultimapNmax 100 --alignEndsType EndToEnd --outSAMprimaryFlag AllBestScore"

########################
# RNAStar Params genome
########################
rnaStarGenome_params = "--runThreadN 8 --outFilterMismatchNoverReadLmax 0.023 --outFilterMatchNmin 18 --outFilterScoreMinOverLread 0  --outFilterMatchNminOverLread 0 --alignIntronMax 1 --outFilterMultimapNmax 100 --alignEndsType EndToEnd --outSAMprimaryFlag AllBestScore --outSAMtype BAM Unsorted"
rnaStarGenome_filterSam = FALSE
rnaStarGenome_outputFormat = "bam"


######################################################
#              Default directories
######################################################

#Name and path of the ncRNA mapping directory
ncRNAmappingDir = file.path(rootDir, "ncRNAmapping")
dir.create(ncRNAmappingDir)

#Name and path of the mapping directory
mappingDir = file.path(rootDir,"mapping")
dir.create(mappingDir)

#Name and path of the Contig Assembly directory
contigAssemblyDir = file.path(rootDir,"contigAssembly")
dir.create(contigAssemblyDir)

#Name and path of the Contig Clustering directory
contigClusterDir = file.path(rootDir,"contigClustering")
dir.create(contigClusterDir)

#Name and path of the Read Counts directory
readCountsDir = file.path(rootDir,"readCounts")
dir.create(readCountsDir)

#Name and path of the file containing the commands that are executed in bash
commandLog = file.path(rootDir, "Commands") 
if(!file.exists(commandLog)){
  file.create(commandLog)
}

#Name and path of the Log file for Commands that are executed in bash
executionLog = file.path(rootDir, "executionLog") 
if(!file.exists(executionLog)){
  file.create(executionLog)
}

#Annotation directory name
annotationDir = file.path(rootDir,"annotation")
dir.create(annotationDir)

#Differential Expression result directory name
diffExpDir = file.path(rootDir,"diffExpAnalysis")
dir.create(diffExpDir)

#Name and Path of the tab separated file containing the sample information
samplesInfo = read.table(file=file.path(rootDir, "samplesInfo.csv"), sep="\t", header=TRUE) #should contain cloumn condition and column sample name
if(!( "sampleName" %in% colnames(samplesInfo) & "condition"%in% colnames(samplesInfo))) {stop("Please provide a sampleInfo file with column condition and column sampleName")}

#Saving the configuration file: DO NOT CHANGE ITS NAME!
save.image(file.path(rootDir,"Configuration.rda") )

####################################################################################################################################
#                                         Additional Infos
#####################################################################################################################################

#############################################
#           Contig Assembly - among samples - MultiIntersectBed
#
#   PERL script is doing the multi intersection for plus and minus strand bed files -> and then merging them. 
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



#############################################
#     Prepare ncRNA index
#
#   Call method
#   RNASeqUtility::createSizeSelectedEnsemblncRNADB
#
#   if one wants to add other annotations next to ensembl please call this method with generateIdx = FALSE
#   Then prepare your files (for example fasta from mirbase) use: PreparencRNAfasta.py 
#   Script can be found by: file.path( system.file( package="RNASeqUtility"), "data", "PreparencRNAfasta.py")
#   This would be an example to change a fasta containing Us to Ts and replacing whiteSpaces in Ids with _ and replacing newLines in the sequences
#   PreparencRNAfasta.py -i testfasta.fa -o testfastaMod.fa --rmSpaceInID T --replacementChar '_' --RNAtoDNA T
#
#   Maybe adding the ncRNA class to the header with sed: 
#   sed -i 's/>/>snoRNA_hacabox_/g' HacaBox.fa
#
#   adding the file to the ENSEMBL ncRNA fasta file
#   cat <PREPARED NCRNA CLASS FASTA> >> <MODIFIED ENSEMBL NCRNA FASTA>
#
#   Now building the genome index with STAR 
#   STAR --runThreadN <NR> --runMode genomeGenerate --genomeDir <DIR>  --genomeFastaFiles <FASTA>
#
#   The resulting path of the index directory is then used for the variable: genomeIndexFilePath_ncRNA
#############################################

#############################################
#     Prepare repeat masker data
#   Download repeats masker bed file from UCSC (Table Browser) for the reuqired organism
# 
#   Optional, if more information is wanted, download from UCSC table browser and select fields from 'Select fields from primary and related table'
#   In this case the table has to be read into R and also changed to GRanges object accordingly
#   Read table
#   rpmsk = read.table( "<PATHTOTABLE>" )
#   convert to GRanges object
#   rmpsk_mm10 = with(rpmsk, GRanges( seq=<COLUMN>, IRanges(<COLUMN>,<COLUMN>), strand=<COLUMN>, ... ))
#   save object
#   save(rpmsk_mm10, file="<PATHTOFILE>")
#############################################


#############################################
#         Prepare annotation for sncRNAannotation
# call sncRNAannotation::createRObject_gtf function, or execute the steps in this function
#############################################

#############################################
#         Prepare infernal covariance model
#
#cmbuild -> use RFam pre-build db
#cmcalibrate [options] cmfile -> use RFam pre-build db
#
#  Start here when using the RFam pre-build db
#Compress the CM database first to accelerate cmscan!
# cmpress Rfam.cm.1_1
# this will generate various files (Rfam.cm.1_1.i1f,...), that have to be put in a directory and the path and infernal base name (in this case Rfam.cm.1_1) to these files has to be provided
# for the variable infernalDB
#############################################