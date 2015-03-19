
args <- commandArgs(TRUE)
if( length(args) == 0 ){
  load( file.path(getwd(),"Configuration.rda") )
} else{
  load(file.path(args[1],"Configuration.rda"))
}

load(file.path(rootDir, "Annotation.rda"))
load(file.path(rootDir, "reports","workspace.rda"))

require(rtracklayer)
require(RNASeqUtility)
require(ggplot2)

coverageAnalysisDir = file.path(rootDir, "reports","readCoverageAnalysis")
dir.create(coverageAnalysisDir)
setwd(coverageAnalysisDir)

##################
#   Get Differentially Expressed candidates from the report and annotate them
##################
resOrdered_tmp = resOrdered
resOrdered_tmp$UID = rownames(resOrdered_tmp)
resOrdered_tmp = as.data.frame(resOrdered_tmp)
resAnnot = merge(resOrdered_tmp, minimalAnnotationDF, by="UID", sort=FALSE)
rm(resOrdered_tmp)
signCand = resAnnot[which(resAnnot$padj < 0.1),]

######################################################
#   Plot candidates ncRNA mappings (not contigs)
######################################################
#First subset to non-contigs
signCand_other = signCand[-grep("contig",signCand$UID),]
#Generate coverage dataframe, containing read counts at each position
roi_ncOtherDFL = generateCoverageDF_ncRNAmappings( signCand_other, genomeCovNcRNAMapping_plusL, genomeCovNcRNAMapping_minusL, samplesInfo$sampleName)
roi_ncOtherDF = do.call(rbind, roi_ncOtherDFL)

#Plot and save the plot to disk
p2 = ggplot(data = roi_ncOtherDF, aes(x=Idx,y = ReadCount, color=Sample)) + geom_density(stat="identity") + 
  facet_wrap( ~ID, ncol=floor( sqrt( length( unique(roi_ncOtherDF$ID) ) )), scales = "free")
ggsave(p2, filename = "otherncRNAsCoverage.png")

#####################################################################
#     Coverage plots with mean normalized read counts (ncRNAs)
#####################################################################
#Generate coverage dataframe, containing read counts at each position
roi_ncOtherDF_cond = RNASeqUtility::splitReadCountCoveragePerCondition(coverageL = roi_ncOtherDFL, deseqds = dds, samplesInfo = samplesInfo)

roi_ncOtherDF_cond_plot = roi_ncOtherDF_cond[,c("Name","condition","ReadMeansNorm","Idx")]
roi_ncOtherDF_cond_plot = unique(roi_ncOtherDF_cond_plot)
p2 = ggplot(data = roi_ncOtherDF_cond_plot, aes(x=Idx,y = ReadMeansNorm, color=condition)) + geom_density(stat="identity") + 
  facet_wrap( ~Name, ncol=floor( sqrt( length( unique(roi_ncOtherDF_cond_plot$Name) ) )), scales = "free")
ggsave(p2, filename = "otherncRNAsCoverageNorm.png")

#################################################################
#     Contig Coverage Plots
#################################################################
#Get all significant candidates that are contigs and not already annotated ncRNAs
signCand_contig = signCand[grep("contig",signCand$UID),]

contigCov = generateCoverageDF_contigs(contigsGR=count_contigsGR, 
                                       contigResTable=signCand_contig, 
                                       bedgraphMapping_plusL=genomeCovMapping_plusL, 
                                       bedgraphMapping_minusL=genomeCovMapping_minusL, 
                                       sampleNames=samplesInfo$sampleName)
contigCovDF = do.call(rbind, contigCov)

#Plot and save the plot to disk
pContig = ggplot(data = contigCovDF, aes(x=Idx,y = ReadCount, color=Sample )) + geom_density(stat="identity") + 
  facet_wrap( ~ID, ncol=floor( sqrt( length( unique(contigCovDF$ID) ) )), scales = "free")
ggsave(pContig, filename = "contigProcessing.png")


#####################################################################
#     Coverage plots with mean normalized read counts (contigs)
#####################################################################
roi_contigCov_cond = RNASeqUtility::splitReadCountCoveragePerCondition(coverageL = contigCov, deseqds = dds, samplesInfo = samplesInfo)
roi_contigCov_cond_plot = roi_contigCov_cond[,c("ID","condition","ReadMeansNorm","Idx")]
roi_contigCov_cond_plot = unique(roi_contigCov_cond_plot)
p2 = ggplot(data = roi_contigCov_cond_plot, aes(x=Idx,y = ReadMeansNorm, color=condition)) + geom_density(stat="identity") + 
  facet_wrap( ~ID, ncol=floor( sqrt( length( unique(roi_contigCov_cond_plot$ID) ) )), scales = "free")
ggsave(p2, filename = "contigProcessingMeanNorm.png")



###################################
#     Contig Coverage Plots with annotation
###################################
# It is also possible to add an annotation to the coverage plot, this means, that the annotation will be taken as reference and therefore, 
# the coverage of an ncRNA can be investigated (this is for example useful for tRNAs to investigate processing)
# In order to get your annotation of interest download for example at UCSC table browser
#Currently this code is commented, since it needs to be modified before running!

#<CODE BEGINS HERE>
# biotypeOfInterest = "tRNA"
# allTRNAs = import( "PATHTOYOURANNOTATION_BEDFILE", asRangedData=FALSE)
# signTRNAs = signCand_contig[which(signCand_contig$biotype == biotypeOfInterest),] #subset to 
# tRNACov = generateCoverageDF_contigs(contigsGR=count_contigsGR, 
#                                      contigResTable=signTRNAs, 
#                                      bedgraphMapping_plusL=genomeCovMapping_plusL, 
#                                      bedgraphMapping_minusL=genomeCovMapping_minusL, 
#                                      sampleNames=samplesInfo$sampleName,
#                                      annotationGR=allTRNAs)
# tRNACovDF = do.call(rbind, tRNACov)
# pTRNA = ggplot(data = tRNACovDF, aes(x=Idx,y = ReadCount, color=Sample )) + geom_density(stat="identity") + 
#   facet_wrap( ~ID, ncol=floor( sqrt( length( unique(tRNACovDF$ID) ) )), scales = "free")
# ggsave(pTRNA, filename = "tRNAProcessing.png")
#<CODE ENDS HERE>

