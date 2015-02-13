args <- commandArgs(TRUE)
if( length(args) == 0 ){
  load( file.path(getwd(),"Configuration.rda") )
} else{
  load(file.path(args[1],"Configuration.rda"))
}

#####################################################################################
#                Differential Expression configuration 
#####################################################################################

load(file.path( rootDir,"Annotation.rda"))#
#extracting a sample fastq subset! sampleFastq.fastq (100.000 entries)
options(stringsAsFactors = FALSE)
setwd(rootDir)

library("edgeR")
library(CLIHelperPackage)
library(RNASeqUtility)
require(rtracklayer)

rownames(count_contigsDF ) = count_contigsDF$contigID
rownames(ensReadsCountDF_clustered ) = ensReadsCountDF_clustered$UID
counts = ensReadsCountDF_clustered[,3:(2+length(samplesInfo$condition))]
counts = rbind(counts, count_contigsDF[,7:dim(count_contigsDF)[2]])
counts_tmp = as.data.frame(apply(counts, 2, as.numeric))
rownames(counts_tmp) = rownames(counts)

if( dim(otherReadsCountDF)[1] > 0 ){
  otherReadsCountDF_tmp = otherReadsCountDF[, which( !colnames(otherReadsCountDF) %in% c("mapStart","mapEnd")  ) ]
  otherReadsCountDF_tmp = otherReadsCountDF_tmp[,-1]
  otherReadsCountDF_tmp = as.data.frame(apply(otherReadsCountDF_tmp, 2, as.numeric))
  rownames(otherReadsCountDF_tmp) = otherReadsCountDF[,1]
  counts_tmp = rbind( counts_tmp, otherReadsCountDF_tmp  )
}

counts = counts_tmp
#GET ALL COUNT Data othe
cpms = cpm(counts)
keep = rowSums(cpms > 1) >=2
counts = counts[keep,]

colnames(counts) = samplesInfo$sampleName
# head(counts[,order(samplesInfo$condition)],5)

#############################
#  Simple Designs...
#############################
# d = DGEList(counts=counts, group=factor(samplesInfo$condition))
# d = calcNormFactors(d)
# plotMDS(d, labels=c("MSA","WT")[sampleInfo$condition], col=c("darkgreen","blue")[sampleInfo$condition])
# d = estimateCommonDisp(d)
# d = estimateTagwiseDisp(d)
# de = exactTest(d, pair = c("MSA", "WT"))
# tt = topTags(de, n=nrow(d))

#############################
#  For More Complex Designs...
#############################
d = DGEList(counts=counts, group=factor(samplesInfo$condition))
d = calcNormFactors(d)
d = estimateCommonDisp(d)

design = model.matrix( ~condition, samplesInfo  )
d = estimateGLMTrendedDisp(d, design)
d = estimateGLMTagwiseDisp(d, design)
f = glmFit(d, design)
de = glmLRT(f, coef=2)

tt = topTags(de, n=nrow(d))
# tt1 = topTags(de1, n=nrow(d))

head(tt$table)
# head(tt1$table)

#mean-variance relationship
#Each dot represents the estimated mean and variance for each gene, with binned variances as well as the trended common dispersion overlaid
png(filename = file.path(diffExpDir,"MeanVar.png"))
plotMeanVar(d, show.tagwise.vars=TRUE, NBline = TRUE)
dev.off()
#plotBCV illustrates the relationship of biological coefficient of variation versus mean log CPM
png(filename = file.path(diffExpDir,"BCV.png"))
plotBCV(d)
dev.off()

png(filename=file.path(diffExpDir,"MDS.png"))
plotMDS(d)
dev.off()

nc = cpm(d, normalized.lib.sizes=TRUE)
rn = rownames(tt$table)
head(nc[rn, order(samplesInfo$condition)], 5)

deg = rn[tt$table$FDR < 0.1]

png(filename=file.path(diffExpDir,"SmearMA.png"))
plotSmear(d, de.tags=deg)
dev.off()
write.table(tt$table, file=file.path(diffExpDir,"toptags_edgeR.csv"), sep="\t")

####################################################################################
#                       DESeq
####################################################################################

library("DESeq")

# libDesign = data.frame("condition"=samplesInfo$condition, "libType"=rep("single-end",length(samplesInfo$condition)))
# cds = newCountDataSet( counts_tmp, samplesInfo$condition)
# counts.filtered = counts[-grep("contig",rownames(counts)),]
# counts.filtered = counts.filtered[apply(counts.filtered,1,min) != 0,]

cds = newCountDataSet( counts, samplesInfo$condition)

cds = estimateSizeFactors(cds)
sizeFactors(cds)

cdsB = estimateDispersions(cds, method="blind")
vsd = varianceStabilizingTransformation(cdsB)
p = plotPCA(vsd,ntop = 100,intgroup =c("condition")  )

png(filename = file.path(diffExpDir,"DESeq_PCA.png"))
p
dev.off()

cds = estimateDispersions(cds)

png(filename = file.path(diffExpDir,"DESeq_DispersionEstimation.png"))
plotDispEsts(cds)
dev.off()

conditionsNbinom = unique(samplesInfo$condition)

res = nbinomTest(cds,condA = conditionsNbinom[1], condB=conditionsNbinom[2])

png(filename = file.path(diffExpDir,"DESeq_MA.png"))
plotMA(res)
dev.off()
# resSig = res[which(res$padj < 0.1),]
# head(resSig[order(resSig$log2FoldChange, decreasing=TRUE, )])

deseqTable = res[ order(res$pval, decreasing=FALSE),]

write.table(deseqTable, file=file.path(diffExpDir,"toptags_DESeq.csv"))

save.image(file.path(rootDir, "DifferentialExpression.rda"))


