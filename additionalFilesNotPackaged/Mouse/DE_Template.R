rootDir = file.path( "/home","simon","PHDStudies","RNA-Seq","IonProton","CavKO_Striessnig","MouseBrain_wo_Ref_extended")#schaffrr
# New Version -> one may load the workspace of the mapping assembly and clustering approach MappingAssemblyClustering.rda
load(file.path( rootDir,"MappingAssemblyClustering.rda"))#
#extracting a sample fastq subset! sampleFastq.fastq (100.000 entries)
options(stringsAsFactors = FALSE)
setwd(rootDir)

library("edgeR")
rownames(count_contigsDF ) = count_contigsDF$contigID
rownames(ensReadsCountDF_clustered ) = ensReadsCountDF_clustered$UID
counts = ensReadsCountDF_clustered[,3:(2+length(samplesInfo$condition))]
counts = rbind(counts, count_contigsDF[,7:dim(count_contigsDF)[2]])
counts_tmp = as.data.frame(apply(counts, 2, as.numeric))
rownames(counts_tmp) = rownames(counts)
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
d = DGEList(counts=counts, group=factor(samplesInfo$condition))
d = calcNormFactors(d)
# plotMDS(d, labels=c("MSA","WT")[sampleInfo$condition], col=c("darkgreen","blue")[sampleInfo$condition])
d = estimateCommonDisp(d)
# d = estimateTagwiseDisp(d)
# de = exactTest(d, pair = c("MSA", "WT"))
# tt = topTags(de, n=nrow(d))

#############################
#  For More Complex Designs...
#############################

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
plotMeanVar(d, show.tagwise.vars=TRUE, NBline = TRUE)
#plotBCV illustrates the relationship of biological coefficient of variation versus mean log CPM
plotBCV(d)
plotMDS(d)

nc = cpm(d, normalized.lib.sizes=TRUE)
rn = rownames(tt$table)
head(nc[rn, order(samplesInfo$condition)], 5)

deg = rn[tt$table$FDR < 0.1]
plotSmear(d, de.tags=deg)

# write.table(tt$table, file="toptags_edgeR.csv")
save.image(file.path(rootDir, "DifferentialExpression.rda"))


