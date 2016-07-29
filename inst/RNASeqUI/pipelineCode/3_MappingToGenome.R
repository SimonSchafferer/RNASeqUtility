print("Mapping to Genome")

options(stringsAsFactors = FALSE)
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
