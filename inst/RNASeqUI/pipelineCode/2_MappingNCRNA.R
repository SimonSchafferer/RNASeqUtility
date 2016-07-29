print("Mapping to ncRNA Genome")

options(stringsAsFactors = FALSE)

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
