
# running the main scripts with arguments: 
# R CMD BATCH --no-save --no-restore '--args /tmp/test' test.R, or Rscript test.R /tmp/test 
# the second command gives only output of R, the first one outputs all code too

options(stringsAsFactors = FALSE)

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

