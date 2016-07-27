print("Running Cutadapt")

options(stringsAsFactors = FALSE)
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
cutAdaptCLI = Cutadapt_CLI(inFilePath=rawDataDir, cliParams =  cutadaptOptions, outputFlag = "_trimmed", 
                           outFilePath = file.path(rootDir,"rawDataTrimmed") )
cutAdatptCLI_cmdRes = generateCommandResult(object = cutAdaptCLI )
#logging
tmpCommandLog = getCommandLog(cutAdatptCLI_cmdRes)

#Execute Cutadapt
cutAdatptCLI_execRes = executeCommandResult(cutAdatptCLI_cmdRes, testing=FALSE)


