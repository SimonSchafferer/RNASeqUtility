
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)

shinyServer(function(input, output, session) {
  sourceDir = file.path(getwd(),"tabPanels")
  pipelineDir = file.path(getwd(),"pipelineCode")
  choiceNone = "none"
  choiceAll = "all"
  customMessages = list(warningOutlier = "Cannot remove anymore outliers!")
  
  cutadaptSourceFN = "1_Cutadapt.R"
  ncRNAMappingSrcFN = "2_MappingNCRNA.R"
  genomeMappingSrcFN = "3_MappingToGenome.R"
  clusteringSrcFN = "4_AssemblyClustering.R"
  annotationSrcFN = "5_Annotation.R"
  diffExpSrcFN = "6_DiffExp.R"
  
  ######################
  #   Samples Info Upload
  ######################
  samplesInfo = reactive({
    rawDat = input$upl_file_csv
    req(rawDat)
    dat = read.csv(rawDat$datapath, header = TRUE,
                   sep = input$sep, quote = input$quote)
    return(dat)
  })
  
  output$upl_samplesInfo <- renderTable({
    req(samplesInfo())
    return(samplesInfo())
  })
  ######################
  #   Run options text respones
  ######################
  output$runScriptsOut <- renderText({ if( input$chooseRunScripts == 1){
    return("Running only Mapping");
  } else if(input$chooseRunScripts == 2){
    return("Running Mapping and annotation");
  } else{
    return("Running Mapping, annotation and \ndifferential expression");
  } })
  
  allOptions = function(){
    
    genomeIndexFilePath = input$config_genomeIndexFilePath
    
    return( list("genomeIndexFilePath"=genomeIndexFilePath) )
    
  };
  
  txtDecl = function(x,y){
    return(paste0(x,"= \"",y,"\""))
  }
  valDecl = function(x,y){
    return(paste0(x,"= ",y))
  }
  dirNameDec = function(x,y){
    return(paste0(x,"= file.path(rootDir,\"",y,"\")\ndir.create(",x,")"))
  }
  fileNameDec = function(x,y){
    return(paste0(x,"= file.path(rootDir,\"",y,"\")\nfile.create(",x,")"))
  }
  
  writeConfigFile = function(){
    cat(paste(
    txtDecl("rootDir",input$config_rootDir),
    txtDecl("genomeIndexFilePath",input$config_genomeIndexFilePath),
    txtDecl("genomeIndexFilePath_ncRNA",input$config_genomeIndexFilePath_ncRNA),
    txtDecl("infernalDB",input$config_infernalDB),
    txtDecl("repeatMaskerDir",input$config_repeatMaskerDir),
    txtDecl("repeatMaskerFN",input$config_repeatMaskerFN),
    txtDecl("organismForAnnotation",input$config_organismForAnnotation),
    valDecl("grouping",input$config_grouping),
    valDecl("withinGroupTH",input$config_withinGroupTH),
    valDecl("read_threshold",input$config_read_threshold),
    valDecl("readOverlap_contig",input$config_readOverlap_contig),
    valDecl("readCompositionIdentity",input$config_readCompositionIdentity),
    paste0("samplesInfo = read.csv(\"",file.path(input$config_rootDir, "SamplesInfo.txt"),"\", header=TRUE, sep=",")"),
    paste0("diffExpFormula = with(samplesInfo,",input$config_diffExpFormula,")"),
    txtDecl("perlPath",input$config_perlPath),
    txtDecl("cutadaptOptions",input$config_cutadaptOptions),
    txtDecl("rnaStarncRNA_params",input$config_rnaStarncRNA_params),
    valDecl("rnaStarncRNA_filterSam",TRUE),
    txtDecl("rnaStarncRNA_outputFormat","sam"),
    txtDecl("rnaStarGenome_params",input$config_rnaStarGenome_params),
    valDecl("rnaStarGenome_filterSam",FALSE),
    txtDecl("rnaStarGenome_outputFormat","bam"),
    dirNameDec("ncRNAmappingDir",input$config_ncRNAmappingDir),
    dirNameDec("mappingDir",input$config_mappingDir),
    dirNameDec("contigAssemblyDir",input$config_contigAssemblyDir),
    dirNameDec("contigClusterDir",input$config_contigClusterDir),
    fileNameDec("commandLog",input$config_commandLog),
    fileNameDec("executionLog",input$config_executionLog),
    dirNameDec("annotationDir",input$config_annotationDir),
    dirNameDec("readCountsDir",input$config_readCountsDir),
    dirNameDec("diffExpReportingDir",input$config_diffExpReportingDir),
    "save.image(file.path(rootDir,\"Configuration.rda\") )",sep="\n"), file=file.path(input$config_rootDir,"Configuration_ncRNA.R"))
  }
  
  
  callValidate = reactive({
    #Installed Programs
#     listOfPrograms = c("samtools --version", "bedtools --version", "cutadapt --version", "STAR --version", "cmsearch -h", "perl --version", "python --version", "awk --version", "cat --version", "sort --version", "blastn -version" )
#     programsNotFound = c("")
#     testPrograms = sapply(listOfPrograms, system)
    #Check Sample Info.txt
#     requiredCols = c("condition","sampleName","pathToFile","filename" )
#     
#     validate(
#       need(dir.exists(input$config_genomeIndexFilePath),"Please specify a valid genome index path directory.")
#     )
#     validate(
#       need(dir.exists(input$config_genomeIndexFilePath_ncRNA),"Please specify a valid ncRNA genome index path directory.")
#     )
#     validate(
#       need(dir.exists(input$config_repeatMaskerDir),"Please specify a valid repeat masker index path directory.")
#     )
#     validate(
#       need( sum(testPrograms) == 0, paste0("Cannot find: ", paste0( sub(" .*","", names( testPrograms[testPrograms > 0]) ), collapse=", ") ))
#     )
#     validate(
#       need(sum( requiredCols %in% colnames(samplesInfo())) == length(requiredCols),"Please provide a sampleInfo file with column: 'condition', 'sampleName', 'pathToFile' and 'filename' ")
#     )
    
    # print("All Ok!")
    
  })
  
  observeEvent(input$runAnalysisBtn, {
    req(samplesInfo())
    write.table(samplesInfo(),file = file.path(input$config_rootDir, "SamplesInfo.txt"), sep=",",col.names = TRUE, row.names=FALSE)
    writeConfigFile()
    runMapping()
    })
  
  runMapping = function(){
    #############################
    #   Run the Whole analysis and indicate the steps with a progress bar, by including sources
    #############################
    withProgress( message = '',min = 0, max = 1, value=0,{
      
      input$config_cutadaptRun
      
      setwd(input$config_rootDir)
      source(file.path(input$config_rootDir, "Configuration_ncRNA.R"), local = TRUE)$value
      #Cutadapt
      incProgress(amount = 0.2, detail = paste("Assembling document"))
      file.copy(from=file.path(pipelineDir,cutadaptSourceFN), to=input$config_rootDir)
      source(file.path(pipelineDir, cutadaptSourceFN), local = TRUE)$value
      
      #Mapping ncRNA
      incProgress(amount = 0.2, detail = paste("Assembling document"))
      file.copy(from=file.path(pipelineDir,ncRNAMappingSrcFN), to=input$config_rootDir)
      source(file.path(pipelineDir, ncRNAMappingSrcFN), local = TRUE)$value
      
      #MappingGenome
      incProgress(amount = 0.2, detail = paste("Assembling document"))
      file.copy(from=file.path(pipelineDir,genomeMappingSrcFN), to=input$config_rootDir)
      source(file.path(pipelineDir, genomeMappingSrcFN), local = TRUE)$value
      
      #Assembly Clustering
      incProgress(amount = 0.2, detail = paste("Assembling document"))
      file.copy(from=file.path(pipelineDir,clusteringSrcFN), to=input$config_rootDir)
      source(file.path(pipelineDir, clusteringSrcFN), local = TRUE)$value
      
      #Annotation
      incProgress(amount = 0.2, detail = paste("Assembling document"))
      file.copy(from=file.path(pipelineDir,annotationSrcFN), to=input$config_rootDir)
      source(file.path(pipelineDir, annotationSrcFN), local = TRUE)$value
      
      #DiffExp
      incProgress(amount = 0.2, detail = paste("Assembling document"))
      file.copy(from=file.path(pipelineDir,diffExpSrcFN), to=input$config_rootDir)
      source(file.path(pipelineDir, diffExpSrcFN), local = TRUE)$value
      
    })
    
  }
  
})
