##########################
# Boxplot Page
##########################

outputDirDescrColSize = 2

fluidPage(
  
  # Application title
  title = ("Pipeline configuration"),
  wellPanel(
  # Sidebar with controls to select the variable to plot against mpg
  # and to specify whether outliers should be included
#   sidebarPanel(width = 6,
    tabsetPanel(
      tabPanel("Directory Settings",
               br(),
               fluidRow(
                 column(3,
                 tags$span("Main Directory where all output files should be placed.")
                 ),
                 column(7,
                 textInput("config_rootDir", "Analysis Base Directory",  value = "/home/schaffrr/RNASeqUtilityTestFiles/Testing")
                 )
                 ),
               br(),
               fluidRow(
                 column(3,
                   tags$span("Path to rna ",
                   tags$a(href = "https://github.com/alexdobin/STAR", "STAR"),
                   " mapper Index file for whole genome (download genome fasta file from UCSC -> twoBit format -> convert with twoBitToFa )")
                 ),
                 column(7,
                  textInput("config_genomeIndexFilePath", "Genome STAR Index Directory",  value = "/home/simon/dbsOfflineUse/HomoSapiens/hg19/rnaStarIndex")
                 )
               ),
               br(),
               fluidRow(
                 column(3,
                   tags$span("Path to rna star mapper Index file for ncRNAs")
                 ), column(7,
                   textInput("config_genomeIndexFilePath_ncRNA", "ncRNA STAR Index Directory",  value = "/home/simon/dbsOfflineUse/HomoSapiens/hg19/ncRNA_ENSEMBL/restr400nt_extended")
                 )
               ),
               br(),
               fluidRow(
                 column(3,
                   tags$span("Path to infernal database (compiled) e.g. Rfam.cm.1_1",
                   tags$a(href = "http://eddylab.org/infernal/", "Infernal"))
                 ), column(7,
                 textInput("config_infernalDB", "Infernal DB File",  value = "/home/simon/dbsOfflineUse/HomoSapiens/hg19/infernal/Rfam.cm.1_1")
                 )
               ),
               br(),
               fluidRow(
                 column(3,
                  tags$span("Repeats masker bed file from UCSC (Table Browser) for the reuqired organism")
                 ), column(6,
                  textInput("config_repeatMaskerDir", "Path to Repeat Masker File Path",  value = "/home/simon/dbsOfflineUse/HomoSapiens/hg19/repeatMskr/")
                 ),
                 column(2,
                        textInput("config_repeatMaskerFN", "Path to Repeat Masker File Name",  value = "rpmsk_hg19.rda")
                 )
               ),
               br(),
               fluidRow(
                 column(3,
                 tags$span("Perl Path (if empty the default system path will be chosen")
                 ), column(7,
                  textInput("config_perlPath", "Perl path",  value = "")
                 )
               ),
               br(),
               fluidRow(
                 column(3,
                 tags$span("Please choose the annotation")
                 ), column(7,
                 selectInput(inputId = "config_organismForAnnotation", label="Organism", choices= c("hg19", "hg38", "mm9", "mm10"))
                 )
               ),
               br(),
                 fluidRow(
                   column(3,
                   tags$span("Sample Info File (must contain following columns: condition, sampleName, pathToFile, filename")
                   ), column(3,
                   fileInput('upl_file_csv', 'Choose Sample Info File to upload',
                             accept = c(
                               'text/csv',
                               'text/comma-separated-values',
                               'text/tab-separated-values',
                               'text/plain',
                               '.csv',
                               '.tsv'
                             )
                       )
                   ),
#                  ),
#                  fluidRow(
                   column(3,
                   radioButtons('sep', 'Separator',
                                c(Comma=',',
                                  Semicolon=';',
                                  Tab='\t'),
                                ',')
                   ),
                   column(3,
                   radioButtons('quote', 'Quote',
                                c(None='',
                                  'Double Quote'='"',
                                  'Single Quote'="'"),
                                '"')
                   )
               ),
               br(),
               fluidRow(
                 column(3,
                  tags$span("Following Model will be used for differential expression analysis")
                 ), column(7,
                 textInput("config_diffExpFormula", "Differential Expression formula", value = "~condition")
                 )
               )
      ),
      tabPanel("Pipeline Options", 
               br(),
               fluidRow(
                 column(3,
                 tags$span("This flag defines if the groups should be considered when intersecting all contigs")
                 ), column(7,
                 checkboxInput("config_grouping", "Grouping considered", value =TRUE)
               )),
               br(),
               fluidRow(
                 column(3,
                 tags$span("This is the threshold for contig assembly within the groups, it defines how many samples of a group do not need to have reads for this contig ")
                 ), column(7,
                 numericInput("config_withinGroupTH","Within Group Threshold",0, min=0, max=30, step=1)
               )),
               br(),
               fluidRow(
                 column(3,
                 tags$span("read threshold for contig assembly ")
                 ), column(7,
                 numericInput("config_read_threshold","Read Threshold",1, min=0, max=30, step=1)
               )),
               br(),
               fluidRow(
                 column(3,
                 tags$span("This parameters defines the number of nt a read has to overlap with the next in order to be clustered (from positive to negative numbers -> APART default: 0) ")
                 ), column(7,
                 numericInput("config_readOverlap_contig","Read Overlap (nt)",-1, min=-100, max=100, step=1)
               )),
               br(),
               fluidRow(
                 column(3,
                 tags$span("This is used for clustering: If 0.95: At least 95% of the reads have to be shared between a contig with the representative contig (highest read count) in order to be clustered, so in this case if there are 5% reads that are different then the contig will be kept next to the representative contig. When this value is set to 1 then contigs will be clustered to the representative contig if they share the same reads. If this is set to 0, then a contig will be deleted if it shares 1 or more reads with an representative contig! ")
                 ), column(7,
                 numericInput("config_readCompositionIdentity","Read Composition Identity",0.95, min=0, max=1, step=0.1)
               ))
      ),
      tabPanel("Cutadapt Options", 
               br(),
               fluidRow(
                 column(3,
                 tags$span("Run Cutadapt")
                 ), column(7,
                 checkboxInput("config_cutadaptRun", "Run Cutadapt", value =TRUE)
               )),
               br(),
               fluidRow(
                 column(3,
                 tags$span("Cutadapt options")
                 ), column(7,
                 textInput("config_cutadaptOptions", "Cutadapt options",  value = "-a ATCACCGACTGCCCATAGAGAGG --minimum-length 16", width=800)
               ))
      ),
      tabPanel("RNAStar Options", 
               br(),
               fluidRow(
                 column(3,
                 tags$span("RNA Star Options for mapping to small RNA genome")
                 ), column(7,
                 textInput("config_rnaStarncRNA_params", "Cutadapt options",  value = "--runThreadN 8 --outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin 16 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --alignIntronMax 1 --outFilterMultimapNmax 100 --outSAMprimaryFlag AllBestScore --outReadsUnmapped Fastx", width=800)
               )),
               br(),
               fluidRow(
                 column(3,
                 tags$span("RNA Star Options for re-mapping to the whole genome")
                 ), column(7,
                 textInput("config_rnaStarGenome_params", "Cutadapt options",  value = "--runThreadN 8 --outFilterMismatchNoverReadLmax 0.023 --outFilterMatchNmin 18 --outFilterScoreMinOverLread 0  --outFilterMatchNminOverLread 0 --alignIntronMax 1 --outFilterMultimapNmax 100 --alignEndsType EndToEnd --outSAMprimaryFlag AllBestScore --outSAMtype BAM Unsorted", width=800)
               ))
      ),
      tabPanel("Output Directories", 
               br(),
               fluidRow(
                 column(outputDirDescrColSize,
                 tags$span("ncRNA Mapping Directory Name")
                 ), column(7,
                 textInput("config_ncRNAmappingDir", "ncRNA Mapping Directory",  value = "ncRNAmapping")
               )),
               br(),
               fluidRow(
                 column(outputDirDescrColSize,
                 tags$span("Re-Mapping Directory Name")
                 ), column(7,
                 textInput("config_mappingDir", "Re-Mapping Directory",  value = "mapping")
               )),
               br(),
               fluidRow(
                 column(outputDirDescrColSize,
                 tags$span("Contig Assembly Directory Name")
                 ), column(7,
                 textInput("config_contigAssemblyDir", "ncRNA Mapping Directory",  value = "contigAssembly")
               )),
               br(),
               fluidRow(
                 column(outputDirDescrColSize,
                 tags$span("Contig Clustering Directory Name")
                 ), column(7,
                 textInput("config_contigClusterDir", "Contig Clustering Directory",  value = "contigClustering")
                )),
               br(),
               fluidRow(
                 column(outputDirDescrColSize,
                 tags$span("Read Count Directory Name")
                 ), column(7,
                 textInput("config_readCountsDir", "Read Count Directory",  value = "readCounts")
               )),
               br(),
               fluidRow(
                 column(outputDirDescrColSize,
                 tags$span("Command Log Filename")
                 ), column(7,
                 textInput("config_commandLog", "Command Log Filename",  value = "Commands")
               )),
               br(),
               fluidRow(
                 column(outputDirDescrColSize,
                 tags$span("Execution Log Filename")
                 ), column(7,
                 textInput("config_executionLog", "Execution Log Filename",  value = "executionLog")
               )),
               br(),
               fluidRow(
                 column(outputDirDescrColSize,
                 tags$span("Annotation Directory Name")
                 ), column(7,
                 textInput("config_annotationDir", "Annotation Directory",  value = "annotation")
               )),
               br(),
               fluidRow(
                 column(outputDirDescrColSize,
                 tags$span("Differential Expression Directory Name")
                 ), column(7,
                 textInput("config_diffExpReportingDir", "DiffExp Directory",  value = "de_analysis_report")
               ))
        )
    )
  )
#   ),
#   
#   mainPanel(
#     dataTableOutput(outputId = 'upl_sampleInfo_data')
#   )
  
)