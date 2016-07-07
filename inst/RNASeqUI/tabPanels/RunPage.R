fluidPage(
  sidebarPanel(width = 3,
    
    radioButtons("chooseRunScripts", label = h3("Run options"),
                 choices = list("Mapping Assembly" = 1, "Annotation" = 2, "Differential Expression" = 3), 
                 selected = 3),
    hr(),
    fluidRow(verbatimTextOutput("runScriptsOut"))
    
  ), 
  mainPanel(
    tableOutput('upl_samplesInfo'),
    hr(),
    actionButton("runAnalysisBtn", "Run Analysis!")
    )
)