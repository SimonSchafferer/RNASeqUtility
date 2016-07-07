
#############
# Structure
############

#rootDirectory



# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#
library(shiny)
sourceDir = "./tabPanels/"

tabsetPanel(
  #First tab featuring the upload of the table
  tabPanel("Configuration", 
           source(file.path(sourceDir, "ConfigurationPageUI.R"))$value
  ),
  tabPanel("Run Pipeline", 
           source(file.path(sourceDir, "RunPage.R"))$value
  )

)