#Need to be started in the rootDir and 
load("Configuration.rda")
load("Annotation.rda")
setwd(rootDir)
fileToKnit = list.files(getwd(), pattern=".Rmd$")
file.copy(file.path(rootDir, fileToKnit), file.path(diffExpDir, fileToKnit) )
library(knitr)
setwd(diffExpDir)
knit2html(input=fileToKnit,output=paste0(sub(".Rmd","",fileToKnit),"_REPORT"), spin(knit = FALSE), force_v1 = TRUE)
