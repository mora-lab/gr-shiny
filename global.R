## Setting up the global environment
## Installing and Loading libraries.

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install(c("devtools", "GenomicRanges", "shiny", "shinythemes", "shinyjs", 
# "tidyr", "ggplot2", "markdown", "shinydashboard", "shinycssloaders", "seq2pathway", "chipenrich"))

library(shinycssloaders)
library(shiny)
library(shinythemes)
library(shinydashboard)
library(tidyr)
library(ggplot2)
library(markdown)
library(shinyjs)
library(devtools)
library(GenomicRanges)
library(seq2pathway)
library(chipenrich)


options(shiny.maxRequestSize = 100*1024^2) ## Setting maximum upload size to 100 Mb.
source("./R/appFunctions.R") ## loading function definitions
system("mkdir ./data/results")