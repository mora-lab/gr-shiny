## Setting up the global environment
## Installing and Loading libraries.

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

requiredPackages <- c("devtools", "GenomicRanges", "shiny", "shinythemes", "shinyjs",
                      "tidyr", "ggplot2", "markdown", "shinydashboard", "shinycssloaders", 
                      "seq2pathway", "chipenrich")

# update previously installed packages
if(requiredPackages %in% installed.packages()[,"Package"]) update.packages(requiredPackages) 


newPackages <- requiredPackages[!(requiredPackages %in% installed.packages()[,"Package"])]
if(length(newPackages)) BiocManager::install(newPackages) # install missing packages

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

## Creating results directory
if(!dir.exists("./data/results")) dir.create("./data/results") else NULL


## Defining variables to represent disease terms

alzheimersDiseasePool <- c()
colorectalCancerPool <- c()
gastricCancerPool <- c()
prostateCancerPool <- c()

## Defining lists to hold user-selected tools and gold-datasets

toolsResults <- c()
diseasePools <- c()


## Defining normalized results for different tools

chipenrichResultsShredded <- list()
broadenrichResultsShredded <- list()
seq2pathwayResultsShredded <- list()


## Our parent list

consolidatedSpecificity <- vector("list", length = length(toolsResults))
consolidatedSensitivity <- vector("list", length = length(toolsResults))
consolidatedPrecision <- vector("list", length = length(toolsResults))
consolidatedPrioritization <- vector("list", length = length(toolsResults))
