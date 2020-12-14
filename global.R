## Author: Shaurya Jauhari
## Last Reviewed: December 14th, 2020.
## Setting up the global environment
## Installing and Loading libraries.

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", dependencies = TRUE)


requiredPackages <- c("devtools", "Biostrings", "Rsamtools", "GenomicAlignments", "polspline", 
                      "GenomicRanges", "shiny", "shinythemes", "shinyjs", "tidyr", "ggplot2", 
                      "markdown", "shinydashboard", "shinycssloaders", "seq2pathway", "chipenrich", 
                      "htmlwidgets", "chipenrich.data", "hash")

newPackages <- requiredPackages[!(requiredPackages %in% installed.packages()[,"Package"])]
if(length(newPackages)) BiocManager::install(newPackages, ask = TRUE) # install missing packages

suppressPackageStartupMessages(library(htmlwidgets))
suppressPackageStartupMessages(library(shinycssloaders))
suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(shinythemes))
suppressPackageStartupMessages(library(shinydashboard))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(markdown))
suppressPackageStartupMessages(library(shinyjs))
suppressPackageStartupMessages(library(devtools))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(seq2pathway))
suppressPackageStartupMessages(library(chipenrich))
suppressPackageStartupMessages(library(chipenrich.data))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(polspline))
suppressPackageStartupMessages(library(hash))

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
