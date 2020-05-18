## This function takes a number (tool index) as input and returns the prioritization data against each given disease.

calculatePrioritization <- function(tool){
  
  forPrioritization <- vector("list", length(ChIPSeqSamples)) # output values to be filled in this list
  
  for (sam in 1:length(ChIPSeqSamples))
  {
    for (dis in 1:length(diseasePools))
    {
      tryCatch({
      ranks <- list()
      x <- 1
      for (i in 1:length(eval(parse(text=(paste0(paste0(toolsResults[tool],"$"), ChIPSeqSamples[sam]))))[[1]]))
        {
          for (j in (eval(parse(text=diseasePools[[dis]]))))
          {
            if((eval(parse(text=(paste0(paste0(toolsResults[tool],"$"), ChIPSeqSamples[sam]))))[[1]])[[i]] == j)
            {
              ranks[[x]]<- i
              x <- x+1
            }
          }
        }
        findPrioritization <- as.double((ranks[[1]]/nrow((eval(parse(text=(paste0(paste0(toolsResults[tool],"$"), 
                                                                                  ChIPSeqSamples[sam]))))[1]))*100))
        forPrioritization[[sam]][dis] <- findPrioritization
      }, error = function(e){}
      , finally = {dis=dis+1}) ## The exception handling moves the pointer to the next index if 'ranks' is empty.
    }
  }
  ## Owing to the discrepnacies in the calculation of prioritization for different samples (intersections may vary), the resultant list of
  ## lists(results) may have different length components. To ensure fullness of the dataframe we shall calculate the maximum length amongst
  ## all (that should technically be the number of diseases), and structure the shorter lists to fit the dimensions of the dataframe by 
  ## explicit zeros. 

return(listToFrame(forPrioritization))
}

## (For GREAT) Note that in the resultant table, the zero entries reflect that there are no matches between the disease terms and the tool 
## output. While the empty cells mean that the output for the particular sample was not available.


calculateSensitivity <- function(tool)
{
  forSensitivity <- vector("list", length(ChIPSeqSamples)) 
  
  ## Classical Sensitivity
  
  for (sam in 1:length(ChIPSeqSamples))
  {
    for (dis in 1:length(diseasePools))
    {
      truePositives <- list()
      trueNegativesIDs <- list()
      falsePositives1IDs <- list()
      falsePositives2IDs <- list()
      falsePositives <- list()
      falseNegatives <- list()
      
      ## Tool results' subsets on the basis of statistical significance.
      greaterThan0.05 <- eval(parse(text=(paste0(paste0(toolsResults[tool],"$"), ChIPSeqSamples[sam]))))[which(eval(parse(text=(paste0(paste0(toolsResults[tool],"$"), ChIPSeqSamples[sam]))))[2] > 0.05),]
      lessThan0.05 <- eval(parse(text=(paste0(paste0(toolsResults[tool],"$"), ChIPSeqSamples[sam]))))[which(eval(parse(text=(paste0(paste0(toolsResults[tool],"$"), ChIPSeqSamples[sam]))))[2] <= 0.05),]
      trueNegativesIDs <- setdiff(greaterThan0.05[[1]], eval(parse(text=diseasePools[dis]))) ## All ids that are there in the tool result with p > 0.05 and absent in the disease pool.
      
      falsePositives1IDs <- intersect(eval(parse(text=diseasePools[dis])),greaterThan0.05[[1]])
      falsePositives2IDs <- setdiff(lessThan0.05[[1]],eval(parse(text=diseasePools[dis])))
      
      falsePositives <- c(falsePositives1IDs,falsePositives2IDs)
      truePositives <- intersect(lessThan0.05[[1]], eval(parse(text=diseasePools[dis])))
      falseNegatives <- intersect(greaterThan0.05[[1]], eval(parse(text=diseasePools[dis])))
      
      ## Results
      
      sensitivity <- length(truePositives)/(length(truePositives)+length(falseNegatives))
      
      # Store results
      
      forSensitivity[[sam]][dis] <- sensitivity
      
    }
  } 
  
  ## Let's return a list for convenience; a list of dataframes.
  
  return(listToFrame(forSensitivity))
}


calculateSpecificity <- function(tool)
{
  
  forSpecificity <- vector("list", length(ChIPSeqSamples)) 
  
  ## Classical Specificity
  
  for (sam in 1:length(ChIPSeqSamples))
  {
    for (dis in 1:length(diseasePools))
    {
      truePositives <- list()
      trueNegativesIDs <- list()
      falsePositives1IDs <- list()
      falsePositives2IDs <- list()
      falsePositives <- list()
      falseNegatives <- list()
      
      ## Tool results' subsets on the basis of statistical significance.
      greaterThan0.05 <- eval(parse(text=(paste0(paste0(toolsResults[tool],"$"), ChIPSeqSamples[sam]))))[which(eval(parse(text=(paste0(paste0(toolsResults[tool],"$"), ChIPSeqSamples[sam]))))[2] > 0.05),]
      lessThan0.05 <- eval(parse(text=(paste0(paste0(toolsResults[tool],"$"), ChIPSeqSamples[sam]))))[which(eval(parse(text=(paste0(paste0(toolsResults[tool],"$"), ChIPSeqSamples[sam]))))[2] <= 0.05),]
      trueNegativesIDs <- setdiff(greaterThan0.05[[1]], eval(parse(text=diseasePools[dis]))) ## All ids that are there in the tool result with p > 0.05 and absent in the disease pool.
      
      falsePositives1IDs <- intersect(eval(parse(text=diseasePools[dis])),greaterThan0.05[[1]])
      falsePositives2IDs <- setdiff(lessThan0.05[[1]],eval(parse(text=diseasePools[dis])))
      
      falsePositives <- c(falsePositives1IDs,falsePositives2IDs)
      truePositives <- intersect(lessThan0.05[[1]], eval(parse(text=diseasePools[dis])))
      falseNegatives <- intersect(greaterThan0.05[[1]], eval(parse(text=diseasePools[dis])))
      
      ## Results
      
      specificity <- length(trueNegativesIDs)/(length(trueNegativesIDs)+length(falsePositives))
      
      # Store results
      
      forSpecificity[[sam]][dis] <- specificity
      
    }
  } 
  
  ## Let's return a list for convenience; a list of dataframes.
  
  return(listToFrame(forSpecificity))
  
}


calculatePrecision <- function(tool)
{
  forPrecision <- vector("list", length(ChIPSeqSamples))
  
  ## Classical Precision
  
  for (sam in 1:length(ChIPSeqSamples))
  {
    for (dis in 1:length(diseasePools))
    {
      truePositives <- list()
      trueNegativesIDs <- list()
      falsePositives1IDs <- list()
      falsePositives2IDs <- list()
      falsePositives <- list()
      falseNegatives <- list()
      
      ## Tool results' subsets on the basis of statistical significance.
      greaterThan0.05 <- eval(parse(text=(paste0(paste0(toolsResults[tool],"$"), ChIPSeqSamples[sam]))))[which(eval(parse(text=(paste0(paste0(toolsResults[tool],"$"), ChIPSeqSamples[sam]))))[2] > 0.05),]
      lessThan0.05 <- eval(parse(text=(paste0(paste0(toolsResults[tool],"$"), ChIPSeqSamples[sam]))))[which(eval(parse(text=(paste0(paste0(toolsResults[tool],"$"), ChIPSeqSamples[sam]))))[2] <= 0.05),]
      trueNegativesIDs <- setdiff(greaterThan0.05[[1]], eval(parse(text=diseasePools[dis]))) ## All ids that are there in the tool result with p > 0.05 and absent in the disease pool.
      
      falsePositives1IDs <- intersect(eval(parse(text=diseasePools[dis])),greaterThan0.05[[1]])
      falsePositives2IDs <- setdiff(lessThan0.05[[1]],eval(parse(text=diseasePools[dis])))
      
      falsePositives <- c(falsePositives1IDs,falsePositives2IDs)
      truePositives <- intersect(lessThan0.05[[1]], eval(parse(text=diseasePools[dis])))
      falseNegatives <- intersect(greaterThan0.05[[1]], eval(parse(text=diseasePools[dis])))
      
      ## Results
      
      precision <- length(truePositives)/(length(truePositives)+length(falsePositives))
      
      # Store results
      
      forPrecision[[sam]][dis] <- precision
      
    }
  } 
  
  ## Let's return a list for convenience; a list of dataframes.
  
  return(listToFrame(forPrecision))
  
}



dataImportClean <- function(loc){
  
    ### Importing the master table ###
    
    setwd(".")

    ## The following code creates a list of samples for which we need to extract the BED files for.
    
    
    ChIPSeqSamples <<- list.files(loc) # Extracting files from the input directory.
    ChIPSeqSamples <<- substr(ChIPSeqSamples,1,nchar(ChIPSeqSamples)-4) # Clipping file extension to retrieve sample names only.
    
    ## Initializing list for storing BED files and the consecutive GRanges objects.
    
    samplesInBED <- list()
    
    for(i in 1:length(ChIPSeqSamples))
    {
      samplesInBED[[i]] <- read.table(paste0(loc,paste0(eval(parse(text="ChIPSeqSamples[i]")),".bed")), sep = "\t", header = FALSE)
      samplesInBED[[i]] <- samplesInBED[[i]][,1:3]
      colnames(samplesInBED[[i]]) <- c("chrom", "start", "end")
      samplesInBED[[i]] <- samplesInBED[[i]][order(samplesInBED[[i]]$chrom),]
      samplesInBED[[i]] <- GRanges(seqnames = samplesInBED[[i]]$chrom, ranges = IRanges(samplesInBED[[i]]$`start`, samplesInBED[[i]]$`end`))
      genome(samplesInBED[[i]]) <- "hg19"
    }
    
    ## Saving BED files as GRanges objects ##
    
    samplesInBED <- GRangesList(samplesInBED)
    names(samplesInBED) <- ChIPSeqSamples
    saveRDS(samplesInBED, "./data/samplesInBED.rds")
    saveRDS(ChIPSeqSamples, "./data/ChIPSeqSamples.rds")
}


diseaseTerms <- function(){
## GO Terms
  
## Import the manually downloaded GO terms (as a text file) from geneontology.org associated with each of the four
## diseases as depicted in the benchmark dataset.

  
gastricCancerGOTerms <- read.table("./data/gastricCancerGOTerms.txt", sep = "\t", header = FALSE)
prostateCancerGOTerms <- read.table("./data/prostateCancerGOTerms.txt", sep = "\t", header = FALSE)
alzheimersDiseaseGOTerms <- read.table("./data/alzheimersDiseaseGOTerms.txt", sep = "\t", header = FALSE)
colorectalCancerGOTerms <- read.table("./data/colorectalCancerGOTerms.txt", sep = "\t", header = FALSE)

## Convert data frame column as character expressions instead of factors.

gastricCancerGOTerms <- data.frame(lapply(gastricCancerGOTerms, as.character), stringsAsFactors=FALSE) 
prostateCancerGOTerms <- data.frame(lapply(prostateCancerGOTerms, as.character), stringsAsFactors=FALSE) 
alzheimersDiseaseGOTerms <- data.frame(lapply(alzheimersDiseaseGOTerms, as.character), stringsAsFactors=FALSE) 
colorectalCancerGOTerms <- data.frame(lapply(colorectalCancerGOTerms, as.character), stringsAsFactors=FALSE) 

## Extract the GO id from the full length string including the id as well as description.

for (i in 1: length(gastricCancerGOTerms$V1))
{
  gastricCancerGOTerms$V1[i] <- substr(gastricCancerGOTerms$V1[i],1,10)
}

for (i in 1: length(colorectalCancerGOTerms$V1))
{
  colorectalCancerGOTerms$V1[i] <- substr(colorectalCancerGOTerms$V1[i],1,10)
}

for (i in 1: length(prostateCancerGOTerms$V1))
{
  prostateCancerGOTerms$V1[i] <- substr(prostateCancerGOTerms$V1[i],1,10)
}

for (i in 1: length(alzheimersDiseaseGOTerms$V1))
{
  alzheimersDiseaseGOTerms$V1[i] <- substr(alzheimersDiseaseGOTerms$V1[i],1,10)
}

## KEGG Terms

## We define a list for KEGG terms for each disease here. 
## The foremost ID is the KEGG term for the disease (self-explanatory variable) and the trailing ids represent it's
## subpathways.

alzheimersDiseaseKEGGTerms <- c("hsa05010", "hsa00190", "hsa04210", "hsa04020")
gastricCancerKEGGTerms <- c("hsa05226", "hsa05120", "hsa04115", "hsa04110", "hsa04310", "hsa04151", "hsa04350", "hsa04520", "hsa04010")
prostateCancerKEGGTerms <- c("hsa05215", "hsa04110", "hsa04210", "hsa04151", "hsa05202", "hsa04115", "hsa04010", "hsa04060", "hsa00140")
colorectalCancerKEGGTerms <- c("hsa05210", "hsa04310", "hsa04110", "hsa04210", "hsa04115", "hsa04151", "hsa04010", "hsa04350", "hsa04012", "hsa04150")


## Merging KEGG and GO terms.
gastricCancerPool <- c(gastricCancerGOTerms[[1]], gastricCancerKEGGTerms)
colorectalCancerPool <- c(colorectalCancerGOTerms[[1]], colorectalCancerKEGGTerms)
alzheimersDiseasePool <- c(alzheimersDiseaseGOTerms[[1]], alzheimersDiseaseKEGGTerms)
prostateCancerPool <- c(prostateCancerGOTerms[[1]], prostateCancerKEGGTerms)


# Saving these objects
saveRDS(alzheimersDiseasePool,"./data/alzheimersDiseasePool")
saveRDS(prostateCancerPool, "./data/prostateCancerPool")
saveRDS(colorectalCancerPool, "./data/colorectalCancerPool")
saveRDS(gastricCancerPool, "./data/gastricCancerPool")

}


frameMe <- function(x)
{
  library(dplyr)
  y <- lapply(x, function(z) z %>% select(Median))
  y <- as.data.frame(y)
  colnames(y) <- toolsResults
  y$Samples <- ChIPSeqSamples
  return(y)
}


listToFrame <- function(listLists)
{
  ## The following chunk stretches the incomplete lists to the maximum length (total number of considered
  ## diseases), by filling the incomplete entries (trailing) with NA.  
  # tempList <- list()
  # for(ind in 1: length(listLists)){tempList <- append(tempList, length(listLists[[ind]]))}
  # maxNum <- max(unlist(tempList))
  healMe <-  function(x) {
    for(ind in 1: length(x))
    {
      if(length(x[[ind]]) < length(diseasePools))
      {
        x[[ind]] <- append(x[[ind]], rep(as.numeric(NA), times = (length(diseasePools)-length(x[[ind]]))))
      }
    }
    return(x)}
  
  listLists <- healMe(listLists)
  
  ## The below lines: (i) convert the list of lists to list a dataframe representing results from a tool,
  ## (ii) annotate the columns, (iii) add columns for median and sample names, and (iv) returns the 
  ## "glorified" dataframe.
  
  finalFrame <- as.data.frame(listLists) # transform to a dataframe
  finalFrame <- as.data.frame(t(finalFrame)) # transpose the data frame as the output from the function is a list
  row.names(finalFrame) <- NULL
  colnames(finalFrame) <- diseasePools
  finalFrame$Median <- apply(finalFrame, 1, median) # median value shall be the basis of plotting the results.
  finalFrame$Samples <- ChIPSeqSamples # add key attribute of sample names. this may be helpful for the purpose of joining dataframes.
  finalFrame[is.na(finalFrame)] <- 0 ## replacing NAs with zero for lists shorther than the maximum length.
  return(finalFrame)
}



plotMetrics <- function(x, metric)
{
  y <- gather(x, Tool, medianValue, -Samples)
  y$Tool <- toupper(y$Tool)
  y$Tool <- substr(y$Tool,1,nchar(y$Tool)-15)
  ggplot(data = y,
         mapping = aes(Tool, medianValue, fill=Tool, na.rm = FALSE)) +
    geom_boxplot(varwidth = TRUE) +
    labs(x= "TOOL", y= toupper(as.character(metric)))
}



## The following function takes the sensitivity and specificity values for the respective tools and draws out an ROC.

rocPlot <- function()
{
  plot(sort(1-plotSpecificity$chipenrichResultsShredded), sort(plotSensitivity$chipenrichResultsShredded), type="b", col=c("red"), 
       xlab="1-Specificity",
       ylab="Sensitivity", main = "Tools for Enrichment of Genomic Regions")
  
  lines(sort(1-plotSpecificity$broadenrichResultsShredded), sort(plotSensitivity$broadenrichResultsShredded), type="b", col=c("green"))
  lines(1-plotSpecificity[order(plotSpecificity$enrichrResultsShredded),][,4], plotSensitivity[order(plotSensitivity$enrichrResultsShredded),][,4], type="b", col=c("blue"))
  lines(1-plotSpecificity[order(plotSpecificity$seq2pathwayResultsShredded),][,3], plotSensitivity[order(plotSensitivity$seq2pathwayResultsShredded),][,3], type="b", col=c("magenta"))
  lines(1-plotSpecificity[order(plotSpecificity$greatResultsShredded),][,5], plotSensitivity[order(plotSensitivity$greatResultsShredded),][,5], type="b", col=c("lightpink"))
  
  legend(0.06, 0.65, legend=c("CHIPENRICH", "BROADENRICH", "ENRICHR", "SEQ2PATHWAY", "GREAT"),
         col=c("red", "green", "blue", "magenta", "lightpink"), lty=1, cex=0.8, box.lty=0)
}



executeSeq2pathway <- function(loc){
  ## Testing individual data in benchmark dataset with each GSA tool package
  ## Also, since several of the tools don't acknowledge the mitochondrial DNA ("chrMT") entries have to be removed from the BED files.
  
  regeneratedSamples <- list()
  
  
  ## This module holds the individual executives for the tools in question
  ## Seq2pathway
  
  seq2pathwayRun <- function(x){
    results <- list()
    results [[i]] <- runseq2pathway(x, genome = "hg19")
    return(results)
  }
  
  
  ## Execution
  
  seq2pathwayResults <- list()
  ChIPSeqSamples <- readRDS("./data/ChIPSeqSamples.rds")
  
  for (i in 1:length(ChIPSeqSamples))
  {
    regeneratedSamples[[i]] <- read_bed(paste0(loc,paste0(eval(parse(text="ChIPSeqSamples[i]")),".bed")))
    seq2pathwayResults[[i]] <- seq2pathwayRun(regeneratedSamples[[i]])
  }
  
  dir.create("./data/results/Seq2pathway")
  saveRDS(seq2pathwayResults, file = "./data/results/Seq2pathway/seq2pathwayResults")
  rm(seq2pathwayResults)
  
  ## Pruning the irrelevant results from Seq2pathway
  
  seq2pathwayResults <- readRDS("./data/results/Seq2pathway/seq2pathwayResults")
  
  for (i in 1:length(seq2pathwayResults))
  {
    seq2pathwayResultsShredded[[i]] <<- rbind(seq2pathwayResults[[i]][[i]]$gene2pathway_result.FET$GO_BP[,c(1,3)],
                                             seq2pathwayResults[[i]][[i]]$gene2pathway_result.FET$GO_CC[,c(1,3)],
                                             seq2pathwayResults[[i]][[i]]$gene2pathway_result.FET$GO_MF[,c(1,3)])
  }
  names(seq2pathwayResultsShredded) <<- as.character(ChIPSeqSamples)
  saveRDS(seq2pathwayResultsShredded, file = "./data/results/Seq2pathway/seq2pathwayResultsShredded")
  toolsResults <<- append(toolsResults, "seq2pathwayResultsShredded")
  rm(seq2pathwayResults)
}

executeChipenrich <- function(loc){
  ## Testing individual data in benchmark dataset with each GSA tool package
  ## Also, since several of the tools don't acknowledge the mitochondrial DNA ("chrMT") entries have to be removed from the BED files.
  
  regeneratedSamples <- list()
  
  ## Chipenrich
  
  chipenrichRun <- function(x){
    results <- list()
    results [[i]] <- chipenrich(peaks = x, out_name = NULL, genesets = c("GOBP", "GOCC", "GOMF", "kegg_pathway"), genome = "hg19", qc_plots = FALSE, n_cores = 1)
    return(results)
  }
  
  chipenrichResults <- list()
  ChIPSeqSamples <- readRDS("./data/ChIPSeqSamples.rds")
  
  for (i in 1:length(ChIPSeqSamples))
  {
    chipenrichResults[[i]] <- chipenrichRun(paste0(loc,paste0(eval(parse(text="ChIPSeqSamples[i]")),".bed")))
  }
  dir.create("./data/results/Chipenrich")
  saveRDS(chipenrichResults, file = "./data/results/Chipenrich/chipenrichResults")
  rm(chipenrichResults)
  
  ## Chipenrich
  
  chipenrichResults <- readRDS("./data/results/Chipenrich/chipenrichResults")
  for(i in 1:length(chipenrichResults))
  {
    chipenrichResultsShredded[[i]] <<- chipenrichResults[[i]][[i]]$results[,c(2,4)]
    
  }
  names(chipenrichResultsShredded) <<- as.character(ChIPSeqSamples)
  saveRDS(chipenrichResultsShredded, file = "./data/results/Chipenrich/chipenrichResultsShredded")
  toolsResults <<- append(toolsResults, "chipenrichResultsShredded")
  rm(chipenrichResults)
}


executeBroadenrich <- function(loc){
  ## Testing individual data in benchmark dataset with each GSA tool package
  ## Also, since several of the tools don't acknowledge the mitochondrial DNA ("chrMT") entries have to be removed from the BED files.
  
  regeneratedSamples <- list()
  
  
  ## Broadenrich
  
  broadenrichRun <- function(x){
    results <- list()
    results [[i]] <- broadenrich(peaks = x, out_name = NULL, genesets = c("GOBP", "GOCC", "GOMF", "kegg_pathway"), genome = "hg19", qc_plots = FALSE, n_cores = 1)
    return(results)
  }
  
  
  broadenrichResults <- list()
  ChIPSeqSamples <- readRDS("./data/ChIPSeqSamples.rds")
  
  for (i in 1:length(ChIPSeqSamples))
  {
    broadenrichResults[[i]] <- broadenrichRun(paste0(loc,paste0(eval(parse(text="ChIPSeqSamples[i]")),".bed")))
  }
  dir.create("./data/results/Broadenrich")
  saveRDS(broadenrichResults, file = "./data/results/Broadenrich/broadenrichResults")
  rm(broadenrichResults)
  
  ## Broadenrich
  
  broadenrichResults <- readRDS("./data/results/Broadenrich/broadenrichResults")
  for(i in 1:length(broadenrichResults))
  {
    broadenrichResultsShredded[[i]] <<- broadenrichResults[[i]][[i]]$results[,c(2,4)]
    
  }
  names(broadenrichResultsShredded) <<- as.character(ChIPSeqSamples)
  saveRDS(broadenrichResultsShredded, file = "./data/results/Broadenrich/broadenrichResultsShredded")
  toolsResults <<- append(toolsResults, "broadenrichResultsShredded") 
  rm(broadenrichResults)
}


## Loading disease definitions

loadcc <- function() { colorectalCancerPool <<- readRDS("./data/colorectalCancerPool"); diseasePools <<- append(diseasePools, "colorectalCancerPool")}
loadpc <- function() { prostateCancerPool <<- readRDS("./data/prostateCancerPool"); diseasePools <<- append(diseasePools, "prostateCancerPool")}
loadgc <- function() { gastricCancerPool <<- readRDS("./data/gastricCancerPool"); diseasePools <<- append(diseasePools, "gastricCancerPool")}
loadad <- function() { alzheimersDiseasePool <<- readRDS("./data/alzheimersDiseasePool"); diseasePools <<- append(diseasePools, "alzheimersDiseasePool")}


## Calculating metrics for all selected tools
## Prioritization

runPrioritization <- function() {
  consolidatedPrioritization <<- lapply(1:length(toolsResults), calculatePrioritization)
  
  ## For easy access, we shall name the list elements by tool names as they are in correspondence.
  
  names(consolidatedPrioritization) <<- toolsResults
  
  ## The resulting prioritization values are as follows. As noticeable, these have been calculated for each sample, 
  ## across given target pathways for each tool.
  
  saveRDS(consolidatedPrioritization, "./data/results/consolidatedPrioritization")
}


## Precision

runPrecision <- function() {
  consolidatedPrecision <<- lapply(1:length(toolsResults), calculatePrecision)
  
  ## For easy access, we shall name the list elements by tool names as they are in correspondence.
  
  names(consolidatedPrecision) <<- toolsResults
  
  ## The resulting prioritization values are as follows. As noticeable, these have been calculated for each sample, 
  ## across given target pathways for each tool.
  
  saveRDS(consolidatedPrecision, "./data/results/consolidatedPrecision")
}


##Sensitivity

runSensitivity <- function() {
  consolidatedSensitivity <<- lapply(1:length(toolsResults), calculateSensitivity)
  
  ## For easy access, we shall name the list elements by tool names as they are in correspondence.
  
  names(consolidatedSensitivity) <<- toolsResults
  
  ## The resulting prioritization values are as follows. As noticeable, these have been calculated for each sample, 
  ## across given target pathways for each tool.
  
  saveRDS(consolidatedSensitivity, "./data/results/consolidatedSensitivity")
}


##Specificity

runSpecificity <- function() {
  consolidatedSpecificity <<- lapply(1:length(toolsResults), calculateSpecificity)
  
  ## For easy access, we shall name the list elements by tool names as they are in correspondence.
  
  names(consolidatedSpecificity) <<- toolsResults
  
  ## The resulting prioritization values are as follows. As noticeable, these have been calculated for each sample, 
  ## across given target pathways for each tool.
  
  saveRDS(consolidatedSpecificity, "./data/results/consolidatedSpecificity")
}