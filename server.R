# Server logic 
shinyServer(function(input, output, session){
  
  source("global.R")
  
  
  ## Display the contents of the table in the main panel.  
  output$contents <- renderTable({
    if (is.null(input$bds))
      return(NULL)
    head(read.table(input$bds$datapath, header = input$header, sep = "\t", quote = ""))
  })
  
  ## Display the summary of the table.
  output$summary <- renderPrint({
    summary(input$bds)
  })
  
  ## Aggregating plots from the radio button inputs from user- tool and disease.
  forPlot <- function(){
    
    
    if (!is.null(input$d) && !is.null(input$t)) ## Check for valid inputs
    {
      
      ## Import results for each metric from external, manually curated files.
      
      samplenames <- readRDS("www/sample_names")
      precision_data <- read.table("www/Precision_Table.txt", sep = "\t", header = TRUE, quote = "")
      prioritization_data <- read.table("www/Prioritization_Table.txt", sep = "\t", header = TRUE, quote = "")
      sensitivity_data <- read.table("www/Sensitivity_Table.txt", sep = "\t", header = TRUE, quote = "")
      specificity_data <- read.table("www/Specificity_Table.txt", sep = "\t", header = TRUE, quote = "")
      
      ## Creating all possible combinations for the tool and disease inputs from the user, 4 diseases * 5 tools, i.e. 20 precisely. 
      ## Colorectal Cancer and Chipenrich
      
      if(input$d == 'cc' && input$t == 'ce')
      {
        db_sn <- sensitivity_data[,which(grepl("Chipenrich_Colorectal", colnames(sensitivity_data)))]
        db_sp <- specificity_data[,which(grepl("Chipenrich_Colorectal", colnames(specificity_data)))]
        db_pr <- prioritization_data[,which(grepl("Chipenrich_Colorectal", colnames(prioritization_data)))]
        db_pn <- precision_data[,which(grepl("Chipenrich_Colorectal", colnames(precision_data)))]
        db <- as.data.frame(cbind(db_sn,db_sp,db_pr,db_pn))
        colnames(db) <- c("Sensitivity", "Specificity", "Precision", "Prioritization")
        db$Samples <- samplenames
        db_gather <- gather(db, Metric, Value, -Samples)
        print(ggplot(data = db_gather,
                     mapping = aes(Metric, Value, fill=Metric, na.rm = FALSE)) +
                geom_boxplot(varwidth = TRUE) +
                labs(x= "Chipenrich", y= "Colorectal Cancer"))
      }
      
      ## Colorectal Cancer and Broadenrich
      if(input$d == 'cc' && input$t == 'be')
      {
        db_sn <- sensitivity_data[,which(grepl("Broadenrich_Colorectal", colnames(sensitivity_data)))]
        db_sp <- specificity_data[,which(grepl("Broadenrich_Colorectal", colnames(specificity_data)))]
        db_pr <- prioritization_data[,which(grepl("Broadenrich_Colorectal", colnames(prioritization_data)))]
        db_pn <- precision_data[,which(grepl("Broadenrich_Colorectal", colnames(precision_data)))]
        db <- as.data.frame(cbind(db_sn,db_sp,db_pr,db_pn))
        colnames(db) <- c("Sensitivity", "Specificity", "Precision", "Prioritization")
        db$Samples <- samplenames
        db_gather <- gather(db, Metric, Value, -Samples)
        print(ggplot(data = db_gather,
                     mapping = aes(Metric, Value, fill=Metric, na.rm = FALSE)) +
                geom_boxplot(varwidth = TRUE) +
                labs(x= "Broadenrich", y= "Colorectal Cancer"))
        
      }
      
      ## Colorectal Cancer and Seq2pathway
      if(input$d == 'cc' && input$t == 'sy')
      {
        db_sn <- sensitivity_data[,which(grepl("Seq2pathway_Colorectal", colnames(sensitivity_data)))]
        db_sp <- specificity_data[,which(grepl("Seq2pathway_Colorectal", colnames(specificity_data)))]
        db_pr <- prioritization_data[,which(grepl("Seq2pathway_Colorectal", colnames(prioritization_data)))]
        db_pn <- precision_data[,which(grepl("Seq2pathway_Colorectal", colnames(precision_data)))]
        db <- as.data.frame(cbind(db_sn,db_sp,db_pr,db_pn))
        colnames(db) <- c("Sensitivity", "Specificity", "Precision", "Prioritization")
        db$Samples <- samplenames
        db_gather <- gather(db, Metric, Value, -Samples)
        print(ggplot(data = db_gather,
                     mapping = aes(Metric, Value, fill=Metric, na.rm = FALSE)) +
                geom_boxplot(varwidth = TRUE) +
                labs(x= "Seq2pathway", y= "Colorectal Cancer"))
      }
      
      ## Colorectal Cancer and Enrichr
      if(input$d == 'cc' && input$t == 'er')
      {
        db_sn <- sensitivity_data[,which(grepl("Enrichr_Colorectal", colnames(sensitivity_data)))]
        db_sp <- specificity_data[,which(grepl("Enrichr_Colorectal", colnames(specificity_data)))]
        db_pr <- prioritization_data[,which(grepl("Enrichr_Colorectal", colnames(prioritization_data)))]
        db_pn <- precision_data[,which(grepl("Enrichr_Colorectal", colnames(precision_data)))]
        db <- as.data.frame(cbind(db_sn,db_sp,db_pr,db_pn))
        colnames(db) <- c("Sensitivity", "Specificity", "Precision", "Prioritization")
        db$Samples <- samplenames
        db_gather <- gather(db, Metric, Value, -Samples)
        print(ggplot(data = db_gather,
                     mapping = aes(Metric, Value, fill=Metric, na.rm = FALSE)) +
                geom_boxplot(varwidth = TRUE) +
                labs(x= "Enrichr", y= "Colorectal Cancer"))
      }
      
      ## Colorectal Cancer and GREAT
      if(input$d == 'cc' && input$t == 'gt')
      {
        db_sn <- sensitivity_data[,which(grepl("GREAT_Colorectal", colnames(sensitivity_data)))]
        db_sp <- specificity_data[,which(grepl("GREAT_Colorectal", colnames(specificity_data)))]
        db_pr <- prioritization_data[,which(grepl("GREAT_Colorectal", colnames(prioritization_data)))]
        db_pn <- precision_data[,which(grepl("GREAT_Colorectal", colnames(precision_data)))]
        db <- as.data.frame(cbind(db_sn,db_sp,db_pr,db_pn))
        colnames(db) <- c("Sensitivity", "Specificity", "Precision", "Prioritization")
        db$Samples <- samplenames
        db_gather <- gather(db, Metric, Value, -Samples)
        print(ggplot(data = db_gather,
                     mapping = aes(Metric, Value, fill=Metric, na.rm = FALSE)) +
                geom_boxplot(varwidth = TRUE) +
                labs(x= "GREAT", y= "Colorectal Cancer"))
      }
      
      ## Gastric Cancer and Chipenrich
      if(input$d == 'gc' && input$t == 'ce')
      {
        db_sn <- sensitivity_data[,which(grepl("Chipenrich_Gastric", colnames(sensitivity_data)))]
        db_sp <- specificity_data[,which(grepl("Chipenrich_Gastric", colnames(specificity_data)))]
        db_pr <- prioritization_data[,which(grepl("Chipenrich_Gastric", colnames(prioritization_data)))]
        db_pn <- precision_data[,which(grepl("Chipenrich_Gastric", colnames(precision_data)))]
        db <- as.data.frame(cbind(db_sn,db_sp,db_pr,db_pn))
        colnames(db) <- c("Sensitivity", "Specificity", "Precision", "Prioritization")
        db$Samples <- samplenames
        db_gather <- gather(db, Metric, Value, -Samples)
        print(ggplot(data = db_gather,
                     mapping = aes(Metric, Value, fill=Metric, na.rm = FALSE)) +
                geom_boxplot(varwidth = TRUE) +
                labs(x= "Chipenrich", y= "Gastric Cancer"))
      }
      
      ## Gastric Cancer and Broadenrich
      if(input$d == 'gc' && input$t == 'be')
      {
        db_sn <- sensitivity_data[,which(grepl("Broadenrich_Gastric", colnames(sensitivity_data)))]
        db_sp <- specificity_data[,which(grepl("Broadenrich_Gastric", colnames(specificity_data)))]
        db_pr <- prioritization_data[,which(grepl("Broadenrich_Gastric", colnames(prioritization_data)))]
        db_pn <- precision_data[,which(grepl("Broadenrich_Gastric", colnames(precision_data)))]
        db <- as.data.frame(cbind(db_sn,db_sp,db_pr,db_pn))
        colnames(db) <- c("Sensitivity", "Specificity", "Precision", "Prioritization")
        db$Samples <- samplenames
        db_gather <- gather(db, Metric, Value, -Samples)
        print(ggplot(data = db_gather,
                     mapping = aes(Metric, Value, fill=Metric, na.rm = FALSE)) +
                geom_boxplot(varwidth = TRUE) +
                labs(x= "Broadenrich", y= "Gastric Cancer"))
      }
      
      ## Gastric Cancer and Seq2pathway
      if(input$d == 'gc' && input$t == 'sy')
      {
        db_sn <- sensitivity_data[,which(grepl("Seq2pathway_Gastric", colnames(sensitivity_data)))]
        db_sp <- specificity_data[,which(grepl("Seq2pathway_Gastric", colnames(specificity_data)))]
        db_pr <- prioritization_data[,which(grepl("Seq2pathway_Gastric", colnames(prioritization_data)))]
        db_pn <- precision_data[,which(grepl("Seq2pathway_Gastric", colnames(precision_data)))]
        db <- as.data.frame(cbind(db_sn,db_sp,db_pr,db_pn))
        colnames(db) <- c("Sensitivity", "Specificity", "Precision", "Prioritization")
        db$Samples <- samplenames
        db_gather <- gather(db, Metric, Value, -Samples)
        print(ggplot(data = db_gather,
                     mapping = aes(Metric, Value, fill=Metric, na.rm = FALSE)) +
                geom_boxplot(varwidth = TRUE) +
                labs(x= "Seq2pathway", y= "Gastric Cancer"))
      }
      
      ## Gastric Cancer and Enrichr
      if(input$d == 'gc' && input$t == 'er')
      {
        db_sn <- sensitivity_data[,which(grepl("Enrichr_Gastric", colnames(sensitivity_data)))]
        db_sp <- specificity_data[,which(grepl("Enrichr_Gastric", colnames(specificity_data)))]
        db_pr <- prioritization_data[,which(grepl("Enrichr_Gastric", colnames(prioritization_data)))]
        db_pn <- precision_data[,which(grepl("Enrichr_Gastric", colnames(precision_data)))]
        db <- as.data.frame(cbind(db_sn,db_sp,db_pr,db_pn))
        colnames(db) <- c("Sensitivity", "Specificity", "Precision", "Prioritization")
        db$Samples <- samplenames
        db_gather <- gather(db, Metric, Value, -Samples)
        print(ggplot(data = db_gather,
                     mapping = aes(Metric, Value, fill=Metric, na.rm = FALSE)) +
                geom_boxplot(varwidth = TRUE) +
                labs(x= "Enrichr", y= "Gastric Cancer"))
      }
      
      ## Gastric Cancer and GREAT
      if(input$d == 'gc' && input$t == 'gt')
      {
        db_sn <- sensitivity_data[,which(grepl("GREAT_Gastric", colnames(sensitivity_data)))]
        db_sp <- specificity_data[,which(grepl("GREAT_Gastric", colnames(specificity_data)))]
        db_pr <- prioritization_data[,which(grepl("GREAT_Gastric", colnames(prioritization_data)))]
        db_pn <- precision_data[,which(grepl("GREAT_Gastric", colnames(precision_data)))]
        db <- as.data.frame(cbind(db_sn,db_sp,db_pr,db_pn))
        colnames(db) <- c("Sensitivity", "Specificity", "Precision", "Prioritization")
        db$Samples <- samplenames
        db_gather <- gather(db, Metric, Value, -Samples)
        print(ggplot(data = db_gather,
                     mapping = aes(Metric, Value, fill=Metric, na.rm = FALSE)) +
                geom_boxplot(varwidth = TRUE) +
                labs(x= "GREAT", y= "Gastric Cancer"))
      }
      
      ## Prostate Cancer and Chipenrich
      if(input$d == 'pc' && input$t == 'ce')
      {
        db_sn <- sensitivity_data[,which(grepl("Chipenrich_Prostate", colnames(sensitivity_data)))]
        db_sp <- specificity_data[,which(grepl("Chipenrich_Prostate", colnames(specificity_data)))]
        db_pr <- prioritization_data[,which(grepl("Chipenrich_Prostate", colnames(prioritization_data)))]
        db_pn <- precision_data[,which(grepl("Chipenrich_Prostate", colnames(precision_data)))]
        db <- as.data.frame(cbind(db_sn,db_sp,db_pr,db_pn))
        colnames(db) <- c("Sensitivity", "Specificity", "Precision", "Prioritization")
        db$Samples <- samplenames
        db_gather <- gather(db, Metric, Value, -Samples)
        print(ggplot(data = db_gather,
                     mapping = aes(Metric, Value, fill=Metric, na.rm = FALSE)) +
                geom_boxplot(varwidth = TRUE) +
                labs(x= "Chipenrich", y= "Prostate Cancer"))
      }
      
      ## Prostate Cancer and Broadenrich
      if(input$d == 'pc' && input$t == 'be')
      {
        db_sn <- sensitivity_data[,which(grepl("Broadenrich_Prostate", colnames(sensitivity_data)))]
        db_sp <- specificity_data[,which(grepl("Broadenrich_Prostate", colnames(specificity_data)))]
        db_pr <- prioritization_data[,which(grepl("Broadenrich_Prostate", colnames(prioritization_data)))]
        db_pn <- precision_data[,which(grepl("Broadenrich_Prostate", colnames(precision_data)))]
        db <- as.data.frame(cbind(db_sn,db_sp,db_pr,db_pn))
        colnames(db) <- c("Sensitivity", "Specificity", "Precision", "Prioritization")
        db$Samples <- samplenames
        db_gather <- gather(db, Metric, Value, -Samples)
        print(ggplot(data = db_gather,
                     mapping = aes(Metric, Value, fill=Metric, na.rm = FALSE)) +
                geom_boxplot(varwidth = TRUE) +
                labs(x= "Broadenrich", y= "Prostate Cancer"))
      }
      
      ## Prostate Cancer and Seq2pathway
      if(input$d == 'pc' && input$t == 'sy')
      {
        db_sn <- sensitivity_data[,which(grepl("Seq2pathway_Prostate", colnames(sensitivity_data)))]
        db_sp <- specificity_data[,which(grepl("Seq2pathway_Prostate", colnames(specificity_data)))]
        db_pr <- prioritization_data[,which(grepl("Seq2pathway_Prostate", colnames(prioritization_data)))]
        db_pn <- precision_data[,which(grepl("Seq2pathway_Prostate", colnames(precision_data)))]
        db <- as.data.frame(cbind(db_sn,db_sp,db_pr,db_pn))
        colnames(db) <- c("Sensitivity", "Specificity", "Precision", "Prioritization")
        db$Samples <- samplenames
        db_gather <- gather(db, Metric, Value, -Samples)
        print(ggplot(data = db_gather,
                     mapping = aes(Metric, Value, fill=Metric, na.rm = FALSE)) +
                geom_boxplot(varwidth = TRUE) +
                labs(x= "Seq2pathway", y= "Prostate Cancer"))
      }
      
      ## Prostate Cancer and Enrichr
      if(input$d == 'pc' && input$t == 'er')
      {
        db_sn <- sensitivity_data[,which(grepl("Enrichr_Prostate", colnames(sensitivity_data)))]
        db_sp <- specificity_data[,which(grepl("Enrichr_Prostate", colnames(specificity_data)))]
        db_pr <- prioritization_data[,which(grepl("Enrichr_Prostate", colnames(prioritization_data)))]
        db_pn <- precision_data[,which(grepl("Enrichr_Prostate", colnames(precision_data)))]
        db <- as.data.frame(cbind(db_sn,db_sp,db_pr,db_pn))
        colnames(db) <- c("Sensitivity", "Specificity", "Precision", "Prioritization")
        db$Samples <- samplenames
        db_gather <- gather(db, Metric, Value, -Samples)
        print(ggplot(data = db_gather,
                     mapping = aes(Metric, Value, fill=Metric, na.rm = FALSE)) +
                geom_boxplot(varwidth = TRUE) +
                labs(x= "Enrichr", y= "Prostate Cancer"))
      }
      
      ## Prostate Cancer and GREAT
      if(input$d == 'pc' && input$t == 'gt')
      {
        db_sn <- sensitivity_data[,which(grepl("GREAT_Prostate", colnames(sensitivity_data)))]
        db_sp <- specificity_data[,which(grepl("GREAT_Prostate", colnames(specificity_data)))]
        db_pr <- prioritization_data[,which(grepl("GREAT_Prostate", colnames(prioritization_data)))]
        db_pn <- precision_data[,which(grepl("GREAT_Prostate", colnames(precision_data)))]
        db <- as.data.frame(cbind(db_sn,db_sp,db_pr,db_pn))
        colnames(db) <- c("Sensitivity", "Specificity", "Precision", "Prioritization")
        db$Samples <- samplenames
        db_gather <- gather(db, Metric, Value, -Samples)
        print(ggplot(data = db_gather,
                     mapping = aes(Metric, Value, fill=Metric, na.rm = FALSE)) +
                geom_boxplot(varwidth = TRUE) +
                labs(x= "GREAT", y= "Prostate Cancer"))
      }
      
      ## Alzheimer's Disease and Chipenrich
      if(input$d == 'ad' && input$t == 'ce')
      {
        db_sn <- sensitivity_data[,which(grepl("Chipenrich_Alzheimer", colnames(sensitivity_data)))]
        db_sp <- specificity_data[,which(grepl("Chipenrich_Alzheimer", colnames(specificity_data)))]
        db_pr <- prioritization_data[,which(grepl("Chipenrich_Alzheimer", colnames(prioritization_data)))]
        db_pn <- precision_data[,which(grepl("Chipenrich_Alzheimer", colnames(precision_data)))]
        db <- as.data.frame(cbind(db_sn,db_sp,db_pr,db_pn))
        colnames(db) <- c("Sensitivity", "Specificity", "Precision", "Prioritization")
        db$Samples <- samplenames
        db_gather <- gather(db, Metric, Value, -Samples)
        print(ggplot(data = db_gather,
                     mapping = aes(Metric, Value, fill=Metric, na.rm = FALSE)) +
                geom_boxplot(varwidth = TRUE) +
                labs(x= "Chipenrich", y= "Alzheimer's Disease"))
      }
      
      ## Alzheimer's Disease and Broadenrich
      if(input$d == 'ad' && input$t == 'be')
      {
        db_sn <- sensitivity_data[,which(grepl("Broadenrich_Alzheimer", colnames(sensitivity_data)))]
        db_sp <- specificity_data[,which(grepl("Broadenrich_Alzheimer", colnames(specificity_data)))]
        db_pr <- prioritization_data[,which(grepl("Broadenrich_Alzheimer", colnames(prioritization_data)))]
        db_pn <- precision_data[,which(grepl("Broadenrich_Alzheimer", colnames(precision_data)))]
        db <- as.data.frame(cbind(db_sn,db_sp,db_pr,db_pn))
        colnames(db) <- c("Sensitivity", "Specificity", "Precision", "Prioritization")
        db$Samples <- samplenames
        db_gather <- gather(db, Metric, Value, -Samples)
        print(ggplot(data = db_gather,
                     mapping = aes(Metric, Value, fill=Metric, na.rm = FALSE)) +
                geom_boxplot(varwidth = TRUE) +
                labs(x= "Broadenrich", y= "Alzheimer's Disease"))
      }
      
      ## Alzheimer's Disease and Seq2pathway
      if(input$d == 'ad' && input$t == 'sy')
      {
        db_sn <- sensitivity_data[,which(grepl("Seq2pathway_Alzheimer", colnames(sensitivity_data)))]
        db_sp <- specificity_data[,which(grepl("Seq2pathway_Alzheimer", colnames(specificity_data)))]
        db_pr <- prioritization_data[,which(grepl("Seq2pathway_Alzheimer", colnames(prioritization_data)))]
        db_pn <- precision_data[,which(grepl("Seq2pathway_Alzheimer", colnames(precision_data)))]
        db <- as.data.frame(cbind(db_sn,db_sp,db_pr,db_pn))
        colnames(db) <- c("Sensitivity", "Specificity", "Precision", "Prioritization")
        db$Samples <- samplenames
        db_gather <- gather(db, Metric, Value, -Samples)
        print(ggplot(data = db_gather,
                     mapping = aes(Metric, Value, fill=Metric, na.rm = FALSE)) +
                geom_boxplot(varwidth = TRUE) +
                labs(x= "Seq2pathway", y= "Alzheimer's Disease"))
      }
      
      ## Alzheimer's Disease and Enrichr
      if(input$d == 'ad' && input$t == 'er')
      {
        db_sn <- sensitivity_data[,which(grepl("Enrichr_Alzheimer", colnames(sensitivity_data)))]
        db_sp <- specificity_data[,which(grepl("Enrichr_Alzheimer", colnames(specificity_data)))]
        db_pr <- prioritization_data[,which(grepl("Enrichr_Alzheimer", colnames(prioritization_data)))]
        db_pn <- precision_data[,which(grepl("Enrichr_Alzheimer", colnames(precision_data)))]
        db <- as.data.frame(cbind(db_sn,db_sp,db_pr,db_pn))
        colnames(db) <- c("Sensitivity", "Specificity", "Precision", "Prioritization")
        db$Samples <- samplenames
        db_gather <- gather(db, Metric, Value, -Samples)
        print(ggplot(data = db_gather,
                     mapping = aes(Metric, Value, fill=Metric, na.rm = FALSE)) +
                geom_boxplot(varwidth = TRUE) +
                labs(x= "Enrichr", y= "Alzheimer's Disease"))
      }
      
      ## Alzheimer's Disease and GREAT
      if(input$d == 'ad' && input$t == 'gt')
      {
        db_sn <- sensitivity_data[,which(grepl("GREAT_Alzheimer", colnames(sensitivity_data)))]
        db_sp <- specificity_data[,which(grepl("GREAT_Alzheimer", colnames(specificity_data)))]
        db_pr <- prioritization_data[,which(grepl("GREAT_Alzheimer", colnames(prioritization_data)))]
        db_pn <- precision_data[,which(grepl("GREAT_Alzheimer", colnames(precision_data)))]
        db <- as.data.frame(cbind(db_sn,db_sp,db_pr,db_pn))
        colnames(db) <- c("Sensitivity", "Specificity", "Precision", "Prioritization")
        db$Samples <- samplenames
        db_gather <- gather(db, Metric, Value, -Samples)
        print(ggplot(data = db_gather,
                     mapping = aes(Metric, Value, fill=Metric, na.rm = FALSE)) +
                geom_boxplot(varwidth = TRUE) +
                labs(x= "GREAT", y= "Alzheimer's Disease"))
      }
    }
  }
  
  
  ## Deploying function for plotting.
  
  output$studyplot <- renderPlot(
    {
      forPlot()
    }, height = 500, width = 900)
  
  
  dataFolder <- NULL ## directory path
  ChIPSeqSamples <- c() ## list of sample names
  
  ## Validating path and confirming upload
  
  observeEvent(input$files, {
    output$upload <- renderText({
        dataFolder <<- input$files  
        dataImportClean(dataFolder)
        return("Samples retrieved successfully.")
      })
  })
  
  ## Executing tools as per user selection 
  
  observeEvent(input$cbt, {
    output$tools <- renderText({
      sapply(input$cbt, do.call, args = list(dataFolder))
      return("Done. Results saved.")
    })
})
  
  
  ## Activating disease definitions
  observeEvent(input$cbd, {
    output$disease <- renderText({
      diseaseTerms()
      sapply(input$cbd, do.call, args = list())
      return("Disease definitions loaded.")
      })
  })
  
  ## Calculating comparison metrics
  observeEvent(input$cbm, {
    output$metric <- renderText({
      sapply(input$cbm, do.call, args = list())
      return("Done. Results saved.")
      })
  })
  
  ## Define plots for each metric
  
  uaPlot <- function(){
    if (!is.null(input$uam)) ## Check for valid inputs
    {
      if(input$uam == 'uasn'){
        consolidatedSensitivity <<- lapply(consolidatedSensitivity, function(x) as.data.frame(x))
        plotMetrics(frameMe(consolidatedSensitivity), "Sensitivity")}
      
      if(input$uam == 'uasp'){
        consolidatedSpecificity <<- lapply(consolidatedSpecificity, function(x) as.data.frame(x))
        plotMetrics(frameMe(consolidatedSpecificity), "Specificity")}
      
      if(input$uam == 'uapr'){
        consolidatedPrecision <<- lapply(consolidatedPrecision, function(x) as.data.frame(x))
        plotMetrics(frameMe(consolidatedPrecision), "Precision")}
      
      if(input$uam == 'uapn'){
        consolidatedPrioritization <<- lapply(consolidatedPrioritization, function(x) as.data.frame(x))
        plotMetrics(frameMe(consolidatedPrioritization), "Prioritization")}
      
    }   
    
  }
  
  
  uaOut <- function(){
    if (!is.null(input$uam)) ## Check for valid inputs
    {
      if(input$uam == 'uasn'){consolidatedSensitivity}
      
      if(input$uam == 'uasp'){consolidatedSpecificity}
      
      if(input$uam == 'uapr'){consolidatedPrecision}
      
      if(input$uam == 'uapr'){consolidatedPrioritization}
      
    }   
  }

  ## Output results of comparison metrics
  
  observeEvent(input$uam, {
    output$toolOut<- renderDataTable({
      uaOut()
    })    
  })
  
  
  ## Plotting results
  
  observeEvent(input$uam, {
    output$uaPlot <- renderPlot(
      {
        if(is.null(input$uam)){
          return(NULL)
        }
        else{
          uaPlot()
        }
        
      }, height = 500, width = 900)

  })
  
})