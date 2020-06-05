# Server logic 
shinyServer(function(input, output, session){
  
 source("global.R")
  
  ## Displaying Master Table for the samples chosne for our study.
  
  output$chart<- renderTable({
    buf <- read.table("./www/samplesMasterTable.txt", header = TRUE, sep = "\t", quote="")
    return(buf[ , c("Experimental_Method", "Cell_Type.Tissue", "GSM", "TF.Histone_Mark", "Disease.Target_Pathway", "PMID")])
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
                labs(x= "Chipenrich | Colorectal Cancer", y= "Metric"))
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
                labs(x= "Broadenrich | Colorectal Cancer", y= "Metric"))
        
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
                labs(x= "Seq2pathway | Colorectal Cancer", y= "Metric"))
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
                labs(x= "Enrichr | Colorectal Cancer", y= "Metric"))
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
                labs(x= "GREAT | Colorectal Cancer", y= "Metric"))
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
                labs(x= "Chipenrich | Gastric Cancer", y= "Metric"))
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
                labs(x= "Broadenrich | Gastric Cancer", y= "Metric"))
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
                labs(x= "Seq2pathway | Gastric Cancer", y= "Metric"))
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
                labs(x= "Enrichr | Gastric Cancer", y= "Metric"))
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
                labs(x= "GREAT | Gastric Cancer", y= "Metric"))
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
                labs(x= "Chipenrich | Prostate Cancer", y= "Metric"))
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
                labs(x= "Broadenrich | Prostate Cancer", y= "Metric"))
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
                labs(x= "Seq2pathway | Prostate Cancer", y= "Metric"))
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
                labs(x= "Enrichr | Prostate Cancer", y= "Metric"))
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
                labs(x= "GREAT | Prostate Cancer", y= "Metric"))
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
                labs(x= "Chipenrich | Alzheimer's Disease", y= "Metric"))
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
                labs(x= "Broadenrich | Alzheimer's Disease", y= "Metric"))
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
                labs(x= "Seq2pathway | Alzheimer's Disease", y= "Metric"))
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
                labs(x= "Enrichr | Alzheimer's Disease", y= "Metric"))
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
                labs(x= "GREAT | Alzheimer's Disease", y= "Metric"))
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
      if(input$uam == 'uasn'){return(frameMe(consolidatedSensitivity))}
      
      if(input$uam == 'uasp'){return(frameMe(consolidatedSpecificity))}
      
      if(input$uam == 'uapr'){return(frameMe(consolidatedPrecision))}
      
      if(input$uam == 'uapn'){return(frameMe(consolidatedPrioritization))}
    }   
  }

  ## Output results of comparison metrics
  
  observeEvent(input$uam, {
    output$toolOut<- renderTable({
      uaOut()
    })    
  })
  
  
  ## Plotting results
  
  observeEvent(input$uam, {
    output$uaPlot <- renderPlot({
          uaPlot()
        }, height = 500, width = 900)

  })
  
})