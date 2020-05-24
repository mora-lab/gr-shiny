shinyUI(fluidPage(theme = shinytheme("sandstone"),
                titlePanel(windowTitle ="GR-Shiny | Mora Lab | GMU",
                  div(column(width = 9, h2("GR-Shiny")), 
                               column(3, img(height = 100, width = 200, src = "labLogo.png")))),
                hr(),
                h5('This application allows analyzing a user defined dataset for ranking various enrichment tools for genomic regions.
   The current suite of tools include GREAT, Enrichr, Chipenrich, Broadenrich, Seq2pathway. However, only the latter 3 are available
   as functions in R and hence form a part of this interface. The user is required to select the dataset, tool, and the metric of
   comparison. The application returns the plot and values as a result.'),
                h5('This interactive module allows for reviewing the results from our current study of benchmark data as we formulated.'),  
                br(),
                tabsetPanel(
                  tabPanel("Review Results from our Study ", fluid = TRUE, 
                           mainPanel(
                             tabsetPanel(type = "pill",
                                         tabPanel("Dataset Reference",
                                                  tags$iframe(style="height:1000px; width:100%; scrolling=yes", 
                                                              src="GSA_ChIP_Seq_Master_Table.pdf")), # Summary
                                         tabPanel("Results", 
                                                  tabsetPanel(type = "pill",
                                                              tabPanel("Plots on relative performance",
                                                                       h5("Please select a Comparison metric from the right panel."),
                                                                       conditionalPanel(condition = "input.m=='sn'", tags$img(src="Sensitivity_ggplot.jpeg", 
                                                                                                                              height="500", 
                                                                                                                              width="900",
                                                                                                                              align="center")),
                                                                       conditionalPanel(condition = "input.m=='sp'", tags$img(src="Specificity_ggplot.jpeg",
                                                                                                                              height="500", 
                                                                                                                              width="900",
                                                                                                                              align="center")),
                                                                       conditionalPanel(condition = "input.m=='pr'", tags$img(src="Prioritization_ggplot.jpeg",
                                                                                                                              height="500", 
                                                                                                                              width="900",
                                                                                                                              align="center")),
                                                                       conditionalPanel(condition = "input.m=='pn'", tags$img(src="Precision_ggplot.jpeg",
                                                                                                                              height="500", 
                                                                                                                              width="900",
                                                                                                                              align="center"))),     # Display plot
                                                              tabPanel("Plots on Overall Performance",
                                                                       h5("Please select a GSA tool and a Gold Standard dataset from the right panel."),
                                                                       plotOutput(outputId = "studyplot", width = "100%") %>% withSpinner(color="#0000FF", type = 4, size = 0.25)),
                                                              tabPanel("ROC Plot",
                                                                       tags$img(src="ROC_Plot.jpeg",
                                                                                height="500",
                                                                                width="900",
                                                                                align="center")))) # ROC Plot
                             )),
                           sidebarLayout("", fluid = TRUE,
                                         sidebarPanel(
                                           radioButtons(inputId = "d", label = "Select a Gold Standard dataset",
                                                        choices = c("Colorectal Cancer"='cc', "Prostate Cancer"='pc',  "Gastric Cancer"='gc', 
                                                                    "Alzheimer's Disease"='ad'),
                                                        selected = 'cc'),
                                           br(),
                                           radioButtons(inputId = "t", label = "Select a GSA tool",
                                                        choices = c("Chipenrich"='ce', "Broadenrich"='be',  "Seq2pathway"='sy', "Enrichr"='er', "GREAT"='gt'),
                                                        selected = 'ce'),
                                           br(),
                                           radioButtons(inputId = "m", label = "Select a Comparison Metric",
                                                        choices = c("Sensitivity"='sn', "Specificity"='sp',  "Prioritization"='pn', "Precision"='pr'),
                                                        selected = 'sn'),
                                           br(),
                                           submitButton("View")
                                         )
                           )),
                  tabPanel("Analyze Data", fluid = TRUE, mainPanel(tabsetPanel(type = "pill",
                                                                               tabPanel("Preview Data", 
                                                                                        dataTableOutput("toolOut") %>% withSpinner(color="#0000FF", type = 4, size = 0.25)),
                                                                               tabPanel("Results", 
                                                                                        plotOutput(outputId = "uaPlot", width = "100%") %>% withSpinner(color="#0000FF", type = 4, size = 0.25))
                                                                               
                  )),
                  sidebarLayout("", fluid = TRUE,
                                sidebarPanel(
                                  textInput('files', '1. Input Path'),
                                  submitButton('Submit'),
                                  
                                  helpText("For details on file structuring, visit the FAQ section."),
                                  
                                  br(),
                                  verbatimTextOutput("upload") %>% withSpinner(color="#0000FF", type = 4, size = 0.25),
                                           
                                  br(),
                                  checkboxGroupInput(inputId = "cbt", label = "2. Select GSA tool(s)",
                                                     choices = c("Chipenrich"='executeChipenrich', "Broadenrich"='executeBroadenrich',  "Seq2pathway"='executeSeq2pathway')),
                                                     
                                  br(),
                                  submitButton("Execute"),
                                  br(),
                                  verbatimTextOutput("tools") %>% withSpinner(color="#0000FF", type = 4, size = 0.25),
                                  br(),
                                  
                                  checkboxGroupInput(inputId = "cbd", label = "3. Select Gold Standard dataset(s)",
                                                     choices = c("Colorectal Cancer"='loadcc', "Prostate Cancer"='loadpc',  "Gastric Cancer"='loadgc',
                                                                 "Alzheimer's Disease"='loadad')),
                                  br(),
                                  submitButton("Load"),
                                  br(),
                                  verbatimTextOutput("disease") %>% withSpinner(color="#0000FF", type = 4, size = 0.25),
                                  br(),
                                  checkboxGroupInput(inputId = "cbm", label = "4. Select Comparison Metric(s)",
                                                     choices = c("Sensitivity"='runSensitivity', "Specificity"='runSpecificity',
                                                                 "Precision"='runPrecision',"Prioritization"='runPrioritization')),
                                  br(),
                                  submitButton("Calculate"),
                                  br(),
                                  verbatimTextOutput("metric") %>% withSpinner(color="#0000FF", type = 4, size = 0.25),
                                  br(),
                                  
                                  radioButtons(inputId = "uam", label = "5. Select a Comparison Metric to view data and plot",
                                               choices = c("Sensitivity"='uasn', "Specificity"='uasp', "Precision"='uapr', "Prioritization"='uapn'),
                                               selected = NULL),
                                  helpText("The data and plots can be viewed at the 'Preview Data' and 'Results' tabs, respectively, on the left."),
                                  br(),
                                  submitButton("View"),
                                  br()
                                  
                  ))),
                  tabPanel("Frequently Asked Questions", fluid = TRUE, mainPanel(tabPanel("", includeHTML("www/faq.html")))),
                  tabPanel("Contact Us", fluid = TRUE, mainPanel(tabPanel("", includeHTML("www/contact.html")))))
))