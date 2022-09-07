library(shiny)
source('EnViBr.R')
library(tidyverse)
library(DT)
library(shinydashboard)

ui <- fluidPage(
# App title ----
  titlePanel("EVB"),
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input: Select a file ----
      fileInput("gene_csv", h4("Choose CSV File:"),
                multiple = FALSE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")),
      helpText("Note: the .csv file is expected to consist of one column named 'genes' with gene symbols."),
      h3(''),
      h4("OR"),
      textAreaInput("genes_text", h4("Paste gene symbols here:"), 
                    value = "Enter gene symbols..."),
      helpText("Note: Input one gene per line."),
      # Horizontal line ----
      tags$hr(),
      checkboxInput("HUGORemap", "Remap gene symbols to HUGO", value = FALSE),
      # Selection of the required enrichment procedures
      checkboxGroupInput("EnrichmentPipeline", label = h3("Select analyses"), 
                         choices = list("MSigDB" = 'msigdb', "Pathways (KEGG, Reactome, GO)" = 'rkg', "Brain regions (ABA+GTEx, Coldcuts) enrichment with brain region visualization" = 'ccuts',"Brain regoions (ABA) - fgsea" = 'aba', "Brain regions (ABA) - ABAEnrichment",
                                        "Tissues (GTEx)" = 'gtex', "Cell types (ABA)"='caba', "Subset compounds (DGIDB, DrugBank, CTD)"='drugs', "Drug ATC category enrichment (DrugBank)" ='atc', 
                                        "Signaling small molecule enrichment (cellinker)" ='sm'),
                         selected = c('ccuts')),
      
      actionButton("RUN", "Run the analysis")
      
      
    ),
    
    
    mainPanel(
      
      # Output: Data file ----
      textOutput("genes"),
      textOutput("HUGO"),
      textOutput("HUGO_res"),
      tags$br(),
      column(3, offset = 0),
      column(3, offset = 0, uiOutput("download_hugo", style = "text-align: center;")),
      column(3, offset = 0),
      #tableOutput("contents"),
      tags$br(),
      tags$hr(),
      #textOutput('status'),
      # conditionalPanel(
      #   condition = "input.breaks == 'custom'",
      #   sliderInput("breakCount", "Break Count", min = 1, max = 50, value = 10)
      # )
      dashboardBody(
        uiOutput("msigdb_ui"),
        uiOutput("rkg_ui")
      )
    )
    )
    )
      


server <- function(input, output, session) {
  showcats=20
  table_opts<-list(dom = 'Bfrtip',
                   #pageLength =  showcats,
                             buttons = 
                               list('colvis', list(
                                 extend = 'collection',
                                 buttons = list(list(extend='csv',
                                                     filename = 'enrichment_table.csv'),
                                                list(extend='excel',
                                                     filename = 'enrichment_table.excel'),
                                                list(extend='pdf',
                                                     filename= 'enrichment_table.pdf')),
                                 text = 'Download'
                               )),
                             scrollX = TRUE
  )
  # cond_genecsv<-eventReactive(input$RUN,{input$gene_csv})
  # cond_genetext<-eventReactive(input$RUN,{input$genes_text})
  # cond_hugo<-eventReactive(input$RUN,{input$HUGORemap})
  # cond_enrichment<-eventReactive(input$RUN, {input$EnrichmentPipeline})
  # input$RUN<-FALSE
  status<-reactiveValues()
  genes<-c()
  observeEvent(input$RUN, {
     
       print(input$gene_csv)
       if (!is.null(input$gene_csv)){
         genes<-read_csv(input$gene_csv$datapath)$genes
       } 
       if (!is.null(input$genes_text) & input$genes_text!=c("Enter gene symbols...")) {
         g<-strsplit(input$genes_text, '\n')[[1]]
         genes<-c(genes, g)
       }
       if (is.null(input$gene_csv) && (is.null(input$genes_text) | input$genes_text==c("Enter gene symbols..."))){
         genes<-c()
         output$genes<-renderText('No genes found in the input. Stopping.')
       }
       
       print(genes)
       if (length(genes)>0){
         output$current<-renderText('Genes read...')
         output$genes<-renderText(paste('Read', length(genes), 'genes.'))
         if (input$HUGORemap==TRUE){
           output$HUGO<-renderText('Remapping genes to HUGO...')
           genes=remap_genes(genes)
           output$HUGO_res<-renderText(paste(length(genes), 'genes were remapped.'))
           output$download_hugo <- renderUI({
             div(style = "display:inline-block;width:0%;", downloadButton("HUGO_remapped_genes", "Download remapped genes", icon = icon("download"), 
                                                                        style = " 
           flex-grow: 1;
           display: inline-block;
           background-color:#999;
           text-decoration: none;
           font-weight: 300;
           border: 1px dash transparent;
           letter-spacing: 0.98pt;
           border-color:#00245d;"))
           })
           output$HUGO_remapped_genes <- downloadHandler(
             filename = function() {
               paste("hugo_genes.txt")
             },
             content = function(file) {
               write_csv(tibble(genes=genes), file)
             }
           )
         }
      
#----MSigDB enrichment module----
         if ('msigdb' %in% input$EnrichmentPipeline){
             print('MDigDB...')
             output$current<-renderText('MSigDB enrichment...')
             msigdb_res<-enrich_msigdb(genes, showcats=showcats)
             
             
             output$msigdb_ui <- renderUI({
               check1 <- 'msigdb' %in% input$EnrichmentPipeline
               if(length(check1)==0){check1 <- F}
               if(check1){
                 tabBox(title = h1("MSigDB enrichment results"),id= "msigdb_tabs", width = 12, #height = "420px",
                        tabPanel("MSigDB enrichment plot", plotOutput("msigdb_output_plot",height = "800px")),
                        tabPanel("MSigDB enrichment table", DTOutput("msigdb_output_df"))
                 )
                 
                 
                 
               }
             
              }
             )
             output$msigdb_output_plot<-renderPlot(msigdb_res$plot)
             output$msigdb_output_df<-renderDT(msigdb_res$df, filter = "top" ,extensions = 'Buttons', options=table_opts, server = FALSE)
             output$current<-renderText('MSigDB enrichment-done.')
         }
         
#----MsigDB enrichment module - end----    
#----Pathway enrichment module----
         if ('rkg' %in% input$EnrichmentPipeline){
           print('KEGG, Reactome, GO...')
           output$current<-renderText('KEGG, Reactome, GO enrichment...')
           pathway_res<-enrich_pathways(genes, showcats=showcats)
           output$rkg_ui <- renderUI({
           check1 <- 'rkg' %in% input$EnrichmentPipeline
           if(length(check1)==0){check1 <- F}
           if(check1){
               tabBox(title = h1("KEGG, Reactome, GO enrichment results"),id= "rkg_tabs", width = 12, #height = "420px",
                      tabPanel("Reactome enrichment plot", plotOutput("reactome_output_plot",height = "1000px")),
                      tabPanel("Reactome enrichment table", DTOutput("reactome_output_df")),
                      tabPanel("Kegg enrichment plot", plotOutput("kegg_output_plot",height = "1000px")),
                      tabPanel("Kegg enrichment table", DTOutput("kegg_output_df")),
                      tabPanel("GO enrichment plot", plotOutput("go_output_plot",height = "1000px")),
                      tabPanel("GO enrichment table", DTOutput("go_output_df"))
               )}
           })
           output$reactome_output_plot<-renderPlot(pathway_res$react_plot)
           output$reactome_output_df<-renderDT(pathway_res$react@result, filter = "top" ,extensions = 'Buttons', options=table_opts, server = FALSE)
           output$kegg_output_plot<-renderPlot(pathway_res$kegg_plot)
           output$kegg_output_df<-renderDT(pathway_res$kegg@result, filter = "top" ,extensions = 'Buttons', options=table_opts, server = FALSE)
           output$go_output_plot<-renderPlot(pathway_res$go_plot)
           output$go_output_df<-renderDT(pathway_res$go@result, filter = "top" ,extensions = 'Buttons', options=table_opts, server = FALSE)
         }
         
         
#----Pathway enrichment module - end----
         
       }
  })
}
shinyApp(ui, server)


