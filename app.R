library(shiny)
source('EnViBr.R')
library(tidyverse)
library(DT)
library(shinydashboard)
#----List of things to improve---
#----1. Debug why cerebellum not shown on enrichment plots----
#----2. Add brain visualization for the standard ABA enrichment----
#----3. Dockerize----


withConsoleRedirect <- function(containerId, expr) {
  # Change type="output" to type="message" to catch stderr
  # (messages, warnings, and errors) instead of stdout.
  txt <- capture.output(results <- expr, type = "output")
  if (length(txt) > 0) {
    insertUI(paste0("#", containerId), where = "beforeEnd",
             ui = paste0(txt, "\n", collapse = "")
    )
  }
  results
}
ui <- fluidPage(
# App title ----
  titlePanel("EVB"),
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input----
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
      tags$hr(),
      checkboxInput("HUGORemap", "Remap gene symbols to HUGO", value = FALSE),
      # Selection of the required enrichment procedures
      checkboxGroupInput("EnrichmentPipeline", label = h3("Select analyses"), 
                         choices = list("MSigDB" = 'msigdb', "Pathways (KEGG, Reactome, GO)" = 'rkg', "Brain regions (ABA, Coldcuts) fgsea enrichment with visualization" = 'ccuts',"Brain regoions (ABA) - fgsea" = 'aba_fgsea', "Brain regions (ABA) - ABAEnrichment" = 'aba_std',
                                        "Tissues (GTEx)" = 'gtex', "Cell types (ABA)"='caba', "Subset compounds (DGIDB, DrugBank, CTD)"='drugs', "Drug target ATC category enrichment (DrugBank)" ='atc', 
                                        "Signaling small molecule analysis (cellinker)" ='sm'),
                         selected = c('ccuts')),
      
      actionButton("RUN", "Run the analysis")
      
      
    ),
    
    
    mainPanel(
      
      # Output ----
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
     # fluidPage(shinyjs::useShinyjs(),
     # pre(id="console")),
      dashboardBody(id='dashbd',
        uiOutput("msigdb_ui"),
        uiOutput("rkg_ui"),
        uiOutput("ccuts_ui"),
        uiOutput("aba_fgsea_ui"),
        uiOutput("aba_std_ui"),
        uiOutput("gtex_ui"),
        uiOutput("caba_ui"),
        uiOutput("drugs_ui"),
        uiOutput('atc_ui'),
        uiOutput('sm_ui')
      )
    )
    )
    )
      


server <- function(input, output, session) {
  showcats=20
  brainplot_hight="700px"
  border_style="border: 10px solid DarkGrey;"
  table_opts<-list(dom = 'Bfrtip',
                   #pageLength =  showcats,
                             buttons = 
                               list('colvis', list(
                                 extend = 'collection',
                                 buttons = list(list(extend='csv',
                                                     filename = 'enrichment_table'),
                                                list(extend='excel',
                                                     filename = 'enrichment_table'),
                                                list(extend='pdf',
                                                     filename= 'enrichment_table')),
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
  
  # global <- reactiveValues(refresh = FALSE)
  # observe({
  #   if(input$RUN) isolate(global$refresh <- TRUE)
  # })
  # output$dt <- renderUI({
  #   if(global$refresh) return()
  # })
  # observeEvent(input$RUN, {
  #   shinyjs::hideElement('gtex_ui')
  #  # shinyjs::showElement("")
  # })
  
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
           showModal(modalDialog("Running MSigDB enrichment...", footer=NULL))
           # withConsoleRedirect("console", {
           #   print('MSigDB...')
           # })
             output$current<-renderText('MSigDB enrichment...')
             msigdb_res<-enrich_msigdb(genes, showcats=showcats)
             
             
             output$msigdb_ui <- renderUI({
               check1 <- 'msigdb' %in% input$EnrichmentPipeline
               if(length(check1)==0){check1 <- F}
               if(check1){
                 fluidRow(style=border_style,tabBox(title = h1("MSigDB enrichment results"),id= "msigdb_tabs", width = 12, #height = "420px",
                        tabPanel("MSigDB enrichment plot", plotOutput("msigdb_output_plot",height = "800px")),
                        tabPanel("MSigDB enrichment table", DTOutput("msigdb_output_df"))
                 ))
                 
                 
                 
               }
             
              }
             )
             output$msigdb_output_plot<-renderPlot(msigdb_res$plot)
             output$msigdb_output_df<-renderDT(msigdb_res$df, filter = "top" ,extensions = 'Buttons', options=table_opts, server = FALSE)
             output$current<-renderText('MSigDB enrichment-done.')
         removeModal()}
         
#----MsigDB enrichment module - end----    
#----Pathway enrichment module----
         if ('rkg' %in% input$EnrichmentPipeline){
           showModal(modalDialog("Running KEGG, Reactome, GO enrichments...", footer=NULL))
           print('KEGG, Reactome, GO...')
           output$current<-renderText('KEGG, Reactome, GO enrichment...')
           pathway_res<-enrich_pathways(genes, showcats=showcats)
           output$rkg_ui <- renderUI({
           check1 <- 'rkg' %in% input$EnrichmentPipeline
           if(length(check1)==0){check1 <- F}
           if(check1){
             fluidRow(style=border_style,tabBox(title = h1("KEGG, Reactome, GO enrichment results"),id= "rkg_tabs", width = 12, #height = "420px",
                      tabPanel("Reactome enrichment plot", plotOutput("reactome_output_plot",height = "1000px")),
                      tabPanel("Reactome enrichment table", DTOutput("reactome_output_df")),
                      tabPanel("Kegg enrichment plot", plotOutput("kegg_output_plot",height = "1000px")),
                      tabPanel("Kegg enrichment table", DTOutput("kegg_output_df")),
                      tabPanel("GO enrichment plot", plotOutput("go_output_plot",height = "1000px")),
                      tabPanel("GO enrichment table", DTOutput("go_output_df"))
               ))}
           })
           output$reactome_output_plot<-renderPlot(pathway_res$react_plot)
           output$reactome_output_df<-renderDT(pathway_res$react@result, filter = "top" ,extensions = 'Buttons', options=table_opts, server = FALSE)
           output$kegg_output_plot<-renderPlot(pathway_res$kegg_plot)
           output$kegg_output_df<-renderDT(pathway_res$kegg@result, filter = "top" ,extensions = 'Buttons', options=table_opts, server = FALSE)
           output$go_output_plot<-renderPlot(pathway_res$go_plot)
           output$go_output_df<-renderDT(pathway_res$go@result, filter = "top" ,extensions = 'Buttons', options=table_opts, server = FALSE)
         removeModal()}
#----Pathway enrichment module - end----
#----Brain region enrichment module with visualization (ABA - fgsea)----
         if ('ccuts' %in% input$EnrichmentPipeline){
           showModal(modalDialog("Running brain region (ABA-fgsea) enrichment with visualizations...", footer=NULL))
           print('Coldcuts region enrichment ...')
           output$current<-renderText('Brain region (ABA) enrichment with visualization...')
           print('Getting expression...')
           expr_df<-get_expression_df(sources=c('aba_ds_adult'))
           enr_results<-enrich_genes_fgsea(genes, expr_df)
           print('Preprocessing enrichment results...')
           print(enr_results)
           enr_results<-preprocess_enrichment_res_forcoldcuts(enr_results, brainonsources=c('ABA_dataset_adult','GTEx'))
           print('Getting the segmentation...')
           seg<-get_seg_default()
           
           print('Getting the assay for the left hemi...')
           seg<-get_var_distribution_brainwise(seg=seg, brainstats=enr_results$left_hemi_minp, assay_name='left_minp')
           print('Getting the assay for the right hemi...')
           seg<-get_var_distribution_brainwise(seg=seg, brainstats=enr_results$right_hemi_minp, assay_name='right_minp')
           print('Getting the assay for both hemi...')
           seg<-get_var_distribution_brainwise(seg=seg, brainstats=enr_results$all_minp, assay_name='all_minp')
           
           print('Getting the assay for the left hemi (medians)...')
           seg<-get_var_distribution_brainwise(seg=seg, brainstats=enr_results$left_hemi_med, assay_name='left_med')
           print('Getting the assay for the right hemi (medains)...')
           seg<-get_var_distribution_brainwise(seg=seg, brainstats=enr_results$right_hemi_med, assay_name='right_med')
           print('Getting the assay for both hemi (medians)...')
           seg<-get_var_distribution_brainwise(seg=seg, brainstats=enr_results$all_med, assay_name='all_med')
           
           
           print('Getting the assay for the left hemi (means)...')
           seg<-get_var_distribution_brainwise(seg=seg, brainstats=enr_results$left_hemi_mean, assay_name='left_mean')
           print('Getting the assay for the right hemi (means)...')
           seg<-get_var_distribution_brainwise(seg=seg, brainstats=enr_results$right_hemi_mean, assay_name='right_mean')
           print('Getting the assay for both hemi (means)...')
           seg<-get_var_distribution_brainwise(seg=seg, brainstats=enr_results$all_mean, assay_name='all_mean')
           
           
           left_nes_plot_minp<-seg_feature_plot(segmentation = seg,
                            assay = "left_minp",
                            feature = "NES",
                            smooth=FALSE,
                            labelsize=6,
                            remove_axes=TRUE)
           right_nes_plot_minp<-seg_feature_plot(segmentation = seg,
                                           assay = "right_minp",
                                           feature = "NES",
                                           smooth=FALSE,
                                           labelsize=6,
                                           remove_axes=TRUE)
           all_nes_plot_minp<-seg_feature_plot(segmentation = seg,
                                            assay = "all_minp",
                                            feature = "NES",
                                            smooth=FALSE,
                                            labelsize=6,
                                            remove_axes=TRUE)
           print('Minp - done.')
           
           left_nes_plot_med<-seg_feature_plot(segmentation = seg,
                                                assay = "left_med",
                                                feature = "NES",
                                                smooth=FALSE,
                                               labelsize=6,
                                               remove_axes=TRUE)
           right_nes_plot_med<-seg_feature_plot(segmentation = seg,
                                                 assay = "right_med",
                                                 feature = "NES",
                                                 smooth=FALSE,
                                                labelsize=6,
                                                remove_axes=TRUE)
           all_nes_plot_med<-seg_feature_plot(segmentation = seg,
                                               assay = "all_med",
                                               feature = "NES",
                                               smooth=FALSE,
                                              labelsize=6,
                                              remove_axes=TRUE)  
           print('Med - done.')
           
           left_nes_plot_mean<-seg_feature_plot(segmentation = seg,
                                               assay = "left_mean",
                                               feature = "NES",
                                               smooth=FALSE,
                                               labelsize=6,
                                               remove_axes=TRUE)
           right_nes_plot_mean<-seg_feature_plot(segmentation = seg,
                                                assay = "right_mean",
                                                feature = "NES",
                                                smooth=FALSE,
                                                labelsize=6,
                                                remove_axes=TRUE)
           all_nes_plot_mean<-seg_feature_plot(segmentation = seg,
                                              assay = "all_mean",
                                              feature = "NES",
                                              smooth=FALSE,
                                              labelsize=6,
                                              remove_axes=TRUE)
           print('Mean - done.')
           output$ccuts_ui <- renderUI({
             check1 <- 'ccuts' %in% input$EnrichmentPipeline
             if(length(check1)==0){check1 <- F}
             if(check1){
               fluidRow(style=border_style,tabBox(title = h1("Brain region enrichment (ABA-fgsea) with visualization results"),id= "ccuts_tabs", width = 12, #height = "420px",
                      tabPanel("Brain enrichment NES plots", 
                               h3('Left hemisphere (aggregated by min(padj)):'),
                               plotOutput("brainreg_output_plot_left",height = brainplot_hight),
                               h3('Right hemisphere (aggregated by min(padj)):'),
                               plotOutput("brainreg_output_plot_right",height = brainplot_hight),
                               h3('Both hemispheres (aggregated by min(padj)):'),
                               plotOutput("brainreg_output_plot_all",height = brainplot_hight),
                               
                               h3('Left hemisphere (aggregated by medians):'),
                               plotOutput("brainreg_output_plot_left_med",height = brainplot_hight),
                               h3('Right hemisphere (aggregated by medians):'),
                               plotOutput("brainreg_output_plot_right_med",height = brainplot_hight),
                               h3('Both hemispheres (aggregated by medians):'),
                               plotOutput("brainreg_output_plot_all_med",height = brainplot_hight),
                      
                               h3('Left hemisphere (aggregated by means):'),
                               plotOutput("brainreg_output_plot_left_mean",height = brainplot_hight),
                               h3('Right hemisphere (aggregated by means):'),
                               plotOutput("brainreg_output_plot_right_mean",height = brainplot_hight),
                               h3('Both hemispheres (aggregated by means):'),
                               plotOutput("brainreg_output_plot_all_mean",height = brainplot_hight)),                  
                      
                      tabPanel("Brain region enrichment (ABA-fgsea) table", DTOutput("brainreg_output_df"))
               ))}
           })
           print('Plotting...')
           output$brainreg_output_plot_left<-renderPlot(plot(left_nes_plot_minp))
           output$brainreg_output_plot_right<-renderPlot(plot(right_nes_plot_minp))
           output$brainreg_output_plot_all<-renderPlot(plot(all_nes_plot_minp))
           
           output$brainreg_output_plot_left_med<-renderPlot(plot(left_nes_plot_med))
           output$brainreg_output_plot_right_med<-renderPlot(plot(right_nes_plot_med))
           output$brainreg_output_plot_all_med<-renderPlot(plot(all_nes_plot_med))
           
           output$brainreg_output_plot_left_mean<-renderPlot(plot(left_nes_plot_mean))
           output$brainreg_output_plot_right_mean<-renderPlot(plot(right_nes_plot_mean))
           output$brainreg_output_plot_all_mean<-renderPlot(plot(all_nes_plot_mean))           
           
           
           
           output$brainreg_output_df<-renderDT(enr_results$full_annotated_enrichment_data %>% 
                                                 select(-pathway), filter = "top" ,extensions = 'Buttons', options=table_opts, server = FALSE)
           removeModal() }
#----Brain region enrichment module with visualization (ABA - fgsea) - end----
#----Brain region enrichment (ABA) - fgsea----
         if ('aba_fgsea' %in% input$EnrichmentPipeline){
           showModal(modalDialog("Running ABA fgsea enrichment...", footer=NULL))
           print('ABA fgsea...')
           output$current<-renderText('ABA fgsea enrichment...')
           expr_df<-get_expression_df(sources=c('aba_ds_adult'))
           enr_results_abafgsea<-enrich_genes_fgsea(genes, expr_df)%>% 
             select(-pathway)
           
           output$aba_fgsea_ui <- renderUI({
             check1 <- 'aba_fgsea' %in% input$EnrichmentPipeline
             if(length(check1)==0){check1 <- F}
             if(check1){
               fluidRow(style=border_style,tabBox(title = h1("ABA brain region enrichment results with fgsea"),id= "aba_fgsea_tabs", width = 12, #height = "420px",
                      tabPanel("ABA brain region enrichment table with fgsea", DTOutput("aba_fgsea_output_df"))
               ))
               
               
               
             }
             
           }
           )
           output$aba_fgsea_output_df<-renderDT(enr_results_abafgsea, filter = "top" ,extensions = 'Buttons', options=table_opts, server = FALSE)
           output$current<-renderText('ABA fgsea enrichment-done.')
         removeModal()}         
#----Brain region enrichment (ABA) - fgsea - end----
#----Brain region enrichment (ABA) - standard----
         if ('aba_std' %in% input$EnrichmentPipeline){
           showModal(modalDialog("Running ABA enrichment...", footer=NULL))
           print('ABA standard...')
           output$current<-renderText('ABA standard enrichment...')
           enr_results_abastd<-aba_enrich_genes(genes)
           enr_results_abastd<-enr_results_abastd$results %>% 
             arrange(-(n_significant), min_FWER)
           output$aba_std_ui <- renderUI({
             check1 <- 'aba_std' %in% input$EnrichmentPipeline
             if(length(check1)==0){check1 <- F}
             if(check1){
               fluidRow(style=border_style,tabBox(title = h1("ABA brain region enrichment results with ABAEnrichment"),id= "aba_std_tabs", width = 12, #height = "420px",
                      tabPanel("ABA brain region enrichment table with with ABAEnrichment", DTOutput("aba_std_output_df"))
               ))
             }
           }
           )
           output$aba_std_output_df<-renderDT(enr_results_abastd, filter = "top" ,extensions = 'Buttons', options=table_opts, server = FALSE)
           output$current<-renderText('ABA fgsea enrichment-done.')
         removeModal()}         
#----Brain region enrichment (ABA) - standard - end----
#----GTEx enrichment----
         if ('gtex' %in% input$EnrichmentPipeline){
           showModal(modalDialog("Running GTEx enrichment...", footer=NULL))
           print('GTEx fgsea...')
           output$current<-renderText('GTEx fgsea enrichment...')
           expr_df<-get_expression_df(sources=c('gtex'))
           enr_results_gt<-enrich_genes_fgsea(genes, expr_df) %>% 
             select(-pathway)
           
           output$gtex_ui <- renderUI({
             check1 <- 'gtex' %in% input$EnrichmentPipeline
             if(length(check1)==0){check1 <- F}
             if(check1){
               fluidRow(style = border_style, tabBox(title = h1("GTEx enrichment results with fgsea"),id= "gtex_tabs", width = 12, #height = "420px",
                      tabPanel("GTEx enrichment table with fgsea", DTOutput("gtex_df"))
               ))
               
               
               
             }
             
           }
           )
           output$gtex_df<-renderDT(enr_results_gt, filter = "top" ,extensions = 'Buttons', options=table_opts, server = FALSE)
           output$current<-renderText('GTEx fgsea enrichment-done.')
         removeModal()}            
#----GTEx enrichment - end----         
#----Cell type enrichment----
         if ('caba' %in% input$EnrichmentPipeline){
           showModal(modalDialog("Running Cell type (ABA) enrichment...", footer=NULL))
           print('Cell types fgsea...')
           output$current<-renderText('Cell types fgsea enrichment...')
           expr_df<-get_expression_df(sources=c('aba_cells'))
           enr_results_caba<-enrich_genes_fgsea(genes, expr_df) %>% 
             select(-pathway)
           
           output$caba_ui <- renderUI({
             check1 <- 'caba' %in% input$EnrichmentPipeline
             if(length(check1)==0){check1 <- F}
             if(check1){
               fluidRow(style=border_style, tabBox(title = h1("Cell type enrichment results with fgsea"),id= "caba_tabs", width = 12, #height = "420px",
                      tabPanel("Cell type enrichment table with fgsea", DTOutput("caba_df"))
               ))
               
               
               
             }
             
           }
           )
           output$caba_df<-renderDT(enr_results_caba, filter = "top" ,extensions = 'Buttons', options=table_opts, server = FALSE)
           output$current<-renderText('Cell types fgsea enrichment-done.')
         removeModal()}            
#----Cell type enrichment - end---- 
#----Compound selection----
        if ('drugs' %in% input$EnrichmentPipeline){
          showModal(modalDialog("Selecting compounds interacting with the genes...", footer=NULL))
           print('Compound selection...')
           output$current<-renderText('Compound selection procedure...')
           drugs_subset<-subset_drugs_bygenes(genes)
           
           
           output$caba_ui <- renderUI({
             check1 <- 'drugs' %in% input$EnrichmentPipeline
             if(length(check1)==0){check1 <- F}
             if(check1){
               fluidRow(style=border_style,tabBox(title = h1("Drug selection results with fgsea"),id= "drug_tabs", width = 12, #height = "420px",
                      tabPanel("Drug-target gene interaction table", DTOutput("drug_df")),
                      tabPanel("Drugs interacting with target genes", DTOutput("unique_drugs_df"))
               ))
               
               
               
             }
             
           }
           )
           output$drug_df<-renderDT(drugs_subset$df %>% 
                                      select(-interaction_present), filter = "top" ,extensions = 'Buttons', options=table_opts, server = FALSE)
           output$unique_drugs_df<-renderDT(tibble(drugs=drugs_subset$drugs), filter = "top" ,extensions = 'Buttons', options=table_opts, server = FALSE)
           output$current<-renderText('Compound selection-done.')
        removeModal()}         
#----Compound selection - end----
#----ATC enrichment----
      if ('atc' %in% input$EnrichmentPipeline){
        showModal(modalDialog("Running ATC target enrichment...", footer=NULL))
           print('ATC enrichment...')
           output$current<-renderText('ATC enrichment...')
           atc_enrichres<-enrich_atc(genes, showcats=showcats)
           
           if (length(atc_enrichres)>0){
             output$atc_ui <- renderUI({
               check1 <- 'atc' %in% input$EnrichmentPipeline
               if(length(check1)==0){check1 <- F}
               if(check1){
                 fluidRow(style=border_style,tabBox(title = h1("ATC categories enrichment results"),id= "atc_tabs", width = 12, #height = "420px",
                        tabPanel("ATC enrichment plot", plotOutput("atc_plot")),
                        tabPanel("ATC enrichment table", DTOutput("atc_df"))
                 ))
                 
                 
                 
               }
               
             }
             )
             #print(atc_enrichres)
             print(atc_enrichres)
             output$atc_plot<-renderPlot(atc_enrichres$plot)
             output$atc_df<-renderDT(atc_enrichres$df, filter = "top" ,extensions = 'Buttons', options=table_opts, server = FALSE)
          
             output$current<-renderText('ATC enrichment-done.')
           } else {
             print('Not enough ATC drug targets for enrichment analysis.')
             output$atc_ui<-renderUI({textOutput('atc_message')})
             output$atc_message<-renderText('Not enough ATC drug targets for enrichment analysis.')
             
           }
      removeModal()}         
#----ATC enrichment - end----
#----Small molecule enrichment----
         if ('sm' %in% input$EnrichmentPipeline){
           showModal(modalDialog("Getting small signaling molecule interactors...", footer=NULL))
           print('Small signaling molecule enrichment...')
           output$current<-renderText('Small signaling molecule enrichment...')
           sm_enrichres<-subset_metabolites_bygenes(genes, showcats=showcats)
           
           if (length(sm_enrichres$df$small_molecule_interactors)>0){
             output$sm_ui <- renderUI({
               check1 <- 'sm' %in% input$EnrichmentPipeline
               if(length(check1)==0){check1 <- F}
               if(check1){
                 fluidRow(style=border_style,tabBox(title = h1("Small signaling molecule enrichment results"),id= "sm_tabs", width = 12, #height = "420px",
                        tabPanel("Small signaling molecule enrichment plot", plotOutput("sm_plot")),
                        tabPanel("Small signaling molecule enrichment table", DTOutput("sm_df")),
                        tabPanel("Small signaling molecule interactions table", DTOutput("sm_inter_df")),
                        tabPanel("Small signaling molecule interactors", DTOutput("sm_comp_df"))
                 ))
                 
                 
                 
               }
               
             }
             )

             output$sm_plot<-renderPlot(sm_enrichres$plot)
             output$sm_df<-renderDT(sm_enrichres$enrichment_df, filter = "top" ,extensions = 'Buttons', options=table_opts, server = FALSE)
             output$sm_inter_df<-renderDT(tibble(small_molecule_interactors=sm_enrichres$df), filter = "top" ,extensions = 'Buttons', options=table_opts, server = FALSE)
             output$sm_comp_df<-renderDT(tibble(small_molecule_interactors=sm_enrichres$small_mol_ligands), filter = "top" ,extensions = 'Buttons', options=table_opts, server = FALSE)
             output$current<-renderText('Small signaling molecule enrichment-done.')
           } else {
             print('Not enough small signaling molecule interactors for enrichment analysis.')
             output$sm_ui<-renderUI({textOutput('sm_message')})
             output$sm_message<-renderText('Not enough small signaling molecule interactors for analysis.')
             
           }
         removeModal()}         
#----Small molecule enrichment - end----         
       }
  })
}

shinyApp(ui, server)


