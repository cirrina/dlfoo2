
# sample_id <- "KIRP_TP_3467"


#' Start Shiny app for plotting taxonomy/annotations on tSNE analyses of GEX, Me etc. Highlights individual TCGA samples. Uses pre-defined data sets loaded from disk.
#' @family tcga
#' @family shiny
#' @family PlotWrap
#' @param sample_id chracter string of TCGA sample or barcode ids
#' @param id_type if sample id or barcode. defaults to dlfoo style sample_id
#' @return a shiny app plot will open
#' @export
tsneShinyIndivSampleEvaluation <- function(my_samples, id_type="sample_id"){
   require(dlfoo2)
   require(shiny)
   require(dplyr)
   require(ggplot2)
   require(survival)

   #data("shape_subtypes", package="dlfoo2")
   #data("color_subtypes", package="dlfoo2")
   #data("color_annotations", package="dlfoo2")

   data("gdc_cn", package="dlfoo2data")
   data("tcga_surv_tab", package="dlfoo2data")
   data("pdata_utc_nt", package="dlfoo2data")
   data("pdata_tcgaPanCan", package="dlfoo2data")
   data("pan_gdc", package="dlfoo2data")
   data("utc_nt", package="dlfoo2data")


   tsne_gex_utc <- readRDS("~/PROJECTS/RCC_ccpRCC_2019/tSNE/tsne_UTCnt_Center_sd05.rds")
   # tsne_me_utc <- readRDS("~/RESOURCES/Methylation_450/Analyses/VMP_CorrGex/utc_me_1126set_VMP_CorrGEX/RCC_1305set_ECR_tSNE/RCC_1305set_ECR_tSNEtSneResults.rds")
   # tsne_mir_utc <- readRDS("~/PROJECTS/RCC_ccpRCC_2019/miRNA/tSNE/tsne_mir_utc_nt_sd05/tsne_mir_utc_nt_sd05tSneResults.rds")
   #tsne_rppa_utc <- readRDS("~/PROJECTS/RCC_ccpRCC_2019/Protein_RPPA/tSNE/tsne_utcnt/tsne_utcnttSneResults.rds")
   ## Me/mir pan instead of utc
   tsne_gex_pancan <- readRDS("~/PROJECTS/RCC_ccpRCC_2019/tSNE/tsne_PanCanGDC_Center_sd1_6211genes.rds")
   tsne_rppa_pancan <- readRDS("~/PROJECTS/RCC_ccpRCC_2019/Protein_RPPA/tSNE/tsne_rppa_PanCan/tsne_rppa_PanCantSneResults.rds")
   tsne_me_pancan <- readRDS("~/PROJECTS/RCC_ccpRCC_2019/Methylation/tSNE/RCC_9739set_ECR_tSNE_p50/RCC_9739set_ECR_tSNE_p50tSneResults.rds")
   tsne_mir_pancan <- readRDS("~/PROJECTS/RCC_ccpRCC_2019/miRNA/tSNE/tsne_mir_panCan_gdc_mir_sd05/tsne_mir_panCan_gdc_mir_sd05tSneResults.rds")

   tsne.res.list <- list(
      tsne_gex = tsne_gex_utc,
      # tsne_me = tsne_me_utc,
      # tsne_mir = tsne_mir_utc,
      tsne_me = tsne_me_pancan,
      tsne_mir = tsne_mir_pancan,
      # tsne_rppa = tsne_rppa_utc,
      tsne_rppa = tsne_rppa_pancan,
      tsne_gex_pancan = tsne_gex_pancan
      )

   pdata <- dlfoo2data::pdata_pan

   tsne_list <- lapply(tsne.res.list, function(x) do.call("rbind", lapply(x, dlfoo2::tsne2df)))
   #tsne_list <- lapply(tsne.res.list, function(x) lapply(x, dlfoo2::tsne2df))
   tsne_list <- lapply(tsne_list, function(x){
      x <- dplyr::left_join(x, pdata)
      return(x)
   })
   #dataset_choices <- names(tsne_list)
   annotation_choices <- colnames(pdata %>% dplyr::select(-contains("barcode"), -sample_id))
   default.annotation <- "taxonomy_published"


   survPlotWrap <- function(y){
                  z <- as.data.frame(tcga_surv_tab)[tcga_surv_tab$patient_barcode %in% y,]
                  z$fu_group <- "selected tumors"
                  dlfoo2::tcgaSurvPlot(z)
                  }


   sample_df <- pdata %>%
      mutate(sample_vec = as.character(taxonomy_published)) %>%
      mutate(sample_vec = if_else(sample_id %in% my_samples, "selected", sample_vec)) %>%
      dplyr::filter(!is.na(sample_vec)) %>%
      dplyr::select(sample_id, sample_vec)
   sample_groups <- unique(sample_df$sample_vec)

   if(length(my_samples)>1) centroid_selected <- apply(exprs(utc_nt)[,sample_df$sample_id[sample_df$sample_vec=="selected"]], 1, mean)
   if(length(my_samples)==1) centroid_selected <- exprs(utc_nt)[,sample_df$sample_id[sample_df$sample_vec=="selected"]]

   gene_deltas <- sapply(sample_groups[!sample_groups=="selected"], function(x){
      x_cent <- apply(exprs(utc_nt[,sample_df$sample_id[sample_df$sample_vec==x]]), 1, mean)
      return(centroid_selected - x_cent)
   })
   gene_deltas <- round(gene_deltas, 2)
   gene_df <- cbind(fData(pan_gdc), gene_deltas)

   # str(gene_table)


   #  Shiny App
## ::::::::::::::::::::::::::::::
   shinyApp(
   ## UI section
   ## ::::::::::::::::::::::::::::::
      ui <- fluidPage(
         tags$hr(),
         h1('tSNE plotter'),

         ## 1: select dataset and perplexes
         # fluidRow(
         #    column(3,
         #       selectInput(inputId = "dataToken", label = "Choose main data set:",
         #          selected = dataset_choices[1],
         #          choices = dataset_choices))),
         fluidRow(
            column(3,
               selectInput(inputId = "annotToken", label = "Choose main annotation:",
                  selected = default.annotation,
                  choices = annotation_choices))),


         # tags$hr(),
         fluidRow(
            # column(2,
            #    uiOutput("perplexity_master")),
            # column(3,
            #    uiOutput("annotation_selection")),
            column(2,
               sliderInput(inputId = "pointSize", label = "Select point size",
                  min=0.5, max=10, ticks = T, step = 0.5, value = 5)),
            column(2,
               sliderInput(inputId = "pointAlpha", label = "Select alpha",
                  min=0, max=1, ticks = T, step = 0.1, value = 0.7))
               ),

         # tags$hr(),
         # h2(textOutput(outputId = 'dataTokenText')),
         #
         # fluidPage(
         #       title = "tSNE Plot",
         #       plotOutput('tsne_plot_master', height=900, width=1200)
         #       ),
         #
         # tags$hr(),
         # fluidPage(
         #       title = "tSNE Plot selected annotation",
         #       plotOutput('tsne_plot_master_annot', height=900, width=1200, brush = "plot_brush_master")
         #       ),
         #
         tags$hr(),
         fluidRow(p(class = 'text-center', downloadButton(outputId = 'master_table_downloadbutton', label = 'Download selected data points'))),

         tags$hr(),
         fluidPage(
            title="Download",
            DT::dataTableOutput('sample_pdata_table')),


         ## GEX UTC
         tags$hr(),
         h2('TCGA UTC GEX Data'),
         fluidRow(
            column(4,
               uiOutput("perplexity_gex_utc")
            )),
         fluidRow(
            column(6,
               title = "tSNE Plot selected annotation",
               plotOutput('tsne_plot_annot_gex_utc', height=600, width=500)
               ),
            column(6,
               title = "tSNE Plot selected samples",
               plotOutput('tsne_plot_selection_gex_utc', height=600, width=500)
               ),
            column(6,
               title = "selected sample",
               plotOutput('tsne_plot_singleSelection_gex_utc', height=600, width=500)
               )
            ),

         ## Methylation
         tags$hr(),
         h2('TCGA UTC Me Data'),
         fluidRow(
            column(4,
               uiOutput("perplexity_me_utc")
            )),
         fluidRow(
            column(6,
               title = "tSNE Plot selected annotation",
               plotOutput('tsne_plot_annot_me_utc', height=600, width=500)
               ),
            column(6,
               title = "tSNE Plot selected samples",
               plotOutput('tsne_plot_selection_me_utc', height=600, width=500)
               ),
            column(6,
               title = "selected sample",
               plotOutput('tsne_plot_singleSelection_me_utc', height=600, width=500)
               )
            ),

         ## miRNA
         tags$hr(),
         h2('TCGA UTC mir Data'),
         fluidRow(
            column(4,
               uiOutput("perplexity_mir_utc")
            )),
         fluidRow(
            column(6,
               title = "tSNE Plot selected annotation",
               plotOutput('tsne_plot_annot_mir_utc', height=600, width=500)
               ),
            column(6,
               title = "tSNE Plot selected samples",
               plotOutput('tsne_plot_selection_mir_utc', height=600, width=500)
               ),
            column(6,
               title = "selected sample",
               plotOutput('tsne_plot_singleSelection_mir_utc', height=600, width=500)
               )
            ),

         ## RPPA
         tags$hr(),
         h2('TCGA UTC RPPA Data'),
         fluidRow(
            column(4,
               uiOutput("perplexity_rppa_utc")
            )),
         fluidRow(
            column(6,
               title = "tSNE Plot selected annotation",
               plotOutput('tsne_plot_annot_rppa_utc', height=600, width=500)
               ),
            column(6,
               title = "tSNE Plot selected samples",
               plotOutput('tsne_plot_selection_rppa_utc', height=600, width=500)
               ),
            column(6,
               title = "selected sample",
               plotOutput('tsne_plot_singleSelection_rppa_utc', height=600, width=500)
               )
            ),

         # PANCAN GEX DATA
         tags$hr(),
         h2('TCGA PanCan GEX Data'),
         fluidRow(
            column(4,
               uiOutput("perplexity_gex_pancan")
               )),
         tags$hr(),

         fluidPage(
               title = "tSNE Plot selected annotation",
               plotOutput('tsne_plot_annot_gex_pancan', height=1100, width=1200,  brush = "plot_brush_pancan")),
         fluidPage(
               title = "tSNE Plot selected samples",
               plotOutput('tsne_plot_selection_gex_pancan', height=1100, width=1200)),
         fluidPage(
            plotOutput('tsne_plot_singleSelection_gex_pancan', height=1100, width=1200
            )),


         fluidPage(
            title="PanCan Selection",
            DT::dataTableOutput('plot_pancan_table_selection')),

         ## Genome plots
         tags$hr(),
         fluidPage(
            uiOutput("SelectIndivSample")
         ),

         fluidPage(
              plotOutput('genome_plot')
           ),


         tags$hr(),
         fluidPage(
            title="Gene Selection",
            DT::dataTableOutput('gene_table')
            ),
         fluidPage(
              plotOutput('tcgaViolinPlot')
           )




         # ## Individual sample plots
         # tags$hr(),
         # h3('GEX urinary tract'),
         # fluidRow(
         #    column(4,
         #       title = h1('GEX UTC'),
         #       plotOutput('tsne_plot_gex_utc', height=600, width=500)
         #       ),
         #    column(4 ,
         #       title = h1('selected sample'),
         #       plotOutput('tsne_plot_singleSelection_gex_utc', height=600, width=500)
         #       )),
         # h3('Methylation urinary tract'),
         # fluidRow(
         #    column(4,
         #       title = h1('Me UTC'),
         #       plotOutput('tsne_plot_me_utc', height=600, width=500)
         #       ),
         #    column(4 ,
         #       title = h1('selected sample'),
         #       plotOutput('tsne_plot_singleSelection_me_utc', height=600, width=500)
         #       )),
         # h3('miRNA urinary tract'),
         # fluidRow(
         #    column(4,
         #       title = h1('mir UTC'),
         #       plotOutput('tsne_plot_mir_utc', height=600, width=500)
         #       ),
         #    column(4 ,
         #       title = h1('selected sample'),
         #       plotOutput('tsne_plot_singleSelection_mir_utc', height=600, width=500)
         #       )),
         # h3('protein RPPA urinary tract'),
         # fluidRow(
         #    column(4,
         #       title = h1('GEX UTC'),
         #       plotOutput('tsne_plot_rppa_utc', height=600, width=500)
         #       ),
         #    column(4 ,
         #       title = h1('selected sample'),
         #       plotOutput('tsne_plot_singleSelection_rppa_utc', height=600, width=500)
         #       )),
         # h3('GEX all TCGA samples'),
         # fluidPage(
         #    plotOutput('tsne_plot_gex_pancan', height=1100, width=1400
         #    )),
         # fluidPage(
         #    plotOutput('tsne_plot_singleSelection_gex_pancan', height=1100, width=1200
         #    ))

         ), # end ui fluid page


         # Define server logic
         server = function(input, output) {

            # Select main annotation

            output$sample_pdata_table <- DT::renderDataTable(
                  subset(pdata, sample_id %in% my_samples), server = FALSE, selection = 'none'
               )
            output$master_table_downloadbutton = downloadHandler(filename = 'pdata_selected_samples.csv', content = function(file) {
                     write.table(subset(pdata, sample_id %in% my_samples), file, sep=";",row.names = F)
                  }) # end download table



            observe({
               x = input$annotToken
               # observe data sets
   perpexity_choices_gex_utc = sort(unique(tsne_list[["tsne_gex"]]$perplexity))
               perpexity_choices_me_utc = sort(unique(tsne_list[["tsne_me"]]$perplexity))
               perpexity_choices_mir_utc = sort(unique(tsne_list[["tsne_mir"]]$perplexity))
               perpexity_choices_rppa_utc = sort(unique(tsne_list[["tsne_rppa"]]$perplexity))
               perpexity_choices_gex_pan = sort(unique(tsne_list[["tsne_gex_pancan"]]$perplexity))

               ## tsne GEX UTC
               output$perplexity_gex_utc <- renderUI({
                    selectInput(inputId = "perplex_gex_utc", label = "Available perplexities:",
                           choices=perpexity_choices_gex_utc, selected = perpexity_choices_gex_utc[3])
                  })
               output$tsne_plot_annot_gex_utc <- renderPlot(
                  dlfoo2::tsnePlot_special(
                     tsne_tab =  subset(tsne_list[["tsne_gex"]], perplexity==input$perplex_gex_utc),
                     annotation.column = input$annotToken,
                     my.size = input$pointSize*0.75, my.alpha = input$pointAlpha, showGuide=F))

               ## tsne Me UTC
               output$perplexity_me_utc <- renderUI({
                    selectInput(inputId = "perplex_me_utc", label = "Available perplexities:",
                           choices=perpexity_choices_me_utc, selected = perpexity_choices_me_utc[3])
                  })
               output$tsne_plot_annot_me_utc <- renderPlot(
                  dlfoo2::tsnePlot_special(
                     tsne_tab =  subset(tsne_list[["tsne_me"]], perplexity==input$perplex_me_utc),
                     annotation.column = input$annotToken,
                     my.size = input$pointSize*0.75, my.alpha = input$pointAlpha, showGuide=F))

               ## tsne mir UTC
               output$perplexity_mir_utc <- renderUI({
                    selectInput(inputId = "perplex_mir_utc", label = "Available perplexities:",
                           choices=perpexity_choices_mir_utc, selected = perpexity_choices_mir_utc[3])
                  })
               output$tsne_plot_annot_mir_utc <- renderPlot(
                  dlfoo2::tsnePlot_special(
                     tsne_tab =  subset(tsne_list[["tsne_mir"]], perplexity==input$perplex_mir_utc),
                     annotation.column = input$annotToken,
                     my.size = input$pointSize*0.75, my.alpha = input$pointAlpha, showGuide=F))

               ## tsne rppa UTC
               output$perplexity_rppa_utc <- renderUI({
                    selectInput(inputId = "perplex_rppa_utc", label = "Available perplexities:",
                           choices=perpexity_choices_rppa_utc, selected = perpexity_choices_rppa_utc[3])
                  })
               output$tsne_plot_annot_rppa_utc <- renderPlot(
                  dlfoo2::tsnePlot_special(
                     tsne_tab =  subset(tsne_list[["tsne_rppa"]], perplexity==input$perplex_rppa_utc),
                     annotation.column = input$annotToken,
                     my.size = input$pointSize*0.75, my.alpha = input$pointAlpha, showGuide=F))


               ## data pancan
               output$perplexity_gex_pancan <- renderUI({
                    selectInput(inputId = "perplex_gex_pancan", label = "Available perplexities:",
                           choices=perpexity_choices_gex_pan, selected = perpexity_choices_gex_pan[3])
                  })
               output$tsne_plot_annot_gex_pancan <- renderPlot(
                  dlfoo2::tsnePlot_special(
                     tsne_tab =  subset(tsne_list[["tsne_gex_pancan"]], perplexity==input$perplex_gex_pancan),
                     annotation.column = input$annotToken,
                     my.size = input$pointSize*0.75, my.alpha = input$pointAlpha, showGuide=F))

               output$plot_pancan_table_selection <- DT::renderDataTable(
                  brushedPoints(subset(tsne_list[["tsne_gex_pancan"]], perplexity==input$perplex_gex_pancan), input$plot_brush_pancan), server = FALSE, selection = 'none'
               )


               # MARK INDIVIDUAL SAMPLES IN TSNE
               output$tsne_plot_selection_gex_utc <- renderPlot(
                  dlfoo2::tsnePlot_sample(
                     tsne_tab = subset(tsne_list[["tsne_gex"]], perplexity==input$perplex_gex_utc),
                     my_samples = my_samples,
                     shape.val = c(4,21), col.val=c("tomato4","gray"), fill.val=c("tomato4","gray"),  size.val=c(2, 0.25), stroke.val=c(1,0.5), alpha.val=c(1,0.75)))
               output$tsne_plot_selection_me_utc <- renderPlot(
                  dlfoo2::tsnePlot_sample(
                     tsne_tab = subset(tsne_list[["tsne_me"]], perplexity==input$perplex_me_utc),
                     my_samples = my_samples,
                     shape.val = c(4,21), col.val=c("tomato4","gray"), fill.val=c("tomato4","gray"),  size.val=c(2, 0.25), stroke.val=c(1,0.5), alpha.val=c(1,0.75)))
               output$tsne_plot_selection_mir_utc <- renderPlot(
                  dlfoo2::tsnePlot_sample(
                     tsne_tab = subset(tsne_list[["tsne_mir"]], perplexity==input$perplex_mir_utc),
                     my_samples = my_samples,
                     shape.val = c(4,21), col.val=c("tomato4","gray"), fill.val=c("tomato4","gray"),  size.val=c(2, 0.25), stroke.val=c(1,0.5), alpha.val=c(1,0.75)))
               output$tsne_plot_selection_rppa_utc <- renderPlot(
                  dlfoo2::tsnePlot_sample(
                     tsne_tab = subset(tsne_list[["tsne_rppa"]], perplexity==input$perplex_rppa_utc),
                     my_samples = my_samples,
                     shape.val = c(4,21), col.val=c("tomato4","gray"), fill.val=c("tomato4","gray"),  size.val=c(2, 0.25), stroke.val=c(1,0.5), alpha.val=c(1,0.75)))
               output$tsne_plot_selection_gex_pancan <- renderPlot(
                  dlfoo2::tsnePlot_sample(
                     tsne_tab = subset(tsne_list[["tsne_gex_pancan"]], perplexity==input$perplex_gex_pancan),
                     my_samples = my_samples,
                     shape.val = c(4,21), col.val=c("tomato4","gray"), fill.val=c("tomato4","gray"),  size.val=c(2, 0.25), stroke.val=c(1,0.5), alpha.val=c(1,0.75)))


               ## Select sample for INDIVIDUAL PLOTS
               output$SelectIndivSample <- renderUI({
                  selectInput(inputId = "indivSample", label = "Genome plot",
                        choices=my_samples, selected = my_samples[1])
               })
            }) # end observe x




         observe({
            x = input$indivSample
            if(length(x)){
               output$selectedSample <- renderText(x)
               ## Render CN plot using dlfoo2::segmentPlotter
               if(x %in% names(gdc_cn)){
                  output$genome_plot <- renderPlot(dlfoo2::segmentPlotter(gdc_cn[[x]], sample_name = x))
                  }else{
                  output$genome_plot <- renderPlot(dlfoo2::ggPlotEmpty("Copy number profile not available"))
                  } # end CN plot

               # Plot individual selected sample in all three data sets

               output$tsne_plot_singleSelection_gex_utc <- renderPlot(
                  dlfoo2::tsnePlot_sample(
                     tsne_tab = subset(tsne_list[["tsne_gex"]], perplexity==input$perplex_gex_utc),
                     my_samples = x,
                     shape.val = c(4,21), col.val=c("tomato4","gray"), fill.val=c("tomato4","gray"),  size.val=c(2, 0.25), stroke.val=c(1,0.5), alpha.val=c(1,0.75)))
               output$tsne_plot_gex_utc <- renderPlot(
                  dlfoo2::tsnePlot_special(
                     tsne_tab = subset(tsne_list[["tsne_gex"]], perplexity==input$perplex_gex_utc),
                     annotation.column = default.annotation,
                     my.size = input$pointSize, my.alpha = input$pointAlpha, showGuide=F)
                  )


               output$tsne_plot_singleSelection_me_utc <- renderPlot(
                  dlfoo2::tsnePlot_sample(
                     tsne_tab = subset(tsne_list[["tsne_me"]], perplexity==input$perplex_me_utc),
                     my_samples = x,
                     shape.val = c(4,21), col.val=c("tomato4","gray"), fill.val=c("tomato4","gray"),  size.val=c(2, 0.25), stroke.val=c(1,0.5), alpha.val=c(1,0.75)))
               output$tsne_plot_me_utc <- renderPlot(
                  dlfoo2::tsnePlot_special(
                     tsne_tab = subset(tsne_list[["tsne_me"]], perplexity==input$perplex_me_utc),
                     annotation.column = default.annotation,
                     my.size = input$pointSize, my.alpha = input$pointAlpha, showGuide=F)
                  )

               output$tsne_plot_singleSelection_mir_utc <- renderPlot(
                  dlfoo2::tsnePlot_sample(
                     tsne_tab = subset(tsne_list[["tsne_mir"]], perplexity==input$perplex_mir_utc),
                     my_samples = x,
                     shape.val = c(4,21), col.val=c("tomato4","gray"), fill.val=c("tomato4","gray"),  size.val=c(2, 0.25), stroke.val=c(1,0.5), alpha.val=c(1,0.75)))
               output$tsne_plot_mir_utc <- renderPlot(
                  dlfoo2::tsnePlot_special(
                     tsne_tab = subset(tsne_list[["tsne_mir"]], perplexity==input$perplex_mir_utc),
                     annotation.column = default.annotation,
                     my.size = input$pointSize, my.alpha = input$pointAlpha, showGuide=F)
                  )

               output$tsne_plot_singleSelection_rppa_utc <- renderPlot(
                  dlfoo2::tsnePlot_sample(
                     tsne_tab = subset(tsne_list[["tsne_rppa"]], perplexity==input$perplex_rppa_utc),
                     my_samples = x,
                     shape.val = c(4,21), col.val=c("tomato4","gray"), fill.val=c("tomato4","gray"),  size.val=c(2, 0.25), stroke.val=c(1,0.5), alpha.val=c(1,0.75)))
               output$tsne_plot_rppa_utc <- renderPlot(
                  dlfoo2::tsnePlot_special(
                     tsne_tab = subset(tsne_list[["tsne_rppa"]], perplexity==input$perplex_rppa_utc),
                     annotation.column = default.annotation,
                     my.size = input$pointSize, my.alpha = input$pointAlpha, showGuide=F)
                  )

            output$tsne_plot_singleSelection_gex_pancan <- renderPlot(
                  dlfoo2::tsnePlot_sample(
                     tsne_tab = subset(tsne_list[["tsne_gex_pancan"]], perplexity==input$perplex_gex_pancan),
                     my_samples = x,
                     shape.val = c(4,21), col.val=c("tomato4","gray"), fill.val=c("tomato4","gray"),  size.val=c(2, 0.25), stroke.val=c(1,0.5), alpha.val=c(1,0.75)))
               output$tsne_plot_gex_pancan <- renderPlot(
                  dlfoo2::tsnePlot_special(
                     tsne_tab = subset(tsne_list[["tsne_gex_pancan"]], perplexity==input$perplex_gex_pancan),
                     annotation.column = default.annotation,
                     my.size = input$pointSize, my.alpha = input$pointAlpha)
                  )
               }# end if observe x input$indivSample
         })

      # Gene Selection table
      output$gene_table <- DT::renderDataTable(
                  gene_df, server = TRUE, selection = 'single')

      observe({
         row_i = input$gene_table_rows_selected
         if(length(row_i)){
            output$tcgaViolinPlot <- renderPlot(
               dlfoo2::geneplotPanTCGA_violin(ensg = gene_df[row_i,]$ENSG,
                  highlight.samples = my_samples, individual.sample=input$indivSample, default.taxonomy="taxonomy_published")
               )}
      })


      ## Clinical data (SurvPlot for selected samples)
         # observe({
         #    x = brushedPoints(gene_table, input$plot_brush)$patient_barcode
         #    if(length(xx)){
         #       output$survPlot <- renderPlot(survPlotWrap(xx)
         #
         #          )
         #
         #       }
         #       })


      }
   ) # end shiny app
} # end function shiny tsne





