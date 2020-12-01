


   # fdata <- readRDS("~/PROJECTS/BRCA_SR_KP/Data/fdata.rds")
   # pdata <- readRDS( "~/PROJECTS/BRCA_SR_KP/Data/pdata.rds")
   # pdata <- pdata[,c("sample", "seqno","time","cells","annot")]
   # rownames(pdata) <- pdata$sample
   # pdata_color <- read.delim("~/PROJECTS/BRCA_SR_KP/Data/pdata_color_key.txt", as.is = T)
   # color_key <- pdata_color$color
   # names(color_key) <- pdata_color$name
   #
   # gene_rcount <- readRDS( "~/PROJECTS/BRCA_SR_KP/Data/gene_rcount.rds")
   # rcount_mat <- as.matrix(gene_rcount[,pdata$sample])
   # rownames(rcount_mat) <- gene_rcount$gene_id
   # str(rcount_mat)
   #
   # sig_list <- readRDS("~/PROJECTS/BRCA_SR_KP/DESeq2_analysis/Significant_genes.rds")
   # names(sig_list)
   # vsd <- readRDS(file="~/PROJECTS/BRCA_SR_KP/DESeq2_analysis/DeSeq2_vst.rds")
   # vsd_mat <- assay(vsd)
   # vsd_mat <- vsd_mat[,pdata$sample]
   # mat <- vsd_mat
   # str(fdata)
   # identical(fdata$gene_id, rownames(vsd_mat))
   # identical(fdata$gene_id, rownames(rcount_mat))
   # rownames(fdata) <- fdata$gene_id
   #
   # ensg <- "ENSG00000188404"
   # symbol <- fdata[ensg,]$hgnc_symbol

   # boxplot_fdataRowSelection_SHINY(mat, fdata, pdata, color.key=color_key, fdata.altName.column = "hgnc_symbol", default.pdata.column="time")
   # eset <- dlfoo2data::pan_gdc
   #   source("~/PROJECTS/BRCA_SR_KP/R_SCRIPTS/R_CAF_Soruce.R")
   #   fdata_shiny <- fdata_res %>% dplyr::select(gene_id, hgnc_symbol, contains("stat"), everything())
   #   identical(rownames(vsd_mat), rownames(fdata_shiny) )
   #mat=vsd_mat, fdata=fdata_shiny, pdata, color.key=color_key, fdata.altName.column = "hgnc_symbol", default.pdata.column="annot")


#' generic shiny boxplots given an expression set (defaults to pan_gdc) or all three of fdata, pdata and a matrix
#' @description Shiny app for plotting annotations (pdata) for genes expression selectable from rows in a matrix (fdata)
#' @family shiny
#' @family boxplot
#' @family ggplot
#' @family violin
#' @alias boxplot_fdataRowSelection_SHINY
#' @param mat gene identifier
#' @param fdata
#' @param pdata sample data frame with sample annotations
#' @param color.key named character vector with sample bxp colors
#' @param fdata.altName.column alternative name column, e.g. hgnc_symbol
#' @param default.pdata.column name of preferred pdata annotations column (not needed)
#' @param ... additional parameters
#' @return a shiny app plot will open
#' @export
boxplot_SHINY_annotation <- function(
   eset = NULL,
   mat=NULL, fdata=NULL, pdata=NULL,
   color.key = dlfoo2::color_subtypes,
   fdata.altName.column=NULL,
   default.pdata.column=NULL
   ){

   require(Biobase)
   require(dlfoo2)
   require(shiny)
   require(dplyr)
   require(ggplot2)
   require(survival)
   require(umap)

   if(!is.null(mat)){
      message("... expression matrix defined. setting eset to null")
      eset = NULL
   }

   if(!is.null(eset)){
      if(class(eset)!="ExpressionSet") stop ("eset must be expression set")
      mat = exprs(eset)
      fdata = fData(eset)
      rownames(mat) <- featureNames(eset)
      rownames(fdata) <- featureNames(eset)
      pdata = pData(eset)
   }


   ## ::: BXP FUNCTION
   ## ----------------
   bxpFoo = function(
      x, mat, fdata, pdata, pdata.column=NULL,
      altName=NA, y.range.add=1, color.key){

      if(length(x)!=1){
         message(x, " only one id is allowed")
         return(NULL)
      }
      gene.i = match(x, rownames(fdata))
      my_df = pdata
      my_df$gex <- mat[gene.i, ]
      my_df$annot <- my_df[,pdata.column]

      range_y <- range(my_df$gex)

      g_theme <- ggplot(my_df) +
         theme(panel.background = element_rect(fill = 'gray95', colour = 'antiquewhite3', size = 0.5)) +
         theme(panel.grid.major =  element_line(colour = 'lightsteelblue2'), panel.grid.minor =  element_line(colour = 'mistyrose3')) +
         theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
         theme(axis.title.y = element_text(size=12)) +
         theme(axis.title.x = element_text(size=12)) +
         theme(axis.text = element_text(size=9)) +
         theme(title = element_text(size=9)) +
         guides(fill=FALSE)+
         theme(strip.text = element_text(size=12)) +
         theme(title = element_text(size=12)) +
         ggtitle(paste(altName, x, sep=" :  "))

      g <- g_theme +
         aes(x=annot, y=gex, fill=annot) +
         labs(y=paste0("expr")) +
         geom_violin(alpha=0.5, width=1, trim = TRUE, scale = "width", adjust = 0.5) +
         geom_boxplot(width=0.2, outlier.colour="grey75", notch = FALSE, notchwidth = .75, alpha = 0.75, colour = "grey50", outlier.size=1) +
         scale_fill_manual(values = color.key) +
         coord_cartesian(ylim =c(range_y[1]-y.range.add, range_y[2]+y.range.add))
      return(g)
   } # end internl bxp function



   ## ::: PREPARE DATA
   ## -----------
   if(!identical(rownames(mat), rownames(fdata))) stop("rownames for matrix and fdata must be idenitcal")
   if(!identical(colnames(mat), rownames(pdata))) stop("colnames for matrix and rownames pdata must be idenitcal")

   default.pdata.column <- if_else(is.null(default.pdata.column), false = default.pdata.column, true = colnames(pdata)[1])
   if(!default.pdata.column %in% colnames(pdata)) stop("pdata object must contain column specified by default.annotation")
   pdata_choices <- colnames(pdata)



   #  ::: Shiny App
   ## -----------------
   shinyApp(

   ## UI section
   ## ................
      ui <- fluidPage(
         tags$hr(),

                  fluidRow(
            # column(2,
            #    uiOutput("perplexity_master")),
            column(3, uiOutput("annotation_selection"))
            #column(2, sliderInput(inputId = "pointSize", label = "Select point size", min=0.5, max=10, ticks = T, step = 0.5, value = 5)),
            #column(2, sliderInput(inputId = "pointAlpha", label = "Select alpha", min=0, max=1, ticks = T, step = 0.1, value = 0.7))
         ),

         #
         # The PAN TCGA BOXPLOT
         tags$hr(),
         fluidPage(
            column(10, plotOutput('bxp'))
         ),

         h1('choose gene'),
         # Show the main table used to select ENSG_1 gene id
         # -- fdataTable_main (dataTableOutput)
         # -- ensg_1 (output)
         fluidRow(column(8, DT::dataTableOutput('fdataTable'))),
         tags$hr(),
         h1(textOutput('SYMBOL_1')),
         h3(textOutput('ENSG_1'))





      ),


      ## Define server logic
      ## ::::::::
      server <- function(input, output){

         output$annotation_selection <- renderUI({
            selectInput(inputId = "annot", label = "Available annotations:",
            choices=pdata_choices, selected = default.pdata.column)
            })


         # fdata data table used to select gene
         output$fdataTable = DT::renderDataTable(fdata, server = TRUE, selection = 'single')
         observe({
            u = input$fdataTable_rows_selected
            # if(length(u)){
            #    output$geneID = renderText({
            #       rownames(fdata)[u]
            #    })
            #    output$altName= renderText({
            #       fdata[u,fdata.altName.column]
            #    })
               output$bxp = renderPlot({
                  bxpFoo(x = rownames(fdata)[u],  mat = mat, fdata = fdata, pdata = pdata, pdata.column = input$annot, altName=fdata[u,fdata.altName.column], y.range.add=1, color.key=color.key)
               })
            }) # end observe
      } # end server
      ) # end shiny app
   }  # end shiny function









   # tsne_res <- readRDS("~/PROJECTS/RCC_ccpRCC/miRNA/tSNE/tsne_mir_utc_nt/tsne_mir_utc_nttSneResults.rds")
   #    tsne_res <- readRDS("~/RESOURCES/Methylation_450/Analyses/VMP_CorrGex/rcc_me450_654set_tumors_ECR/RCC_790set_ECR_tSNE/RCC_790set_ECR_tSNEtSneResults.rds")
   # str(tSNEres)
   # tsne_res <- readRDS("~/PROJECTS/RCC_ccpRCC/tSNE/tsne_RCC_BLCA_Center_sd05.rds")


   #   results.object <- readRDS(file = "~/PROJECTS/RCC_ccpRCC_2019/UMAP/mRNA/umap_RCC_NT_Center_sd1_n15_pearson.rds")

#' scatterPoint_annotationPlot_SHINY
#' @description Shiny app for plotting taxonomy/annotations on generic 2D scatter plot results objects, as PCA umap, tsne etc.
#' @family shiny
#' @family tsne
#' @family umap
#' @family PCA
#' @family PlotWrap
#' @param results.object results object. list (tSNEres) from run.tsne-analysis, or umap (umap object), or PCA (prcomp object)
#' @param pdata
#' @param ... additional parameters as used in the scatterPoint_annotationPlot function
#' @return a shiny app plot will open
#' @export
scatterPoint_annotationPlot_SHINY <- function(
   results.object, pdata=NULL, default.annotation="tax_simp", ...
   ){

   require(dlfoo2)
   require(shiny)
   require(dplyr)
   require(ggplot2)
   require(survival)
   require(umap)


   if(is.null(pdata)) pdata <- dlfoo2::pdata_panCanFull
   if(!"sample_id" %in% colnames(pdata)) stop("pdata object must contain 'sample_id' column")
   if(!default.annotation %in% colnames(pdata)) stop("pdata object must contain column specified by default.annotation")


    if(class(results.object)=="umap"){
      message("Found UMAP results")
      res_df <- dlfoo2::umap2df(results.object)
      plot_type <- "umap"
   }

   if(class(results.object)=="prcomp"){
      message("Found PCA results")
      plot_type <- "pca"
      res_df <- pca2df(results.object)
   }

   # if(class(results.object)=="list" & !is.null(results.object$Y) ){
   #    message("Found list with Y matrix slot ... ")
   #    message("... assuming x is a tsne reults list. Begin plotting")
   #    plot_type <- "tsne"
   #    res_df <- tsne2df(results.object)
   # }
   if(class(results.object)=="data.frame"){
      message("Found matrix ... ")
      message("... will plot the first 2 columns as dim1 and dim2")
      plot_type <- "df"
      res_df <- results.object
    }

   if(is.null(plot_type)) stop("results.object does not seem to be a prcomp, umap or tsne reusults list")


   ## plot_df: Join with pdata - check columns/ids
   u <- match(res_df$sample_id, pdata$sample_id)
   uu <- !is.na(u)
   if(all(!uu)) stop("cannot match sample_id for results.object (rownames) with sample_id column in pdata. if data frame you must have a sample_id column")

   message("... found ",length(which(uu==T))," matching saple_id in pdata out of ", length(which(uu==F)))
   plot_df <- dplyr::left_join(res_df, pdata[u[uu],])
   rownames(plot_df) <- plot_df$sample_id

   ## Set annotation
   # plot_df$annotation <- plot_df[, default.annotation]
   # plot_df$annotation[plot_df$annotation==""] <- NA

   annotation_choices <- colnames(pdata %>% dplyr::select(-contains("barcode"), -sample_id))
   #default.annotation <- "tax_simp"

   # survPlotWrap <- function(y){
   #                z <- as.data.frame(tcga_surv_tab)[tcga_surv_tab$patient_barcode %in% y,]
   #                z$fu_group <- "selected tumors"
   #                dlfoo2::tcgaSurvPlot(z)
   #                }




#  Shiny App
## ::::::::::::::::::::::::::::::
   shinyApp(

   ## UI section
   ## ::::::::::::::::::::::::::::::
      ui <- fluidPage(
         tags$hr(),

         h1('Annotation plotter'),

         fluidRow(
            # column(2,
            #    uiOutput("perplexity_master")),
            column(3, uiOutput("annotation_selection")),
            column(2, sliderInput(inputId = "pointSize", label = "Select point size", min=0.5, max=10, ticks = T, step = 0.5, value = 5)),
            column(2, sliderInput(inputId = "pointAlpha", label = "Select alpha", min=0, max=1, ticks = T, step = 0.1, value = 0.7))
         ),

         tags$hr(),

         fluidPage(title = "scatter Plot selected annotation",
               plotOutput('scatter_plot_annot', height=1200, width=1200, brush = "plot_brush")
         ),

         tags$hr(),

         fluidRow(p(class = 'text-center', downloadButton(outputId = 'table_downloadbutton', label = 'Download selected data points'))
            ),

         tags$hr(),

         fluidPage(title="Download", DT::dataTableOutput('plot_table_selection')),


         ## Genome plots
         tags$hr(),

         fluidPage(uiOutput("SelectGenomeSample")),

         fluidPage(plotOutput('genome_plot'))
      ),


      ## Define server logic
      ## ::::::::
      server <- function(input, output){

            output$annotation_selection <- renderUI({
                    selectInput(inputId = "annot", label = "Available annotations:",
                           choices=annotation_choices, selected = default.annotation)
            })

            output$scatter_plot_annot <- renderPlot(
               scatterPoint_annotationPlot(plot_df, pdata = pdata, annotation.column = input$annot, my.size=input$pointSize, my.alpha =input$pointAlpha , ...)
            )


            output$plot_table_selection <- DT::renderDataTable(
                  brushedPoints(plot_df, input$plot_brush)
            )

            # download the filtered data
            output$table_downloadbutton = downloadHandler(filename = 'selected_points.csv', content = function(file) {
                     write.table(brushedPoints(plot_df, input$plot_brush, file, sep=";",row.names = F))
            }) # end download table


         # Obseve selected samples
            observe({
               x = brushedPoints(plot_df,  input$plot_brush)$sample_id
               if(length(x)){
                  ## CN DATA - Select sample for CN plot (dropdown list)
                  output$SelectGenomeSample <- renderUI({
                  selectInput(inputId = "genomeSample", label = "Genome plot",
                        choices=x, selected = x[1])
                  })
               } # end if observe x brushed
            })

         observe({
            x = input$genomeSample
            if(length(x)){
               output$selectedSample <- renderText(x)
               ## Render CN plot using dlfoo2::segmentPlotter
               if(x %in% names(gdc_cn)){
                  output$genome_plot <- renderPlot(dlfoo2::segmentPlotter(dlfoo2data::gdc_cn[[x]], sample_name = x))
                  }else{
                  output$genome_plot <- renderPlot(dlfoo2::ggPlotEmpty("Copy number profile not available"))
                  } # end CN plot
            }# end if observe x input$genomeSample
         })
})}



# my_df <- pdata %>% mutate(beta=x.me[x,pdata$sample_id,"beta"])


   # fdata <- readRDS("~/PROJECTS/BRCA_SR_KP/Data/fdata.rds")
   # pdata <- readRDS( "~/PROJECTS/BRCA_SR_KP/Data/pdata.rds")
   # pdata <- pdata[,c("sample", "seqno","time","cells","annot")]
   # rownames(pdata) <- pdata$sample
   # pdata_color <- read.delim("~/PROJECTS/BRCA_SR_KP/Data/pdata_color_key.txt", as.is = T)
   # color_key <- pdata_color$color
   # names(color_key) <- pdata_color$name
   #
   # gene_rcount <- readRDS( "~/PROJECTS/BRCA_SR_KP/Data/gene_rcount.rds")
   # rcount_mat <- as.matrix(gene_rcount[,pdata$sample])
   # rownames(rcount_mat) <- gene_rcount$gene_id
   # str(rcount_mat)
   #
   # sig_list <- readRDS("~/PROJECTS/BRCA_SR_KP/DESeq2_analysis/Significant_genes.rds")
   # names(sig_list)
   # vsd <- readRDS(file="~/PROJECTS/BRCA_SR_KP/DESeq2_analysis/DeSeq2_vst.rds")
   # vsd_mat <- assay(vsd)
   # vsd_mat <- vsd_mat[,pdata$sample]
   # mat <- vsd_mat
   # str(fdata)
   # identical(fdata$gene_id, rownames(vsd_mat))
   # identical(fdata$gene_id, rownames(rcount_mat))
   # rownames(fdata) <- fdata$gene_id
   #
   # ensg <- "ENSG00000188404"
   # symbol <- fdata[ensg,]$hgnc_symbol

   # boxplot_fdataRowSelection_SHINY(mat, fdata, pdata, color.key=color_key, fdata.altName.column = "hgnc_symbol", default.pdata.column="time")
   # eset <- dlfoo2data::pan_gdc
   #   source("~/PROJECTS/BRCA_SR_KP/R_SCRIPTS/R_CAF_Soruce.R")
   #   fdata_shiny <- fdata_res %>% dplyr::select(gene_id, hgnc_symbol, contains("stat"), everything())
   #   identical(rownames(vsd_mat), rownames(fdata_shiny) )
   #mat=vsd_mat, fdata=fdata_shiny, pdata, color.key=color_key, fdata.altName.column = "hgnc_symbol", default.pdata.column="annot")


    ##GEX
      # pdata <- data.frame(dlfoo2::pdata_panCanFull %>%  dplyr::filter(!duplicated(sample_id)) %>% dplyr::filter(sample_id %in% sampleNames(dlfoo2data::rcc_nt)) %>% dplyr::filter(grepl("ccRCC|pRCC|chRCC|CIMP|ccpRCC|mesRCC|chONC|pONC", tax_simp)))
      # rownames(pdata) <- pdata$sample_id
      #
      # my_mat <- exprs(dlfoo2data::rcc_nt[,pdata$sample_id])
      # str(my_mat)
      # my_fdata <- fData(dlfoo2data::rcc_nt)
      # identical(rownames(my_fdata), rownames(my_mat))
      # identical(rownames(pdata), colnames(my_mat))
#       #
# shinyPlot_annotation(
#        mat= my_mat,
#          fdata=my_fdata,
#          pdata = pdata,
#          fdata.altName.column = "SYMBOL",
#          default.pdata.column = "tax_simp2",
#          color.key = color_subtypes_gex)
#



#' shinyPlot_annotation of two variables. Matrix provides Dim1 whereas dataframe matched with matrix colnames provides Dim2 and annotations..
#' @description Shiny app for plotting taxonomy/annotations on 2 generic variables, e.g. methylation vs purity
#' @description Matrix with matched fdata and pdata. matrix/fdata supplies Dim1, pdata dim2 and annotations
#' @description variation of the Shiny app for plotting taxonomy/annotations on generic 2D scatter plot results objects, as PCA umap, tsne etc.
#' @family shiny
#' @family scatterplot
#' @family scatterpoint
#' @param mat gene identifier. should
#' @param fdata
#' @param pdata sample data frame with sample annotations
#' @param color.key named character vector with sample bxp colors
#' @param fdata.altName.column alternative name column, e.g. hgnc_symbol
#' @param default.pdata.column name of preferred pdata annotations column (not needed)
#' @param ... additional parameters
#' @return a shiny app plot will open
#' @export
shinyPlot_annotation <- function(
   eset = NULL,
   mat=NULL, fdata=NULL, pdata=NULL,
   color.key = dlfoo2::color_subtypes,
   fdata.altName.column=NULL,
   default.xaxis.column=NULL,
   default.pdata.column=NULL,
   default.annot.column=NULL,
   x.lim=NULL, y.lim=NULL
   ){

   require(Biobase)
   require(dlfoo2)
   require(shiny)
   require(dplyr)
   require(ggplot2)
   require(survival)
   require(umap)

   if(!is.null(mat)){
      message("... expression matrix defined. setting eset to null")
      eset = NULL
   }

   if(!is.null(eset)){
      if(class(eset)!="ExpressionSet") stop ("eset must be expression set")
      mat = exprs(eset)
      fdata = fData(eset)
      rownames(mat) <- featureNames(eset)
      rownames(fdata) <- featureNames(eset)
      pdata = pData(eset)
   }



   ## ::: BXP FUNCTION
   ## ----------------
   plotFoo = function(
      x, mat, fdata, pdata, pdata.column=NULL, annot.column=NULL,
      altName=NA, y.range.add=1, color.key, xlim=x.lim, ylim=y.lim){

      if(length(x)!=1){
         message(x, " only one id is allowed")
         return(NULL)
      }
      gene.i = match(x, rownames(fdata))
      my_df = pdata
      my_df$dim1 <- mat[gene.i, ]
      my_df$dim2 <- my_df[,pdata.column]
      my_df$annot <- my_df[,annot.column]

      range_y <- range(my_df$dim1)

      g_theme <- ggplot(my_df) +
         theme(panel.background = element_rect(fill = 'gray95', colour = 'antiquewhite3', size = 0.5)) +
         theme(panel.grid.major =  element_line(colour = 'lightsteelblue2'), panel.grid.minor =  element_line(colour = 'mistyrose3')) +
         theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
         theme(axis.title.y = element_text(size=12)) +
         theme(axis.title.x = element_text(size=12)) +
         theme(axis.text = element_text(size=9)) +
         theme(title = element_text(size=9)) +
         guides(fill=FALSE)+
         theme(strip.text = element_text(size=12)) +
         theme(title = element_text(size=12)) +
         ggtitle(paste(altName, x, sep=" :  "))

      if(!is.numeric(my_df$dim2)){
      g <- g_theme +
         aes(x=dim2, y=dim1, fill=dim2) +
         #labs(y=paste0("expr")) +
         geom_violin(alpha=0.5, width=1, trim = TRUE, scale = "width", adjust = 0.5) +
         geom_boxplot(width=0.2, outlier.colour="grey75", notch = FALSE, notchwidth = .75, alpha = 0.75, colour = "grey50", outlier.size=1) +
         scale_fill_manual(values = color.key) +
         coord_cartesian(ylim =c(range_y[1]-y.range.add, range_y[2]+y.range.add))
      }

      if(is.numeric(my_df$dim2)){
         g <- g_theme +
         aes(x=dim2, y=dim1) +
         # labs(y=paste0("expr")) +
         geom_point(aes(color=annot), stroke=1, size=1, alpha=0.75) +
         scale_color_manual(values=color_subtypes, na.value=NA) +
         guides(color = guide_legend(override.aes = list(size=7)))
      }
   if(!is.null(xlim)) g <- g+ coord_cartesian(xlim=xlim)
   if(!is.null(ylim)) g <- g+ coord_cartesian(ylim=ylim)

      return(g)
   } # end internl bxp function






   ## ::: PREPARE DATA
   ## -----------
   if(!identical(rownames(mat), rownames(fdata))) stop("rownames for matrix and fdata must be idenitcal")
   if(!identical(colnames(mat), rownames(pdata))) stop("colnames for matrix and rownames pdata must be idenitcal")

   default.pdata.column <- if_else(is.null(default.pdata.column), false = default.pdata.column, true = colnames(pdata)[1])
   if(!default.pdata.column %in% colnames(pdata)) stop("pdata object must contain column specified by default.annotation")
   pdata_choices <- colnames(pdata)

   uu <- unlist(sapply(pdata, function(x) !is.numeric(x)))
   if(any(uu)){
      annot_choices <- colnames(pdata)[uu]
      default.annot.column <- if_else(is.null(default.annot.column), false = default.annot.column, true = annot_choices[1])
      if(!default.annot.column %in% annot_choices) stop("pdata object must contain column specified by default.annotation and be numeric")
   }
   if(!any(uu)){
      annot_choices <- NULL
      default.annot.column <- NULL
   }

   #  ::: Shiny App
   ## -----------------
   shinyApp(

   ## UI section
   ## ................
      ui <- fluidPage(
         tags$hr(),

                  fluidRow(
            # column(2,
            #    uiOutput("perplexity_master")),
            column(3, uiOutput("dim2_selection")),
            column(5, uiOutput("annotation_selection"))
            #column(2, sliderInput(inputId = "pointSize", label = "Select point size", min=0.5, max=10, ticks = T, step = 0.5, value = 5)),
            #column(2, sliderInput(inputId = "pointAlpha", label = "Select alpha", min=0, max=1, ticks = T, step = 0.1, value = 0.7))
         ),

         #
         # The PAN TCGA BOXPLOT
         tags$hr(),
         fluidPage(
            column(10, plotOutput('myPlot'))
         ),

         h1('choose gene'),
         # Show the main table used to select ENSG_1 gene id
         # -- fdataTable_main (dataTableOutput)
         # -- ensg_1 (output)
         fluidRow(column(8, DT::dataTableOutput('fdataTable'))),
         tags$hr(),
         h1(textOutput('SYMBOL_1')),
         h3(textOutput('ENSG_1'))
      ),


      ## Define server logic
      ## ::::::::
      server <- function(input, output){

         output$dim2_selection <- renderUI({
            selectInput(inputId = "dim2", label = "Available x-lim varaibles:",
            choices=pdata_choices, selected = default.pdata.column)
            })
         output$annotation_selection <- renderUI({
            selectInput(inputId = "annot", label = "Available annotations if scatterplot:",
            choices=annot_choices, selected = default.annot.column)
            })


         # fdata data table used to select gene
         output$fdataTable = DT::renderDataTable(fdata, server = TRUE, selection = 'single')
         observe({
            u = input$fdataTable_rows_selected
            output$myPlot = renderPlot({
               plotFoo(x = rownames(fdata)[u],  mat = mat, fdata = fdata, pdata = pdata, pdata.column=input$dim2, annot.column = input$annot, altName=fdata[u,fdata.altName.column], y.range.add=1, color.key=color.key)
            })


            }) # end observe
      } # end server
      ) # end shiny app
   }  # end shiny function











