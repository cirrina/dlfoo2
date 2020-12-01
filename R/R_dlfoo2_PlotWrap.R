## Scatterpoints/rankplot with annotations in pdata. variant of the scatterPoint_annotationPlot
#' plotViolinScatterRankBarGeneric
#' @description Generic plot for tSNE (list from tsneWrapper analysis), umap (umap object), or PCA (prcomp object) but also for generic data frames and for plotting Ranks only
#' @description If data.frama supplied. Then:
#' @description - Bar plot: Only when x is data.frame and both dim1 non-numeric characrer/factors. Primary annotation supplied in pdata/pdata.column.
#' @description - Violin: only dim1 numeric OR if eset, gene.id (but not gene.id_2) supplied
#' @description - Rank plot: as for Violin but rank.plot set to TRUE,
#' @description - Scatter plot: dim1 and dim2 numeric OR if eset, gene.id and gene.id.2 supplied
#' @family umap
#' @family tsne
#' @family pca
#' @family Violin
#' @family RankPlot
#' @family Barplot
#' @family ScatterPlot
#' @family expressonset
#' @family PlotWrap
#' @param x results object. list (tSNEres) from run.tsne-analysis, or umap (umap object), or PCA (prcomp object). Can also be a data.frame (sample_id, dim1, dim2) or an expressionSet
#' @param gene.id supplied if x is an ExpressionSet and if scatter or rankplot
#' @param gene.id.2 supplied if x is an ExpressionSet and if scatterplot
#' @param gene.id.signature supplied if x is expression set and if to create a score for violin plot over multiple genes
#' @param pdata Optional. Pairs 'sample_id' with the "annotation' used for group-colors. If not provided an 'annotation' column must be present in x. data frame with, minimally, columns 'sample_id' and a column for annotation (defined by 'annotation.column').
#' @param pdata.column should be a column in the pdata object to use for annotation (colors).
#' @param violin.2nd.group.column specify a columnname if to split violin groups into secondary groups
#' @param drop.na.annotation if to drop NA values from annotation column. defaults TRUE. May be set to FALSE for e.g. barplots to get a better representation of the data.
#' @param my.samples if to highligt samples - character vector of sample names
#' @param plot.title Title of plot
#' @param x.title
#' @param y.title
#' @param axis_names_n if to add number of samples to sample groups (x-axis names). Note: not feasible if downstream facet_wrap
#' @param color.key what color key to use for annotation. Detaults to dlfoo2::color_subtypes (but could be set to e.g. dlfoo2::color_annotations)
#' @param do.gradient set to FALSE if to manually override numeric check for gradient or manual fil
#' @param color.key.gradient what gradient.
#' @param gradient.nomalize TRUE or FALSE if to z normalize annotation value (x-mean(x)) sd(x)
#' @param gradient.squish.probs if to cap squish gradient to remove effects on color scale of outliers. vector of tow percent values. Uses quantile-probs. Defaults to c(0.05, 0.95),i.e. cap scale at 5 and 95 of values
#' @param stripchart.plot set to TRUE if to do stripchart instead of scatterplot. Provided as dim1 (values) and dim2 (groupings).
#' @param rank.plot set to TRUE if to do Rank plot instead of Violin. Defaults to False. Possile if dim1 onöy is  supplied (not dim1 and dim2, then scatterplot).
#' @param rank.decreasing set to True if Rank plot and if rank decerasing on x-scale
#' @param barplot.y.freq TRUE if barplot y-axis should be frequency percent, 0-100 (not total N)
#' @param my.size point size
#' @param my.alpha point alpha
#' @param my.stroke stroke width
#' @param my.stroke.color stroke line colors
#' @param do.boxplot set to true if boxplot insead of violin
#' @param showGuide if to show guide
#' @param coord.fixed if coord_fixed ggplot param
#' @param swap.axes if to swap x and y axes
#' @param xlim.cap  CAPs x - values (dim2) if any outsidde
#' @param xlim.rescale passed on to axis. vector of 2. set upper/lower to NA if only rescaling one value
#' @param ylim.cap  CAPs x - values (dim2) if any outsidde
#' @param ylim.rescale passed on to axis. vector of 2. set upper/lower to NA if only rescaling one value
#' @param zero.lines if to higlight h- and -vlines through origo
#' @param shape.key if to plot annotations in different shapes (connected to annotation)
#' @param zoom.coord list of x and y coordinates to specify/zoom plot coordinates
#' @return a ggplot object
#' @export
plotViolinScatterRankBarGeneric <- function(
   x,
   gene.id=NULL, gene.id.2=NULL, gene.id.signature = NULL,
   pdata=NULL, pdata.column=NULL,
   plot.title=NULL, x.title=NULL, y.title=NULL,
   axis_names_n=T,
   rank.plot = F,
   stripchart.plot =F,
   violin.2nd.group.column = NULL,
   drop.na.annotation = T,
   barplot.y.freq = T,
   rank.decreasing = F,
   my.size=2, my.alpha=0.7, my.stroke=0.05, my.stroke.color="gray10",
   color.key=NULL,
   do.gradient=T, pseudo.gradient=T, gradient.sat = c(-1,1),
   palette.gradient=rev(dlfoo2::palette_gradients[["red2white2blue"]]), gadient.na.col=NA, gradient.normalize=F, gradient.squish.probs=c(0.05, 0.95),
   my.samples=NULL,
   do.boxplot = F,
   swap.axes =F,
   showGuide=F, coord.fixed=F, xlim.cap=NULL, ylim.cap=NULL, xlim.rescale=NULL, ylim.rescale=NULL, zero.lines=F,
   shape.key=NULL, zoom.coord=NULL
   ){


   require(umap)
   require(Rtsne)
   require(dplyr)
   require(tidyverse)
   require(tibble)
   require(ggrepel)
   # outputPlotly=T
   #  if plotly text, add  text =  paste("Point: ", annotation, "\n", sample_id, "\n", annotation2)

   plot_type <- NULL

   ## if umap/tsne/pca then ScatterPlot
   ## -------------
   if(class(x)=="umap"){
      message("Found UMAP results. Begin plotting ... ")
      x_df <- umap2df(x)
      plot_type <- "scatterPlot"
   }

   if(class(x)=="prcomp"){
      message("Found PCA results. Begin plotting ... ")
      plot_type <- "scatterPlot"
      x_df <- pca2df(x)
   }

    if(class(x)=="list"){
     if(!is.null(x$Y)){
      message("Found list with Y matrix slot  ... ")
      message("... assuming x is a tsne reults list. Begin plotting")
      plot_type <- "scatterPlot"
      x_df <- tsne2df(x)
     }}



   ## Expression set
   ## -------------
   if(class(x)=="ExpressionSet"){
      message("expressionset provided")

      x_samples <- intersect(unique(pdata$sample_id), sampleNames(x))


      gene_i = match(gene.id, featureNames(x))

      if(!is.null(gene.id.signature)){ # if multiple genes - gene signature
         gene.id.signature <- unique(gene.id.signature[!is.na(gene.id.signature)])
         gene_i = match(gene.id.signature, featureNames(x))

      }
      if(all(is.na(gene_i))){
         message(gene.id, " not found in feature data for pan_gdc")
         return(NULL)
      }

      if(!is.null(gene.id.2)){
         # if gene_2 is supplied both dim1 and dim2 has to be numeric and set plot to scatterPlot
         plot_type <- "scatterPlot"
         gene_ii = match(gene.id.2, featureNames(x))
         if(is.na(gene_i)){
            message(gene.id.2, " not found in feature data for pan_gdc")
            return(NULL)
         }

         x_df <- data.frame(sample_id=sampleNames(x),
            dim1 = as.numeric(exprs(x[gene_i, x_samples])),
            dim2 = as.numeric(exprs(x[gene_ii, x_samples]))
            )
         x_df <- x_df %>% dplyr::filter(!is.na(dim1)) %>% dplyr::filter(!is.na(dim2))
      }else{
         # if not dim2 is supplied then TWO possible plot types (ViolinPlot or rankPlot) defined by 'rank.plot' param
         ## NOTE: stripchart.plot=T will override rank.plot=T
         if(rank.plot){
            plot_type <- "rankPlot"
            message("... setting plot type to rankPlot")
            } # both dim1 and dim2 has to be neumeric}

         if(!rank.plot) {
            plot_type <- "violinPlot"
            message("... setting plot type to violinplot")
            }

         ## create pdata_plot with matching sample_id as the eset sample names
         if(length(gene_i)==1){ # if only one gene
            x_df <- data.frame(sample_id=x_samples, dim1 = as.numeric(exprs(x[gene_i, x_samples])))
            message("...removing NA-values")
            x_df <- x_df %>% dplyr::filter(!is.na(dim1))
            }
         if(length(gene_i)>1){
            message("\n ... multiple genes supplied - creting gene score based on mean value")
            es_sig <- esetCenterGenes(x[gene_i,])
            x_df <- data.frame(sample_id=x_samples, dim1 = as.numeric(apply(exprs(es_sig[,x_samples]), 2, mean, na.rm=T)))
            message("...removing NA-values")
            x_df <- x_df %>% dplyr::filter(!is.na(dim1))
         }

      } # end if else


   } # end if x is expression set




   ## Data frame
   ## -------------
   if(class(x) == "data.frame"){
      message("Found data frame ... looking for sample_id and dim1 as columns")
      if(!all(c("sample_id","dim1") %in% colnames(x))) stop()

      message("...removing NA-values")
      x_df <- x %>% dplyr::filter(!is.na(dim1)) %>% dplyr::filter(!is.na(sample_id))
      x_samples <- unique(as.character(x_df$sample_id))
      x_df <- x_df[match(x_samples, x_df$sample_id), ]

      if("dim2" %in% colnames(x_df)){
         message(" ... found dim2")
         # if dim2 is supplied both dim1 and dim2 has to be numeric and set plot to scatterPlot
         # OR if stripchart.plot=T then allow di2 as not numeric
         if(stripchart.plot){
            plot_type <- "stripChart"
            x_df$dim2 <- as.character(x_df$dim2)
            }else{
               if(is.numeric(x$dim2)) plot_type <- "scatterPlot"
               if(!is.numeric(x$dim2)) stop("dim2 must be numeric. If barPlot then supply binning annotation (y-axis) through the pdata and pdata.colums parameters while dim1 should be the primary grouping annotation (x-axis) ")
            }
         message(" ... setting plot type to:  ",plot_type)
         x_df <- x_df %>% dplyr::filter(!is.na(dim2))

      }else{
         # if  not dim2 is supplied then 3 possible plot types
         if(is.numeric(x_df$dim1) && rank.plot) plot_type <- "rankPlot"
         if(is.numeric(x_df$dim1) && !rank.plot) plot_type <- "violinPlot"
         if(!is.numeric(x_df$dim1) && !is.numeric(pdata[,pdata.column])){
            plot_type <- "barPlot"
            if(!is.factor(x_df$dim1)) x_df$dim1 <- factor(x_df$dim1)
            x_df <- droplevels(x_df)
            if(is.null(pdata)) stop("pdata must be provided for barPlot")
            }
         message(plot_type)
      } # end if else
   } # end if dataframe



   # if rank plot then create dim2 (x-axis) rank value
   if(rank.plot && plot_type=="rankPlot"){
         x_df$dim2 <- rank(x_df$dim1, ties.method = "random")
         if(rank.decreasing) x_df$dim2 <- rank((-1)*x_df$dim1, ties.method = "random")
      }

   # stop if no x
   if(is.null(plot_type)) stop("x must be an umap, tsne, expressionSet or a data frame with sample_id and dim1 as columns ")



   ## If no pdata provided - then a annotation column must be found in x
   ## ------------
   if(is.null(pdata)){
      message("... no pdata provided\n ... ... checking annotation column in x")
      if(!"annotation" %in% colnames(x_df)) stop("x must contain 'annotation' if no pdata is provided !")
      my_df <- x_df
      ## Check annotation
      if(!is.factor(my_df$annotation) && !is.numeric(my_df$annotation)){
         my_df$annotation <- factor(my_df$annotation)
         }
   }

   ## check pdata - add annotation (group - colors) to my_df
   ## ------------
   ## Check x and pdata
   if(!is.null(pdata)){
      message("...  pdata provided\n ... ... checking annotations in x")
      message("... ... integrating pdata into x")
      if(!"sample_id" %in% colnames(pdata)) stop("pdata object must contain 'sample_id' column")
      if(!(pdata.column %in% colnames(pdata))) stop("pdata.column not in selected pdata object")

      x_samples <-  unique(intersect(pdata$sample_id, x_df$sample_id))
      my_df <- as.data.frame(
            pdata[match(x_samples, pdata$sample_id),] %>%
            mutate(annotation = UQ(rlang::sym(pdata.column))) %>%
            #dplyr::select(!!!c("sample_id", pdata.column), everything()))
            dplyr::select(c(sample_id, annotation, everything())))
      #colnames(my_df)[2] <- "annotation"
      # str(my_df)

      ## Check annotation
      if(!is.factor(my_df$annotation) && !is.numeric(my_df$annotation)){
         my_df$annotation <- factor(my_df$annotation)
         }
      suppressWarnings(my_df <- x_df %>% dplyr::filter(sample_id %in% x_samples) %>% dplyr::left_join(my_df, by="sample_id"))
      # Drop NA from annotation colums (default)
      if(drop.na.annotation) my_df <- my_df %>% dplyr::filter(!is.na(annotation))
      if(!drop.na.annotation){
         temp_levels <- levels(my_df$annotation)
         my_df <- my_df %>%
            mutate(annotation=as.character(annotation)) %>%
            mutate(annotation = if_else(is.na(annotation), "_NA_", annotation)) %>%
            mutate(annotation = factor(annotation, levels = c(temp_levels, "_NA_")))
         }
      my_df <- my_df %>% droplevels() %>% arrange(annotation)
      # str(my_df)

      } # end integrate pdata


   ## Check for Inf values
   if(any(my_df$dim1 %in% c("Inf","-Inf"))){
      message("!! WARNING: 'Inf' values found -possibly after log - replacing thse with max/min of non Inf values")
      my_df <- my_df %>% dplyr::mutate(dim1 = if_else(dim1=="Inf", max(my_df$dim1[!my_df$dim1=="Inf"]), dim1))
      my_df <- my_df %>% dplyr::mutate(dim1 = if_else(dim1=="-Inf", min(my_df$dim1[!my_df$dim1=="-Inf"]), dim1))
   }

   if("dim2" %in% colnames(my_df)){
      if( any(my_df$dim2 %in% c("Inf","-Inf"))){
      message("Inf values found - replacing thse w max/min")
      my_df <- my_df %>% dplyr::mutate(dim2 = if_else(dim2=="Inf", max(my_df$dim2[!my_df$dim2=="Inf"]), dim2))
      my_df <- my_df %>% dplyr::mutate(dim2 = if_else(dim2=="-Inf", min(my_df$dim2[!my_df$dim2=="-Inf"]), dim2))
      }
   }



   ## Add manual samples
   ## ----------
      # if manually added my.samples Then label these
         ## ??? working?
      if(!is.null(my.samples)){
         if(length(which(my.samples %in% my_df$sample_id)) > 0) {
            my_df$plot_text <- NA
            my_df$plot_text[my_df$sample_id %in% my.samples] <- as.character(my_df$sample_id[my_df$sample_id %in% my.samples] )
         }
         if(length(which(my.samples %in% my_df$sample_id)) < 1) {
            message("cant find any of the supplied sample_ids")
            my.samples <- NULL}
      }

   # if fill manual or fill gradient
   if(!is.numeric(my_df$annotation)){
      message("... ... annotation not numeric. setting fill to discrete colors defined by 'color.key'")
      do.fill = "manual"
   }

   ## if to do gradient
   ## -------------
   if(is.numeric(my_df$annotation) && do.gradient && !pseudo.gradient){
      message("... ... annotation is numeric. setting fill to gradient as defined by 'color.key.gradient'")
      do.fill = "gradient"
      if(gradient.normalize) my_df$annotation <- my_df$annotation - mean(my_df$annotation, na.rm=T) / sd(my_df$annotation, na.rm=T)
      #if(cap.outliers)
      #color.key.gradient = palette_gradientRamp(x = my_df$annotation, my.pal = palette.gradient, palette.sat = gradient.sat, na.col = gadient.na.col)
      color.key.gradient = palette.gradient
      my.gradient.limits <- quantile(my_df$annotation, probs=gradient.squish.probs, na.rm=T)
      }


   ## Color.key - Use Pseudo gradient
   ## ----------------------
      # fix for problems with gradient_fill in plotly - define a pseudigradient instaed
   if(pseudo.gradient==T && do.fill == "gradient"){
      message(" ... setting upp pseudo-gradient - manual_fill")
      if(gradient.normalize) my_df$annotation <- my_df$annotation - mean(my_df$annotation, na.rm=T) / sd(my_df$annotation, na.rm=T)
      ## gradient colors is defined here by creating a manual palette for each row
      my_df$annotation_bin <- paste0("annotation",1:nrow(my_df))
      # color_key_pseudo_gradient
      my_df <- my_df %>% dplyr::rename(annotation_num = annotation, annotation = annotation_bin)
      color.key <- dlfoo2::heatcol_creator2(my_df$annotation_num, my.pal = rev(palette.gradient), palette.sat = gradient.sat)$col
      names(color.key) <- my_df$annotation
      do.fill = "manual"
      } # end pseudo.gradient


   ## if do.shape
   do.shape = F
   if(!is.null(shape.key)){
      if(any(levels(my_df$annotation) %in% names(shape.key))){
      message("... ... setting do.shape to TRUE")
      do.shape = T
      shape_key <- shapeKeyFix(shape.key,  as.character(levels(my_df$annotation)))
      shape_key <- shape_key[order(match(names(shape_key), levels(my_df$annotation)))]
      }}

   # zoom coordinates (?? outdated)
   do.zoom = FALSE
   if(!is.null(zoom.coord)){
      message("... ... 'zoom.chord' found. setting do.zoom to T")
      do.zoom=T
      if(class(zoom.coord)!="list") stop("zoom.coord must be list with named character slots x and y character slots of lengtht 2")
      if(length(zoom.coord)!=2) stop("zoom.coord must be list with named character slots x and y character slots of lengtht 2")
      if(!identical(sort(names(zoom.coord)), c("x","y"))) stop("zoom.coord must be list with named character slots x and y character slots of lengtht 2")
      }

   ## if capping x/y values
   ## ------------------
      #!!! note that x is dim2 and y dim1 for historic reasons
   if(!is.null(xlim.cap)){
      if(!is.numeric(xlim.cap) | length(xlim.cap)!=2) stop("xlim.cap is not nueric vector of 2 values")
      my_df$dim2[my_df$dim2 < min(xlim.cap)] <- min(xlim.cap)
      my_df$dim2[my_df$dim2 > max(xlim.cap)] <- max(xlim.cap)
      message("x-axis (dim2) is capped to: ", xlim.cap[1], " and ", xlim.cap[2])
   }
   if(!is.null(ylim.cap)){
       if(!is.numeric(ylim.cap) | length(ylim.cap)!=2) stop("xlim.cap is not nueric vector of 2 values")
      my_df$dim1[my_df$dim1 < min(ylim.cap)] <- min(ylim.cap)
      my_df$dim1[my_df$dim1 > max(ylim.cap)] <- max(ylim.cap)
      message("y-axis (dim1) is capped to: ", ylim.cap[1], " and ", ylim.cap[2])
      }



   ## color.key & shape keys (named vectors w colors)
   ## ---------
   if(!drop.na.annotation) color.key['_NA_'] <- NA
   color_key_x <-  colorKeyFix(color.key, as.character(levels(factor(my_df$annotation))))

   if(!is.null(violin.2nd.group.column) & plot_type=="violinPlot"){
      if(!c(violin.2nd.group.column %in% colnames(my_df))) stop(paste(violin.2nd.group.column, "not among pdata columns"))
      color_key_x <-  colorKeyFix(color.key, as.character(levels(factor(my_df[,violin.2nd.group.column]))))
   }

   ## if swap axes
   ## ------------
   if(swap.axes){
      my_df <- my_df %>%
         dplyr::mutate(dim2_temp=dim1) %>%
         dplyr::mutate(dim1_temp=dim2) %>%
         dplyr::select(-dim1, -dim2) %>%
         dplyr::rename(dim1=dim1_temp, dim2=dim2_temp) %>%
         dplyr::select(sample_id, dim1, dim2, everything())

      x.title_temp <- y.title
      y.title <- x.title
      x.title <- x.title_temp
      }


   ##    Plot ggplot
   ## :::::::::::::::::::
   # if(is.null(plot.title)) plot.title = pdata.column

   ## scatterPlot or rankPlot
   ##
   message("... plotting")
   if(plot_type %in% c("scatterPlot","rankPlot")){

       g <- ggplot(my_df) +
            aes(x=dim2, y=dim1)
            # text =  paste("Point: ", annotation, "\n", sample_id, "\n", annotation2) :: if plotly

       if(do.shape & do.fill == "manual"){
         g <- g +
            geom_point(aes(fill=annotation, shape=annotation), colour=my.stroke.color, stroke=my.stroke, size=my.size, alpha=my.alpha) +
            scale_fill_manual(values=color_key_x, na.value=NA) +
            scale_shape_manual(values = as.numeric(shape_key), na.value = 4) +
            guides(color = guide_legend(override.aes = list(size=7)))
      }

      if(!do.shape & do.fill == "manual"){
          g <- g +
            geom_point(aes(fill=annotation), shape=21, stroke=my.stroke, colour=my.stroke.color, size=my.size, alpha=my.alpha) +
            scale_fill_manual(values=color_key_x, na.value=NA)+
            guides(fill=guide_legend(title=pdata.column))
            }
      if(!do.shape & do.fill == "gradient"){
            g <- g +
            geom_point(aes(fill=annotation), shape=21, stroke=my.stroke, colour=my.stroke.color, size=my.size, alpha=my.alpha) +
            # scale_fill_gradientn(colours=color.key.gradient, na.value=NA) +
            scale_fill_gradientn(colours=color.key.gradient, na.value=NA, limits=my.gradient.limits, oob=scales::squish) +
            guides(fill=guide_legend(title=pdata.column))
      }

      ## add dim1 and dim2 lines for selected sample name(s)
      if(!is_null(my.samples)){
         # g <- g + geom_label_repel(aes(label = plot_text),  box.padding   = 0.35, point.padding = 0.5, segment.color = 'grey50')
         g <- g +
            geom_hline(yintercept = data.frame(my_df %>% dplyr::filter(!is.na(plot_text)))$dim1) +
            geom_vline(xintercept = data.frame(my_df %>% dplyr::filter(!is.na(plot_text)))$dim2)
         }

   } # end if scatter/rankplot

   ## Stripchart
    if(plot_type %in% c("stripChart")){

      my_sizes <- my_df %>% group_by(dim2) %>% summarize(num=n())
      if(axis_names_n){ axis_names <- paste0(my_sizes$dim2, ", ", "n=", my_sizes$num)
         }else{
         axis_names <- paste0(my_sizes$dim2)
      }
      # if(stripchart.box){
      #    g <- ggplot(my_df) +
      #       message(" ... stripchart.box=T, adding box")
      #       aes(x=dim2, y=dim1) +
      #       geom_boxplot() +
      #       theme(legend.position='none') +
      #       geom_jitter(aes(fill=annotation), position=position_jitter(0.2), shape=21, stroke=my.stroke, colour=my.stroke.color, size=my.size, alpha=my.alpha) +
      #       scale_fill_manual(values = color_key_x) +
      #       scale_x_discrete( labels = axis_names)
      # }

      g <- ggplot(my_df) +
         aes(x=dim2, y=dim1, fill=annotation) +
         theme(legend.position='none') +
         geom_jitter(position=position_jitter(0.2), shape=21, stroke=my.stroke, colour=my.stroke.color, size=my.size, alpha=my.alpha) +
         scale_fill_manual(values = color_key_x) +
         scale_x_discrete( labels = axis_names)

   }


   # Vioin
   if(plot_type %in% c("violinPlot")){

      my_sizes <- my_df %>% group_by(annotation) %>% summarize(num=n())
      if(axis_names_n){ axis_names <- paste0(my_sizes$annotation, ", ", "n=", my_sizes$num)
         }else{
         axis_names <- paste0(my_sizes$annotation)
      }

      # axis_num <- my_sizes$num
      # my_df %>% left_join(my_sizes, by="annotation") %>% mutate(myaxis = paste0(dim1, "\n", "n=", x.group.num))

      if(is.null(violin.2nd.group.column)){
         if(!do.boxplot){
            g <- ggplot(my_df) +
               aes(x=annotation, y=dim1, fill=annotation) +
               # text = paste("Group: ", annotation)
               theme(legend.position='none') +
               geom_violin(alpha=0.5, width=1, trim = TRUE, scale = "width", adjust = 0.5) +
               geom_boxplot(width=0.2, outlier.colour="grey25", notch = FALSE, notchwidth = .75, alpha = my.alpha, colour = "grey10", outlier.size=0.3) +
               scale_fill_manual(values = color_key_x) +
               scale_x_discrete( labels = axis_names)
         }else{
            g <- ggplot(my_df) +
               aes(x=annotation, y=dim1, fill=annotation) +
               # text = paste("Group: ", annotation)
               theme(legend.position='none') +
               geom_boxplot(outlier.colour="grey25", notch = FALSE, notchwidth = .75, alpha = my.alpha, colour = "grey10", outlier.size=0.3) +
               scale_fill_manual(values = color_key_x) +
               scale_x_discrete( labels = axis_names)
            }
         }
      if(!is.null(violin.2nd.group.column)){
         message(" ... ... viollin () plot with 2nd grouping variable")
         g <- ggplot(my_df) +
            aes(x=annotation, y=dim1, fill = UQ(rlang::sym(violin.2nd.group.column))) +
            theme(legend.position='none') +
            # geom_violin(alpha=0.5, width=1, trim = TRUE, scale = "width", adjust = 0.5) +
            geom_boxplot(outlier.colour="grey25", notch = FALSE, notchwidth = .75, alpha = my.alpha, colour = "grey10", outlier.size=0.3) +
            scale_fill_manual(values = color_key_x) +
            scale_x_discrete(labels = axis_names)
      }
   }

   # barplot
   if(plot_type %in% c("barPlot")){

      ftab <- my_df %>% group_by(dim1, annotation) %>% summarise(n=n()) %>% mutate(f=n/sum(n))
      # ftab$annotation <- factor(ftab$annotation, levels=rev(ftab$annotation))
      # ftab <- ftab %>% arrange(!annotation)
      my_sizes <- my_df %>% group_by(dim1) %>% summarize(num=n())
      axis_names <- paste0(my_sizes$dim1, ", ", "n=", my_sizes$num)

      if(barplot.y.freq == T ){
         # factor(annotation, levels = c(NA, as.character(levels(my_df$annotation))), exclude = NULL)

         g <- ggplot(ftab, aes(fill=forcats::fct_rev(annotation), x=dim1, y=n)) +
            geom_bar(position="fill", stat="identity") +
            scale_fill_manual(values =color_key_x) +
            scale_x_discrete(labels = axis_names)
      }

      if(barplot.y.freq == F ){
         g <- ggplot(ftab, aes(fill=forcats::fct_rev(annotation), x=dim1, y=n)) +
            geom_col() +
            scale_fill_manual(values =color_key_x) +
            scale_x_discrete(plabels = axis_names)
      }
         # text = paste("Group: ", dim1, "\n annotation: ", annotation)


   } # end barplot

   ## plot stuff
   g <- g + ggtitle(label = plot.title) +
            labs(x=x.title, y=y.title, colour='') +
            theme(axis.title.y = element_text(size=9)) +
            theme(axis.title.x = element_text(size=9)) +
            theme(legend.text = element_text(size=9)) +
            theme(title = element_text(size=9)) +
            theme(axis.text.x=element_text(angle = 45, hjust = 1))

   ## if rescaling of x/y
   if(!is.null(xlim.rescale)){
      if(!is.numeric(xlim.rescale[!is.na(xlim.rescale)]) | length(xlim.rescale)!=2) stop("xlim.rescale is not numeric vector of 2 values or NA adn one value ")
      g <- g + xlim(xlim.rescale)
      }
   if(!is.null(ylim.rescale)){
      if(!is.numeric(ylim.rescale[!is.na(ylim.rescale)]) | length(ylim.rescale)!=2) stop("ylim.rescale is not numeric vector of 2 values or NA adn one value ")
      g <- g + ylim(ylim.rescale)
      }

   # Guide
   if(showGuide==F){
      g <- g + guides(fill=FALSE, shape=FALSE, color=FALSE) +
               theme(legend.position='none')
   }

   if(coord.fixed) {g = g + coord_fixed(ratio = 1)}

   if(do.zoom){g = g + coord_cartesian(xlim = zoom.coord$x, ylim = zoom.coord$y)}

   if(zero.lines) g <- g + geom_hline(yintercept=0, linetype="dashed", color="black", size=0.5) + geom_vline(xintercept=0, linetype="dashed", color="black", size=0.5)
      message("... done")

   return(g)


} # end violinScatterRankBarPlotGeneric








#' @description Wrapper for ggplot heatmap suing geom_raster
#' @description NOTE! cannot export thse as vetor graphics!! Need geom_rect??
#' @description Input (x) must be a matrix with n sample columns and m feature/gene rows
#' @description Note: since plotly cant cope with plotting pdata (character/factors) in geom_raster and custom colorings - groups have to be defined/marked creatively. here by adding vertical line using 'pdata.column.vert.lines'
#' @description also: plotly.js does not (yet) support horizontal legend items. You can track progress here: https://github.com/plotly/plotly.js/issues/53
#' @family Heatmap
#' @family PlotWrap
#' @param my_mat matrix for heatmap with n samples and n features
#' @param pdata if to add sample information to the ggplot object (used for downstream ggplotly tootips)
#' @param pdata.column.vert.lines A pdata column that define groups. if to draw vertical lines that separate groups. Will only work if x is sorted properly and if cluster.samples==F.
#' @param pdata.column.vert.lines.lwd
#' @param normalize.rows if to z normalize rows (divide by sd)
#' @param palette.gradient palette for gradient coloring
#' @param cluster.features if to cluster rows/features
#' @param cluster.samples if to cluster samples/columns
#' @param strip.background if to leave backround blank
#' @param col.labels if to draw ticks and labes for x axis
#' @param nrow.max.for.labels the max rows for x if to plot y-labels
#' @param gradient.squish.probs how palette is squished
#' @param gradient.limits if to manually set the gradient limits (not by squish values in matrix)
#' @param x.title
#' @param y.title
#' @param plot.title
#' @param showGuide defaults to F
#' @param ... parameters passed on to function violinScatterRankBarPlotGeneric
#' @return a ggplot object. Z vill represent the matix values
#' @export
heatmapGeomRaster <- function(
   my_mat =NULL,
   pdata = NULL,
   pdata.column.vert.lines = NULL,
   pdata.column.vert.lines.lwd = 0.5,
   #color.key.pdata = NULL, # not suppored since ggploty cant understand scale_fill_manual
   #pdata.heatmap = F, # not suppored since ggploty cant understand scale_fill_manual
   normalize.rows = F,
   cluster.features = T,
   cluster.samples = T,
   nrow.max.for.labels = 50,
   col.labels = F,
   strip.background = F,
   plot.title=NULL,
   x.title=NULL, y.title=NULL,
   palette.gradient = rev(c(rev(RColorBrewer::brewer.pal(9,"YlOrRd")[1:9]), rep("white",3), RColorBrewer::brewer.pal(9,"Blues")[1:9])),
   gradient.limits = NULL,
   gradient.squish.probs=c(0.025, 0.975),
   showGuide=F,  ...
   ){

   x <- my_mat
   gadient.na.col=NA
   if(!is.matrix(x)) stop("x must be amatrix")

   ## match up columns (sample ids) with pdata
   ## -------------------------------
      # pdata must include sample_id column that match colnames of x
   if(!is.null(pdata)){
      if(!"sample_id" %in% colnames(pdata)) stop("pdata must include 'sample_id' column with values matching colnames of x")
      pdata <- pdata %>% dplyr::filter(sample_id %in% colnames(x))
      pdata <- pdata[match(colnames(x), pdata$sample_id),]
      if(ncol(pdata)>10) message("WARNING: pdata contains >10 columns - are these really needed for downstream plotly?")
      #if(!is.null(pdata.column.plotly.text)){
      #   if(!pdata.column.plotly.text %in% colnames(pdata)) stop(paste0("pdata.column.plotly.text not in colnames of pdata:  ", pdata.column.plotly.text))
         # pdata$plotly_text <- as.character( pdata[,pdata.column.plotly.text])
      #}
   }

   ## if normalize.rows
   if(normalize.rows){
      ## remove genes that do not vary (with sd=0). Not compatible with hca clustering.
      ## -------------------------
      u <- apply(x, 1, sd)
      if(length(which(u==0))>0){
         message(" .. WARNING: genes with zero standard deviation found. Removing these")
         x <- x[!u==0,]
      }
      message("... normalizing rows")
      x <- t(apply(x, 1, function(y) y/sd(y)))
   }

   ## cluster rows and columns
   ## -----------------------
   if(cluster.features){
      ## remove genes that do not vary (with sd=0). Not compatible with hca clustering.
      ## -------------------------
      u <- apply(x, 1, sd)
      if(length(which(u==0))>0){
         message(" .. WARNING: genes with zero standard deviation found. Removing these")
         x <- x[!u==0,]
      }
      message("... clustering features")
      hca_rows <- dlfoo2::hca(t(x), plot.dendogram = F)
      x <- x[hca_rows$hclust$order, ]
   }
   if(cluster.samples){
      message("... clustering samples")
      hca_cols <- dlfoo2::hca(x, plot.dendogram = F)
      x <- x[, hca_cols$hclust$order]
      pdata <- pdata[hca_cols$hclust$order, ]
   }



   ## Matrix as heatmap
   ## :::::::::::::::::::::::::::::
   x_df <- x %>%
      as_tibble() %>%
      mutate(feature_id = rownames(x)) %>%
      gather(key="sample_id", value="Z", -feature_id) %>%
      mutate(Y=factor(feature_id, levels=unique(rownames(x)))) %>%
      mutate(X=factor(sample_id, levels=unique(colnames(x))))
   x_df$Y <- as.numeric(x_df$Y)
   x_df$X <- as.numeric(x_df$X)

   # if to add pdata group to plotly text
   # if(!is.null(pdata.column.plotly.text)){
   #    u <- match(x_df$sample_id, pdata$sample_id)
   #    x_df <- x_df %>% mutate(plotly_text = pdata[,pdata.column.plotly.text][u])
   # }

   ## add pdata to x_df
   ## -----------------
   if(!is.null(pdata)){
      x_df <- data.frame(x_df) %>% dplyr::left_join(pdata, by="sample_id")
      x_df <- as_tibble(x_df)
   }


   ## set up color gradient
   ## -------------
      message("... ... setting fill to gradient as defined by 'palette.gradient'")
      #if(gradient.normalize) my_df$annotation <- my_df$annotation - mean(my_df$annotation, na.rm=T) / sd(my_df$annotation, na.rm=T)
      #if(cap.outliers)
      #color.key.gradient = palette_gradientRamp(x = my_df$annotation, my.pal = palette.gradient, palette.sat = gradient.sat, na.col = gadient.na.col)
      color.key.gradient = palette.gradient
      my.gradient.limits <- quantile(x_df$Z, probs=gradient.squish.probs, na.rm=T)
      if(!is.null(gradient.limits)) my.gradient.limits = gradient.limits

   ## Axis labels
   ## -------
   y_labs <- unique(x_df$feature_id)
   y.n <- length(y_labs)


   ## ggPlot Heatmap
   ## ---------------
   message("... plotting")
   g_mat <-  ggplot(data.frame(x_df), aes(x=X, y=Y, fill= Z)) +
      geom_raster() +
      scale_fill_gradientn(colours=color.key.gradient, na.value=NA, limits=my.gradient.limits, oob=scales::squish)

   if(!col.labels){
      g_mat <- g_mat +
         theme(axis.ticks.x = element_blank(),
               axis.text.x = element_blank()
               )
      }


   if(strip.background){
      g_mat <- g_mat + theme_bw() +
         theme(panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            # axis.line = element_line(colour = "black")
            axis.line = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank()
            )
      }



   g_mat <- g_mat +
      ggtitle(label = plot.title) +
      labs(x=x.title, y=y.title, colour='black') +
      #theme(axis.title.y = element_text(size=9)) +
      #theme(axis.title.x = element_text(size=9)) +
      theme(legend.text = element_text(size=9)) +
      theme(title = element_text(size=9))
      # theme(axis.text.x=element_text(angle = 45, hjust = 1))
      # theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
      #theme(panel.background = element_rect(fill = NA, colour = NA)) +
      #theme(panel.grid.major =  element_line(colour = NA), panel.grid.minor =  element_line(colour = NA)) +
      #theme(strip.text = element_text(size=9)) +
      #labs(x="", y="" )

   if(nrow(x) < nrow.max.for.labels) g_mat <- g_mat + scale_y_continuous(limits=c(0.5, y.n+0.5), breaks=c(1:y.n), labels=as.character(y_labs), expand=c(0,0))

   ## if pdata.column.vert.lines, then identify groups, plot vert lines and add x-labels
   if(!is.null(pdata.column.vert.lines)){
      if(!pdata.column.vert.lines %in% colnames(pdata)) stop("'pdata.column.vert.lines' among pdata columns:  ", pdata.column.vert.lines)
      if(!is.factor(pdata[,pdata.column.vert.lines])) stop("'pdata.column.vert.lines' must be factor:  ", pdata.column.vert.lines)
      my_groups <- levels(pdata[,pdata.column.vert.lines])
      if(any(is.na(my_groups))) stop("'pdata.column.vert.lines' column must not contain NA values")
      my_x_lines <- unlist(lapply(my_groups, function(x){ which(pdata[,pdata.column.vert.lines]==x)[1]}))
      my_x_lines <- c(my_x_lines-0.5, nrow(pdata)+0.5)
      my_x_mid <- my_x_lines[-1] - c(my_x_lines[-1] - my_x_lines[-length(my_x_lines)])/2
      x.n <- length(my_x_mid)
      g_mat <- g_mat + geom_vline(xintercept = my_x_lines, lwd=pdata.column.vert.lines.lwd, col="black")
      g_mat <- g_mat + scale_x_continuous(breaks =my_x_mid, labels=as.character(my_groups), expand=c(0,0))

   }


    if(showGuide==F) g_mat <- g_mat +
         guides(fill=FALSE, shape=FALSE, color=FALSE) +
         theme(legend.position='none')

      #if(is.null(pdata.column.plotly.text)) g_mat <- g_mat + (aes(text=paste0(sample_id, "\n",feature_id,"\n", round(Z,2))))
      #if(!is.null(pdata.column.plotly.text)) g_mat <- g_mat + (aes(text=paste0(sample_id, "\n",feature_id,"\n", round(Z,2), "\n", plotly_text)))
      # ggplotly(g_mat)



   return(g_mat)
} # end function heatmapGeomRaster




#' @description Plot annotations in heatmap format
#' @description Input should be a data frame with columns representing annotations to be plotted
#' @description a color key (named list of color character vectors) must beprovided for each annotation
#' @description The data frame shpuld be sorted to desired layout
#' @description the function will automatically determine if annotations are numeric or binary (character or factor)
#' @family heatmap
#' @family Bubbleblot
#' @family annotations
#' @family PlotWrap
#' @param annot.tab data frame of annotations to be plotted. must include 'sample_id' column as well. if numeric create pseudo gradients are created
#' @param color.keys named list of colors. if non-numeric annotation, the color vector must be a named vector with unique values of that annotation
#' @param color.keys.sat named list containing saturation limits for continous color keys. if not present the function sets top/bottom 5 percentile as limits
#' @param split.column if to split up into multiple panels baesd on one of the annotations. If so the function returns a list of gg-plots.
#' @param reverse.y if to reverse y axis order
#' @param reverse.x if to reverse y axis order
#' @param size.ylabels
#' @param y.title
#' @param plot.title
#' @param my.stroke
#' @param my.stroke.color
#' @return a ggplot object
#' @export

heatmapAnnotations <- function(
      annot.tab,
      color.keys,
      color.keys.sat=NULL,
      split.column = NULL,
      y.title = NULL,
      plot.title=NULL,
      size.ylabels = 3,
      # score.title ="score"
      # size.title ="size"
      ## shoild be pre-sorted!
      # border.col = NA
      reverse.x=F,
      reverse.y=T
) {

   if(!"sample_id" %in% colnames(annot.tab)) stop("annot.tab must include 'sample_id' column")
   my_annots <- colnames(annot.tab %>% dplyr::select(-sample_id))
   if(!all(my_annots %in% names(color.keys))) stop("not all columns in plot_df are present in color.keys")
   my_annot_classes <- lapply(my_annots, function(x){
     class(annot.tab[,x] )
   })
   names(my_annot_classes) <- my_annots


   ## first create annot list to query each individual annotation and assign colors (and pseudo color gradients)
   annot_list <- lapply(my_annots, function(x) annot.tab %>% dplyr::select(c("sample_id",x)))
   names(annot_list) <- my_annots

   annot_list <- lapply(my_annots, function(x){
      y <- annot_list[[x]]
      if(class(y[,x]) %in% c("character","factor")){
         y$col <- color.keys[[x]][as.character(y[,x])]
      }else{ # expect numeric and create gradient
         if(!is.null(color.keys.sat[[x]])){
            gradient.limits <- color.keys.sat[[x]]
         }else{
         gradient.squish.probs=c(0.05, 0.95)
         gradient.limits <- quantile(as.numeric(y[,x]), probs=gradient.squish.probs, na.rm=T)
         }
         color.key.gradient <- as.character(color.keys[[x]])

         y$col <- dlfoo2::heatcol_creator2(as.numeric(y[,x]), my.pal = color.key.gradient, palette.sat = gradient.limits, na.col = NA)$col
      }
      # colnames(y) <- c("sample_id","annot","col")

      yy <- y %>% dplyr::select(sample_id, col)
      yy$annot <- x
      return(yy)
   })
   names(annot_list) <- my_annots

   annot_melt <- dplyr::bind_rows(annot_list)
   annot_melt$size <- 1
   annot_melt$pch <- 15
   annot_melt$alpha <- 1

   ## Create x (sample_id) and y (annot) values
   annot_melt <- annot_melt %>%
      mutate(y=factor(annot, levels=my_annots)) %>%
      mutate(x=factor(sample_id, levels=unique(sample_id)))


   if(reverse.x) annot_melt <- annot_melt %>% mutate(x = factor(x, levels=rev(levels(x))))
   if(reverse.y) annot_melt <- annot_melt %>% mutate(y = factor(y, levels=rev(levels(y))))
   y_labs <- levels(annot_melt$y)
   y.n <- length(y_labs)

   ## If to produce multiple plots based on one of the annotations
   ## -----------------------------------
      ## e.g. if annotation groups are of very different sizes - split and scale outside
      ## split.column = "tax_simp2"
   if(is.null(split.column)) plot_list <- list(annot_melt)
   if(!is.null(split.column)){
      annot_melt <- annot_melt %>%  dplyr::left_join(annot.tab %>% dplyr::select(sample_id, split.column))
      if(!split.column %in% colnames(annot.tab)) stop()
      # my_splits <- as.character(unique(annot_melt[,split.column]))
      plot_list <- split(annot_melt[,], annot_melt[,split.column])
      # str(plot_list)
   }

   gg_list <- list()

   ## loop all data frames in plot_list
   for(i in 1:length(plot_list)){
      my_melt <- plot_list[[i]] %>% droplevels()
      x_labs <- levels(my_melt$x)
      x.n <- length(x_labs)
      # color_key_x
      my_melt$annot_color_id <- paste0("annot_", 1:nrow(my_melt))
      color_key_x <- my_melt$col
      names(color_key_x) <- my_melt$annot_color_id

         # str(my_melt)
      ## Plot rects
      ## ------------------
       g <- ggplot(my_melt) +
            aes(y=as.numeric(y), x=as.numeric(x), fill=annot_color_id, color=annot_color_id) +
            guides(fill=FALSE, shape=FALSE, color=FALSE) +
            theme(legend.position='none') +
            #scale_color_gradientn(name=score.title, colours = rev(gradient.pal), values=c(0,0.5,1)) +
            geom_point(alpha=1, pch=15, size=3) +
            scale_color_manual(values=color_key_x)

      #g


      g <- g +
         ggtitle(label = plot.title) +
         labs(y=y.title, colour='') +
         theme(axis.title.y = element_text(size=size.ylabels)) +
         theme(axis.title.x = element_text(size=9)) +
         theme(legend.text = element_text(size=6)) +
         theme(title = element_text(size=9)) +
         # theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
         theme(axis.text.x=element_blank()) +
         theme(axis.ticks.x=element_blank()) +
         # theme(panel.grid.major =  element_line(colour = 'lightsteelblue2'), panel.grid.minor =  element_line(colour = 'mistyrose3')) +
         theme(panel.grid.major =  element_blank(), panel.grid.minor =  element_blank()) +
         # theme(panel.background = element_rect(fill = 'mintcream', colour = NA)) +
         # theme(strip.text = element_text(size=5)) +
         #theme(legend.title =  element_text(size=9, face="bold")) +
         #theme(legend.direction = "horizontal", legend.box = "vertical") +
         # guides(size=guide_legend(title=size.title, direction="horizontal")) +
         labs(x="", y="" ) +
         scale_x_continuous(limits=c(0.5, x.n+0.5), breaks=c(1:x.n), labels=dlfoo2::label_wrap_mod(as.character(x_labs)), expand=c(0,0)) +
         scale_y_continuous(limits=c(0.5, y.n+0.5), breaks=c(1:y.n), labels=dlfoo2::label_wrap_mod(as.character(y_labs)), expand=c(0,0))
      # g

       #g
       gg_list[[i]] <- g
   } # end i all plot_lists

   if(length(gg_list)==1) gg_list <- gg_list[[1]]
   if(length(gg_list)>=1) message(" *** Output is a list of ggplot objects defined by 'split.column' *** ")
   return(gg_list)
   } # end functionheatmapAnnotations



#' @description Wrapper (shiny) to plot pre-calculated EA enrichment results for pair-wise comparisons as a table - bubble graph
#' @description Uses violinScatterRankBarPlotGeneric to plot a scatter plot
#' @description The coefficients c1 and c2 are dim1 and dim2 (factors - x/y positions)
#' @description p.adjust is translated inot annot - which then is used for a annot_bin and psedo-gradient
#' @family violinScatterRankBarPlotGeneric
#' @family Entrichment analysis
#' @family Bubbleblot
#' @family ScatterPlot
#' @family dotplot
#' @family PlotWrap
#' @param ea.tab results table from enrichment analysis. Must contaoin coef1 and coef2 columns (as factors for choosing plot order)
#' @param ea.term name of the one term to plot
#' @param gradient.pal palette for gradient (p.value) coloring
#' @param gradient.sat
#' @param reverse.y if to reverse y axis order
#' @param reverse.x if to reverse y axis order
#' @param force.levels.c1 if to force order on c1 axis
#' @param force.levels.c2 if to force order on c2 axis
#' @param drop.levels if to drop levels from c1/c2 when these are factors (if not to plot empty resutls)
#' @param x.title
#' @param y.title
#' @param plot.title
#' @param my.stroke
#' @param my.stroke.color
#' @param my.size
#' @param my.alpha
#' @param showGuide defaults to F
#' @param coord.fixed defaults to T
#' @param ... parameters passed on to function violinScatterRankBarPlotGeneric
#' @return a ggplot object
#' @export
enrichmentCoocPlotrWrap <- function(
   ea.tab=NULL,
   ea.term=NULL,
   gradient.pal = c(rev(RColorBrewer::brewer.pal(9,"YlOrRd")[1:9]), rep("white",3), RColorBrewer::brewer.pal(9,"BuGn")[1:9]),
   gradient.sat = 5,
   reverse.y=F, reverse.x=F,
   force.levels.c1 = NULL,
   force.levels.c2 = NULL,
   drop.levels=F,
   x.title=NULL, y.title=NULL,
   plot.title = NULL,
   my.stroke=1, my.stroke.color=NA, my.size=1, my.alpha=1,
   showGuide=F, coord.fixed=T,...
      ){

      # outputPlotly=F
      pseudo.gradient=T
      if(is.null(plot.title)) plot.title=ea.term
      term_df <- ea.tab %>% dplyr::filter(ID==ea.term)
      if(nrow(term_df)<1) stop()

      # remove if duplicated terms, that is if both upreg AND downreg for one comparison - keep the most significant
      term_df <- term_df %>% arrange(p.adjust) %>%
         dplyr::filter(!duplicated(paste(ontology, contrast, Description)))

      ## Here
         # coef1 - Y / dim 1
         # coef2 - X / dim 2
         # p.adjust - score//color//annotation
         # Cluster - will shift sign for p.adjust to tranlate to palette
         # GeneRatio // possibly size of circles
      term_df$size <- sapply(term_df$GeneRatio, function(x) eval(parse(text=x)))
      if(!is.null(force.levels.c1)){
         term_df$coef1 <- factor(term_df$coef1, levels=force.levels.c1)
      }
      if(!is.null(force.levels.c2)){
         term_df$coef2<- factor(term_df$coef2, levels=force.levels.c2)
      }

      if(!is.factor(term_df$coef1)) term_df$coef1 <- factor(term_df$coef1)
      if(!is.factor(term_df$coef2)) term_df$coef2 <- factor(term_df$coef2)
      if(reverse.y) term_df <- term_df %>% mutate(coef1 = factor(coef1, levels=rev(levels(coef1))))
      if(reverse.x) term_df <- term_df %>% mutate(coef2 = factor(coef2, levels=rev(levels(coef2))))

      term_df <- term_df %>%
         mutate(dim1=as.numeric(coef1), dim2=as.numeric(coef2), annot=p.adjust) %>%
         mutate(annot= (-1)*log10(annot)) %>%
         mutate(annot=if_else(Cluster=="Dnreg", annot*-1, annot)) %>%
         dplyr::select(dim1, dim2, annot, size, coef1, coef2, p.adjust)

      if(drop.levels) term_df <- droplevels(term_df)

      ## gradient colors is defined here by creating a manual palette for each row
      term_df$annot_bin <- paste0("annot",1:nrow(term_df))
      gradient_sat <- c(-1*gradient.sat, gradient.sat)
      color_key_pseudo_gradient <- dlfoo2::heatcol_creator2(term_df$annot, my.pal = rev(gradient.pal), palette.sat = gradient_sat)$col
      names(color_key_pseudo_gradient) <- term_df$annot_bin

      my_df <- term_df %>%
         arrange(dim2, dim1)

      x_labs <- levels(my_df$coef2)
      y_labs <- levels(my_df$coef1)
      x.n <- length(x_labs)
      y.n <- length(y_labs)

      # START EA CoOc PLOT
      ## ................
         message("... plotting")

         # If plotly then replace gradient with manual pseudogradient
         g <- ggplot(my_df) +
                  aes(x=dim2, y=dim1)
            # if plotly: text =  paste(coef1, " vs ", coef2, "\np.adj: ", signif(p.adjust, 2) )


         if(pseudo.gradient) g <- g +
                  geom_point(aes(color=annot_bin), size=my.size, alpha=my.alpha) +
                  scale_color_manual(values=color_key_pseudo_gradient, na.value=NA)


         ## plot stuff
         g <- g + ggtitle(label = plot.title) +
                  labs(x=x.title, y=y.title, colour='') +
                  theme(axis.title.y = element_text(size=9)) +
                  theme(axis.title.x = element_text(size=9)) +
                  theme(legend.text = element_text(size=9)) +
                  theme(title = element_text(size=9)) +
                  theme(axis.text.x=element_text(angle = 45, hjust = 1)) +

                  theme(panel.grid.major =  element_line(colour = 'lightsteelblue2'), panel.grid.minor =  element_line(colour = 'mistyrose3')) +
                  theme(panel.background = element_rect(fill = 'mintcream', colour = NA)) +

                  theme(panel.grid.major =  element_line(colour = NA), panel.grid.minor =  element_line(colour = 'black')) +
                  theme(strip.text = element_text(size=9)) +
                  theme(legend.title =  element_text(size=9, face="bold")) +
                  theme(legend.direction = "horizontal", legend.box = "vertical") +
                  # guides(size=guide_legend(title=size.title, direction="horizontal")) +
                  labs(x="", y="" ) +
                  scale_x_continuous(limits=c(0.5, x.n+0.5), breaks=c(1:x.n), labels=dlfoo2::label_wrap_mod(as.character(x_labs)), expand=c(0,0)) +
                  scale_y_continuous(limits=c(0.5, y.n+0.5), breaks=c(1:y.n), labels=dlfoo2::label_wrap_mod(as.character(y_labs)), expand=c(0,0))


           if(showGuide==F) g <- g +
               guides(fill=FALSE, shape=FALSE, color=FALSE) +
               theme(legend.position='none')

            if(coord.fixed) g <- g +
               coord_fixed(ratio = 1)

        return(g)

} # exportenrichmentCoocPlotrWrap





#' @description Wrapper for shiny to plot FGSEA enrichment results for pair-wise comparisons as a table - bubble graph
#' @description Uses violinScatterRankBarPlotGeneric to plot a scatter plot
#' @description The coefficients coef1 and coef2 are dim1 and dim2 (factors - x/y positions)
#' @description p.adjust is translated inot annot - which then is used for a annot_bin and psedo-gradient
#' @family violinScatterRankBarPlotGeneric
#' @family Entrichment analysis
#' @family fgsea
#' @family Bubbleblot
#' @family ScatterPlot
#' @family dotplot
#' @family PlotWrap
#' @param ea.tab results table from enrichment analysis. Must contaoin coef1 and coef2 columns (as factors for choosing plot order)
#' @param ea.term name of the one term to plot
#' @param gradient.pal palette for gradient (p.value) coloring
#' @param gradient.sat
#' @param reverse.y if to reverse y axis order
#' @param reverse.x if to reverse y axis order
#' @param force.levels.c1 if to force order on c1 axis
#' @param force.levels.c2 if to force order on c2 axis
#' @param drop.levels if to drop levels from c1/c2 when these are factors (if not to plot empty resutls)
#' @param x.title
#' @param y.title
#' @param plot.title
#' @param my.stroke
#' @param my.stroke.color
#' @param my.size
#' @param my.alpha
#' @param showGuide defaults to F
#' @param coord.fixed defaults to T
#' @param ... parameters passed on to function violinScatterRankBarPlotGeneric
#' @return a ggplot object
#' @export
fgseaCoocPlotrWrap <- function(
   ea.tab=NULL,
   ea.term=NULL,
   gradient.pal = c(rev(RColorBrewer::brewer.pal(9,"YlOrRd")[1:9]), rep("white",3), RColorBrewer::brewer.pal(9,"BuGn")[1:9]),
   gradient.sat = 5,
   plot.title=NULL,
   # border.col = NA
   reverse.y=F, reverse.x=F,
   force.levels.c1 = NULL,
   force.levels.c2 = NULL,
   drop.levels=F,
   x.title=NULL, y.title=NULL,
   my.stroke=1, my.stroke.color=NA, my.size=1, my.alpha=1,
   showGuide=F, coord.fixed=T, ...
   ){

      pseudo.gradient=T
      if(is.null(plot.title)) plot.title=ea.term
      term_df <- ea.tab %>% dplyr::filter(pathway==ea.term)
      if(nrow(term_df)<1) stop()

      # remove if duplicated terms, that is if both upreg AND downreg for one comparison - keep the most significant
      term_df <- term_df %>% arrange(padj) %>%
         dplyr::filter(!duplicated(paste(ontology, contrast, pathway)))

      ## Here
         # coef1 - Y / dim 1
         # coef2 - X / dim 2
         # p.adjust - score//color//annotation
         # Cluster - will shift sign for p.adjust to tranlate to palette
         # GeneRatio // possibly size of circles
      #term_df$size <- sapply(term_df$GeneRatio, function(x) eval(parse(text=x)))
      if(!is.null(force.levels.c1)){
         term_df$coef1 <- factor(term_df$coef1, levels=force.levels.c1)
      }
      if(!is.null(force.levels.c2)){
         term_df$coef2<- factor(term_df$coef2, levels=force.levels.c2)
      }

      if(!is.factor(term_df$coef1)) term_df$coef1 <- factor(term_df$coef1)
      if(!is.factor(term_df$coef2)) term_df$coef2 <- factor(term_df$coef2)
      if(reverse.y) term_df <- term_df %>% mutate(coef1 = factor(coef1, levels=rev(levels(coef1))))
      if(reverse.x) term_df <- term_df %>% mutate(coef2 = factor(coef2, levels=rev(levels(coef2))))

      term_df <- term_df %>%
         mutate(dim1=as.numeric(coef1), dim2=as.numeric(coef2), annot=NES) %>%
         #mutate(annot= (-1)*log10(annot)) %>%
         #mutate(annot=if_else(regulation=="Dnreg", annot*-1, annot)) %>%
         dplyr::select(dim1, dim2, annot, size, coef1, coef2, padj, NES)

      if(drop.levels) term_df <- droplevels(term_df)

      ## gradient colors is defined here by creating a manual palette for each row
      term_df$annot_bin <- paste0("annot",1:nrow(term_df))
      gradient_sat <- c(-1*gradient.sat, gradient.sat)
      color_key_pseudo_gradient <- dlfoo2::heatcol_creator2(term_df$annot, my.pal = rev(gradient.pal), palette.sat = gradient_sat)$col
      names(color_key_pseudo_gradient) <- term_df$annot_bin

      my_df <- term_df %>%
         arrange(dim2, dim1)

      x_labs <- levels(my_df$coef2)
      y_labs <- levels(my_df$coef1)
      x.n <- length(x_labs)
      y.n <- length(y_labs)

      # START EA CoOc PLOT
      ## ................
         message("... plotting")

         g <- ggplot(my_df) +
               aes(x=dim2, y=dim1)
         #  if plotly text =  paste(coef1, " vs ", coef2, "\nNES: ", signif(NES, 2), "\npadj",signif(padj,2) )

         if(pseudo.gradient) g <- g +
                  geom_point(aes(color=annot_bin), size=my.size, alpha=my.alpha) +
                  scale_color_manual(values=color_key_pseudo_gradient, na.value=NA)

         ## plot stuff
         g <- g + ggtitle(label = plot.title) +
                  labs(x=x.title, y=y.title, colour='') +
                  theme(axis.title.y = element_text(size=9)) +
                  theme(axis.title.x = element_text(size=9)) +
                  theme(legend.text = element_text(size=9)) +
                  theme(title = element_text(size=9)) +
                  theme(axis.text.x=element_text(angle = 45, hjust = 1)) +

                  theme(panel.grid.major =  element_line(colour = 'lightsteelblue2'), panel.grid.minor =  element_line(colour = 'mistyrose3')) +
                  theme(panel.background = element_rect(fill = 'mintcream', colour = NA)) +

                  theme(panel.grid.major =  element_line(colour = NA), panel.grid.minor =  element_line(colour = 'black')) +
                  theme(strip.text = element_text(size=9)) +
                  theme(legend.title =  element_text(size=9, face="bold")) +
                  theme(legend.direction = "horizontal", legend.box = "vertical") +
                  # guides(size=guide_legend(title=size.title, direction="horizontal")) +
                  labs(x="", y="" ) +
                  scale_x_continuous(limits=c(0.5, x.n+0.5), breaks=c(1:x.n), labels=dlfoo2::label_wrap_mod(as.character(x_labs)), expand=c(0,0)) +
                  scale_y_continuous(limits=c(0.5, y.n+0.5), breaks=c(1:y.n), labels=dlfoo2::label_wrap_mod(as.character(y_labs)), expand=c(0,0))


        if(showGuide==F) g <- g +
            guides(fill=FALSE, shape=FALSE, color=FALSE) +
            theme(legend.position='none')

         if(coord.fixed) g <- g +
            coord_fixed(ratio = 1)

         return(g)

} #  endfgseaCoocPlotrWrap


#' @description Wrapper for shiny to plot gene llimma pvalues for all pair-wise comparisons as a table - bubble graph
#' @description The coefficients c1 and c2 are dim1 and dim2 (factors - x/y positions)
#' @description p.value is log10 and translated into annot - which then is used for a annot_bin and psedo-gradient
#' @family Limma
#' @family Bubbleblot
#' @family ScatterPlot
#' @family dotplot
#' @family PlotWrap
#' @param limma.res.tab results table from enrichment analysis. Must contaoin c1 and c2 columns (as factors for choosing plot order)
#' @param gene.id ENSG of the gene
#' @param gradient.pal palette for gradient (p.value) coloring
#' @param gradient.sat
#' @param reverse.y if to reverse y axis order
#' @param reverse.x if to reverse y axis order
#' @param force.levels.c1 if to force order on c1 axis
#' @param force.levels.c2 if to force order on c2 axis
#' @param drop.levels if to drop levels from c1/c2 when these are factors (if not to plot empty resutls)
#' @param x.title
#' @param y.title
#' @param plot.title
#' @param my.stroke
#' @param my.stroke.color
#' @param my.size
#' @param my.alpha
#' @param showGuide defaults to F
#' @param coord.fixed defaults to T
#' @param ... parameters passed on
#' @return a ggplot object
#' @export
limmaGeneCoOcPlot <- function(
   limma.res.tab=NULL,
   gene.id=NULL,
   gradient.pal = c(rev(RColorBrewer::brewer.pal(9,"YlOrRd")[1:9]), rep("white",3), RColorBrewer::brewer.pal(9,"BuGn")[1:9]),
   gradient.sat = 5,
   plot.title=NULL,
   # border.col = NA
   reverse.y=F, reverse.x=F,
   force.levels.c1 = NULL,
   force.levels.c2 = NULL,
   drop.levels=F,
   x.title=NULL, y.title=NULL,
   my.stroke=1, my.stroke.color=NA, my.size=1, my.alpha=1,
   showGuide=F, coord.fixed=T, ...
   ){

      pseudo.gradient=T
      if(is.null(plot.title)) plot.title=gene.id
      limma_df <- limma.res.tab %>%
         dplyr::filter(ENSG==gene.id & gene_padj_BH_call!="ns")

      if(nrow(limma_df)<1) return(ggPlotEmpty())


      limma_df <- limma_df %>% arrange(p.value)


      ## Here
         # c1 - Y / dim 1
         # c2 - X / dim 2
         # plog - score//color//annotation
         # gene_padj_BH_call - will shift sign for p.adjust to tranlate to palette


      if(!is.null(force.levels.c1)){
         limma_df$c1 <- factor(limma_df$c1, levels=force.levels.c1)
      }
      if(!is.null(force.levels.c2)){
         limma_df$c2 <- factor(limma_df$c2, levels=force.levels.c2)
      }

      if(!is.factor(limma_df$c1)) limma_df$c1 <- factor(limma_df$c1)
      if(!is.factor(limma_df$c2)) limma_df$c2 <- factor(limma_df$c2)
      if(reverse.y) limma_df <- limma_df %>% mutate(c1 = factor(c1, levels=rev(levels(c1))))
      if(reverse.x) limma_df <- limma_df %>% mutate(c2 = factor(c2, levels=rev(levels(c2))))

      limma_df <- limma_df %>%
         mutate(dim1=as.numeric(c1), dim2=as.numeric(c2)) %>%
         dplyr::mutate(plog = -1*log10(p.value)) %>%
         mutate(annot = if_else(gene_padj_BH_call=="Downreg", plog*-1, plog)) %>%
         dplyr::select(dim1, dim2, annot,  c1, c2, p.value, plog, lfc, gene_padj_BH_call)

      if(drop.levels) limma_df <- droplevels(limma_df)

      ## gradient colors is defined here by creating a manual palette for each row
      limma_df$annot_bin <- paste0("annot",1:nrow(limma_df))
      gradient_sat <- c(-1*gradient.sat, gradient.sat)
      color_key_pseudo_gradient <- dlfoo2::heatcol_creator2(limma_df$annot, my.pal = rev(gradient.pal), palette.sat = gradient_sat)$col
      names(color_key_pseudo_gradient) <- limma_df$annot_bin

      my_df <- limma_df %>%
         arrange(dim2, dim1)

      x_labs <- levels(my_df$c2)
      y_labs <- levels(my_df$c1)
      x.n <- length(x_labs)
      y.n <- length(y_labs)

      # START EA CoOc PLOT
      ## ................
         message("... plotting")

         g <- ggplot(my_df) +
               aes(x=dim2, y=dim1)
         #  if plotly text =  paste(c1, " vs ", c2, "\nNES: ", signif(NES, 2), "\np.value",signif(p.value,2) )

         if(pseudo.gradient) g <- g +
                  geom_point(aes(color=annot_bin), size=my.size, alpha=my.alpha) +
                  scale_color_manual(values=color_key_pseudo_gradient, na.value=NA)

         ## plot stuff
         g <- g + ggtitle(label = plot.title) +
                  labs(x=x.title, y=y.title, colour='') +
                  theme(axis.title.y = element_text(size=9)) +
                  theme(axis.title.x = element_text(size=9)) +
                  theme(legend.text = element_text(size=9)) +
                  theme(title = element_text(size=9)) +
                  theme(axis.text.x=element_text(angle = 45, hjust = 1)) +

                  theme(panel.grid.major =  element_line(colour = 'lightsteelblue2'), panel.grid.minor =  element_line(colour = 'mistyrose3')) +
                  theme(panel.background = element_rect(fill = 'mintcream', colour = NA)) +

                  theme(panel.grid.major =  element_line(colour = NA), panel.grid.minor =  element_line(colour = 'black')) +
                  theme(strip.text = element_text(size=9)) +
                  theme(legend.title =  element_text(size=9, face="bold")) +
                  theme(legend.direction = "horizontal", legend.box = "vertical") +
                  # guides(size=guide_legend(title=size.title, direction="horizontal")) +
                  labs(x="", y="" ) +
                  scale_x_continuous(limits=c(0.5, x.n+0.5), breaks=c(1:x.n), labels=dlfoo2::label_wrap_mod(as.character(x_labs)), expand=c(0,0)) +
                  scale_y_continuous(limits=c(0.5, y.n+0.5), breaks=c(1:y.n), labels=dlfoo2::label_wrap_mod(as.character(y_labs)), expand=c(0,0))


        if(showGuide==F) g <- g +
            guides(fill=FALSE, shape=FALSE, color=FALSE) +
            theme(legend.position='none')

         if(coord.fixed) g <- g +
            coord_fixed(ratio = 1)

         return(g)

} #  limmaGeneCoOcPlot





#' Enrichemnt plot in gsea-line format. One line/vertical bar for each gene
#' Genes are provided as a list of named vectors vit values (different groups)
#' gene values are given a rank and plotted (x-axis). on y-axis the different groups are ordered (acornding to mean rank)
#' rank lines are colored and according to color.key provided
#' @family plot wraps
#' @family gsea
#' @param x Data frame with, at least:  one 'group', one 'value', and one 'gene_id' column
#' @param gene.id.column specify name of gene.id column
#' @param value.column specify name of group column
#' @param group.column specify name of value column - used to create RANKS
#' @param n_genes the max rank (x-value)
#' @param rank.descending TRUE if leftmost should be the highest x value
#' @param pathway.genes Vector of genes that belong to a specific pathway to plot
#' @param pathway.name I to add statistics
#' @param reverse.y if to reverse y axis
#' @param color.gradient.on.rank if to color each group (row) on a group color a gradient based on the rank value. Defaults to TRUE
#' @param color.key color key for group colors if gradient is FALSE
#' @param ggplotply.prepare if output is for ggplotly, then set this to TRUE to get a ggploty gradient (ggplotly uses fill_ instead of color_)
#' @return
#' @export
geneRankPlot <- function(
   x,
   gene.id.column = "SYMBOL",
   group.column = "contrast",
   value.column = "t.value",
   rank.descending = T,
   n_genes=NULL,
   pathway.genes=NULL,
   pathway.name=NULL,
   reverse.y=T,
   color.gradient.on.rank = T,
   ggplotply.prepare = F,
   color.key=dlfoo2::color_annotations, ...
   ){

   if(is.null(n_genes)) stop("You must specify n_genes - this sets length of x-axis")
   if(is.null(pathway.genes)) stop("You must specify pathway.genes - character vecor of gene identifiers")
   #if(!("sample_id" %in% colnames(x)) stop("each entity must have a sample_id dummy name")
    plot.title = pathway.name

    if(ggplotply.prepare==F) message("\n ... note!!! plot colors not optimized for ggplotly - if so set to 'ggplotply.prepare=T'")


   # POINT SIZE and SHAPE:: create wighted sizes for the points - (alos arrange accornding to this)
   # n_genes <- length(rank.list[[1]])

   my_df <- as_tibble(x) %>%
      mutate(gene_id = UQ(rlang::sym(gene.id.column))) %>%
      mutate(plot.group = UQ(rlang::sym(group.column))) %>%
      mutate(rank.value = UQ(rlang::sym(value.column))) %>%
      dplyr::mutate(orig.value = rank.value) %>%
      group_by(plot.group) %>%
      dplyr::mutate(rank.value = rank(-1*rank.value)) %>%
      mutate(point_size = abs(n_genes/2 - rank.value)) %>%
      mutate(point_size = (point_size/(n_genes/2))^3) %>%
      #mutate(point_alpha = (point_size/(n_genes/2))) %>%
      mutate(point_alpha  = 1) %>%
      mutate(point_shape = if_else(rank.value <= n_genes/2, "up",  "dn")) %>%
      arrange(point_size) %>%
      dplyr::filter(gene_id %in% pathway.genes) %>%
      droplevels()



   ##  Rank the grupus on mean orig value
   groups_ranked <- my_df %>% group_by(plot.group) %>% summarise(val_mean = mean(orig.value))
   groups_ranked <- unique(groups_ranked$plot.group)[order(groups_ranked$val_mean, decreasing = T)]

   my_df <- data.frame(my_df) %>% mutate(plot.group=factor(plot.group, levels=groups_ranked))
   if(reverse.y) my_df <- my_df %>% mutate(plot.group = factor(plot.group, levels=rev(levels(plot.group))))



   ## Color.key - IF Pseudo gradient
   ## ----------------------
      # fix for problems with gradient_fill in plotly - define a pseudigradient instaed
      if(color.gradient.on.rank==T){
         message(" ... setting upp pseudo-gradient - manual_fill")
         # if(gradient.normalize) my_df$annotation <- my_df$annotation - mean(my_df$annotation, na.rm=T) / sd(my_df$annotation, na.rm=T)
         ## gradient colors is defined here by creating a manual palette for each row
         palette.gradient <- rev(c(rev(RColorBrewer::brewer.pal(9,"YlOrRd")[c(3,3,3:7,7,8,9)]), rep("#FFEDA0",5), rep("white",3), rep("#DEEBF7",5), RColorBrewer::brewer.pal(9,"Blues")[c(3,3,3:7,7,8,9)]))

         my_df$value_bin <- paste0("annotation", 1:nrow(my_df))
         # color_key_pseudo_gradient
         # my_df <- my_df %>% dplyr::rename(annotation_num = annotation, annotation = annotation_bin)
         color.key <- dlfoo2::heatcol_creator2(my_df$rank.value, my.pal = rev(palette.gradient), palette.sat = c(1,n_genes))$col
         names(color.key) <- my_df$value_bin
         color_key_x <-  color.key
         } # end pseudo.gradient


      ## color.key & shape keys (named vectors w colors)
      ## ---------
    if(color.gradient.on.rank==F){
      color_key_x <-  colorKeyFix(color.key, as.character(levels(my_df$plot.group)))
    }
      shape_key_x <- c(up=3, dn=3)


   ##    Plot ggplot
   ## :::::::::::::::::::
      # if(is.null(plot.title)) plot.title = pdata.column

   if(!color.gradient.on.rank & !ggplotply.prepare ) {
      g <- ggplot(my_df) +
         aes(y=plot.group, x=rank.value, size=point_size, shape=point_shape, color=as.character(plot.group), fill=as.character(plot.group))   +
         #coord_flip() +
         #geom_violin(fill="white", width=1, trim = TRUE, scale = "width", adjust = 0.5) +
         coord_cartesian(xlim = c(0,n_genes)) +
         geom_point(stroke=0.5, alpha=0.8) +
         guides(fill=FALSE, shape=FALSE, color=FALSE) + theme(legend.position='none') +
         scale_color_manual(values=color_key_x, na.value=NA) +
         scale_shape_manual(values=shape_key_x)}

    if(!color.gradient.on.rank & ggplotply.prepare ) {
      g <- ggplot(my_df) +
         aes(y=plot.group, x=rank.value, size=point_size, shape=point_shape, color=as.character(plot.group), fill=as.character(plot.group))   +
         #coord_flip() +
         #geom_violin(fill="white", width=1, trim = TRUE, scale = "width", adjust = 0.5) +
         coord_cartesian(xlim = c(0,n_genes)) +
         geom_point(stroke=0.5, alpha=0.8) +
         guides(fill=FALSE, shape=FALSE, color=FALSE) + theme(legend.position='none') +
         scale_fill_manual(values=color_key_x, na.value=NA) +
         scale_shape_manual(values=shape_key_x)
      }


   if(color.gradient.on.rank & ggplotply.prepare) {
      message("\n .. plotting pseudo graient - ggplotply.prepare option TRUE")
      g <- ggplot(my_df) +
         aes(y=plot.group, x=rank.value, size=point_size, shape=point_shape, fill=as.character(value_bin), color=as.character(value_bin))   +
         #coord_flip() +
         #geom_violin(fill="white", width=1, trim = TRUE, scale = "width", adjust = 0.5) +
         coord_cartesian(xlim = c(0,n_genes)) +
         geom_point(aes(fill=plot.group), stroke=0.5, alpha=0.8) +
         guides(fill=FALSE, shape=FALSE, color=FALSE) + theme(legend.position='none') +
         scale_fill_manual(values=color_key_x, na.value=NA) +
         scale_shape_manual(values=shape_key_x)
      }

   if(color.gradient.on.rank & !ggplotply.prepare) {
      message("\n .. plotting pseudo graient - ggplotply.prepare option FALSE")
      g <- ggplot(my_df) +
         aes(y=plot.group, x=rank.value, size=point_size, shape=point_shape, color=value_bin, fill=value_bin)   +
         #coord_flip() +
         #geom_violin(fill="white", width=1, trim = TRUE, scale = "width", adjust = 0.5) +
         coord_cartesian(xlim = c(0,n_genes)) +
         geom_point(aes(fill=plot.group), stroke=0.5, alpha=0.8) +
         guides(fill=FALSE, shape=FALSE, color=FALSE) + theme(legend.position='none') +
         scale_color_manual(values=color_key_x, na.value=NA) +
         scale_shape_manual(values=shape_key_x)
      }


   ## plot stuff
   g <- g + ggtitle(label = plot.title) +
            labs(x="gene rank",  colour='') +
            theme(axis.title.y = element_text(size=7)) +
            theme(axis.title.x = element_text(size=9)) +
            theme(legend.text = element_text(size=9)) +
            theme(title = element_text(size=9))
            #theme(axis.text.x=element_text(angle = 45, hjust = 1))

   # Warning: The 'plotly_click' event tied a source ID of 'A' is not registered. In order to obtain this event data,
   # please add `event_register(p, 'plotly_click')` to the plot (`p`) that you wish to obtain event data from.
   # ggplotly(g)

   return(g)


} # end  geneRankPlot










      # plot.name = "ChordPlot_ReTax_ConsensusClustering_Methylation_k25",
      # plot.dir = "~/PROJECTS/RCC_ccpRCC_2019/reTaxonomy/Figures",
      # cc.k = 25,
      # cc.dir = "~/PROJECTS/RCC_ccpRCC_2019/ConsensusClust/ConsensusCluster_rcc_t_me_589set_2065genes_reps.1000_maxK.25_innerLinkage.ward_sd.0.01/"
      #




   # source("~/PROJECTS/RCC_ccpRCC_2019//R_SCRIPTS/R_ccpRCC_2019_SOURCE.R")
  # cc_key <- unlist(cc_grouping_list$cc_gex$cc_groups_col)
  #  names(cc_key) <- unlist(cc_grouping_list$cc_gex$cc_groups)
  #  color_key <- c(color_subtypes, cc_key)
  #
  #
  #     plot.name = "ChordPlot_ReTax_ConsensusClustering_GEX_k25_cc880set",
  #     plot.dir = "~/PROJECTS/RCC_ccpRCC_2019/reTaxonomy/ChordPlots_reTaxonomy/",
  #     cc.k = 25,
  #     cc.dir = "~/PROJECTS/RCC_ccpRCC_2019/ConsensusClust/ConsensusCluster_rcc_t_5000_880set_5000genes_reps.1000_maxK.25_innerLinkage.ward_sd.rank.5000/",
  #     dend.color.key = cc_key,
  #     dend.rev.clusters <- c(rep(F,24),T)
  #

#' chordPlot_ccHCA_to_reTax
#' @description  specialized function alluvial chordPlot re-taxonomy w consensus HCA
#' @family chordplot
#' @family PlotWrap
#' @family circlize
#' @param plot.name name of file to plot
#' @param plot.dir directory to plot in
#' @param cc.dir directory in which consensus clust results are - parsed using dlfoo2::ConsensusCluster_process
#' @param cc.k number of consensusclusters (cannot be higher than performed in cc.dir)
#' @param group2_levels factor levels for the re-taxonomy. specifies order of groups in plot
#' @param group_names Title of plot
#' @param scale_factors what color key to use for annotation. Detaults to dlfoo2::color_subtypes (but could be set to e.g. dlfoo2::color_annotations)
#' @param do.gradient set to FALSE if to manually override numeric check for gradient or manual fil
#' @param color_key what color key for subgroups.
#' @param dend.color.key if a specific color key is to be used for dendogram
#' @param dend.rev.clusters if to flip clusters.
#' @return plot pdf to file
#' @export

chordPlot_ccHCA_to_reTax <- function(
   plot.name = "temp",
   plot.dir = "~/PROJECTS/RCC_ccpRCC_2019/reTaxonomy/Figures/",
   cc.dir = "~/PROJECTS/RCC_ccpRCC_2019/ConsensusClust/ConsensusCluster_rcc_t_me_589set_2065genes_reps.1000_maxK.25_innerLinkage.ward_sd.0.01/",
   cc.k = 25,
   group2_levels = c("outlier","metanephric",  "mtscRCC", "dCIMP", "mesRCC", "ccpRCC","chONC", "pONC"), # "Uro"
   group_names = c("cc","retax"),
   scale_factors = c(1, 3),
   color_key=NULL,
   dend.rev.clusters = NULL,
   dend.color.key = NULL
){



   # source("~/PROJECTS/RCC_ccpRCC_2019//R_SCRIPTS/R_ccpRCC_2019_SOURCE.R")

   require(reshape2)
   require(dplyr)
   options(dplyr.width = Inf)
   require(tibble)
   require(tidyverse)
   require(tidyr)
   require(circlize)
   require(dendextend)

   if(is.null(color_key)) color_key <- color_subtypes
   if(is.null(dend.color.key)) dend.color.key <- color_annotations
   color_key <- c(color_key, dend.color.key)


   ## DENDROGRAM
   ## :::::::::::::
   message(".. fetching consensusCluster results")
   my_cc <- ConsensusCluster_process(cc.dir, k=cc.k)
      # identical(pdata$sample_id, my_cc$sample.names)
      #pdata <- dlfoo2data::pdata_tcgaPanCan %>% dplyr::filter(platform == "HumanMethylation450") %>%  dplyr::filter(cancer.type %in% c("KIRC","KIRP","KICH"))
   pdata <- dlfoo2::pdata_panCanFull
   u <- match(my_cc$sample.names, pdata$sample_id)
   pdata <-  pdata[u, ]
   pdata <- pdata %>% mutate(cc = paste0("cc_",my_cc$cc))
   pdata <- droplevels(pdata)
   # table(pdata$taxonomy_published)
   # table(pdata$tax_simp)
   # table(pdata$cc, pdata$tax_simp)

   if(!all(group2_levels %in% unique(pdata$tax_simp))){
      message("... not all retaxonomy levels present. correcting group2_levels vector")
      group2_levels <- group2_levels[group2_levels %in% unique(pdata$tax_simp)]
      message("... ...",paste(group2_levels, cat=", "))
   }


   # str(pdata) # 589 obs. of  70 variables:
   message("... ",nrow(pdata)," samples")  # 582
   message("... ",sum(table(pdata$taxonomy_published))," samples annotated in taxonomy published")  # 582
   message("... ",sum(table(pdata$tax_simp))," samples annotated in re-taxonomy simplified")  # 582
   # sum(table(pdata$tax_simp)) # 584


   # GROUPS
   # :::::::::

   group1_levels <- paste0("cc_",unique(my_cc$cc[my_cc$i]))
   group1_levels <- factor(group1_levels, levels=group1_levels)

   dend <- as.dendrogram(my_cc$consensus.tree)
   dend <- dend %>%
         dendextend::set("branches_lwd", value = 1.5) %>%
         dendextend::set("branches_k_color", k=cc.k, value= dend.color.key[as.character(group1_levels)])

   ## flip clusters
   if(!is.null(dend.rev.clusters)){
      stopifnot(length(dend.rev.clusters)==cc.k)
      dend <- dlfoo2::dend_rev_clusters(dend, rev.clusters = dend.rev.clusters)
      # re-structure my_cc object since sample order is changed
         # my_cc2 <- my_cc
         new_i <- match(labels(dend), my_cc$sample.names)
         my_cc$i <- new_i
   }
   # plot(dend %>% set())


   dend_list <- list(cc=dend, retax=NULL)
      #group2_levels <- c("outlier","metanephric",  "mtscRCC", "dCIMP", "mesRCC", "ccpRCC","chONC", "pONC") # "Uro"
      #group_names <- c("cc","retax")
   group_levels <- list(
         cc=group1_levels,
         retax=group2_levels
      )
      #scale_factors <- c(1, 3)


   # SUB-SEGMENTS
   # ::::::::::::
   seg_list <- list()
   df1 <- pdata %>%
      mutate(cc = paste0("cc_",my_cc$cc)) %>%
      dplyr::rename(sector=cc) %>%
      dplyr::select(sector) %>%
      mutate(sector = factor(sector, levels=group1_levels)) %>%
      group_by(sector) %>%
      arrange(sector) %>%
      summarise(n=n()) %>%
      mutate(group=factor("cc"))
   df2 <-  pdata %>%
      #dplyr::filter(!is.na(taxonomy_published)) %>%
      dplyr::rename(sector=tax_simp) %>%
      dplyr::filter(sector %in% group2_levels) %>%
      #mutate_all(.funs = as.character) %>%
      dplyr::select(sector) %>%
      mutate(sector = factor(sector, levels=group2_levels)) %>%
      group_by(sector) %>%
      arrange(sector) %>%
      summarise(n=n()*scale_factors[2])  %>%
      mutate(group=factor("retax"))


      ##  ANNOTATIONS track 1 - (cc-rectangles)
      ## ---------------------
      seg_list$track1$segments <- list(
            cc=as.data.frame(df1),
            retax=as.data.frame(df2)
            )
      seg_list$track1$colors <- lapply(seg_list$track1$segments, function(x){
         if(!is.null(x)) unlist(apply(x, 1, function(y) {rep(color_key[as.character(y[1])], y[2])}))
      })



      ##  ANNOTATIONS track 2 - (color_vec/color_tab)
      ## ---------------------
      ## re-ordered from pdata!!
      annot_levels <- c("CC-e.2",   "CC-e.1",   "CC-e.3", "mixed", "Ch-e" ,  "P-e.1a",   "P-e.1b"   ,"P-e.2", "P.CIMP-e")
      df3 <- pdata %>%
            rename(annot = taxonomy_published) %>%
            #dplyr::filter(!is.na(annot))
            dplyr::select(annot) %>%
            mutate(annot = factor(annot, levels=annot_levels)) %>%
            group_by(annot) %>%
            arrange(annot) %>%
            summarise(n=n()) %>%
            mutate(group=factor("cc"))


      my_df <- data.frame(pdata %>% dplyr::select(taxonomy_published), row.names = pdata$sample_id)
      seg_list$track2$colors$cc <- sapply(my_df[labels(dend_list[["cc"]]),], function(x){
          u <- match(as.character(x), names(color_key))
          return(color_subtypes[u])
       })

      seg_list$track2$segments <- list(
         cc=as.data.frame(df3),
         #cc=rectMatrix(names(seg_list$track2$colors$cc))
         retax=as.data.frame(df2)
         )

      seg_list$track2$colors$retax <- unlist(
         apply(
            seg_list$track2$segments$retax, 1, function(y) {
               rep(color_key[as.character(y[1])], y[2])
               }))


   ##    X-LIMS
   ## :::::::::::
   xlim_list <- list()
   xlim_list$track1 <- lapply(seg_list$track1$segments, function(x) rectMatrix(x$n))
   xlim_list$track2 <- lapply(seg_list$track2$segments, function(x) rectMatrix(x$n))



   ## LINKS
   ## ::::::::::
   links_df <- pdata %>%
      mutate(cc = paste0("cc_",my_cc$cc)) %>%
      arrange(order(match(my_cc$sample.names[my_cc$i], pdata$sample_id))) %>%
      #dplyr::rename(s1=cc, s2=tax_simp) %>%
      mutate(s1=1:nrow(pdata)) %>%
      mutate(s2 = tax_simp) %>%
      dplyr::select(s1, s2, tax_simp) %>%
      mutate_all(.funs = as.character) %>%
      mutate(n=1)
      #dplyr::filter(s2%in%group2_levels)

   links_df <- links2sectors(links_df, scale.x2 = scale_factors[2])
   # links_df <- links_df[!apply(links_df, 1, function(x) x[1]==x[2]),]
   links_df <- links_df[links_df$s2 %in% group2_levels,]
   links_df$s1x1 <- links_df$order-1
   links_df$s1x2 <- links_df$order
   links_df$taxonomy_published <- links_df$s1
   links_df$tax_simp <- links_df$s2
   u <- sapply(links_df$s2, function(x) {
      match(x, seg_list$track2$segments$retax$sector) })
   uu <- as.numeric(xlim_list$track2$retax[,1])[u]
   links_df$s2x1 <- links_df$s2x1+uu
   links_df$s2x2 <- links_df$s2x2+uu
   links_df$s1 <- "cc"
   links_df$s2 <- "retax"




   ## PLOT TO FILE
   ## ::::::::::::::::
   my.file <- paste0(file.path(plot.dir, plot.name),".pdf")
   pdf(my.file, width = 12, height = 12, useDingbats = F)

   ## START PAR
   circos.par(cell.padding = c(0, 0, 0, 0), gap.degree = 15, start.degree = (-45),  "track.height" = 0.2)
   circos.initialize(factors = group_names, xlim = matrix(unlist(lapply(xlim_list$track1, range)), nrow=2, byrow=T))

      # Title
      title(plot.name)

      ## add DENDOGRAM
      max_height = attr(dend_list[["cc"]], "height")

      circos.track(ylim = c(0, max_height), bg.border = NA, track.height = 0.3,
       panel.fun = function(x, y) {
        sector.index = get.cell.meta.data("sector.index")
           dend = dend_list[[sector.index]]
           if(sector.index=="cc") circos.dendrogram(dend, max_height = max_height, facing = "inside")
      })



      ## add HCA cluster rects (cc)
      circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.1,
         panel.fun = function(x, y){
            xlim = CELL_META$xlim
            ylim = CELL_META$ylim
            m = xlim_list$track1[[CELL_META$sector.index]]

            if(CELL_META$sector.index=="cc") circos.rect(xleft = m[,1], ybottom = 0, xright = m[,2], ytop = 1,
                     col = color_key[as.character((seg_list$track1$segments[[CELL_META$sector.index]]$sector))],

                     border = "black")
            if(CELL_META$sector.index=="cc") circos.text(x = m[,3], y = 0.5, labels = as.character(seg_list$track1$segments[[CELL_META$sector.index]]$sector), facing = "bending.inside", cex=0.5, col="white")
            if(CELL_META$sector.index=="retax") circos.rect(xleft = m[,1], ybottom = 0, xright = m[,2], ytop = 1,
                     col = NA, border = NA)
            })


      ## Rectangles track2
      circos.track(ylim = c(0,1), bg.border = NA, track.height=0.175, panel.fun = function(x, y) {
         l = CELL_META$sector.index
         if(l=="cc"){
            m = xlim_list$track2[[CELL_META$sector.index]]
            circos.rect(xleft = c(0:(max(xlim_list$track2$cc)-1)), ybottom = 0, xright = 1:max(xlim_list$track2$cc), ytop = 1,
                  col = seg_list$track2$colors[["cc"]], border = NA)}
         if(l=="retax"){
            m = xlim_list$track2$retax
            circos.rect(xleft = m[,1], ybottom = 0, xright = m[,2], ytop = 1,
                  col = unique(seg_list$track2$colors[[l]]), border = "black")
                  circos.text(x = m[,3], y = 0.5, labels = seg_list$track2$segments[["retax"]]$sector, facing = "clockwise", cex=0.6, col="black")
                  circos.text(x = m[,3], y = 0, labels = as.character(seg_list$track2$segments[["retax"]]$n/scale_factors[2]), facing = "inside", cex=0.5, col="black", adj=c(0,1.5))
                  }
         })

          ## Links
      for(i in 1:nrow(links_df)){
         #links_df
         circos.link(
            sector.index1 = links_df$s1[i],
            point1 = c(links_df$s1x1[i], links_df$s1x2[i]),
            sector.index2 = links_df$s2[i],
            point2 = c(links_df$s2x1[i], links_df$s2x2[i]),
            col= color_subtypes[links_df$tax_simp[i]])
      }

      circos.clear()

     dev.off()
} # end function chord plot consensus clusters





#' chordPlot_ccHCA_to_reTax - Plot HCA tree for Consensus Clustering
#' Here one can shoose to plot with Chord iagram or not.
#' Choose what annotations and track heights usong annotation.columns
#' @description  specialized function alluvial chordPlot re-taxonomy w consensus HCA
#' @family chordplot
#' @family PlotWrap
#' @family circlize
#' @param pdata what pdata to use for plot
#' @param annotation.columns what pdata columns to plot as annotation segments
#' @param predefined.cc.column if the loaded cc already has been put in a specific pdata colum, e.g. cc_gex
#' @param chord.column name of pdata column if to plot add a second grup, pie section, to plot a chord plot.
#' @param chord.levels factor levels for the re-taxonomy. specifies order of groups in plot
#' @param scale_factors if to scale size of the chord part
#' @param cc.dir directory in which consensus clust results are - parsed using dlfoo2::ConsensusCluster_process
#' @param cc.k number of consensusclusters (cannot be higher than performed in cc.dir)
#' @param plot.name name of file to plot
#' @param plot.dir directory to plot in
#' @param do.gradient set to FALSE if to manually override numeric check for gradient or manual fil
#' @param color.key what color key for subgroups.
#' @param dend.color.key if a specific color key is to be used for dendogram
#' @param dend.rev.clusters if to flip clusters.
#' @param annotation.track.height what plot heights for the annotation tracks. integer vector of same length as annotation.columns
#' @return plot pdf to file
#' @export

circosPlot_cc_hca_tree <- function(
   pdata = NULL,
   #group.names = c("cc","hdbcl","taxonomy")
   annotation.columns = c("hdbcl_gex", "taxonomy_published"),
   predefined.cc.column = "cc_gex",
   plot.name = "ChordPlot_GEX_cc_hdb_chen",
   plot.dir = "~/PROJECTS/RCC_ccpRCC_2019/reTaxonomy/2_Evaluate_cc_hdb/CircosPlots/",
   cc.dir = "~/PROJECTS/RCC_ccpRCC_2019/ConsensusClust/ConsensusCluster_rcc_t_880set_3177genes_reps.1000_maxK.25_innerLinkage.ward.D2_sd.1/",
   cc.k = 25,
   chord.column = NULL,
   chord.levels = c("outlier","metanephric",  "mtscRCC", "dCIMP", "mesRCC", "ccpRCC","chONC", "pONC"), # "Uro"
   scale_factors = c(1, 3),
   color.key=NULL,
   dend.rev.clusters = NULL,
   dend.color.key = NULL,
   annotation.track.height = NULL
){

   ## Groups is if to splir the cicle into >1 parts.
   ## Sectrors are the pie slices of each group.
   ## Tracks are the different rows of the graph (and per sector and group)
   # Here

   require(reshape2)
   require(dplyr)
   options(dplyr.width = Inf)
   require(tibble)
   require(tidyverse)
   require(tidyr)
   require(circlize)
   require(dendextend)

   if(is.null(color.key)) color.key <- color_subtypes
   if(is.null(dend.color.key)) dend.color.key <- color_annotations
   color_key <- c(color.key, dend.color.key, color_subtypes)
   annot.i <- match(annotation.columns, colnames(pdata))
   if(any(is.na(annot.i))) stop("not all annotation columns in pdata")
   if(is.null(annotation.track.height)) annotation.track.height <- rep(0.15, length(annotation.columns))
   if(!is.null(predefined.cc.column)) {
         u <- match(predefined.cc.column, colnames(pdata))
         if(any(is.na(u))) stop("predefined.cc.column not present in pdata chord column")}

   ## HCA tree: Fetch dendrogram & create my.k annatations
   ## :::::::::::::::::::::::::::
      ## Track 1 is here always "cc" which is set to "main"
      ## If a pre-defined cc track is supplied (predefined.cc.column) then this is set as the main
   message(".. fetching consensusCluster results")
   my_cc <- ConsensusCluster_process(cc.dir, k=cc.k)
   u <- match(my_cc$sample.names, pdata$sample_id)
   if(any(is.na(u))) stop("not all consensus cluster names in pdata")
   pdata <-  pdata[u, ]

   ## create cc-column for main cc annotation. if predefined.cc.column is supplied, rename this to 'cc'
   if(is.null(predefined.cc.column)) pdata <- pdata %>% mutate(cc = factor(paste0("cc_",my_cc$cc)))
   if(!is.null(predefined.cc.column)) pdata <- pdata %>% dplyr::rename(cc = !!predefined.cc.column)
   pdata <- droplevels(pdata)
   message("... ",nrow(pdata)," samples")  # 582



   # GROUPS
   # :::::::::
      # groups here define the two different plot parts of thegraph, i.e. a HCA/CC/TAX and a TAX part used for chord diagram
      # only extended and used if chord column i specified
      # Groups correspond to each Pie-section of the circle. Multiple groups only if to draw e.g. Chord diagram

   group_names <- c("main")
   if(!is.null(chord.column)){
      if(is.na(match(chord.column, colnames(pdata)))) stop ("chord.column not found in pdata")
      group_names <- c("main","chord")
      if(!is.null(chord.levels))
         u <- match(chord.levels, unique(as.character(pdata %>% dplyr::pull(chord.column))))
         if(any(is.na(u))) stop("not all chord.levels present in pdata chord column")
      }

   ## note for dend colors - cc has to be ordered as factor to get colors appear proper order
   dend <- as.dendrogram(my_cc$consensus.tree)
   dend <- dend %>%
         dendextend::set("branches_lwd", value = 1.5) %>%
         dendextend::set("branches_k_color", k=cc.k, value = color_key[as.character(levels(pdata$cc))])

   ## if rev/flip clusters !! not sure how this reflects of annotations ... has to be double checked
   if(!is.null(dend.rev.clusters)){
      stopifnot(length(dend.rev.clusters)==cc.k)
      dend <- dlfoo2::dend_rev_clusters(dend, rev.clusters = dend.rev.clusters)
      new_i <- match(labels(dend), my_cc$sample.names)
      my_cc$i <- new_i
      u <- match(my_cc$sample.names, pdata$sample_id)
      pdata <-  pdata[u, ]
   }
   dend_list <- list(main=dend, chord=NULL)



   # SUB-SEGMENTS/ROWS (multiple annotation rows. defined by annotation.columns)
   # ::::::::::::
      # 'Tracks' corresponds to annotations/rows
      # Note!! in the first row (main track) consensus clusters as defined by the dendrogram will be plotted as defined by cc or predefined.cc.column

   group_df_list <- list()

   group_df_list[["main"]] <- pdata %>%
      dplyr::rename(sector=cc) %>%
      dplyr::select(sector) %>%
      group_by(sector) %>%
      arrange(sector) %>%
      summarise(n=n()) %>%
      mutate(group=factor("main"))

   # add df of sectors for second group (chord.column)
   if(!is.null(chord.column)){
      group_df_list[["chord"]] <- pdata %>%
         dplyr::rename(sector:=!!chord.column) %>%
         dplyr::filter(sector %in% chord.levels) %>%
            #mutate_all(.funs = as.character) %>%
         dplyr::select(sector) %>%
         mutate(sector = factor(sector, levels=chord.levels)) %>%
         group_by(sector) %>%
         arrange(sector) %>%
         summarise(n=n()*scale_factors[2])  %>%
         mutate(group=factor(chord.column))
   }



   ##  Tracks/Annotations . color bars/rectangles for consensus clusters
   ## ----------------------------------------
      # step through each row/track and create df for each of the groups
      # track list [[1]] is reserved for the cc's
   track_list <- list()
   track_list[[1]] <- list()
   track_list[[1]]$segments <- lapply(names(group_df_list), function(x){
            as.data.frame(group_df_list[[x]])
            })
   names(track_list[[1]]$segments) <- names(group_df_list)
   track_list[[1]]$colors <- lapply(track_list[[1]]$segments, function(x){
      if(!is.null(x)) unlist(apply(x, 1, function(y) {rep(color_key[as.character(y[1])], y[2])}))
   })



   ## add/loop additional tracks/annotations (color_vec/color_tab)
   ## ---------------------
      # annot_levels <- c("CC-e.2",   "CC-e.1",   "CC-e.3", "mixed", "Ch-e" ,  "P-e.1a",   "P-e.1b"   ,"P-e.2", "P.CIMP-e")
   n_tot_annotations <- length(annotation.columns)+1
   for(i in 1:length(annotation.columns)){
      ii <- i+1
      track_list[[ii]] <- list()
      df_i <- pdata %>%
            rename(annot:= !!annotation.columns[i]) %>%
            #dplyr::filter(!is.na(annot))
            dplyr::select(annot) %>%
            mutate(annot = as.character(annot)) %>%
            group_by(annot) %>%
            arrange(annot) %>%
            summarise(n=n()) %>%
            mutate(group=factor("main"))

      my_df <- data.frame(pdata %>% dplyr::select(!!annotation.columns[i]), row.names = pdata$sample_id)
      track_list[[ii]][["colors"]][["main"]] <- sapply(my_df[labels(dend_list[["main"]]),], function(x){
             u <- match(as.character(x), names(color_key))
             return(color_key[u])
          })
      track_list[[ii]][["segments"]][["main"]] <- as.data.frame(df_i)
      } # end annotation.columns




   ##    X-LIMS
   ## :::::::::::
      #  define coordinates for the different segments/sections in each group and for each track / row
   xlim_list <- vector(mode="list", length=length(track_list))
   for(i in 1:length(track_list)){
      xlim_list[[i]] <- lapply(track_list[[i]][["segments"]], function(x) dlfoo2::rectMatrix(x$n))
   }



   ## LINKS
   ## ::::::::::
   links_df <- data.frame()
   if(!is.null(chord.column)){
      links_df <- pdata %>%
         # mutate(cc = paste0("cc_",my_cc$cc)) %>%
         arrange(order(match(my_cc$sample.names[my_cc$i], pdata$sample_id))) %>%
         #dplyr::rename(s1=cc, s2=tax_simp) %>%
         mutate(s1=1:nrow(pdata)) %>%
         # mutate(s2 = chord.column) %>%
         dplyr::rename(s2 = !!chord.column) %>%
         dplyr::select(s1, s2) %>%
         mutate_all(.funs = as.character) %>%
         mutate(n=1)
         #dplyr::filter(s2%in%chord.levels)

      links_df <- links2sectors(links_df, scale.x2 = scale_factors[2])
      # links_df <- links_df[!apply(links_df, 1, function(x) x[1]==x[2]),]
      links_df <- links_df[links_df$s2 %in% chord.levels,]
      links_df$s1x1 <- links_df$order-1
      links_df$s1x2 <- links_df$order

      #links_df$taxonomy_published <- links_df$s1
      #links_df$tax_simp <- links_df$s2

      u <- sapply(links_df$s2, function(x) {
         match(x, track_list[[1]]$segments$chord$sector)
         })
      uu <- as.numeric(xlim_list[[1]]$chord[,1])[u]

      links_df$s2x1 <- links_df$s2x1+uu
      links_df$s2x2 <- links_df$s2x2+uu
      links_df$s1 <- "main"
      links_df$s2 <- "chord"
   }



   ## PLOT TO FILE
   ## ::::::::::::::::
   circlize::circos.clear()
   my.file <- paste0(file.path(plot.dir, plot.name),".pdf")
   pdf(my.file, width = 12, height = 12, useDingbats = F)

   ## START PAR
   # circos.clear()
   circos.par(cell.padding = c(0, 0, 0, 0), gap.degree = 15, start.degree = (-45),  "track.height" = 0.2)
   circos.initialize(factors = group_names, xlim = matrix(unlist(lapply(xlim_list[[1]], range)), nrow=length(group_names), byrow=T))

      # Title
      title(plot.name)

      ## add DENDOGRAM
      max_height = attr(dend_list[["main"]], "height")

      circos.track(ylim = c(0, max_height), bg.border = NA, track.height = 0.3,
       panel.fun = function(x, y) {
        sector.index = get.cell.meta.data("sector.index")
           dend = dend_list[[sector.index]]
           if(sector.index=="main") circos.dendrogram(dend, max_height = max_height, facing = "inside")
      })



      ## add HCA cluster rects (cc)
      circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.1,
         panel.fun = function(x, y){
            xlim = CELL_META$xlim
            ylim = CELL_META$ylim
            m = xlim_list[[1]][[CELL_META$sector.index]]
            # here plot only the cc's as these are the main plot

            if(CELL_META$sector.index=="main") circos.rect(xleft = m[,1], ybottom = 0, xright = m[,2], ytop = 1,
                     col = color_key[as.character((track_list[[1]][["segments"]][[CELL_META$sector.index]]$sector))],
                     border = "black")

            if(CELL_META$sector.index=="main") circos.text(x = m[,3], y = 0.5, labels = as.character(track_list[[1]][["segments"]][[CELL_META$sector.index]]$sector), facing = "bending.inside", cex=0.5, col="white")

            if(CELL_META$sector.index=="chord") circos.rect(xleft = m[,1], ybottom = 0, xright = m[,2], ytop = 1,
                     col = NA, border = NA)
            })



      ## add second track and chord if available
      if(n_tot_annotations>1){
         for(i in 1:length(annotation.columns)){
            ii <- i+1

            circos.track(ylim = c(0,1), bg.border = NA, track.height=annotation.track.height[i], panel.fun = function(x, y) {
               if(CELL_META$sector.index=="main"){
                  m = xlim_list[[ii]][["main"]]
                  circos.rect(xleft = c(0:(max(xlim_list[[ii]]$main)-1)), ybottom = 0, xright = 1:max(xlim_list[[2]]$main), ytop = 1,
                        col = as.character(track_list[[ii]][["colors"]]$main), border = NA)
               }

               if(CELL_META$sector.index=="chord" & i==1){
                  m = xlim_list[[1]][["chord"]]
                  circos.rect(xleft = m[,1], ybottom = 0, xright = m[,2], ytop = 1,
                        col = unique(track_list[[1]]$colors[["chord"]]), border = "black")
                        circos.text(x = m[,3], y = 0.5, labels = track_list[[1]]$segments[["chord"]]$sector, facing = "clockwise", cex=0.6, col="black")
                        circos.text(x = m[,3], y = 0, labels = as.character(track_list[[1]]$segments[["chord"]]$n/scale_factors[2]), facing = "inside", cex=0.5, col="black", adj=c(0,1.5))
               }
               # if(CELL_META$sector.index=="chord" & i!=1){
               #    m = xlim_list[[1]][["chord"]]
               #    circos.rect(xleft = m[,1], ybottom = 0, xright = m[,2], ytop = 1, col = NA, border = NA)
               #    }
               })
            } # end i all tracks/annotation columns
         } # end n_tot_annots

         ## Links
      for(i in 1:nrow(links_df)){
         #links_df
         circos.link(
            sector.index1 = links_df$s1[i],
            point1 = c(links_df$s1x1[i], links_df$s1x2[i]),
            sector.index2 = links_df$s2[i],
            point2 = c(links_df$s2x1[i], links_df$s2x2[i]),
            col= color_subtypes[links_df$tax_simp[i]])
      }

      circos.clear()

     dev.off()
} # end function chord plot consensus clusters


         # source("~/PROJECTS/RCC_ccpRCC_2019/R_SCRIPTS/R_ccpRCC_2019_SOURCE.R", echo=TRUE)
         # setwd("~/PROJECTS/RCC_ccpRCC_2019/DNA_mutations/summary_plots/taxonomy_groups/")
         # require(maftools)
         # require(circlize)
         # #pdata <- pdata_panCanFull
         #
         # table(duplicated(names(gdc_cn)))
         # rcc_maf <- readRDS("~/RESOURCES/TCGA_PANCAN/Processed_data/MAF_RCC_tumors_708set_mc3_v028.rds")
         # ignore_genes <- c("TTN","MUC16","CSMD3")
         #
         #
         #
         # pdata_maf <- getClinicalData(rcc_maf)
         # pdata_maf <- pdata_maf %>% dplyr::select("sample_id","Tumor_Sample_Barcode","Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del" ,"In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Splice_Site", "Translation_Start_Site", "total")
         # pdata_maf <- pdata_panCanGex %>% dplyr::filter(sample_id %in% as.character(pdata_maf$sample_id)) %>% dplyr::left_join(pdata_maf, by="sample_id") %>% dplyr::filter(grepl("ccRCC|pRCC|chRCC|CIMP|ccpRCC|mesRCC|chONC|pONC", tax_simp)) %>% droplevels()
         # pdata_cn <- pdata_panCanFull %>% dplyr::filter(platform == "Genome_Wide_SNP_6") %>% dplyr::filter(!duplicated(sample_id)) %>% dplyr::filter(sample_id %in% names(gdc_cn)) %>% dplyr::filter(grepl("ccRCC|pRCC|chRCC|CIMP|ccpRCC|mesRCC|chONC|pONC", tax_simp)) %>% droplevels()
         #
         # pdata_maf
         # pdata_cn # A tibble: 847 x 133
         # my_groups <- levels(pdata_cn$tax_simp2)
         #
         #
         # ## Gene lists for all top genes
         # ## -------------------
         # #  getSampleSummary(rcc_maf)
         # # gene_list[["rcc"]] <- getGeneSummary(subsetMaf(rcc_maf, tsb=pdata_maf$Tumor_Sample_Barcode))[1:20, ]
         # gene_list <- lapply(my_groups, function(x, n.genes=15){
         #    my_barcodes <- as.character(pdata_maf %>% dplyr::filter(tax_simp2==x) %>% pull(Tumor_Sample_Barcode))
         #    y <- getGeneSummary(subsetMaf(rcc_maf, tsb=my_barcodes, ))
         #    y <- y[!y$Hugo_Symbol%in%ignore_genes, ]
         #    # yy <- dlfoo2::featureGetBM(y$Hugo_Symbol[1:n.genes], query.type = "gene_symbol")
         #    u <- match(y$Hugo_Symbol[1:n.genes], fData(rcc_t)$SYMBOL)
         #    u <- u[!is.na(u)]
         #    yy <- fData(rcc_t)[u,] %>% dplyr::select(chr, start, end, SYMBOL) %>% rename(Hugo_Symbol=SYMBOL)
         #    yy <- dplyr::left_join(yy, y, by="Hugo_Symbol")
         #    yy <- yy %>% mutate(percent=total/length(my_barcodes)*100) %>% mutate(percent=round(percent,1))
         #    yy <- yy %>% mutate(label=paste0(Hugo_Symbol," (",percent," %)"))
         #    return(yy)
         # })
         # names(gene_list) <- my_groups
         #
         #
         # ## Loop CN data thorugh all groups
         # for(i in 1:length(my_groups)){
         #    pdf(file.path("~/PROJECTS/RCC_ccpRCC_2019/CopyNumber/GenomomeHeatmaps/Circlize_Tax_Simp2", paste0("CopyNumber_Circle_",my_groups[i],".pdf")))
         #       my_cn <- gdc_cn[as.character(pdata_cn %>% dplyr::filter(tax_simp2 == my_groups[i]) %>% pull(sample_id))]
         #       dlfoo2::genomeCircleSegments(x=my_cn, plot.name = my_groups[i], genome.version = "hg38", value.column = "score", label.bed = gene_list[my_groups[i]])
         #    dev.off()
         # }
         #
         # i <- 1
         # my_cn <- gdc_cn[as.character(pdata_cn %>% dplyr::filter(tax_simp2 == my_groups[i]) %>% pull(sample_id))]
         # x <- my_cn
         # plot.name = my_groups[i]
         # genome.version = "hg38"
         # value.column = "score"
         # label.bed = gene_list[[my_groups[i]]]


#' Wrapper for Plotting (Multiple) genome tracks in Circlize -- One sample multiple tracks.
#' @description wrapper for circlize function described at https://jokergoo.github.io/circlize_book/book/modes-of-input.html
#' @description version 2.0. accepts multiple track types
#' @description  Plot genome data from BED/genomic ranges objects.
#' @description Input can be GenomicRanges object (if one 1 sample) or list of genomic ranges objects OR a BED style data.frame. Columns must include "ID" if multiple samples, "chr" if bed.
#' @family circlize
#' @family PlotWrap
#' @family copy number
#' @family genomic ranges
#' @param x GenomicRanges or bed style data.frame  (if one 1 sample) or list of genomic ranges objects OR a BED style data.frames (if multiple samples or if multiple tracks). Chromosomes should be in 'chrX' format (not 'X')
#' @param value.column name of the column in bed/GRanges that will define segment colors ('score' is default)
#' @param plot.mode If x is a list and if multiple tracks (not multiple samples for a single track) is to be plotted. If multiple tracks set to "multiple.track.types". Default is "single.track.type".
#' @param track.types What type of plot connected to each of the list objects. "rect" or "line.bar" supported. Character vector of length(x) if multiple and different tracks.
#' @param track.height Numeric vector. If different track.heights then this should be the length of x.
#' @param my.chr vector with what chromosomes to plot
#' @param genome.version Defaults to hg38. accepts multiple species
#' @param sort.pos what chromosome position to sort on (chr, start, end)
#' @param palette list of cahacter vectors. the palette(s) to use for the different tracks. If tracks are not numeic, then names of palette will be used to define colors.
#' @param palette.sat palette sat
#' @param plot.name name of plot
#' @param my.lwd line width for segments
#' @param label.bed bed style data.frame. if to label genes/regions in plot. Must include 'chr' (in chrX format), 'start' 'end' and 'label' column. May include 'color' column if user-defined colors.
#' @param label.side
#' @param label.cex
#' @param track.height
#' @param ideogram.plotType for circos.initializeWithIdeogram function. One, tow or all of "ideogram", "axis", "labels"
#' @param ideogram.track.height Height of the track which contains "axis" and "labels"
#' @param ideogram.ideogram.height Height of the ideogram track
#' @return genomeCircleSegments
#' @export
genomeCircleSegments <- function(
      x=my_cn,
      plot.mode = c("single.track.type"),
      plot.name = "",
      track.types = c("rect"),
      value.column = c("score"),
      genome.version = "hg38",
      track.height = NULL,
      palette = list(rev(c(rev(RColorBrewer::brewer.pal(9,"YlOrRd")[c(3:9,9)]), rep("white",1), RColorBrewer::brewer.pal(9,"Blues")[2:9]))),
      palette.sat = list(c(-0.75, 0.75)),
      my.chr = paste0("chr", c(1:22,"X")),
      sort.pos = c("chr3","10141635"),
      my.lwd = NULL,
      label.bed = NULL ,
      label.side="inside",
      label.cex=0.8,
      ideogram.plotType=c("ideogram", "axis", "labels"),
      ideogram.ideogram.height = convert_height(2, "mm")
){
   require(reshape2)
   require(dplyr)
   options(dplyr.width = Inf)
   require(tibble)
   require(tidyverse)
   require(tidyr)
   require(circlize)
   require(rtracklayer)
   require(genoset)
   require(GenomicRanges)
   require(maftools)

   ## Input data Single BED or GRanges (x is not a list)
   ## ::::::::::::::::::::::::::::::::
   if(!is.list(x)){
      x <- list(x)
      message("x is not a list - forcing plot.type to 'single.segment.type'")
      plot.mode == "single.segment.type"
   } #


   ## Input data List of BED or GRanges

   ## ::::::::::::::::::::::::
   if(class(x) %in% c("list","GRangesList")){
      ## If multiple segments (x is a list) - two possible plot modes - 'single.track.type' or 'multiple.track.types'
      # plot.mode will decide if single.track.type ( one tracl - multiple rows - with same settings) or multiple.track.types (multiple tracks with or wothout multiple samples)
      ## Check what type of plot - if multiple GRanges/tables (x is a list) then either plot all with same settings, or continue with lists of bed value.column, palette and palette.sat
      plot.mode = match.arg(plot.mode, choices = c("single.track.type", "multiple.track.types"),  several.ok = F)

      if(plot.mode == "single.track.type"){
         message("... 'plot.mode': ", plot.mode)
         message(" ... .... assuming all segments in x to be plotted with same value.column, palette and palette.sat ")
         message(" ... .... set plot.mode to 'multiple.segment.types' if to plot multiple segments with different value.column, palette and palette.sat ")

         u <- unique(lapply(x, function(y) class(y) %in% c("GRanges","data.frame")))
         if(u==T){
            message( "... .... ... x is a list of GRanges. generating bed_list from these. Assuming seqnames are chromsomes in chrX format.")
            bed <- dplyr::bind_rows(lapply(x, as.data.frame)) %>%
               dplyr::rename(chr=seqnames)
            bed_list <- list(bed)
            }
         u <- unique(lapply(x, function(y) class(y)=="data.frame"))
         if(u==T){
            message(" ... ... ... x is a list of Data Frames. generating bed_list from these. Assuming seqnames are chromsomes in chrX format.")
            x <- do.call("rbind", x)
            bed_list <- list(x)
         }
         # if(!is.list(value.column)) list(value.column)
         if(!is.list(palette)) palette <- list(palette)
         if(!is.list(palette.sat)) stop("palettee.sat must be provided as list")

      } # end if multiple samples

      ## if multiple.segment.types, i.e. single sample wit different types of segments to plot
      ## here it is possible to supply a mixed list of data.frames (bed-style) and GRanges
      if(plot.mode == "multiple.track.types"){
         message("... 'plot.mode': ", plot.mode)
         message(" ... .... checking if segments in x to be plotted with different value.column, palette and palette.sat. (if so, these have to be provided as lists of the same length as x)")
         message(" ... .... ... ")
         # if(!is.list(value.column)) value.column = list(value.column)
         if(!is.list(palette)) palette = list(palette)
         if(!is.list(palette.sat)) stop("palettee.sat must be provided as list")

         # if(!is.list(palette.sat)) palette.sat = list(palette.sat)
         if(!length(x)==length(value.column)){
            message("... ... ... Warning: value.column is not a list with length(x). Setting all value columns to: ", value.column[1] )
            value.column.temp <- vector(mode="character", length = length(x))
            value.column.temp[1:length(x)] <- value.column[1]
            value.column <- value.column.temp
            }
         if(!length(x)==length(palette)){
            message("... ... ... Warning: palette is not a list with length(x). Setting all value columns to: ", palette[[1]] )
            palette.temp <- vector(mode="list", length = length(x))
            palette.temp[1:length(x)] <- palette[[1]]
            palette <- palette.temp
            }
         if(!length(x)==length(palette.sat)){
            message("... ... ... Warning: palette.sat is not a list with length(x). Setting all value columns to: ", paste(palette.sat[[1]], collapse = ", ") )
            palette.sat.temp <- vector(mode="list", length = length(x))
            palette.sat <- lapply(palette.sat.temp, function(x) palette.sat[[1]])
         }

         bed_list <- list()
         for(i in 1:length(x)){
            if(class(x[[i]]) == "GRanges"){
               message( "... .... ... converting GRanges to bed-style data.frame. Assuming seqnames are chromsomes in chrX format.")
               bed_list[[i]] <- as.data.frame(x[[i]]) %>%
                  dplyr::rename(chr=seqnames)
               }
            if(class(x[[i]]) == "data.frame"){
               message(" ... ... ... found Data Frame. Assuming bed-style and that seqnames are chromsomes in chrX format.")
               bed_list[[i]] <- x[[i]]
            }
         }
      } # end if mode multiple.track.types
   } # end if x is a list



   ## Check data integrity (geneomic position column names)
   ## -----------------------------------
   lapply(bed_list, function(x){
      if(!all(sapply(c("chr","start","end","ID"), function(y)
         y %in% colnames(x)))) {
         stop("colnames must include: 'chr', 'start', 'end', and 'ID'")
         }
   })
   bed_list <- lapply(bed_list, function(x){
      x$chr <- as.character(x$chr)
      return(x)
   })



   ## Check data integrity (data value.column)
   ## -----------------------------------
   if(!is.character(value.column)) stop("value column should be a cahracer vector matching value columns to color in bed/GRanges objects")
   if(length(bed_list)>1){
      if(length(value.column)==length(x)) {
         message("... length(x) 'value.column's provided. using as is")
      }
      if(length(value.column)!=length(x)) {
         value.column <- rep(value.column[1], length(bed_list))
         message("The number of supplied 'value.column's do not match length of x. All 'value columns' set to:  " , value.column)
      }
   }

   for(i in 1:length(bed_list)){
      if(!value.column[i] %in% colnames(bed_list[[i]])) stop(paste0(value.column[i],": colnames must include column as defined by 'value.column'"))
       bed_list[[i]] <- bed_list[[i]] %>% dplyr::mutate(score = UQ(rlang::sym(value.column[i])))
   }

   ## Check track.types
   ## -----------------
   if(length(bed_list)>1){
      if(length(track.types)==length(x)) {
         message("... multiple 'track.types' provided. using as is")
      }
      if(length(track.types)!=length(x)) {
         track.types <- rep(track.types[1], length(bed_list))
         message("The number of supplied track.types do not match length of x. All 'track types' set to:  " , track.types)
      }
   }
   if(!all(track.types %in% c("rect","line.bar"))) stop("'track.types' must be a single character vector of length 1 or length(x). Allowed values are 'rect' or 'line.bar'")

   # generate segment colors
   # :::::::::::::::::
   message("... generating segment colorings from palettes")
   for( i in 1:length(bed_list)){
      ## If a given palette has names (named character vector) - then assume score values should be matched with these names
      ## i.e. in this case the score will be treated as a character vector (not numeric). The palette.sat will thus be ignored.
      if(is.character(bed_list[[i]]$score)) if(is.null(names(palette[[i]]))) stop(" ... ... No Names found for given palette. Must be named character vector if value column is character: ", value.column[i])
      if(!is.null(names(palette[[i]]))){
         message(" ... ... Names found for given palette. Setting score values to character and matched with these palette names: ", value.column[i])
         bed_list[[i]]$color <- palette[[i]][as.character(bed_list[[i]]$score)]
      }else{
         message(" ... ... generating ramped colors using palette and palette.sat ")
         bed_list[[i]]$color <- heatcol_creator2(bed_list[[i]]$score, palette.sat = palette.sat[[i]], my.pal = palette[[i]])$col
      }
   }



   # Sort ID's  & define Y positons from ID (defined by ID in sorting_scores data frame )
   # :::::::::::::::::::::::::::
      ## Only applicabe if multiple ID's present in individual beds (within the bed_list - in particular if single.segment.type and multiple samples)

   n.samples_list <- list()
   for(i in 1:length(bed_list)){
      if(!(sort.pos[1] %in% bed_list[[i]]$chr)){
         message("sort.pos not in data")
         sort.pos <- c(bed_list[[i]]$chr[1], bed_list[[i]]$start[1])
      }
      sorting_df <- bed_list[[i]] %>% filter(chr==sort.pos[1])
      sorting_gr <- GenomicRanges::makeGRangesFromDataFrame(sorting_df)
      elementMetadata(sorting_gr) <- DataFrame(score=as.numeric(sorting_df$score), ID=sorting_df$ID)
      gr_list <- split(sorting_gr, sorting_gr$ID)
      gr_query = GRanges(seqnames=sort.pos[1], ranges=IRanges(start=as.numeric(sort.pos[2]), end=as.numeric(sort.pos[2])))

      sorting_scores <- lapply(gr_list, function(y) {
         u <- GenomicRanges::nearest(gr_query, y)
         return(data.frame(score=y[u]$score, ID=y[u]$ID))
      })
      sorting_scores <- do.call("rbind", sorting_scores)
      sorting_scores <- dplyr::left_join(data.frame(ID=unique(bed_list[[i]]$ID)), sorting_scores)
         # mutate(score = if_else(is.na(score), 0, score))   #fix if any IDs not present - set score to 0
      sorting_scores <- sorting_scores %>% arrange(-sorting_scores$score,sorting_scores$ID) # sort on score, then sample_id
      bed_list[[i]] <- bed_list[[i]] %>% mutate(y = factor(ID, levels=as.character(sorting_scores$ID))) %>% mutate(y = as.numeric(y))
      my_ylim <- c(0, nrow(sorting_scores)+1)

      n.samples_list[[i]] <- nrow(sorting_scores)
   }


   ## Use only speifed chromosomes
   ## :::::::::::::::::::::::::
   for(i in 1:length(bed_list)){
      bed_list[[i]] <- bed_list[[i]] %>% dplyr::filter(chr %in% my.chr) %>% arrange(y, abs(as.numeric(score)))
   }

   # generate segment y-coordinates (value column) (if point ... or line?)
   # :::::::::::::::::
   for(i in 1:length(bed_list)){
      if(track.types[i] == "line.bar"){
         bed_list[[i]]$value1 <- 0.45
      }
   }


   # region_list: create region for each chromosome
   # :::::::::
   region_list <- list()
   for(i in 1:length(bed_list)){
      region_list[[i]] <- split(bed_list[[i]], bed_list[[i]]$chr)
      #length(region_list)
   }

   # my_lwd
   ## -------
   if(is.null(my.lwd)){
      # my.lwd = 100/n.samples
      my.lwd = 100/sum(unlist(n.samples_list))
   }

   # Label regions/genes
   if(!is.null(label.bed)){
      if(!c("color" %in% colnames(label.bed))) {
         label.bed <- label.bed %>% mutate(color="black")
      }
      label.bed <- label.bed %>% dplyr::filter(chr %in% my.chr) %>% mutate(chr = as.character(chr))
      label.bed_list <- split(label.bed, label.bed$chr)
   }

   # define track height if null or if multiple tracks
   # ----------
   if(!is.null(track.height)) {
      my.track.height=track.height
   }
   if(is.null(track.height) & plot.mode=="single.track.type") my.track.height <- 0.3
   if(is.null(track.height) & plot.mode=="multiple.track.types") my.track.height <- 0.1

   if(plot.mode=="multiple.track.types"){
      if(length(my.track.height) != length(bed_list)){
         message("... only one track hight provieded for ", length(bed_list), " tracks. Setting to default track height (0.1 or 0.3)")
         my.track.height <- rep(0.1, length(bed_list)) # /length(bed_list)
      }else{
         message(" ... individual track heights provided. using as is.")
      }
   }


   ## PLOT ::::::::
   ## ::::::::::::::::::
   message("... PLOTTNG")
   circos.initializeWithIdeogram(species = genome.version, chromosome.index = my.chr, plotType = ideogram.plotType, ideogram.height = ideogram.ideogram.height)
   text(0, 0, plot.name, cex = 1)
   for(ii in 1:length(bed_list)){
      track.types.i <- track.types[ii]

      if(track.types.i == "rect"){
         my_region_list <- region_list[[ii]]
         track.hight.i <- my.track.height[ii]
         circos.genomicTrack(my_region_list, ylim = my_ylim, track.height = track.hight.i, panel.fun = function(region, value, ...) {
               i = getI(...)
               circos.genomicRect(
                  region, value,
                  ybottom=my_region_list[[CELL_META$sector.index]]$y-0.85,
                  ytop=my_region_list[[CELL_META$sector.index]]$y+0.85,
                  lty=0, border=NA,
                  col=my_region_list[[CELL_META$sector.index]]$color, ... )
         })
      }

      if(track.types.i == "line.bar"){
         my_region_list <- lapply(region_list[[ii]], function(x){
           return(x[,c("chr","start","end","value1","color")])
         })
         track.hight.i <- my.track.height[ii]
         circos.genomicTrack(
            my_region_list, ylim = c(-0.5,0.5),
            track.height = track.hight.i,
            panel.fun = function(region, value, ...) {
               # i = getI(...)
               circos.genomicLines(
                  region, value, type="h",
                  col=my_region_list[[CELL_META$sector.index]]$color)
                 #  ...)
            })

         }

   } # end all i

   if(!is.null(label.bed)){
      message(" ... ... adding labels deined in label.bed")
      circos.genomicLabels(label.bed, labels.column = "label", side = label.side, col=label.bed$color, cex=label.cex)
   }

   circos.clear()

} # end function genomeCircleSegments



#' Wrapper for Plotting Copy number prfiles of Chrs in Circlize. each sample represented by a segment
#' @description wrapper for circlize function described at https://jokergoo.github.io/circlize_book/book/modes-of-input.html
#' @description  Plot genome data from BED/genomic ranges objects.
#' @description Input can be GenomicRanges object (if one 1 sample) or list of genomic ranges objects OR a BED style data.frame. Columns must include "ID" if multiple samples, "chr" if bed.
#' @family circlize
#' @family PlotWrap
#' @family copy number
#' @family genomic ranges
#' @param plot.name name of file to plot
#' @param plot.dir directory to plot in
# #' @param plot.quartz if to plot to quartz insted
#' @param x GenomicRanges object (if one 1 sample) or list of genomic ranges objects OR a BED style data.frame
#' @param value.column name of column in bed/granges to define copy numnber colors (score)
#' @param my.chr vector with what chromosomes to plot
#' @param genome.version Defaults to hg38. accepts multiple species
# #' @param sort.pos what chromosome position to sort on (chr, start, end)
#' @param palette palette
#' @param palette.sat palette sat
#' @param my.lwd line width for segments
#' @param label.bed bed style data.frame. if to label genes/regions in plot. Must include 'label' column. May include 'color' column if user-defined colors.
#' @param label.side
#' @param label.cex
#' @param track.height
#' @return plot circle
#' @export

genomeCircleSegments_legacy <- function(
   x,
   plot.name = "temp",
   #plot.dir = NULL,
   # plot.quartz = FALSE,
   my.chr = paste0("chr", c(1:22,"X")), genome.version = "hg38",
   sort.pos = c("chr3","10141635"),
   value.column = "score",
   palette = rev(c(rev(RColorBrewer::brewer.pal(9,"YlOrRd")[2:9]), rep("white",3), RColorBrewer::brewer.pal(9,"Blues")[2:9])),
   palette.sat = c(-0.75, 0.75),
   my.lwd = NULL,
   label.bed = NULL, label.side="inside", label.cex=0.8,
   track.height = 0.3
){


   require(reshape2)
   require(dplyr)
   options(dplyr.width = Inf)
   require(tibble)
   require(tidyverse)
   require(tidyr)
   require(circlize)
   require(rtracklayer)
   require(genoset)
   require(GenomicRanges)
   require(maftools)


   ## Input data BED or GRanges
   ## ::::::::::::::::::::::::
   if(class(x) %in% c("list","GRangesList")){
      u <- unique(lapply(x, function(y) class(y)=="GRanges"))
      if(u==T){
         message("x is a list of GRanges. generating bed from these. Assuming seqnames are chromsomes in chrX format.")
         bed <- do.call("rbind", lapply(x, as.data.frame)) %>% dplyr::rename(chr=seqnames)
         }
      u <- unique(lapply(x, function(y) class(y)=="data.frame"))
      if(u==T){
         x <- do.call("rbind", x)
      }
   }
   if(class(x)=="GRanges"){
      "x is GRanges object. generating bed from these. Assuming seqnames are chromsomes in chrX format."
      bed <- as.data.frame(x) %>% dplyr::rename(chr=seqnames)
   }
   if(class(x)=="data.frame"){
      bed <- x
      if(!c("ID" %in% colnames(bed))){
         "no ID column in x. Assuming only one sample exists"
         bed$ID <- 1
      }
   }
   if(class(bed)!="data.frame") stop("x cannot be converted to bed")


   # check column names
   if(!all(sapply(c("chr","start","end","ID",value.column), function(x)  x %in% colnames(bed)))) {message("colnames must include: chr, start, stop, ID and value column defined by 'value.column'")}
   if(value.column != "score") bed <- bed %>% dplyr::rename(score=value.column)
   bed$chr <- as.character(bed$chr)

   if(!is.null(label.bed))
   if(!all(sapply(c("chr","start","end","label"), function(x)  x %in% colnames(label.bed)))) {message("colnames for label.bed must include: chr, start, stop, and label")}

   message("... ",nrow(bed)," segments")
   message("... ",length(unique(bed$ID))," IDs")  #


   # generate segment colors
   # :::::::::::::::::
   bed$color <- heatcol_creator2(bed$score, palette.sat = palette.sat, my.pal = palette )$col

   # Sort & define Y positons from ID (defined by ID in sorting_scores data frame )
   # :::::::::::::::::::::::::::
   if(!(sort.pos[1] %in% bed$chr)){
      message("sort.pos not in data")
      sort.pos <- c(bed$chr[1], bed$start[1])
   }
   sorting_df <- bed %>% filter(chr==sort.pos[1])
   sorting_gr <- GenomicRanges::makeGRangesFromDataFrame(sorting_df)
   elementMetadata(sorting_gr) <- DataFrame(score=sorting_df$score, ID=sorting_df$ID)
   gr_list <- split(sorting_gr, sorting_gr$ID)
   gr_query = GRanges(seqnames=sort.pos[1], ranges=IRanges(start=as.numeric(sort.pos[2]), end=as.numeric(sort.pos[2])))

   sorting_scores <- lapply(gr_list, function(y) {
      u <- GenomicRanges::nearest(gr_query, y)
      return(data.frame(score=y[u]$score, ID=y[u]$ID))
   })
   sorting_scores <- do.call("rbind", sorting_scores)
   sorting_scores <- dplyr::left_join(data.frame(ID=unique(bed$ID)), sorting_scores)
      # mutate(score = if_else(is.na(score), 0, score))   #fix if any IDs not present - set score to 0
   sorting_scores <- sorting_scores %>% arrange(-sorting_scores$score,sorting_scores$ID) # sort on score, then sample_id
   bed <- bed %>% mutate(y = factor(ID, levels=as.character(sorting_scores$ID))) %>% mutate(y = as.numeric(y))
   my_ylim <- c(0, nrow(sorting_scores)+1)
   n.samples <-  nrow(sorting_scores)

   ## Use only speifed chromosomes
   bed <- bed %>% dplyr::filter(chr %in% my.chr) %>% arrange(abs(score))


   # region_list: create region for each chromosome
   # :::::::::
   region_list <- split(bed, bed$chr)
      #length(region_list)

   # my_lwd
   ## -------
   if(is.null(my.lwd)){
      my.lwd = 100/n.samples
   }

   # Label regions/genes
   if(!is.null(label.bed)){
      if(!c("color" %in% colnames(label.bed))) {
         label.bed <- label.bed %>% mutate(color="black")
      }
      label.bed <- label.bed %>% dplyr::filter(chr %in% my.chr) %>% mutate(chr = as.character(chr))
      label.bed_list <- split(label.bed, label.bed$chr)
   }


   ## PLOT ::::::::
   circos.initializeWithIdeogram(species = genome.version, chromosome.index = my.chr)
   text(0, 0, plot.name, cex = 1)

      # circos.genomicTrack(region_list, ylim = my_ylim, track.height = track.height, panel.fun = function(region, value, ...) {
      #       i = getI(...)
      #       circos.genomicLines(region, value=region_list[[CELL_META$sector.index]]$y, type = "segment", lwd = my.lwd, border=NA, col=region_list[[CELL_META$sector.index]]$color, ... )
      # })

      circos.genomicTrack(region_list, ylim = my_ylim, track.height = track.height, panel.fun = function(region, value, ...) {
            i = getI(...)
            circos.genomicRect(region, value, ybottom=region_list[[CELL_META$sector.index]]$y-0.5,
               ytop=region_list[[CELL_META$sector.index]]$y+0.5, lty=0,
               border=NA, col=region_list[[CELL_META$sector.index]]$color, ... )
      })


      if(!is.null(label.bed)){
         circos.genomicLabels(label.bed, labels.column = "label", side = label.side, col=label.bed$color, cex=label.cex)
      }

   circos.clear()

} # end function genomeCircleSegments





  #  my_res <- readRDS("~/OneDrive - Lund University/PROJECTS_od/RCC_VHPT/DESeq2/EnrichmentAnalysis/fGSEA/fgsea_Apeglm_lfc_all_combos_smythies_nperm1e6/fgsea_Apeglm_lfc_all_combos_smythies_nperm1e6_fgsea_allOntologies.rds")
  #  my_pathways <- readRDS(file="~/RESOURCES/MSigDB/msigdb_v7.0_GMTs/smythies_hif_canonical.symbols_pathwayList.rds")
  #   my_files <- list.files("~/OneDrive - Lund University/PROJECTS_od/RCC_VHPT/DESeq2/GSEA/gsea_rnk/", full.names = T, pattern = "Apeglm")
  #  rank_list <- lapply(my_files, function(x){
  #   y <- read.delim(x, header = T, sep="\t")
  #   yy <- y$rank.out
  #   names(yy) <- as.character(y$id.out)
  #     return(yy)
  #  })
  # names(rank_list) <- gsub(".*Apeglm_","",my_files)
  # names(rank_list) <- gsub("_lfc.*","",names(rank_list))
  #
  #  str(rank_list)

   # pathway_list <- my_pathways
   # rnk_list <- rank_list
   # fgsea_res <- my_res [[1]]

#
#
# #' adaption of the plotGseaTable & plotEnrichment {fgsea} from fgsea package
# #' @description faceted waterfall plot of multiple GSEAs over multiple samples.
# #' @description  adaption of the plotGseaTable from fgsea package https://bioconductor.org/packages/release/bioc/manuals/fgsea/man/fgsea.pdf
# #' @family enrichment
# #' @family plotwrap
# #' @family msig
# #' @family gsea
# #' @family fgsea
# #' @param run.pathwats gene set(s) to plot. character vector. warns if >25.
# #' @param pathway_list list of pathways to plot. All these pathwasys will be plotted and for all sample rank lists. Should be in the form of 'pathways' for fgsea analysis, i.e. a list with character vectors of genes within wah pathway. If null, the entire msig db will be loaded and list flattened.
# #' @param rnk_list list of named integer ranks for a series of samples or experiments. Each vector should have names corresponding to ids used in the pathways.
# #' @param fgsea_res fgsea results object obtained from using the rnk_list and pathways from which pathway_list is subset from
# #' @return ggplot object
# #' @export
# #' @examples
# #' \dontrun{
# #' ## load pathwys from MSIg and create subseted pathway list before run
# #' fgsea_res_tab <- read.delim("~/PROJECTS/RCC_cc_eos_HN/Limma/Enrichment/fgsea_nperm1e6/fgsea_nperm1e6_fgsea_Ontologies_significant_Plot.txt", as.is=T)
# #'   pathway_tab <- fgsea_res_tab[fgsea_res_plot$plot!="", ]
# #'   run.pathways <- pathway_tab$pathway
# #'
# #'   pathway_files <- list.files("~/RESOURCES/MSigDB/msigdb_v7.0_GMTs/", pattern="symbols_pathwayList.rds", full.names = T)
# #'   pathway_list <- lapply(pathway_files, readRDS)
# #'   str(pathway_list)
# #'   pathway_list <- unlist(pathway_list, recursive = F)
# #'   names(pathway_list)
# #'   ## select from pathway list
# #'   pathway_list <- pathway_list[my_pathways]
# #'
# #'   # load fgsea res
# #'   fgsea_res <- readRDS("~/PROJECTS/RCC_cc_eos_HN/Limma/Enrichment/fgsea_nperm1e6/fgsea_nperm1e6_fgsea_allOntologies.rds")
# #'   str(fgsea_res)
# #'
# #'   # load Rnk tav
# #'   rank_tab <- read.delim(file="~/PROJECTS/RCC_cc_eos_HN/Limma/eos_vs_cc_Paired4v4_fullData16238.txt",as.is=T)
# #'   y <- rank_tab$logFC
# #'   names(y) <- rank_tab$SYMBOL
# #'   table(duplicated(y))
# #'   #    FALSE  TRUE
# #'   # 15012  1226
# #'   y <- y[!duplicated(y)]
# #'   table(duplicated(names(y)))
# #'   y <- y[!duplicated(names(y))]
# #'   rnk_list <- list()
# #'   rnk_list[["cc_eos"]] <- y
# #' }
#
# gseaPlot <- function(
#    plot.pathways=NULL,
#    pathway_list=NULL,
#    rnk_list,
#    fgsea_res,
#    score.pal = c(rev(RColorBrewer::brewer.pal(9,"YlOrRd")[1:9]),rep("#FFFFCC",2),rep("white",5), rep("#FFFFCC",2), RColorBrewer::brewer.pal(9,"PuBu")[1:9])
# ){
#    require(dplyr)
#    require(tibble)
#    require(tidyr)
#    require(ggplot2)
#    require(fgsea)
#
#
#
#    fgsea_run_names <- unique(unlist(lapply(fgsea_res, function(x) names(x))))
#    if(!all(names(rnk_list) %in% fgsea_run_names)) stop("not all results in fgsea_results object present in rnk_list. Names must match ... ")
#    if(!length(plot.pathways)) stop("you must supply valid plot.pathways as character vector and present in pathway_list")
#
#    if(is.null(pathway_list)){
#       sig.dir <- "~/RESOURCES/MSigDB/msigdb_v7.0_GMTs/"
#       message("... loading pathways from local repository: ", sig.dir)
#       #pathway_files <- list.files(sig.dir, pattern="symbols_pathwayList.rds", full.names = T)
#       pathway_list <- lapply("/Users/david/RESOURCES/MSigDB/msigdb_v7.0_GMTs//msigdb.v7.0.symbols_pathwayList.rds", readRDS)
#       #str(pathway_list)
#       pathway_list <- unlist(pathway_list, recursive = F)
#       names(pathway_list)
#       ## select from pathway list
#       }
#    if(!all(plot.pathways %in% names(pathway_list))) stop("not all supplied plot.pathways present in pathway_list")
#    pathway_list <- pathway_list[plot.pathways]
#
#    run_names <- names(rnk_list) ## will define what runs/experiments that will be run and in what order
#
#
#
#
#    ## flatten fgsea_res into data frame instad of list
#    ## ----
#    fgsea_tab <- lapply(names(fgsea_res), function(x){
#
#       lapply(names(fgsea_res[[x]]), function(y, xx=x){
#          z <- fgsea_res[[xx]][[y]]
#          z$ontology <- xx
#          z$condition <- y
#          return(z)
#       })
#    })
#    fgsea_tab <- unlist(fgsea_tab, recursive = F)
#    fgsea_tab <- dplyr::bind_rows(fgsea_tab)
#    fgsea_tab <- fgsea_tab[match(plot.pathways, fgsea_tab$pathway),]
#
#
#    ## check all rnk_lists & Convert to tables with rnk_coordinates & pathway names
#    ## ----
#    rnk_lengths <- lapply(rnk_list, function(x) length(x))
#    if(!all(rnk_lengths==rnk_lengths[[1]])) stop("not all rnk lists are the same length", rnk_lengths)
#    my.xlim <- c(0,rnk_lengths[[1]])
#
#    # ## pathway_list names defines what pathways, then all genes in these pathwyas will be queried in all rnk_lists
#    rnk_tab <- lapply(plot.pathways, function(x){
#       my_list <- lapply(names(rnk_list), function(y, xx=x){
#          my_rnk <- rnk_list[[y]][order(rnk_list[[y]], decreasing=T)]
#          #if(length(which(is.na(my_rnk)))>1) message("warning: NA values among rnk values")
#          #if(length(which(duplicated(my_rnk)))>1) message("warning: multiple duplicates in rnk list. may produce unreproducible results")
#          u <- match(pathway_list[[xx]], names(my_rnk))
#          u
#
#          d <- data.frame(
#             rnk.x = sort(u[!is.na(u)], decreasing = F),
#             rnk.y = sort(my_rnk[u[!is.na(u)]], decreasing=T),
#             pathway=xx,
#             condition=y
#             )
#          return(d)
#          })
#       names(my_list) <- names(rnk_list)
#       return(my_list)
#       # return(dplyr::bind_rows(rnk_tab))
#    })
#    # names(rnk_tab) <- plot.pathways
#
#    ## Create plot_tab
#    # pathway_stats <- lapply(names(pathway_list), function(x){
#    #    my_list <- lapply(names(rnk_list), function(y, xx=x){
#    #       dd <- fgsea_res[[y]][match(xx, fgsea_res[[y]]$pathway),c("pathway","pval","padj","NES")]
#    #       dd$condition <- y
#    #       return(dd)
#    #    })
#    #    names(my_list) <- names(rnk_list)
#    #    return(dplyr::bind_rows(my_list))
#    # })
#    # pathway_rnk_tab <-  bind_rows(pathway_rnk)
#    # pathway_stats_tab <- bind_rows(pathway_stats_tab)
#
#    ## combine rnk_tab and
#    rnk_tab <- dplyr::bind_rows(unlist(rnk_tab, recursive = F))
#
#    rnk_tab <- rnk_tab %>%
#
#       mutate(y.direction =  if_else(rnk.y>0, "up", "dn")) %>%
#       mutate(pathway_condition = factor(paste(pathway,condition, sep="_"))) %>% arrange(pathway, condition, rnk.x)
#
#    #
#    #
#    # ggplot(pathway_rnk_tab, aes(x=rnk.x, y=rnk.y, fill=y.direction, color=y.direction)) +
#    #    scale_fill_discrete(guide="none") +
#    #    scale_color_discrete(guide="none") +
#    #    geom_bar(stat="identity", width=0.7, position = position_dodge(width=0.4)) +
#    #    scale_fill_manual(aes(fill=y.direction))
#    #
#    # color_key <- c("red","blue")
#    # names(color_key) <- c("up","dn")
#    #
#    #  ggplot(pathway_rnk_tab, aes(x=rnk.x, y=rnk.y, fill=y.direction, color=y.direction)) +
#    #    geom_bar(stat="identity", width=0.1, position = position_dodge(width=0.05), alpha=0.5) +
#    #    scale_fill_manual(values = color_key) +
#    #    scale_color_manual(values = color_key)
#    #
#    #
#    #
#    #
#    # ggplot(pathway_rnk_tab %>% arrange(rnk.x), aes(x=rnk.x, y=rnk.y, color=rnk.x)) +
#    #    coord_cartesian(xlim=my.xlim) +
#    #    geom_bar(stat="identity", width=0.001, alpha=0.25) +
#    #    scale_color_gradientn(guide="none",colours = rev(score.pal), limits=c(0,my.xlim[2])) +
#    #    facet_wrap(pathway~condition, scales="free", nrow = length(names(my_pathways)))
#    #    #scale_color_gradientn(guide="none",colours = rev(score.pal), limits=c(0,my.xlim[2]))
#    #
#
# ticksSize <- 0.7
#
#       ggplot(rnk_tab, aes(x = rnk.x, y = rnk.y)) +
#          facet_wrap(~pathway_condition) +
#          coord_cartesian(xlim=my.xlim, ylim =  c(-ceiling(max(rnk_tab$rnk.y)),ceiling(max(rnk_tab$rnk.y)))) +
#          # geom_point(color = "green", size = 0.1) +
#          geom_point(aes(color = rnk.y), size = 0.5) +
#
#
#          geom_hline(yintercept = max(rnk_tab$rnk.y), colour = "red", linetype = "dashed") +
#          geom_hline(yintercept = min(rnk_tab$rnk.y), colour = "red", linetype = "dashed") +
#          geom_hline(yintercept = 0, colour = "black") + geom_line(color = "light blue") +
#          theme_bw() +
#          geom_segment(mapping = aes(x = rnk.x, y = 0, xend = rnk.x, yend = rnk.y, color=rnk.y), size = ticksSize) +
#          scale_color_gradientn(guide="none",colours = rev(score.pal)) +
#          theme(panel.border = element_blank(), panel.grid.minor = element_blank()) +
#          labs(x = "rank", y = "enrichment score")
#
#
#     #        geom_segment(data = data.frame(x = pathway), mapping = aes(x = x,
#     #        y = -diff/2, xend = x, yend = diff/2), size = ticksSize) +
#
#
#
#    return(g)
#   #ggsave(g, filename = "Correlation/TCGA/EA_plot_fgsea_tcga.pdf", width=18, height=10, useDingbats = F)
#   }
#
#






# ea_pathways <- read.delim("~/PROJECTS/BRCA_SR_KP/DESeq2_analysis/EnrichmentAnalysis/fGSEA/fgsea_caf_e2_all_combos/fgsea_caf_e2_all_combos_fgsea_Ontologies_significant_PLOT.txt", as.is=T)
#    table(ea_pathways$PLOT)
#    ea_pathways <- unique(ea_pathways$pathway[ea_pathways$PLOT!=""])
#    ea_tab <- read.delim("~/PROJECTS/BRCA_SR_KP/DESeq2_analysis/EnrichmentAnalysis/fGSEA/fgsea_caf_e2_all_combos/fgsea_caf_e2_all_combos_fgsea_Ontologies_all.txt", as.is=T)
#    ea_tab <- ea_tab[ea_tab$pathway %in% ea_pathways,]
#    str(ea_tab)

   # ea.tab <- as_tibble(ea_tab) %>% dplyr::rename(group=tissue, term=pathway, score=NES, size.n=size) %>% mutate(size=abs(score)) %>% mutate(border=if_else(padj<0.05, "black", "white"))




#' Generic plot for enrichment analyses
#' @description input is dataframe with, group, term, n, and p columns
#' @description  https://bioconductor.org/packages/release/bioc/manuals/fgsea/man/fgsea.pdf
#' @description plot order is determined by clustering - if not then set term column to factor and let levels decide order
#' @description v2.0 20191016. fix pbug with color scales. Also move to plotting NES instead of pvalue ... and NES for ALL not only the significant.
#' @family enrichment
#' @family plotwrap
#' @family msig
#' @family gsea
#' @family fgsea
#' @param ea.tab, data frame with enrichment analysis to plot.
#' @param score.title name of score value, e.g. p-value
#' @param size.title name of size value, e.g. p-value
#' @param border.col set to NA. If column named border is icluded then this column will be used
#' @param score.pal, color palette for p-val. if negative sizes (n).
#' @param order.by.clustering.x, let order of terms/groups be decided by clustering. Otherwise defined by levels
#' @param order.by.clustering.y, let order of terms/groups be decided by clustering. Otherwise defined by levels
#' @param reverse.x if to reverse x-axis
#' @param reverse.y if to reverse y-axis
#' @return ggplot object
#' @export
enrichmentPlot_bubble <- function(
      ea.tab=NULL,
      score.title ="score", size.title ="size",
      order.by.clustering.x=T, order.by.clustering.y=T,
      border.col = NA,
      #score.pal = c(rev(RColorBrewer::brewer.pal(9,"YlOrRd")[1:9]), RColorBrewer::brewer.pal(9,"PuBuGn")[1:9]),
      #score.pal = c(rev(RColorBrewer::brewer.pal(9,"Reds")[1:9]),"white", RColorBrewer::brewer.pal(9,"Blues")[1:9]),
      score.pal = c(rev(RColorBrewer::brewer.pal(9,"YlOrRd")[1:9]),rep("#FFFFCC",2),rep("white",5), rep("#FFFFCC",2), RColorBrewer::brewer.pal(9,"YlGn")[1:9]) ,
      reverse.x=F, reverse.y=F
   ){
   require(dplyr)
   require(tibble)
   require(tidyr)
   require(ggplot2)


   ea.tab <- data.frame(ea.tab)
   if(!all(sapply(c("group","term","size","score"), function(x)  x %in% colnames(ea.tab)))) {message("colnames must include: group, term,size, score")}

   if(!is.factor(ea.tab$term)) ea.tab$term <- factor(ea.tab$term)
   if(!is.factor(ea.tab$group)) ea.tab$group <- factor(ea.tab$group)
   if(!c("border" %in% colnames(ea.tab))) ea.tab$border <- border.col

   ea.tab <- droplevels(ea.tab)
   my_terms <- levels(ea.tab$term)
   my_groups <- levels(ea.tab$group)

   my_tab <- as_tibble(ea.tab) %>%
      arrange(group) %>%
      mutate(x0 = as.numeric(term)) %>%
      mutate(y0 = as.numeric(group))



   # if to order rows/cols by cluster analysis
   if(any(order.by.clustering.x, order.by.clustering.y)){
      my_mat <- matrix(nrow=max(my_tab$y0), ncol=max(my_tab$x0), data=0, dimnames=list(levels(my_tab$group), levels(my_tab$term)))
      for(i in 1:nrow(my_tab)){
         my_mat[my_tab$y0[i], my_tab$x0[i]] <- my_tab$score[i]
      }

      mydist=function(c) {amap::Dist(c,method="pearson")}
      myclust=function(c) {hclust(c,method="ward.D2")}
      hca.y <- my_mat %>% mydist %>% myclust %>% as.dendrogram
      hca.x <- t(my_mat) %>% mydist %>% myclust %>% as.dendrogram

      if(order.by.clustering.x){
        my_tab <- my_tab %>%  mutate(term = factor(term, levels=labels(hca.x)))
      }
      # if(!order.by.clustering.x){
      #   my_tab <- my_tab %>%  mutate(group = factor(group, levels=my_terms))
      # }
      if(order.by.clustering.y){
        my_tab <- my_tab %>%  mutate(group = factor(group, levels=labels(hca.y)))
      }

      # if(!order.by.clustering.y){
      #   my_tab <- my_tab %>%  mutate(group = factor(group, levels=my_groups)) ## ?!?!?"?!"??
      # }
      my_tab <- my_tab %>% arrange(group, term, -size)
   }

   if(reverse.x) my_tab <- my_tab %>% mutate(term = factor(term, levels=rev(levels(term)))) %>% arrange(group, term)
   if(reverse.y) my_tab <- my_tab %>% mutate(group = factor(group, levels=rev(levels(group)))) %>% arrange(group, term)

   # str(unique(my_tab$group))

   x_labs <- levels(my_tab$term)
   y_labs <- levels(my_tab$group)
   x.n <- length(x_labs)
   y.n <- length(y_labs)
   geom_point_color = my_tab$score
   #geom_point_border = my_tab$border

   g <-  ggplot(data = my_tab,
         aes(y=as.numeric(group), x=as.numeric(term), size=abs(as.numeric(size)))) +
         coord_fixed() +
         theme(panel.grid.major =  element_line(colour = 'lightsteelblue2'), panel.grid.minor =  element_line(colour = 'mistyrose3')) +
         theme(panel.background = element_rect(fill = 'mintcream', colour = NA)) +
         theme(panel.grid.major =  element_line(colour = NA), panel.grid.minor =  element_line(colour = 'black')) +
         theme(axis.title.y = element_text(size=12)) +
         theme(axis.title.x = element_text(size=12)) +
         theme(axis.text = element_text(size=9)) +
         theme(title = element_text(size=9)) +
         theme(strip.text = element_text(size=12)) +
         theme(title = element_text(size=12)) +
         theme(axis.text.y=element_text(angle = 0, hjust = 0, size = 12)) +
         theme(axis.text.x=element_text(angle = 45, hjust = 1, size = 10)) +
         theme(legend.text =  element_text(size=10)) +
         theme(legend.title =  element_text(size=10, face="bold")) +
         theme(legend.direction = "horizontal", legend.box = "vertical") +
         guides(size=guide_legend(title=size.title, direction="horizontal")) +
         labs(x="", y="" ) +
         scale_x_continuous(limits=c(0.5, x.n+0.5), breaks=c(1:x.n), labels=dlfoo2::label_wrap_mod(as.character(x_labs)), expand=c(0,0)) +
         scale_y_continuous(limits=c(0.5, y.n+0.5), breaks=c(1:y.n), labels=dlfoo2::label_wrap_mod(as.character(y_labs)), expand=c(0,0)) +
         scale_color_gradientn(name=score.title, colours = rev(score.pal), values=c(0,0.5,1)) +
         geom_point(aes(colour=my_tab$score), alpha=1) +
         scale_size_continuous(range=c(3,10)) +
         geom_point(aes(), fill=NA, colour=my_tab$border, alpha=1, pch=21, stroke=1)

         #scale_color_gradientn(name="-log10*(p-value)", colours = rev(score_pal), values=c(0,0.1,1)) +
         #scale_size_continuous(range=c(4,10)) +


   return(g)
  #ggsave(g, filename = "Correlation/TCGA/EA_plot_fgsea_tcga.pdf", width=18, height=10, useDingbats = F)
  }





#' Generic plot for plottng co-occurrence freuency matrix using Bubble blot.
#' @description Framction of y-axis group members (defined by my.columns[1]) that fall within x-axis groups, i.e. Sum of each plot row fraction equals 1.
#' @description v1.0 20191025.
#' @family Co-occurrence
#' @family plotwrap
#' @param my.tab, data frame with two grouping values to compare overlap. Should be factors to set order
#' @param my.columns two columns with group assignments that should be compared in frequency matrix. Y-axis and X-axis, respectively
#' @param color.key color key to match to primary column (y-axis)
#' @param order.by.clustering.x, let order of terms/groups be decided by clustering. Otherwise defined by levels
#' @param order.by.clustering.y, let order of terms/groups be decided by clustering. Otherwise defined by levels
#' @param reverse.x if to reverse x-axis
#' @param reverse.y if to reverse y-axis
#' @return ggplot object
#' @export
coOccurence_bubble <- function(
      my.tab=NULL,
      my.columns =NULL, color.key ="size",
      title.text=NULL,
      order.by.clustering.x=T, order.by.clustering.y=T,
      border.col = NA,
      #score.pal = c(rev(RColorBrewer::brewer.pal(9,"YlOrRd")[1:9]),rep("#FFFFCC",2),rep("white",5), rep("#FFFFCC",2), RColorBrewer::brewer.pal(9,"YlGn")[1:9]) ,
      reverse.x=F, reverse.y=F
   ){
   require(dplyr)
   require(tibble)
   require(tidyr)
   require(ggplot2)

   if(!all(sapply(my.columns, function(x)  x %in% colnames(my.tab)))) {message("colnames must include: ", paste(my.columns, sep = "  "))}

   my.tab <- droplevels(data.frame(my.tab))
   my.tab[,my.columns[1]] <- factor(my.tab[,my.columns[1]])
   my.tab[,my.columns[2]] <- factor(my.tab[,my.columns[2]])
   str(my.tab)

   my_ftab <- ftable(my.tab[,my.columns])
   my_ftab <- as.data.frame(t(apply(my_ftab, 1, function(x) round(x/sum(x)*100,1))))


   my_ftab <- my_ftab %>% mutate(y = levels(my.tab[,my.columns[1]]))
   my_ptab <- melt(my_ftab, id.vars = "y", variable.name = "x", value.name = "score")
   my_ptab <- my_ptab %>%
      mutate(y = factor(y, levels=levels(my.tab[,my.columns[1]]))) %>%
      mutate(x = factor(x, levels=levels(my.tab[,my.columns[2]]))) %>%
      mutate(score = score/100) %>%
      dplyr::filter(score>0) %>%
      mutate(score_scaled = (score+0.5)/1.5)

   x.n <- length(levels(my_ptab$x))
   y.n <- length(levels(my_ptab$y))

    # if to order rows/cols by cluster analysis
   if(any(order.by.clustering.x, order.by.clustering.y)){
      my_mat <- matrix(nrow=y.n, ncol=x.n, data=0, dimnames=list(levels(my_ptab$y), levels(my_ptab$x)))
      for(i in 1:nrow(my_ptab)){
         my_mat[as.numeric(my_ptab$y)[i], as.numeric(my_ptab$x[i])] <- my_ptab$score[i]
      }

      mydist=function(c) {amap::Dist(c,method="pearson")}
      myclust=function(c) {hclust(c,method="ward.D2")}
      hca.y <- my_mat %>% mydist %>% myclust %>% as.dendrogram
      hca.x <- t(my_mat) %>% mydist %>% myclust %>% as.dendrogram

      if(order.by.clustering.x){
        my_ptab <- my_ptab %>%  mutate(x = factor(x, levels=labels(hca.x)))
      }
      # if(!order.by.clustering.x){
      #   my_tab <- my_tab %>%  mutate(group = factor(group, levels=my_terms))
      # }
      if(order.by.clustering.y){
        my_ptab <- my_ptab %>%  mutate(y = factor(y, levels=labels(hca.y)))
      }

      # if(!order.by.clustering.y){
      #   my_tab <- my_tab %>%  mutate(group = factor(group, levels=my_groups)) ## ?!?!?"?!"??
      # }
      my_ptab <- my_ptab %>% arrange(y,x,-score)
   }

   if(reverse.x) my_ptab <- my_ptab %>% mutate(x = factor(x, levels=rev(levels(x)))) %>% arrange(y, x)
   if(reverse.y) my_ptab <- my_ptab %>% mutate(y = factor(y, levels=rev(levels(y)))) %>% arrange(y, x)



   x_labs <- levels(my_ptab$x)
   y_labs <- levels(my_ptab$y)
   # geom_point_color = my_tab$score
   #geom_point_border = my_tab$border

   my_color_key <- color.key[levels(my_ptab$y)]

   g <- ggplot(data = my_ptab,
         aes(y=as.numeric(y), x=as.numeric(x), size=abs(as.numeric(score)))) +
         coord_fixed() +
         theme(panel.grid.major =  element_line(colour = 'lightsteelblue2'), panel.grid.minor =  element_line(colour = 'mistyrose3')) +
         theme(panel.background = element_rect(fill = 'mintcream', colour = NA)) +
         theme(panel.grid.major =  element_line(colour = NA), panel.grid.minor =  element_line(colour = 'black')) +
         theme(axis.title.y = element_text(size=12)) +
         theme(axis.title.x = element_text(size=12)) +
         theme(axis.text = element_text(size=9)) +
         theme(title = element_text(size=9)) +
         theme(strip.text = element_text(size=12)) +
         theme(title = element_text(size=12)) +
         theme(axis.text.y=element_text(angle = 0, hjust = 0, size = 12)) +
         theme(axis.text.x=element_text(angle = 45, hjust = 1, size = 10)) +
         theme(legend.text =  element_text(size=10)) +
         theme(legend.title =  element_text(size=10, face="bold")) +
         theme(legend.direction = "horizontal", legend.box = "vertical") +
         ggtitle(label=title.text) +
         # guides(size=guide_legend(title=size.title, direction="horizontal")) +
         labs(x="", y="" ) +
         scale_x_continuous(limits=c(0.5, x.n+0.5), breaks=c(1:x.n), labels=dlfoo2::label_wrap_mod(as.character(x_labs)), expand=c(0,0)) +
         scale_y_continuous(limits=c(0.5, y.n+0.5), breaks=c(1:y.n), labels=dlfoo2::label_wrap_mod(as.character(y_labs)), expand=c(0,0)) +
         # scale_color_gradientn(name=score.title, colours = rev(score.pal), values=c(0,0.5,1)) +
         #geom_point(aes(colour=as.character(my_ptab$y), alpha=my_ptab$score_scaled)) +
         geom_point(aes(colour=as.character(y), alpha=score_scaled)) +
         scale_size_continuous(range=c(3,10)) +
         scale_color_manual(values=my_color_key, na.value=NA) +
         #geom_point(aes(), fill=NA, colour="black", alpha=my_ptab$score_scaled, pch=21, stroke=1)
         geom_point(aes(), fill=NA, colour="black", alpha=my_ptab$score_scaled, pch=21, stroke=1)

   return(g)

} # end coOccurrence bubble



#' function tsnePlot
#' @description Plots RCC and BLCA taxonomy annotations
#' @family tsne
#' @family PlotWrap
#' @aliases tsnePlotTax
#' @param tSNEres results list (tSNEres) from run.tsne-analysis
#' @param pdata defaults to pdata_panCanGex
#' @param annotation.column should be a column in the pdata object to for annotation.
#' @param perplex what perplexity to use. if null then plott all perplexes available in results table
#' @param my.samples if to highligt samples - character vector of sample names
#' @param plot.title Title of plot
#' @param color_key what color key to use for annotation. Detaults to dlfoo2::color_subtypes (but could be set to e.g. dlfoo2::color_annotations)
#' @param do.gradient set to FALSE if to manually override numeric check for gradient or manual fil
#' @param color.key.gradient what gradient.
#' @param gradient.nomalize TRUE or FALSE if to z normalize annotation value (x-mean(x)) sd(x)
#' @param gradient.squish.probs if to cap squish gradient to remove effects on color scale of outliers. vector of tow percent values. Uses quantile-probs. Defaults to c(0.05, 0.95),i.e. cap scale at 5 and 95 of values
#' @param my.size point size
#' @param my.alpha point alpha
#' @param my.stroke stroke width
#' @param my.stroke.color stroke line colors
#' @param showGuide if to show guide
#' @param coord.fixed if coord_fixed ggplot param
#' @param shape.key if to plot annotations in different shapes (connected to annotation)
#' @param zoom.coord list of x and y coordinates to specify/zoom plot coordinates
#' @return a ggplot object
#' @export
tsnePlot <- function(tSNEres, pdata=NULL, annotation.column="taxonomy_published", my_samples=NULL, plot.title="",
   perplex = NULL, my.size=2, my.alpha=0.7, my.stroke=0.05, my.stroke.color="gray10",
   color.key.annotation=dlfoo2::color_subtypes,
   do.gradient=T, palette.gradient=rev(dlfoo2::palette_gradients[["red2white2blue"]]), gadient.na.col=NA, gradient.normalize=F, gradient.squish.probs=c(0.05, 0.95),
   showGuide=T, coord.fixed=T, shape.key=dlfoo2::shape_subtypes, zoom.coord=NULL
   ){

   require(Rtsne)
   require(dplyr)
   require(tidyverse)
   require(tibble)
   require(ggrepel)

   if(is.null(pdata)) pdata <- dlfoo2data::pdata_panCanGex
   if(!(annotation.column %in% colnames(pdata))) stop("annotation.column not in selected pdata object")


   plot_df <- do.call("rbind", lapply(tSNEres, dlfoo2::tsne2df))
   available.perplexes <- sort(unique(plot_df$perplexity))
   if(!is.null(perplex)) if(!(perplex %in% available.perplexes)){stop("perplex not among available perplexes ",paste(available.perplexes,collapse = " "))}

   plot_df <- dplyr::left_join(plot_df, pdata)
   plot_df$annotation <- plot_df[, annotation.column]
   plot_df$annotation[plot_df$annotation==""] <- NA

   if(!is.null(my_samples)){
      if(length(which(my_samples %in% plot_df$sample_id)) > 0) {
         plot_df$plot_text <- NA
         plot_df$plot_text[plot_df$sample_id %in% my_samples] <- as.character(plot_df$sample_id[plot_df$sample_id %in% my_samples] )
      }
      if(length(which(my_samples %in% plot_df$sample_id)) < 1) {
         message("cant find any of the supplied sample_ids")
         my_samples <- NULL}

   }

   # if plot all perplexes
   if(is.null(perplex)){
      my_df = as_tibble(plot_df)
      do.perplex.all=T
   }else{
      my_df = as_tibble(subset(plot_df, perplexity==perplex))
      do.perplex.all=F
   }

   # if fill manual or fill gradient
   if(!is.numeric(my_df$annotation)){
      do.fill = "manual"
      my_df <- my_df %>% mutate(annotation = factor(annotation)) %>% arrange(annotation)
      color.key <- color.key.annotation[sort(levels(my_df$annotation))]
      }
   if(is.numeric(my_df$annotation) && do.gradient){
      do.fill = "gradient"
      if(gradient.normalize) my_df$annotation <- my_df$annotation - mean(my_df$annotation, na.rm=T) / sd(my_df$annotation, na.rm=T)
      #if(cap.outliers)
      #color.key.gradient = palette_gradientRamp(x = my_df$annotation, my.pal = palette.gradient, palette.sat = gradient.sat, na.col = gadient.na.col)
      color.key.gradient = palette.gradient
      my.gradient.limits <- quantile(my_df$annotation, probs=gradient.squish.probs, na.rm=T)
      }

   ## if do.shape
   do.shape = F
   if(class(my_df$annotation)=="character") my_df$annotation <- factor(my_df$annotation)
   if(any(levels(my_df$annotation) %in% names(shape.key))){
      do.shape = T
      #shape.key <- shape.key[sort(levels(my_df$annotation))]
      shape.key <- shape.key[(levels(my_df$annotation))]
      }

   # zoom coordinates
   do.zoom = FALSE
   if(!is.null(zoom.coord)){
      do.zoom=T
      if(class(zoom.coord)!="list") stop("zoom.coord must be list with named character slots x and y character slots of lengtht 2")
      if(length(zoom.coord)!=2) stop("zoom.coord must be list with named character slots x and y character slots of lengtht 2")
      if(!identical(sort(names(zoom.coord)), c("x","y"))) stop("zoom.coord must be list with named character slots x and y character slots of lengtht 2")
      }


   ## Plot standard taxonomy
   ## ----------------------
      g <- ggplot(my_df) +
            aes(x=dim1, y=dim2) +
            ggtitle((plot.title), label = annotation.column) +
            labs(x="tSNE dim1", y="tSNE dim2", colour='') +
            theme(axis.title.y = element_text(size=24)) +
            theme(axis.title.x = element_text(size=24)) +
            theme(legend.text = element_text(size=18)) +
            theme(title = element_text(size=24))



         if(do.shape & do.fill == "manual"){
            g <- g +
               geom_point(aes(fill=annotation, shape=annotation), colour=my.stroke.color, stroke=my.stroke, size=my.size, alpha=my.alpha) +
               scale_fill_manual(values=color.key, na.value=NA) +
               scale_shape_manual(values = as.numeric(shape.key), na.value = 4) +
               guides(color = guide_legend(override.aes = list(size=7)))
         }

         if(!do.shape & do.fill == "manual"){
             g <- g +
               geom_point(aes(fill=annotation), shape=21, stroke=my.stroke, colour=my.stroke.color, size=my.size, alpha=my.alpha, ) +
               scale_fill_manual(values=color.key, na.value=NA)+
               guides(fill=guide_legend(title=annotation.column))
               }
         if(!do.shape & do.fill == "gradient"){
               g <- g +
               geom_point(aes(fill=annotation), shape=21, stroke=my.stroke, colour=my.stroke.color, size=my.size, alpha=my.alpha) +
               # scale_fill_gradientn(colours=color.key.gradient, na.value=NA) +
               scale_fill_gradientn(colours=color.key.gradient, na.value=NA, limits=my.gradient.limits, oob=scales::squish) +
               guides(fill=guide_legend(title=annotation.column))
         }

         if(do.perplex.all){
            g <- g + facet_wrap(~perplexity, scales = "free")
         }

      if(showGuide==F){
            g <- g + guides(fill=FALSE, shape=FALSE, color=FALSE)
         }

      if(!is_null(my_samples)){
            g <- g + geom_label_repel(aes(label = plot_text),
                           box.padding   = 0.35,
                           point.padding = 0.5,
                           segment.color = 'grey50')
            }
      if(coord.fixed & !do.perplex.all) {g <- g + coord_fixed(ratio = 1)}

      if(do.zoom){
            g <- g + coord_cartesian(xlim = zoom.coord$x, ylim = zoom.coord$y)
               }

   return(g)
   } # end tsnePlot




#' Scatterpoint with annotations
#' @description Generic plot for tSNE (list from tsneWrapper analysis), umap (umap object), or PCA (prcomp object)
#' @family umap
#' @family tsne
#' @family pca
#' @family PlotWrap
#' @param x results object. list (tSNEres) from run.tsne-analysis, or umap (umap object), or PCA (prcomp object)
#' @param pdata data frame with, minimally, columns 'sample_id' and a column for annotation (defined by 'annotation.column'). defaults to pdata_panCanGex
#' @param annotation.column should be a column in the pdata object to for annotation.
#' @param plot.dims what 2Ds plot. PCA has multiple dims and sometimes umap. defaults to c(1,2)
#' @param my.samples if to highligt samples - character vector of sample names
#' @param plot.title Title of plot
#' @param color.key.annotation what color key to use for annotation. Detaults to dlfoo2::color_subtypes (but could be set to e.g. dlfoo2::color_annotations)
#' @param do.gradient set to FALSE if to manually override numeric check for gradient or manual fil
#' @param color.key.gradient what gradient.
#' @param gradient.nomalize TRUE or FALSE if to z normalize annotation value (x-mean(x)) sd(x)
#' @param gradient.squish.probs if to cap squish gradient to remove effects on color scale of outliers. vector of tow percent values. Uses quantile-probs. Defaults to c(0.05, 0.95),i.e. cap scale at 5 and 95 of values
#' @param my.size point size
#' @param my.alpha point alpha
#' @param my.stroke stroke width
#' @param my.stroke.color stroke line colors
#' @param showGuide if to show guide
#' @param coord.fixed if coord_fixed ggplot param
#' @param shape.key if to plot annotations in different shapes (connected to annotation)
#' @param zoom.coord list of x and y coordinates to specify/zoom plot coordinates
#' @return a ggplot object
#' @export
scatterPoint_annotationPlot <- function(x, pdata=NULL, annotation.column="tax_simp", my_samples=NULL, plot.title=NULL,
   my.size=2, my.alpha=0.7, my.stroke=0.05, my.stroke.color="gray10", plot.dims=c(1,2),
   # color.key.annotation=dlfoo2::color_subtypes,
   color.key.annotation=NULL,
   do.gradient=T, palette.gradient=rev(dlfoo2::palette_gradients[["red2white2blue"]]), gadient.na.col=NA, gradient.normalize=F, gradient.squish.probs=c(0.05, 0.95),
   showGuide=T, coord.fixed=T, shape.key=dlfoo2::shape_subtypes, zoom.coord=NULL
   ){

   require(umap)
   require(Rtsne)
   require(dplyr)
   require(tidyverse)
   require(tibble)
   require(ggrepel)

   if(is.null(pdata)) pdata <- dlfoo2data::pdata_panCanGex
   if(!"sample_id" %in% colnames(pdata)) stop("pdata object must contain 'sample_id' column")
   if(!(annotation.column %in% colnames(pdata))) stop("annotation.column not in selected pdata object")

   plot_type <- NULL


   if(class(x)=="umap"){
      message("Found UMAP results. Begin plotting ... ")
      res_df <- umap2df(x)
      plot_type <- "umap"
   }

   if(class(x)=="prcomp"){
      message("Found PCA results. Begin plotting ... ")
      plot_type <- "pca"
      res_df <- pca2df(x)
   }

    if(class(x)=="list"){
     if(!is.null(x$Y)){
      message("Found list with Y matrix slot  ... ")
      message("... assuming x is a tsne reults list. Begin plotting")
      plot_type <- "tsne"
      res_df <- tsne2df(x)
     }}

   if(class(x) == "data.frame"){
      message("Found data frame ... ")
      message("Lokking for  sample_id dim1, dim2 as columns")
      if(!all(c("sample_id","dim1","dim2") %in% colnames(x))) stop()
      plot_type <- "dataframe"
      res_df <- x[,c("sample_id","dim1","dim2")]
   }

   if(is.null(plot_type)) stop("x does not seem to be a prcomp, umap or tsne reusults list")


   ## plot_df: Join with pdata - check columns/ids
   u <- match(res_df$sample_id, pdata$sample_id)
   uu <- !is.na(u)
   if(all(!uu)) stop("cannot match sample_id for x (rownames) with sample_id column in pdata")

   message("... found ",length(which(uu==T))," matching saple_id in pdata out of ", length(which(uu==F)))
   plot_df <- dplyr::left_join(res_df, pdata)
   plot_df$annotation <- plot_df[, annotation.column]
   plot_df$annotation[plot_df$annotation==""] <- NA



   if(!is.null(my_samples)){
      if(length(which(my_samples %in% plot_df$sample_id)) > 0) {
         plot_df$plot_text <- NA
         plot_df$plot_text[plot_df$sample_id %in% my_samples] <- as.character(plot_df$sample_id[plot_df$sample_id %in% my_samples] )
      }
      if(length(which(my_samples %in% plot_df$sample_id)) < 1) {
         message("cant find any of the supplied sample_ids")
         my_samples <- NULL}
   }

   ## plot_df to tibble my_df
   my_df = as_tibble(plot_df)

   # if fill manual or fill gradient
   if(!is.numeric(my_df$annotation)){
      message("... ... annotation not numeric. setting fill to discrete colors defined by 'color.key.annotation'")
      do.fill = "manual"
      my_df <- my_df %>% mutate(annotation = factor(annotation)) %>% arrange(annotation)
      color.key <- color.key.annotation[sort(levels(my_df$annotation))]
      }
   if(is.numeric(my_df$annotation) && do.gradient){
      message("... ... annotation is numeric. setting fill to gradient as defined by 'color.key.gradient'")
      do.fill = "gradient"
      if(gradient.normalize) my_df$annotation <- my_df$annotation - mean(my_df$annotation, na.rm=T) / sd(my_df$annotation, na.rm=T)
      #if(cap.outliers)
      #color.key.gradient = palette_gradientRamp(x = my_df$annotation, my.pal = palette.gradient, palette.sat = gradient.sat, na.col = gadient.na.col)
      color.key.gradient = palette.gradient
      my.gradient.limits <- quantile(my_df$annotation, probs=gradient.squish.probs, na.rm=T)
      }

   ## if do.shape
   do.shape = F
   if(class(my_df$annotation)=="character") my_df$annotation <- factor(my_df$annotation)
   if(any(levels(my_df$annotation) %in% names(shape.key))){
      message("... ... setting do.shape to TRUE")
      do.shape = T
      #shape.key <- shape.key[sort(levels(my_df$annotation))]
      shape.key <- shape.key[(levels(my_df$annotation))]
      }

   # zoom coordinates
   do.zoom = FALSE
   if(!is.null(zoom.coord)){
      message("... ... 'zoom.chord' found. setting do.zoom to T")
      do.zoom=T
      if(class(zoom.coord)!="list") stop("zoom.coord must be list with named character slots x and y character slots of lengtht 2")
      if(length(zoom.coord)!=2) stop("zoom.coord must be list with named character slots x and y character slots of lengtht 2")
      if(!identical(sort(names(zoom.coord)), c("x","y"))) stop("zoom.coord must be list with named character slots x and y character slots of lengtht 2")
      }


   ##    Plot ggplot
   ## :::::::::::::::::::
      g <- ggplot(my_df) +
            aes(x=dim1, y=dim2) +
            ggtitle((plot.title), label = annotation.column) +
            labs(x="dim1", y="dim2", colour='') +
            theme(axis.title.y = element_text(size=24)) +
            theme(axis.title.x = element_text(size=24)) +
            theme(legend.text = element_text(size=18)) +
            theme(title = element_text(size=24))



         if(do.shape & do.fill == "manual"){
            g <- g +
               geom_point(aes(fill=annotation, shape=annotation), colour=my.stroke.color, stroke=my.stroke, size=my.size, alpha=my.alpha) +
               scale_fill_manual(values=color.key, na.value=NA) +
               scale_shape_manual(values = as.numeric(shape.key), na.value = 4) +
               guides(color = guide_legend(override.aes = list(size=7)))
         }

         if(!do.shape & do.fill == "manual"){
             g <- g +
               geom_point(aes(fill=annotation), shape=21, stroke=my.stroke, colour=my.stroke.color, size=my.size, alpha=my.alpha, ) +
               scale_fill_manual(values=color.key, na.value=NA)+
               guides(fill=guide_legend(title=annotation.column))
               }
         if(!do.shape & do.fill == "gradient"){
               g <- g +
               geom_point(aes(fill=annotation), shape=21, stroke=my.stroke, colour=my.stroke.color, size=my.size, alpha=my.alpha) +
               # scale_fill_gradientn(colours=color.key.gradient, na.value=NA) +
               scale_fill_gradientn(colours=color.key.gradient, na.value=NA, limits=my.gradient.limits, oob=scales::squish) +
               guides(fill=guide_legend(title=annotation.column))
         }


      if(showGuide==F){
            g <- g + guides(fill=FALSE, shape=FALSE, color=FALSE)
         }

      if(!is_null(my_samples)){
            g <- g + geom_label_repel(aes(label = plot_text),
                           box.padding   = 0.35,
                           point.padding = 0.5,
                           segment.color = 'grey50')
            }
      if(coord.fixed) {g = g + coord_fixed(ratio = 1)}

      if(do.zoom){g = g + coord_cartesian(xlim = zoom.coord$x, ylim = zoom.coord$y)}

   return(g)
   } # end scatterPoint_pdata









#' plot segments from granges file (typically Copy-number) using vanilla R-style plotting (lines)
#' @family PlotWrap
#' @family genomics
#' @family old functions
#' @family granges
#' @param x granges object
#' @param sample_name vector of sample ids
#' @param all_seqlengths if to run all available seqlengths (chromosomes)
#' @param cutoff cutoff to flag alterations
#' @param snp_cut number of snpt to flag alteration
#' @param lwd line width
#' @param mark.genes genomic ranges object. start postiion sets ablines (not completed)
#' @return plots to quartz
#' @export
segmentPlotter<-function(x, sample_name=NULL, all_seqlengths = F, cutoff=c(-0.2,0.2), snp_cut=150, lwd=2.5, mark.genes = NULL){
   require(GenomicRanges)
   require(RColorBrewer)
  	#if(!exists("centromere.limit.hg19")) load("~/RESOURCES/Annotation Files/centromere.limits.Rdata")
   data("centromere.limit.hg19", package = "dlfoo2")

      centromere.limit <- centromere.limit.hg19[!centromere.limit.hg19$chromosome%in%c(24),]
   centromere.limit$chromosome <- toChr(centromere.limit$chromosome)
   rownames(centromere.limit) <- centromere.limit$chromosome

   ## Check structure of x
	if(class(x)!="GRanges"){stop("x is not an 'GRanges' Object")} # Class = RangedData
   if(!length(seqlengths(x))){ stop("seqlengths must be defined")}
   my_ylim <- c(-2,2)
   score(x)[score(x) < my_ylim[1]]  <- my_ylim[1]
   score(x)[score(x) > my_ylim[2]]  <- my_ylim[2]

   min_seg <- min(score(x[x$n.snp>snp_cut]))
   if(min_seg>cutoff[1]) min_seg<-cutoff[1]
   max_seg <- max(score(x[x$n.snp>snp_cut]))
   if(max_seg<cutoff[2]) max_seg<-cutoff[2]

   if(all_seqlengths) {chroms = as.character(seqlevels(x))}
   if(!all_seqlengths) {chroms = as.character(unique(seqnames(x)))}
   my_seqlengths <- seqlengths(x)
   my_seqlengths <- my_seqlengths[chroms]
   genome_add <- data.frame(wgpos=c(0,cumsum(as.numeric(my_seqlengths))[-length(my_seqlengths)]), row.names = names(my_seqlengths))

   centromere.limit <- centromere.limit[chroms, -c(1)]
   #centromere.limit$genome <- genome_add

   # x <- dl_genomicPosCum(x, all_seqlengths = F)
   x1 <- as.numeric(start(x)) + as.numeric(genome_add[as.character(seqnames(x)), ])
   x2 <- as.numeric(end(x)) + as.numeric(genome_add[as.character(seqnames(x)), ])

   if(!all_seqlengths) my_xlim <- c(0, rev(cumsum( as.numeric(seqlengths(x)[as.character(unique(seqnames(x)))])) )[1])
   if(all_seqlengths) my_xlim <- c(0, rev(cumsum(as.numeric(seqlengths(x))))[1])

	plot.new()
	plot.window(xlim= my_xlim, ylim= my_ylim, xlab="", xaxt="n", yaxt="n", ylab="")
	rect(xleft = my_xlim[1], xright = my_xlim[2], ybottom = my_ylim[1], ytop = my_ylim[2], col="gray98", border = F)
	#rect(xleft = my_xlim[1], xright = my_xlim[2], ybottom = cutoff[2], ytop = cutoff[1], col="gray92", border = F)
	   p.args <- list(xleft=as.numeric(genome_add[,1]), xright=as.numeric(centromere.limit$maxParm+genome_add[,1]), ybottom=cutoff[1], ytop=cutoff[2],
	      col="gray85", border = NA, lwd=0.5)
	   do.call("rect", args = p.args)
	   q.args <- list(xleft=as.numeric(centromere.limit$minQarm+genome_add[,1]), xright=as.numeric(centromere.limit$maxQarm+genome_add[,1]), ybottom=cutoff[1], ytop=cutoff[2],
	      col="gray92", border = NA, lwd=0.5)
	   do.call("rect", args = q.args)

	if(min(score(x))<cutoff[1]) rect(xleft = my_xlim[1], xright = my_xlim[2], ybottom = min_seg, ytop = cutoff[1], col="lightcyan", border = F)
	if(max(score(x))>cutoff[2]) rect(xleft = my_xlim[1], xright = my_xlim[2], ybottom = cutoff[2], ytop = max_seg, col="mistyrose1", border = F)

	title(ylab="log2 CN", main=sample_name)


	abline(h=seq(-10,10,0.5),lty=3,lwd=0.75)
	abline(h=c(cutoff[1], cutoff[2]),lty=2, lwd=0.5)
	abline(v=c(genome_add[,1], my_xlim[2]), col="gray")
	axis(2,at=seq(-10,10,1),labels=seq(-10,10,1),las=1)
	#abline(h=c(0.25, -0.25), col=c("darkorange3","steelblue3"), lty=1, lwd=0.75)
   #abline(h=c(0), col="olivedrab1")

   lines_df <- data.frame(x1,x2, y1=score(x), y2=score(x), col="black", stringsAsFactors = F)
   my.pal <- rev(brewer.pal(11,"RdBu")[c(2:10)])
   my.pal <- c("dodgerblue2","wheat3","tomato")
   lines_df$col <- dlfoo2::heatcol_creator2(lines_df$y1, palette.sat = c(-0.5,0.5), my.pal = my.pal)$col

	plot.lines<-function(x){
		lines(c(x[1], x[2]), c(x[3], x[4]), col=x[5], lwd=lwd)
		}
	apply(lines_df, 1, plot.lines)

	axis(1, at=c(genome_add[,1]+seqlengths(x)[chroms]/2), labels=chroms, las=2, tick = "")
   abline(v=as.numeric(genome_add[,1]), lwd=1.5,col="gray")




	} ## END FUNCTION PLOT SEGMENTS





#' Plot Methylation Beta heatmap - standard
#' @description Version 1.0 20190124
#' @description Function to perform a standard HCA and heatmap plot using RCC (and bladder) taxonomy. Saves .png, .pdf and cluster trees (.rds) to file.
#' @family PlotWrap
#' @family heatmap
#' @family methylation
#' @family genoset
#' @param x.me GenoSet object with 'beta' slot
#' @param runName name of plot/analyses - subfolder
#' @param runFolder directory where to plot/save objects
#' @param clusterFeatures_k number of clusters for CpG feature tree
#' @param clusterSamples_k number of clusters for sample tree
#' @param amapDist Supplied to amap::Dist function. Can be one or more of distance methods. "euclidean", "maximum", "manhattan", "canberra", "binary", "pearson", "abspearson", "correlation", "abscorrelation", "spearman" or "kendall". Any unambiguous substring can be given.
#' @param palNo What palette (dl_heat)
#' @param palSat saturation levels for palette, defaults to '(0,1)'
#' @return plots to file
#' @export
#'
meHeatmapRCC <- function(
   x.me, runName, runFolder,
   clusterFeatures_k = 25,
   clusterSamples_k = 25,
   amapDist = c("pearson","euclidean","manhattan","correlation"),
   palNo=3, palSat=c(0, 1)
   ){

   require(dplyr)
   options(dplyr.width = Inf)
   require(tibble)
   require(tidyverse)
   require(tidyr)
   require(ggplot2)
   require(alluvial)
   require(ggalluvial)
   require(dendextend)
   require(RColorBrewer)
   if(class(x.me)!="GenoSet") stop("gset is not of class 'GenoSet'")
   amapDist <- match.arg(amapDist, choices = c("euclidean", "maximum", "manhattan", "canberra", "binary", "pearson", "abspearson", "correlation", "abscorrelation", "spearman", "kendall"), several.ok = T)

   #gsSinkGenosetDescription()
   mainPlotDir <- file.path(runFolder, runName)
   if(!file.exists(mainPlotDir)) dir.create(mainPlotDir, recursive = T)
   message("... plotting beta methylation values for:", paste(dim(x.me), c(" Features and "," sammples in folder:\n ")), mainPlotDir)

   ## Source ccpRCC source to get sample pData
   #source("~/PROJECTS/RCC_ccpRCC//R_SCRIPTS/R_ccpRCC_SOURCE.R")
   #message("... sourcing R_ccpRCC_SOURCE.R to get sample info - pdata_panCanGex")


   ## Select genes & samples for tree
   ## ------------------------s
   # my_es <- es_me450k[var_cpg, pdata_me$sample_id ] # 7562 features, 701 samples
   data("pdata_panCanGex", package="dlfoo2")
   pdata <- pdata_panCanGex

   message("... Retrieving curated sample data table for urinary tract TCGA samples 'pdata_panCanGex'")
   u <- match(colnames(x.me), pdata$sample_id)
   if(length(u[!is.na(u)])!=length(colnames(x.me))){
      if(all(is.na(u))) stop("no samples (colnames) are found in pdata_panCanGex")
      message("... ... removing samples not present in pdata file")
      message("... ... ... from: ", length(u))
      message("... ... ... to: ", length(which(!is.na(u))))
      x.me <- x.me[,!is.na(u)]
      u<-u[!is.na(u)]
      }
   pdata <- pdata[u,]
   pdata <- pdata %>% arrange(taxonomy_simplified, taxonomy_dl)
   x.me <- x.me[,pdata$sample_id]
   stopifnot(identical(pdata$sample_id, colnames(x.me)))
   my_mat <- x.me[,,"beta"]
   fdata <- rowRanges(x.me)


   for(i in 1:length(amapDist)){
      ## Create subfolder for Distanace method
      myDist <- amapDist[i]
      distPlotDir <- file.path(mainPlotDir, myDist)
      if(!dir.exists(distPlotDir)) dir.create(distPlotDir, recursive = T)
      message(" <- Using amap::Dist measure ", myDist, " ->")

      ## HCA - Features
      ## ---------------------------------------
      message("... HCA of features/rows")
      my_kf <- clusterFeatures_k
      my_colors <- rep(palette.dl.annot[[19]],5)[1:my_kf]
      mydist=function(c) {amap::Dist(c,method = myDist)}
      myclust=function(c) {hclust(c,method="ward.D2")}
      hca_rows <- my_mat %>% mydist %>% myclust %>% as.dendrogram

      message("... ... Plotting HCA tree for rows to file ")
      myFile <- paste0(distPlotDir, "/", runName,"_hcaTree.pdf")
      pdf(myFile, width=11.692, height=8.267)
         hca_rows %>%
            dendextend::set("labels_cex", value=0.6)  %>%
               #set("labels", value=my.labels) %>%
               dendextend::set("labels", value=NA) %>%
               dendextend::set("branches_k_color", value=my_colors, k=my_kf)  %>%
               dendextend::set("branches_lwd", value = 1) %>%
               plot(yaxt="n", horiz=T)
         dev.off()

      my_k_classes <- dendextend::cutree(hca_rows, k = my_kf, order_clusters_as_data = FALSE)
      #dend.list.split[[my.group]] <- cbind(fdata[names(my.k.classes),], my.k.classes)
      my_tab <- dplyr::left_join(data.frame(FeatureID=names(my_k_classes), row_cl=paste0("cl_",my_k_classes), stringsAsFactors = F),
         as.data.frame(fdata), by="FeatureID")
      my_tab$row_cl_col <- as.factor(my_tab$row_cl)
      levels(my_tab$row_cl_col) <- my_colors
      my_tab$row_cl_col  <- as.character(my_tab$row_cl_col)
      head(my_tab)
      myFile <- paste0(distPlotDir, "/", runName,"_hcaRows.txt")
      write.table(my_tab, file=myFile, sep="\t", row.names=F, quote=F)
   #
   # # Cluster Columns (Samples - ALL samples)
   # ## ----------------
      message("... HCA of samples/columns")
      my_ks <- clusterSamples_k
      my_colors <- rep(palette.dl.annot[[19]],5)[1:my_ks]
      #hca_rows <- hca(t(my_mat), plot=F)$hclust
      mydist=function(c) {amap::Dist(c,method=myDist)}
      myclust=function(c) {hclust(c,method="ward.D2")}
      hca_cols <- t(my_mat) %>% mydist %>% myclust %>% as.dendrogram
      my_k_classes_cols <- dendextend::cutree(hca_cols, k = my_ks, order_clusters_as_data = FALSE)

   ## Alluvial graph
   ## --------------
      message("... ... Plotting Alluvial for HCA sample tree")
      hca_k_order <- c(paste0("cl_",1:my_ks))
      my_t <- pdata %>%
         dplyr::select(sample_id, taxonomy_simplified) %>%
         mutate(HCA_k=paste0("cl_",my_k_classes_cols[pdata$sample_id])) %>%
         mutate(HCA_k = factor(HCA_k, levels=hca_k_order))

      #ggalluvial::is_alluvia_form(my_t, axes = 1:6, silent = TRUE)
      #my_t_melt <- reshape2::melt(my_t, id.vars=c("sample_id","taxonomy"), measure.vars="HCA_k")
      my_a <- my_t %>% group_by(taxonomy_simplified, HCA_k) %>%
         summarise(n = n())
      myFile <- paste0(distPlotDir, "/", runName,"_TaxAlluvial.pdf")
      pdf(myFile)
         alluvial(my_a[,c("taxonomy_simplified","HCA_k")], freq=my_a$n,
            cex=0.5,
           #gap.width = 0.25,
            col= color_subtypes[match(my_a$taxonomy_simplified, names(color_subtypes))]
            #ordering = list(
            #   match(my_a$taxonomy_simplified, my_a$taxonomy_simplified),
            #   match(my_a$HCA_k, hca_k_order))
         )
         dev.off()

      myFile <- paste0(distPlotDir, "/", runName,"_hcaColumns.txt")
      write.table(my_t, file=myFile, sep="\t", row.names=F, quote=F)



      ##    Row annotations (rcol, heatmap.3)
      ## ------------------------------
      rownames(my_tab) <- my_tab$FeatureID
      row_annot <- data.frame(row_cl=my_tab$row_cl_col)
      str(row_annot)
      rownames(row_annot) <- rownames(my_tab)
      rcol <- as.matrix(t(row_annot))[,rownames(my_mat)]
      if(class(rcol)=="character"){
         rcol=matrix(data=rcol, nrow=1, dimnames = list(c(),names(rcol)))
         }
      stopifnot(identical(rownames(my_mat), colnames(rcol)))


      ##    column annotations (ccol, heatmap.3)
      ## ------------------------------
      data('color_subtypes', package='dlfoo2')
      data('shape_subtypes', package='dlfoo2')
      identical(colnames(my_mat), pdata$sample_id)
      annot <- pdata$taxonomy_simplified
      levels(annot) <- color_subtypes[match(levels(annot), names(color_subtypes))]

      annot1 <- as.factor(pdata$taxonomy_dl)
      levels(annot1) <- color_subtypes[match(levels(annot1), names(color_subtypes))]

      annot2 <- as.factor(my_k_classes_cols[colnames(my_mat)])
      levels(annot2) <- my_colors[1:length(levels(annot2))]

      ccol <- cbind(as.character(annot), as.character(annot1), as.character(annot2))
      colnames(ccol) <- c("taxonomy_simp","taxonomy","cluster")
      rownames(ccol) <- colnames(my_mat)
      head(ccol)


      ## Plot PDFs - heatmap.3
      ## ---------------------------
      #my_heatpal <- palette_dl_heatmap2(my_mat, pal.vec= colorRamps::blue2yellow(21), palette.sat = c(0, 1))
      my_heatpal <- palette_dl_heatmap2(my_mat, pal.no = palNo, palette.sat = palSat)

      myPdf <- paste(runName,"_heatmap2way_",paste(dim(my_mat),collapse = "x"),".pdf", sep = "")
      myFile <- file.path(distPlotDir, myPdf)
     # my_file <- "HCA_Me_change/HCA_Me_change_NvsT/HCA_Me_N_vs_T_heatmap_PearsonWards_9984x281.pdf"
      pdf(file=myFile, height = 9, width = 7)
         heatmap.3(my_mat,
            col=my_heatpal$legend$col,
            breaks=my_heatpal$legend$legend,
            trace=c("none"),
            notecol = NA,
            Rowv = as.dendrogram(hca_rows %>% set("branches_k_color", value=my_colors, k=my_kf)),
            Colv= as.dendrogram(hca_cols %>% set("branches_k_color", value=my_colors[1:my_ks], k=my_ks)),
            dendrogram = "both",
            labRow = NA,
            labCol = NA,

            RowSideColors = rcol,
            ColSideColors = ccol,
            RowSideColorsSize=5
            #cex=0.7
            )
         dev.off()
      myPdf <- paste(runName,"_heatmap2way_",paste(dim(my_mat),collapse = "x"),".png", sep = "")
      myFile <- file.path(distPlotDir, myPdf)
     # my_file <- "HCA_Me_change/HCA_Me_change_NvsT/HCA_Me_N_vs_T_heatmap_PearsonWards_9984x281.pdf"
      png(file=myFile, width = 596, height =  842)
         heatmap.3(my_mat,
            col=my_heatpal$legend$col,
            breaks=my_heatpal$legend$legend,
            trace=c("none"),
            notecol = NA,
            Rowv = as.dendrogram(hca_rows %>% set("branches_k_color", value=my_colors, k=my_kf)),
            Colv= as.dendrogram(hca_cols %>% set("branches_k_color", value=my_colors[1:my_ks], k=my_ks)),
            dendrogram = "both",
            labRow = NA,
            labCol = NA,

            RowSideColors = rcol,
            ColSideColors = ccol,
            RowSideColorsSize=5
            #cex=0.7
            )
         dev.off()


      ## 2-way
      ## ------
      myPdf <- paste(runName,"_heatmap1way_",paste(dim(my_mat),collapse = "x"),".pdf", sep = "")
      myFile <- file.path(distPlotDir, myPdf)
         pdf(file=myFile, height = 9, width = 7)
         heatmap.3(my_mat,
            col=my_heatpal$legend$col,
            breaks=my_heatpal$legend$legend,
            trace=c("none"),
            notecol = NA,
            Rowv = as.dendrogram(hca_rows %>% set("branches_k_color", value=my_colors, k=my_kf)),
            Colv= NA,
            dendrogram = "row",
            labRow = NA,
            labCol = NA,

            RowSideColors = rcol,
            ColSideColors = ccol,
            RowSideColorsSize=5
            #cex=0.7
            )
         dev.off()
      myPdf <- paste(runName,"_heatmap1way_",paste(dim(my_mat),collapse = "x"),".png", sep = "")
      myFile <- file.path(distPlotDir, myPdf)
         png(file=myFile,  width = 596, height =  842)
         heatmap.3(my_mat,
            col=my_heatpal$legend$col,
            breaks=my_heatpal$legend$legend,
            trace=c("none"),
            notecol = NA,
            Rowv = as.dendrogram(hca_rows %>% set("branches_k_color", value=my_colors, k=my_kf)),
            Colv= NA,
            dendrogram = "row",
            labRow = NA,
            labCol = NA,

            RowSideColors = rcol,
            ColSideColors = ccol,
            RowSideColorsSize=5
            #cex=0.7

            )
         dev.off()


      } #end amap dist i

}


#' Plot Methylation Beta heatmap with matched GEX heatmao
#' @description Version 1.0 20190124
#' @description Function to perform a standard HCA and heatmap plot of meatched me/gex data using RCC (and bladder) taxonomy. Saves .png, .pdf and cluster trees (.rds) to file. Requires genoset with ECR tab as input
#' @description ECR tab describes Me-Gex pairs. And lists correlaions, Extension (and replaces meHeatmapRCC but adds also a GEX heatmap
#' @family plot wraps
#' @family heatmap
#' @family methylation
#' @family genoset
#' @family eset
#' @param x.me methylation GenoSet with 'beta' slot
#' @param x.eset expressionSet with gex data
#' @param runName name of plot/analyses - subfolder
#' @param runFolder directory where to plot/save objects
#' @param clusterFeatures_k number of clusters for CpG feature tree
#' @param clusterSamples_k number of clusters for sample tree
#' @param amapDist Supplied to amap::Dist function. Can be one or more of distance methods. "euclidean", "maximum", "manhattan", "canberra", "binary", "pearson", "abspearson", "correlation", "abscorrelation", "spearman" or "kendall". Any unambiguous substring can be given.
#' @param palNo What palette (dl_heat)
#' @param palNoGex numer (id) of palette to use. uses dlfoo2::palette_dl_heatmap2
#' @param palSat saturation levels for beta levels. defaults to '0,1'
#' @param paSatGex saturatuon levels for gene expression. defaults to '0,3'
#' @return plots to file
#' @export
#'
me2gexHeatmapRCC <- function(
   x.me, x.gex,
   runName, runFolder,
   clusterFeatures_k = 25,
   clusterSamples_k = 25,
   amapDist = c("pearson","euclidean","manhattan","correlation"),
   palNo=3, palSat=c(0, 1), palNoGex=1, palSatGex=c(-3,3)
   ){

   require(dplyr)
   options(dplyr.width = Inf)
   require(tibble)
   require(tidyverse)
   require(tidyr)
   require(ggplot2)
   require(alluvial)
   require(ggalluvial)
   require(RColorBrewer)

   if(class(x.me)!="GenoSet") stop("gset is not of class 'GenoSet'")
   if(class(x.gex)!="ExpressionSet") stop("gset is not of class 'ExpressionSet'")
   if(is.null(metadata(x.me)[["meFilterECR_tab"]])) stop("x.me GenoSet must contain metadata slot 'meFilterECR_tab'")
   if(!identical(colnames(x.me), sampleNames(x.gex))) stop("x.me coluns and sampleNames in x.gex must be identical")
   amapDist <- match.arg(amapDist, choices = c("euclidean", "maximum", "manhattan", "canberra", "binary", "pearson", "abspearson", "correlation", "abscorrelation", "spearman", "kendall"), several.ok = T)


   #gsSinkGenosetDescription()
   mainPlotDir <- file.path(runFolder, runName)
   if(!file.exists(mainPlotDir)) dir.create(mainPlotDir, recursive = T)
   message("... plotting beta methylation values for:", paste(dim(x.me), c(" Features and"," Sammples"), " in folder: ", mainPlotDir))

   # ## Source ccpRCC source to get sample pData
   # #message("... sourcing R_ccpRCC_SOURCE.R to get sample info - pdata_panCanGex")
   # # source("~/PROJECTS/RCC_ccpRCC//R_SCRIPTS/R_ccpRCC_SOURCE.R")
   # data("pdata_panCanGex")
   #
   # ## Select genes & samples for tree
   # ## ------------------------s
   # u <- match(colnames(x.me),pdata_panCanGex$sample_id)
   # if(length(u[!is.na(u)])!=length(colnames(x.me))) stop("not all samples (colnames) are in pdata_panCanGex")
   # pdata <- pdata_panCanGex[u,]
   # stopifnot(identical(pdata$sample_id, colnames(x.me)))
   #
   # ## SAMPLE INFO and ORDER
   # ## ------------------
   # # table(pdata$taxonomy_dl)
   # my_groups <- c("crtx","crtmed","med","infl1","infl2","ccpRCC","Ch-e","chONC", "pONC","CC-e.2","CC-e.1","CC-e.3",
   #    "P-e.1a","P-e.1b","P-e.2","P.CIMP-e","CIMP_o","rcc_UroMesInf",
   #    "outlier",
   #    "uro","Urobasal","GU","Basal/SCCL","Mes-Inf","Inf_other","SC/NE","blca_normal")
   #
   # pdata <- pdata %>%
   #    mutate(taxonomy_dl = factor(taxonomy_dl, levels=my_groups)) %>%
   #    mutate(taxonomy_dl = droplevels(taxonomy_dl)) %>%
   #    arrange(taxonomy_dl) %>%
   #    mutate(taxonomy_simplified = factor(taxonomy_simplified, levels = c("ccpRCC","chRCC","chONC","pONC","ccRCC","pRCC_e1","pRCC_e2","CIMP_p","CIMP_o","rcc_UroMesInf","uro","outlier")))
   # x.me <- x.me[,pdata$sample_id]


   data("pdata_panCanGex", package="dlfoo2")
   pdata <- pdata_panCanGex

   message("... Retrieving curated sample data table for urinary tract TCGA samples 'pdata_panCanGex'")
   u <- match(colnames(x.me), pdata$sample_id)
   if(length(u[!is.na(u)])!=length(colnames(x.me))){
      if(all(is.na(u))) stop("no samples (colnames) are found in pdata_panCanGex")
      message("... ... removing samples not present in pdata file")
      message("... ... ... from: ", length(u))
      message("... ... ... to: ", length(which(!is.na(u))))
      x.me <- x.me[,!is.na(u)]
      u<-u[!is.na(u)]
      }
   pdata <- pdata[u,]
   pdata <- pdata %>% arrange(taxonomy_simplified, taxonomy_dl)
   x.me <- x.me[,pdata$sample_id]
   stopifnot(identical(pdata$sample_id, colnames(x.me)))
   my_mat <- x.me[,,"beta"]
   fdata <- rowRanges(x.me)
   x.gex <- x.gex[,pdata$sample_id]


   ## Get ECR table used to split pos/neg correlated Me Features
   ## -------------------
   ecr_tab <- metadata(x.me)[["meFilterECR_tab"]]


   ## Loop all distance measures
   ## --------------------
   for(i in 1:length(amapDist)){
      ## Create subfolder for Distanace method
      myDist <- amapDist[i]
      distPlotDir <- file.path(mainPlotDir, myDist)
      if(!dir.exists(distPlotDir)) dir.create(distPlotDir, recursive = T)
      message(" <- Using amap::Dist measure ", myDist, " ->")

      ## loop Pos/neg ECR pairs separately
      for(ii in c("pos","neg")){
            x.me.ii <- x.me[unique(ecr_tab %>% dplyr::filter(ECR==ii) %>% pull(FeatureID))]
            my_mat <- x.me.ii[,,"beta"]
            fdata <- rowRanges(x.me.ii)

         ## HCA - Features
         ## ---------------------------------------
            message("... HCA of features/rows")
            my_kf <- clusterFeatures_k
            my_colors <- rep(palette.dl.annot[[19]],5)[1:my_kf]
            mydist=function(c) {amap::Dist(c,method = myDist)}
            myclust=function(c) {hclust(c,method="ward.D2")}
            hca_rows <- my_mat %>% mydist %>% myclust %>% as.dendrogram

            message("... ... Plotting HCA tree for rows to file ")
            myFile <- paste0(distPlotDir, "/", runName,"_hcaTree_ECR_",ii,"_k",my_kf,".pdf")
            pdf(myFile, width=11.692, height=8.267)
               hca_rows %>%
                  dendextend::set("labels_cex", value=0.6)  %>%
                     #set("labels", value=my.labels) %>%
                     dendextend::set("labels", value=NA) %>%
                     dendextend::set("branches_k_color", value=my_colors, k=my_kf)  %>%
                     dendextend::set("branches_lwd", value = 1) %>%
                     plot(yaxt="n", horiz=T)
               dev.off()

            my_k_classes <- dendextend::cutree(hca_rows, k = my_kf, order_clusters_as_data = FALSE)
            #dend.list.split[[my.group]] <- cbind(fdata[names(my.k.classes),], my.k.classes)
            my_tab <- dplyr::left_join(data.frame(FeatureID=names(my_k_classes), row_cl=paste0("cl_",my_k_classes), stringsAsFactors = F),
               as.data.frame(fdata), by="FeatureID")
            my_tab$row_cl_col <- as.factor(my_tab$row_cl)
            levels(my_tab$row_cl_col) <- my_colors
            my_tab$row_cl_col  <- as.character(my_tab$row_cl_col)
            myFile <- paste0(distPlotDir, "/", runName,"_hcaRows_ECR_",ii,"_k",my_kf,".txt")
            write.table(my_tab, file=myFile, sep="\t", row.names=F, quote=F)

         # # Cluster Columns (Samples - ALL samples)
         # ## ----------------
            message("... HCA of samples/columns")
            my_ks <- clusterSamples_k
            my_colors <- rep(palette.dl.annot[[19]],5)[1:my_ks]
            #hca_rows <- hca(t(my_mat), plot=F)$hclust
            mydist=function(c) {amap::Dist(c,method=myDist)}
            myclust=function(c) {hclust(c,method="ward.D2")}
            hca_cols <- t(my_mat) %>% mydist %>% myclust %>% as.dendrogram
            my_k_classes_cols <- dendextend::cutree(hca_cols, k = my_ks, order_clusters_as_data = FALSE)



         ## Alluvial graph
         ## --------------
            message("... ... Plotting Alluvial for HCA sample tree")
            hca_k_order <- c(paste0("cl_",1:my_ks))
            my_t <- pdata %>%
               dplyr::select(sample_id, taxonomy_dl) %>%
               mutate(HCA_k = paste0("cl_", my_k_classes_cols[pdata$sample_id])) %>%
               mutate(HCA_k = factor(HCA_k, levels=hca_k_order))

            my_a <- my_t %>% group_by(taxonomy_dl, HCA_k) %>%
               summarise(n = n())
            myFile <- paste0(distPlotDir, "/", runName,"_TaxAlluvial_ECR_",ii,"_k",my_ks,".pdf")
            pdf(myFile)
               alluvial(my_a[,c("taxonomy_dl","HCA_k")], freq=my_a$n,
                  cex=0.5,
                 #gap.width = 0.25,
                  col= color_subtypes[match(my_a$taxonomy_dl, names(color_subtypes))]
                  #ordering = list(
                  #   match(my_a$taxonomy_simplified, my_a$taxonomy_simplified),
                  #   match(my_a$HCA_k, hca_k_order))
               )
               dev.off()

            myFile <- paste0(distPlotDir, "/", runName,"_hcaColumns_ECR_",ii,"_k",my_ks,".txt")
            write.table(my_t, file=myFile, sep="\t", row.names=F, quote=F)



         ##    Row annotations (rcol, heatmap.3)
         ## ------------------------------
         rownames(my_tab) <- my_tab$FeatureID
         row_annot <- data.frame(row_cl=my_tab$row_cl_col)
         str(row_annot)
         rownames(row_annot) <- rownames(my_tab)
         rcol <- as.matrix(t(row_annot))[,rownames(my_mat)]
         if(class(rcol)=="character"){
            rcol=matrix(data=rcol, nrow=1, dimnames = list(c(),names(rcol)))
            }
         stopifnot(identical(rownames(my_mat), colnames(rcol)))


         ##    column annotations (ccol, heatmap.3)
         ## ------------------------------
         data("color_subtypes", package="dlfoo2")
         identical(colnames(my_mat), pdata$sample_id)
         annot <- as.factor(pdata$taxonomy_simplified)
         levels(annot) <- color_subtypes[match(levels(annot), names(color_subtypes))]

         annot1 <- as.factor(pdata$taxonomy_dl)
         levels(annot1) <- color_subtypes[match(levels(annot1), names(color_subtypes))]

         annot2 <- as.factor(my_k_classes_cols[colnames(my_mat)])
         levels(annot2) <- my_colors[1:length(levels(annot2))]

         ccol <- cbind(as.character(annot), as.character(annot1), as.character(annot2))
         colnames(ccol) <- c("taxonomy_simp","taxonomy","cluster")
         rownames(ccol) <- colnames(my_mat)
         head(ccol)


         ## Plot PDFs - heatmap.3
         ## ---------------------------
         #my_heatpal <- palette_dl_heatmap2(my_mat, pal.vec= colorRamps::blue2yellow(21), palette.sat = c(0, 1))
         my_heatpal <- palette_dl_heatmap2(my_mat, pal.no = palNo, palette.sat = palSat)

         myPdf <- paste(runName,"_heatmap2way_",paste(dim(my_mat),collapse = "x"),"_ECR_",ii,"_k",my_ks,".pdf", sep = "")
         myFile <- file.path(distPlotDir, myPdf)
        # my_file <- "HCA_Me_change/HCA_Me_change_NvsT/HCA_Me_N_vs_T_heatmap_PearsonWards_9984x281.pdf"
         pdf(file=myFile, height = 9, width = 7)
            heatmap.3(my_mat,
               col=my_heatpal$legend$col,
               breaks=my_heatpal$legend$legend,
               trace=c("none"),
               notecol = NA,
               Rowv = as.dendrogram(hca_rows %>% set("branches_k_color", value=my_colors, k=my_kf)),
               Colv= as.dendrogram(hca_cols %>% set("branches_k_color", value=my_colors[1:my_ks], k=my_ks)),
               dendrogram = "both",
               labRow = NA,
               labCol = NA,

               RowSideColors = rcol,
               ColSideColors = ccol,
               RowSideColorsSize=5
               #cex=0.7
               )
            dev.off()

         png(file= gsub(".pdf",".png",myFile), width = 596, height =  842)
            heatmap.3(my_mat,
               col=my_heatpal$legend$col,
               breaks=my_heatpal$legend$legend,
               trace=c("none"),
               notecol = NA,
               Rowv = as.dendrogram(hca_rows %>% set("branches_k_color", value=my_colors, k=my_kf)),
               Colv= as.dendrogram(hca_cols %>% set("branches_k_color", value=my_colors[1:my_ks], k=my_ks)),
               dendrogram = "both",
               labRow = NA,
               labCol = NA,

               RowSideColors = rcol,
               ColSideColors = ccol,
               RowSideColorsSize=5
               #cex=0.7
               )
            dev.off()


         ## 1-way
         ## ------
         myPdf <- paste(runName,"_heatmap1way_",paste(dim(my_mat),collapse = "x"),"_ECR_",ii,"_k",my_ks,".pdf", sep = "")
         myFile <- file.path(distPlotDir, myPdf)
            pdf(file=myFile, height = 9, width = 7)
            heatmap.3(my_mat,
               col=my_heatpal$legend$col,
               breaks=my_heatpal$legend$legend,
               trace=c("none"),
               notecol = NA,
               Rowv = as.dendrogram(hca_rows %>% set("branches_k_color", value=my_colors, k=my_kf)),
               Colv= NA,
               dendrogram = "row",
               labRow = NA,
               labCol = NA,

               RowSideColors = rcol,
               ColSideColors = ccol,
               RowSideColorsSize=5
               #cex=0.7
               )
            dev.off()
         png(file=gsub(".pdf",".png",myFile),  width = 596, height =  842)
            heatmap.3(my_mat,
               col=my_heatpal$legend$col,
               breaks=my_heatpal$legend$legend,
               trace=c("none"),
               notecol = NA,
               Rowv = as.dendrogram(hca_rows %>% set("branches_k_color", value=my_colors, k=my_kf)),
               Colv= NA,
               dendrogram = "row",
               labRow = NA,
               labCol = NA,

               RowSideColors = rcol,
               ColSideColors = ccol,
               RowSideColorsSize=5
               #cex=0.7

               )
            dev.off()


         ## Plot GEX heatmaps - 1:n rows & n:n
         ## -----------------

            #  All Me-Gex Pairs for the pos/neg
         my_tab_gex <- left_join(my_tab, ecr_tab %>% dplyr::filter(ECR==ii), by="FeatureID")
         myFile <- paste0(distPlotDir, "/", runName,"_hcaRows_GEX_ECR_",ii,"_k",my_kf,".txt")
         write.table(my_tab_gex, file=myFile, sep="\t", row.names=F, quote=F)

            #str(my_tab_gex) # data.frame':	1836 obs
         my_mat_gex <- exprs(x.gex)[my_tab_gex$ENSG,]
         my_mat_gex <- t(apply(my_mat_gex, 1, function(x) x-mean(x)))
            #str(my_mat_gex)

         # rownames(my_tab_gex) <- my_tab_gex$ENSG
         row_annot <- data.frame(row_cl=my_tab_gex$row_cl_col)
         str(row_annot)
         rownames(my_mat_gex) <- rownames(row_annot)
         #rcol <- as.matrix(t(row_annot))[,rownames(my_tab_gex)]
         rcol <- as.matrix(t(row_annot))
         if(class(rcol)=="character"){
            rcol=matrix(data=rcol, nrow=1, dimnames = list(c(),names(rcol)))
         }
         colnames(rcol) <- rownames(row_annot)
         #stopifnot(identical(rownames(my_mat_gex), colnames(rcol)))


         ## Plot PDFs - heatmap.3
         ## ---------------------------
         #my_heatpal <- palette_dl_heatmap2(my_mat, pal.vec= colorRamps::blue2yellow(21), palette.sat = c(0, 1))
         my_heatpal <- palette_dl_heatmap2(my_mat_gex, pal.no = palNoGex, palette.sat = palSatGex)

         myPdf <- paste(runName,"_heatmap2way_GEX_",paste(dim(my_mat_gex),collapse = "x"),"_ECR_",ii,"_k",my_ks,".pdf", sep = "")
         myFile <- file.path(distPlotDir, myPdf)
         pdf(file=myFile, height = 9, width = 7)
            heatmap.3(my_mat_gex,
               col=my_heatpal$legend$col,
               breaks=my_heatpal$legend$legend,
               trace=c("none"),
               notecol = NA,
               Rowv = NA,
               Colv= 8,
               dendrogram = "col",
               labRow = NA,
               labCol = NA,

               RowSideColors = rcol,
               ColSideColors = ccol,
               RowSideColorsSize=5
               #cex=0.7
               )
            dev.off()

         png(file= gsub(".pdf",".png",myFile), width = 596, height =  842)
            heatmap.3(my_mat_gex,
               col=my_heatpal$legend$col,
               breaks=my_heatpal$legend$legend,
               trace=c("none"),
               notecol = NA,
               Rowv = NA,
               Colv= as.dendrogram(hca_cols %>% set("branches_k_color", value=my_colors[1:my_ks], k=my_ks)),
               dendrogram = "col",
               labRow = NA,
               labCol = NA,

               RowSideColors = rcol,
               ColSideColors = ccol,
               RowSideColorsSize=5
               #cex=0.7
               )
            dev.off()


         ## 1-way
         ## ------
         myPdf <- paste(runName,"_heatmap1way_GEX_",paste(dim(my_mat),collapse = "x"),"_ECR_",ii,"_k",my_ks,".pdf", sep = "")
         myFile <- file.path(distPlotDir, myPdf)
            pdf(file=myFile, height = 9, width = 7)
            heatmap.3(my_mat_gex,
               col=my_heatpal$legend$col,
               breaks=my_heatpal$legend$legend,
               trace=c("none"),
               notecol = NA,
               Rowv = NA,
               Colv= NA,
               dendrogram = "none",
               labRow = NA,
               labCol = NA,

               RowSideColors = rcol,
               ColSideColors = ccol,
               RowSideColorsSize=5
               #cex=0.7
               )
            dev.off()
         png(file=gsub(".pdf",".png",myFile),  width = 596, height =  842)
            heatmap.3(my_mat_gex,
               col=my_heatpal$legend$col,
               breaks=my_heatpal$legend$legend,
               trace=c("none"),
               notecol = NA,
               Rowv = NA,
               Colv= NA,
               dendrogram = "none",
               labRow = NA,
               labCol = NA,

               RowSideColors = rcol,
               ColSideColors = ccol,
               RowSideColorsSize=5
               #cex=0.7

               )
            dev.off()


         } #end amap dist i
      } # end ii pos/neg Me-Gex corrs
   }









#' function Plot GEX eset heatmaps to file using standard Taxonomy.
#' @description Performs HCA and heatmap plot using different distances
#' @description Version 1.0 20190129
#' @family plot wraps
#' @family heatmap
#' @family eset
#' @family gene expression
#' @aliases gexHeatmapRCC
#' @param x.eset gene expression eSet with expression values to be plotted. Input shpuld be Sorted and Filtered!
#' @param runName name of plot/analyses - subfolder
#' @param runFolder directory where to plot/save objects
#' @param clusterFeatures_k How many subclusters should the gene HCA tree be split into
#' @param clusterSamples_k  How many subclusters should the sample HCA tree be split into
#' @param amapDist Supplied to amap::Dist function. Can be one or more of distance methods. "euclidean", "maximum", "manhattan", "canberra", "binary", "pearson", "abspearson", "correlation", "abscorrelation", "spearman" or "kendall". Any unambiguous substring can be given.
#' @param palNo What palette (dl_heat)
#' @param palSat saturation levels for palette. defaults to -3,3
#' @return heatmaps to file
#' @export
#'
esetHeatmapRCC <- function(
   x.eset, runName, runFolder,
   clusterFeatures_k = 25,
   clusterSamples_k = 25,
   amapDist = c("pearson","euclidean","manhattan","correlation"),
   palNo=3, palSat=c(-3, 3)
   ){

   message("... This function requires ccpRCC SOURCE file")


   require(dplyr)
   options(dplyr.width = Inf)
   require(tibble)
   require(tidyverse)
   require(tidyr)
   require(ggplot2)

   require(alluvial)
   require(ggalluvial)

   if(class(x.eset)!="ExpressionSet") stop("gset is not of class 'ExpressionSet'")
   amapDist <- match.arg(amapDist, choices = c("euclidean", "maximum", "manhattan", "canberra", "binary", "pearson", "abspearson", "correlation", "abscorrelation", "spearman", "kendall"), several.ok = T)


   #gsSinkGenosetDescription()
   mainPlotDir <- file.path(runFolder, runName)
   if(!file.exists(mainPlotDir)) dir.create(mainPlotDir, recursive = T)
   message("... plotting GEX values for:", paste(dim(x.eset), c(" Features and "," Sammples")), " \nin folder: ", mainPlotDir)

   ## Source ccpRCC source to get sample pData
   #message("... sourcing R_ccpRCC_SOURCE.R to get sample info - pdata_panCanGex")
   #source("~/PROJECTS/RCC_ccpRCC//R_SCRIPTS/R_ccpRCC_SOURCE.R")

   # ## Select genes & samples for tree
   # ## ------------------------s
   # # my_es <- es_me450k[var_cpg, pdata_me$sample_id ] # 7562 features, 701 samples
   # u <- match(sampleNames(x.eset),pdata_panCanGex$sample_id)
   #
   # if(length(u[!is.na(u)])!=length(colnames(x.eset))) stop("not all samples (colnames) are in pdata_panCanGex")
   # pdata <- pdata_panCanGex[u,]
   # pdata <- pdata %>% mutate(taxonomy_simplified = factor(taxonomy_simplified))
   # stopifnot(identical(pdata$sample_id, sampleNames(x.eset)))
   #
   # ## SAMPLE INFO and ORDER
   # table(pdata$taxonomy_simplified)
   #
   # my_groups <- c("crtx","crtmed","med","infl1","infl2","ccpRCC","Ch-e","chONC", "pONC","CC-e.2","CC-e.1","CC-e.3",
   #    "P-e.1a","P-e.1b","P-e.2","P.CIMP-e","CIMP_o","rcc_UroMesInf",
   #    "outlier",
   #    "uro","Urobasal","GU","Basal/SCCL","Mes-Inf","Inf_other","SC/NE","blca_normal")
   #
   # pdata <- pdata %>%
   #    mutate(taxonomy_dl = factor(taxonomy_dl, levels=my_groups)) %>%
   #    arrange(taxonomy_dl) %>%
   #    mutate(taxonomy_simplified = factor(taxonomy_simplified, levels = c("ccpRCC","chRCC","chONC","pONC","ccRCC","pRCC_e1","pRCC_e2","CIMP_p","CIMP_o","rcc_UroMesInf","uro","outlier")))
   #
   #
   # ## fdata and probe info
   # ## -------------------------
   # fdata <- fData(x.eset)
   # x.eset <- x.eset[,pdata$sample_id]
   # ## Row Center /normalize eset
   # my_mat <- exprs(x.eset)
   # my_mat <- t(apply(my_mat, 1, function(x) x-mean(x)))

   data("pdata_panCanGex", package="dlfoo2")
   pdata <- pdata_panCanGex

   message("... Retrieving curated sample data table for urinary tract TCGA samples 'pdata_panCanGex'")
   u <- match(colnames(x.eset), pdata$sample_id)
   if(length(u[!is.na(u)])!=length(colnames(x.eset))){
      if(all(is.na(u))) stop("no samples (colnames) are found in pdata_panCanGex")
      message("... ... removing samples not present in pdata file")
      message("... ... ... from: ", length(u))
      message("... ... ... to: ", length(which(!is.na(u))))
      x.eset <- x.eset[,!is.na(u)]
      u<-u[!is.na(u)]
      }
   pdata <- pdata[u,]
   pdata <- pdata %>% arrange(taxonomy_simplified, taxonomy_dl)
   x.eset <- x.eset[,pdata$sample_id]
   stopifnot(identical(pdata$sample_id, colnames(x.eset)))
   my_mat <- exprs(x.eset)
   my_mat <- t(apply(my_mat, 1, function(x) x-mean(x)))
   fdata <- fData(x.eset)




   for(i in 1:length(amapDist)){
      ## Create subfolder for Distanace method
      myDist <- amapDist[i]
      distPlotDir <- file.path(mainPlotDir, myDist)
      if(!dir.exists(distPlotDir)) dir.create(distPlotDir, recursive = T)
      message(" <- Using amap::Dist measure ", myDist, " ->")

      ## HCA - Features
      ## ---------------------------------------
         message("... HCA of features/rows")
         my_kf <- clusterFeatures_k
         my_colors <- rep(palette.dl.annot[[19]],5)[1:my_kf]
         mydist=function(c) {amap::Dist(c,method = myDist)}
         myclust=function(c) {hclust(c,method="ward.D2")}
         hca_rows <- my_mat %>% mydist %>% myclust %>% as.dendrogram

         message("... ... Plotting HCA tree for rows to file ")
         myFile <- paste0(distPlotDir, "/", runName,"_hcaTree.pdf")
         pdf(myFile, width=11.692, height=8.267)
            hca_rows %>%
               dendextend::set("labels_cex", value=0.6)  %>%
                  #set("labels", value=my.labels) %>%
                  dendextend::set("labels", value=NA) %>%
                  dendextend::set("branches_k_color", value=my_colors, k=my_kf)  %>%
                  dendextend::set("branches_lwd", value = 1) %>%
                  plot(yaxt="n", horiz=T)
            dev.off()

         my_k_classes <- dendextend::cutree(hca_rows, k = my_kf, order_clusters_as_data = FALSE)
         str(my_k_classes)
         #dend.list.split[[my.group]] <- cbind(fdata[names(my.k.classes),], my.k.classes)
         my_tab <- dplyr::left_join(data.frame(ENSG=names(my_k_classes), row_cl=paste0("cl_",my_k_classes), stringsAsFactors = F),
            as.data.frame(fdata), by="ENSG")
         str(my_tab)
         my_tab$row_cl_col <- as.factor(my_tab$row_cl)
         levels(my_tab$row_cl_col) <- my_colors
         my_tab$row_cl_col  <- as.character(my_tab$row_cl_col)
         head(my_tab)
         myFile <- paste0(distPlotDir, "/", runName,"_hcaRows.txt")
         write.table(my_tab, file=myFile, sep="\t", row.names=F, quote=F)
   #
      # # Cluster Columns (Samples - ALL samples)
      # ## ----------------
         message("... HCA of samples/columns")
         my_ks <- clusterSamples_k
         my_colors <- rep(palette.dl.annot[[19]],5)[1:my_ks]
         #hca_rows <- hca(t(my_mat), plot=F)$hclust
         mydist=function(c) {amap::Dist(c,method=myDist)}
         myclust=function(c) {hclust(c,method="ward.D2")}
         hca_cols <- t(my_mat) %>% mydist %>% myclust %>% as.dendrogram
         my_k_classes_cols <- dendextend::cutree(hca_cols, k = my_ks, order_clusters_as_data = FALSE)

      ## Alluvial graph
      ## --------------
         message("... ... Plotting Alluvial for HCA sample tree")
         hca_k_order <- c(paste0("cl_",1:my_ks))
         my_t <- pdata %>%
            dplyr::select(sample_id, taxonomy_simplified) %>%
            mutate(HCA_k=paste0("cl_",my_k_classes_cols[pdata$sample_id])) %>%
            mutate(HCA_k = factor(HCA_k, levels=hca_k_order))

         #ggalluvial::is_alluvia_form(my_t, axes = 1:6, silent = TRUE)
         #my_t_melt <- reshape2::melt(my_t, id.vars=c("sample_id","taxonomy"), measure.vars="HCA_k")
         my_a <- my_t %>% group_by(taxonomy_simplified, HCA_k) %>%
            summarise(n = n())
         myFile <- paste0(distPlotDir, "/", runName,"_TaxAlluvial.pdf")
         pdf(myFile)
            alluvial(my_a[,c("taxonomy_simplified","HCA_k")], freq=my_a$n,
               cex=0.5,
              #gap.width = 0.25,
               col= color_subtypes[match(my_a$taxonomy_simplified, names(color_subtypes))]
               #ordering = list(
               #   match(my_a$taxonomy_simplified, my_a$taxonomy_simplified),
               #   match(my_a$HCA_k, hca_k_order))
            )
            dev.off()

         myFile <- paste0(distPlotDir, "/", runName,"_hcaColumns.txt")
         write.table(my_t, file=myFile, sep="\t", row.names=F, quote=F)



      ##    Row annotations (rcol, heatmap.3)
      ## ------------------------------
      rownames(my_tab) <- my_tab$ENSG
      row_annot <- data.frame(row_cl=my_tab$row_cl_col)
      str(row_annot)
      rownames(row_annot) <- rownames(my_tab)
      rcol <- as.matrix(t(row_annot))[,rownames(my_mat)]
      if(class(rcol)=="character"){
         rcol=matrix(data=rcol, nrow=1, dimnames = list(c(),names(rcol)))
         }
      stopifnot(identical(rownames(my_mat), colnames(rcol)))


      ##    column annotations (ccol, heatmap.3)
      ## ------------------------------
      data("color_subtypes", package="dlfoo2")
      identical(colnames(my_mat), pdata$sample_id)
      annot <- as.factor(pdata$taxonomy_simplified)
      levels(annot) <- color_subtypes[match(levels(annot), names(color_subtypes))]

      annot1 <- as.factor(pdata$taxonomy_dl)
      levels(annot1) <- color_subtypes[match(levels(annot1), names(color_subtypes))]

      annot2 <- as.factor(my_k_classes_cols[colnames(my_mat)])
      levels(annot2) <- my_colors[1:length(levels(annot2))]

      ccol <- cbind(as.character(annot), as.character(annot1), as.character(annot2))
      colnames(ccol) <- c("taxonomy_simp","taxonomy","cluster")
      rownames(ccol) <- colnames(my_mat)
      head(ccol)


      ## Plot PDFs - heatmap.3
      ## ---------------------------
      #my_heatpal <- palette_dl_heatmap2(my_mat, pal.vec= colorRamps::blue2yellow(21), palette.sat = c(0, 1))
      my_heatpal <- dlfoo2::palette_dl_heatmap2(my_mat, pal.no = palNo, palette.sat = palSat)

      myPdf <- paste(runName,"_heatmap2way_",paste(dim(my_mat),collapse = "x"),".pdf", sep = "")
      myFile <- file.path(distPlotDir, myPdf)
     # my_file <- "HCA_Me_change/HCA_Me_change_NvsT/HCA_Me_N_vs_T_heatmap_PearsonWards_9984x281.pdf"
      pdf(file=myFile, height = 9, width = 7)
         dlfoo2::heatmap.3(my_mat,
            col=my_heatpal$legend$col,
            breaks=my_heatpal$legend$legend,
            trace=c("none"),
            notecol = NA,
            Rowv = as.dendrogram(hca_rows %>% set("branches_k_color", value=my_colors, k=my_kf)),
            Colv= as.dendrogram(hca_cols %>% set("branches_k_color", value=my_colors[1:my_ks], k=my_ks)),
            dendrogram = "both",
            labRow = NA,
            labCol = NA,

            RowSideColors = rcol,
            ColSideColors = ccol,
            RowSideColorsSize=5
            #cex=0.7
            )
         dev.off()
      myPdf <- paste(runName,"_heatmap2way_",paste(dim(my_mat),collapse = "x"),".png", sep = "")
      myFile <- file.path(distPlotDir, myPdf)
     # my_file <- "HCA_Me_change/HCA_Me_change_NvsT/HCA_Me_N_vs_T_heatmap_PearsonWards_9984x281.pdf"
      png(file=myFile, width = 596, height =  842)
         dlfoo2::heatmap.3(my_mat,
            col=my_heatpal$legend$col,
            breaks=my_heatpal$legend$legend,
            trace=c("none"),
            notecol = NA,
            Rowv = as.dendrogram(hca_rows %>% set("branches_k_color", value=my_colors, k=my_kf)),
            Colv= as.dendrogram(hca_cols %>% set("branches_k_color", value=my_colors[1:my_ks], k=my_ks)),
            dendrogram = "both",
            labRow = NA,
            labCol = NA,

            RowSideColors = rcol,
            ColSideColors = ccol,
            RowSideColorsSize=5
            #cex=0.7
            )
         dev.off()


      ## 2-way
      ## ------
      myPdf <- paste(runName,"_heatmap1way_",paste(dim(my_mat),collapse = "x"),".pdf", sep = "")
      myFile <- file.path(distPlotDir, myPdf)
         pdf(file=myFile, height = 9, width = 7)
         dlfoo2::heatmap.3(my_mat,
            col=my_heatpal$legend$col,
            breaks=my_heatpal$legend$legend,
            trace=c("none"),
            notecol = NA,
            Rowv = as.dendrogram(hca_rows %>% set("branches_k_color", value=my_colors, k=my_kf)),
            Colv= NA,
            dendrogram = "row",
            labRow = NA,
            labCol = NA,

            RowSideColors = rcol,
            ColSideColors = ccol,
            RowSideColorsSize=5
            #cex=0.7
            )
         dev.off()
      myPdf <- paste(runName,"_heatmap1way_",paste(dim(my_mat),collapse = "x"),".png", sep = "")
      myFile <- file.path(distPlotDir, myPdf)
         png(file=myFile,  width = 596, height =  842)
         dlfoo2::heatmap.3(my_mat,
            col=my_heatpal$legend$col,
            breaks=my_heatpal$legend$legend,
            trace=c("none"),
            notecol = NA,
            Rowv = as.dendrogram(hca_rows %>% set("branches_k_color", value=my_colors, k=my_kf)),
            Colv= NA,
            dendrogram = "row",
            labRow = NA,
            labCol = NA,

            RowSideColors = rcol,
            ColSideColors = ccol,
            RowSideColorsSize=5
            #cex=0.7

            )
         dev.off()
      } #end amap dist i
   }







#' function Plot GEX eset heatmaps to file - for signatures
#' @description Version 1.0 20190129
#' @family plot wraps
#' @family heatmap
#' @family eset
#' @family gene expression
#' @family signature
#' @param x.eset filtered gene expression eSet with expression values to be plotted. Input shpuld be Sorted and Filtered!
#' @param signature.object object aproduct from featureGetBM analysis. A 'signature' slot should have been added.
#' @param signature.annotation.column defaults to 'cell.type.tax'. The annotation to break down the signature. Can be NULL if no subdivistion of the signature,
#' @param signature.id.type id type to match signature - gex eset. degfaults to 'hgnc_symbol'
#' @param eset.id.type SYMBOL default
#' @param clusterRows if to cluster signature rows (and each of the sub-signatuers)
#' @param runName name of plot/analyses - subfolder
#' @param runFolder directory where to plot/save objects
#' @param amapDist Supplied to amap::Dist function. Can be ONE of following distance methods. "euclidean", "maximum", "manhattan", "canberra", "binary", "pearson", "abspearson", "correlation", "abscorrelation", "spearman" or "kendall". Any unambiguous substring can be given.
#' @param palNo What palette (dl_heat)
#' @param palSat saturation levels for palette. detfaults to '-3,3'
#' @return heatmaps to file
#' @export
#'
esetGexSignatures <- function(
   x.eset,
   signature.object, signature.annotation.column='cell.type.tax',
   signature.id.type='hgnc_symbol', eset.id.type = "SYMBOL",
   clusterRows =TRUE,
   sd.cut = 0.25,
   runName, runFolder,
   amapDist = c("pearson"),
   palNo=3, palSat=c(-3, 3)
   ){


   require(dplyr)
   options(dplyr.width = Inf)
   require(tibble)
   require(tidyverse)
   require(tidyr)
   require(ggplot2)

   require(alluvial)
   require(ggalluvial)

   if(class(x.eset)!="ExpressionSet") stop("gset is not of class 'ExpressionSet'")
   if(!signature.annotation.column %in% colnames(signature.object)) stop("'signature.annotation.column' not present in 'signature.object'")
   if(!signature.id.type %in% colnames(signature.object)) stop("'signature.id.type' not present in 'signature.object'")
   if(!eset.id.type %in% colnames(fData(x.eset))) stop("'eset.id.type' not in fdata for eset")

   amapDist <- match.arg(amapDist, choices = c("euclidean", "maximum", "manhattan", "canberra", "binary", "pearson", "abspearson", "correlation", "abscorrelation", "spearman", "kendall"),
         several.ok = F)
   signature.object <- signature.object %>% arrange(cell_)

   #gsSinkGenosetDescription()
   mainPlotDir <- file.path(runFolder, runName)
   if(!file.exists(mainPlotDir)) dir.create(mainPlotDir, recursive = T)
   message("   ---> plotting GEX signature(s) in folder: ", mainPlotDir)


   data("pdata_panCanGex", package="dlfoo2")
   pdata <- pdata_panCanGex


   message("... Retrieving curated sample data table for urinary tract TCGA samples 'pdata_panCanGex'")
   u <- match(colnames(x.eset), pdata$sample_id)
   if(length(u[!is.na(u)])!=length(colnames(x.eset))){
      if(all(is.na(u))) stop("no samples (colnames) are found in pdata_panCanGex")
      message("... ... removing samples not present in pdata file")
      message("... ... ... from: ", length(u))
      message("... ... ... to: ", length(which(!is.na(u))))
      x.eset <- x.eset[,!is.na(u)]
      u<-u[!is.na(u)]
      }
   pdata <- pdata[u,]
   pdata <- pdata %>% arrange(taxonomy_simplified, taxonomy_dl)
   x.eset <- x.eset[,pdata$sample_id]
   stopifnot(identical(pdata$sample_id, colnames(x.eset)))
   my_mat <- exprs(x.eset)
   message("normalizing matrix & applying sd filter")
   x.eset <- dlfoo2::esetVarianceFilter(x.eset, cutoff = sd.cut)
   my_mat <- t(apply(my_mat, 1, function(x) {
      y <- x-mean(x)
      #y/sd(y)
      }
      ))
   fdata <- fData(x.eset)





   ## Loop all levels of the signature grouping variable
   ## ----------------------------
   sig_annot <- as.factor(data.frame(signature.object)[,signature.annotation.column])
   sig_groups <- levels(sig_annot)

   mat_list <- list()
   fdata_list <- list()
   for(i in 1:length(sig_groups)){
      group_df <- data.frame(signature.object)
      group_ids <- group_df[which(group_df[,signature.annotation.column]==sig_groups[i]), signature.id.type]
      group_ids <- group_ids[!is.na(group_ids)]
      group_ids <- group_ids[!group_ids==""]

      group_ids <- group_ids[!duplicated(group_ids)] # remove duplicates
      group_ids <- intersect(group_ids, fdata[,eset.id.type])
      u <- match(group_ids, fdata[,eset.id.type])
      if(length(u)==0) next(i)
      if(length(u)>1) group_mat <- my_mat[u,]
      if(length(u)==1) {
         group_mat <- matrix(nrow=1, data=my_mat[u,])
         colnames(group_mat) <- colnames(my_mat)
      }
      rownames(group_mat) <- fdata[u,eset.id.type]

      if(clusterRows & nrow(group_mat)>2){
         myDist <- amapDist
         #distPlotDir <- file.path(mainPlotDir, myDist)
         #if(!dir.exists(distPlotDir)) dir.create(distPlotDir, recursive = T)
         #message(" <- Using amap::Dist measure ", myDist, " ->")
         mydist=function(c) {amap::Dist(c,method = myDist)}
         myclust=function(c) {stats::hclust(c,method="ward.D2")}
         hca_rows <- group_mat %>% mydist %>% myclust %>% as.dendrogram
         group_mat <- group_mat[labels(hca_rows), ]
         }
      group_fdata <- group_df[match(rownames(group_mat), group_df[,signature.id.type]),]

      rownames(group_fdata) <- apply(group_fdata[,c(signature.annotation.column, signature.id.type)],1, function(x) paste0(x, collapse =  "_"))
      rownames(group_mat) <- rownames(group_fdata)
      mat_list[[i]] <- group_mat
      fdata_list[[i]] <- group_fdata
      } # end i


   plot_mat <- do.call("rbind", mat_list)
   plot_fdata <- do.call("rbind",fdata_list)



   ##    Row annotations (rcol, heatmap.3)
   ## ------------------------------
   data("color_celltypes", package="dlfoo2")

   annot <- plot_fdata[,signature.annotation.column]
   levels(annot) <- color_celltypes$color[match(levels(annot), color_celltypes$cell.type.tax)]
   rcol=matrix(data=annot, nrow=1, dimnames = list(c(),rownames(plot_fdata)))
   stopifnot(identical(rownames(plot_fdata), colnames(rcol)))


   ##    column annotations (ccol, heatmap.3)
   ## ------------------------------
   data("color_subtypes", package="dlfoo2")
   stopifnot(identical(colnames(my_mat), pdata$sample_id))
   annot <- as.factor(pdata$taxonomy_simplified)
   levels(annot) <- color_subtypes[match(levels(annot), names(color_subtypes))]

   annot1 <- as.factor(pdata$taxonomy_dl)
   levels(annot1) <- color_subtypes[match(levels(annot1), names(color_subtypes))]

   ccol <- cbind(as.character(annot), as.character(annot1))
   colnames(ccol) <- c("taxonomy_simp","taxonomy")
   rownames(ccol) <- colnames(my_mat)
   head(ccol)


   ## Plot PDFs - heatmap.3
   ## ---------------------------
   #my_heatpal <- palette_dl_heatmap2(my_mat, pal.vec= colorRamps::blue2yellow(21), palette.sat = c(0, 1))
   my_heatpal <- dlfoo2::palette_dl_heatmap2(my_mat, pal.no = palNo, palette.sat = palSat)

   ## Plot Heatmap
      dlfoo2::heatmap.3(
            plot_mat,
            col=my_heatpal$legend$col,
            breaks=my_heatpal$legend$legend,
            trace=c("none"),
            notecol = NA,
            Rowv = NA,
            Colv= NA,
            dendrogram = "none",
            labRow = NA,
            labCol = NA,

            RowSideColors = rcol,
            ColSideColors = ccol,
            RowSideColorsSize=5
            #cex=0.7
            )

   } # end -> esetGexSignatures



#
#    fdata <- readRDS("~/PROJECTS/BRCA_SR_KP/Data/fdata.rds")
#    pdata <- readRDS( "~/PROJECTS/BRCA_SR_KP/Data/pdata.rds")
#    pdata <- pdata[,c("sample", "seqno","time","cells","annot")]
#
#    sig_list <- readRDS("~/PROJECTS/BRCA_SR_KP/DESeq2_analysis/Significant_genes.rds")
#    vst <- readRDS(file="~/PROJECTS/BRCA_SR_KP/DESeq2_analysis/DeSeq2_vst.rds")
#
#    my_mat <- assay(vst)
#    my_mat <- my_mat[,pdata$sample]
#    str(fdata)
#    identical(fdata$gene_id, rownames(my_mat))
#    rownames(fdata) <- fdata$gene_id
#
#    names(sig_list)
#    my_run <- "MCF7_E2_24h_cells_CAF_vs_noCAF"
#    u <- match(sig_list[[my_run]]$gene_id, rownames(my_mat))
#    x <- my_mat[u, ]
#    fdata <- fdata[u,]
#    fdata$regulation <- if_else(sig_list[[my_run]]$stat<0, "neg","pos")
#    str(x)
#    str(fdata)
#    identical(rownames(x), rownames(fdata))
#
#    runName <- "heatmap_temp"
#    runFolder <-     "~/PROJECTS/BRCA_SR_KP/DESeq2_analysis"
#    pdata_color <- read.delim("~/PROJECTS/BRCA_SR_KP/Data/pdata_color_key.txt", as.is = T)
#    color.key.pdata <- pdata_color$color
#    names(color.key.pdata) <- pdata_color$name
#    color.key.genes = color.key.pdata
# pdata.id.column = "sample"
# fdata.id.column = "gene_id"


#' plot heatmap - genereic function
#' @description Version 1.0 20191018
#' @family plot wraps
#' @family heatmap
#' @family gene expression
#' @param x expression matrix. Input shpuld be Sorted and Filtered!
#' @param center.rows if to center rows
#' @param pdata data frame. sample colData with sample phenotype labels. Should contain 1 column with sample_id that match colnames for pdata/coldata. remaining columns will be plotted. The order of pdata will determine the sample-order to plot
#' @param pdata.id.column name of primary sample id column in pdata file. unique names used as colnames in matrix
#' @param color.key.pdata color key for annottaion of samples
#' @param fdata feature info if to plot gene-groups
#' @param fdata.annot.columns what fdata annotations to plot
#' @param color.key.genes color key for annotation of genes in signature
#' @param runName name of plot/analyses
#' @param runFolder name of directory
#' @param amapDist Supplied to amap::Dist function. Can be ONE of following distance methods. "euclidean", "maximum", "manhattan", "canberra", "binary", "pearson", "abspearson", "correlation", "abscorrelation", "spearman" or "kendall". Any unambiguous substring can be given.
#' @param pal.heatmap What palette
#' @param pal.sat saturation levels for palette. detfaults to '-3,3'
#' @return heatmaps to file
#' @export
#'
plotHeatmap <- function(
   x, center.rows=T, cluster.rows =T,
   pdata=NULL, pdata.id.column="sample_id",
   fdata=NULL, fdata.annot.columns = NULL,
   color.key.pdata,
   color.key.genes,
   pal.heatmap= rev(c(rev(RColorBrewer::brewer.pal(9,"YlOrRd")[1:9]),rep("#FFFFCC",1),rep("white",1), rep("#FFFFCC",1), RColorBrewer::brewer.pal(9,"YlGn")[1:9]) ),
   pal.sat=c(-3, 3),
   # amapDist = c("pearson"),
   trace="none",
   notecol = NA,
   dendrogram = "none",
   Rowv = NA, Colv= NA,
   labRow = NA, labCol = NA, ...
   ){

   require(dplyr)
   options(dplyr.width = Inf)
   require(tibble)
   require(tidyverse)
   require(tidyr)
   require(ggplot2)



   # amapDist <- match.arg(amapDist, choices = c("euclidean", "maximum", "manhattan", "canberra", "binary", "pearson", "abspearson", "correlation", "abscorrelation", "spearman", "kendall"),
   #      several.ok = F)

   # mainPlotDir <- file.path(runFolder, runName)
   # if(!file.exists(mainPlotDir)) dir.create(mainPlotDir, recursive = T)
   # message("   ---> plotting GEX signature(s) in folder: ", mainPlotDir)


   ## Pdata stuff
   if(!is.null(pdata)){
      pdata <- data.frame(pdata)
      message("... matching matrix with supplided pdata")
      if(!pdata.id.column %in% colnames(pdata)) stop("'sample_id' not in pdata colnames. sample_ids must be same as colnames in matrix")
      pdata[,pdata.id.column] <- as.character(pdata[,pdata.id.column]) # if factor
      u <- match(colnames(x), pdata[,pdata.id.column])

      if(length(u[!is.na(u)])!=length(colnames(x))){
         if(all(is.na(u))) stop("no samples (colnames) are found in pdata")
         message("... ... removing samples not present in pdata file")
         message("... ... ... from: ", length(u))
         message("... ... ... to: ", length(which(!is.na(u))))
         x <- x[,!is.na(u)]
         u<-u[!is.na(u)]
         }
      x <- x[,pdata[,pdata.id.column]]
      stopifnot(identical(pdata[,pdata.id.column], colnames(x)))
   }
   # str(x)


   # fdata stuff
   if(!is.null(fdata)){
      message("... matching matrix with supplided fdata")
      stopifnot(identical(rownames(fdata), rownames(x)))
      #if(!fdata.id.column %in% colnames(fdata)) stop("'sample_id' not in pdata colnames. sample_ids must be same as colnames in matrix")
      #fdata[,fdata.id.column] <- as.character(fdata[,fdata.id.column]) # if factor
   }
   if(is.null(fdata)){
      fdata <- data.frame(gene_id=rownames(x))
   }
   #str(fdata)


   ## generate matrix
   message("... generate row-centered matrix for clustering and plot")
   if(center.rows){
      x <- t(apply(x, 1, function(y) {
         y <- y-mean(y)
      }
      ))
   }


   # ## cluster rows (non row-centered data)
   # if(cluster.rows){
   #    message("...clustering rows")
   #    myDist <- amapDist
   #    mydist=function(c) {amap::Dist(c,method = myDist)}
   #    myclust=function(c) {stats::hclust(c,method="ward.D2")}
   #    hca_rows <- x %>% mydist %>% myclust %>% as.dendrogram
   #    x <- x[labels(hca_rows), ]
   #    fdata <- fdata[labels(hca_rows), ]
   # }

   ##Fix fdata order and rownames
   # fdata <- fdata[match(rownames(x), fdata[,fdata.id.column]), ]
   # if(!is.null(fdata.rownames.column)){
   #    if(!c(fdata.rownames.column %in% colnames(fdata))) stop("fdata.rownames.column not in fdata")
   #    y <- fdata[,fdata.rownames.column]
   #    if(any(y=="")|any(is.na(y))|any(duplicated(y))) y <- paste(y, 1:length(y), sep = "__")
   #    rownames(fdata) <- y
   #    rownames(x) <- y
   # }



   ##    Row annotations - rowside (rcol, heatmap.3)
   ## ------------------------------
   if(!is.null(fdata.annot.columns)){
      fdata_annot_list <- as.list(fdata)
      u<-match(fdata.annot.columns, names(fdata_annot_list))
      fdata_annot_list <- fdata_annot_list[u[!is.na(u)]]
      if(length(fdata_annot_list)>0){
         fdata_annot_list <- lapply(fdata_annot_list, function(y){
            y <-factor(y)
            if(is.null(color.key.genes)) levels(y) <- dlfoo2::palette.dl.annot[[19]][1:length(levels(y))]
            if(!is.null(color.key.genes)) y<-recode_factor(y, !!!color.key.genes) # levels(y) <- color.key.genes[match(levels(y), names(color.key.genes))]
         })
      rcol=matrix(data=unlist(fdata_annot_list), nrow=length(fdata_annot_list), dimnames = list(names(fdata_annot_list), rownames(x) ))
      # stopifnot(identical(rownames(plot_fdata), colnames(rcol)))
      }
   }


   ##    Row annotations - colside (ccol, heatmap.3)
   ## ------------------------------
   stopifnot(identical(colnames(x), pdata[,pdata.id.column]))
   pdata_annot_list <- as.list(pdata)
   u<-match(c(pdata.id.column), names(pdata_annot_list))
   pdata_annot_list <- pdata_annot_list[-u]
   if(length(pdata_annot_list)>0){
      pdata_annot_list <- lapply(pdata_annot_list, function(y){
         y <- factor(y)
         if(is.null(color.key.pdata)) levels(y) <- dlfoo2::palette.dl.annot[[19]][1:length(levels(y))]
         if(!is.null(color.key.pdata)) y <-recode_factor(y, !!!color.key.pdata) #color.key.pdata[match(levels(y), names(color.key.pdata))]
         return(y)
      })
      ccol=matrix(data=unlist(pdata_annot_list), ncol=length(pdata_annot_list), dimnames = list(colnames(x), names(pdata_annot_list)))
   }



   ## Plot PDFs - heatmap.3
   ## ---------------------------
   #my_heatpal <- palette_dl_heatmap2(my_mat, pal.vec= colorRamps::blue2yellow(21), palette.sat = c(0, 1))
   my_heatpal <- dlfoo2::palette_dl_heatmap2(x, pal.vec = pal.heatmap ,  palette.sat = pal.sat)


   ## Plot Heatmap
   # my_pars <- list(
   #          x = x,
   #          col=my_heatpal$legend$col,
   #          breaks=my_heatpal$legend$legend,
   #          trace=c("none"),
   #          notecol = NA,
   #          Rowv = NA,
   #          Colv= NA,
   #          dendrogram = "none",
   #          labRow = NA,
   #          labCol = NA,
   #          )
   #if(exists("ccol")) my_pars$ColSideColors = ccol
   #if(exists("rcol")) my_pars$RowSideColors = rcol
   #do.call(dlfoo2::heatmap.3, args=my_pars)

   if(!exists("ccol")) ccol <- NULL
   if(!exists("rcol")) rcol <- NULL

   heatmap.3(x = x, col= my_heatpal$legend$col, breaks=my_heatpal$legend$legend,
            ColSideColors=ccol, RowSideColors=rcol,
            trace=trace,
            Rowv = Rowv,
            notecol = notecol,
            Colv= Colv,
            dendrogram = dendrogram,
            labRow = labRow,
            labCol = labCol, ...
            )

   } # plot egeneric heatmap






#' Generate a hetmap.
#' @description Uses code from devtools::source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
#' @family plot functions
#' @param x a matrix containg expression values
#' @param ... additional stuff
#' @return plot to quartz
#' @export

heatmap.3 <- function(x,
                      Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
                      distfun = dist,
                      hclustfun = hclust,
                      dendrogram = c("both","row", "column", "none"),
                      symm = FALSE,
                      scale = c("none","row", "column"),
                      na.rm = TRUE,
                      revC = identical(Colv,"Rowv"),
                      add.expr,
                      breaks,
                      symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",
                      col = "heat.colors",
                      colsep,
                      rowsep,
                      sepcolor = "white",
                      sepwidth = c(0.05, 0.05),
                      cellnote,
                      notecex = 1,
                      notecol = "cyan",
                      na.color = par("bg"),
                      trace = c("none", "column","row", "both"),
                      tracecol = "cyan",
                      hline = median(breaks),
                      vline = median(breaks),
                      linecol = tracecol,
                      margins = c(5,5),
                      ColSideColors = NULL,
                      RowSideColors = NULL,
                      side.height.fraction=0.3,
                      cexRow = 0.2 + 1/log10(nr),
                      cexCol = 0.2 + 1/log10(nc),
                      labRow = NULL,
                      labCol = NULL,
                      key = TRUE,
                      keysize = 1.5,
                      density.info = c("none", "histogram", "density"),
                      denscol = tracecol,
                      symkey = max(x < 0, na.rm = TRUE) || symbreaks,
                      densadj = 0.25,
                      main = NULL,
                      xlab = NULL,
                      ylab = NULL,
                      lmat = NULL,
                      lhei = NULL,
                      lwid = NULL,
                      ColSideColorsSize = 1,
                      RowSideColorsSize = 1,
                      KeyValueName="Value",...){

    invalid <- function (x) {
      if (missing(x) || is.null(x) || length(x) == 0)
          return(TRUE)
      if (is.list(x))
          return(all(sapply(x, invalid)))
      else if (is.vector(x))
          return(all(is.na(x)))
      else return(FALSE)
    }

    x <- as.matrix(x)
    scale01 <- function(x, low = min(x), high = max(x)) {
        x <- (x - low)/(high - low)
        x
    }
    retval <- list()
    scale <- if (symm && missing(scale))
        "none"
    else match.arg(scale)
    dendrogram <- match.arg(dendrogram)
    trace <- match.arg(trace)
    density.info <- match.arg(density.info)
    if (length(col) == 1 && is.character(col))
        col <- get(col, mode = "function")
    if (!missing(breaks) && (scale != "none"))
        warning("Using scale=\"row\" or scale=\"column\" when breaks are",
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
    if (is.null(Rowv) || is.na(Rowv))
        Rowv <- FALSE
    if (is.null(Colv) || is.na(Colv))
        Colv <- FALSE
    else if (Colv == "Rowv" && !isTRUE(Rowv))
        Colv <- FALSE
    if (length(di <- dim(x)) != 2 || !is.numeric(x))
        stop("`x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1)
        stop("`x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2)
        stop("`margins' must be a numeric vector of length 2")
    if (missing(cellnote))
        cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
    if (!inherits(Rowv, "dendrogram")) {
        if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
            c("both", "row"))) {
            if (is.logical(Colv) && (Colv))
                dendrogram <- "column"
            else dedrogram <- "none"
            warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
                dendrogram, "'. Omitting row dendogram.")
        }
    }
    if (!inherits(Colv, "dendrogram")) {
        if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
            c("both", "column"))) {
            if (is.logical(Rowv) && (Rowv))
                dendrogram <- "row"
            else dendrogram <- "none"
            warning("Discrepancy: Colv is FALSE, while dendrogram is `",
                dendrogram, "'. Omitting column dendogram.")
        }
    }
    if (inherits(Rowv, "dendrogram")) {
        ddr <- Rowv
        rowInd <- order.dendrogram(ddr)
    }
    else if (is.integer(Rowv)) {
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd))
            stop("row dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Rowv)) {
        Rowv <- rowMeans(x, na.rm = na.rm)
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd))
            stop("row dendrogram ordering gave index of wrong length")
    }
    else {
        rowInd <- nr:1
    }
    if (inherits(Colv, "dendrogram")) {
        ddc <- Colv
        colInd <- order.dendrogram(ddc)
    }
    else if (identical(Colv, "Rowv")) {
        if (nr != nc)
            stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
        if (exists("ddr")) {
            ddc <- ddr
            colInd <- order.dendrogram(ddc)
        }
        else colInd <- rowInd
    }
    else if (is.integer(Colv)) {
        hcc <- hclustfun(distfun(if (symm)
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd))
            stop("column dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Colv)) {
        Colv <- colMeans(x, na.rm = na.rm)
        hcc <- hclustfun(distfun(if (symm)
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd))
            stop("column dendrogram ordering gave index of wrong length")
    }
    else {
        colInd <- 1:nc
    }
    retval$rowInd <- rowInd
    retval$colInd <- colInd
    retval$call <- match.call()
    x <- x[rowInd, colInd]
    x.unscaled <- x
    cellnote <- cellnote[rowInd, colInd]
    if (is.null(labRow))
        labRow <- if (is.null(rownames(x)))
            (1:nr)[rowInd]
        else rownames(x)
    else labRow <- labRow[rowInd]
    if (is.null(labCol))
        labCol <- if (is.null(colnames(x)))
            (1:nc)[colInd]
        else colnames(x)
    else labCol <- labCol[colInd]
    if (scale == "row") {
        retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
        x <- sweep(x, 1, rm)
        retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
        retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
        x <- sweep(x, 2, rm)
        retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/")
    }
    if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
        if (missing(col) || is.function(col))
            breaks <- 16
        else breaks <- length(col) + 1
    }
    if (length(breaks) == 1) {
        if (!symbreaks)
            breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                length = breaks)
        else {
            extreme <- max(abs(x), na.rm = TRUE)
            breaks <- seq(-extreme, extreme, length = breaks)
        }
    }
    nbr <- length(breaks)
    ncol <- length(breaks) - 1
    if (class(col) == "function")
        col <- col(ncol)
    min.breaks <- min(breaks)
    max.breaks <- max(breaks)
    x[x < min.breaks] <- min.breaks
    x[x > max.breaks] <- max.breaks
    if (missing(lhei) || is.null(lhei))
        lhei <- c(keysize, 4)
    if (missing(lwid) || is.null(lwid))
        lwid <- c(keysize, 4)
    if (missing(lmat) || is.null(lmat)) {
        lmat <- rbind(4:3, 2:1)

        #if (!missing(ColSideColors)) {
        if (!is.null(ColSideColors)) {
           #if (!is.matrix(ColSideColors))
           #stop("'ColSideColors' must be a matrix")
            if (!is.character(ColSideColors) || nrow(ColSideColors) != nc)
                stop("'ColSideColors' must be a matrix of nrow(x) rows")
            lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
            #lhei <- c(lhei[1], 0.2, lhei[2])
             lhei=c(lhei[1], side.height.fraction*ColSideColorsSize/2, lhei[2])
        }

        #if (!missing(RowSideColors)) {
        if (!is.null(RowSideColors)) {
            #if (!is.matrix(RowSideColors))
            #stop("'RowSideColors' must be a matrix")
            if (!is.character(RowSideColors) || ncol(RowSideColors) != nr)
                stop("'RowSideColors' must be a matrix of ncol(x) columns")
            lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[,2] + 1)
            #lwid <- c(lwid[1], 0.2, lwid[2])
            lwid <- c(lwid[1], side.height.fraction*RowSideColorsSize/2, lwid[2])
        }
        lmat[is.na(lmat)] <- 0
    }

    if (length(lhei) != nrow(lmat))
        stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
    if (length(lwid) != ncol(lmat))
        stop("lwid must have length = ncol(lmat) =", ncol(lmat))
    op <- par(no.readonly = TRUE)
    on.exit(par(op))

    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)

    if (!is.null(RowSideColors)) {
        if (!is.matrix(RowSideColors)){
                par(mar = c(margins[1], 0, 0, 0.5))
                image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
        } else {
            par(mar = c(margins[1], 0, 0, 0.5))
            rsc = t(RowSideColors[,rowInd, drop=F])
            rsc.colors = matrix()
            rsc.names = names(table(rsc))
            rsc.i = 1
            for (rsc.name in rsc.names) {
                rsc.colors[rsc.i] = rsc.name
                rsc[rsc == rsc.name] = rsc.i
                rsc.i = rsc.i + 1
            }
            rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
            image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
            if (length(rownames(RowSideColors)) > 0) {
                axis(1, 0:(dim(rsc)[2] - 1)/max(1,(dim(rsc)[2] - 1)), rownames(RowSideColors), las = 2, tick = FALSE)
            }
        }
    }

    if (!is.null(ColSideColors)) {

        if (!is.matrix(ColSideColors)){
            par(mar = c(0.5, 0, 0, margins[2]))
            image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
        } else {
            par(mar = c(0.5, 0, 0, margins[2]))
            csc = ColSideColors[colInd, , drop=F]
            csc.colors = matrix()
            csc.names = names(table(csc))
            csc.i = 1
            for (csc.name in csc.names) {
                csc.colors[csc.i] = csc.name
                csc[csc == csc.name] = csc.i
                csc.i = csc.i + 1
            }
            csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
            image(csc, col = as.vector(csc.colors), axes = FALSE)
            if (length(colnames(ColSideColors)) > 0) {
                axis(2, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), colnames(ColSideColors), las = 2, tick = FALSE)
            }
        }
    }

    par(mar = c(margins[1], 0, 0, margins[2]))
    x <- t(x)
    cellnote <- t(cellnote)
    if (revC) {
        iy <- nr:1
        if (exists("ddr"))
            ddr <- rev(ddr)
        x <- x[, iy]
        cellnote <- cellnote[, iy]
    }
    else iy <- 1:nr
    image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
    retval$carpet <- x
    if (exists("ddr"))
        retval$rowDendrogram <- ddr
    if (exists("ddc"))
        retval$colDendrogram <- ddc
    retval$breaks <- breaks
    retval$col <- col
    if (!invalid(na.color) & any(is.na(x))) { # load library(gplots)
        mmat <- ifelse(is.na(x), 1, NA)
        image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
            col = na.color, add = TRUE)
    }
    axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
        cex.axis = cexCol)
    if (!is.null(xlab))
        mtext(xlab, side = 1, line = margins[1] - 1.25)
    axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
        cex.axis = cexRow)
    if (!is.null(ylab))
        mtext(ylab, side = 4, line = margins[2] - 1.25)
    if (!missing(add.expr))
        eval(substitute(add.expr))
    if (!missing(colsep))
        for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, length(csep)), xright = csep + 0.5 + sepwidth[1], ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    if (!missing(rowsep))
        for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    min.scale <- min(breaks)
    max.scale <- max(breaks)
    x.scaled <- scale01(t(x), min.scale, max.scale)
    if (trace %in% c("both", "column")) {
        retval$vline <- vline
        vline.vals <- scale01(vline, min.scale, max.scale)
        for (i in colInd) {
            if (!is.null(vline)) {
                abline(v = i - 0.5 + vline.vals, col = linecol,
                  lty = 2)
            }
            xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
            xv <- c(xv[1], xv)
            yv <- 1:length(xv) - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (trace %in% c("both", "row")) {
        retval$hline <- hline
        hline.vals <- scale01(hline, min.scale, max.scale)
        for (i in rowInd) {
            if (!is.null(hline)) {
                abline(h = i + hline, col = linecol, lty = 2)
            }
            yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
            yv <- rev(c(yv[1], yv))
            xv <- length(yv):1 - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (!missing(cellnote))
        text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
            col = notecol, cex = notecex)
    par(mar = c(margins[1], 0, 0, 0))
    if (dendrogram %in% c("both", "row")) {
        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    }
    else plot.new()
    par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
    if (dendrogram %in% c("both", "column")) {
        plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    }
    else plot.new()
    if (!is.null(main))
        title(main, cex.main = 1.5 * op[["cex.main"]])
    if (key) {
        par(mar = c(5, 4, 2, 1), cex = 0.75)
        tmpbreaks <- breaks
        if (symkey) {
            max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
            min.raw <- -max.raw
            tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
            tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
        }
        else {
            min.raw <- min(x, na.rm = TRUE)
            max.raw <- max(x, na.rm = TRUE)
        }

        z <- seq(min.raw, max.raw, length = length(col))
        image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
            xaxt = "n", yaxt = "n")
        par(usr = c(0, 1, 0, 1))
        lv <- pretty(breaks)
        xv <- scale01(as.numeric(lv), min.raw, max.raw)
        axis(1, at = xv, labels = lv)
        if (scale == "row")
            mtext(side = 1, "Row Z-Score", line = 2)
        else if (scale == "column")
            mtext(side = 1, "Column Z-Score", line = 2)
        else mtext(side = 1, KeyValueName, line = 2)
        if (density.info == "density") {
            dens <- density(x, adjust = densadj, na.rm = TRUE)
            omit <- dens$x < min(breaks) | dens$x > max(breaks)
            dens$x <- dens$x[-omit]
            dens$y <- dens$y[-omit]
            dens$x <- scale01(dens$x, min.raw, max.raw)
            lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
                lwd = 1)
            axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
            title("Color Key\nand Density Plot")
            par(cex = 0.5)
            mtext(side = 2, "Density", line = 2)
        }
        else if (density.info == "histogram") {
            h <- hist(x, plot = FALSE, breaks = breaks)
            hx <- scale01(breaks, min.raw, max.raw)
            hy <- c(h$counts, h$counts[length(h$counts)])
            lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
                col = denscol)
            axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
            title("Color Key\nand Histogram")
            par(cex = 0.5)
            mtext(side = 2, "Count", line = 2)
        }
        else title("Color Key")
    }
    else plot.new()
    retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
        high = retval$breaks[-1], color = retval$col)
    invisible(retval)
}



#' boxplot single gene epxression in TCGA urinary tract (UTC) data
#' @family plot wraps
#' @family tcga
#' @family gene expression
#' @family eset
#' @param ensg A ensemble gene identifier (homo sapiens)
#' @param highligt.samples NULL or character vector of sampple_id (e.g. KIRP_TP_5883). If NULL, no action. Set sample names to hihlight expression for. if >1 sample, a separate group is created in the bxp
#' @param individual.sample if to highlight one addtiontional individual sample (often within the highlight group above)
#' @param default.taxonomy what taxonomy to use
#' @return a plot
#' @export
#'
geneplotPanTCGA_violin = function(ensg=NULL, highlight.samples = NULL, individual.sample=NULL, default.taxonomy="taxonomy_published", ...){
   require(gplots)
   require(ggplot2)
   require(tidyr)
   require(tidyverse)
   require(dplyr)
   data("pan_gdc", package="dlfoo2data")
   data('pdata_panCanGex', package='dlfoo2data')
   data('color_subtypes', package='dlfoo2')

   message("boxplot gene epxression in TCGA pan_gdc data")

   if(length(ensg)!=1){
      message(ensg, " only one id is allowed")
      return(NULL)
   }

   gene.i = match(ensg, featureNames(pan_gdc))

   if(is.na(gene.i)){
      message(ensg, " not found in feature data for pan_gdc")
      return(NULL)
   }

   pdata <- pData(pan_gdc)
   if(!is.null(highlight.samples) & !all(highlight.samples %in% pdata$sample_id)) stop("ERROR : Not all highlighted samples in pdata for pdan_gdc")
   if(!is.null(individual.sample) & !all(individual.sample %in% pdata$sample_id)) stop("ERROR : Not all highlighted samples in pdata for pdan_gdc")


   my_df <- pdata[,c("sample_id","sample_type2", default.taxonomy)]
   my_df$taxonomy <- my_df[,default.taxonomy]
   my_df <- my_df %>%
      mutate(taxonomy = recode(factor(taxonomy),
         "crtx"="kidney_n", "crtmed"="kidney_n", "med"="kidney_n", "infl1"="kidney_n", "infl2"="kidney_n",
         "CC-e.2"="ccRCC_e2", "CC-e.1"="ccRCC_e1", "CC-e.3"="ccRCC_e3",
         "P-e.1a"="pRCC_e1a","P-e.1b"="pRCC_e1b", "P-e.2"="pRCC_e2",
         "P.CIMP-e"="CIMP", "Ch-e"="chRCC","mixed"="mixed",
         "Basal/SCCL" = "blca_Basal_SCCL", "GU"="blca_GU", "Inf_other"="blca_InfOther", "Mes-Inf"="blca_MesInf", "SC/NE"="blca_SC_NE", "Urobasal"="blca_Urobasal",
         "blca_normal"="bladder_n")) %>%
      mutate(taxonomy = as.character(taxonomy)) %>%
      mutate(taxonomy = if_else(is.na(taxonomy), sample_type2, taxonomy)) %>%
      mutate(taxonomy = if_else(sample_id %in% highlight.samples, "highlighted", taxonomy)) %>%
      mutate(taxonomy = factor(taxonomy, levels=
            c("highlighted","kidney_n","ccRCC_e2", "ccRCC_e1","ccRCC_e3", "pRCC_e1a", "pRCC_e1b", "pRCC_e2",
         "CIMP", "chRCC","mixed",
         "blca_Basal_SCCL", "blca_GU", "blca_InfOther", "blca_MesInf", "blca_SC_NE", "blca_Urobasal",
         "bladder_n", unique(pdata$sample_type2[-grep("BLCA|KIRC|KICH|KIRP",pdata$sample_type2)]))))
   table(my_df$taxonomy)

   my_df = my_df %>%
      mutate(gex=exprs(pan_gdc)[gene.i,])

   ## PLOT CONDENSED GROUPS
   g <- ggplot(my_df) +
      aes(x=taxonomy, y=gex, fill=taxonomy) +
      theme(panel.background = element_rect(fill = 'gray95', colour = 'antiquewhite3', size = 0.5)) +
      theme(panel.grid.major =  element_line(colour = 'lightsteelblue2'), panel.grid.minor =  element_line(colour = 'mistyrose3')) +
      theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
      labs(y=paste0("log2 FPKM"), x="") +
      guides(fill=FALSE) +
      theme(axis.title.y = element_text(size=12)) +
      theme(axis.title.x = element_text(size=12)) +
      theme(axis.text = element_text(size=9)) +
      theme(title = element_text(size=9)) +
      theme(strip.text = element_text(size=12)) +
      theme(title = element_text(size=12))
   g1 <- g +
      geom_violin(alpha=0.5, width=1, trim = TRUE, scale = "width", adjust = 0.5) +
      geom_boxplot(width=0.2, outlier.colour="grey25", notch = FALSE, notchwidth = .75, alpha = 0.75, colour = "grey10", outlier.size=0.3) +
      scale_fill_manual(values = color_subtypes)

   if(!is.null(highlight.samples)){
      sample_mean <- mean(subset(my_df,  taxonomy=="highlighted")$gex)
      g1 <- g1 + geom_hline(yintercept=sample_mean, linetype="dashed", color = "blue", size=2)
      }
   if(!is.null(individual.sample)){
      sample_gex <- my_df$gex[my_df$sample_id==individual.sample]
      g1 <- g1 + geom_hline(yintercept=sample_gex, linetype="solid", color = "red", size=1)
   }

   return(g1)
} # end geneplotTCGA_violin





#' wrapper for VlnPlot in Seurat package
#' @family plot wraps
#' @family seurat
#' @aliases VlnPlot.dl
#' @param gene.vec vector with gene symbols
#' @param orientation, prortrait (P) or Landscape (L)
#' @return plot to quartz
#' @export
VlnPlot_seurat_dl <- function(gene.vec, color.vec, x.text.size=12, y.title.size=12, ...){

   require(ggplot2)
   require(Seurat)
   pl <-  lapply(gene.vec, function(x){
      plx <- Seurat::VlnPlot(object = park,
                             features.plot = x, cols.use = color.vec,
                             y.log=T,
                             single.legend = T, nCol=1, remove.legend = T,
                             return.plotlist=T, do.return=F, ...)

      plx <- plx +
         theme(panel.background = element_rect(fill = 'gray95', colour = "gray50", linetype = "solid", size=0.25)) +
         ggtitle(label = "") +
         labs(x="", y=x) +
         theme(axis.title.y = element_text(angle=0, size=y.title.size)) +
         theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0, size=x.text.size)) +
         theme(plot.title = element_blank()) +
         theme(plot.margin = unit(c(0, 1, 0, 2.5), "cm"))
   }) # end lapply

   x.text <- ggplot_build(pl[[4]])$layout$panel_params[[1]]$x.labels
   x.n <- length(ggplot_build(pl[[4]])$layout$panel_params[[1]]$x.labels)

   pl[1] <- pl[1] %>%
      map(~ .x +
             #scale_x_discrete(position = "top", limits=c(0.5,x.n+0.5), breaks=c(1:x.n), labels=x.text)
             scale_x_discrete(position = "top") +
             theme(axis.text.x = element_text(vjust= -0.5))
      )

   pl[2:(length(pl)-1)] <- pl[2:(length(pl)-1)] %>%
      map(~ .x + theme(
         axis.text.x=element_blank(),
         axis.title.x=element_blank(),
         axis.ticks.x=element_blank()
      )
      )

   return(pl)
} # end VlnPlot.dl











#' boxplot single gene epxression in TCGA urinary tract (UTC) data
#' @family plot wraps
#' @family tcga
#' @family gene expression
#' @family eset
#' @param ensg A ensemble gene identifier (homo sapiens)
#' @return a plot
#' @export
#'
geneplotTcgaRccTax = function(ensg=NULL, ...){
   require(gplots)
   require(ggplot2)
   require(tidyr)
   require(tidyverse)
   require(dplyr)
   message("boxplot gene epxression in TCGA urinary tract (UTC) data")

   if(length(ensg)!=1){
      message(ensg, " only one id is allowed")
      return(NULL)
   }
   gene.i = match(ensg, featureNames(utc_nt))
   if(is.na(gene.i)){
      message(ensg, " not found in feature data for utc_nt.")
      return(NULL)
   }

   ## PAN UTC data
   ## --------------
   data('pdata_panCanGex', package='dlfoo2')
   pdata <- pdata_panCanGex  %>% dplyr::filter(sample_id%in%sampleNames(utc_nt))

   tax_vals <- c("kidney_n","ccpRCC","ccRCC","pRCC_e1","pRCC_e2","CIMP_p","CIMP_o","chRCC","chONC","pONC",
"rcc_UroMesInf","uro","outlier","bladder_t","bladder_n")

   my_df <- pdata %>%
      dplyr::select(sample_id, taxonomy_simplified, taxonomy_combo) %>%
      mutate(taxonomy_1 = recode(taxonomy_simplified, "kidney_n"="kidney_n", "ccpRCC"="RCC", "ccRCC"="RCC", "pRCC_e1"="pRCC_e1","pRCC_e2"="pRCC_e2","CIMP_p"="RCC","CIMP_o"="RCC","chRCC"="RCC","chONC"="RCC","pONC"="RCC",
"rcc_UroMesInf"="RCC","uro"="BLCA", "outlier"="RCC","bladder_t"="BLCA","bladder_n"="bladder_n")) %>%
      mutate(taxonomy_2 = recode(taxonomy_simplified, "kidney_n"="kidney_n", "ccpRCC"="ccpRCC", "ccRCC"="ccRCC", "pRCC_e1"="pRCC_e1","pRCC_e2"="pRCC_e2","CIMP_p"="CIMP","CIMP_o"="CIMP","chRCC"="chRCC", .default=NA_character_)) %>%
      mutate(taxonomy_2 = factor(taxonomy_2, levels=c("kidney_n", "ccRCC","pRCC_e1", "pRCC_e2", "CIMP","chRCC", "ccpRCC"))) %>%
      mutate(taxonomy_3 = if_else(as.character(taxonomy_2)=="ccRCC", as.character(taxonomy_combo), as.character(taxonomy_2))) %>%
      mutate(taxonomy_3 = factor(taxonomy_3, levels= c("kidney_n","ccpRCC", "CC-e.2","CC-e.1","CC-e.3","pRCC_e1","pRCC_e2","CIMP","chRCC")))

   my_df = my_df %>%
      mutate(gex=exprs(utc_nt)[gene.i,])


   ## PLOT CONDENSED GROUPS
   g <- ggplot(my_df) +
      aes(x=taxonomy_simplified, y=gex, fill=taxonomy_simplified) +
      theme(panel.background = element_rect(fill = 'gray95', colour = 'antiquewhite3', size = 0.5)) +
      theme(panel.grid.major =  element_line(colour = 'lightsteelblue2'), panel.grid.minor =  element_line(colour = 'mistyrose3')) +
      theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
      labs(y=paste0("log2 FPKM"), x="") +
      guides(fill=FALSE) +
      theme(axis.title.y = element_text(size=12)) +
      theme(axis.title.x = element_text(size=12)) +
      theme(axis.text = element_text(size=9)) +
      theme(title = element_text(size=9)) +
      theme(strip.text = element_text(size=12)) +
      theme(title = element_text(size=12))
   g1 <- g +
      geom_violin(alpha=0.5, width=1, trim = TRUE, scale = "width", adjust = 0.5) +
      geom_boxplot(width=0.2, outlier.colour="grey25", notch = FALSE, notchwidth = .75, alpha = 0.75, colour = "grey10", outlier.size=0.3) +
      scale_fill_manual(values = color_subtypes)

   g <- ggplot(my_df %>% filter(!is.na(taxonomy_3)) %>% filter(!taxonomy_3=="ccpRCC")) +
      aes(x=taxonomy_3, y=gex, fill=taxonomy_3) +
      theme(panel.background = element_rect(fill = 'gray95', colour = 'antiquewhite3', size = 0.5)) +
      theme(panel.grid.major =  element_line(colour = 'lightsteelblue2'), panel.grid.minor =  element_line(colour = 'mistyrose3')) +
      theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
      labs(y=paste0("log2 FPKM"), x="") +
      guides(fill=FALSE) +
      theme(axis.title.y = element_text(size=12)) +
      theme(axis.title.x = element_text(size=12)) +
      theme(axis.text = element_text(size=9)) +
      theme(title = element_text(size=9)) +
      theme(strip.text = element_text(size=12)) +
      theme(title = element_text(size=12))
   g2 <- g +
      geom_violin(alpha=0.5, width=1, trim = TRUE, scale = "width", adjust = 0.5) +
      geom_boxplot(width=0.2, outlier.colour="grey25", notch = FALSE, notchwidth = .75, alpha = 0.75, colour = "grey10", outlier.size=0.3) +
      scale_fill_manual(values = color_subtypes)

   return(list(g1,g2))

}





#' boxplot gene epxression in TCGA PAN tissues
#' @family plot wraps
#' @family tcga
#' @family gene expression
#' @family eset
#' @param ensg Ensemble gene id
#' @return list of ggplot objects
#' @export
#'
geneplotTcgaPan = function(ensg){
   require(gplots)
   require(ggplot2)
   require(tidyr)
   require(tidyverse)
   require(dplyr)
   message("boxplot gene epxression in TCGA PAN tissues")


   if(length(ensg)!=1){
      message(ensg, " only one id is allowed")
      return(NULL)
   }
   gene.i = match(ensg, featureNames(pan_gdc))
   if(is.na(gene.i)){
      message(ensg, " not found in feature data for TCGA ")
      return(NULL)
   }

   my_df <- pdata_gdc %>% dplyr::select(sample_id, sample_type3)
   stopifnot(identical(pdata_gdc$sample_id, sampleNames(pan_gdc)))

   my_df = my_df %>%
      mutate(gex=exprs(pan_gdc)[gene.i,])

   g <- ggplot(my_df) +
      aes(x=sample_type3, y=gex, fill=sample_type3) +
      theme(panel.background = element_rect(fill = 'gray95', colour = 'antiquewhite3', size = 0.5)) +
      theme(panel.grid.major =  element_line(colour = 'lightsteelblue2'), panel.grid.minor =  element_line(colour = 'mistyrose3')) +
      theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
      labs(y=paste0("log2 FPKM")) +
      theme(axis.title.y = element_text(size=12)) +
      theme(axis.title.x = element_text(size=12)) +
      theme(axis.text = element_text(size=9)) +
      theme(title = element_text(size=9)) +
      guides(fill=FALSE)+
      theme(strip.text = element_text(size=12)) +
      theme(title = element_text(size=12))
   g1 <- g +
       geom_violin(alpha=0.5, width=1, trim = TRUE, scale = "width", adjust = 0.5) +
       geom_boxplot(width=0.2, outlier.colour="grey25", notch = FALSE, notchwidth = .75, alpha = 0.75, colour = "grey50", outlier.size=0.3) +
      scale_fill_manual(values = color_subtypes)

   b <- my_df %>%
         group_by(sample_type3) %>%
         summarise(group_median=median(gex)) %>%
         arrange(-group_median)
   my_df2 <- my_df %>%
         mutate(sample_type3 = factor(as.character(sample_type3), levels=b$sample_type3))
   g <- ggplot(my_df2) +
         aes(x=sample_type3, y=gex, fill=sample_type3) +
         theme(panel.background = element_rect(fill = 'gray95', colour = 'antiquewhite3', size = 0.5)) +
         theme(panel.grid.major =  element_line(colour = 'lightsteelblue2'), panel.grid.minor =  element_line(colour = 'mistyrose3')) +
         theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
         labs(y=paste0("log2 FPKM ")) +
         theme(axis.title.y = element_text(size=12)) +
         theme(axis.title.x = element_text(size=12)) +
         theme(axis.text = element_text(size=9)) +
         theme(title = element_text(size=9)) +
         guides(fill=FALSE)+
         theme(strip.text = element_text(size=12)) +
         theme(title = element_text(size=12))
      g2 <- g +
          geom_violin(alpha=0.5, width=1, trim = TRUE, scale = "width", adjust = 0.5) +
          geom_boxplot(width=0.2, outlier.colour="grey25", notch = FALSE, notchwidth = .75, alpha = 0.75, colour = "black", outlier.size=0.3) +
         scale_fill_manual(values = color_subtypes)

      return(list(g1,g2))
} # end bxpTcgaPan




#' boxplot gene epxression in GTEX data set
#' @family plot wraps
#' @family gtex
#' @family gene expression
#' @family eset
#' @param ensg one ENSEMBL gene id
#' @return ggplot object
#' @export
geneplotGtex = function(ensg){
   require(gplots)
   require(ggplot2)
   require(tidyr)
   require(tidyverse)
   require(dplyr)
   message("boxplot gene epxression in GTEX data set")
   es_gtex_8555set = dltcgaDataGDC::es_gtex_8555set

   if(length(ensg)!=1){
      message(ensg, " only one id is allowed")
      return(NULL)}

   gene.i = match(ensg, featureNames(es_gtex_8555set))
   if(is.na(gene.i)){
      message(ensg, " not found in feature data for GTEX")
      return(NULL)
   }

   pdata <- as.tibble(pData(es_gtex_8555set))
   #table(pdata$SMTS)
   #head(fData(es_gtex_8555set))

   my_df <- pdata %>%
      mutate(gex = exprs(es_gtex_8555set)[gene.i, ])

   b <- my_df %>%
      group_by(SMTS) %>%
      summarise(group_median=median(gex)) %>%
      arrange(-group_median)
   my_df2 <- my_df %>%
      mutate(sample_type3 = factor(as.character(SMTS), levels=b$SMTS))

    g <- ggplot(my_df2) +
      aes(x=sample_type3, y=gex, fill=sample_type3) +
      theme(panel.background = element_rect(fill = 'gray95', colour = 'antiquewhite3', size = 0.5)) +
      theme(panel.grid.major =  element_line(colour = 'lightsteelblue2'), panel.grid.minor =  element_line(colour = 'mistyrose3')) +
      theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
      labs(y=paste0("log2 FPKM"), x="") +
      theme(axis.title.y = element_text(size=12)) +
      theme(axis.title.x = element_text(size=12)) +
      theme(axis.text = element_text(size=9)) +
      theme(title = element_text(size=9)) +
      guides(fill=FALSE)+
      theme(strip.text = element_text(size=12)) +
      theme(title = element_text(size=12))

   g1 <- g +
       geom_violin(alpha=0.5, width=1, trim = TRUE, scale = "width", adjust = 0.5) +
       geom_boxplot(width=0.2, outlier.colour="grey25", notch = FALSE, notchwidth = .75, alpha = 0.75, colour = "black", outlier.size=0.3) +
      scale_fill_manual(values = color_subtypes)
   return(g1)
}


#' correlation gene epxression in GTEX data set
#' @family analysis
#' @family gene expression
#' @family eset
#' @family gtex
#' @param ensg ENSEMBL id
#' @param tissue what tossie types to include in analysis. Defautls to "AdiposeTissue","AdrenalGland","BloodVessel","Bladder","Brain","Breast","Blood","Skin","CervixUteri","Colon","Esophagus","FallopianTube","Heart","Kidney","Liver","Lung","SalivaryGland","Muscle","Nerve","Ovary","Pancreas","Pituitary","Prostate","SmallIntestine","Spleen","Stomach","Testis","Thyroid","Uterus","Vagina"
#' @return table with correlations
#' @export
corrGtex = function(ensg, tissue=c("AdiposeTissue","AdrenalGland","BloodVessel","Bladder","Brain","Breast","Blood","Skin","CervixUteri","Colon","Esophagus","FallopianTube","Heart","Kidney","Liver","Lung","SalivaryGland","Muscle","Nerve","Ovary","Pancreas","Pituitary","Prostate","SmallIntestine","Spleen","Stomach","Testis","Thyroid","Uterus","Vagina")){
   require(tidyr)
   require(tidyverse)
   require(dplyr)
   message("correlation gene epxression in GTEX data set")
   tissue = match.arg(tissue, several.ok = T)
   # es_gtex_8555set = dltcgaDataGDC::es_gtex_8555set


   if(length(ensg)!=1){
      message(ensg, " only one id is allowed")
      return(NULL)}

   gene.i = match(ensg, featureNames(es_gtex_8555set))
   if(is.na(gene.i)){
      message(ensg, " not found in feature data for GTEX")
      return(NULL)
   }

   pdata <- as.tibble(pData(es_gtex_8555set))
   #unique(pdata$SMTS)
   my_df <- pdata %>%
      filter(SMTS %in% tissue)
   my_mat <- exprs(es_gtex_8555set[,my_df$SAMPID])
   my_df <- my_df %>%
      mutate(gex = my_mat[gene.i, ])

   my_cor <- apply(my_mat, 1, function(x) cor(x, my_df$gex))
   my_cor[is.na(my_cor)] <- 0
   my_cor_df <- as.tibble(data_frame(ENSG=featureNames(es_gtex_8555set), SYMBOL=fData(es_gtex_8555set)$SYMBOL, cor=my_cor)) %>% arrange(-cor)

   return(my_cor_df)
}





#' correlations in CCLE data set
#' @family analysis
#' @family gene expression
#' @family eset
#' @family ccle
#' @param ensg ENSEMBL gene id
#' @param tissue tissue types to include in analysis. Defualts to "prostate","stomach","urinary_tract","central_nervous_system","ovary","haematopoietic_and_lymphoid_tissue","kidney","thyroid","skin","soft_tissue","salivary_gland","lung","bone","pleura","endometrium","pancreas","breast","upper_aerodigestive_tract","large_intestine","autonomic_ganglia","oesophagus","fibroblast","liver","biliary_tract","small_intestine","cervix","haematopoietic_and_lymphoid","buccal","placenta","eye","adrenal_cortex"
#' @param geneSubset A vector of ensg gene ids to limit the analysis
#' @return correlation data frame
#' @export
#'
corrCcle = function(ensg, geneSubset=NULL, tissue.site=c("prostate","stomach","urinary_tract","central_nervous_system","ovary","haematopoietic_and_lymphoid_tissue","kidney","thyroid","skin","soft_tissue","salivary_gland","lung","bone","pleura","endometrium","pancreas","breast","upper_aerodigestive_tract","large_intestine","autonomic_ganglia","oesophagus","fibroblast","liver","biliary_tract","small_intestine","cervix","haematopoietic_and_lymphoid","buccal","placenta","eye","adrenal_cortex")){
   require(tidyr)
   require(tidyverse)
   require(dplyr)
   message("correlation gene epxression in CCLE data set")
   tissue.site = match.arg(tissue.site, several.ok = T)
   my_es = readRDS(file="/Volumes/MacPro2TB/RESOURCES/CCLE/Curated/ccle_1156set_rpkm_DepMap_18q3_20180718_log2.rds")

   #table(pData(es_ccle)$site.primary)
   if(length(ensg)!=1){
      message(ensg, " only one id is allowed")
      return(NULL)}
   if(!is.null(geneSubset)){
      my_subset <- intersect(featureNames(my_es), c(geneSubset,ensg))
      if(length(my_subset)==0){
         message(" no subset genes not identified. Check that ENSGs are input as a character vector")
         return(NULL)}
      my_es <- my_es[my_subset,]
   }
   gene.i = match(ensg, featureNames(my_es))
   if(is.na(gene.i)){
      message(ensg, " not found in feature data for ccle")
      return(NULL)
   }

   pdata <- as.tibble(pData(my_es))
   #unique(pdata$SMTS)
   my_df <- pdata %>%
      dplyr::filter(site.primary %in% tissue.site)

   my_mat <- exprs(my_es[,my_df$id])
   my_df <- my_df %>%
      mutate(gex = my_mat[gene.i, ])

   my_cor <- apply(my_mat, 1, function(x) cor(x, my_df$gex))
   my_cor[is.na(my_cor)] <- 0
   my_cor_df <- as.tibble(data_frame(ENSG=featureNames(my_es), SYMBOL=fData(my_es)$SYMBOL, cor=my_cor)) %>% arrange(-cor)

   return(my_cor_df)
}


#' boxplot gene epxression in CCLE data set
#' @family plot wraps
#' @family gene expression
#' @family eset
#' @family ccle
#' @param ensg ENSEMBL gene id
#' @return ggplot object
#' @export
#'
geneplotCcle1156 = function(ensg){
   require(Biobase)
   require(gplots)
   require(ggplot2)
   require(tidyr)
   require(tidyverse)
   require(dplyr)
   message("boxplot gene epxression in CCLE data set")
   ## TWO DATASETS TO CHOOSE FROM
   #ccle_1019set = dltcgaDataGDC::ccle_1019set
   my_es = readRDS("~/RESOURCES/CCLE/CURATED/ccle_1156set_rpkm_20180718_log2.rds")

   if(length(ensg)!=1){
      message(ensg, " only one id is allowed")
      return(NULL)}

   gene.i = match(ensg, featureNames(my_es))
   if(is.na(gene.i)){
      message(ensg, "not found in feature data for GTEX")
      return(NULL)
   }

   my_df <- as.tibble(pData(my_es)) %>%
   mutate(gex = exprs(my_es)[gene.i, ])

   b <- my_df %>%
      group_by(tissue) %>%
      summarise(group_median=median(gex)) %>%
      arrange(-group_median)
   my_df2 <- my_df %>%
      mutate(sample_type3 = factor(as.character(tissue), levels=b$tissue))

    g <- ggplot(my_df2) +
      aes(x=sample_type3, y=gex, fill=sample_type3) +
      theme(panel.background = element_rect(fill = 'gray95', colour = 'antiquewhite3', size = 0.5)) +
      theme(panel.grid.major =  element_line(colour = 'lightsteelblue2'), panel.grid.minor =  element_line(colour = 'mistyrose3')) +
      theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
      labs(y=paste0("log2 FPKM"), x="") +
       guides(fill=FALSE)+
      theme(axis.title.y = element_text(size=12)) +
      theme(axis.title.x = element_text(size=12)) +
      theme(axis.text = element_text(size=9)) +
      theme(title = element_text(size=9)) +
      theme(strip.text = element_text(size=12)) +
      theme(title = element_text(size=12))
   g1 <- g +
       geom_violin(alpha=0.5, width=1, trim = TRUE, scale = "width", adjust = 0.5) +
       geom_boxplot(width=0.2, outlier.colour="grey25", notch = FALSE, notchwidth = .75, alpha = 0.75, colour = "black", outlier.size=0.3) +
      scale_fill_manual(values = color_subtypes)

   ## Rankplot?
   # g <- ggplot(my_df2) +
   #    aes(x=sample_type3, y=gex, color=sample_type3) +
   #    theme(panel.background = element_rect(fill = 'gray95', colour = 'antiquewhite3', size = 0.5)) +
   #    theme(panel.grid.major =  element_line(colour = 'lightsteelblue2'), panel.grid.minor =  element_line(colour = 'mistyrose3')) +
   #    theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
   #    labs(y=paste0("log2 FPKM"), x="") +
   #    # guides(fill=FALSE)+
   #    theme(axis.title.y = element_text(size=12)) +
   #    theme(axis.title.x = element_text(size=12)) +
   #    theme(axis.text = element_text(size=9)) +
   #    theme(title = element_text(size=9)) +
   #    theme(strip.text = element_text(size=12)) +
   #    theme(title = element_text(size=12))
   #
   # g +
   #    aes(x=rank(-gex)) +
   #    geom_point() +
   #    scale_color_manual(values = color_subtypes)
   return(g1)

} # end bxp ccle




#' barplot gene epxression in Knepper Rat nephron data
#' @family plot wraps
#' @family gene expression
#' @family eset
#' @family knepper
#' @param entrez Entrez id for Rattus Norvegicus
#' @return ggplot object
#' @export
geneplotKnepper = function(entrez){
   require(gplots)
   require(ggplot2)
   require(tidyr)
   require(tidyverse)
   require(dplyr)
   message("barplot gene epxression in Knepper Rat nephron data")
   es.knepper.rnaseq = dltcgaDataGDC::es.knepper.rnaseq

   if(length(entrez)!=1){
      message(entrez, "only one id is allowed")
      return(NULL)}

   gene.i = match(entrez, featureNames(es.knepper.rnaseq))
   if(is.na(gene.i)){
      message(entrez, " not found in feature data for GTEX")
      return(NULL)}

   pdata = pData(es.knepper.rnaseq)
   pdata$sample_type3 = color_tab_celtypes$tax.name[match(pdata$sample_type, color_tab_celtypes$name)]
   pdata$sample_type3 <- factor(pdata$sample_type3, levels=as.character(unique(color_tab_celtypes$tax.name)))

   my_df <- as.tibble(pdata) %>%
      dplyr::select(geo_accession, name, tissue, sample_type, sample_type3) %>%
      mutate(gex = exprs(es.knepper.rnaseq)[gene.i, ]) %>%
      mutate(sample_type3 = droplevels(sample_type3))

   g <- ggplot(my_df) +
      aes(x=sample_type3, y=gex, fill=sample_type3) +
      theme(panel.background = element_rect(fill = 'gray95', colour = 'antiquewhite3', size = 0.5)) +
      theme(panel.grid.major =  element_line(colour = 'lightsteelblue2'), panel.grid.minor =  element_line(colour = 'mistyrose3')) +
      theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
      labs(y=paste0("log2 gex"), x="") +
      guides(fill=FALSE)+
      theme(axis.title.y = element_text(size=12)) +
      theme(axis.title.x = element_text(size=12)) +
      theme(axis.text = element_text(size=9)) +
      theme(title = element_text(size=9)) +
      theme(strip.text = element_text(size=12)) +
      theme(title = element_text(size=12)) +
      expand_limits(y=0)

   if(max(my_df$gex)<2) g = g + scale_y_continuous(limits=c(0,2))
   g1 <- g +
       #geom_violin(alpha=0.5, width=1, trim = TRUE, scale = "width", adjust = 0.5) +
      geom_boxplot( outlier.colour="grey25", notch = FALSE, alpha = 0.75, colour = "black", outlier.size=1) +
      scale_fill_manual(values = color_celltypes)
   return(g1)

   } # end barplot knepper




#' plot clinical data TCGA PAN tissues
#' @family plot wraps
#' @family gene expression
#' @family eset
#' @family tcga
#' @family clinical
#' @param ensg Ensembl ID Human, e.g. from featureGetBM
#' @return ggplot object
#' @export
#'
geneplotTcgaClinical = function(ensg){
   require(gplots)
   require(ggplot2)
   require(tidyr)
   require(tidyverse)
   require(dplyr)
   require(survminer)
   require(survival)
   message("plot clinical data TCGA PAN tissues")


   if(length(ensg)!=1){
      message(ensg, " only one id is allowed")
      return(NULL)
   }
   gene.i = match(ensg, featureNames(pan_gdc))
   if(is.na(gene.i)){
      message(ensg, " not found in feature data for TCGA ")
      return(NULL)
   }

   my_df <- pdata_gdc
   my_df <- my_df %>%
      mutate(gex=exprs(pan_gdc)[gene.i,]) %>%
      filter(sample_type=="Primary Tumor") # A tibble: 9,528 x 20
   my_es <- pan_gdc[,my_df$sample_id] # assayData: 19676 features, 2007 samples


   ## Match up w tcga pan clinical data
   #clinical_pan <- readRDS(file = "/Users/david/RESOURCES/MSig_DL_Samples/pantcga_clinical.rds")
   data('clinical_pan', package='dlfoo2')
   table(clinical_pan$vital_status)
   clinical_pan <- clinical_pan %>%
      rename(patient_barcode = bcr_patient_barcode) %>%
      filter(patient_barcode %in% my_df$patient_barcode) %>%
      filter(vital_status != "[Discrepancy]")

    my_patients <- intersect(clinical_pan$patient_barcode, my_df$patient_barcode)

   my_df <- my_df %>%
      dplyr::select(-c(gender, vital_status, days_to_birth, days_to_death, days_to_last_followup, follow_up_version)) %>%
      left_join(clinical_pan, by="patient_barcode")

   surv_tab <- my_df %>%
               mutate(sample_type3 = factor(sample_type3)) %>%
               arrange(sample_type3) %>%
               dplyr::select(sample_id, patient_barcode, sample_type3, gender, vital_status, days_to_death, days_to_last_followup,
                  pathologic_T,pathologic_M, clinical_M, pathologic_N, clinical_stage, clinical_T, clinical_N, pathologic_stage,
                  neoplasm_histologic_grade, nuclear_grade_III_IV, gex) %>%
               group_by(sample_type3) %>%
                filter(!is.na(vital_status)) %>%
               mutate(fu_status = as.numeric(recode(vital_status, "Alive"="0", "Dead"="1"))) %>%
               filter(!all(is.na(days_to_death) & is.na(days_to_last_followup))) %>%
               mutate(fu_time = as.numeric(if_else(fu_status=="0", days_to_last_followup, days_to_death)))

   ## Survplot - Kaplan
   g_surv = lapply(levels(surv_tab$sample_type3), function(x){
      surv_tab2 = surv_tab %>%
         filter(sample_type3 %in% x) %>%
         mutate(gex.bin = findInterval(gex, vec=quantile(gex, probs=c(0,0.5,1)), rightmost.closed = T)) %>%
         mutate(gex.bin = recode(gex.bin, "1"="low","2"="high")) %>%
         mutate(gex.bin = factor(gex.bin, levels=c("low","high")))

      fit <- survfit(Surv(time = fu_time, event = fu_status) ~ gex.bin, data = surv_tab2)

      ggsurv <- ggsurvplot(fit, data = surv_tab2,
         pval = TRUE, conf.int = TRUE,
         palette = c("#2E9FDF", "chocolate3"), xlab = x, break.time.by = 365,
         conf.int.style = "step", size=1.5,
         ggtheme = theme(
            axis.title=element_text(size=14),
            axis.text.y=element_text(size=12),
            axis.text.x=element_text(angle = 45, hjust = 1, size=12),
            text = element_text(size = 12))
         )
      return(ggsurv)
   })
   names(g_surv) <- levels(surv_tab$sample_type3)


   ## Boxplots
   # Stage
   stage_tab <- surv_tab %>%
      filter(!is.na(pathologic_stage)) %>%
      filter(!pathologic_stage%in%c("[Discrepancy]","[Unknown]")) %>%
      ungroup() %>%
      mutate(sample_type3 = droplevels(sample_type3))

   g <- ggplot(stage_tab) +
      aes(x=pathologic_stage, y=gex, fill=pathologic_stage) +
      theme(panel.background = element_rect(fill = 'gray95', colour = 'antiquewhite3', size = 0.5)) +
      theme(panel.grid.major =  element_line(colour = 'lightsteelblue2'), panel.grid.minor =  element_line(colour = 'mistyrose3')) +
      labs(y=paste0("log2 FPKM")) +
      guides(fill=FALSE) +
      theme(axis.title.y = element_text(size=12)) +
      theme(axis.title.x = element_text(size=12)) +
      theme(axis.text = element_text(size=9)) +
      theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
      theme(title = element_text(size=9)) +
      theme(strip.text = element_text(size=12)) +
      theme(title = element_text(size=12)) +
      facet_wrap(~sample_type3, scales = "free_x", drop=T, nrow=3)
   g_stage <- g +
      geom_violin(alpha=0.5, width=1, trim = TRUE, scale = "width", adjust = 0.5) +
      geom_boxplot(width=0.2, outlier.colour="grey25", notch = FALSE, notchwidth = .75, alpha = 0.75, colour = "grey50", outlier.size=0.3) +
      scale_fill_manual(values = rainbow(20))
      #RColorBrewer::brewer.pal(9,"Greens")[-c(1:2)]


   # Grade neoplasm_histologic_grade
    grade_tab <- surv_tab %>%
      filter(!is.na(neoplasm_histologic_grade)) %>%
      filter(!neoplasm_histologic_grade%in%c("[Discrepancy]","[Unknown]")) %>%
      ungroup() %>%
      mutate(sample_type3 = droplevels(sample_type3))

    g_colors <- c(rev(RColorBrewer::brewer.pal(5, "Spectral")), "gray","darkorange1","darkolivegreen1")
    names(g_colors) = sort(unique(grade_tab$neoplasm_histologic_grade))
   g <- ggplot(grade_tab) +
      aes(x=neoplasm_histologic_grade, y=gex, fill=neoplasm_histologic_grade) +
      theme(panel.background = element_rect(fill = 'gray95', colour = 'antiquewhite3', size = 0.5)) +
      theme(panel.grid.major =  element_line(colour = 'lightsteelblue2'), panel.grid.minor =  element_line(colour = 'mistyrose3')) +
      labs(y=paste0("log2 FPKM")) +
      guides(fill=FALSE) +
      theme(axis.title.y = element_text(size=12)) +
      theme(axis.title.x = element_text(size=12)) +
      theme(axis.text = element_text(size=9)) +
      theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
      theme(title = element_text(size=9)) +
      theme(strip.text = element_text(size=12)) +
      theme(title = element_text(size=12)) +
      facet_wrap(~sample_type3, scales = "free_x", drop=T, nrow=3)
   g_grade <- g +
      geom_violin(alpha=0.7, width=1, trim = TRUE, scale = "width", adjust = 0.5) +
      geom_boxplot(width=0.2, outlier.colour="grey25", notch = FALSE, notchwidth = .75, alpha = 0.85, colour = "grey50", outlier.size=0.3) +
      scale_fill_manual(values = g_colors)
      #RColorBrewer::brewer.pal(9,"Greens")[-c(1:2)]

  return(list(g_surv, g_stage, g_grade))
}


      #  pdata <- pdata_panCanFull %>% dplyr::filter(!duplicated(sample_id)) %>% dplyr::filter(grepl("ccRCC|pRCC|chRCC|CIMP|ccpRCC|mesRCC|chONC|pONC", tax_simp)) %>% droplevels()
   # surv.tab = pdata
   # colnames(pdata)
   # survPlot_pdata(pdata = as.data.frame(pdata),  surv.bin.column = "OS", surv.time.column = "OS.time", annot.column = "tax_simp2", title.main="OS tax_simp2", color.key = color_subtypes)


#' Generic for plotting survival data from sample table separated on defined groups
#' @family plot wraps
#' @family clinical
#' @family survival
#' @family kaplan
#' @param pdata sample table
#' @param surv.bin.column colunmn name fopr surv bin
#' @param surv.time.column column in pdata for surv time
#' @param annot.column name of annotations to use
#' @param title plot title
#' @param color.key what color key
#' @param conf.int if to add confindence intervals
#' @return ggplot object
#' @export
#'
survPlot_pdata = function(
   surv.tab,
   surv.bin.column =  NULL,
   surv.time.column = NULL,
   annot.column = NULL,
   color.key = dlfoo2::color_subtypes,
   conf.int = TRUE,
   title.main=""
   ){
   require(dplyr)
   require(gplots)
   require(ggplot2)
   require(tidyr)
   require(tidyverse)
   require(dplyr)
   require(survminer)
   require(survival)

   x <- data.frame(surv.tab)
   vars = c(surv_bin=surv.bin.column, surv_time=surv.time.column, annot=annot.column)
   x <- x %>% dplyr::rename(!!vars)

   str(x)
   class(x)
   x$surv_bin <- as.numeric(x$surv_bin)
   x$surv_time <- as.numeric(x$surv_time)


   message("... removing NA rows. from: ", nrow(x))
   x <- x[!is.na(x$surv_bin), ]
   x <- x[!is.na(x$surv_time), ]
   x <- x[!is.na(x$annot), ]
   message("... ... to: ", nrow(x))

   ## create fit object
   fit <- survfit(Surv(time = surv_time, event = surv_bin) ~ x$annot, data = x)
   str(fit)
   names(fit$strata) <- gsub(".*=","",names(fit$strata))

   # plot.
   ggsurv <- ggsurvplot(fit, data = x,
         pval = TRUE, conf.int = conf.int,
         palette = color.key, xlab = "time", break.time.by = 365,
         conf.int.style = "step", size=1.5,
         ggtheme = theme(
            axis.title=element_text(size=14),
            axis.text.y=element_text(size=12),
            axis.text.x=element_text(angle = 45, hjust = 1, size=12),
            text = element_text(size = 12))
         )
   # ggsurv
   return(ggsurv)
}





#' plot survival data TCGA PAN tissues. Select samples using my_samples and ensg for gene
#' @family plot wraps
#' @family gene expression
#' @family eset
#' @family tcga
#' @family clinical
#' @param ensg Ensembl ID Human, e.g. from featureGetBM
#' @param my_samples sample names
#' @return ggplot object
#' @export
#'
tcgaSurvPlot2 = function(ensg, my_samples=NULL){
   require(gplots)
   require(ggplot2)
   require(tidyr)
   require(tidyverse)
   require(dplyr)
   require(survminer)
   require(survival)
   message("plot clinical data TCGA PAN tissues")


   if(length(ensg)!=1){
      message(ensg, " only one id is allowed")
      return(NULL)
   }
   data('pan_gdc', package='dlfoo2data')
   gene.i = match(ensg, featureNames(pan_gdc))
   if(is.na(gene.i)){
      message(ensg, " not found in feature data for TCGA ")
      return(NULL)
   }

   data('pdata_panCanFull', package='dlfoo2data')
   my_df <- pdata_pan
   my_df <- my_df %>%
      mutate(gex=exprs(pan_gdc)[gene.i,]) %>%
      filter(sample_type=="Primary Tumor") # A tibble: 9,528 x 20
   if(!is.null(my_samples)) my_df <- my_df %>% dplyr::filter(sample_id %in% my_samples)

         surv_tab2 = my_df %>%
         # filter(sample_type3 %in% x) %>%
         mutate(gex.bin = findInterval(gex, vec=quantile(gex, probs=c(0,0.5,1)), rightmost.closed = T)) %>%
         mutate(gex.bin = recode(gex.bin, "1"="low","2"="high")) %>%
         mutate(gex.bin = factor(gex.bin, levels=c("low","high")))

      fit <- survfit(Surv(time = fu_time, event = fu_status) ~ gex.bin, data = surv_tab2)

      ggsurv <- ggsurvplot(fit, data = surv_tab2,
         pval = TRUE, conf.int = TRUE,
         palette = c("#2E9FDF", "chocolate3"), xlab = "", break.time.by = 365,
         conf.int.style = "step", size=1.5,
         ggtheme = theme(
            axis.title=element_text(size=14),
            axis.text.y=element_text(size=12),
            axis.text.x=element_text(angle = 45, hjust = 1, size=12),
            text = element_text(size = 12))
         )

   return(ggsurv)
}



#' plot clinical data TCGA PAN tissues
#' @family plot wraps
#' @family tcga
#' @family clinical
#' @family survival
#' @param surv_tab data frame (or tibble) with patient_barcode, fu_time, fu_status, and fu_group
#' @return ggplot object
#' @export
tcgaSurvPlot <- function(surv_tab){

   require(survminer)
   require(survival)

   if(!all(c("patient_barcode","fu_group","fu_status","fu_time") %in% colnames(surv_tab))) stop()


   fit <- survfit(Surv(time = fu_time, event = fu_status) ~ fu_group, data = surv_tab)

   # png(file=paste0("Survplots/",my_tissues[i],".png"), width = 596, height=596)

      ggsurv <- ggsurvplot(fit, data = surv_tab,
         #pval = TRUE, conf.int = TRUE,
         # palette = c("#2E9FDF", "chocolate3"),
         xlab = "Time in days", break.time.by = 365,
         conf.int.style = "step", size=1.5,
         ggtheme = theme(
            axis.title=element_text(size=14),
            axis.text.y=element_text(size=12),
            axis.text.x=element_text(angle = 45, hjust = 1, size=12),
            text = element_text(size = 12)))
      return(ggsurv)
      #dev.off()
      }







