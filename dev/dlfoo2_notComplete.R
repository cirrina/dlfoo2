


#' Generic plot for co-occurence matrix
#' @description input is dataframe with, group, term, n, and p columns
#' @description  https://bioconductor.org/packages/release/bioc/manuals/fgsea/man/fgsea.pdf
#' @description plot order is determined by clustering - if not then set term column to factor and let levels decide order
#' @family enrichment
#' @family plotwrap
#' @family msig
#' @family gsea
#' @family fgsea
#' @param ea.tab, data frame with enrichment analysis to plot.
#' @param pval, p value
#' @param size, size factor for bubbles, can correspond to n or NES etc.. Allowed can negative, e.g. NES results from gsea (fgsea) and then be colored accordingly
#' @param pval.pal, color palette for p-val. if negative sizes (n).
#' @param order.by.clustering, let order of terms/groups be decided by clustering
#' @param reverse.x if to reverse x-axis
#' @param reverse.y if to reverse y-axis
#' @return ggplot object
#' @export
coOccurrence_bubble <- function(
      ea.tab=NULL, size, pval,
      order.by.clustering=T,
      pval.pal = c(rev(RColorBrewer::brewer.pal(9,"YlOrRd")[3:9]), RColorBrewer::brewer.pal(9,"PuBuGn")[3:9]),
      reverse.x=F, reverse.y=F
   ){
   require(tibble)
   require(tidyr)
   require(ggplot2)



my.k <- 25
my.cc <- dlfoo::ConsensusCluster_process(cc.dir = "~/PROJECTS/RCC_ccpRCC/ConsensusClust/ConsensusCluster_utc_nt_1435set_reps.1000_maxK.25_innerLinkage.ward_sd.rank.5000/", k = my.k)
pdata <- pData(utc_nt)
identical(pdata$sample_id, my.cc$sample.names)

#png(file="~/PROJECTS/RCC_ccpRCC/ConsensusClust/CoOccurence_TAX_Consensus_1435_k25.png",width = 842, height = 842)
pdf(file="~/PROJECTS/RCC_ccpRCC/ConsensusClust/CoOccurence_TAX_Consensus_1435_k25.pdf", width=11.692*2, height=11.692)

   table(my.cc$cc, pdata$taxonomy)
   my_tab <- ftable(paste0("cc-",my.cc$cc), pdata$taxonomy)
   my_tab <- as.matrix(my_tab)
   my_tab <- as.data.frame(my_tab, stringsAsFactors=F)
   #rownames(my_tab) <- paste0("cc-",rownames())
   # my_tab <- my_tab[c(4,8,2,5,6,1,3,7,9),c(2,1,3:9)]
   my_tab <- t(my_tab)
   my_tab <- my_tab[-1,]
   my_tab <- my_tab[c(2,3,4,16,17,18,19,5,14,7,6,10,11,5,12,21,8,1,20,9,13,15),]


   par(mfrow=c(1,2))
   my.mat <- as.matrix(apply(my_tab, 1, function(x) round(x/sum(x)*100,0)))
   heatcol.mat <-  heatcol_creator(my.mat, pal = "white.green", palette.sat = c(0,100))
   color.mat <- heatcol.mat$col
   par(mar=c(12,6,6,2))
   plot.new()
   plot.window(xlim = c(0, ncol(color.mat)), ylim = c(0, nrow(color.mat)))
   rect.args <- list(
            xleft = rep( seq(0, ncol(color.mat)-1), nrow(color.mat)),
            ybottom = rep(seq(nrow(color.mat)-1,0), each=ncol(color.mat)),
            xright = rep( seq(1, ncol(color.mat)), nrow(color.mat)),
            ytop = rep(seq(nrow(color.mat),1), each=ncol(color.mat)),
            col = as.character(t(color.mat)),
            border = NA)
   do.call("rect", args = rect.args)
   rect(xleft = 0, ybottom = 0, ytop = nrow(color.mat), xright = ncol(color.mat), lwd=0.5, border = "gray")
   axis(side=1, labels = colnames(my.mat), at=1:ncol(my.mat)-0.5, las=2)
   axis(side=2, labels = rownames(my.mat), at=nrow(my.mat):1-0.5, las=2)
   for(i in 1:ncol(my.mat)){
   do.call(text, args = list(x=rep(i,ncol(my.mat))-0.5, y=nrow(my.mat):1-0.5, labels=my.mat[,i], col=c("black","white")[findInterval(my.mat[,i], vec = c(0,66,100), all.inside = T) ]))
   }

   my.mat <- as.matrix(apply(my_tab, 2, function(x) round(x/sum(x)*100,0)))
   my.mat <- t(my.mat)
   heatcol.mat <-  heatcol_creator(my.mat, pal = "white.green", palette.sat = c(0,100))
   color.mat <- heatcol.mat$col
   plot.new()
   plot.window(xlim = c(0, ncol(color.mat)), ylim = c(0, nrow(color.mat)))
   rect.args <- list(
            xleft = rep( seq(0, ncol(color.mat)-1), nrow(color.mat)),
            ybottom = rep(seq(nrow(color.mat)-1,0), each=ncol(color.mat)),
            xright = rep( seq(1, ncol(color.mat)), nrow(color.mat)),
            ytop = rep(seq(nrow(color.mat),1), each=ncol(color.mat)),
            col = as.character(t(color.mat)),
            border = NA)
   do.call("rect", args = rect.args)
   rect(xleft = 0, ybottom = 0, ytop = nrow(color.mat), xright = ncol(color.mat), lwd=0.5, border = "gray")
   axis(side=1, labels = colnames(my.mat), at=1:ncol(my.mat)-0.5, las=2)
   axis(side=2, labels = rownames(my.mat), at=nrow(my.mat):1-0.5, las=2)
   for(i in 1:ncol(my.mat)){
   do.call(text, args = list(x=rep(i,ncol(my.mat))-0.5, y=nrow(my.mat):1-0.5, labels=my.mat[,i], col=c("black","white")[findInterval(my.mat[,i], vec = c(0,66,100), all.inside = T) ]))
   }
dev.off()

   ea.tab <- data.frame(ea.tab)
   if(!all(sapply(c("group","term","size","pval"), function(x)  x %in% colnames(ea.tab)))) {message("colnames must include: group, term, n, padj")}

   if(!is.factor(ea.tab$term)) ea.tab$term <- factor(ea.tab$term)
   if(!is.factor(ea.tab$group)) ea.tab$group <- factor(ea.tab$group)

   ea.tab <- droplevels(ea.tab)
   my_terms <- levels(ea.tab$term)
   my_groups <- levels(ea.tab$group)

   my_tab <- as_tibble(ea.tab) %>%
      mutate(logp = -log10(pval)) %>%
      mutate(logp.sat = if_else(logp>10, 10, logp)) %>%
      mutate(logp.sat = if_else(size<0, (-1)*logp.sat, logp.sat)) %>%
      arrange(group) %>%
      mutate(x0 = as.numeric(term)) %>%
      mutate(y0 = as.numeric(group))



   # if to order rows/cols by cluster analysis
   if(order.by.clustering){
      my_mat <- matrix(nrow=max(my_tab$y0), ncol=max(my_tab$x0), data=0, dimnames=list(levels(my_tab$group), levels(my_tab$term)))
      for(i in 1:nrow(my_tab)){
         my_mat[my_tab$y0[i], my_tab$x0[i]] <- my_tab$logp.sat[i]
      }

      mydist=function(c) {amap::Dist(c,method="pearson")}
      myclust=function(c) {hclust(c,method="ward.D2")}
      hca_rows <- my_mat %>% mydist %>% myclust %>% as.dendrogram
      hca_cols <- t(my_mat) %>% mydist %>% myclust %>% as.dendrogram

      my_tab <- my_tab %>%
         mutate(group = factor(group, levels=labels(hca_rows))) %>%
         mutate(term = factor(term, levels=labels(hca_cols))) %>%
         arrange(group, term)
   }

   if(reverse.x) my_tab <- my_tab %>% mutate(term = factor(term, levels=rev(levels(term)))) %>% arrange(group, term)
   if(reverse.y) my_tab <- my_tab %>% mutate(group = factor(group, levels=rev(levels(group)))) %>% arrange(group, term)


   str(unique(my_tab$group))

   x_labs <- levels(my_tab$term)
   y_labs <- levels(my_tab$group)
   x.n <- length(x_labs)
   y.n <- length(y_labs)
   geom_point_color = my_tab$logp.sat

   g <- ggplot(data = my_tab,
         aes(y=as.numeric(group), x=as.numeric(term), size=abs(as.numeric(size)))) +
         coord_fixed() +
         theme(panel.grid.major =  element_line(colour = 'lightsteelblue2'), panel.grid.minor =  element_line(colour = 'mistyrose3')) +
         theme(panel.background = element_rect(fill = 'mintcream', colour = NA)) +
         theme(panel.grid.major =  element_line(colour = NA), panel.grid.minor =  element_line(colour = 'black')) +
         # theme(text=element_text(family = "Arial")) +
         labs(x="", y="" ) +
         scale_x_continuous(limits=c(0.5, x.n+0.5), breaks=c(1:x.n), labels=label_wrap_mod(as.character(x_labs)), expand=c(0,0)) +
         scale_y_continuous(limits=c(0.5, y.n+0.5), breaks=c(1:y.n), labels=label_wrap_mod(as.character(y_labs)), expand=c(0,0)) +
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
         #scale_color_gradientn(name="-log10*(p-value)", colours = rev(pval_pal), values=c(0,0.1,1)) +
         scale_color_gradientn(name="-log10*(p-value)", colours = rev(pval.pal), values=c(0,0.5,1)) +
         geom_point(aes(colour=my_tab$logp.sat), alpha=1) +
         #scale_size_continuous(range=c(4,10)) +
         scale_size_continuous(range=c(3,10)) +
         guides(size=guide_legend(title="relative size", direction="horizontal"))

  return(g)
  #ggsave(g, filename = "Correlation/TCGA/EA_plot_fgsea_tcga.pdf", width=18, height=10, useDingbats = F)



}





#' Generate ramped color palettes as color input for heatmaps and ggplot gradients - values and ramps described by vectors with llength n
#' This should be replaced by plalette_gradients and palette_gradientRamp
#' @family plot functions
#' @family misc
#' @family older
#' @param x a matrix containg expression values


#devtools::use_vignette("dlfoo2_vignette")

# devtools::use_vignette("my-vignette")
# setwd("~/R_PKG/dlfoo2/")
# proj_set(rstudioapi::getActiveProject())
# require(usethis)
# setwd(proj_get())
# usethis::use_vignette("~/R_PKG/dlfoo2/")

# ● Your working directory is not the same as the active usethis project.
#   To set working directory to the project: `setwd(proj_get())`
#   To activate project in working directory: `proj_set(getwd())`
# ● Your active RStudio Project is not the same as the active usethis project.
#   To set usethis project to RStudio Project: `proj_set(rstudioapi::getActiveProject())`

usethis::use_vignette("dlfoo2", title = "dlfoo2")



#' #' Wrapper for umap analysis
#' #' will run analysis using multiple knns
#' #' umap Aanalysis using \code{\link[umap]{umap}} function
#' #' @family umap
#' #' @param x matrix with expression data
#' #' @return list of umap results
#' #' @export
#' #'
#' umapWrapper<-function(
#'     x=NULL, #matrix with expression data
#'     theta=0.5,
#'     max_iter=5000,
#'     outfile=NULL,
#'     floor_value=NULL,
#'     z_score=F,
#'     log=F,
#'     variance=NULL,
#'     perplex=c(5:50)
#'    ){
#'
#'    require(Rtsne)
#'    if(log){
#'       x<-log(x,2)
#'       x[which(x<0)] <-0
#'    }
#'
#'     if(!is.null(floor_value)){
#'       x[which(x<floor_value)]<-floor_value
#'     }
#'
#'    if(!is.null(variance)){
#'          yy<-apply(x,1,function(y)var(y))
#'          ord<-order(yy, decreasing = T)
#'          ord<-ord[c(1:(length(ord)*variance))]
#'          x<-x[ord,]
#'        }
#'    if(z_score){
#'          x<-t(apply(x,1,function(y)((y)-mean(y))/sd(y)))
#'          x[which(is.na(x))]<-0
#'     }
#'
#'     number<-dim(x)[1]
#'     #How many samples are in the matrix
#'     #min_perplex<-400 # As suggested by author, minimal perplexity=5
#'     #max_perplex<-401#round((number-1)/3,digits=0) #Set the max perplexity to (number-1)/3
#'     #if(max_perplex>50){max_perplex=50} #set a maximum perplexity to try
#'
#'     ## res<-list(perplex)
#'     res <- vector(mode = "list", length=length(perplex))
#'     names(res) <- paste0("perplex_",perplex)
#'
#'     cat("#T-SNE result file from Rtsne\n", file=outfile, append=FALSE)
#'     for(i in 1:length(perplex)) {
#'       set.seed(42) #set the seed for every iteration
#'       res[[i]] <- Rtsne::Rtsne(as.matrix(x), theta=theta, max_iter=max_iter, perplexity=perplex[i], check_duplicates = F) #Perform Rtsne function
#'       rownames(res[[i]]$Y)<-rownames(x) #Add the samplename names
#'
#'       #Header and result, cat sets the first col
#'       cat("\n#H:perplexity_",perplex[i],"\t", file=outfile, append=TRUE)
#'       write.table(res[[i]]$Y, file=outfile, col.names = F, quote=F, sep="\t", append=T)
#'     }
#'     #rtsne.data<-list(min_perplex=min(perplex), max_perplex=max(perplex), res=res)
#'     #return(rtsne.data)
#'     return(res)
#'    }





# hca_tab <- read.delim(file="~/PROJECTS/RCC_ccpRCC/Limma/HCA_Limma_Genes/hca_ccpRCC_Limma_v3_B_n15_fc4_EuclWard_tree_seed3.txt", as.is=T)
# gene.list <- split(hca_tab[,"ENSG"], hca_tab[,"k_name"])
# my_genes <- gene.list[[9]]
# tssGR <- genes2tssEnsg(ensgs=ensgs, genome.version = 'hg19')
# hifGR <- readRDS("~/RESOURCES/ChIPseq/ProcessedData/hif_smythies/hif_smythies_GRangesList.rds")
# myGR <- hifGR[[1]]
# plotLinesGR <- hifGR[c(1,6)]


#' \code{tssLinePlot}
#' NOT COMPLETE !!
#' @param tssGR A vector of >1 gene ensembl ids
#' @param upstream
#' @param downstream
#' @param plotLinesGR
#'
tssLinePlot = function(
   tssGR = NULL,
   upstream = 5e5,
   downstream = 5e5,
   plotLinesGR = NULL,
   plotMarksGR = NULL
   ){
   require(biomaRt)
   require(GenomicRanges)
   require(genoset)
   require(dplyr)
   require(tidyr)
   require(tidyverse)
   message("NOT COMPLETE!!!!!!!")

   # Get genome locations
   class(tssGR)
   class(plotLinesGR)
   #if(class(plotLinesGR)!="CompressedGRangesList") plotLinesGR<-GRangesList(plotLinesGR)

   ## message
   windowGR <- promoters(tssGR, upstream=5e4, downstream=5e4+1)
   tssPointGR <- promoters(tssGR, upstream=0, downstream=1)

   #findOverlaps(plotLinesGR[[1]], gr)@from

   ## Create (list of) data frames with line coordinates
   markerList <- lapply(names(plotLinesGR), function(z){
      #x <- plotLinesGR
      x <- plotLinesGR[[z]]
      my_overlaps <- findOverlaps(windowGR, x)
      if(length(my_overlaps)==0) return()
      y <- x[my_overlaps@to]
      my_df <- data.frame(
         name=names(tssPointGR[my_overlaps@from]),
         x1=start(y)-start(tssPointGR[my_overlaps@from]),
         x2=end(y)-start(tssPointGR[my_overlaps@from]),
         stringsAsFactors = F)

      my_df$marker <- z
      if("score" %in% colnames(elementMetadata(y))) my_df$score = y$score
      my_df <- as.tibble(dplyr::left_join(data.frame(name=names(tssPointGR), stringsAsFactors = F), my_df))
      my_df <- my_df %>% arrange(x1)
      my_df <- my_df %>% mutate(y1 = factor(name, levels = my_df$name[!duplicated(my_df$name)])) %>% mutate(y1 = as.integer(y1))
      my_df <- my_df %>% rowwise %>% mutate(xmid = mean(c(x1,x2)))
      return(my_df)
      })


   ggplot(markerList[[1]]) +
      geom_segment(aes(x=x1, xend=x2, y=y1, yend=y1))

   ggplot(markerList[[36]]) +
      geom_linerange(aes(ymin=x1, ymax=x2, x=y1)) +
      geom_point(aes(x=y1, y=xmid), size=2, shape=21,  fill="red", alpha=5/10) +

      coord_flip()



}


## NOT COMPLETED !!
dl_genomicPosCum <- function(x, all_seqlengths=F){
   if(!exists("cytoband.hg19")) load("~/RESOURCES/Annotation Files/cytoBand.Rdata")

   # calculate the cumulative genomic position
   stopifnot(class(x)=="GRanges")
   if(all_seqlengths) {chroms = as.character(seqlevels(x))}
   if(!all_seqlengths) {chroms = as.character(unique(seqnames(x)))}

   # if no seqlengths provided
   if(!length(seqlengths(x))){
      cat("\n ... no seqlengths provided - trying to addasdasd OT IMP")
   }

   # if seqlengths are supplied in GRanges object
   if(length(seqlengths(x))){
      my_seqlengths <- seqlengths(x)
      my_seqlengths <- my_seqlengths[chroms]
      # genome_add <- data.frame(wgpos=c(0,cumsum(as.numeric(seqlengths(x)))[-length(seqlengths(x))]), row.names = as.character(seqlevels(x)))
      genome_add <- data.frame(wgpos=c(0,cumsum(as.numeric(my_seqlengths))[-length(my_seqlengths)]), row.names = names(my_seqlengths))

      x$cum_start <- as.numeric(start(x)) + as.numeric(genome_add[as.character(seqnames(x)), ])
      x$cum_end <- as.numeric(start(x)) + as.numeric(genome_add[as.character(seqnames(x)), ])
   } # end  if seqlengths are provided
   return(x)
}


# ==============================================================================
#  Functions for parsing and filtering Illumina Methylation 450k data
# ==============================================================================




#' ATAC peak regions - overlap, function eaGeneList
#' Version 1.0 20190129
#' @param x.me
#' @param runName name of plot/analyses - subfolder
#' @param runFolder directory where to plot/save objects
#'
atacPeakOverlap <- function(
   x.me, runName, runFolder
   ){
   message("NOT COMPLETE")
   }




