







#' tsneME
#' @description  wrapper for RtSNE
#' @description  Version 1.0 20190213
#' @description tSNE Aanalysis using \code{\link[Rtsne]{Rtsne}} function in rtsne package.
#' @family tsne
#' @family analysis
#' @param x.me Genoset object, methylation beta values as x.me beta
#' @param runName name of plot/analyses - subfolder
#' @param runFolder higher level directory where to plot/save objects
#' @param perplex vector of integers. Suggested max is n-1 3
#' @param doPlotTaxonomy if to load rcc phenotype from file & plot
#' @param theta Rtsne theta value
#' @param max_iter Rtsne max_iter value
#' @param floor_value Rtsne floor_value
#' @param z_score if to use z-scores (False default)
#' @param log if to log matrix(False default)
#' @param variance numeric between 0 and 1 describing the ratio of Features to keep (0.8 keeps the 80 percent highest vatying Features)
#' @return write tables and tSNEres to file
#' @export
tsneMe <- function(x.me, runName, runFolder,
   perplex = seq(5, 50, by = 5), doPlotTaxonomy=F,
   theta=0.5, max_iter=5000, floor_value=NULL, z_score=F,
   log=F, variance=NULL){

   # data('pdata_utc_nt', package="dlfoo2data")

   if(class(x.me)!="GenoSet") stop("gset is not of class 'GenoSet'")
   x <- x.me[,,"beta"]
   if(max(perplex) > ncol(x)/3-1) stop("max perplex should not exeed n/3-1")

   myDir <- file.path(runFolder, runName)
   if(dir.exists(myDir)) stop("Run already exists. Set alternative 'runName' or delete/move old run from disk")
   dir.create(myDir)

    if(log){
      x<-log(x,2)
      x[which(x<0)]<-0
      }
    if(!is.null(floor_value)){
      x[which(x<floor_value)]<-floor_value
      }
    if(!is.null(variance)){
      yy<-apply(x,1,function(y)var(y))
      ord<-order(yy, decreasing = T)
      ord<-ord[c(1:(length(ord)*variance))]
      x<-x[ord,]
      }
    if(z_score){
      x<-t(apply(x,1,function(y)((y)-mean(y))/sd(y)))
      x[which(is.na(x))]<-0
      }

    x <- t(x)
    number <- dim(x)[1]
    #How many samples are in the matrix
    res <- vector(mode = "list", length=length(perplex))
    names(res) <- paste0("perplex_",perplex)

    outfile <- file.path(myDir, paste0(runName,"_tSneResults.txt"))
    message("... #T-SNE result file from Rtsne\n", file=outfile, append=FALSE)
    cat("... #T-SNE result file from Rtsne\n", file=outfile, append=FALSE)

    for(i in 1:length(perplex)) {
      set.seed(42)
       #set the seed for every iteration
      res[[i]] <- Rtsne::Rtsne(as.matrix(x), theta=theta, max_iter=max_iter, perplexity=perplex[i], check_duplicates = F)
      #Perform Rtsne function
      rownames(res[[i]]$Y)<-rownames(x)
      #Add the samplename names

      #Header and result, cat sets the first col
      cat("\n#H:perplexity_",perplex[i],"\t", file=outfile, append=TRUE)
      write.table(res[[i]]$Y, file=outfile, col.names = F, quote=F, sep="\t", append=T)

      } # end i all peplexes

     ## Plot standard taxonomy
   ## ----------------------
   if(doPlotTaxonomy){
      g_tsne <- dlfoo2::tsnePlot(res)
      pdf(gsub(".txt","TaxAllPerplex.pdf",outfile), width=18, height = 12, useDingbats = F)
         print(g_tsne)
      dev.off() } # end if plot

    saveRDS(res, file=gsub(".txt",".rds",outfile))
    write.table(data.frame(FeatureID=rownames(x)), file=gsub(".txt","_FeatureIDs.txt",outfile), row.names = F, quote=F)

   } # end me tSNE



#' function tsneGex wrapper for RtSNE
#' @description Version 1.0 20190418
#' @description tSNE Aanalysis using \code{\link[Rtsne]{Rtsne}} function in stats package.
#' @family tsne
#' @family analysis
#' @param eset expressionSet object
#' @param runName name of plot/analyses - subfolder
#' @param runFolder higher level directory where to plot/save objects
#' @param perplex vector of integers. Suggested max is n-1  3
#' @param doPlotTaxonomy if to load rcc phenotype from file & plot
#' @param theta Rtsne theta value
#' @param max_iter Rtsne max_iter value
#' @param floor_value Rtsne floor_value
#' @param z_score if to use z-scores (False default)
#' @param log if to log matrix(False default)
#' @param variance numeric between 0 and 1 describing the ratio of Features to keep (0.8 keeps the 80 percent highest vatying Features)
#' @return write tables and tSNEres to file
#' @export
tsneGex <- function(eset, perplex = seq(5, 50, by = 5), runName, runFolder,
   doPlotTaxonomy=F, theta=0.5, max_iter=5000, floor_value=NULL, z_score=F,
   log=F, variance=NULL){
   require(Rtsne)

   # data('pdata_utc_nt', package="dlfoo2data")

   if(class(eset)!="ExpressionSet") stop("eset not of class 'GenoSet'")
   x <- exprs(eset)
   if(max(perplex) > ncol(x)/3-1) stop("max perplex should not exeed n/3-1")

   myDir <- file.path(runFolder, runName)
   if(dir.exists(myDir)) stop("Run already exists. Set alternative 'runName' or delete/move old run from disk")
   dir.create(myDir)

    if(log){
      x<-log(x,2)
      x[which(x<0)]<-0
      }
    if(!is.null(floor_value)){
      x[which(x<floor_value)]<-floor_value
      }
    if(!is.null(variance)){
      yy<-apply(x,1,function(y)var(y))
      ord<-order(yy, decreasing = T)
      ord<-ord[c(1:(length(ord)*variance))]
      x<-x[ord,]
      }
    if(z_score){
      x<-t(apply(x,1,function(y)((y)-mean(y))/sd(y)))
      x[which(is.na(x))]<-0
      }

    x <- t(x)
    number <- dim(x)[1]
    #How many samples are in the matrix
    res <- vector(mode = "list", length=length(perplex))
    names(res) <- paste0("perplex_",perplex)

    outfile <- file.path(myDir, paste0(runName,"_tSneResults.txt"))
    message("... #T-SNE result file from Rtsne\n", file=outfile, append=FALSE)
    cat("... #T-SNE result file from Rtsne\n", file=outfile, append=FALSE)

    for(i in 1:length(perplex)) {
      set.seed(42)
       #set the seed for every iteration
      res[[i]] <- Rtsne::Rtsne(as.matrix(x), theta=theta, max_iter=max_iter, perplexity=perplex[i], check_duplicates = F)
      #Perform Rtsne function
      rownames(res[[i]]$Y)<-rownames(x)
      #Add the samplename names

      #Header and result, cat sets the first col
      cat("\n#H:perplexity_",perplex[i],"\t", file=outfile, append=TRUE)
      write.table(res[[i]]$Y, file=outfile, col.names = F, quote=F, sep="\t", append=T)

      } # end i all peplexes

   saveRDS(res, file=gsub(".txt",".rds",outfile))
   write.table(data.frame(FeatureID=rownames(x)), file=gsub(".txt","_FeatureIDs.txt",outfile), row.names = F, quote=F)

   ## Plot standard taxonomy
   ## ----------------------
   if(doPlotTaxonomy){
      g_tsne <- dlfoo2::tsnePlot(res)
      pdf(gsub(".txt","TaxAllPerplex.pdf",outfile), width=18, height = 12, useDingbats = F)
         print(g_tsne)
      dev.off()
      } # end if plot
} # end me tSNE







#' produce a data frame from a tSNE results list
#' @description Version 1.0 20190213
#' @family tsne
#' @param tSNEres results list (tSNEres) from run.tsne-analysis
#' @return data frame
#' @export
#'
tsne2df <- function(tSNEres){
   x <- tSNEres
   data.frame(sample_id=rownames(x$Y), dim1=x$Y[,1], dim2=x$Y[,2], perplexity=x$perplexity, theta=x$theta, row.names=rownames(x$Y))
}



#' produce a data frame from an umap results object
#' @description Version 1.0 20190916
#' @family umap
#' @param x results object from umap analysis
#' @return data frame
#' @export
#'
umap2df <- function(x){
      stopifnot(class(x)=="umap")
      ndim <- ncol(x$layout)
      dim_mat <- data.frame(x$layout)
      colnames(dim_mat) <- paste0("dim",1:ndim)
      res_df <- cbind(data.frame(sample_id = rownames(x$layout), stringsAsFactors = F), dim_mat)
      return(res_df)
      }



#'  produce a data frame from a PCA results object
#' @description Version 1.0 20190916
#' @family pca
#' @param x results object from pca analysis (prcomp)
#' @return data frame
#' @export
#'
pca2df <- function(x){
      stopifnot(class(x)=="prcomp")
      ndim <- ncol(x$x)
      dim_mat <- x$x
      colnames(dim_mat) <- paste0("dim",1:ndim)
      res_df <- cbind(data.frame(sample_id = rownames(x$x), stringsAsFactors = F), dim_mat)
      return(res_df)
   }


#' Wrapper for HCA
#' @description Perform Hierarchical Cluster Aanalysis using \code{\link[stats]{hclust}} function in stats package.
#' @family hca
#' @family transformation
#' @param x, a matrix to cluster samples or genes
#' @param cor.method, can be one of 'pearson', 'kendall', 'spearman', 'jacard', or 'euclidean'
#' @param linkage, default is 'ward.d'. as defined by the \code{\link[stats]{hclust}} function in stats package.
#' @param plot.dendogram, if to plot HCA tree. Default is TRUE.
#' @param sample.names, a vector that will replace sample names.
#' @param my.plot.args a list with graphical arguments.
#' @return a list with two slots: $hr: reults object of class hclust and $d, the distance matrix
#' @export
hca <- function(x, cor.method = "pearson", linkage = "ward.D", plot.dendogram = TRUE, sample.names = NULL,
   my.plot.args = list(cex = 0.5, main = NULL), new = FALSE, mar = c(4,4, 2, 1), hang = -1) {

   requireNamespace("stats")

   if (linkage == "ward") {
      linkage <- "ward.D"
   }

   if (class(x) == "ExpressionSet") {
      if (is.null(sample.names)) {
         sample.names = sampleNames(x)
      }
      x <- exprs(x)
   }

   ## Performs HCA and subsequent tree plotting It also outputs an object (list) with: $dist the originaldistance matrix $hr hclust object with e.g. $order
   y <- x
   if (cor.method %in% c("pearson", "kendall", "spearman")) {
      c <- cor(y, method = cor.method, use = "complete.obs")
      d <- as.dist(1 - c)
   }

   if (cor.method == "jaccard") {
      library(prabclus)
      c <- jaccard(y)
      d <- as.dist(c)
   }

   if (cor.method == "euclidean") {
      d <- dist(t(y))

   }

   hr <- hclust(d, method = linkage, members = NULL)
   if (!is.null(sample.names)) {
      hr$labels = sample.names
   }
   # if(plot.dendogram){plot(hr, hang=hang)}
   if (plot.dendogram) {
      my.plot.args$x <- hr
      my.plot.args$hang <- hang
      do.call("plot", args = my.plot.args)
   }
   out.list <- list()
   out.list$hclust <- hr
   out.list$dist.mat <- d
   return(out.list)
}  # end function hca



#' dendextend package utils
#' @description Classify a dendogram into n=k subclusters using \code{\link[dendextend]{cutree}}. Can be used for input to reorder branches in a dendogram
#' @family hca
#' @family transformation
#' @param dend, an object of class \code{\link[dendextend]{dendrogram}}
#' @param k, the number of subclusters to be defined
#' @return an integer vector with the same length as terminal nodes and that defines subcluster groups
#' @export
dend_classify <- function(dend, k){
   requireNamespace("dendextend")
   my.vec <- dendextend::cutree(tree = dend, k = k)
   my.vec <- my.vec[match(labels(dend), names(my.vec))]
   my.fac <- as.factor(my.vec)
   levels(my.fac)[unique(my.vec)] <- 1:length(unique(my.vec))
   my.vec.out <- as.integer(as.character(my.fac))
   names(my.vec.out) <- names(my.vec)
   return(my.vec.out)
}

#' dendextend package utils
#' @description Reorder clusters within a dendrogram. Uses the \code{\link[dendextend]{rotate}} and \code{\link[dlfoo]{dend_classify}} functions.
#' @family hca
#' @family transformation
#' @param dend, an object of class \code{\link[stats]{dendrogram}}
#' @param cluster.order, index vector that specify i) the number (k) of subclusters (length of vector), and ii) what subcluster that should come first, second ... etc
#' @return an object of class \code{\link[stats]{dendrogram}}
#' @export
dend_reorder <- function(dend, cluster.order=c(2,1)){
   requireNamespace("dendextend")
   my.k <- length(cluster.order)
   my.classes <- dend_classify(dend, k=my.k)

   names.new.order <- unlist(sapply(cluster.order, function(x) names(my.classes)[my.classes==x]))

   dend.reordered <- dend %>% rotate(names.new.order)
   return(dend.reordered)
}# end function dl.dend.reorder


#' dendextend package utils
#' @description dendextend package utils
#' @describeIn  Reverses the order of nodes within a specified cluster. Uses the \code{\link[dendextend]{rotate}}  and \code{\link[dlfoo]{dend_classify}} functions.
#' @family hca
#' @family transformation
#' @param dend, an object of class \code{\link[stats]{dendrogram}}
#' @param rev.clusters, vector of T or F values. Will define i) how many clusters (k) and ii) what subclusters that should be reversed.
#' @return an object of class \code{\link[stats]{dendrogram}}
#' @export
dend_rev_clusters <- function(dend, rev.clusters=c(T,F,F)){
   requireNamespace("dendextend")
   # k is defined from rev.clusters
   # cluster.order specify
   #     - the number (k) of subclusters (length of vector)
   #     - AND which subcluster that should be reversed
   my.k <- length(rev.clusters)
   my.classes <- dend_classify(dend, k=my.k)

   rev.i <- which(rev.clusters==T)

   my.classes.rev <- my.classes

   for(i in rev.i){
      names(my.classes.rev)[my.classes==i] <- rev(names(my.classes)[my.classes==i])
   }
   dend.reordered <- dend %>% rotate(names(my.classes.rev))
   return(dend.reordered)
}# end function dl.dend.reorder



# gs_rcc <-   readRDS(file="/Volumes/Crab2000/RESOURCES/Methylation_Infinium/datasets/gdac_rcc/gdac_rcc_867set_norm.hg38_vmfAtac_6968")
#x.me <- gs_rcc[,,"beta"] # dim: 485512 636
#   str(x.me)



#' The consensusCluster-dl function
#' @description  Wrapper for ConsensusClusterPlus analysis \code{\link[ConsensusClusterPlus]{ccRun}}
#' @family hca
#' @family analysis
#' @param d expressionSet
#' @param run.name index vector that specify i) the number (k) of subclusters (length of vector), and ii) what subcluster that should come first, second ... etc
#' @param results.dir name of directory where analysis is to be saved
#' @return consensus cluster results
#' @export


ConsensusClusterPlus.dl <- function (
   d = NULL,
   run.name = NULL,
   results.dir = NULL,
   maxK = 3,
   reps = 10, # 10,
   innerLinkage = "ward", # "average"
   finalLinkage = "ward", # "average",
   distance = "pearson",
   var.method = 'sd',
   var.cut = 1,

   pItem = 0.8,
   pFeature = 1,
   clusterAlg = "hc",
   #title =  NULL # '~/PROJECTS/RCC_TCGA_2014/ConsensuClust/' # "untitled_consensus_cluster"

   ml = NULL,
   tmyPal = NULL,
   seed = NULL,
   plot = "png",
   writeTable = FALSE,
   weightsItem = NULL,
   weightsFeature = NULL,
   verbose = F,
   corUse = "everything"
){
   require(ConsensusClusterPlus)

   stopifnot(!is.null(results.dir))

   if(class(d) == 'ExpressionSet'){
      if(!is.null(var.cut)){ d <- esetVarianceFilter(d, method = var.method, cutoff = var.cut)}
      cat('\n...... eset with',length(sampleNames(d)),'samples.',length(featureNames(d)),'genes')
      d = t(apply(exprs(d), 1, function(y) y-mean(y)))
      cat('\n...... centering each row on mean value \n')
   } # if d is an eset

   cat('\n\n... START ConsensusClustering of sammples')
   stopifnot(!is.null(results.dir))

   # Define name (title)
   my.dir <- paste(
      'ConsensusCluster_', run.name, '_',
      ncol(d),'set','_',
      nrow(d),'genes',
      '_reps.',reps,
      '_maxK.', maxK,
      '_innerLinkage.', innerLinkage,sep='')
   if(!is.null(var.cut)) my.dir <- paste0(my.dir, '_',var.method,'.',var.cut)
   title <- paste0(results.dir,'/' ,my.dir, '/')
   #if(file.exists(out.dir)) stop('Out dir already exists! Remove or use other name')
   #if(!file.exists(out.dir)) dir.create(out.dir)
   dir.create(title)
   cat('\n ... Created output directory:', title)

   # Start actual ConsensusClustering
   ## ----------------------
   if (is.null(seed) == TRUE) {seed = timeSeed = as.numeric(Sys.time())} #  is.null(seed)
   set.seed(seed)

   if (is.null(ml) == TRUE) {
      if (!class(d) %in% c("dist", "matrix", "ExpressionSet")) {stop("d must be a matrix, distance object or ExpressionSet (eset object)")}
      if (inherits(d, "dist")) {
         if (is.null(attr(d, "method"))) {
            attr(d, "method") <- distance <- "unknown - user-specified"
         }
         if (is.null(distance) || (distance != attr(d, "method"))) {
            distance <- attr(d, "method")
         }
         if ((!is.null(pFeature)) && (pFeature < 1)) {
            message("Cannot use the pFeatures parameter when specifying a distance matrix as the data object\n")
            pFeature <- 1
         }
         if (!is.null(weightsFeature)) {
            message("Cannot use the weightsFeature parameter when specifying a distance matrix as the data object\n")
            weightsFeature <- NULL
         }
         if (clusterAlg == "km") {
            message("Note: k-means will cluster the distance matrix you provided.  This is similar to kmdist option when suppling a data matrix")
         }
      }else{
         if (is.null(distance)) {distance <- "pearson"}
      }
      if ((clusterAlg == "km") && inherits(distance, "character") &&
          (distance != "euclidean")) {
         message("Note: The km (kmeans) option only supports a euclidean distance metric when supplying a data matrix.  If you want to cluster a distance matrix using k-means use the 'kmdist' option, or use a different algorithm such as 'hc' or 'pam'.  Changing distance to euclidean")
         distance <- "euclidean"}

      if (inherits(d, "ExpressionSet")) {d <- exprs(d)}
      # ml <- ccRun(d = d, maxK = maxK, repCount = reps, diss = inherits(d,
      #             "dist"), pItem = pItem, pFeature = pFeature, innerLinkage = innerLinkage,
      #             clusterAlg = clusterAlg, weightsFeature = weightsFeature,
      #             weightsItem = weightsItem, distance = distance, verbose = verbose,
      #             corUse = corUse)

      ## Perform ccRun function
      ## -----------------------
      ml <- ConsensusClusterPlus:::ccRun(d = d, maxK = maxK, repCount = reps, diss = inherits(d,
                                                                                              "dist"), pItem = pItem, pFeature = pFeature, innerLinkage = innerLinkage,
                                         clusterAlg = clusterAlg, weightsFeature = weightsFeature,
                                         weightsItem = weightsItem, distance = distance, verbose = verbose,
                                         corUse = corUse)

   } # end is.null(ml) == TRUE (end first section)

   # Start process results
   # --------------------
   res = list()
   if ((is.null(plot) == FALSE | writeTable) & !file.exists(paste(title,
                                                                  sep = ""))) {
      dir.create(paste(title, sep = ""))} # end create dir

   log <- matrix(ncol = 2, byrow = T, c("title", title, "maxK",
                                        maxK, "input matrix rows", ifelse(inherits(d, "matrix"),
                                                                          nrow(d), "dist-mat"), "input matrix columns", ifelse(inherits(d,
                                                                                                                                        "matrix"), ncol(d), ncol(as.matrix(d))), "number of bootstraps",
                                        reps, "item subsampling proportion", pItem, "feature subsampling proportion",
                                        ifelse(is.null(pFeature), 1, pFeature), "cluster algorithm",
                                        clusterAlg, "inner linkage type", innerLinkage, "final linkage type",
                                        finalLinkage, "correlation method", distance, "plot",
                                        if (is.null(plot)) NA else plot, "seed", if (is.null(seed)) NA else seed))
   colnames(log) = c("argument", "value")


   if (writeTable) {
      write.csv(file = paste(title, "/", title, ".log.csv",
                             sep = ""), log, row.names = F)}


   ## Plotting
   ## -----------------
   if (is.null(plot)) {
   } else if (plot == "pngBMP") {
      bitmap(paste(title, "/", "consensus%03d.png", sep = ""))
   } else if (plot == "png") {
      png(paste(title, "/", "consensus%03d.png", sep = ""))
   } else if (plot == "pdf") {
      pdf(onefile = TRUE, paste(title, "/", "consensus.pdf",
                                sep = ""))
   } else if (plot == "ps") {
      postscript(onefile = TRUE, paste(title, "/", "consensus.ps",
                                       sep = ""))
   }


   # Plot consensus matrix legend
   # -----------------------------
   colorList = list()
   colorM = rbind()
   thisPal <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C",
                "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6",
                "#6A3D9A", "#FFFF99", "#B15928", "#bd18ea", "#2ef4ca",
                "#f4cced", "#f4cc03", "#05188a", "#e5a25a", "#06f106",
                "#85848f", "#000000", "#076f25", "#93cd7f", "#4d0776",
                "#ffffff")
   colBreaks = NA
   if (is.null(tmyPal) == TRUE) {
      colBreaks = 10
      tmyPal <- ConsensusClusterPlus:::myPal(colBreaks)
   }else {
      colBreaks = length(tmyPal)
   }
   sc = cbind(seq(0, 1, by = 1/(colBreaks)))
   rownames(sc) = sc[, 1]
   sc = cbind(sc, sc)
   heatmap(sc, Colv = NA, Rowv = NA, symm = FALSE, scale = "none",
           col = tmyPal, na.rm = TRUE, labRow = rownames(sc), labCol = F,
           main = "consensus matrix legend")


   # Plot Consensus Matrixes
   # --------------------
   my.maxK <- 2:maxK

   for (tk in my.maxK) {
      if (verbose) {
         message(paste("consensus ", tk))}

      fm = ml[[tk]]
      hc = hclust(as.dist(1 - fm), method = finalLinkage)
      message("clustered")
      ct = cutree(hc, tk)
      names(ct) = colnames(d)
      if (class(d) == "dist") {
         names(ct) = colnames(as.matrix(d))}
      c = fm
      colorList <- ConsensusClusterPlus:::setClusterColors(res[[tk - 1]][[3]], ct,
                                                           thisPal, colorList)
      pc = c
      pc = pc[hc$order, ]
      if (!is.null(plot) && plot == "pngBMP") {
         pc = pc[, hc$order]
         pc = rbind(pc, 0)
         oc = colorList[[1]][hc$order]
         heatmap(pc, Colv = NA, Rowv = NA, symm = FALSE, scale = "none",
                 col = tmyPal, na.rm = TRUE, labRow = F, labCol = F,
                 mar = c(5, 5), main = paste("consensus matrix k=",
                                             tk, sep = ""), ColSideCol = oc)
      } else {
         pc = rbind(pc, 0)
         heatmap(pc, Colv = as.dendrogram(hc), Rowv = NA,
                 symm = FALSE, scale = "none", col = tmyPal, na.rm = TRUE,
                 labRow = F, labCol = F, mar = c(5, 5), main = paste("consensus matrix k=",
                                                                     tk, sep = ""), ColSideCol = colorList[[1]])
      }

      legend("topright", legend = unique(ct), fill = unique(colorList[[1]]),
             horiz = FALSE)
      res[[tk]] = list(consensusMatrix = c, consensusTree = hc,
                       consensusClass = ct, ml = ml[[tk]], clrs = colorList)
      colorM = rbind(colorM, colorList[[1]])
   } # end for (tk in 2:maxK)
   res[[1]] = colorM


   # Plot CDF
   # -----------
   ConsensusClusterPlus:::CDF(ml)
   ConsensusClusterPlus:::clusterTrackingPlot(colorM[, res[[length(res)]]$consensusTree$order])
   if (is.null(plot) == FALSE) {
      dev.off()}


   # If writeTable
   ## ---------------
   if (writeTable) {
      for (i in 2:length(res)) {
         write.csv(file = paste(title, "/", title, ".k=",
                                i, ".consensusMatrix.csv", sep = ""), res[[i]]$consensusMatrix)
         write.table(file = paste(title, "/", title, ".k=",
                                  i, ".consensusClass.csv", sep = ""), res[[i]]$consensusClass,
                     col.names = F, sep = ",")}}

   # Return Results
   ## ---------------
   cc.results <- res
   save(cc.results, file = paste0(title, 'cc.results.Rdata'))
   return(res)
}



#' wrapper that process resuts from ConsensusClusterPlus analyses
#' @description Function used to process output data from \code{\link[ConsensusClusterPlus]{ConsensusClusterPlus}} function.
#' @family hca
#' @family analysis
#' @param cc.dir, directory created by ConsensusClusterPlus analysis containing results. File must be called "cc.results.Rdata".
#' @param k, define the number of consensus cluster groups
#' @return a list with 5 slots:
#' \itemize{
#'   \item cc: The absolute expression values over all genes present for the signature is averaged. Negatively correlated genes are weighted by (-1)
#'   \item i: Genes in the signature are treated as an expression vector based on their regulation (a series of +1 and -1, respectively). For each sample, the correlation between its expression values and the signature vector will represent the pathway score.
#'   \item consensus.tree: The consensus clutering tree
#'   \item mat: the consensus clust matrix representing co-clutering
#'   \item sample.names: sample names
#' }
#' @export

ConsensusCluster_process <- function(
   cc.dir = NULL,
   k = 2
){
   # get cc.results
   requireNamespace("ConsensusClusterPlus")

   cc = get(load( paste0(cc.dir,'cc.results.Rdata')))

   cc.k <- list()
   cc.k$mat <- cc[[k]]$consensusMatrix
   cc.k$sample.names <- names(cc[[k]][["consensusClass"]])
   dimnames(cc.k$mat) <- list(cc.k$sample.names, cc.k$sample.names)
   cc.k$cc <- as.character(cc[[k]][["consensusClass"]])
   cc.k$i = cc[[k]]$consensusTree$order

   cc.k$consensus.tree<-cc[[k]][["consensusTree"]]
   cc.k$consensus.tree$labels = names(cc[[k]][["consensusClass"]])

   return(cc.k)
} # end plot consensus clusters






   # source("~/PROJECTS/BRCA_SR_KP/R_SCRIPTS/R_CAF_Soruce.R")
   #    my_dir <- "~/PROJECTS/BRCA_SR_KP/DESeq2_analysis/EnrichmentAnalysis/"
   #
   #
   #    rank_list <- lapply(results_list, function(x){
   #          my_tab <- data.frame(x@listData)
   #          my_tab$gene_id <- x@rownames
   #          my_tab <- dplyr::left_join(my_tab, fdata, by="gene_id")
   #          head(my_tab)
   #          my_tab$stat[is.na(my_tab$stat)] <- 0
   #          my_tab <- my_tab[!duplicated(my_tab$hgnc_symbol),]
   #          y <- my_tab$stat
   #          names(y) <- my_tab$hgnc_symbol
   #          return(y)
   #       })
   #    str(rank_list)
   # gene.list <- lapply(rank_list, names)
   ## run.folder = "~/PROJECTS/BRCA_SR_KP/DESeq2_analysis/EnrichmentAnalysis/ClusterProfiler/", run.name="ClusterProfiler_caf_e2_all_combos"

#gene.list <-lapply(readRDS( file="/Volumes/Crab2000/PROJECTS/RCC_ccpRCC_2019//Limma/rcc_nt/Limma_taxSimp2_fit2_top250_list.rds")[[1]], function(x) x$SYMBOL)
#run.folder <- "/Volumes/Crab2000/PROJECTS/RCC_ccpRCC_2019//Limma/rcc_nt/EA_ClusterProfiler"
#run.name <- "test"

#' enrichmentAnalysis_ClusterProfiler v2
#' @description wrapper that perform Enrichment analysis using ClusterProfiler and DOSE packages
#' @description v2 that includes v7 of Msig 20190923
#' @family enrichment
#' @family gene ontology
#' @family go
#' @param gene.list, directory created by ConsensusClusterPlus analysis containing results. File must be called "cc.results.Rdata".
#' @param run.folder, name of higher level directory
#' @param run.name, name of run (directory created)
#' @param my.ontologies what ontologies to run
#' @param entrez.or.symbol if to run on suymbol or entrezgene
#' @param min.genes.in.term, minimum included genes in sigificant term
#' @return write results to utput directory: EnricherResults.txt: tab delimited text file with all significant results from all terms and all gene lists. EnricherResults.rds: same as .txt but saved as .rds. EnricherResults_comparecluster.rds. ResultsSummary.txt: Results Summary on how many significant terms per gene list were obtained.
#' @export
enrichmentAnalysis_ClusterProfiler <- function(
   gene.list = NULL,
   run.folder = NULL,
   run.name = NULL,
   entrez.or.symbol = "symbol",
   my_ontologies = c("h.all","c2.cp.biocarta","c2.cp.kegg","c2.cp.reactome","c3.tft","c5.bp","c5.cc","c5.mf"),
   #term.names =  c("c2","c3.tft","c6","c7","c5.cc","c5.bp","c5.mf"),
   min.genes.in.term = 3
   #plot.results = T
   ){

   if(!is.list(gene.list)) stop("gene.list must be list vith character vectors of gene symbols")

   # myDir <- file.path(runFolder, runName)
   myDir <- run.folder
   if(!dir.exists(myDir)) dir.create(myDir)
   myDir <- file.path(myDir,run.name)
   if(dir.exists(myDir)) stop("run already exists - delete old run folder")
   dir.create(myDir)


   # Run/Output objects
      # enrich.list = vector(mode="list", length=length(term.names))
      # names(enrich.list) = term.names
      # enrich.cmp = enrich.list
   require(clusterProfiler)
   require(DOSE)

   # Load signatures (mSig etc)
   # message(" ... loading mSig gene signatures...")
   # if(!exists("ch")) ch <- read.gmt("/Volumes/Crab2000/RESOURCES/MSigDB/msigdb_v6.0_GMTs/h.all.v6.0.symbols.gmt")
   # if(!exists("c1")) c1 <- read.gmt("/Volumes/Crab2000/RESOURCES/MSigDB/msigdb_v6.0_GMTs/c1.all.v6.0.symbols.gmt")
   # if(!exists("c2.bc")) c2.bc <- read.gmt("/Volumes/Crab2000/RESOURCES/MSigDB/msigdb_v6.0_GMTs/c2.cp.biocarta.v6.0.symbols.gmt")
   # if(!exists("c2.kegg")) c2.kegg <- read.gmt("/Volumes/Crab2000/RESOURCES/MSigDB/msigdb_v6.0_GMTs/c2.cp.kegg.v6.0.symbols.gmt")
   # if(!exists("c2.cgp")) c2.cgp <- read.gmt("/Volumes/Crab2000/RESOURCES/MSigDB/msigdb_v6.0_GMTs/c2.cgp.v6.0.symbols.gmt")
   # if(!exists("c2.rea") )c2.rea <- read.gmt("/Volumes/Crab2000/RESOURCES/MSigDB/msigdb_v6.0_GMTs/c2.cp.reactome.v6.0.symbols.gmt")
   # if(!exists("c3.tft")) c3.tft <- read.gmt("/Volumes/Crab2000/RESOURCES/MSigDB/msigdb_v6.0_GMTs/c3.tft.v6.0.symbols.gmt")
   # if(!exists("c6")) c6 <- read.gmt("/Volumes/Crab2000/RESOURCES/MSigDB/msigdb_v6.0_GMTs/c6.all.v6.0.symbols.gmt")
   # if(!exists("c7")) c7 <- read.gmt("/Volumes/Crab2000/RESOURCES/MSigDB/msigdb_v6.0_GMTs/c7.all.v6.0.symbols.gmt")
   # if(!exists("c5.cc")) c5.cc <- read.gmt("/Volumes/Crab2000/RESOURCES/MSigDB/msigdb_v6.0_GMTs/c5.cc.v6.0.symbols.gmt")
   # if(!exists("c5.bp")) c5.bp <- read.gmt("/Volumes/Crab2000/RESOURCES/MSigDB/msigdb_v6.0_GMTs/c5.bp.v6.0.symbols.gmt")
   # if(!exists("c5.mf")) c5.mf <- read.gmt("/Volumes/Crab2000/RESOURCES/MSigDB/msigdb_v6.0_GMTs/c5.mf.v6.0.symbols.gmt")
   # c2 <- do.call("rbind", list(c2.bc, c2.kegg, c2.rea, ch))

   message("... start")
   # pathway_files <- list.files("/Volumes/Crab2000/RESOURCES/MSigDB/msigdb_v6.0_GMTs/", pattern="entrez_pathwayList.rds", full.names = T)
   #if(entrez.or.symbol == "entrez") pathway_files <- list.files("/Volumes/Crab2000/RESOURCES/MSigDB/msigdb_v7.0_GMTs/", pattern="entrez_pathwayList.rds", full.names = T)
   #if(entrez.or.symbol == "symbol") pathway_files <- list.files("/Volumes/Crab2000/RESOURCES/MSigDB/msigdb_v7.0_GMTs/", pattern="symbols_pathwayList.rds", full.names = T)
   if(entrez.or.symbol == "entrez") pathway_files <- list.files("/Volumes/Crab2000/RESOURCES/MSigDB/msigdb_v7.0_GMTs/", pattern="entrez.gmt", full.names = T)
   if(entrez.or.symbol == "symbol") pathway_files <- list.files("/Volumes/Crab2000/RESOURCES/MSigDB/msigdb_v7.0_GMTs/", pattern="symbols.gmt", full.names = T)


   pathway_files <- pathway_files[apply(sapply(my_ontologies, function(x) grepl(x, pathway_files)),1, any)]

   term.names <-  gsub("[.]v7.*","",gsub(".*[//]","",pathway_files))
   enrich.list = vector(mode="list", length=length(term.names))
   names(enrich.list) = term.names
   names(pathway_files) = term.names
   enrich.cmp = enrich.list

   ## internal functions
   enrich2compare <- function(enrichResult.list){
         # Create a compareClusterResult object from a collection of enrichResult objects
         u <- unlist(lapply(enrichResult.list, function(x) is.null(x)))
         enrichResult.list <- enrichResult.list[!u]
         my.results <- lapply(enrichResult.list, function(x) x@result)
         my.tab <- do.call("rbind", my.results)
         my.tab <- cbind( data.frame(Cluster=c(rep(names(my.results), times=unlist(lapply(my.results, nrow)))), my.tab) )
         if(nrow(my.tab)==0) return(c())
         rownames(my.tab) <- c(1:nrow(my.tab))
         my.clusters <- lapply(enrichResult.list, function(x) x@gene)
         my.compare <- new("compareClusterResult", compareClusterResult = my.tab, geneClusters = my.clusters)
         return(my.compare)
         }
   foo.trim.res <- function(res){
         u <- nchar(res@compareClusterResult$Description)>40
         if(!all(u)) return(res)
         res@compareClusterResult$Description[u] <- paste(substr(res@compareClusterResult$Description[u],1,40),"...")
         res@compareClusterResult <- res@compareClusterResult[res@compareClusterResult$Count>2,]
         rownames(res@compareClusterResult) <- 1:nrow(res@compareClusterResult)
         return(res)
         }


   ## Loop all msig terms
   for(i in 1:length(term.names)){
      term.name <- term.names[i]
      message("... ... performing EA using mSig term: ", term.name)
      #my.term <- eval(parse(text=term.name))
      #my.term <- readRDS(pathway_files[term.name])
      my.term <- read.gmt(pathway_files[term.name])

      # if(!any(unlist(lapply(gene.list, function(x){sapply(x, function(y) any(unlist(y %in% my.term$gene)))})))) message("no genes from list present in ", term.name)

      enrich.list[[term.name]] <- lapply(gene.list, function(x) enricher(x, TERM2GENE=my.term, pvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 5, qvalueCutoff = 0.2))
      enrich.cmp[[term.name]] <- enrich2compare(enrich.list[[term.name]])
      my_res <- foo.trim.res(enrich.cmp[[term.name]])
      #    if(plot.results){
      #       pdf(file=paste0(myDir,"/",run.name,"_", term.name, ".pdf"), width=8.267, height=11.692, useDingbats = F)
      #          plot(enrich.cmp[[2]], type="dot", title=term.name, showCategory=50)
      #       dev.off()}
      #
      # }
   }


   message("... ... DONE")
   ####
   my.tab <- do.call("rbind", lapply(enrich.cmp, function(x) x@compareClusterResult))
   my.tab <- my.tab[-which(my.tab$Count < min.genes.in.term), ]
   my.tab$ontology <- gsub("[.].*","",rownames(my.tab))
   str(my.tab)
   summaryResults <- table(my.tab$Cluster, my.tab$ontology)

   message("... Saving results to file")
   write.table(my.tab, file=paste0(myDir,"/",run.name,"_EnricherResults.txt"), row.names=F, quote=F, sep="\t")
   saveRDS(enrich.list, file=paste0(myDir,"/",run.name,"_EnricherResults.rds"))
   saveRDS(enrich.cmp, file=paste0(myDir,"/",run.name,"_EnricherResults_comparecluster.rds"))
   write.table(summaryResults, file=paste0(myDir,"/",run.name,"_Summary_SiginficantTerms.txt"), row.names=T, quote=F, sep="\t")

} # end function erichment analysis



#' create a rank-list as input for GSEA
#' @aliases gsea.rnk.object
#' @family enrichment
#' @family gene ontology
#' @family msig
#' @family gsea
#' @param gene.identifiers, list of named numeric vectors. One for each entity (sample/group etc). Each numeric vector is named with gene identifiers. no duplicated (duplicated names allowed)
#' @param rank.values, name of higher level directory
#' @param invert.rank, name of run (directory created)
#' @param collapse.method, what mSig or other signatues to load
#' @param convert.input.to.rank fgsea param min size
#' @param use.absolute.for.triming fgsea param max size
#' @export
gsea_RankList <- function(
         gene.identifiers, rank.values, invert.rank=FALSE,
         collapse.method="highest rank", convert.input.to.rank=FALSE,
         use.absolute.for.triming=TRUE){

      if(length(gene.identifiers)!=length(rank.values)) stop()
      u <- gene.identifiers %in% c("",NA)
      if(length(which(u)==T)>0) message("removing ", length(which(u)==T), " blank or NA values from string of gene identifiers" )
      gene.identifiers <- gene.identifiers[!u]
      rank.values <- rank.values[!u]

   	## collapse.methods
   	  #	"highest rank"  (the highest ranked identifier is chosen)
   	  #	"median"
   	    #	collapse.median.method=	"rank" or "input value" !! Not working!!!

   	#	rank values: 	vector with values to interpret as rank; for GSEA, the highest value is top scored (or a pvalue)
   	#	invert.rank (TRUE/FALSE);

   	#	transform.pval.to.rank; pval transformed with 1-log(p) to "invert" rank
   	#	use.absolute.for.triming; the absolute rank-value (eg p-value or log-diff) is used to identify the highest scoring replicate
   	gene.rank<-rank.values
   	u<-order(gene.rank, decreasing=T)
   	gene.rank<-gene.rank[u]
   	gene.identifiers<-gene.identifiers[u]
   	rm(u)

     # Invert rank
   	if(invert.rank){
   	    gene.rank<-gene.rank*(-1)
   	    u<-order(gene.rank, decreasing=T)
   	    gene.rank<-gene.rank[u]
   	    gene.identifiers<-gene.identifiers[u]
   	    } # end if invert rank

   	 if(convert.input.to.rank){
   	    gene.rank<-c(length(gene.rank):1)
   	   } # end if covert to rank

     # rank.output
   	 rank.output<-gene.rank

     # if use absolute for trimming
   	if(use.absolute.for.triming){
   	    gene.rank<-abs(gene.rank)
   	    u<-order(gene.rank, decreasing=T)
   	    gene.rank<-gene.rank[u]
   	    gene.identifiers<-gene.identifiers[u]
   	    rank.output<-rank.output[u]
   	    rm(u)	}


   	## find unique identifiers (remove null and NA)
   	  unique.ids<-unique(gene.identifiers)
       gene.rank<-gene.rank[!is.na(gene.identifiers)]
       rank.output<-rank.output[!is.na(gene.identifiers)]
       gene.identifiers<-gene.identifiers[!is.na(gene.identifiers)]

   	  gene.rank<-gene.rank[!unique.ids==""]
   	  rank.output<-rank.output[!unique.ids==""]
   	  gene.identifiers<-gene.identifiers[!unique.ids==""]


     ## If duplicated reporters ->  collapse on highest - or - median values
       if( collapse.method=="highest rank"){
   	      u<-duplicated(gene.identifiers)
   	      id.out<-gene.identifiers[!u]
   	      rank.output<-rank.output[!u]
           } # end collapse highest rank

   	  if( collapse.method=="median"){
   	    unique.ids<-unique(gene.identifiers)
         id.out<-c()
   	    rank.out<-c()
         for(i in 1:length(unique.ids)){
   	      gene.i<-unique.ids[i]
   	      ii<-which(gene.identifiers==gene.i)
   	      id.out<-c(id.out, gene.i)
   	      rank.out<-c(rank.out, median(rank.output[ii]))
   	       }
         rank.output<-rank.out
   	      } # End collapse median

     ## reorder data to rank.output
   	  u<-order(rank.output, decreasing=T)
   	  rank.output<-rank.output[u]
   	  id.out<-id.out[u]


     # create output object
   	  my.out<-data.frame(id.out, stringsAsFactors=FALSE)
   	  my.out$rank.out<-rank.output

   	return(my.out)
	} ## end function .rnk object



# my_dir <- "~/PROJECTS/BRCA_SR_KP/DESeq2_analysis"
#
#    my_files <- list.files(my_dir, pattern = "_resultsTable_")
#    rank_list <- lapply(my_files, function(x){
#          my_tab <- read.delim(file.path(my_dir,x), as.is=T)
#          my_tab$stat[is.na(my_tab$stat)] <- 0
#          my_tab <- my_tab[!duplicated(my_tab$hgnc_symbol),]
#          y <- my_tab$stat
#          names(y) <- my_tab$hgnc_symbol
#          return(y)
#       })
#    str(rank_list)
#
#    x <- gsub("DeSeq2_resultsTable_","",my_files)
#    x <- gsub(".txt","",x)
#    x <- gsub("MCF7_CAF_E2_","exp",x)
#    names(rank_list) <- x
# rank.list <- rank_list
# entrez.or.symbol="symbol", run.folder = "~/PROJECTS/BRCA_SR_KP/DESeq2_analysis/EnrichmentAnalysis/fGSEA/", run.name="fgsea_caf_e2_all_combos"


# limma_allDeg_list <- readRDS(file="/Volumes/Crab2000/PROJECTS/RCC_ccpRCC_2019//Limma/rcc_nt/Limma_taxSimp2_fit2_allDEG_list.rds")
# my_dir <- "/Volumes/Crab2000/PROJECTS/RCC_ccpRCC_2019//Limma/rcc_nt/EA_fgsea_ClusterProfiler_allDeg/"
# my_fit <- readRDS(file="/Volumes/Crab2000/PROJECTS/RCC_ccpRCC_2019//Limma/rcc_nt/Limma_taxSimp2_fit2.rds")
#
# for(my_contrast in colnames(my_fit)){
#       message("\n....... ", my_contrast)
#
#       t_vec <- my_fit$t[,my_contrast]
#       names(t_vec) <- my_fit$genes$SYMBOL
#       t_list <- list()
#       t_list[[my_contrast]] <- t_vec
#
#       # dir.create(my_dir)
#       enrichmentAnalysis_fgsea(rank.list = t_list,
#          run.folder = my_dir, run.name = my_contrast)
#
#    }#}
# my_contrast <- colnames(my_fit)[1]

# rank.list = t_list, pathway.list.gmt=kidney_pathways
   # limma_allDeg_list <- readRDS(file="/Volumes/Crab2000/PROJECTS/RCC_ccpRCC_2019//Limma/rcc_nt/Limma_taxSimp2_fit2_allDEG_list.rds")
   # my_dir <- "/Volumes/Crab2000/PROJECTS/RCC_ccpRCC_2019//Limma/rcc_nt//EA_fgsea_KidneySignatures/"
   # my_fit <- readRDS(file="/Volumes/Crab2000/PROJECTS/RCC_ccpRCC_2019//Limma/rcc_nt/Limma_taxSimp2_fit2.rds")
   # kidney_pathways <- readRDS("/Volumes/Crab2000/RESOURCES/MSig_DL/KIDNEY_COLLECTION/msig_kidney_list_gmt.rds")
   #
   # my_contrast <- colnames(my_fit)[1]
   # message("\n....... ", my_contrast)
   #       t_vec <- my_fit$t[,my_contrast]
   #       names(t_vec) <- my_fit$genes$SYMBOL
   #       t_list <- list()
   #       t_list[[my_contrast]] <- t_vec
   # run.folder = my_dir, run.name = my_contrast)

#' perform Enrichment analysis fgsea using rank values
#' @description https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/06_Gene_set_testing.nb.html
#' @description https://bioconductor.org/packages/release/bioc/manuals/fgsea/man/fgsea.pdf
#' @description updated 20190916. Saveing complete list as well as significant ontologies
#' @family enrichment
#' @family gene ontology
#' @family msig
#' @family gsea
#' @family fgsea
#' @param rank.list, list of named numeric vectors. One for each entity (sample/group etc). Each numeric vector is named with gene identifiers. no duplicated (duplicated names allowed)
#' @param run.folder, name of higher level directory
#' @param run.name, name of run (directory created)
#' @param entrez.or.symbol what identifiers to use
#' @param my.ontologies, what mSig or other signatues to load
#' @param minSize fgsea param min size
#' @param maxSize fgsea param max size
#' @param nperm fgsea param nperm
#' @return write results to utput directory
#' \itemize{
#'   \item _fgsea_allOntologies: tab delimited text file with all significant ontologies.
#' }
#' @export
enrichmentAnalysis_fgsea <- function(
   rank.list = NULL,
   run.folder=NULL,
   run.name=NULL,
   pathway.list.gmt = NULL,
   my_ontologies = c("h.all","c2.cp.biocarta","c2.cp.kegg","c2.cp.reactome","c3.tft","c5.bp","c5.cc","c5.mf"),
   entrez.or.symbol = c("symbol"),
   minSize=10, maxSize = 500, nperm=1e6, minLeadingEdge=3, adjp.filter=0.05, leadingEdgeMin=4
   ){
   require(tibble)
   require(tidyr)
   require(dplyr)
   require(fgsea)
   my_ranks <- rank.list

   myDir <- file.path(run.folder, run.name)
   #myDir <- run.folder
   if(!dir.exists(myDir)){ dir.create(myDir)}
   # myDir <- file.path(myDir,run.name)
   # if(dir.exists(myDir)) stop("run already exists - delete old run folder")


   # Run fgsea for all selected tissues and selected ontology subsets
   # Clean upp resutls to tab format and keep padj<0.05

   message("... start")
   if(is.null(pathway.list.gmt)){
      message("... no list of patwhays provided, loading MSig GMTs")
      # pathway_files <- list.files("/Volumes/Crab2000/RESOURCES/MSigDB/msigdb_v6.0_GMTs/", pattern="entrez_pathwayList.rds", full.names = T)
      if(entrez.or.symbol == "entrez") pathway_files <- list.files("/Volumes/Crab2000/RESOURCES/MSigDB/msigdb_v7.0_GMTs/", pattern="entrez_pathwayList.rds", full.names = T)
      if(entrez.or.symbol == "symbol") pathway_files <- list.files("/Volumes/Crab2000/RESOURCES/MSigDB/msigdb_v7.0_GMTs/", pattern="symbols_pathwayList.rds", full.names = T)
      pathway_files <- pathway_files[apply(sapply(my_ontologies, function(x) grepl(x, pathway_files)),1, any)]

      fgsea_res_list <- lapply(pathway_files, function(x){
         my_pathway <- readRDS(file = x)
         message("... running fgsea pathway ", x)
         lapply(my_ranks, function(y) {
            my_f <- fgsea(y, pathways = my_pathway, minSize=minSize, maxSize = maxSize, nperm=nperm)
            return(my_f)
            #return(my_f[padj<0.05,])
            })
         })
      names(fgsea_res_list) <- gsub("[.]v7.*","",gsub(".*[//]","",pathway_files))

   }
   if(!is.null(pathway.list.gmt)){
      message("... list of pathways entries provided using this")

      fgsea_res_list <- list(lapply(my_ranks, function(y) {
            my_f <- fgsea(y, pathways = pathway.list.gmt, minSize=minSize, maxSize = maxSize, nperm=nperm)
            return(my_f)
            #return(my_f[padj<0.05,])
            }))
      names(fgsea_res_list) <- run.name
      }

   # str(fgsea_res_list)
      # fgsea_res_list <- readRDS("~/OneDrive - Lund University/PROJECTS_od/RCC_VHPT/DESeq2/EnrichmentAnalysis/fGSEA/fgsea_all_combos_nperm1e6/fgsea_all_combos_nperm1e6_fgsea_allOntologies.rds")


   ## process full lists
   ## ::::::::::::::::
      apa <- fgsea_res_list
      apa <- lapply(apa, function(x){
         y <- unlist(lapply(x, nrow))
         x[!(y==0)]
         })
     # str(apa)

      apa <- lapply(apa, function(x){
         y <- names(x)
         lapply(y, function(y, z=x){
            zz <-z[[y]]
            zz$tissue <- y
            return(zz)
         })
         })
    #  str(apa)

      apa2 <- lapply(names(apa), function(x){
         y <- apa[[x]]
         # yy <- do.call("rbind", y)
         yy <- dplyr::bind_rows(y)
         yy$ontology <- x
         return(yy)
         })
      #str(apa2)
   u <- lapply(apa2, function(x) length(x))
   apa2 <- apa2[(unlist(u)>1)]
   apa2 <- as_tibble(dplyr::bind_rows(apa2))
   #str(apa2)
   #apa2[,1]
   leading_edge <- unlist(lapply(apa2[,"leadingEdge"]$leadingEdge, function(x) paste0(x, collapse="|")))
   leading_edge_n <- unlist(lapply(apa2[,"leadingEdge"]$leadingEdge, length))
   apa2 <- apa2 %>% mutate(leading_edge = leading_edge) %>% mutate(leading_edge_n = leading_edge_n)

   #str(leading_edge)
   u<-which(colnames(apa2)=="leadingEdge")
   #str(apa2[,-u])

   write.table(apa2[,-u], file=paste0(myDir,"/",run.name,"_fgsea_Ontologies_all_table.txt"), sep="\t", quote=F, row.names=F)
   saveRDS(fgsea_res_list, file=paste0(myDir,"/",run.name,"_fgsea_allOntologies.rds"))

   #Hmisc::describe(apa2$padj)
   uu <- apa2$padj < adjp.filter & apa2$leading_edge_n>=leadingEdgeMin
   table(uu)
   write.table(apa2[uu,-u], file=paste0(myDir,"/",run.name,"_fgsea_Ontologies_significant.txt"), sep="\t", quote=F, row.names=F)

   }





