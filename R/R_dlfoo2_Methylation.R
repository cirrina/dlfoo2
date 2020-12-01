
#' meMergeGenoSetsBeta
#' @description Version 1.0 20190209
#' @family merge
#' @family methylation
#' @family genoset
#' @family genomics
#' @param gs_list list of objects of class genoset
#' @return save rds to file
#' @export
#'
meMergeGenoSetsBeta <- function(gs_list){
      require(genoset)
      require(dplyr)
      message("\n  merging GenoSets!!!")

      stopifnot(all(unlist(lapply(gs_list, class))=="GenoSet"))
      if(!all(unlist(lapply(gs_list, function(x) identical(rownames(x), rownames(gs_list[[1]])))))) stop("not all rownames (featurenames) identical")
      if(!all(unlist(lapply(gs_list, function(x) identical(genome(x), genome(gs_list[[1]])))))) stop("not all genomes are identical")

      mat_list <- lapply(gs_list, function(x){x[,,"beta"]})
      pdata_list <- lapply(gs_list, function(x){colData(x)})
      fdata_list <- lapply(gs_list, function(x){elementMetadata(x)})
      sample_list <- lapply(gs_list, function(x) colnames(x))

      my_mat <- do.call("cbind", mat_list)
      colnames(my_mat) <- unlist(lapply(gs_list, colnames))

      my_list <- lapply(pdata_list, as.data.frame)
      pdata <- dplyr::bind_rows(my_list)
      pdata$sample_id <- unlist(lapply(gs_list, colnames))
      rownames(pdata) <- pdata$sample_id
      pdata <- dplyr::select(pdata, c(sample_id, everything()))

      pdata <- DataFrame(pdata)


      # fdata <- rowRanges(gs_list[[1]])

      fgr <- rowRanges(gs_list[[1]])
      #fgr <- GRanges(seqnames=fdata$seqnames, ranges=IRanges(start=fdata$start, end=fdata$end), FeatureID=fdata$FeatureID)
      #names(fgr) <- rownames(fdata$FeatureID)
      #genome(fgr) <- "hg19"
      stopifnot(identical(colnames(my_mat), rownames(pdata)))
      stopifnot(identical(names(rowRanges(gs_list[[1]])), rownames(my_mat)))

      gs <- GenoSet(rowRanges = fgr, assays=list(beta=as.matrix(my_mat)), colData = pdata)
      print(gs)
      message("\n  DONE!!!")
      return(gs)
   } # end meMergeGenoSets





#' Normalize Illumina Methylation Ringner/Staaf
#' @description function meNormalizeBeta
#' @description  Version 1.0 20190108
#' @family methylation
#' @family genoset
#' @family genomics
#' @family transformation
#' @param x object of class GenoSet
#' @param runName Fie dir-Name to use forthis run, e.g. rcc_me450.
#' @param runFolder directory where processed genoset object should be saved to
#' @param returnObject if to return object instead of write to disk
#' @param projectLabel character vecto. name of project
#' @param probe.annotations.file file name for probe annotations that inlcude 'INFINIUM_DESIGN_TYPE' column: type I and II ILMN probe sets
#' @return genoset
#' @export
meNormalizeBeta <- function(x, runName, runFolder, returnObject=F,
   projectLabel=NULL, probe.annotations.file='~/RESOURCES/Methylation_Infinium/Annotations/Probe.manifest_masked_intersect_450.epic_hg38.rds'){
   require(genoset)
   if(class(x)!="GenoSet") stop("gset is not of class 'GenoSet'")
   if(is.null(colnames(x))) stop("no colnames, sample names, set for genoset. please fix!! ")

   # myDir <- file.path(runFolder, runName)
   myDir <- file.path(runFolder)
   # outFile <- paste0(myDir, "/", runName, ".rds")
   outFile <- paste0(myDir, "/", runName, ".rds")
   if(file.exists(outFile) & !returnObject) stop("Outfile already exists. Delete or change run name")
   #if(dir.exists(myDir) & !returnObject) stop("Run already exists. Set new name or delete/move old run")
   if(!is.null(myDir)) dir.create(myDir)

      find.max<-function(x,span=1,method="strict"){
      # x = density vector or series vector with values
      # span = number of items on each side of the tested point. 1 = window of 3, i-1,i,i+1
      #function supposes that the data is correctly sorted,
      # and that the data is relatively smooth
      maxima<-c()
      for(t in 1:length(x)){
         span.low<-t-span
         span.high<-t+span
         if(method=="strict"){
            if( (span.low<0) | (span.high>length(x))){
               #not valid for a check as the window spans outside the vector indexes
            }else{
               if(which.max(x[span.low:span.high])==(span+1)) maxima<-c(maxima,t) #this is a local maxima
            }}}
      maxima.order<-order(x[maxima],decreasing=TRUE)
      return(maxima[maxima.order])
   } # end funtion find.max

   ## select peaks function
   select.peaks<-function(m){
      peaks<-c(2,-1)
      for(t in 1:nrow(m)){
         if((m$x[t]<0.4) & (m$x[t]<peaks[1]) & (m$y[t]>(0.05*m$y[1]))) {peaks[1]=m$x[t]}
         if((m$x[t]>0.6) & (m$x[t]>peaks[2]) & (m$y[t]>(0.05*m$y[1]))) {peaks[2]=m$x[t]}
      } # inf for t
      if(peaks[1]==2) {peaks[1]=0}
      if(peaks[2]==-1) {peaks[2]=1}
      return(peaks)
   } # end function select peaks





   # Log
   logFile <- file.path(myDir, "log.txt")
   cat(file = logFile, "PARAMS\n")
   cat(file = logFile, append = T, "runName:", runName, "\n")
   cat(file = logFile, append = T, "runFolder:", runFolder,"\n")
   cat(file = logFile, append = T, "Processd data directory:", myDir,"\n")
   cat(file = logFile, append = T, "Probe info file:", probe.annotations.file,"\n")

   cat(file = logFile, append = T, "\nDATA\n")
   cat(file = logFile, append = T, "dim input genoSet: ", dim(x),"\n")
   # if(length(metadata(x))>0){
   #    lapply(names(metadata(x)), function(y){
   #       cat(file = logFile, append = T, y,": ",metadata(x)[[y]],"\n")
   #    })
   # }
   cat(file = logFile, append=T, "\n------\n")

   #


   message("... reading probe annotations. matching genoset to probes in file")
   #if(grepl("rdata",probe.annotations.file, ignore.case = T)) probe.annotations <- get(load(probe.annotations.file))
   if(grepl("rds",probe.annotations.file, ignore.case = T)) probe.annotations <- readRDS(probe.annotations.file)
   probe.annotations <- as.data.frame(elementMetadata(probe.annotations))
   # rownames(probe.annotations)

   #
   # u <- rowSums(x[,,"beta"]==0)
   # which(u==1)


   my_probes <- intersect(rownames(x), rownames(probe.annotations))
   cat(file = logFile, append = T, "n features intersect with probe info file: ", length(my_probes),"\n")
   x <- x[my_probes]
   probe.annotations <- probe.annotations[my_probes,]
      #str(probe.annotations)

   message("... Defining Type I/II probes")
   probes_I = rownames(probe.annotations[(which(probe.annotations$designType=="I")),])
   probes_II = rownames(probe.annotations[(which(probe.annotations$designType=="II")),])

   data = x[,,"beta"]
   dim(data)
   #dimnames(data)


      # head(data[1000:2000,])
   data_I=data[which(!is.na(match(rownames(data),probes_I))),]
   data_II=data[which(!is.na(match(rownames(data),probes_II))),]

   peaks <- data.frame(unmet_I=rep(NA,ncol(data)),met_I=rep(NA,ncol(data)),unmet_II=rep(NA,ncol(data)),met_II=rep(NA,ncol(data)))
   rownames(peaks)=colnames(data)

   plotFile = file.path(myDir, 'peaks_preNorm.pdf')
   pdf(plotFile)
   message("... plotting peaks pre normalization")
   for(sample in 1:ncol(data)) {
      cat("loop", sample, "\n")
      d_I=density(data_I[,sample],kernel="epanechnikov",bw=0.02,from=0,to=1,na.rm=T)
      d_II=density(data_II[,sample],kernel="epanechnikov",bw=0.02,from=0,to=1,na.rm=T)
      #d_Ir=density(data_Ir[,sample],kernel="epanechnikov",bw=0.02,from=0,to=1,na.rm=T)
      #d_Ig=density(data_Ig[,sample],kernel="epanechnikov",bw=0.02,from=0,to=1,na.rm=T)

      plot(d_I,col="blue",main="")
      lines(d_II,col="orange")
      # lines(d_Ir,col="red")
      #lines(d_Ig,col="green")
      title(colnames(data)[sample])
      max_I<-find.max(d_I$y,span=3);
      max_I=cbind(max_I,d_I$x[max_I])
      max_I=cbind(max_I,d_I$y[max_I])
      colnames(max_I)=c("index","x","y")
      max_II<-find.max(d_II$y,span=3);
      max_II=cbind(max_II,d_II$x[max_II])
      max_II=cbind(max_II,d_II$y[max_II])
      colnames(max_II)=c("index","x","y")
      peaks_I=select.peaks(as.data.frame(max_I))
      peaks_II=select.peaks(as.data.frame(max_II))
      abline(v=peaks_I,col="blue")
      abline(v=peaks_II,col="orange")
      peaks[colnames(data)[sample],]=append(peaks_I,peaks_II)
   } # end sample loop
   dev.off()


   ## Normalize!
   #####################
   message("... normalizing")
   data.norm=as.matrix(data)
   rows_I=which(!is.na(match(rownames(data.norm),probes_I)))
   rows_II=which(!is.na(match(rownames(data.norm),probes_II)))
   tmp=as.matrix(data) # str(tmp)
   for(sample in 1:ncol(data.norm)) {
      data.norm[rows_I,sample]=(tmp[rows_I,sample]-peaks[sample,1])/(peaks[sample,2]-peaks[sample,1])
      data.norm[rows_II,sample]=(tmp[rows_II,sample]-peaks[sample,3])/(peaks[sample,4]-peaks[sample,3])
      } # head(data.norm[1:10,], 20)

   rm(tmp)
   gc()

   data.norm=as.data.frame(data.norm)
   data.norm[]=lapply(data.norm, function(x) ifelse (x<0,0,x))
   data.norm[]=lapply(data.norm, function(x) ifelse (x>1,1,x))

   data.norm_I=data.norm[which(!is.na(match(rownames(data.norm),probes_I))),]
   data.norm_II=data.norm[which(!is.na(match(rownames(data.norm),probes_II))),] # str(data.norm_II)

   plotFile = file.path(myDir, 'peaks_postNorm.pdf')
   #file.name = paste(eset.dir, eset.name,'_peaks2.pdf',sep='')
   pdf(plotFile)
   op<-par(mfrow=c(2,1))
   for(sample in 1:ncol(data)) {

      d_I=density(data_I[,sample],kernel="epanechnikov",bw=0.02,from=0,to=1,na.rm=T)
      d_II=density(data_II[,sample],kernel="epanechnikov",bw=0.02,from=0,to=1,na.rm=T)
      d.norm_I=density(data.norm_I[,sample],kernel="epanechnikov",bw=0.02,from=0,to=1,na.rm=T)
      d.norm_II=density(data.norm_II[,sample],kernel="epanechnikov",bw=0.02,from=0,to=1,na.rm=T)

      plot(d_I,col="blue",main="")
      lines(d_II,col="orange")
      abline(v=peaks[sample,1:2],col="blue")
      abline(v=peaks[sample,3:4],col="orange")
      title(colnames(data.norm)[sample])

      plot(d.norm_I,col="blue",main="")
      lines(d.norm_II,col="orange")

      }
   par(op)
   dev.off()

   stopifnot(all(dim(data) == dim(data.norm)))
   message("... creating normalized genoset")
   gs <- x # str(x)
   gs[,,"beta"] <- as.matrix(data.norm) # str(as.matrix(data.norm))
      # str(gs[,,"beta"])
      # str(data.norm)

   message("... adding project metadata")
   if(!is.null(projectLabel)) metadata(gs)[["projectLabel"]] <- projectLabel
   metadata(gs)[["meNormalizeBeta"]] <- timestamp()
   cat(file = logFile, append=T, "OUTPUT\n")
   cat(file = logFile, append=T, "\n------\n")
   cat(file = logFile, append = T, "dim normalized genoSet:", dim(gs) ,"\n")
   # lapply(names(metadata(gs)), function(y){
   #       cat(file = logFile, append = T, y,": ",metadata(x)[[y]],"\n")
   #    })
   print(gs)
   message("... saving normalized genoSet to file")

   if(!returnObject){
      saveRDS(gs, outFile)
      cat(file = logFile, append = T, "normalized genoSet saved to:", myDir ,"\n")
      message("... DONE!")
      }
   if(returnObject) return(gs)

} # end function meNormalizeBeta











#' function meFilterBetaNa
#' @description Version 1.0 20190108
#' @family genoset
#' @family methylation
#' @family genomics
#' @family filter
#' @param x object of class GenoSet with assay data slot for 'beta' values
#' @param na.filter numeric. ratio of maximum allowed na-values
#' @return genoset
#' @export
#'
meFilterBetaNa <- function(x, na.filter=0.5){
   require(genoset)
   if(class(x)!="GenoSet") stop("gset is not of class 'GenoSet'")

   u <- apply(x[,,"beta"], 1, function(x) length(which(is.na(x)))/length(x) )
   uu <- u<=na.filter
   print(table(uu))

   x <- x[uu]
   metadata(x)[["na.filter"]] <- paste0("na.filter, max NA-vals: ", na.filter*100,"%; ",
      length(which(uu==F))," features removed; ",
      length(which(!uu==F))," features kept; ",timestamp())
   message(metadata(x)[["na.filter"]])
   print(x)
   return(x)
} # end function meFilterBetaNa




#


#' VMF - Varying Methylated Features
#' @description FIlter GenoSet for Varying Methylated features
#' @family genoset
#' @family genomics
#' @family filter
#' @family methylation
#' @param x Genoset with methylation beta values as assay data slot
#' @param vmfQantile What quantiles to use for VMF variance filter.
#' @param vmfBetaCut What beta cutoff to use for quantile VMF varition cutoff
#' @param naFilter MAX ratio NA values that are allowed per feature
#' @param projectLabel character vector; name of project
#' @return genoSet object with fitered features (VMF)
#'@examples
#'\dontrun{
#' gs_rcc <- readRDS(file="/Volumes/MacPro2TB/RESOURCES/Methylation_Infinium/ProcessedData/gdac_rcc/gdac_rcc_867set_norm.hg38.rds")
#' x <- gs_rcc[,sample(1:ncol(gs_rcc), size=100, replace=F)]
#'}
#' @export
#'
meFilterVMF <- function(x,
   #runName="Temp", runFolder="/Volumes/MacPro2TB/RESOURCES/Methylation_450/Analyses/VMP",
   projectLabel = NULL,
   vmfQuantile = c(0.05,0.95), vmfBetaCut=0.1,
   naFilter = 0
   ){
   require(reshape2)
   require(dplyr)
   options(dplyr.width = Inf)
   require(tibble)
   require(tidyverse)
   require(tidyr)
   require(GenomicRanges)
   require(GenomicFeatures)


   require(genoset)
   if(class(x)!="GenoSet") stop("gset is not of class 'GenoSet'")

   ## filter NA
   if(!is.null(naFilter)){
      message("... performing optional NA filter")
      x <- dlfoo2::meFilterBetaNa(x, na.filter = naFilter)
      print(metadata(x)[["na.filter"]])
   }

   ## set features to genome order
   message("... matching up genome order")
   x <- toGenomeOrder(x)

   ## Filter to remove non-varying features
   ## --------------
   message("... filtering Variable Methylated Features (VMF) based on percentiles & beta-cut")


   memat <- x[,,"beta"]
   myFlag <- apply(memat, 1, function(x){
      y<-quantile(x, probs=vmfQuantile, na.rm=T)
      abs(y[1]-y[2])>=vmfBetaCut
   })
   # table(myFlag)


   memat <- memat[myFlag,]
   # message("... Filter Variable Methylated Features (VMF). Kept: ", length(which(myFlag==T)), " features. Removed: ", length(which(myFlag==F)))

   ## create filtered genoset and return
   ## ----
   message("... returning genoset object")
   fdata <- rowRanges(x[myFlag,])
   stopifnot(identical(names(fdata), rownames(memat)))

   gs <- GenoSet(rowRanges = fdata, assays=list(beta=as.matrix(memat)), colData = colData(x), metadata = metadata(x))
   if(!is.null(projectLabel)) metadata(gs)[["projectLabel"]] <- projectLabel
   metadata(gs)[["meFilterVMF"]] <- paste0("Filter Variable Methylated Features (VMF). Kept: ",
         length(which(myFlag==T)), " features. Removed: ", length(which(myFlag==F)),
         "features. Params: vmfBetaCut= ", vmfBetaCut, "; vmfQuantile= ",paste(vmfQuantile,collapse=","),
         ". ", timestamp())
   print(gs)
   message(metadata(gs)[["meFilterVMF"]])
   message("... DONE!")
   return(gs)
} # end function MeFilterVMF




#' meMergeCpGs - merge adjacent CpGs
#' @family methylation
#' @family granges
#' @family genoset
#' @family genomics
#' @family transformation
#' @family filter
#' @param x Genoset with methylation beta values as assay data slot
#' @param mergeCpGwindow window size for CpG (mean) merging to decrease redundancy
#' @param projectLabel character vector; name of project
#' @param naFilter MAX ratio NA values that are allowed per feature
#' @return genoset
#' @export
#'
meMergeCpGs <- function(x,
   projectLabel = NULL, mergeCpGwindow=50, naFilter = NULL
   ){
   require(GenomicRanges)
   require(genoset)
   if(class(x)!="GenoSet") stop("gset is not of class 'GenoSet'")

## filter NA
   if(!is.null(naFilter)){
      message("... performing additional NA filter")
      x <- meFilterBetaNa(x, na.filter = 0)
      print(metadata(x)[["na.filter"]])
      }

   ## Match up with probe info & order according to
   message("... matching up genome order")
   x <- toGenomeOrder(x)


   ## Merge adjacent CpGs
   ## ----------------
   message("... merging beta values for adjacent Features; window = ", mergeCpGwindow)

   meProbes <- rowRanges(x)
   meProbesWind <- GRanges(seqnames = seqnames(meProbes), ranges = IRanges(start=start(meProbes)-round(mergeCpGwindow/2,0), end = end(meProbes)+round(mergeCpGwindow/2,0)))
   stopifnot(all(width(meProbesWind)==mergeCpGwindow+1))
   meProbesReduced <- GenomicRanges::reduce(meProbesWind) ##tror det är “reduce” som kollapsar överlappande ranges till en enda range
   # meProbesReduced # GRanges object with 199088 ranges and 0 metadata columns:

   my_overlaps <- findOverlaps(meProbes, meProbesReduced)
   o <- split(my_overlaps@from, my_overlaps@to)

   memat <- x[,,"beta"]
   mematMerge <- lapply(o, function(y){
      y <- unlist(y)
      if(length(y)>1) {return(apply(memat[y,],2,mean, na.rm=T))
         }
      else{
         return(memat[y,])
      }
   })
   mematMerge <- do.call("rbind", mematMerge)
   # str(mematMerge) #  'data.frame':	26585 obs. of  20 variables:
   mematMerge <- as.matrix(mematMerge)

   myProbes <- lapply(o, function(x){
      x <- unlist(x)
      paste(meProbes$FeatureID[x], collapse = "_")
   })

   pgr <- meProbesReduced
   #pgr$FeatureID <- unlist(myProbes)
   pgr <- GRanges(seqnames = seqnames(meProbesReduced), ranges = IRanges(start=mid(meProbesReduced), end=mid(meProbesReduced)), FeatureID=unlist(myProbes))
   names(pgr) <- pgr$FeatureID
   rownames(mematMerge) <- as.character(unlist(myProbes))
   colnames(mematMerge)
   colnames(x)

   gs <- GenoSet(rowRanges = pgr, assays=list(beta=mematMerge), colData = colData(x), metadata = metadata(x))
   if(!is.null(projectLabel)) metadata(gs)[["projectLabel"]] <- projectLabel
   metadata(gs)[["meMergeCpGs"]] <- paste0("Merge beta for adjacent CpGs. From: ",
         nrow(x), " features. To: ", nrow(gs),
         " features. Params: mergeCpG window= ", mergeCpGwindow,
         ".  ", timestamp())
   print(gs)
   message(metadata(gs)[["meMergeCpGs"]])
   message("... DONE!")
   return(gs)
   } # end function merge CpGs


#' Create logfile with genoset info
#' @family genoset
#' @family info
#' @param file.name full file path of .rds genoset object
#' @export
#'
gsSinkGenosetDescription <- function(file.name){
   require(genoset)
   x <- readRDS(file.name)
   myFile <- gsub(".rds","_gsDescription.txt",file.name)

   sink(file=myFile)
      print(x)
   sink(file=NULL)
   cat(file = myFile, append=T, "\n\n--- METADATA ---\n")

   sink(file=myFile, append=T)
      lapply(names(x@metadata), function(y){
         if(length((x@metadata[[y]])) == 1){
            print(y)
            print("------------")
            print(x@metadata[[y]])
            }
         if(length((x@metadata[[y]])) > 1){
            print(y)
            print("------------")
            head(x@metadata[[y]])
            }
         })
   sink(file=NULL)
   return()
   } # end function merge CpGs



#' pearson correlation between CpG beta and GEX
#' @description Version 1.0 20190110
#' @description Function:  Correlation Me beta vs GEX
#' @family methylation
#' @family granges
#' @family genoset
#' @family genomics
#' @family analysis
#' @family filter
#' @param x.me object of class GenoSet. Contains assayData slot 'beta' with beta values. Sample names (colnames) must match sampleNames for x.gex
#' @param x.gex object of class expressionSet. Contains gene expression data. sampleNames for x.gex must match sample names (colnames) for x.me
#' @param runName Fie/dir-Name to use forthis run, e.g. rcc_me450.
#' @param runFolder Higher level directory where processed genoset object should be saved to
#' @param gexSDcut Variance filter cutoff for GEX data. Defaults to 0.5.
#' @param corWindow window surrounding CpG to include genes for correlation analysis
#' @return saving table of correlations to file
#' @export
#'
meAnalysisCorrGex <- function(x.me,x.gex,
   runName, runFolder,
   gexSDcut=0.25,
   corWindow=1e6
   ){
   require(GenomicRanges)
   require(genoset)
   require(rtracklayer)
   require(biomaRt)

   require(reshape2)
   require(dplyr)
   options(dplyr.width = Inf)
   require(tibble)
   require(tidyverse)
   require(tidyr)


   if(class(x.me)!="GenoSet") stop("gset is not of class 'GenoSet'")
   if(class(x.gex)!="ExpressionSet") stop("gset is not of class 'ExpressionSet'")
   if(!identical(sampleNames(x.gex), colnames(x.me))) stop("sample names in genoset and expression set are not identical")

   #source("~/PROJECTS/RCC_ccpRCC//R_SCRIPTS/R_ccpRCC_SOURCE.R")
   #setwd("~/PROJECTS/RCC_ccpRCC/global_me450/")
   #fdata_me <- as.data.frame(readRDS("Data/probes_450k_ccpStudy.rds"))
   # pdata_me
   myDir <- file.path(runFolder, runName)
   if(dir.exists(myDir)) stop("Run already exists. Set alternative 'runName' or delete/move old run from disk")
   dir.create(myDir)

   #  Log
   # ------
   logFile <- file.path(myDir, "AnalysisLog.txt")
   cat(file = logFile, "PARAMS\n-------\n")
   cat(file = logFile, append = T, "runName:", runName, "\n")
   cat(file = logFile, append = T, "runFolder:", runFolder,"\n")
   cat(file = logFile, append = T, "Processd data directory:", myDir,"\n")
   cat(file = logFile, append = T, "corWindow: ", corWindow,"\n")
   cat(file = logFile, append = T, "gexSDcut: ", gexSDcut,"\n")
   cat(file = logFile, append = T, "\n\n")
   cat(file = logFile, append=T, "\n------\n")
   sink(file=logFile, append=T)
      cat("Genoset:\n")
      print(x.me)
      cat("\n")
      print(x.me@metadata)
      cat("\nExpressionSet:\n")
      print(x.gex)
   sink(file=NULL)
   cat(file = logFile, append=T, "\n------\n")


   ## GEX data - variance filter
   ## -------------------
   message("... variance filter of GEX data. SD-cutoff = ", gexSDcut)
   cat(file = logFile, append = T, "\n\n------\nVariance filter of GEX data. SD-cutoff = ", gexSDcut)
   u <- apply(exprs(x.gex), 1, function(x) sd(x)>=gexSDcut)
   #table(u)
   message("... ... n features in GEX eset pre variance filter: ", dim(x.gex))
   x.gex <- x.gex[u,]
   message("... ... n features in GEX eset post variance filter: ", dim(x.gex))
   cat(file = logFile, append = T, "\nSD filter gex removed: ", length(which(u==F)), "SD filter gex kept: ", length(which(u==T)), "\n")
   print(x.gex)
   sink(file=logFile, append=T)
      print(x.gex)
   sink(file=NULL)


   ## Genomic positions for TCGA gex data
   ## ----------------------------
   message("... acquire genomic positions for gex features using biomart host 'grch37.ensembl.org'")

   cat(file = logFile, append = T, "\n\n------\nGenomic positions, genes GEX data.")
   cat(file = logFile, append = T, "\n fething chr locations for ensembl_gene_id using biomart host 'grch37.ensembl.org")

   message("... ... 1: use 'ensembl_gene_id' to obtain chr locations")
   fdata.gex <- fData(x.gex)

   mart_human = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org",
                      path="/biomart/martservice",dataset="hsapiens_gene_ensembl")
   # listAttributes(mart_human)
   my_attributes <- c("ensembl_gene_id" ,"start_position", "end_position", "chromosome_name","strand", "hgnc_symbol", "entrezgene_id")
   if(!("ENSG" %in% colnames(fdata.gex))) stop("ensembl_gene_id column in gex fData must be named 'ENSG'")
   probe_df <- getBM(attributes=my_attributes, filters="ensembl_gene_id",  values=fdata.gex$ENSG, mart=mart_human, uniqueRows = F)
   #str(probe_df) #
   #as_tibble(probe_df) %>% filter(hgnc_symbol=="HNF1A")
   message("... ... ... grch37 (hg19) chr locations obtained for ", length(unique(probe_df$ensembl_gene_id)), " ensembl_gene_ids")
   cat(file = logFile, append = T, "\n... grch37 (hg19) chr locations obtained for ", length(unique(probe_df$ensembl_gene_id)), " ensembl_gene_ids")

   ## for genes that not are found - e.g. HNF1B - redo using ENTREZ as value
   message("... ... 2: for 'ensembl_gene_id' whose chr locations not available in biomart, try obtain using entrezgene id (genes which ensemble ids have been repalced in grch38)")
   cat(file = logFile, append = T, "\n some 'ensembl_gene_id' have been repalced in grch38")
   cat(file = logFile, append = T, "\n fetching chr locations for these using entrezgene id's")


   if(!("ENTREZID" %in% colnames(fdata.gex))) stop("ensembl_gene_id column in gex fData must be named 'ENTREZID'")
   my_entrez <- as_tibble(as.data.frame(fdata.gex)) %>% dplyr::filter(!(ENSG %in% probe_df$ensembl_gene_id)) %>% dplyr::pull(ENTREZID)
   probe_df2 <- getBM(attributes=my_attributes, filters="entrezgene_id",  values=my_entrez, mart=mart_human, uniqueRows = F)
   str(probe_df2) #'data.frame':	221 obs. of  7 variables:

   ## get the corresponding hg38 ensemble id for these entrez-derived probes
   my_ensgs <- fdata.gex$ENSG[match(probe_df2$entrezgene_id, fdata.gex$ENTREZID)]
   #str(my_ensgs) # chr [1:221]
   probe_df2$ensembl_gene_id <- my_ensgs
   probe_df2 <- probe_df2[which(probe_df2$chromosome_name %in% c(1:22,"X")),]
   probe_df <- rbind(probe_df, probe_df2)
   #table(duplicated(probe_df$ensembl_gene_id))
   probe_df <- probe_df[!duplicated(probe_df$ensembl_gene_id),]
   my_probes <- intersect(fdata.gex$ENSG, probe_df$ensembl_gene_id)
   str(my_probes)
   fdata.gex.hg19 <- as_tibble(probe_df) %>%
      dplyr::rename(ENSG=ensembl_gene_id) %>%
      dplyr::filter(ENSG %in% my_probes) %>%
      left_join(as_tibble(fdata.gex) %>% dplyr::select(ENSG, SYMBOL, ENTREZID), by="ENSG")

   #fdata.gex.hg19 <- dplyr::left_join(probe_df, fdata.gex, by.x="ensembl_gene_id", by.y="ENSG")
   # str(my_probes) # chr [1:15504]
   #fdata.gex <- fdata.gex[fdata.gex$ENSG %in% my_probes,]
   # fdata.gex[fdata.gex$SYMBOL=="HNF1B",]
   #fdata.gex <- fdata.gex[my_probes,] #
   #probe_df <- probe_df[probe_df$ensembl_gene_id %in% my_probes,]
   #str(probe_df)
   #str(fdata.gex)

   # gex_gr <- GRanges(seqnames = paste0("chr", fdata.gex$chr), strand=fdata.gex$strand,
   #       ranges=IRanges(start=fdata.gex$start, end=fdata.gex$end, names = fdata.gex$ENSG),
   #       ENSG=fdata.gex$ENSG, SYMBOL=fdata.gex$SYMBOL, ENTREZID=fdata.gex$ENTREZID)
   gex_gr <- GRanges(seqnames = paste0("chr", fdata.gex.hg19$chromosome_name), strand=fdata.gex.hg19$strand,
          ranges=IRanges(start=fdata.gex.hg19$start_position, end=fdata.gex.hg19$end_position, names = fdata.gex.hg19$ENSG),
          ENSG=fdata.gex.hg19$ENSG, SYMBOL=fdata.gex.hg19$SYMBOL, ENTREZID=fdata.gex.hg19$ENTREZID)


   genome(gex_gr) <- "hg19"
   gex_gr <- toGenomeOrder(gex_gr)
   message("... ... ... grch37 (hg19) chr locations obtained for a total of ", length(unique(probe_df$ensembl_gene_id)), " genes")
   cat(file = logFile, append = T, "\n... grch37 (hg19) chr locations obtained for a total of ", length(unique(probe_df$ensembl_gene_id)), " genes")

   x.gex <- x.gex[gex_gr$ENSG,]
   stopifnot(identical(featureNames(x.gex), names(gex_gr)))
   myFile <- paste0(myDir,"/",runName,"_gex_genoset_hg19.rds")
   message("... ... ... saving matched gex eset (hg19) to: ", myFile)
   cat(file = logFile, append = T, "\n... saving matched gex eset (hg19) to: ", myFile)
   gs.gex <- GenoSet(rowRanges = gex_gr, assays=list(fpkm=exprs(x.gex)), colData = pData(x.gex))
   if(!is.null(metadata(x.me)[["projectLabel"]])) metadata(gs.gex)[["projectLabel"]] <- metadata(x.me)[["projectLabel"]]
   metadata(gs.gex)[["info"]] <- paste0("Matched gene expression data for gex-me450 correlations. Gene locations (hg19) otained using biomart 'grch37.ensembl.org' ", timestamp())
   sink(file=logFile, append=T)
      print(gs.gex)
   sink(file=NULL)
   print(gs.gex)
   saveRDS(gs.gex, file=myFile)

   rm(fdata.gex, probe_df, probe_df2, u, myFile, my_ensgs, my_entrez, my_probes)


   ## START CORRELATIONS
   ## --------------
   message("Begin correlations - me vs. gex")
   cat(file = logFile, append = T, "\n\n------\nCalculating pearson correlations between CpGs and genes within given window size")

   memat <- as.matrix(x.me[,,"beta"])
   p_gr <- rowRanges(x.me)
   gexmat <- gs.gex[,,"fpkm"]
   gex_gr <- rowRanges(gs.gex)

   stopifnot(identical(colnames(memat), colnames(gexmat)))

   ## Correlations - VMPs to genes
   # Window size overlap with TSS
   cor_list <- list()
   message("... looping all ",length(p_gr)," methylation FeatureIDs. Calculating correlation to all genes within window size of ", corWindow, " bp")
   for(i in 1:length(p_gr)){
      if(i %in% seq(1,length(p_gr)+50,50)) message("...",i," of ", length(p_gr))
      gr_i <- promoters(p_gr[i], upstream = round(corWindow/2,0), downstream = round(corWindow/2,0))
      my_gex <- gex_gr[findOverlaps(gr_i,  promoters(gex_gr, upstream=0, downstream=1))@to]

      if(length(my_gex)<1) next(i)
      if(length(my_gex)>1){
         my_mat <- gexmat[my_gex$ENSG,]
         #my_cors <- apply(my_mat, 1, function(x) cor(x, memat[gr_i$FeatureID,]))
         my_cors <- apply(my_mat, 1, function(x) unlist(cor.test(x, memat[gr_i$FeatureID,])[c("estimate","p.value")]))
         my_cors <- data.frame(FeatureID=gr_i$FeatureID, ENSG=my_gex$ENSG, cor_r=my_cors[1,], cor_p=my_cors[2,], stringsAsFactors = F, row.names = 1:ncol(my_cors))
         cor_list[[gr_i$FeatureID]] <- my_cors
         next(i)
      }
      if(length(my_gex)==1){
         #my_cors <- cor(gexmat[my_gex$ENSG,], memat[gr_i$FeatureID,])
         my_cors <- unlist(cor.test(gexmat[my_gex$ENSG,], memat[gr_i$FeatureID,])[c("estimate","p.value")])
         #my_cors <- data.frame(FeatureID=gr_i$FeatureID, ENSG=my_gex$ENSG, cor=my_cors, stringsAsFactors = F, row.names = 1:length(my_cors))
         my_cors <- data.frame(FeatureID=gr_i$FeatureID, ENSG=my_gex$ENSG, cor_r=my_cors[1], cor_p=my_cors[2], stringsAsFactors = F, row.names = 1)
         #my_cors$distancetoTSS <- ChIPpeakAnno::annotatePeakInBatch(p_gr[i], AnnotationData = my_gex)$distancetoFeature
         cor_list[[gr_i$FeatureID]] <- my_cors
         next(i)
      }
   }
   # length(cor_list) - length(gex_gr)
   # [1] 15454, 50 genes without any VMP within 1Mbp
   #b <- boxplot(unlist(lapply(cor_list, nrow)))
   #str(b)
   message("... DONE! Correlations between ", length(cor_list)," CpGs and gex was performed")
   cat(file = logFile, append = T, "\n... summary of methylation CpGs and number of gene correlations:\n")
   sink(file=logFile, append=T)
      Hmisc::describe(unlist(lapply(cor_list, nrow)))
   sink(file=NULL)

   myFile <- paste0(myDir,"/",runName,"_corrList_temp.rds")
   message("... saving list of correlations")
   saveRDS(cor_list, myFile)

   cor_tab  <- do.call("rbind", cor_list)
   myFile <- paste0(myDir,"/",runName,"_corrTable.rds")
   message("... saving table of correlations. Total number of correlations: ", nrow(cor_tab))
   saveRDS(cor_tab, myFile)
   cat(file = logFile, append = T, "\n... saving table of correlations to: ", myFile)
   cat(file = logFile, append = T, "\n... Total number of me-gex correlations: ", nrow(cor_tab))
} # end function








#' randomize correlation pairs CpG and GEX
#' @description Version 1.0 20190110
#' @description  Function:  Correlation Me beta vs GEX - RANDOM SAMPLING
#' @family methylation
#' @family granges
#' @family genoset
#' @family genomics
#' @family analysis
#' @family filter
#' @param x.me object of class GenoSet. Contains assayData slot 'beta' with beta values. Sample names (colnames) must match sampleNames for x.gex
#' @param x.gex object of class expressionSet. Contains gene expression data. sampleNames for x.gex must match sample names (colnames) for x.me
#' @param runName Fie/dir-Name to use forthis run, e.g. rcc_me450.
#' @param runFolder Higher level directory where processed genoset object should be saved to
#' @param gexSDcut Variance filter cutoff for GEX data. Defaults to 0.5.
#' @param corWindow window surrounding CpG to include genes for correlation analysis
#' @param mySeed Seed number for random sampling
#' @param nRandomPairs Number of random Me-GEX pairs to be sampled
#' @return saving table of correlations to file
#' @export
#'
meAnalysisCorrGexRandom <- function(x.me,x.gex,
   runName, runFolder,
   gexSDcut=0.25,
   corWindow=1e6,
   mySeed = 100, nRandomPairs=5e6
   ){
   require(GenomicRanges)
   require(genoset)
   require(rtracklayer)
   require(biomaRt)

   require(reshape2)
   require(dplyr)
   options(dplyr.width = Inf)
   require(tibble)
   require(tidyverse)
   require(tidyr)


   if(class(x.me)!="GenoSet") stop("gset is not of class 'GenoSet'")
   if(class(x.gex)!="ExpressionSet") stop("gset is not of class 'ExpressionSet'")
   if(!identical(sampleNames(x.gex), colnames(x.me))) stop("sample names in genoset and expression set are not identical")

   #source("~/PROJECTS/RCC_ccpRCC//R_SCRIPTS/R_ccpRCC_SOURCE.R")
   #setwd("~/PROJECTS/RCC_ccpRCC/global_me450/")
   #fdata_me <- as.data.frame(readRDS("Data/probes_450k_ccpStudy.rds"))
   # pdata_me
   myDir <- file.path(runFolder, runName)
   if(dir.exists(myDir)) stop("Run already exists. Set alternative 'runName' or delete/move old run from disk")
   dir.create(myDir)

   #  Log
   # ------
   logFile <- file.path(myDir, "AnalysisLog.txt")
   cat(file = logFile, "PARAMS\n-------\n")
   cat(file = logFile, append = T, "runName:", runName, "\n")
   cat(file = logFile, append = T, "runFolder:", runFolder,"\n")
   cat(file = logFile, append = T, "Processd data directory:", myDir,"\n")
   cat(file = logFile, append = T, "corWindow: ", corWindow,"\n")
   cat(file = logFile, append = T, "gexSDcut: ", gexSDcut,"\n")
   cat(file = logFile, append = T, "mySeed: ", mySeed,"\n")
   cat(file = logFile, append = T, "nRandomPairs: ", nRandomPairs,"\n")

   cat(file = logFile, append = T, "\n\n")
   cat(file = logFile, append=T, "\n------\n")
   sink(file=logFile, append=T)
      cat("Genoset:\n")
      print(x.me)
      cat("\n")
      print(x.me@metadata)
      cat("\nExpressionSet:\n")
      print(x.gex)
   sink(file=NULL)
   cat(file = logFile, append=T, "\n------\n")


   ## GEX data - variance filter
   ## -------------------
   message("... variance filter of GEX data. SD-cutoff = ", gexSDcut)
   cat(file = logFile, append = T, "\n\n------\nVariance filter of GEX data. SD-cutoff = ", gexSDcut)
   u <- apply(exprs(x.gex), 1, function(x) sd(x)>=gexSDcut)
   #table(u)
   message("... ... n features in GEX eset pre variance filter: ", dim(x.gex))
   x.gex <- x.gex[u,]
   message("... ... n features in GEX eset post variance filter: ", dim(x.gex))
   cat(file = logFile, append = T, "\nSD filter gex removed: ", length(which(u==F)), "SD filter gex kept: ", length(which(u==T)), "\n")
   print(x.gex)
   sink(file=logFile, append=T)
      print(x.gex)
   sink(file=NULL)


   ## Genomic positions for TCGA gex data
   ## ----------------------------
   message("... acquire genomic positions for gex features using biomart host 'grch37.ensembl.org'")

   cat(file = logFile, append = T, "\n\n------\nGenomic positions, genes GEX data.")
   cat(file = logFile, append = T, "\n fething chr locations for ensembl_gene_id using biomart host 'grch37.ensembl.org")

   message("... ... 1: use 'ensembl_gene_id' to obtain chr locations")
   fdata.gex <- fData(x.gex)

   mart_human = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org",
                      path="/biomart/martservice",dataset="hsapiens_gene_ensembl")
   # listAttributes(mart_human)
   # my_attributes <- c("ensembl_gene_id" ,"start_position", "end_position", "chromosome_name","strand", "hgnc_symbol", "entrezgene")
   my_attributes <- c("ensembl_gene_id" ,"start_position", "end_position", "chromosome_name","strand", "hgnc_symbol","entrezgene_id")
   if(!("ENSG" %in% colnames(fdata.gex))) stop("ensembl_gene_id column in gex fData must be named 'ENSG'")
   probe_df <- getBM(attributes=my_attributes, filters="ensembl_gene_id",  values=fdata.gex$ENSG, mart=mart_human, uniqueRows = F)

   #str(probe_df) #
   #as_tibble(probe_df) %>% filter(hgnc_symbol=="HNF1A")
   message("... ... ... grch37 (hg19) chr locations obtained for ", length(unique(probe_df$ensembl_gene_id)), " ensembl_gene_ids")
   cat(file = logFile, append = T, "\n... grch37 (hg19) chr locations obtained for ", length(unique(probe_df$ensembl_gene_id)), " ensembl_gene_ids")

   ## for genes that not are found - e.g. HNF1B - redo using ENTREZ as value
   message("... ... 2: for 'ensembl_gene_id' whose chr locations not available in biomart, try obtain using entrezgene id (genes which ensemble ids have been repalced in grch38)")
   cat(file = logFile, append = T, "\n some 'ensembl_gene_id' have been repalced in grch38")
   cat(file = logFile, append = T, "\n fetching chr locations for these using entrezgene id's")


   if(!("ENTREZID" %in% colnames(fdata.gex))) stop("ensembl_gene_id column in gex fData must be named 'ENTREZID'")
   my_entrez <- as_tibble(as.data.frame(fdata.gex)) %>% dplyr::filter(!(ENSG %in% probe_df$ensembl_gene_id)) %>% dplyr::pull(ENTREZID)
   probe_df2 <- getBM(attributes=my_attributes, filters="entrezgene_id",  values=my_entrez, mart=mart_human, uniqueRows = F)
   str(probe_df2) #'data.frame':	221 obs. of  7 variables:

   ## get the corresponding hg38 ensemble id for these entrez-derived probes
   my_ensgs <- fdata.gex$ENSG[match(probe_df2$entrezgene_id, fdata.gex$ENTREZID)]
   #str(my_ensgs) # chr [1:221]
   probe_df2$ensembl_gene_id <- my_ensgs
   probe_df2 <- probe_df2[which(probe_df2$chromosome_name %in% c(1:22,"X")),]
   probe_df <- rbind(probe_df, probe_df2)
   #table(duplicated(probe_df$ensembl_gene_id))
   probe_df <- probe_df[!duplicated(probe_df$ensembl_gene_id),]
   my_probes <- intersect(fdata.gex$ENSG, probe_df$ensembl_gene_id)
   str(my_probes)
   fdata.gex.hg19 <- as_tibble(probe_df) %>%
      dplyr::rename(ENSG=ensembl_gene_id) %>%
      dplyr::filter(ENSG %in% my_probes) %>%
      left_join(as_tibble(fdata.gex) %>% dplyr::select(ENSG, SYMBOL, ENTREZID), by="ENSG")

   #fdata.gex.hg19 <- dplyr::left_join(probe_df, fdata.gex, by.x="ensembl_gene_id", by.y="ENSG")
   # str(my_probes) # chr [1:15504]
   #fdata.gex <- fdata.gex[fdata.gex$ENSG %in% my_probes,]
   # fdata.gex[fdata.gex$SYMBOL=="HNF1B",]
   #fdata.gex <- fdata.gex[my_probes,] #
   #probe_df <- probe_df[probe_df$ensembl_gene_id %in% my_probes,]
   #str(probe_df)
   #str(fdata.gex)

   # gex_gr <- GRanges(seqnames = paste0("chr", fdata.gex$chr), strand=fdata.gex$strand,
   #       ranges=IRanges(start=fdata.gex$start, end=fdata.gex$end, names = fdata.gex$ENSG),
   #       ENSG=fdata.gex$ENSG, SYMBOL=fdata.gex$SYMBOL, ENTREZID=fdata.gex$ENTREZID)
   gex_gr <- GRanges(seqnames = paste0("chr", fdata.gex.hg19$chromosome_name), strand=fdata.gex.hg19$strand,
          ranges=IRanges(start=fdata.gex.hg19$start_position, end=fdata.gex.hg19$end_position, names = fdata.gex.hg19$ENSG),
          ENSG=fdata.gex.hg19$ENSG, SYMBOL=fdata.gex.hg19$SYMBOL, ENTREZID=fdata.gex.hg19$ENTREZID)


   genome(gex_gr) <- "hg19"
   gex_gr <- toGenomeOrder(gex_gr)
   message("... ... ... grch37 (hg19) chr locations obtained for a total of ", length(unique(probe_df$ensembl_gene_id)), " genes")
   cat(file = logFile, append = T, "\n... grch37 (hg19) chr locations obtained for a total of ", length(unique(probe_df$ensembl_gene_id)), " genes")

   x.gex <- x.gex[gex_gr$ENSG,]
   stopifnot(identical(featureNames(x.gex), names(gex_gr)))
   myFile <- paste0(myDir,"/",runName,"_gex_genoset_hg19.rds")
   message("... ... ... saving matched gex eset (hg19) to: ", myFile)
   cat(file = logFile, append = T, "\n... saving matched gex eset (hg19) to: ", myFile)
   gs.gex <- GenoSet(rowRanges = gex_gr, assays=list(fpkm=exprs(x.gex)), colData = pData(x.gex))
   if(!is.null(metadata(x.me)[["projectLabel"]])) metadata(gs.gex)[["projectLabel"]] <- metadata(x.me)[["projectLabel"]]
   metadata(gs.gex)[["info"]] <- paste0("Matched gene expression data for gex-me450 correlations. Gene locations (hg19) otained using biomart 'grch37.ensembl.org' ", timestamp())
   sink(file=logFile, append=T)
      print(gs.gex)
   sink(file=NULL)
   print(gs.gex)
   saveRDS(gs.gex, file=myFile)

   rm(fdata.gex, probe_df, probe_df2, u, myFile, my_ensgs, my_entrez, my_probes)



   ## START CORRELATIONS
   ## --------------
   message("Begin correlations - me vs. gex")
   cat(file = logFile, append = T, "\n\n------\nCalculating pearson correlations between CpGs and genes within given window size")

   memat <- x.me[,,"beta"]
   p_gr <- rowRanges(x.me)
   gexmat <- gs.gex[,,"fpkm"]
   gex_gr <- rowRanges(gs.gex)

   stopifnot(identical(colnames(memat), colnames(gexmat)))

   #message("Generating random pairs of me-gex features, n = ", nRandomPairs)
   set.seed(mySeed)
   #randomPairTab <- sapply(1:nRandomPairs, function(x) c(sample(rownames(memat), 1), sample(rownames(gexmat), 1)))
   #str(randomPairTab)

   ## Correlations - VMPs to genes
   cor_list <- list()
   message("... looping all ",nRandomPairs," random samplings. Calculating correlation of me-ge feature pair")
   # for(i in 1:length(p_gr)){
   for(i in 1:nRandomPairs){
      if(i %in% seq(1,nRandomPairs+50,50)) message("...",i," of ", nRandomPairs)
      random_i <- sample(1:length(p_gr), 1)
      gr_i <- p_gr[random_i]
      my_gex <- gex_gr[sample(1:length(gex_gr), 1)]
      my_cors <- unlist(cor.test(gexmat[my_gex$ENSG,], memat[gr_i$FeatureID,])[c("estimate","p.value")])
      my_cors <- data.frame(FeatureID=gr_i$FeatureID, ENSG=my_gex$ENSG, cor_r=my_cors[1], cor_p=my_cors[2], stringsAsFactors = F, row.names = 1)
      cor_list[[i]] <- my_cors
      next(i)
   }


      # length(cor_list) - length(gex_gr)
   # [1] 15454, 50 genes without any VMP within 1Mbp
   #b <- boxplot(unlist(lapply(cor_list, nrow)))
   #str(b)
   message("... DONE! Correlations between ", length(cor_list)," CpGs and gex was performed")
   cat(file = logFile, append = T, "\n... summary of methylation CpGs and number of gene correlations:\n")
   sink(file=logFile, append=T)
      Hmisc::describe(unlist(lapply(cor_list, nrow)))
   sink(file=NULL)

   myFile <- paste0(myDir,"/",runName,"_corrList_temp.rds")
   message("... saving list of correlations")
   saveRDS(cor_list, myFile)

   cor_tab  <- do.call("rbind", cor_list)
   myFile <- paste0(myDir,"/",runName,"_corrTable.rds")
   message("... saving table of correlations. Total number of correlations: ", nrow(cor_tab))
   saveRDS(cor_tab, myFile)
   cat(file = logFile, append = T, "\n... saving table of correlations to: ", myFile)
   cat(file = logFile, append = T, "\n... Total number of me-gex correlations: ", nrow(cor_tab))
} # end function







#' meAnalysisCorTableStats
#' @description Version 1.0 20190115
#' @description Uses methylation genoset from a Me-GEX Correlation Analysis of VMPs \code{\link[dlfoo2]{meAnalysisCorrGex}}
#' @description Used to select signficant Expression Correlated Refions (ECRs). Additional inpt must be a run with random generated me::gex correlation pairs used to define sgnificance. from \code{\link[dlfoo2]{meAnalysisCorrGexRandom}}
#' @family methylation
#' @family granges
#' @family genoset
#' @family genomics
#' @family filter
#' @family stats
#' @param meAnalysisCorrGexFolder Path to directory where Me-GEX Correlation Analysis was performed
#' @param meAnalysisCorrGexFolderRandom Folder path to random subsampling data.
#' @param x.me methylation genoset with 'beta' values
#' @param densAdjust Adjustment value for density plots
#' @param myQuantiles Quantiles for random correlation analyses to define correlation cutoffs.
#' @param returnPlots if to retrun list of ggplot obejects
#' @return Returns tables and pdfs to file
#' @export
meAnalysisCorTableStats <- function(
   meAnalysisCorrGexFolder,
   meAnalysisCorrGexFolderRandom,
   x.me,
   densAdjust = 1/5,
   myQuantiles=c(0.001, 0.999),
   returnPlots=FALSE
   ){
   require(GenomicRanges)
   require(genoset)
   require(rtracklayer)
   require(biomaRt)
   require(gridExtra)

   require(reshape2)
   require(dplyr)
   options(dplyr.width = Inf)
   require(tibble)
   require(tidyverse)
   require(tidyr)


   myDir <- meAnalysisCorrGexFolder
   plotDir <- file.path(myDir,"Figures")
   if(!file.exists(plotDir)) dir.create(plotDir)

   #  Log
   # ------
   logFile <- file.path(myDir, "Stats_MeGexCor.txt")
   cat(file = logFile, "PARAMS\n-------\n")
   cat(file = logFile, append = T, "meAnalysisCorrGexFolder:", meAnalysisCorrGexFolder, "\n")
   cat(file = logFile, append = T, "meAnalysisCorrGexFolderRandom:", meAnalysisCorrGexFolderRandom,"\n")

   cat(file = logFile, append = T, "\n\n")
   cat(file = logFile, append=T, "\n------\n")
   sink(file=logFile, append=T)
      cat("Genoset:\n")
      print(x.me)
   sink(file=NULL)
   cat(file = logFile, append=T, "\n------\n")


   message("... reading data & corrTable ")
   x.me
   x.gex <- readRDS(file=list.files(myDir, pattern="hg19.rds", full.names = T))
   cor_tab <- readRDS(file=list.files(myDir, pattern="corrTable.rds", full.names = T))
   message("... ... nrow = ", nrow(cor_tab))

   #
   message("... calculating tss distance")
   cor_tab <- as_tibble(cor_tab) %>%
         dplyr::left_join(as_tibble(rowRanges(x.me))  %>% dplyr::mutate(pos=round((start+end)/2,0)) %>% dplyr::select(-strand, -width, -end, -start), by="FeatureID") %>%
         dplyr::left_join(as_tibble(rowRanges(x.gex)) %>% dplyr::select(-width, -seqnames), by="ENSG") %>%
         mutate(tss = if_else(strand=="+",start,end)) %>%
         mutate(tss.dist = if_else(strand=="+", pos-tss, tss-pos))


   # Density - Random subset correlations
   message("... loading random data - define cutoffs")
   cor_tab_random <- as_tibble(readRDS(file=list.files(meAnalysisCorrGexFolderRandom, pattern="corrTable.rds", full.names = T)))
   message("... ... nrow = ", nrow(cor_tab_random))

   dens.random <- density(cor_tab_random$cor_r, adjust = densAdjust)
   dens.random <- data.frame(x=dens.random$x, y=dens.random$y)
   myQ <- quantile(cor_tab_random$cor_r, probs=myQuantiles)
   dens.random$quant <- factor(findInterval(dens.random$x, myQ))
   levels(dens.random$quant) <- c("lowerQ","midQ","upperQ")

   # Density cor Analysis
   dens.cor <- density(cor_tab$cor_r, adjust = densAdjust)
   dens.cor <- data.frame(x=dens.cor$x, y=dens.cor$y)
   dens.cor$quant <- factor(findInterval(dens.cor$x, myQ))
   levels(dens.cor$quant) <- c("lowerQ","midQ","upperQ")
   my_df <- cor_tab %>%
      mutate(ECR = factor(findInterval(cor_r, myQ))) %>%
      mutate(ECR = factor(ECR, labels=c("neg","noCor","pos")))
   myQ.labels <- my_df %>% group_by(ECR) %>% summarize(n=n()) %>% mutate(pct=round(n/sum(n)*100,3))
   myQ.labels <- apply(myQ.labels, 1, function(x) paste0("bin = ",x[1],"\nn = ", x[2], "\n",x[3]," %"))


   myColors <- c("dodgerblue2","orangered2","cornsilk")
   names(myColors) <- c("neg","pos","noCor")


   message(".. ... \n... ... PLOTTING\n... ...")

   # main scatter, tss vs cor
   gg.scatter <- ggplot(my_df,
      aes(x=tss.dist, y=cor_r)) +
      theme(axis.title.y = element_text(size=8)) +
      guides(fill=FALSE)+
      ylim(c(-1,1)) +
      geom_point(aes(fill=ECR), size=0.5,  alpha=0.5, pch=21, stroke=0.1) +
      scale_fill_manual(values=myColors)



   # Density plots Random and
   gg.dens <-  ggplot(dens.cor, aes(x,y)) +
         theme(panel.background = element_rect(fill = 'gray95', colour = 'antiquewhite3', size = 0.5)) +
         theme(panel.grid.major =  element_line(colour = 'lightsteelblue2'), panel.grid.minor =  element_line(colour = 'mistyrose3')) +
         # ggtitle(label="ECRs - Expression Correlated VMPs - Fixed Threshold r>=abs(0.5)") +
         xlab("Pearson correlation") +
         ylab("density") +
         theme(axis.text.x=element_text(angle = 45, hjust = 1), axis.title.x = element_text(size=12)) +
         theme(axis.title.y = element_text(size=12)) +
         theme(axis.text = element_text(size=9), title = element_text(size=9)) +
         xlim(c(-1,1)) +
         geom_line() +
         guides(fill=FALSE) +
         geom_ribbon(aes(ymin=0, ymax=y, fill=quant), alpha=0.9) +
         scale_fill_manual(values = c("dodgerblue2","ivory","orangered2")) +
         geom_vline(xintercept = myQ, lwd=0.5, lty=3) +
         geom_vline(xintercept = 0, lwd=0.5, lty=3) +
         annotate("text", x = c(-0.75,0,0.75), y = 0.5, size=4, label = myQ.labels) +
         annotate("text", x = c(-0.75,0.75), y = 2, size=4, label = paste0("ECR cut: ",round(myQ,2))) +
         annotate("text", x = c(0), y = 2.5, size=4, label = paste0("Quantiles in random data:\n",paste(myQuantiles, collapse = ", ")))  +
         coord_flip()



   # Violin plot
   gg.vio <- ggplot(my_df)+
      aes(x=ECR, y=tss.dist, fill=ECR) +
      theme(panel.background = element_rect(fill = 'gray95', colour = 'antiquewhite3', size = 0.5)) +
      theme(panel.grid.major =  element_line(colour = 'lightsteelblue2'), panel.grid.minor =  element_line(colour = 'mistyrose3')) +
      #geom_violin(width=1) +
      #scale_fill_manual(values = color_key) +
      theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
      theme(axis.title.y = element_text(size=8)) +
      theme(axis.title.x = element_text(size=12)) +
      theme(axis.text = element_text(size=9)) +
      theme(title = element_text(size=9)) +
      theme(strip.text = element_text(size=12)) +
      theme(title = element_text(size=12)) +
      geom_violin(alpha=0.5, width=1, trim = TRUE, scale = "width", adjust = 0.2) +
      geom_boxplot(width=0.2, outlier.colour="grey50", notch = FALSE, notchwidth = .75, alpha = 0.75, colour = "black", outlier.size=0.3) +
      scale_fill_manual(values=myColors) +
      guides(fill=FALSE) +
      coord_flip()



   ## Cor over threshold within 1kb of tss
   my_df <- my_df %>%
      mutate(within1k = if_else(abs(tss.dist)<=1000, "1kb", ">1kb"))
   my_freq <- my_df %>% group_by(ECR, within1k) %>% summarize(n=n()) %>% mutate(pct=round(n/sum(n)*100,2))

   gg.pie <- ggplot(my_freq, aes('', pct, fill=within1k)) +
      ggtitle("% probes within 1kb of TSS (per ECR-group)") +
      facet_grid(~ECR) +
      geom_col(position = 'fill') +
      geom_label(aes(label = pct), position = position_fill(vjust = 0.5)) +
      coord_polar(theta = 'y') +
      scale_fill_manual(values = c("beige", "darkorange3"))


   gx <- gridExtra::arrangeGrob(grobs = list(gg.scatter, gg.dens, gg.vio, gg.pie), ncol=2, nrow=2, widths=c(8, 6), heights=c(8, 4))

   ggsave(gx, filename = file.path(plotDir, "DensityPlot_Correlations_FixedThreshold.png"), width=12, height=12)
   # Save data frame corTable
   saveRDS(my_df, file=gsub("corrTable.rds","corrTableTSS.rds",list.files(myDir, pattern="corrTable.rds", full.names = T)))

   if(returnPlots) return(list(gg.scatter, gg.dens, gg.vio, gg.pie))
} # end function
