


#' meParseGDACv2
#' @description Version 2.0 20190209 (the initial version got lost (!!)
#' @description initial parsing of methylation data downloaded from GDAC database
#' @family import
#' @family methylation
#' @family genoset
#' @family genomics
#' @param fileName file name. tab separated beta values
#' @param runName Fie dir-Name to use forthis run, e.g. rcc_me450.
#' @param runFolder Higher level directory to store processed genoset object.
#' @param projectLabel project lebel to be added through metadata(genoSet) - projectLabel
#' @return save rds to file
#' @export
#'
meParseGDACv2 <- function(fileName, runName, projectLabel=NULL,
   runFolder="/Volumes/MacPro2TB/RESOURCES/Methylation_Infinium/ProcessedData/GDAC_humanmethylation450/"){

   require(data.table)
   require(GenomicRanges)
   require(genoset)
   require(GenomeInfoDb)

   # Define and create output processed data dir
   myDir <- file.path(runFolder, runName)
   if(dir.exists(myDir)) stop("Run already exists. Set new name or delete/move old run")
   dir.create(myDir)

   # define samples and columns
   message("... begin parsing Me-file")

   my_samples <- as.character(read.table(fileName, sep="\t", nrows = 1, colClasses = "character"))
   my_samples <- unique(my_samples)[-1]

   my_columns <- as.character(read.table(fileName, sep="\t", nrows = 1, colClasses = "character", skip=1))
   #my_columns_r3 <- as.character(read.table(fileName, sep="\t", nrows = 1, colClasses = "character", skip=2))
   str(my_columns)
   #my_columns <- as.character(apply(read.table(gdacFile, sep="\t", nrows = 2, skip=0, colClasses = "character"), 2, function(x) paste(x, collapse = " ")))

   # probe_id <- read.table(fileName, sep="\t", header=F, skip=2, colClasses = c("character", rep("NULL", length(my_columns)-1)) )
   #3apa <- data.table::fread(gdacFile, drop = c(2:length(my_columns)), data.table = F)
   # Log

   logFile <- file.path(myDir, "log.txt")
   cat(file = logFile, "PARAMS\n")
   cat(file = logFile, append = T, "gdacFile:", fileName,"\n")
   cat(file = logFile, append = T, "Processd data directory:", myDir,"\n")

   cat(file = logFile, append = T, "\nDATA\n")
   cat(file = logFile, append = T, "n samples in gdac file:", length(my_samples) ,"\n")
   cat(file = logFile, append = T, "n columns in gdac file:", length(my_columns) ,"\n")

   cat(file = logFile, append=T, "\n------\n")

   # Read data matrix (beta values)
   message("... reading Beta values from file")
   u <- grep("Beta_value", my_columns)



   # beta_samples <- gsub(".AVG_Beta","",my_columns[u])
   #stopifnot(identical(beta_samples, my_samples))

   my_df <- data.table::fread(fileName, select = u,  data.table = F, skip=2, header = F, colClasses = "numeric")
   my_mat <- as.matrix(my_df)
   message("... rounding Beta to 4 decimals")
   str(my_mat)
   my_mat <- round(my_mat, 4)
   head(my_mat)
   str(my_mat)
   #str(beta_samples)
   cat(file = logFile, append = T, "dim beta matrix:", dim(my_mat) ,"\n")
   colnames(my_mat) <- my_samples


   #head(my_mat[1:5, 1:5])
   # read pribe info
   uu <- which(my_columns=="Composite Element REF")
   stopifnot(length(uu)==1)
   fdata_df <- data.table::fread(fileName, select = uu,  data.table = F, skip=1, header = T, colClasses = "character")
   cat(file = logFile, append = T, "dim feature id table:", dim(fdata_df) ,"\n")
   colnames(fdata_df)[1] <- c("FeatureID")
   rownames(fdata_df) <- fdata_df$FeatureID
   rownames(my_mat) <- fdata_df$FeatureID


   message("... Update Feature chr positions")
   message(" .... ... read FeatureID info from file")
   fdata450 <- get(load("~/RESOURCES/Methylation_Infinium/Annotations/probe.annotations.Methylation.Ringner.20031112.gRanges.Rdata"))
   fdata <- dplyr::left_join(fdata_df, as.data.frame(fdata450), by="FeatureID")
   #str(fdata2) # 'data.frame':	485577 obs. of  22 variables:
   # str(fdata) # 'data.frame':	485577 obs. of  5 variables:
   rownames(fdata) <- fdata$FeatureID
   stopifnot(identical(rownames(fdata), rownames(my_mat)))
   message("...removing NA chromosome rows")
   fdata  <- fdata[!is.na(fdata$seqnames), ]

   message("... building genoSet")
   cat(file = logFile, append = T, "\nOutput GenoSet\n")

   stopifnot(identical(my_samples, colnames(my_mat)))
   fgr <- GRanges(seqnames=fdata$seqnames, ranges=IRanges(start=fdata$start, end=fdata$end), FeatureID=fdata$FeatureID)
   names(fgr) <- rownames(fdata$FeatureID)
   genome(fgr) <- "hg19"

   gs <- GenoSet(rowRanges = fgr, assays=list(beta=my_mat[fgr$FeatureID, ]), colData = data.frame(sample_id=my_samples, stringsAsFactors = F, row.names=my_samples))
   cat(file = logFile, append = T, "dim initial genoSet:", dim(gs) ,"\n")



   message("... keep only features that annotated to chromosomes")
   gs <- gs[GenomeInfoDb::seqnames(gs)%in%paste0("chr",c(1:22,"X","Y"))]
   gs <- genoset::toGenomeOrder(gs)
   GenomeInfoDb::seqlevels(gs) <- GenomeInfoDb::seqlevelsInUse(gs)
   cat(file = logFile, append = T, "dim final genoSet, chr1-chrY:", dim(gs) ,"\n")



   message("... adding project metadata")
   metadata(gs)[["projectLabel"]] <- projectLabel
   metadata(gs)[["meParseBeadStudio"]] <- timestamp()
   metadata(gs)[["fileName"]] <- fileName
   cat(file = logFile, append = T, "\nprojectLabelt:",  metadata(gs)[["projectLabel"]] ,"\n")
   cat(file = logFile, append = T, "meParseGDAC:",  metadata(gs)[["meParseGDAC"]] ,"\n")

   print(gs)
   message("... saving genoSet to file")
   outFile <- paste0(myDir, "/", runName, ".rds")
   saveRDS(gs, outFile)
   cat(file = logFile, append = T, "final genoSet saved to:", myDir ,"\n")


   ##  save sample names to logfile
   write.table(data.frame(gdac_name=my_samples), file = gsub("log","log_gdac_sampeNames",logFile), row.names=F, sep="\t", quote=F)

   }


#' meParseInhouseBeadStudio - parse inhouse generated beadstudio to genoSet
#' @description Version 1.0 20190108
#' @description parse beadstudio to genoSet
#' @family import
#' @family methylation
#' @family genoset
#' @family genomics
#' @param fileName file name. tab separated beta values
#' @param runName Fie dir-Name to use forthis run, e.g. rcc_me450.
#' @param runFolder Higher level directory to store processed genoset object.
#' @param projectLabel project lebel to be added as metadata(genoSet)
#' @return save genoset to file
#' @export
#'
meParseInhouseBeadStudio <- function(fileName, runName, projectLabel=NULL, runFolder="/Volumes/MacPro2TB/RESOURCES/Methylation_Infinium/ProcessedData/"){
   require(data.table)
   require(GenomicRanges)
   require(genoset)

   # find uniwue GDAC data file in dir
      # gdacFile <- list.files(gdacDir, pattern = "data.txt", full.names = T)
   stopifnot(length(fileName)==1)

   # Define and create output processed data dir
   myDir <- file.path(runFolder, runName)
   if(dir.exists(myDir)) stop("Run already exists. Set new name or delete/move old run")
   dir.create(myDir)

   # define samples and columns
   message("... begin parsing beadstudio file")

   my_columns <- as.character(read.table(fileName, sep="\t", nrows = 1, colClasses = "character"))
   u <- grep("Beta", my_columns)
   my_samples <- my_columns[u]
   my_ids <- my_columns[-u]

    # Log
   logFile <- file.path(myDir, "log.txt")
   cat(file = logFile, "PARAMS\n")
   cat(file = logFile, append = T, "fileName:", fileName,"\n")
   cat(file = logFile, append = T, "Processd data directory:", myDir,"\n")

   cat(file = logFile, append = T, "\nDATA\n")
   cat(file = logFile, append = T, "n samples in gdac file:", length(my_samples),"\n")
   #cat(file = logFile, append = T, "n columns in gdac file:", length(my_columns) ,"\n")

   cat(file = logFile, append=T, "\n------\n")

   # Read data matrix (beta values)
   my_df <- data.table::fread(fileName, select = u,  data.table = F, header = T, colClasses = "numeric")
   my_mat <- as.matrix(my_df)
   message("... rounding Beta to 4 decimals")
   my_mat <- round(my_mat, 4)
   cat(file = logFile, append = T, "dim beta matrix:", dim(my_mat) ,"\n")
   colnames(my_mat) <- my_samples


   #head(my_mat[1:5, 1:5])
   # read pribe info
   fdata_df <- data.table::fread(fileName, select = which(my_columns%in%my_ids),  data.table = F, header = T, colClasses = "character")
   #head(fdata)
   cat(file = logFile, append = T, "dim feature id table:", dim(fdata_df) ,"\n")
   colnames(fdata_df)[1] <- "FeatureID"
   rownames(fdata_df) <- fdata_df$FeatureID
   rownames(my_mat) <- fdata_df$FeatureID

   message("... Update Feature chr positions")
   message(" .... ... read FeatureID info from file")
   fdata450 <- get(load("~/RESOURCES/Methylation_Infinium/Annotations/probe.annotations.Methylation.Ringner.20031112.gRanges.Rdata"))
   fdata <- dplyr::left_join(fdata_df, as.data.frame(fdata450), by="FeatureID")
   #str(fdata2) # 'data.frame':	485577 obs. of  22 variables:
   # str(fdata) # 'data.frame':	485577 obs. of  5 variables:
   rownames(fdata) <- fdata$FeatureID
   stopifnot(identical(rownames(fdata), rownames(my_mat)))
   message("...removing NA chromosome rows")
   fdata  <- fdata[!is.na(fdata$seqnames), ]

   message("... building genoSet")
   cat(file = logFile, append = T, "\nOutput GenoSet\n")

   stopifnot(identical(my_samples, colnames(my_mat)))
   fgr <- GRanges(seqnames=fdata$seqnames, ranges=IRanges(start=fdata$start, end=fdata$end), FeatureID=fdata$FeatureID)
   names(fgr) <- rownames(fdata$FeatureID)
   genome(fgr) <- "hg19"

   gs <- GenoSet(rowRanges = fgr, assays=list(beta=my_mat[fgr$FeatureID, ]), colData = data.frame(sample_id=my_samples, stringsAsFactors = F, row.names=my_samples))
   cat(file = logFile, append = T, "dim initial genoSet:", dim(gs) ,"\n")

   message("... keep only features that annotated to chromosomes")
   gs <- gs[seqnames(gs)%in%paste0("chr",c(1:22,"X","Y"))]
   gs <- toGenomeOrder(gs)
   seqlevels(gs) <- seqlevelsInUse(gs)
   cat(file = logFile, append = T, "dim final genoSet, chr1-chrY:", dim(gs) ,"\n")

   message("... adding project metadata")
   metadata(gs)[["projectLabel"]] <- projectLabel
   metadata(gs)[["meParseBeadStudio"]] <- timestamp()
   metadata(gs)[["fileName"]] <- fileName
   cat(file = logFile, append = T, "\nprojectLabelt:",  metadata(gs)[["projectLabel"]] ,"\n")
   cat(file = logFile, append = T, "meParseGDAC:",  metadata(gs)[["meParseGDAC"]] ,"\n")

   print(gs)
   message("... saving genoSet to file")
   outFile <- paste0(myDir, "/", runName, ".rds")
   saveRDS(gs, outFile)
   cat(file = logFile, append = T, "final genoSet saved to:", myDir ,"\n")
} # end meParseInhouseBeadStudio



#' meParseBeadStudio
#' @description Version 1.0 20190108
#' @family import
#' @family methylation
#' @family genoset
#' @family genomics
#' @param fileName file name. tab separated beta values
#' @param runName Fie dir-Name to use forthis run, e.g. rcc_me450.
#' @param runFolder Higher level directory to store processed genoset object.
#' @param projectLabel project lebel to be added as metadata(genoSet)
#' @return save rds to file
#' @export
#'
meParseBeadStudio <- function(fileName, runName, projectLabel=NULL, runFolder="/Volumes/MacPro2TB/RESOURCES/Methylation_Infinium/ProcessedData/"){

   require(data.table)
   require(GenomicRanges)
   require(genoset)

   # Define and create output processed data dir
   myDir <- file.path(runFolder, runName)
   if(dir.exists(myDir)) stop("Run already exists. Set new name or delete/move old run")
   dir.create(myDir)

   # define samples and columns
   message("... begin parsing Me-file")
   my_columns <- as.character(read.table(fileName, sep="\t", nrows = 1, colClasses = "character"))
   #my_samples <- unique(my_samples)[-1]
   #my_columns <- as.character(apply(read.table(gdacFile, sep="\t", nrows = 2, skip=0, colClasses = "character"), 2, function(x) paste(x, collapse = " ")))

         # probe_id <- read.table(gdacFile, sep="\t", header=F, skip=2, colClasses = c("character", rep("NULL", length(my_samples))) )
      #3apa <- data.table::fread(gdacFile, drop = c(2:length(my_columns)), data.table = F)
   # Log
   logFile <- file.path(myDir, "log.txt")
   cat(file = logFile, "PARAMS\n")
   cat(file = logFile, append = T, "gdacFile:", fileName,"\n")
   cat(file = logFile, append = T, "Processd data directory:", myDir,"\n")

   cat(file = logFile, append = T, "\nDATA\n")
   cat(file = logFile, append = T, "n columns in gdac file:", length(my_columns) ,"\n")

   cat(file = logFile, append=T, "\n------\n")

   # Read data matrix (beta values)
   message("... reading Beta values from file")
   u <- grep("Beta", my_columns)
   beta_samples <- gsub(".AVG_Beta","",my_columns[u])
   #stopifnot(identical(beta_samples, my_samples))
   my_df <- data.table::fread(fileName, select = u,  data.table = F, skip=1, header = F, colClasses = "numeric")
   my_mat <- as.matrix(my_df)
   message("... rounding Beta to 4 decimals")
   my_mat <- round(my_mat, 4)
   head(my_mat)
   str(my_mat)
   str(beta_samples)
   cat(file = logFile, append = T, "dim beta matrix:", dim(my_mat) ,"\n")
   colnames(my_mat) <- beta_samples


   #head(my_mat[1:5, 1:5])
   # read pribe info
   uu <- which(my_columns=="TargetID")
   fdata_df <- data.table::fread(fileName, select = uu,  data.table = F, skip=0, header = T, colClasses = "character")
   #head(fdata)
   #colnames(fdata)
   cat(file = logFile, append = T, "dim feature id table:", dim(fdata_df) ,"\n")
   colnames(fdata_df)[1] <- c("FeatureID")
   rownames(fdata_df) <- fdata_df$FeatureID
   rownames(my_mat) <- fdata_df$FeatureID

   message("... Update Feature chr positions")
   message(" .... ... read FeatureID info from file")
   fdata450 <- get(load("~/RESOURCES/Methylation_Infinium/Annotations/probe.annotations.Methylation.Ringner.20031112.gRanges.Rdata"))
   fdata <- dplyr::left_join(fdata_df, as.data.frame(fdata450), by="FeatureID")
   #str(fdata2) # 'data.frame':	485577 obs. of  22 variables:
   # str(fdata) # 'data.frame':	485577 obs. of  5 variables:
   rownames(fdata) <- fdata$FeatureID
   stopifnot(identical(rownames(fdata), rownames(my_mat)))
   message("...removing NA chromosome rows")
   fdata  <- fdata[!is.na(fdata$seqnames), ]

   message("... building genoSet")
   cat(file = logFile, append = T, "\nOutput GenoSet\n")

   stopifnot(identical(beta_samples, colnames(my_mat)))
   fgr <- GRanges(seqnames=fdata$seqnames, ranges=IRanges(start=fdata$start, end=fdata$end), FeatureID=fdata$FeatureID)
   names(fgr) <- rownames(fdata$FeatureID)
   genome(fgr) <- "hg19"

   gs <- GenoSet(rowRanges = fgr, assays=list(beta=my_mat[fgr$FeatureID, ]), colData = data.frame(sample_id=beta_samples, stringsAsFactors = F, row.names=beta_samples))
   cat(file = logFile, append = T, "dim initial genoSet:", dim(gs) ,"\n")

   message("... keep only features that annotated to chromosomes")
   gs <- gs[seqnames(gs)%in%paste0("chr",c(1:22,"X","Y"))]
   gs <- toGenomeOrder(gs)
   seqlevels(gs) <- seqlevelsInUse(gs)
   cat(file = logFile, append = T, "dim final genoSet, chr1-chrY:", dim(gs) ,"\n")

   message("... adding project metadata")
   metadata(gs)[["projectLabel"]] <- projectLabel
   metadata(gs)[["meParseBeadStudio"]] <- timestamp()
   metadata(gs)[["fileName"]] <- fileName
   cat(file = logFile, append = T, "\nprojectLabelt:",  metadata(gs)[["projectLabel"]] ,"\n")
   cat(file = logFile, append = T, "meParseGDAC:",  metadata(gs)[["meParseGDAC"]] ,"\n")

   print(gs)
   message("... saving genoSet to file")
   outFile <- paste0(myDir, "/", runName, ".rds")
   saveRDS(gs, outFile)
   cat(file = logFile, append = T, "final genoSet saved to:", myDir ,"\n")
   }




#' get.HPA.images
#' @description Wrapper for importing Human Protein atlas images
#' @family import
#' @family human protein atlas
#' @param symbol gene symbol
#' @param run.name name of this run. needed to generate file names
#' @param hpa.filePaths.file file where HPA filepaths is stored. defaults to ~/RESOURCES/HPR/proteinatlas_v13_2011127.xml_all_tissues_normal_filePaths.Rdata
#' @return save images to disk
#' @export

get.HPA.images <- function(
   symbol = NULL,
   run.name = NULL,
   hpa.filePaths.file = '~/RESOURCES/HPR/proteinatlas_v13_2011127.xml_all_tissues_normal_filePaths.Rdata',
   tissue.type = 'kidney',
   hpa.images.dir = '/Volumes/MacPro2TB//RAW DATA/HPA/HPA_IMAGES/', # default dir to store all images downloaded from HPA url
   out.dir = '/Volumes//MacPro2TB//RAW DATA/HPA/HPA_images_geneSignatures/',
   copy.to.out.dir = F,
   plot.to.file = T
){
   cat('\n\n ... START FUNCTION get.HPA.images\n... ... ... ... ... ... ... ... ... ... ... ... ...\n\n')

   find.closest<-function(x, data){
      # x = input value
      # data = data vector to find the closest value in
      y<-data-x
      closest.i<-which(abs(y)==(min(abs(y))))
      return(closest.i)
   }
   # [1] "adrenal gland"     "appendix"          "bone marrow"       "breast"            "bronchus"          "cerebellum"        "cerebral cortex"   "cervix, uterine"
   # [9] "colon"             "duodenum"          "endometrium 1"     "endometrium 2"     "epididymis"        "esophagus"         "fallopian tube"    "gallbladder"
   # [17] "heart muscle"      "hippocampus"       "kidney"            "lateral ventricle" "liver"             "lung"              "lymph node"        "nasopharynx"
   # [25] "oral mucosa"       "ovary"             "pancreas"          "parathyroid gland" "placenta"          "prostate"          "rectum"            "salivary gland"
   # [33] "seminal vesicle"   "skeletal muscle"   "skin 1"            "skin 2"            "small intestine"   "smooth muscle"     "soft tissue 1"     "soft tissue 2"
   # [41] "spleen"            "stomach 1"         "stomach 2"         "testis"            "thyroid gland"     "tonsil"            "urinary bladder"   "vagina"


   require(jpeg)
   #hpr.reliability <- get(load("/Users/david/PROJECTS/HPR/HPA_reliability_normal_tissue_v13.csv.Rdata"))
   ab.validations <- get(load(file='~/RESOURCES/HPR/proteinatlas_v13_2011127.xml_HPA.reliability.AB.level.Rdata'))
   ab.validations <- ab.validations[-which(duplicated(ab.validations$AB)), ]



   if(is.null(run.name)) stop('you must provid valid run name')

   cat('\n\n ... loading Robject with HPA urls for all ABs and all Tissues: ... ', hpa.filePaths.file)
   hpr.files <- get(load(hpa.filePaths.file))
   match.arg(arg = tissue.type, choices = unique(hpr.files$tissue_type))
   hpr.files <- hpr.files[hpr.files$tissue_type==tissue.type,]

   stopifnot(file.exists(hpa.images.dir))
   stopifnot(file.exists(out.dir))
   my.dir <- paste0(out.dir, '/', run.name,'/')
   if(copy.to.out.dir) dir.create(my.dir)

   if(plot.to.file){
      my.count <- 1
      plot.no <- 1
      plot.file.name <- paste0(out.dir,run.name,'_',plot.no,'.jpg')

      tiff(file=plot.file.name, height=9,width=6,unit="in",res=1280,compression="jpeg")
      par(mfrow=c(9,6))
   } # end plot to file

   # First remove genes not present in list
   u <- match(symbol, hpr.files$SYMBOL)
   cat('\n ... ', length(which(is.na(u))), ' symbols not present in hpa filepaths; removing')
   symbol <- symbol[!is.na(u)]
   symbol <- symbol[!is.na(symbol)]


   for(i in 1:length(symbol)){
      symbol.i <- which(hpr.files$SYMBOL == symbol[i])

      for(ii in 1:length(symbol.i)){
         cat('\n\n... processing gene: ',symbol[i],' ... ', my.count)

         if(plot.to.file & my.count %in% seq(55, by = 54, length.out = 1000) ){
            dev.off()
            plot.no <- plot.no+1
            plot.file.name <- paste0(out.dir,run.name,'_',plot.no,'.jpg')
            tiff(file=plot.file.name,height=9,width=6,unit="in",res=1280,compression="jpeg")
            par(mfrow=c(9,6))
         } # end plot to file

         # download destfile to local HPA_IMAGES repository - if it does not exist
         destfile <- paste(hpa.images.dir, hpr.files$SYMBOL[symbol.i[ii]], '_', hpr.files$AB[symbol.i[ii]], '_',hpr.files$patient.id[symbol.i[ii]], '.jpg', sep='')
         if(file.exists(destfile)) cat('\n ... found HPA gene jpeg in local repository')
         if(!file.exists(destfile)) download.file(url=hpr.files$url[symbol.i[ii]], destfile=destfile)
         if(file.info(destfile)$size == 0){
            my.count <- my.count+1
            cat('\n\n... ... ',destfile,' is of size 0')
            next(ii)}

         # if to copy to out.dir
         if(copy.to.out.dir){
            destfile2 <- paste(my.dir, hpr.files$SYMBOL[symbol.i[ii]], '_', hpr.files$AB[symbol.i[ii]], '_',hpr.files$patient.id[symbol.i[ii]], '.jpg', sep='')
            file.copy(destfile, destfile2)}

         # if plot to file
         if(plot.to.file){
            cat('\n... plotting to file: ',plot.file.name)
            my.count <- my.count+1
            par(mar=c(0,0,0,0))
            plot(1:1280,1:1280,type="n",bty="n",xaxt="n",yaxt="n",xlab="",ylab="",xaxs = "i",yaxs = "i")
            img <- readJPEG(destfile, native=T)
            rasterImage(img,0,0,1280,1280)
            box(lwd=0.5)

            plot.caption <- paste(hpr.files$SYMBOL[symbol.i[ii]])
            text(x = 1280/2, y = 1280-50, labels = hpr.files$SYMBOL[symbol.i[ii]], col='red', cex=0.5)
            text(x = 1280/2, y = 1280-120, labels = paste0(hpr.files$AB[symbol.i[ii]],"  pid: ",hpr.files$patient.id[symbol.i[ii]]), col='black', cex=0.35)

            ab.info <- ab.validations[which(ab.validations$AB==hpr.files$AB[symbol.i[ii]]),2:4]
            foo<-function(x){
               if(is.na(x)) return('white')
               if(x=='supportive') return('green')
               if(x=='uncertain') return('orange')
               if(x=='non-supportive') return('red')
            }
            my.color <- sapply(ab.info, foo)
            legend(x = 50, y=200, legend = names(my.color), col = my.color, pch = 19, cex=0.2, box.lwd = 0.3)
         }

      } # end ii

      if(plot.to.file & !(length(symbol.i)%in% seq(from = 3,by = 3,length.out = 900))){
         u<-find.closest(length(symbol.i), seq(from = 3,by = 3,length.out = 900))
         uu<-seq(from = 3,by = 3,length.out = 900)[u]-length(symbol.i)
         if(uu==(-1)) uu <-2

         for(j in 1:uu){
            cat('\n ... ',i,' .. ',ii,' .. ',j)
            par(mar=c(0,0,0,0))
            plot(1:1280,1:1280,type="n",bty="n",xaxt="n",yaxt="n",xlab="",ylab="",xaxs = "i",yaxs = "i")
            my.count <- my.count+1
         } # end j
      } # end if not 3 images

   } # end i

   dev.off()
} # end function	get.HPA.images

   #

      # geo.id <- "GSE126441"
      # runName    <- "GPL13534_linehan_cimp"
   # projectLabel <- runName
   # runFolder <- "/Volumes/MacPro2TB/RESOURCES/Methylation_Infinium/Raw_data/"

    #     geo.id = "GSE105260", runName = "GSE105260_Nam_rcc", runFolder = "/Volumes/MacPro2TB/RESOURCES/Methylation_Infinium/Raw_data/", projectLabel = "GSE105260_Nam_rcc"


   # geo.id = "GSE35069", runName = "GSE35069_Reinius_Blood", runFolder = "/Volumes/MacPro2TB/RESOURCES/Methylation_Infinium/Raw_data/", projectLabel = "GSE35069_Reinius_Blood"
   # geo.id = "GSE113501", runName = "GSE113501_Evelonn_rcc", runFolder = "/Volumes/MacPro2TB/RESOURCES/Methylation_Infinium/Raw_data/", projectLabel = "GSE113501_Evelonn_rcc"


#' Download and parse GEO .soft methylation data into GenoSet
#' @description Wrapper for importing soft data from GEO using the GEOquery package
#' @description Methylation 450 data only
#' @description will download the full soft file and thereafter extract info
#' @family GEO
#' @family methylation
#' @family genoset
#' @family import
#' @param geo.id GEO GSE id
#' @param runName Fie dir-Name to use forthis run, e.g. rcc_me450.
#' @param runFolder Higher level directory to store processed genoset object
#' @param projectLabel  project lebel to be added as metadata(genoSet)
#' @param NormalizeBeta if to normalize beta infimium probe sets ringner/staaf style
#' @param platform specifiec what platform to extract. defaults to 450k GPL13534
#' @return save tables to disk
#' @export
meGetGEOsoft_me450 <- function(
   geo.id=NULL, runName=NULL, runFolder="/Volumes/MacPro2TB/RESOURCES/Methylation_Infinium/Raw_data/", projectLabel=NULL,
   NormalizeBeta = F, platform = "GPL13534"
){
   require(GEOquery)
   require(data.table)
   require(GenomicRanges)
   require(genoset)
   require(Biobase)

   # Define and create output processed data dir
   myDir <- file.path(runFolder, runName)
   if(dir.exists(myDir)) stop("Run already exists. Set new name or delete/move old run")
   dir.create(myDir)



   ## Download .soft file using getGEO
   gse <- getGEO(geo.id, GSEMatrix = FALSE)
   # show(gse)


   ## check platform
   stopifnot(platform %in% gse@header$platform_id)

   if(length(gse@header$platform_id)!=1){
      message("... NOTE!!!\n... >1 platforms available. Subsetting and using only the specified one")
      u <- lapply(gse@gsms, function(x){
         x@header$platform_id == platform
         })
      message(" ... ... from: ",length(gse@gsms))
      gse@gsms <- gse@gsms[unlist(u)]
      message(" ... ... to: ", length(gse@gsms))
   }
         #if(length(gse@header$platform_id)==1) u <- rep(TRUE, length(gse@gsms))


   ## get a reference for platform
   id_ref <- gse@gsms[[1]]@dataTable@table$ID_REF
   sample_names <- names(gse@gsms)

   ## beta values
   ## ::::::::
   message("... Extracting beta values. rounding Beta to 4 decimals")
   beta_list <- lapply(gse@gsms, function(x){
      stopifnot(identical(x@dataTable@table$ID_REF, id_ref))
      return(round(x@dataTable@table$VALUE, 4))
   })
      #str(beta_list)
   message("... ... creating Beta matrix")
   beta_mat <- matrix(ncol=length(sample_names), nrow=length(id_ref), data=unlist(beta_list), byrow = F, dimnames = list(id_ref, sample_names))
      #str(beta_mat)

   ## feature data
   ## ::::::::
   message("... Update Feaature chr positions")
   message(" .... ... read FeatureID info from file")
   fdata450 <- get(load("~/RESOURCES/Methylation_Infinium/Annotations/probe.annotations.Methylation.Ringner.20031112.gRanges.Rdata"))
   fdata <- dplyr::left_join(data.frame(FeatureID=id_ref, stringsAsFactors = F), as.data.frame(fdata450), by="FeatureID")
   #str(fdata2) # 'data.frame':	485577 obs. of  22 variables:
   # str(fdata) # 'data.frame':	485577 obs. of  5 variables:
   rownames(fdata) <- fdata$FeatureID
   stopifnot(identical(rownames(fdata), rownames(beta_mat)))
   message("...removing NA chromosome rows")
   fdata  <- fdata[!is.na(fdata$seqnames), ]


   ## Build GenoSet
   ## ::::::::::::::
   message("... building genoSet")
   # cat(file = logFile, append = T, "\nOutput GenoSet\n")

   fgr <- GenomicRanges::GRanges(seqnames=fdata$seqnames, ranges=IRanges(start=fdata$start, end=fdata$end), FeatureID=fdata$FeatureID)
   names(fgr) <- rownames(fdata$FeatureID)
   genome(fgr) <- "hg19"

   stopifnot(identical(rownames(beta_mat[fgr$FeatureID, ]), fgr$FeatureID))
   stopifnot(identical(sample_names, colnames(beta_mat)))

   gs <- genoset::GenoSet(rowRanges = fgr, assays=list(beta=beta_mat[fgr$FeatureID, ]), colData = data.frame(sample_id=sample_names, stringsAsFactors = F, row.names=sample_names))
      # cat(file = logFile, append = T, "dim initial genoSet:", dim(gs) ,"\n")
      #str(gs[,,"beta"])

   ## Filter GenoSet
   ## ::::::::::::::
   message("... keep only features that annotated to chromosomes")
   gs <- gs[GenomicRanges::seqnames(gs)%in%paste0("chr",c(1:22,"X","Y"))]
   gs <- genoset::toGenomeOrder(gs)
   GenomeInfoDb::seqlevels(gs) <- GenomeInfoDb::seqlevelsInUse(gs)
   # cat(file = logFile, append = T, "dim final genoSet, chr1-chrY:", dim(gs) ,"\n")



   ## update sample coldata
   ## :::::::::::::::
   message("... extracting GSM sample meta data for colData slot")
   coldata_list <- lapply(gse@gsms, function(x){
      xx <- x@header

      u <- sapply(xx, length)
      y <- data.frame(xx[!u>1])

      if(length(which(u>1))>0){

         xxx <- xx[u>1]
         xxx <- unlist(xxx)
         y2 <- data.frame(t(data.frame(xxx)))
         y <- cbind(y, y2)
      }
      rownames(y) <- y$geo_accession
      return(y)
   })
   coldata_df <- dplyr::rbind_list(coldata_list)
   coldata_df <- data.frame(coldata_df)
   # coldata_df <- do.call("rbind", coldata_list)
   #str(coldata_df)

   rownames(coldata_df) <- names(gse@gsms)
   stopifnot(identical(colnames(gs), rownames(coldata_df)))
   colData(gs) <- DataFrame(coldata_df)



   # IF EXTRA SAMPLE DATA TABLE EXISTS
   df1 <- getGSEDataTables(geo.id)
   if(length(df1)!=0) {
         message("... NOTE!: extra sample info exists, saving RDS file")
         outFile <- paste0(myDir, "/", runName, "_SampleDataTable.rds")
         saveRDS(df1, outFile)
         }




   ## add metadata
   ## ::::::::::::::::
   message("... adding project metadata")

   metadata(gs)[["projectLabel"]] <- projectLabel
   metadata(gs)[["me_getGEO_methylation450"]] <- timestamp()
   metadata(gs)[["fileName"]] <- fileName
   #cat(file = logFile, append = T, "\nprojectLabelt:",  metadata(gs)[["projectLabel"]] ,"\n")
   #cat(file = logFile, append = T, "meParseGDAC:",  metadata(gs)[["meParseGDAC"]] ,"\n")
   metadata(gs)[["gse_header"]] <- gse@header


   ## Save to file
   ## ::::::::::
   print(gs)
   message("... saving genoSet to file")
   outFile <- paste0(myDir, "/", runName, ".rds")
   saveRDS(gs, outFile)
   outFile <- paste0(myDir, "/", runName, "_sampleTable.txt")
   write.table(colData(gs), file = outFile, sep="\t", row.names=F, quote=F)



   ## IF normalize beta (run meNormalizeBeta)
   ## :::::
   if(NormalizeBeta){
      dlfoo2::meNormalizeBeta(
         x = gs, runName = paste0(runName,"_norm"), runFolder = runFolder, projectLabel = projectLabel
      )
   }
} # me_getGEO_methylation450







      # pdata <- read.delim("/Volumes/MacPro2TB/RESOURCES/Methylation_Infinium/Raw_data/GSE122126_moss_IDAT/GSE122126-GPL21145_series_matrix_sampleTable.txt")
      # pdata <- pdata %>%
      #    mutate(Basename = gsub(".*suppl[/]GSM","GSM",as.character(supplementary_file_1))) %>%
      #    mutate(Basename = gsub("_Grn.*","",Basename)) %>%
      #    dplyr::filter(sample.type!="cfDNA") %>% dplyr::filter(disease.state=="Normal")
      #
      # idatFolder <- "/Volumes/MacPro2TB/RESOURCES/Methylation_Infinium/Raw_data/GSE122126_moss_IDAT/"

#' Download and parse GEO .soft methylation data into GenoSet
#' @description Wrapper for importing soft data from GEO using the GEOquery package
#' @description Methylation 450 data only
#' @description will download the full soft file and thereafter extract info
#' @family GEO
#' @family methylation
#' @family genoset
#' @family import
#' @param pdata sample sheet with idat file names. used as pdata.
#' @param pdata.sampleName.column column to use for sample names. must be unique
#' @param idatFolder top level folder for IDATs
#' @param runName Fie dir-Name to use forthis run, e.g. rcc_me450.
#' @param runFolder Higher level directory to store processed genoset object
#' @param projectLabel  project lebel to be added as metadata(genoSet)
#' @return save tables to disk
#' @export
meParseIDAT <- function(
   pdata=NULL, pdata.sampleName.column="geo_accession", idatFolder=NULL, runName=NULL, runFolder="/Volumes/MacPro2TB/RESOURCES/Methylation_Infinium/Raw_data/",
   projectLabel=NULL, NormalizeBeta = F
){
   require(minfi)
   require(GEOquery)
   require(data.table)
   require(GenomicRanges)
   require(genoset)
   require(Biobase)
   require(dplyr)
   require( tidyr)
   # browseVignettes("minfi")
   # require(GEOmetadb)

   # Define and create output processed data dir
   myDir <- file.path(runFolder, runName)
   if(dir.exists(myDir)) stop("Run already exists. Set new name or delete/move old run")
   dir.create(myDir)


   if(!("Basename" %in% colnames(pdata))) stop("pdata must have IDAT Basename as column. should be IDAT full filename without e.g. _Grn.idat.gz")
   sample_names <- pdata[,"geo_accession"]
   rownames(pdata) <- sample_names
   if(is.null(sample_names)) stop("pdata.sampleName.column problem")
   if(any(duplicated(sample_names))) stop("pdata.sampleName.column problem")

   ## Read Methylation IDATs
   setwd(idatFolder)
   RGset <- minfi::read.metharray.exp(targets = pdata, force = T)
   MSet <- preprocessIllumina(RGset)
   # str(getBeta(MSet))

   ## convert MethylSet to genoset
    #  class(MSet)
   beta_mat <- round(getBeta(MSet),4)
   str(beta_mat)
   str(pdata)
   if(ncol(beta_mat)!=nrow(pdata)) message("\n NOTE!! NOT ALL IDAT LOADED")
   pdata <- pdata[match(colnames(beta_mat), pdata$Basename), ]
   stopifnot(identical(colnames(beta_mat), pdata$Basename))
   colnames(beta_mat) <- sample_names

   ## feature data
      ## ::::::::
      message("... Update Feaature chr positions")
      message(" .... ... read FeatureID info from file - USE EPIC array matched hg38 annotation from Zhou et al")
      if(annotation(MSet)["array"]=="IlluminaHumanMethylationEPIC") fdata <- readRDS("~/RESOURCES/Methylation_Infinium/Annotations/Zhou_Infinium_Annotations/EPIC.hg38.manifest.rds")
      if(annotation(MSet)["array"]=="IlluminaHumanMethylation450k") fdata <- readRDS("~/RESOURCES/Methylation_Infinium/Annotations/Zhou_Infinium_Annotations/hm450.hg38.manifest.rds")
      my_probes <- intersect(rownames(beta_mat), names(fdata))
      fdata <- fdata[my_probes]
      genome(fdata) <- "hg38"

      ## Build GenoSet
      ## ::::::::::::::
      message("... building genoSet")
        gs <- genoset::GenoSet(rowRanges = fdata, assays=list(beta=beta_mat[names(fdata), ]), colData = data.frame(sample_id=sample_names, stringsAsFactors = F, row.names=sample_names))
      gs <- genoset::toGenomeOrder(gs)
      GenomeInfoDb::seqlevels(gs) <- GenomeInfoDb::seqlevelsInUse(gs)

      message("... adding colData slot")
      coldata_df <- pdata
      stopifnot(identical(colnames(gs), rownames(coldata_df)))
      colData(gs) <- DataFrame(coldata_df)

   ## add metadata
   ## ::::::::::::::::
   message("... adding project metadata")
   metadata(gs)[["projectLabel"]] <- projectLabel
   metadata(gs)[["minfi.IDAT"]] <- timestamp()
   metadata(gs)[["IDAT.dir"]] <- getwd()

   ## SAVE

   saveRDS(gs, file.path(myDir,paste0(runName,".rds")))

}






   # require(genoset)
   #my_file <- "/Volumes/MacPro2TB/RESOURCES/Methylation_Infinium/Raw_data/GSE35069_Reinius_Blood_norn/GSE35069_Reinius_Blood_norn.rds"
   # gs <- readRDS(my_file)
   # gs
   # pdata <- colData(gs)

#' small script for polisihing pdata seriesmatrix from geo
#' @description Wrapper for importing soft data from GEO using the GEOquery package
#' @description Methylation 450 data only
#' @description will download the full soft file and thereafter extract info
#' @family GEO
#' @family series matrix
#' @family genoset
#' @family import
#' @param pdata pdata from geo series matrix entry
#' @export
geoSeriesMatrixPdataPolish <- function(
   pdata=NULL
){

   pdata <- as.data.frame(pdata, stringsAsFactors=F)
   my_colnames <- colnames(pdata)
   if(any(grepl("characteristics_ch1", my_colnames))){
      my_colnames <- gsub("characteristics_ch1*.","characteristics_ch1",my_colnames)
      #my_colnames <- gsub("description*.","description",my_colnames)
   }
   if(any(grepl("Sample_", my_colnames))) my_colnames <- gsub("*.Sample_","",colnames(pdata))
   # head(pdata)

   u <- grep("characteristics_ch1", my_colnames)
   if(length(u)){
      my_colnames2 <- my_colnames
      for(i in 1:length(u)){
         pdata[,u[i]] <- as.character(pdata[,u[i]])
         uu <- strsplit(pdata[,u[i]], ": ")
         ui <- grep(": ", pdata[,u[i]])[1]
         uu[!grepl(": ", pdata[,u[i]] )] <- ""

         if(length(uu[[ui[1]]]) >1 ){
            my_colnames2[u[i]] <- uu[[ui]][1]
            pdata[,u[i]] <- unlist(lapply(uu, function(x) paste(x[-1], collapse = ": ")))
         }
      }
   my_colnames2 <- gsub(" |[/]",".",my_colnames2)
   if(any(duplicated(my_colnames2))){
      my_dups <- unique(my_colnames2[duplicated(my_colnames2)])
      for(i in 1:length(my_dups)){
         u <- which(my_colnames2 %in% my_dups[i])
         my_colnames2[u] <- paste(my_dups[i], 1:length(u), sep="_")
      }
   }
   colnames(pdata) <- my_colnames2
   }

   ## polish special characters
   colnames(pdata) <- gsub("[(]|[)]","" ,colnames(pdata))
   colnames(pdata) <- gsub("[-]","." ,colnames(pdata))
   colnames(pdata) <- gsub(" ","" ,colnames(pdata))
   colnames(pdata) <- gsub("[.][.]","." ,colnames(pdata))

   return(pdata)
}


      # file.name="/Volumes/MacPro2TB/RESOURCES/Methylation_Infinium/Raw_data/GSE122126_moss_IDAT/GSE122126-GPL13534_series_matrix.txt"

# file.name = "/Volumes/MacPro2TB/RESOURCES/ChIPseq/zanconato_yaptaz_2018_GSE102409/GSE102409_series_matrix.txt"

#' parse GEO series matrix to get sampel info sheet
#' @description Wrapper for importing soft data from GEO using the GEOquery package
#' @description Methylation 450 data only
#' @description will download the full soft file and thereafter extract info
#' @family GEO
#' @family series matrix
#' @family genoset
#' @family import
#' @param geo.id GEO GSE id
#' @param runName Fie dir-Name to use forthis run, e.g. rcc_me450.
#' @param runFolder Higher level directory to store processed genoset object
#' @return save tables to disk
#' @export
parseGeoSeriesMatrixSampleFile <- function(
   file.name=NULL
){
   u <- grep("!Sample_title", readLines(file.name))
   uu <- grep("!series_matrix", readLines(file.name))[1]
   my_tab <- read.delim(file.name, skip=u-1, , header = F, nrows = uu-u)
   my_colnames <- gsub("!Sample_","",my_tab[,1])
   pdata <- data.frame(t(my_tab[,-1]), stringsAsFactors = F)

   # u <- which(my_colnames=="characteristics_ch1")
   # if(length(u)){
   #    my_colnames2 <- my_colnames
   #    for(i in 1:length(u)){
   #       uu <- strsplit(pdata[,u[i]], ": ")
   #       ui <- grep(": ", pdata[,u[i]])[1]
   #       uu[!grepl(": ", pdata[,u[i]])] <- ""
   #
   #       if(length(uu[[1]])>1){
   #          my_colnames2[i] <- uu[[ui]][1]
   #          pdata[,u[i]] <- unlist(lapply(uu, function(x) paste(x[-1], collapse = ": ")))
   #       }
   #    }
   # }
   # my_colnames2 <- gsub(" |[/]",".",my_colnames2)
   # if(any(duplicated(my_colnames2))){
   #    my_dups <- unique(my_colnames2[duplicated(my_colnames2)])
   #    for(i in 1:length(my_dups)){
   #       u <- which(my_colnames2 %in% my_dups[i])
   #       my_colnames2[u] <- paste(my_dups[i], 1:length(u), sep="_")
   #    }
   # }
   #colnames(pdata) <- my_colnames2

   colnames(pdata) <- my_colnames
   pdata <- geoSeriesMatrixPdataPolish(pdata)
   file.name.out <- gsub(".txt|.csv", "_sampleTable.txt", file.name)
   message("... writing sampleTable to disk:", file.name.out)
   write.table(pdata, file=file.name.out, row.names = F, quote=F, sep="\t")
}
