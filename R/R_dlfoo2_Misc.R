#' copy a vetor to memory
#' @param x vecor to be copied
#' @param sep what delim
#' @return copied to memory to be pasted elsewhere
#' @export
copyToMemory <- function(x, sep="\t"){
   clip <- pipe("pbcopy", "w")
   write.table(x, file=clip, sep = '\t', row.names = FALSE)
   close(clip)
}





#' splits a tcga barcode into its defined codes
#' @description  Version 1.0 20190108
#' @family misc
#' @family tcga
#' @param x TCGA aliquot barcode id ('-' or '_' separated)
#' @return barcode character vector of length 1
#' @export
#'
tcgaBarcodeSplitter <- function(x){
   if(nchar(x)!=28) {stop('Aliquot barcode does not have 28 digits')}
   bc.df <- data.frame(
      bc.patient = substr(x, 9, 12),
      bc.type = substr(x, 14, 15),
      bc.vial = substr(x, 16, 16),
      bc.portion = substr(x, 18, 19),
      bc.analyte =  substr(x, 20, 20),
      bc.plate =  substr(x, 22, 25),
      bc.center = substr(x, 27, 28)
      )
   return(bc.df)
} # end function tcgaBarcodeSplitter


#' splits a tcga barcode into its defined codes
#' @description Version 1.0 20190506
#' @family misc
#' @family tcga
#' @param x TCGA aliquot barcode id ('-' or '_' separated)
#' @return barcode character vector of length 1
#' @export
#'
tcgaBarcodeSplitter_v2 <- function(x){
   # if(nchar(x)!=28) {stop('Aliquot barcode does not have 28 digits')}
   y <- strsplit(x, split = "-" )[[1]]
   bc_tab <- data.frame(
      patient_barcode = NA,
      sample_barcode = NA,
      vial_barcode = NA,
      portion_barcode = NA,
      analyte_barcode = NA
   )

   if(length(y)>=3) bc_tab$patient_barcode = paste(y[1:3], collapse = "-")
   if(length(y)>=4){
      yy <- strsplit(y[4], split="")[[1]]
      bc_tab$sample_barcode = paste(c(bc_tab$patient_barcode, paste0(yy[1:2], collapse = "")), collapse = "-")
      if(length(yy)>=3) bc_tab$vial_barcode = paste(c(bc_tab$patient_barcode, paste0(yy, collapse = "")), collapse = "-")
   }
   if(length(y)>=5){
      yy <- strsplit(y[5], split="")[[1]]
      bc_tab$portion_barcode = paste(c(paste(c(bc_tab$patient_barcode, y[4]), collapse = "-"),
               paste0(yy[1:2], collapse = "")), collapse = "-")
      if(length(yy)>=3) bc_tab$analyte_barcode = paste(c(paste(c(bc_tab$patient_barcode, y[4]), collapse = "-"),
               paste0(yy, collapse = "")), collapse = "-")
      }
   return(bc_tab
      )
} # end function tcgaBarcodeSplitter


#' Create short inhouse name for TCGA and TARGET samples
#' @description  Shortnames for TCGA samples using barcode and project
#' @family misc
#' @family tcga
#' @param barcode, the TCGA aliquot barcode id, e.g. TCGA-KL-8323-01A-21R-2315-07
#' @param project, the TCGA project id, e.g. BRCA
#' @param program, TCGA or TARGET initiative
#' @return character string of short barcode name
#' @export
tcga_shortname <- function(barcode, project, program="TCGA"){
   sample_types <- list(
      s01=c("TP","Primary solid Tumor", "Tumor"),
      s02=c("TR",	"Recurrent Solid Tumor", "Tumor"),
      s03=c("TB",	"Primary Blood Derived Cancer - Peripheral Blood", "Tumor"),
      s04=c("TRBM",	"Recurrent Blood Derived Cancer - Bone Marrow", "Tumor"),
      s05=c("TAP",	"Additional - New Primary", "Tumor"),
      s06=c("TM",	"Metastatic", "Tumor"),
      s07=c("TAM",	"Additional Metastatic", "Tumor"),
      s08=c("THOC","Human Tumor Original Cells", "Tumor"),
      s09=c("TBM",	"Primary Blood Derived Cancer - Bone Marrow", "Tumor"),
      s10=c("NB",	"Blood Derived Normal", "Normal"),
      s11=c("NT",	"Solid Tissue Normal", "Normal"),
      s12=c("NBC",	"Buccal Cell Normal"),
      s13=c("NEBV",	"EBV Immortalized Normal"),
      s14=c("NBM,	Bone Marrow Normal", "Normal"),
      s20=c("CELLC",	"Control Analyte", "Normal"),
      s40=c("TRB",	"Recurrent Blood Derived Cancer - Peripheral Blood", "Tumor"),
      s50=c("CELL",	"Cell Lines", "Tumor"),
      s60=c("XP",	"Primary Xenograft Tissue", "Tumor"),
      s61=c("XCL",	"Cell Line Derived Xenograft Tissue", "Tumor")
   )
   if(nchar(barcode)<=15) {stop(' barcode does not have 28 digits')}

   if(program=="TCGA"){
      patient = substr(barcode, 9,12)
      sample.type = substr(barcode, 14,15)
      sample.vial = substr(barcode, 16,16)
      portion = substr(barcode, 18,19)
      sample.type <- sample_types[[paste0("s",sample.type)]][1]
   }

   if(program=="TARGET"){
      patient = substr(barcode, 11, 16)
      sample.type = substr(barcode, 18,19)
      sample.vial = substr(barcode, 19,19)
      portion = substr(barcode, 21,22)
      sample.type <- sample_types[[paste0("s",sample.type)]][1]
   }
   sample.id <- paste(project, sample.type, patient, sep="_")
   return(sample.id)
} # end sub-function create.tcga.shortnames



#' Get individual sample pdata
#' @description Uses dlfoo2::pdata_tcgaPanCan pData table to get sample info
#' @description Incuding the sample_id (tcga shortname) for TCGA and TARGET samples: project + sample type + patient
#' @description This Id can be used to map samples between different platsforms - not sensitive to e.g. DNA/RNA type
#' @description Wrapper that uses tcga_shortname function as well as the PanCan curated sample lists
#' @description NOTE 1: Samples excluded by PanCan will be given no name
#' @description NOTE 2: Replicate id's may be generated here
#' @family misc
#' @family tcga
#' @param aliquot.barcode, the TCGA aliquot barcode id, e.g. TCGA-KL-8323-01A-21R-2315-07
#' @param platform  what platform from choices; "Genome_Wide_SNP_6","HumanMethylation27","HumanMethylation450","IlluminaGA_miRNASeq","IlluminaGA_RNASeqV2","IlluminaHiSeq_DNASeqC","IlluminaHiSeq_miRNASeq","IlluminaHiSeq_RNASeqV2","MDA_RPPA_Core")
#' @return character string of short barcode name
#' @export
tcgaPdataPanCan <- function(
      aliquot.barcode,
      platform = NULL
   ){
   require(dplyr)
   require(tidyverse)
   my_platforms <- c("Genome_Wide_SNP_6","HumanMethylation27","HumanMethylation450","IlluminaGA_miRNASeq","IlluminaGA_RNASeqV2","IlluminaHiSeq_DNASeqC","IlluminaHiSeq_miRNASeq","IlluminaHiSeq_RNASeqV2","MDA_RPPA_Core")
   if(is.null(platform)) stop("platform must be defined")
   platform = match.arg(platform, choices=my_platforms, several.ok = TRUE)

   x <- aliquot.barcode
   x <- unlist(lapply(x, function(y){gsub("[.]","-",y)}))

   if(class(x) != "character") stop("aliquot.barcode must be a character vector")

   pancan_df <- dlfoo2::pdata_tcgaPanCan
   pancan_df <- pancan_df[pancan_df$platform %in% platform,]

   #u <- match(x, pancan_df$aliquot_barcode)
   #uu <- pancan_df$aliquot_barcode[is.na(u)]
   pancan_df <- pancan_df[pancan_df$aliquot_barcode %in% x, ]
   str(pancan_df)

   if(nrow(pancan_df)<1) stop("no matching aliquot ids in TCGA PanCan file")
   message(" ... identified ", length(unique(pancan_df$aliquot_barcode)), " of ", length(x), " aliquot ids in TCGA PanCan file")
   #message(" ... NOTE!: ",length(unique(duplicated(pancan_df$aliquot_barcode))), " aliquot ids are found on duplicate rows in TCGA PanCan file")
   message(" ... returning data frame with ", nrow(pancan_df), " rows")
   # pancan_df <- pancan_df[u[!is.na(u)],]
   my_df <- as.data.frame(pancan_df %>% dplyr::select(sample_id, everything()))
   table(my_df$Do_not_use)
   #message("... returning dataframe with sample_id aliquot_barcode with", nrow(my_df), " rows")
   return(my_df)
}





#' Correlation matrix to table
#' @description Uses a correlation martrix to create a table with correlations for all possible combinations
#' @family old functions
#' @family transformation
#' @family misc
#' @param cor.mat, a correlation matrix
#' @return a data frame with all combinations of correlations
#' @export
cormat2indiv <- function(cor.mat){
   r.n <- nrow(cor.mat)
   r.i <- c(1:(r.n-1))
   c.i <- sapply(r.i, function(x, y=(r.n)) seq(from=x+1, to = y))
   c.length <- unlist(lapply(c.i, length))
   cor1 <- rep(r.i, c.length)
   cor2 <- unlist(c.i)
   c.mat <- matrix(nrow=length(cor1), data=c(cor1, cor2), byrow=F)

   my.cor <- apply(c.mat, 1, function(x) cor.mat[x[1],x[2]])
   out.tab <- data.frame(name1=colnames(cor.mat)[c.mat[,1]], name2=colnames(cor.mat)[c.mat[,2]], cor=my.cor, stringsAsFactors = F)
   return(out.tab)
      }

