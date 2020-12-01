






# hif_canonical_smythies <- readRDS("~/RESOURCES/ChIPseq/Smythies_hif_2018_GSE120887/ProcessedData/hif_smythies/hif_smythies_GRangesList_Canonical.rds")
#    names(hif_canonical_smythies)
#
#   x <- hif_canonical_smythies[["HKC8_HIF2a_hypo05_kdWT"]]
#

#' Wrapper Liftover using AnnotaitonHub package
#' Using AnnotaionHub  \code{\link[AnnotationHub]{AnnotationHub}} to perform liftover between genome bulds.
#' @family methylation
#' @family granges
#' @family genoset
#' @family genomics
#' @family transformation
#' @param x A \code{\link[GenomicRangges]{GRanges}} object. x must have genome specified. all sequnces must be from one unique genome.
#' @param to what genome build to transfer to. Default is hg19
#' @param collapse.to.unique If output should contain uniqe ids (per row) or if multiple rows are allowed for out. The non-collapsed table will be stored in metadata slot
#' @param keep.chromosomes what chrosomosomes to keep. set to NULL if to keep all chromosomes
#' @return GRanges object with genomic start/stop for supplied identifiers (or all if 'all' option used)
#' @export
#'
liftoverGRanges = function(
   x = NULL,
   id_type = c("ensembl_gene_id","hgnc_symbol", "entrezgene"),
   to = c('hg19'),
   keep.chromosomes = c(1:22, "X")
   ){
   require(AnnotationHub)
   stopifnot(class(x)=="GRanges")
   genome_choices <- c('hg19','hg38')
   if(length(unique(genome(x)))!=1) stop("no or mere than one genome specified for x")
   from = match.arg(unique(genome(x)), choices = genome_choices, several.ok = F)
   message(" --> chaining from genome: ", from)
   to = match.arg(to, choices=genome_choices, several.ok = F)
   if(from==to) stop("Argument 'to' genome cannot be the same as genome for x")

   message("... connecting to AnnotationHub")
   ah <- AnnotationHub()
   message("... ... querying ChainFiles")
   ahub.chain <- subset(ah, rdataclass == "ChainFile" & species == "Homo sapiens")
   #query(ahub.chain, c(to, from))
   #chain <- ahub.chain[ahub.chain$title == "hg38ToHg19.over.chain.gz"]
   u <- query(ahub.chain, paste(from,to, sep="To"))
   message("... ... using chain: ", names(u), "; ", u$title)
   chain <- ahub.chain[[names(u)]]
   # chain <- chain[[1]]
   xlo <- liftOver(x, chain)

   message("... ... n ranges pre liftOver: ", length(x))
   message("... ... n ranges post liftOver: ", length(unlist(xlo)))

   xlo <- unlist(xlo) # GRanges object with 564380
   return(xlo)
   }



#' expand ranges up/downstream of a GRanges object
#' From https://support.bioconductor.org/p/78652/
#' @family genomics
#' @family granges
#' @family transformation
#' @param x GRanges object
#' @param upstream no bases upstream to expand
#' @param downstream no bases downstream to expand
#' @return GRanges object with expanded ranges from start stop as defined
#' @export
#'
expandRanges <- function(x, upstream=0, downstream=0) {
      require(GenomicRanges)
      if(class(x)!= "GRanges") stop()
      if (any(strand(x) == "*"))
           warning("'*' ranges were treated as '+'")
       on_plus <- strand(x) == "+" | strand(x) == "*"
       new_start <- start(x) - ifelse(on_plus, upstream, downstream)
       new_end <- end(x) + ifelse(on_plus, downstream, upstream)
       ranges(x) <- IRanges(new_start, new_end)
       trim(x)
       return(x)
   }


#' Load local UCSC refGene object corresponding to the genome of input genomic Ranges object
#' Load local UCSC refGene object corresponding to the genome of input genomic Ranges object. Uses the genome slot
#' @family old functions
#' @family genomics
#' @family granges
#' @param gr.obejct a genomic ranges object whose matching refgene to be loaded
#' @return a granges object
#' @export
#'
   loadRefgene <- function(gr.obejct){
      refgene.list <- list(
                        mm8 = "/Users/david/R/R_FUNCTIONS_AND_REFERENCE_FILES/refgene.mm8.transcripts.20151026.gr.Rdata",
                        mm9 = "/Users/david/R/R_FUNCTIONS_AND_REFERENCE_FILES/refgene.mm9.transcripts.20151021.gr.Rdata",
                        mm10 = "/Users/david/R/R_FUNCTIONS_AND_REFERENCE_FILES/refgene.mm10.transcripts.20151021.gr.Rdata",
                        hg19 = "/Users/david/R/R_FUNCTIONS_AND_REFERENCE_FILES/refgene.hg19.transcripts.20151021.gr.Rdata",
                        hg18 = "/Users/david/R/R_FUNCTIONS_AND_REFERENCE_FILES/refgene.hg18.transcripts.20151026.gr.Rdata",
                        rn6 = "/Users/david/R/R_FUNCTIONS_AND_REFERENCE_FILES/refgene.rn6.transcripts.20151116.gr.Rdata"
                           )

      my.genome <- unique(genome(gr.obejct))
      if(length(my.genome)!=1) stop("cannot determine genome")
      u <- match(my.genome, names(refgene.list))
      if(length(u)<1) stop(paste("no refegene file for genome",my.genome,"defined in function"))
      cat("\n ... loading refgene of genome ",my.genome," from local file:", refgene.list[[u]])
      refgene <- get(load(file=refgene.list[[u]]))
      return(refgene)
   }



#' Extract chip calls in the promoter of genes
#' Wrapper using GenomicRanges::promoters to extract chip peaks from a genomic ranges object located in gene promoters (obtained from UCSC refgene object).
#' @family old functions
#' @family genomics
#' @family granges
#' @family bed
#' @param my.bed a bed-style genomicRanges object with peak positions
#' @param upstream,number of bases upstream TSS
#' @param downstream, .. of TSS
#' @param refgene, what local refgene gRanges object to use
#' @return a genomic ranges object of the chip peaks contained in gene promoters
#' @export
chip_peak2gene <- function(my.bed=NULL, upstream=2000, downstream=200,  refgene=NULL){
   my.bed$input.index <- 1:length(my.bed)

   if(is.null(refgene)) refgene <- dlfoo::loadRefgene(my.bed)
   cat("\n ... extracting promoter regions from refgene using ",upstream,"bases upstream and", downstream,"bases downstream TSS")
   refgene.prom <- GenomicRanges::promoters(refgene, upstream=upstream, downstream = downstream)
   my.overlaps <- as.matrix(findOverlaps(my.bed, refgene.prom, minoverlap=1))

   prom.gr <- refgene[my.overlaps[,2],]
   elementMetadata(prom.gr) <- cbind(elementMetadata(prom.gr), elementMetadata(my.bed[my.overlaps[,1],]))
   return(prom.gr)
    }


#' From a gRanges object, get the max peak value from all peaks within one gene promoter
#' @family old functions
#' @family genomics
#' @family granges
#' @family bed
#' @family chip
#' @param peak2gene.gr, a bed-style genomicRanges object with peak positions from chip_peak2gene function
#' @param value.column, names of chip peak value column
#' @param refgene.ref.column, what refgene column to match on
#' @param refgene.cols, what refgene columns to output (condense)
#' @param refgene, refgene, what local refgene gRanges object to use
#' @return a genomic ranges object with one row per queried refgene entry
#' @export
chip_peak2gene_maxValue <- function(peak2gene.gr, value.column=NULL, refgene.ref.column="GENEID", refgene.cols=c("GENEID","SYMBOL"), refgene=NULL){
      if(is.null(refgene)) refgene <- loadRefgene(peak2gene.gr)
      cat("\n ... Fetching max values for the supplied value columns collapsed on ", refgene.ref.column)
      u <- match(refgene.ref.column, colnames(elementMetadata(peak2gene.gr)))
      if(length(u)!=1 | is.na(u)) stop("cant find refgene column to collapse on")
      u<-match(value.column, colnames(elementMetadata(peak2gene.gr)))
      if(length(u)!=1 | is.na(u)) stop("cant find supplied value column in elementMetadata")
      apa <- split(elementMetadata(peak2gene.gr), elementMetadata(peak2gene.gr)[,"GENEID"])
      my.vals <- lapply(apa, function(x) max(x[,value.column], na.rm=T))
      u <- match(names(my.vals), elementMetadata(refgene)[,refgene.ref.column])
      apa2 <- cbind(as.data.frame(refgene[u,])[,refgene.cols], unlist(my.vals))
      colnames(apa2)[length(refgene.cols)+1] <- value.column
      return(apa2)
}

#' From a gene/promoter-centered gRanges object - get 1) the number of peaks per promoter and 2) the max peak value from all peaks within one gene promoter
#' @family old functions
#' @family genomics
#' @family granges
#' @family chip
#' @param peak2gene.gr, a bed-style genomicRanges object with peak positions from chip_peak2gene function
#' @param value.column, names of chip peak value column
#' @param refgene.ref.column, what refgene column to match on
#' @param refgene.cols, what refgene columns to output (condense)
#' @param refgene, refgene, what local refgene gRanges object to use
#' @return a data frame with one row per queried refgene entry
#' @export
chip_genePeaks2occupancyTab <- function(peak2gene.gr, value.column=NULL, refgene.ref.column="GENEID", refgene.cols=c("GENEID","SYMBOL"), refgene=NULL){
      if(is.null(refgene)) refgene <- loadRefgene(peak2gene.gr)
      cat("\n ... Fetching nHits & max values for the supplied value columns collapsed on ", refgene.ref.column)
            u <- match(refgene.ref.column, colnames(elementMetadata(peak2gene.gr)))
            if(length(u)!=1 | is.na(u)) stop("cant find refgene column to collapse on")
            u<-match(value.column, colnames(elementMetadata(peak2gene.gr)))
            if(length(u)!=1 | is.na(u)) stop("cant find supplied value column in elementMetadata")
            apa <- S4Vectors::split(GenomicRanges::elementMetadata(peak2gene.gr), GenomicRanges::elementMetadata(peak2gene.gr)[,"GENEID"])
            my.vals <- lapply(apa, function(x) c(nrow(x), max(x[,value.column], na.rm=T)))
            my.vals.mat <- matrix(nrow=length(my.vals), ncol=2, data=unlist(my.vals), byrow=T, dimnames = list(names(my.vals), c("nHits","maxVal")))
            colnames(my.vals.mat) <- c("nHits","maxVal")
            u <- match(names(my.vals), GenomicRanges::elementMetadata(refgene)[,refgene.ref.column])
            occupance.tab <- cbind(as.data.frame(refgene[u,])[,refgene.cols], my.vals.mat)
            return(occupance.tab)
      }






#' Wrapper to Create a Biobase ExpressionSet
#' Create a Biobase \code{\link[Biobase]{ExpressionSet}} from a data matrix, fdata and pddta objects
#' @family eset
#' @family gene expression
#' @param data, a matrix containg expression values
#' @param pdata, data frame with phenotype data (pdata), i.e. sample information.
#' @param fdata,  data frame with feature data (fdata), i.e. reporter/probe information
#' @param sample.names.col, index (number) of what column from pdata to use as sampleNames
#' @param feature.names.col, index (number) of what column from fdata to use as featureNames
#' @param platform, what platform is used - the annotation slot in the ExpressionSet object.
#' @param data_processing, additional information on data processing.
#' @return an \code{\link[Biobase]{ExpressionSet}} object
#' @export
esetCreate<-function(data, pdata, fdata, sample.names.col=1, feature.names.col=1, platform=c(), data_processing=NULL){
	requireNamespace("Biobase")

	rownames(pdata)<-as.character(pdata[,sample.names.col])
	p<-new("AnnotatedDataFrame", data=pdata, varMetadata=data.frame(row.names=names(pdata), labelDescription=rep(NA, length(names(pdata)))))

	rownames(fdata)<-as.character(fdata[,feature.names.col])

	f<-new("AnnotatedDataFrame", data=fdata, varMetadata=data.frame(row.names=names(fdata), labelDescription=rep(NA, length(names(fdata)))))

	rownames(data)<-as.character(fdata[,feature.names.col])
	colnames(data)<-as.character(pdata[,sample.names.col])

	dp<-new("MIAME")
	if(!is.null(data_processing)){
		preproc(dp)<-as.list(data_processing)
				}
	data <- as.matrix(data)
	my.eset<-new("ExpressionSet", phenoData=p, featureData=f, exprs=data, experimentData=dp)
	#annotation(my.eset)<-platform

	print(show(my.eset))
	return(my.eset)
} ## end create.eset





#' Combine two or more ExpressionSet sets into one object
#' Pdata will be generated from annotations common to all esets (remainig annotations will be lost). Identical featureNames for all eSets are needed
#' @aliases combine.esets.dl
#' @family eset
#' @family gene expression
#' @family transformation
#' @param eset.names, specify the names of objects of class eset to be merged separated by comma
#' @return an object of class ExpressionSet
#' @export
esetCombine <- function(eset.names) {
    requireNamespace("Biobase")
    requireNamespace("methods")

    # check if duplicated objects
    u <- duplicated(eset.names)
    if (any(u != F)) {
        cat("\n ... ERROR: duplicated object(s)", eset.names[!u])
        return()
    }
    # check if objects exists
    u <- unlist(lapply(eset.names, function(x) {
        exists(x)
    }))
    if (any(u == F)) {
        cat("\n ... object", eset.names[!u], "not found in workspace")
        return()
    }

    # create a list of provided esets
    eset.list <- lapply(eset.names, get)
    u <- lapply(eset.list, class)
    if (any(unlist(u) != "ExpressionSet")) {
        cat("\n ... ERROR: object", eset.names[unlist(u) != "ExpressionSet"], "is not an ExpressionSet")
        return()
    }

    # create expr.data, pdata and fdata.
    featurenames <- lapply(eset.list, Biobase::featureNames)
    u <- lapply(featurenames[-1], function(x) identical(featurenames[[1]], x))
    if (any(u == F)) {
        cat("\n ... ERROR: not all esets have identical featureNames")
        return()
    }
    fdata <- Biobase::fData(eset.list[[1]])

    # expression.data
    my.data <- lapply(eset.list, Biobase::exprs)
    my.data <- do.call("cbind", my.data)

    # pdata : get columns present in all data sets
    col.fun <- function(x) {
        o <- paste("colnames(Biobase::pData(", x, "))", sep = "")
        return(eval(parse(text = o)))
    }
    col.list = sapply(eset.names, col.fun)
    if(class(col.list)=="matrix") col.list <- as.list(as.data.frame(col.list[,1:ncol(col.list)], stringsAsFactors=F))
    col.intersect = Reduce(intersect, col.list)
    # col.intersect = col.intersect[!col.intersect=='days_to_new_tumor_event_after_initial_treatment']
    pdata.list <- lapply(eset.list, Biobase::pData)
    pdata.reduce.fun <- function(x, col.intersect) {
        u <- match(col.intersect, colnames(x))
        x <- x[, u]
    }
    pdata.list <- lapply(pdata.list, pdata.reduce.fun, col.intersect = col.intersect)
    pdata <- do.call("rbind", pdata.list)
    sample.ids <- unlist(lapply(eset.list, Biobase::sampleNames))
    rownames(pdata) <- sample.ids
    colnames(my.data) <- sample.ids

    ## create eset
    p <- methods::new("AnnotatedDataFrame", data = pdata, varMetadata = data.frame(row.names = names(pdata), labelDescription = rep(NA, length(names(pdata)))))
    f <- methods::new("AnnotatedDataFrame", data = fdata, varMetadata = data.frame(row.names = names(fdata), labelDescription = rep(NA, length(names(fdata)))))

    dp <- methods::new("MIAME")
    # if(!is.null(data_processing)){ preproc(dp)<-as.list(data_processing) }
    my.eset <- methods::new("ExpressionSet", phenoData = p, featureData = f, exprs = my.data, experimentData = dp)
    # annotation(my.eset)<-platform

    print(methods::show(my.eset))
    return(my.eset)

}  # end function combine.esets.dl




#' Calculate sample group means for each row within an eset
#' Returns a matrix with mean values for each sample group and for each row (gene). Colums will be named after the annotation provide.
#' @family eset
#' @family gene expression
#' @family transformation
#' @param x, an object of class ExpressionSet
#' @param annot, a vector that defines sample groups
#' @param log2.data, if to log2 the data
#' @return a matrix with mean values for each sample group and for each row (gene). Colums will be named after the annotation provide and rows after featureNames of the eset
#' @export
esetGroupMeans <- function(x, annot, log2_data = F) {
    requireNamespace("Biobase")
   requireNamespace("BiocGenerics")
    if (class(x) == "ExpressionSet") {
        y <- Biobase::exprs(x)
    } else {
        y <- x
    }
    if (!log2_data) {
        print(paste("No additional Log2 transformation performed - check that data is logged"))
    }
    if (log2_data) {
        print(paste("Log2 transforming eset values"))
        y <- log(y, 2)
    }
    annot <- as.factor(as.character(annot))
    annot_i_list <- sapply(levels(annot), function(a) which(annot == a))

    if (class(annot_i_list) == "matrix")
        annot_i_list <- lapply(apply(annot_i_list, 2, as.list), unlist)

    mean_vec_list <- lapply(annot_i_list, function(z, d = y) apply(d[, z], 1, mean))
    group_mat <- sapply(mean_vec_list, BiocGenerics::rbind)
    rownames(group_mat) <- rownames(x)
    return(group_mat)
}



#' Center features/rows in an eset
#' Each row/feature in an eset will be re-rentered based on mean median or midpoint values. Can be performed on only a subset of samples.
#' @family eset
#' @family gene expression
#' @family transformation
#' @param eset, an object of class ExpressionSet
#' @param measure, can be 'mean','median', or 'max.min' (i.e., the midpoint between min and max values)
#' @param subset, specify sample names if measure should be calculated from only a subset of samples. Only these samples will get mean/median of 0.
#' @return an ExpressionSet where each row now is centered to the mean, median or the max.min
#' @export
esetCenterGenes <- function(eset, measure = c("mean", "median", "max.min"), subset = NULL) {
    requireNamespace("Biobase")
    measure <- match.arg(arg = measure, choices = c("mean", "median", "max.min"))
    cat("\n\n ... centering genes using ", measure, " method")
    if (class(eset) %in% c("LumiBatch", "ExpressionSet")) {
        if (is.null(subset)) {
            subset <- c(1:ncol(exprs(eset)))
        }
        if (measure %in% c("mean", "median"))
            apa <- esApply(eset, 1, function(x, m, s) x - eval(parse(text = m))(x[s], na.rm = T), m = measure, s = subset)
        if (measure %in% c("max.min"))
            apa <- esApply(eset, 1, function(x, s) x - c(max(x[s], na.rm = T) + min(x[s], na.rm = T))/2, s = subset)
        exprs(eset) <- t(apa)
        return(eset)
    } else {
        stop("object is not an eset!")
    }

}  # end center genes



# eset <- dlfoo2::pan_gdc_rppa

#' Variance filter an eset
#' Returns a matrix with mean values for each sample group and for each row (gene). Colums will be named after the annotation provide.
#' @family eset
#' @family gene expression
#' @family filter
#' @param eset, an object of class ExpressionSet
#' @param method, sd (standard deviation), IQR, or
#' @param cutoff, specify the cutoff (sd, number of genes) to be applied
#' @return a filtered eset where only features that meet the variance filter criteria are kept
#' @export
esetVarianceFilter <- function(eset, method = "sd", cutoff = 1) {
    requireNamespace("Biobase")
    requireNamespace("BiocGenerics")
    requireNamespace("stats")

    stopifnot(class(eset) %in% c("LumiBatch", "ExpressionSet"))
    stopifnot(method %in% c("sd", "IQR", "sd.rank"))
    match.arg(method, choices = c("sd", "IQR", "sd.rank"))

    if(length(which(is.na(exprs(eset))))>=1) message("\n  WARNING: NAs present in data matrix: ", length(which(is.na(exprs(eset)))), " out of ", length(exprs(eset))," values (",round(length(which(is.na(exprs(eset))))/length(exprs(eset))*100,2), "%)\n")

    if (method == "sd.rank") {
      o <- apply(Biobase::exprs(eset), 1, function(x) sd(x, na.rm = T))
      oo <- order(o, decreasing = T, na.last = F)
      eset <- eset[oo, ][1:cutoff]
      cat("\n ... ... Returning eset with top ", cutoff, " sd genes")
      return(eset)
    }

    cat("\n ...", method, " filter eSet using cutoff ", cutoff, sep = "")
    o <- apply(exprs(eset), 1, function(x) do.call(what = method, args = list(x, na.rm=T))) > cutoff

    if (length(which(o == T)) == nrow(exprs(eset))) {
        cat("\n ... ... All features passed threshold")
        return(eset)
    }
    if (length(which(o == F)) == nrow(exprs(eset))) {
        cat("\n ... ... No features passed threshold")
        return(eset[o, ])
    }
    cat("\n ... ... ", length(which(o == T)), " features remaining", sep = "")
    return(eset[o, ])
}  # end esetVarianceFilter





#' Merge (collapse) rows of an eset on replicate probes
#' Returns an eset where replicate probes are merged based on mean, median or highest sd probe
#' @family eset
#' @family gene expression
#' @family transformation
#' @family old functions
#' @param eset, an object of class ExpressionSet
#' @param fdata.column, name of fData column in eset that specify what reporters to merge on
#' @param method, "mean" or "highest.sd"
#' @return ExpressionSet for which replicate reporters are merged. fData columns are merged using pipe (|) when needed
#' @export
esetMergeRows <- function(eset = NULL, fdata.column = NULL, method='mean'){
   requireNamespace("Biobase")

   stopifnot(class(eset)=='ExpressionSet')
   fdata <- fData(eset)
   fdata.c <- colnames(fdata)
   fdata.n <- match.arg(arg = fdata.column, choices = fdata.c)
   fdata.i <- match(fdata.n, fdata.c)
   my.c <- colnames(fdata)[!fdata.c==fdata.n]


   # get unique ids
   ref.ids <- as.character(unique(fdata[,fdata.i]))
   ref.ids <- ref.ids[!is.na(ref.ids)]
   ref.table <- data.frame(row.names = ref.ids)
   ref.table[,1] <- ref.ids
   colnames(ref.table) <- fdata.n

   # data
   my.mat <- exprs(eset)

   index.list<-sapply(ref.ids, function(x, y=fdata[,fdata.i]) which(y==x))
   data.list <- lapply(index.list, function(x, my.data=my.mat) my.data[x,])

   if(method=='mean') {
      collapse.data.mean <- function(x){
       #if(class(x)=='numeric') return(x)
       #if(class(x)=='matrix') return(apply(x, 2, mean, na.rm=T))
         if(class(x)=='matrix') return(apply(x, 2, mean, na.rm=T))
         else(return(x))}
      cat('\n ... merged rows on mean vaue - NA values removed')
      data.list.merged <- lapply(data.list, collapse.data.mean)
   } # end if collapse mean

   if(method=='highest.sd') {
      collapse.data.sd <- function(x){
         if(class(x)=='matrix'){
            foo = apply(x,1,sd)
            if(all(foo==max(foo))) return(x[1,])
            if(!all(foo==max(foo))) return(x[which(foo==max(foo)),])
            } # end if x is matrix
         else(return(x))
         }
      data.list.merged <- lapply(data.list, collapse.data.sd)
      } # end if collapse highest sd

     if(method=='sum') {
      collapse.data.sum <- function(x){
         if(class(x)=='matrix') return(apply(x, 2, sum, na.rm=T))
         else(return(x))}
      cat('\n ... merged rows by sum - NA values removed')
      data.list.merged <- lapply(data.list, collapse.data.sum)
      } # end if collapse sum


   stopifnot(identical(names(data.list.merged), rownames(ref.table)))
   my.mat.merge <- matrix(nrow=nrow(ref.table), ncol=ncol(my.mat), data=unlist(data.list.merged), byrow=T)
   rownames(my.mat.merge) <- rownames(ref.table)
   colnames(my.mat.merge) <- colnames(my.mat)

   cat('\n ... ... collapsing fdata using pipe')
   for(i in 1:length(my.c)){
         col.i <- my.c[i]
         col.vals <- fdata[,col.i]
         col.vals.collapsed <- lapply(index.list, function(x, y=col.vals) paste(unique(col.vals[x]), collapse='|'))
         col.vals.collapsed <-  gsub(pattern = '[|]NA', '', col.vals.collapsed)
         col.vals.collapsed[col.vals.collapsed=='[|]'] <- ""
         col.vals.collapsed[col.vals.collapsed==""] <- NA
         substrRight <- function(x){
            if(is.na(x)) return(NA)
            if(substr(x, nchar(x), nchar(x)) == "|") return(substr(x, 1, nchar(x)-1))
            if(substr(x, 1, 1) == "|") return(substr(x, 2, nchar(x)))
            return(x)
            }
         col.vals.collapsed <- sapply(col.vals.collapsed, substrRight)
         ref.table[,(i+1)] <- col.vals.collapsed
         colnames(ref.table)[(i+1)] <- col.i
         }

   new.eset <- esetCreate(data = my.mat.merge, fdata = ref.table, pdata = pData(eset), feature.names.col = 1)
   sampleNames(new.eset) <- sampleNames(eset)
   return(new.eset)
   } # end function esetMergeRows




#' Convert epression data for eset to rank values
#' Converting epression data for eset to rank values using \code{\link[Biobase]{rank}} function in BiocGenerics
#' @family eset
#' @family gene expression
#' @family transformation
#' @family normalization
#' @param eset, an object of class ExpressionSet
#' @param ties.method, passed on to function \code{\link[Biobase]{rank}}
#' @return ExpressionSet for which replicate values (per sample) are concerted to ranks
#' @export
esetRankValues <- function(eset, ties.method = "min"){
   requireNamespace("BiocGenerics")
   cat('\n ... ... Converting epression data for eset to rank values')
   cat('\n ... ... .. ties method: min')
   cat('\n ... ... .. NAs are kept NA')
      rank.foo <- function(x){
         rank(x, ties.method = 'min', na.last=NA)
         }
   exprs(eset) <- apply(exprs(eset), 2,  rank.foo)
   return(eset)
}




#' Create a matrix based on input entrez_ids
#' Use an eset to extract available features. Extract features as defined by vector.
#' If not available in eset - then this row is set to NA
#' @family eset
#' @family gene expression
#' @family transformation
#' @param entrez, a vector containing entrez ids to generate matrix from.
#' @param eset, an ExpressionSet
#' @return a matrix where rows represent entrez ids as defined by input parameter.
#' @export
create.dummy.matrix <- function(entrez, eset){
         u <- grep("ENTREZ_ID", colnames(fData(eset)))
         if(length(u)==1) colnames(fData(eset))[u] <- "ENTREZID"
         my.mat <- matrix(nrow=length(entrez), ncol=ncol(exprs(eset)))
         colnames(my.mat) <- sampleNames(eset)
         rownames(my.mat) <- entrez

         entrez.d <- fData(eset)$ENTREZID[fData(eset)$ENTREZID %in% entrez]
         stopifnot(length(which(duplicated(entrez.d)))==0)
         u <- which(fData(eset)$ENTREZID %in% entrez.d)
         eset <- eset[u,]
         u <- match(fData(eset)$ENTREZID, rownames(my.mat))
         my.mat[u,] <- exprs(eset)
         return(my.mat)
      }





#' Switch from enembl to UCSC stype chromosome nomenclature
#' @family genomics
#' @family old functions
#' @family misc
#' @param chr numeric vector
#' @return character vector
#' @export
#'
toChr<-function(chr){
	if(!is.numeric(chr)){stop("input is not numeric")}
	chr<-paste("chr",chr,sep="")
	if(length(which(chr=="chr23"))>0){chr[which(chr=="chr23")]<-"chrX"}
	if(length(which(chr=="chr24"))>0){chr[which(chr=="chr24")]<-"chrX"}
	return(chr)
}


#' Switch from enembl to UCSC stype chromosome nomenclature
#' @family genomics
#' @family old functions
#' @family misc
#' @param chr object of class GenoSet
#' @return numeric vector
#' @export
#'
fromChr<-function(chr){
	#if(!is.numeric(chr)){stop("input is not numeric")}
	chr<-sub("chr","",chr)
	if(length(which(chr=="X"))>0){chr[which(chr=="X")]<-"23"}
	if(length(which(chr=="Y"))>0){chr[which(chr=="Y")]<-"24"}
	chr<-as.numeric(chr)
	return(chr)
	}


#' get full (accumulated) genomic position
#' @family old functions
#' @family genomics
#' @family misc
#' @param pos vector of genomic position(s)
#' @param chrom vector of chromosomes
#' @param force.full.genome.positions if to add virual (cumulative) genome positions
#' @param add.space if to add space in between chromosomes
#' @return gemomic locations
#' @export
genomicPosition <- function(pos, chrom, force.full.genome.positions=F, add.space=0){

   centromere.limit <- readRDS(file="/Users/david/R/R_FUNCTIONS_AND_REFERENCE_FILES/centromere_limits_grch38.rds")

	x.start<-as.integer(pos)
	x.chr<-as.character(chrom)
	x.chr<-fromChr(x.chr)

	include.chromosomes<-unique(x.chr)

	if(force.full.genome.positions==F){
		addPos<-function(chromosome, include.chromosomes=include.chromosomes){
			include.chromosomes<-sort(include.chromosomes)
			u<-match(chromosome, include.chromosomes)-1
			if(u==0){
				add.pos<-0
				}else{
					uu<-include.chromosomes[1:u]
					add.pos<-sum(centromere.limit$maxQarm[uu])+u*add.space
					}
				return(add.pos)
				} # end addpos function
			x.out<-x.start+sapply(x.chr, addPos, unique(x.chr))
			}

	if(force.full.genome.positions){
		addPos<-function(chromosome){
			if(chromosome==1){
				add.pos<-0
			}else{
				add.pos<-sum(centromere.limit$maxQarm[1:(chromosome-1)])+(chromosome-1)*add.space
				}
			} # end addpos function
			x.out<-x.start+sapply(x.chr, addPos)
		}
	return(x.out)
	}

