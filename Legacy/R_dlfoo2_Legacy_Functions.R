

#' enrichmentAnalysis_ClusterProfiler v1
#' @description wrapper that perform Enrichment analysis using ClusterProfiler and DOSE packages
#' @family enrichment
#' @family gene ontology
#' @family go
#' @param gene.list, directory created by ConsensusClusterPlus analysis containing results. File must be called "cc.results.Rdata".
#' @param run.folder, name of higher level directory
#' @param run.name, name of run (directory created)
#' @param term.names, what mSig or other signatues to load
#' @param min.genes.in.term, minimum included genes in sigificant term
#' @return write results to utput directory: EnricherResults.txt: tab delimited text file with all significant results from all terms and all gene lists. EnricherResults.rds: same as .txt but saved as .rds. EnricherResults_comparecluster.rds. ResultsSummary.txt: Results Summary on how many significant terms per gene list were obtained.
#' @export
enrichmentAnalysis_ClusterProfiler <- function(
   gene.list = NULL,
   run.folder = NULL,
   run.name = NULL,
   term.names =  c("c2","c3.tft","c6","c7","c5.cc","c5.bp","c5.mf"),
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
   enrich.list = vector(mode="list", length=length(term.names))
   names(enrich.list) = term.names
   enrich.cmp = enrich.list


   require(clusterProfiler)
   require(DOSE)

   # Load signatures (mSig etc)
   message(" ... loading mSig gene signatures...")
   if(!exists("ch")) ch <- read.gmt("~/RESOURCES/MSigDB/msigdb_v6.0_GMTs/h.all.v6.0.symbols.gmt")
   if(!exists("c1")) c1 <- read.gmt("~/RESOURCES/MSigDB/msigdb_v6.0_GMTs/c1.all.v6.0.symbols.gmt")
   if(!exists("c2.bc")) c2.bc <- read.gmt("~/RESOURCES/MSigDB/msigdb_v6.0_GMTs/c2.cp.biocarta.v6.0.symbols.gmt")
   if(!exists("c2.kegg")) c2.kegg <- read.gmt("~/RESOURCES/MSigDB/msigdb_v6.0_GMTs/c2.cp.kegg.v6.0.symbols.gmt")
   if(!exists("c2.cgp")) c2.cgp <- read.gmt("~/RESOURCES/MSigDB/msigdb_v6.0_GMTs/c2.cgp.v6.0.symbols.gmt")
   if(!exists("c2.rea") )c2.rea <- read.gmt("~/RESOURCES/MSigDB/msigdb_v6.0_GMTs/c2.cp.reactome.v6.0.symbols.gmt")
   if(!exists("c3.tft")) c3.tft <- read.gmt("~/RESOURCES/MSigDB/msigdb_v6.0_GMTs/c3.tft.v6.0.symbols.gmt")
   if(!exists("c6")) c6 <- read.gmt("~/RESOURCES/MSigDB/msigdb_v6.0_GMTs/c6.all.v6.0.symbols.gmt")
   if(!exists("c7")) c7 <- read.gmt("~/RESOURCES/MSigDB/msigdb_v6.0_GMTs/c7.all.v6.0.symbols.gmt")
   if(!exists("c5.cc")) c5.cc <- read.gmt("~/RESOURCES/MSigDB/msigdb_v6.0_GMTs/c5.cc.v6.0.symbols.gmt")
   if(!exists("c5.bp")) c5.bp <- read.gmt("~/RESOURCES/MSigDB/msigdb_v6.0_GMTs/c5.bp.v6.0.symbols.gmt")
   if(!exists("c5.mf")) c5.mf <- read.gmt("~/RESOURCES/MSigDB/msigdb_v6.0_GMTs/c5.mf.v6.0.symbols.gmt")
   c2 <- do.call("rbind", list(c2.bc, c2.kegg, c2.rea, ch))

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
         res@compareClusterResult$Description[u] <- paste(substr(res@compareClusterResult$Description[u],1,40),"...")
         res@compareClusterResult <- res@compareClusterResult[res@compareClusterResult$Count>2,]
         rownames(res@compareClusterResult) <- 1:nrow(res@compareClusterResult)
         return(res)
         }


   ## Loop all msig terms
   for(i in 1:length(term.names)){
      term.name <- term.names[i]
      message("... ... performing EA using mSig term: ", term.name)
      my.term <- eval(parse(text=term.name))
      if(!any(unlist(lapply(gene.list, function(x){sapply(x, function(y) any(unlist(y %in% my.term$gene)))})))) message("no genes from list present in ", term.name)
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



#' OLD function - but the featureIdObject is used in som plot functions
#' Use function \code{\link[dlfoo2]{featureGetBM}} instead
#' Use one gene/feature id to create object with multiple feature annotations, e.g. symbol, entrez, ENSG, and for multiple species
#' The featureIdObject can be used in downstream analyses
#' @family feature annotation
#' @family old functions
#' @param id A gene identifier (default is hgnc_gene symbol)
#' @param idType What type of identifer: c("hgnc_symbol","entrezgene_id","ensembl_gene_id")
#' @param Species What species (hs, mm, or rn)
#' @return object list with slots
#' @export
featureIdObject = function(id, idType=c("hgnc_symbol","entrezgene","ensembl_gene_id"), Species=c("hs","mm","rn")){
   idType = match.arg(idType)
   Species = match.arg(Species)
   require(biomaRt)
   featureId <- list(id=id,  idType=idType, species=Species,hs=NULL, mm=NULL, rn=NULL, geneLoc=NULL)

   message("Fetching gene identifiers for ",id," using bomaRt")
   mart_hs = useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
   mart_mm = useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl")
   mart_rn = useMart("ENSEMBL_MART_ENSEMBL", dataset ="rnorvegicus_gene_ensembl")
   mart.attributes.hs = c("hgnc_symbol","entrezgene","ensembl_gene_id")
   mart.attributes.mm = c("mgi_symbol","entrezgene","ensembl_gene_id")
   mart.attributes.rn = c("rgd_symbol","entrezgene","ensembl_gene_id")

   if(featureId$species=="hs"){
      message("... homo sapiens")

      hsBM = getBM(attributes=mart.attributes.hs, values=featureId$id, filter=featureId$idType, mart=mart_hs)
      if(nrow(hsBM)<1){stop("Cannot find ",id," in Ensembl database. Please supply a valid hgnc_symbol")}
      if(nrow(hsBM)>1){
         mult_idTypes = colnames(hsBM)[apply(hsBM, 2, function(x) length(unique(x))>1)]
         message(" ... note! mutiple matches (1:n) found for:", paste(mult_idTypes))
         message(" ... ... coercing these with pipes '|'")
         hsBM = apply(hsBM, 2, function(y){paste(unique(y), collapse = "|")})
         hsBM = data.frame(hgnc_symbol=gsub(" ","", hsBM[[1]]), entrezgene=gsub(" ","", hsBM[[2]]), ensembl_gene_id=gsub(" ","", hsBM[[3]]), stringsAsFactors = F)
         }
      featureId$hs = hsBM

      ## Get homologs for Mouse
      message("... obtaining mouse homolog(s)")
      mm_homolog = getBM(attributes="mmusculus_homolog_ensembl_gene", values=featureId$id, filter=featureId$idType, mart=mart_hs)
      if(nrow(mm_homolog)<1){
         message("... ... no ensembl mouse homolog found")
         #getBM(attributes="mmusculus_homolog_ensembl_gene", values=featureId$hs$ensembl_gene_id, filter="ensembl_gene_id", mart=mart_hs)
         featureId$mm <- data.frame(mgi_symbol=NA,entrezgene=NA,ensembl_gene_id=NA)
         }
      if(nrow(mm_homolog)>1){
         message("... ... multiple mouse homologs identified (1:n): ",mm_homolog)
      }
      if(nrow(mm_homolog)>0){
         featureId$mm = getBM(attributes=mart.attributes.mm, values=as.character(mm_homolog$mmusculus_homolog_ensembl_gene), filter="ensembl_gene_id", mart=mart_mm)
         print(featureId$mm)
         }

      ## Get homologs for Rat
      message("... obtaining rat homolog(s)")
      rn_homolog = getBM(attributes="rnorvegicus_homolog_ensembl_gene", values=featureId$id, filter=featureId$idType, mart=mart_hs)
      if(nrow(rn_homolog)<1){
         message("... ... no ensembl rat homolog found")
         #getBM(attributes="mmusculus_homolog_ensembl_gene", values=featureId$hs$ensembl_gene_id, filter="ensembl_gene_id", mart=mart_hs)
         featureId$rn <- data.frame(rgd_symbol=NA,entrezgene=NA,ensembl_gene_id=NA)
         }
      if(nrow(rn_homolog)>1){
         message("... ... multiple rat homologs identified (1:n): ", rn_homolog)
      }
      if(nrow(rn_homolog)>0){
         featureId$rn = getBM(attributes=mart.attributes.rn, values=as.character(rn_homolog$rnorvegicus_homolog_ensembl_gene), filter="ensembl_gene_id", mart=mart_rn)
         print(featureId$rn)
      }
   } # end if hs
   return(featureId)
}


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
   my_ontologies = c("h.all","c2.cp.biocarta","c2.cp.kegg","c2.cp.reactome","c3.tft","c5.bp","c5.cc","c5.mf"),
   entrez.or.symbol = c("symbol"),
   minSize=10, maxSize = 500, nperm=1000
   ){
   require(tibble)
   require(tidyr)
   require(fgsea)
   my_ranks <- rank.list

   # myDir <- file.path(runFolder, runName)
   myDir <- run.folder
   if(!dir.exists(myDir)) dir.create(myDir)
   myDir <- file.path(myDir,run.name)
   if(dir.exists(myDir)) stop("run already exists - delete old run folder")
   dir.create(myDir)


   # Run fgsea for all selected tissues and selected ontology subsets
   # Clean upp resutls to tab format and keep padj<0.05

   message("... start")
   # pathway_files <- list.files("~/RESOURCES/MSigDB/msigdb_v6.0_GMTs/", pattern="entrez_pathwayList.rds", full.names = T)
   if(entrez.or.symbol == "entrez") pathway_files <- list.files("~/RESOURCES/MSigDB/msigdb_v7.0_GMTs/", pattern="entrez_pathwayList.rds", full.names = T)
   if(entrez.or.symbol == "symbol") pathway_files <- list.files("~/RESOURCES/MSigDB/msigdb_v7.0_GMTs/", pattern="symbols_pathwayList.rds", full.names = T)

   pathway_files <- pathway_files[apply(sapply(my_ontologies, function(x) grepl(x, pathway_files)),1, any)]

   fgsea_res_list <- lapply(pathway_files, function(x){
      my_pathway <- readRDS(file = x)
      message("... running fgsea pathway ", x)
      lapply(my_ranks, function(y) {
         my_f <- fgsea(y, pathways = my_pathway, minSize=10, maxSize = 500, nperm=1000)
         return(my_f)
         #return(my_f[padj<0.05,])
         })
      })
   # str(fgsea_res_list)
   names(fgsea_res_list) <- gsub("[.]v7.*","",gsub(".*[//]","",pathway_files))

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
         yy <- do.call("rbind", y)
         yy$ontology <- x
         return(yy)
         })
      #str(apa2)
   u <- lapply(apa2, function(x) length(x))
   apa2 <- apa2[(unlist(u)>1)]
   apa2 <- as_tibble(do.call("rbind", apa2))
   #str(apa2)
   #apa2[,1]
   leading_edge <- unlist(lapply(apa2[,"leadingEdge"]$leadingEdge, function(x) paste0(x, collapse="|")))
   #str(leading_edge)
   u<-which(colnames(apa2)=="leadingEdge")
   #str(apa2[,-u])

   write.table(cbind(apa2[,-u],leading_edge), file=paste0(myDir,"/",run.name,"_fgsea_Ontologies_all.txt"), sep="\t", quote=F, row.names=F)
   saveRDS(fgsea_res_list, file=paste0(myDir,"/",run.name,"_fgsea_allOntologies.rds"))

   Hmisc::describe(apa2$padj)
   uu <- apa2$padj<0.05
   table(uu)
   write.table(cbind(apa2[uu,-u],leading_edge[uu]), file=paste0(myDir,"/",run.name,"_fgsea_Ontologies_significant.txt"), sep="\t", quote=F, row.names=F)




   #
   #    ## create filtered small list with significant terms
   #    ## :::::::::::::::
   #    fgsea_res_list_small <- lapply(pathway_files, function(x){
   #       my_pathway <- readRDS(file = x)
   #       message("... running fgsea pathway (ledingEdge): ", x)
   #       lapply(my_ranks, function(y) {
   #          my_f <- fgsea(y, pathways = my_pathway, minSize=10, maxSize = 500, nperm=1000)
   #          return(as.data.frame(my_f[padj<0.05,]))
   #          })
   #       })
   #    names(fgsea_res_list_small) <- gsub("[.]v7.*","",gsub(".*[//]","",pathway_files))
   #
   #    #str(fgsea_res_list_small[[1]])
   #    #nrow(fgsea_res_list_small[[1]][[1]])
   #
   #    apa <- fgsea_res_list_small
   #    apa <- lapply(apa, function(x){
   #       y <- unlist(lapply(x, nrow))
   #       x[!(y==0)]
   #       })
   #    #str(apa)
   #
   #    apa <- lapply(apa, function(x){
   #       y <- names(x)
   #       lapply(y, function(y, z=x){
   #          zz <-z[[y]]
   #          zz$tissue <- y
   #          return(zz)
   #       })
   #       })
   #    #str(apa)
   #
   #
   #       apa2 <- lapply(names(apa), function(x){
   #       y <- apa[[x]]
   #       yy <- do.call("rbind", y)
   #       yy$ontology <- x
   #       return(yy)
   #       })
   #    #str(apa2)
   #    u <- lapply(apa2, function(x) length(x))
   #    apa2 <- apa2[(unlist(u)>1)]
   #
   #    apa2 <- as_tibble(do.call("rbind", apa2))
   #    apa2 # A tibble: 12,129 x 9
   #
   #    leading_edge <- unlist(lapply(apa2[,"leadingEdge"]$leadingEdge, function(x) paste0(x, collapse="|")))
   #    #str(leading_edge)
   #    u<-which(colnames(apa2)=="leadingEdge")
   #    write.table(cbind(apa2[,-u], leading_edge), file=paste0(myDir,"/",run.name,"_fgsea_Ontologies_significant.txt"), sep="\t", quote=F, row.names=F)

   }





#' OLD Function to start Shiny app from input tSNE results list.
#' @family tsne
#' @family shiny
#' @family PlotWrap
#' @family OLD
#' @param tSNEres results list obect from tsneWrapper analyeis (as also used in e.g. dlfoo2::tsneMe wrappers)
#' @return a shiny app plot will open
#' @export
tsneShinyPlot <- function(tSNEresList){
   require(dlfoo2)
   require(shiny)
   require(dplyr)
   require(ggplot2)
   require(survival)

   data("gdc_cn", package="dlfoo2")
   #data("tcga_surv_tab", package="dlfoo2")
   #data("pdata_utc_nt", package="dlfoo2")
   data("shape_subtypes", package="dlfoo2")
   data("color_subtypes", package="dlfoo2")
   pdata <- pdata_utc_nt
   tsne_df <- do.call("rbind", lapply(tSNEres, tsne2df))


   annot_columns <- c("project","sample_type","sample_type2","gender",
         "vital_status","person_neoplasm_cancer_status",
         "included.in.initial.TCGA.marker.paper","explicitly.removed.from.marker.paper",
         "EOS_davis_2014","SARC_davis_2014",
         "FGFR2_IIIb_Zhao","uro","histology_davis_2014",
         "taxonomy_combo","taxonomy_chen",
         "taxonomy_uro","VHL_silenced","TFE3.TFEB.fusion",
         "SWI.SNF.module","PI3K.module","TP53","VHL",
         "SETD2","BAP1","PBRM1","KDM5C",
         "NF2","PTEN","ARID1A","MICALCL",
         "STAG2","SLC1A3","CDKN1A","MTOR",
         "MET","SMARCB1","TCEB1","NFE2L2",
         "PIK3CA","MLL3","cancer.type",
         "taxonomy_dl","taxonomy_simplified","taxonomy_limma",
         "taxonomy_dl2","atacSeq","cc_1435set_k25","cc_882_set_k15",
         "cc_14_set_k2")





   tsne_annotations <- pdata[,c("sample_id", "patient_barcode", annot_columns)]
   tsne_df <- do.call("rbind", lapply(tSNEres, tsne2df))
   tsne_df <- dplyr::left_join(tsne_df, tsne_annotations)
   tsne_df <- reshape2::melt(tsne_df, measure.vars=annot_columns, variable.name="annotation", value.name="annotation_value")
   str(tsne_df)
   tsne_df <- dplyr::left_join(tsne_df, pdata[,c("sample_id","taxonomy_dl")])



   # table(tsne_df$annotation_value)
   # tsne_df$taxonomy_dl2 <- droplevels(tsne_df$taxonomy_dl2)

   tsne_df_surv <- tsne_df[-grep("-NT",tsne_df$annotation_value),]


   available.annotations <- as.character(sort(unique(tsne_df$annotation)))
   available.perplexes <- sort(unique(tsne_df$perplexity))

   tsnePlotWrap <- function(my_df, my_perplex, my_annot, pointsize, Alpha){
         my_df$annotation_value = as.factor(my_df$annotation_value)
         shape_key <- shape_subtypes[levels(my_df$annotation_value)]

         if(my_annot %in% c("cancer.type","project","sample_type","sample_type2") | grepl("taxonomy",my_annot)){
            g <-  ggplot(my_df,
                  aes(x=dim1, y=dim2)) +
                  ggtitle(paste0("tSNE-plot, perplexity=", my_perplex), label = my_annot) +
                  guides(colour = guide_legend(override.aes = list(size=7))) +
                  guides(fill=guide_legend(title=my_annot)) +
                  labs(x="tSNE dim1", y="tSNE dim2", colour='') +
                  theme(axis.title.y = element_text(size=24)) +
                  theme(axis.title.x = element_text(size=24)) +
                  theme(legend.text = element_text(size=18)) +
                  theme(title = element_text(size=24))

            g <- g +
                  geom_point(aes(fill=annotation_value, shape=annotation_value), colour="gray20",stroke=0.25, size=pointsize, alpha=Alpha) +
                  scale_fill_manual(values=color_subtypes, na.value="black") +
                  scale_shape_manual(values = as.numeric(shape_key), na.value = 4)
         }else{
            my_pal <- c(dlfoo2::palette.dl.annot[[19]], rainbow(15))[sample(length(levels(my_df$annotation_value)), replace=F)]
            names(my_pal) <- levels(my_df$annotation_value)
              g <- g +
                  geom_point(aes(fill=annotation_value), shape=21, stroke=0.25, size=pointsize, alpha=Alpha) +
                  scale_fill_manual(values=my_pal, na.value="white")


            }
         return(g)
      }

   survPlotWrap <- function(y){
                  z <- as.data.frame(tcga_surv_tab)[tcga_surv_tab$patient_barcode %in% y,]
                  z$fu_group <- "selected tumors"
                  dlfoo2::tcgaSurvPlot(z)
                  #surv.obj<-Surv(z$fu_time, z$fu_status)
                  #survfit.obj<-survfit(surv.obj~z$fu_group)
                  #plot(survfit.obj)
                  }


   #  Shiny App
   # -------------------------------
   shinyApp(
      # UI section
      ## ---------
      ui <- fluidPage(

         tags$hr(),
         h1('tSNE plotter'),

         fluidRow(
            column(3,
               uiOutput("annotation_selection")
               ),

            column(3,
               uiOutput("perplexity_selection")
               ),

            column(3,
               sliderInput(inputId = "pointSize", label = "Select point size",
                  min=0.5, max=10, ticks = T, step = 0.5, value = 5)
             ),
            column(3,
               sliderInput(inputId = "pointAlpha", label = "Select alpha",
                  min=0, max=1, ticks = T, step = 0.1, value = 0.7)
             )
            ),

         tags$hr(),
         fluidPage(
               title = "tSNE Plot",
               plotOutput('tsne_plot', height=900, width=1200, brush = "plot_brush")
               ),

         tags$hr(),
         fluidRow(p(class = 'text-center', downloadButton(outputId = 'table_download', label = 'Download selected data points'))),

         tags$hr(),
         fluidPage(
            title="Download",
            DT::dataTableOutput('plot_table_selection')),


         ## Genome plots
         tags$hr(),
         fluidPage(
            uiOutput("SelectGenomeSample")
         ),

         fluidPage(
           # textOutput('bajs'),
               plotOutput('genomePlot')

            ),


         tags$hr(),
         h1("Survival plot - selected samples"),
         fluidPage(
            plotOutput('survPlot')
         )

         ), # end ui fluid page


         # Define server logic
         server = function(input, output) {

            output$annotation_selection <- renderUI({
                    selectInput(inputId = "annot", label = "Available annotations:",
                           choices=available.annotations, selected = "taxonomy_dl")
                  })

            output$perplexity_selection <- renderUI({
                    selectInput(inputId = "perplex", label = "Available perplexities:",
                           choices=available.perplexes, selected = available.perplexes[1])
                  })

            output$tsne_plot <- renderPlot(
               tsnePlotWrap(my_df = subset(tsne_df, perplexity==input$perplex & annotation==input$annot), my_perplex = input$perplex, my_annot = input$annot, pointsize = input$pointSize, Alpha = input$pointAlpha)
               )

            output$plot_table_selection <- DT::renderDataTable(
               brushedPoints(subset(tsne_df, perplexity==input$perplex & annotation==input$annot), input$plot_brush), server = FALSE, selection = 'none'
            )
            # download the filtered data
            output$table_download = downloadHandler(filename = 'selected_points.csv', content = function(file) {
                  write.table(brushedPoints(subset(tsne_df, perplexity==input$perplex & annotation==input$annot), input$plot_brush), file, sep=";",row.names = F)
               })





      ## CN DATA
      ## ----------------
         ## Select sample for CN plot (dropdown list)
         observe({
            x = brushedPoints(subset(tsne_df,  perplexity==input$perplex & annotation==input$annot), input$plot_brush)$sample_id
            if(length(x)){
            output$SelectGenomeSample <- renderUI({
               selectInput(inputId = "genomeSample", label = "Genome plot",
                     choices=x, selected = x[1])
               })}
         })
         ## Render CN plot using dlfoo2::segmentPlotter
         observe({
            x = input$genomeSample
            if(length(x)){
               #y <- renderText(x)
               if(x %in% names(gdc_cn)){
                  output$genomePlot <- renderPlot(dlfoo2::segmentPlotter(gdc_cn[[x]], sample_name = x))
               }else{
                  output$genomePlot <- renderPlot(dlfoo2::ggPlotEmpty("Copy number profile not available"))
               }
               }
         })

      ## Clinical data (SurvPlot for selected samples)
      observe({
         xx=brushedPoints(subset(tsne_df,  perplexity==input$perplex & annotation==input$annot), input$plot_brush)$patient_barcode
         if(length(xx)){
            output$survPlot <- renderPlot(survPlotWrap(xx)

               )

            }
            })


      }
   ) # end shiny app
} # end function shiny tsne




#' OLD function. Start Shiny app from Input list of tSNE results - plot taxonomy etc
#' @family tsne
#' @family shiny
#' @family PlotWrap
#' @param tSNEresList List containing >1 results list obects from tsneWrapper analyeis (as also used in e.g. dlfoo2::tsneMe wrappers)
#' @return a shiny app plot will open
#' @export
tsneShinyPlotMultiple <- function(tSNEres){
   require(dlfoo2)
   require(shiny)
   require(dplyr)
   require(ggplot2)
   require(survival)



   data("gdc_cn", package="dlfoo2")
   data("tcga_surv_tab", package="dlfoo2")
   data("pdata_utc_nt", package="dlfoo2")
   data("shape_subtypes", package="dlfoo2")
   data("color_subtypes", package="dlfoo2")
   data("color_annotations", package="dlfoo2")

   # define what columns to use from pdata as annotation columns
   annot_columns <- c("project","sample_type","sample_type2","gender",
         "vital_status","person_neoplasm_cancer_status",
         "included.in.initial.TCGA.marker.paper","explicitly.removed.from.marker.paper",
         "EOS_davis_2014","SARC_davis_2014",
         "FGFR2_IIIb_Zhao","uro","histology_davis_2014",
         "taxonomy_combo","taxonomy_chen",
         "taxonomy_uro","VHL_silenced","TFE3.TFEB.fusion",
         "SWI.SNF.module","PI3K.module","TP53","VHL",
         "SETD2","BAP1","PBRM1","KDM5C",
         "NF2","PTEN","ARID1A","MICALCL",
         "STAG2","SLC1A3","CDKN1A","MTOR",
         "MET","SMARCB1","TCEB1","NFE2L2",
         "PIK3CA","MLL3","cancer.type",
         "taxonomy_dl","taxonomy_simplified","taxonomy_limma",
         "taxonomy_dl2","atacSeq","cc_1435set_k25","cc_882_set_k15",
         "cc_14_set_k2")

   pdata <- pdata_utc_nt
   tsne_annotations <- pdata[,c("sample_id", "patient_barcode", annot_columns)]

   tsne_list <- lapply(tSNEresList, function(x) do.call("rbind", lapply(x, dlfoo2::tsne2df)))


   available.datasets <- as.character(names(tsne_list))


   ## Select what pdata column to annotate on
   # annot_list <- lapply(annot_columns, function(x){
   #    pdata_annot <- pdata[,c("sample_id", "patient_barcode",x)]
   #    colnames(pdata_annot) <- c("sample_id", "patient_barcode","annotation")
   #    if(!is.factor(pdata_annot$annotation)) pdata_annot$annotation <- as.factor(pdata_annot$annotation)
   #    })


   #tsne_df <- do.call("rbind", lapply(tSNEres, tsne2df))
   tsne_list <- lapply(tsne_list, function(x){
      x <- dplyr::left_join(x, tsne_annotations)
      x <- reshape2::melt(x, measure.vars=annot_columns, variable.name="annotation", value.name="annotation_value")
      x <- dplyr::left_join(x, pdata[,c("sample_id","taxonomy_dl")])
      return(x)
   })

   # tsne_df_surv <- tsne_df[-grep("-NT",tsne_df$annotation_value),]
   # available.annotations <- as.character(sort(unique(tsne_df$annotation)))


   tsnePlotWrap <- function(my_df, my_perplex, my_annot, pointsize, Alpha, showGuide=T){
         my_df$annotation_value = as.factor(my_df$annotation_value)
         shape_key <- shape_subtypes[levels(my_df$annotation_value)]
         my_df <- my_df[rev(order(my_df$annotation_value)),]

         g <-  ggplot(my_df,
                  aes(x=dim1, y=dim2)) +
                  ggtitle(paste0("tSNE-plot, perplexity=", my_perplex), label = my_annot) +
                  guides(color = guide_legend(override.aes = list(size=7))) +
                  guides(fill=guide_legend(title=my_annot)) +
                  labs(x="tSNE dim1", y="tSNE dim2", colour='') +
                  theme(axis.title.y = element_text(size=24)) +
                  theme(axis.title.x = element_text(size=24)) +
                  theme(legend.text = element_text(size=18)) +
                  theme(title = element_text(size=24))
         if(my_annot %in% c("cancer.type","project","sample_type","sample_type2") | grepl("taxonomy",my_annot)){
            g <- g +
                  geom_point(aes(fill=annotation_value, shape=annotation_value), colour="gray20",stroke=0.25, size=pointsize, alpha=Alpha) +
                  scale_fill_manual(values=color_subtypes, na.value="black") +
                  scale_shape_manual(values = as.numeric(shape_key), na.value = 4)
            }else{
              g <- g +
                  geom_point(aes(fill=annotation_value), shape=21, stroke=0.25, size=pointsize, alpha=Alpha) +
                  scale_fill_manual(values=color_annotations, na.value=NA)
            }
         if(showGuide==F){
            g <- g + guides(fill=FALSE, shape=FALSE, color=FALSE)

         }
         return(g)
      }

   tsnePlotWrapSampleHighlight <- function(my_df, my_samples, my_perplex, pointsize, Alpha, showGuide=T){
         my_df$annotation_value = NA
         my_df$annotation_value[my_df$sample_id %in% my_samples] <- "selected"
         my_df$annotation_value <- as.factor(my_df$annotation_value)
         my_df <- my_df[rev(order(my_df$annotation_value)),]

         g <-  ggplot(my_df,
               aes(x=dim1, y=dim2)) +
               ggtitle(paste0("tSNE-plot, perplexity=", my_perplex), label = "") +
               guides(color = guide_legend(override.aes = list(size=7))) +
               guides(fill=guide_legend(title="selected samples")) +
               labs(x="tSNE dim1", y="tSNE dim2", colour='') +
               theme(axis.title.y = element_text(size=24)) +
               theme(axis.title.x = element_text(size=24)) +
               theme(legend.text = element_text(size=18)) +
               theme(title = element_text(size=24)) +
               geom_point(aes(colour=annotation_value, shape=annotation_value), stroke=2, size=pointsize, alpha=Alpha) +
               scale_color_manual(values="tomato4", na.value="gray75") +
               scale_shape_manual(values = 4, na.value = 21) +
               guides(fill=FALSE, shape=FALSE, color=FALSE)
         return(g)
      }




   survPlotWrap <- function(y){
                  z <- as.data.frame(tcga_surv_tab)[tcga_surv_tab$patient_barcode %in% y,]
                  z$fu_group <- "selected tumors"
                  dlfoo2::tcgaSurvPlot(z)
                  #surv.obj<-Surv(z$fu_time, z$fu_status)
                  #survfit.obj<-survfit(surv.obj~z$fu_group)
                  #plot(survfit.obj)
                  }


   #  Shiny App
   # -------------------------------
   shinyApp(
      # UI section
      ## ---------
      ui <- fluidPage(

         tags$hr(),
         h1('tSNE plotter'),

         ## 1: select dataset and perplexes

         fluidRow(
            column(3,
               selectInput(inputId = "dataToken", label = "Choose main tsne plot:",
               selected = names(tsne_list)[1],
               choices = names(tsne_list))),
               #uiOutput("dataset_selection")),
            column(2,
               uiOutput("perplexity_main")),
            column(3,
               uiOutput("annotation_selection")),
            column(2,
               sliderInput(inputId = "pointSize", label = "Select point size",
                  min=0.5, max=10, ticks = T, step = 0.5, value = 5)),
            column(2,
               sliderInput(inputId = "pointAlpha", label = "Select alpha",
                  min=0, max=1, ticks = T, step = 0.1, value = 0.7))
            ),


         tags$hr(),
         h2(textOutput(outputId = 'dataTokenText')),

         fluidPage(
               title = "tSNE Plot",
               plotOutput('tsne_plot_tax', height=900, width=1200)
               ),

         tags$hr(),
         fluidPage(
               title = "tSNE Plot selected annotation",
               plotOutput('tsne_plot_annot', height=900, width=1200, brush = "plot_brush")
               ),

         tags$hr(),
         fluidRow(p(class = 'text-center', downloadButton(outputId = 'table_download', label = 'Download selected data points'))),

         tags$hr(),
         fluidPage(
            title="Download",
            DT::dataTableOutput('plot_table_selection')),

         # Second data set
         tags$hr(),
         h2(textOutput(outputId = 'dataTokenText2')),

         fluidRow(
            column(4,
               uiOutput("perplexity_selection2")
               )),
         tags$hr(),

         fluidRow(
            column(6,
               title = "tSNE Plot selected annotation",
               plotOutput('tsne_plot_annot2', height=600, width=500)
               ),
            column(6 ,
               title = "tSNE Plot selected samples",
               plotOutput('tsne_plot_samples2', height=600, width=500)
               )),

         ## Third data set
         tags$hr(),
         h2(textOutput(outputId = 'dataTokenText3')),
         fluidRow(
            column(4,
               uiOutput("perplexity_selection3")
            )),
         fluidRow(
            column(6,
               title = "tSNE Plot selected annotation",
               plotOutput('tsne_plot_annot3', height=600, width=500)
               ),
            column(6,
               title = "tSNE Plot selected samples",
               plotOutput('tsne_plot_samples3', height=600, width=500)
               )),



         # Genome plots
         tags$hr(),
         fluidPage(
            uiOutput("SelectGenomeSample")
         ),

         fluidPage(
          # textOutput(renderText(output$dataToken3)),
              plotOutput('genomePlot')
           ),


         tags$hr(),
         fluidRow(
            column(4,textOutput('tsne_plot_indiv1_name')),
            column(4,textOutput('tsne_plot_indiv2_name')),
            column(4,textOutput('tsne_plot_indiv3_name'))),



         fluidRow(
            column(4,
               title = h1(textOutput('tsne_plot_indiv1_name')),
               plotOutput('tsne_plot_tax1', height=500, width=400)
               ),
            column(4 ,
               title = textOutput('tsne_plot_indiv2_name'),
               plotOutput('tsne_plot_tax2', height=500, width=400)
               ),
            column(4 ,
               title = textOutput('tsne_plot_indiv3_name'),
               plotOutput('tsne_plot_tax3', height=500, width=400)
            )),

         fluidRow(
            column(4,
               title = textOutput('tsne_plot_indiv1_name'),
               plotOutput('tsne_plot_indiv1', height=500, width=400)
               ),
            column(4 ,
               title = textOutput('tsne_plot_indiv2_name'),
               plotOutput('tsne_plot_indiv2', height=500, width=400)
               ),
            column(4 ,
               title = textOutput('tsne_plot_indiv3_name'),
               plotOutput('tsne_plot_indiv3', height=500, width=400)

            ))


         #
         # tags$hr(),
         # h1("Survival plot - selected samples"),
         # fluidPage(
         #    plotOutput('survPlot')
         # )

         ), # end ui fluid page


         # Define server logic
         server = function(input, output) {

            # Select main dataset - output$mainDataNo
            #output$dataset_selection <- renderUI({
            #   selectInput(inputId = "dataName", label = "Available tSNEs", choices=available.annotations, selected=available.annotations[1])
            #})


            # When dataset is changed - create main data table - tsne_df
            observe({
               x = input$dataToken

               output$dataTokenText <- renderText(x)
               available.annotations <- as.character(sort(unique(tsne_list[[input$dataToken]]$annotation)))

               output$annotation_selection <- renderUI({
                    selectInput(inputId = "annot", label = "Available annotations:",
                           choices=available.annotations, selected = "taxonomy_dl")
                  })

               available.perplexes <- sort(unique(tsne_list[[input$dataToken]]$perplexity))
               output$perplexity_main <- renderUI({
                    selectInput(inputId = "perplex", label = "Available perplexities:",
                           choices=available.perplexes, selected = available.perplexes[3])
                  })

               output$tsne_plot_tax <- renderPlot(
                  tsnePlotWrap(my_df = subset(tsne_list[[input$dataToken]], perplexity==input$perplex & annotation=="taxonomy_dl"),
                     my_perplex = input$perplex, my_annot = "taxonomy_dl", pointsize = input$pointSize, Alpha = input$pointAlpha)
               )

               output$tsne_plot_annot <- renderPlot(
                  tsnePlotWrap(my_df = subset(tsne_list[[input$dataToken]], perplexity==input$perplex & annotation==input$annot),
                     my_perplex = input$perplex, my_annot = input$annot, pointsize = input$pointSize, Alpha = input$pointAlpha)
               )

               output$plot_table_selection <- DT::renderDataTable(
                  brushedPoints(subset(tsne_list[[input$dataToken]], perplexity==input$perplex & annotation==input$annot), input$plot_brush), server = FALSE, selection = 'none'
               )

               # download the filtered data
               output$table_download = downloadHandler(filename = 'selected_points.csv', content = function(file) {
                     write.table(brushedPoints(subset(tsne_list[[input$dataToken]], perplexity==input$perplex & annotation==input$annot), input$plot_brush), file, sep=";",row.names = F)
               })

               # observe main vs secondary data sets
               secondary.datasets = names(tsne_list)[!names(tsne_list) %in% input$dataToken]
               output$dataTokenText2 <- renderText(secondary.datasets[1])
               output$dataTokenText3 <- renderText(secondary.datasets[2])
               secondary.perplexes1 = sort(unique(tsne_list[[secondary.datasets[1]]]$perplexity))
               secondary.perplexes2 = sort(unique(tsne_list[[secondary.datasets[2]]]$perplexity))

               ## data #2 (secondary.datasets 1)
               output$perplexity_selection2 <- renderUI({
                    selectInput(inputId = "perplex2", label = "Available perplexities:",
                           choices=secondary.perplexes1, selected = secondary.perplexes1[3])
                  })
               output$tsne_plot_annot2 <- renderPlot(
                  tsnePlotWrap(my_df = subset(tsne_list[[secondary.datasets[1]]], perplexity==input$perplex2 & annotation==input$annot),
                     my_perplex = input$perplex2, my_annot = input$annot, pointsize = input$pointSize*0.75, Alpha = input$pointAlpha, showGuide=F))


               ## data #3 (secondary.datasets 2)
               output$perplexity_selection3 <- renderUI({
                    selectInput(inputId = "perplex3", label = "Available perplexities:",
                           choices=secondary.perplexes2, selected = secondary.perplexes2[3])
                  })
               output$tsne_plot_annot3 <- renderPlot(
                  tsnePlotWrap(my_df = subset(tsne_list[[secondary.datasets[2]]], perplexity==input$perplex3 & annotation==input$annot),
                     my_perplex = input$perplex3, my_annot = input$annot, pointsize = input$pointSize*0.75, Alpha = input$pointAlpha, showGuide=F))


            })



         # Obseve selected samples
         observe({
            x = brushedPoints(subset(tsne_list[[input$dataToken]],  perplexity==input$perplex & annotation==input$annot), input$plot_brush)$sample_id
            if(length(x)){
               secondary.datasets.sampleselect = names(tsne_list)[!names(tsne_list) %in% input$dataToken]
               output$tsne_plot_samples2 <- renderPlot(
                  tsnePlotWrapSampleHighlight(my_df = subset(tsne_list[[secondary.datasets.sampleselect[1]]], perplexity==input$perplex2 & annotation==input$annot),
                     my_perplex = input$perplex2, my_samples = x, pointsize = input$pointSize*0.75, Alpha = input$pointAlpha, showGuide=F))
               output$tsne_plot_samples3 <- renderPlot(
                  tsnePlotWrapSampleHighlight(my_df = subset(tsne_list[[secondary.datasets.sampleselect[2]]], perplexity==input$perplex3 & annotation==input$annot),
                     my_perplex = input$perplex3, my_samples = x, pointsize = input$pointSize*0.75, Alpha = input$pointAlpha, showGuide=F))
               ## CN DATA - Select sample for CN plot (dropdown list)
               output$SelectGenomeSample <- renderUI({
                  selectInput(inputId = "genomeSample", label = "Genome plot",
                        choices=x, selected = x[1])
               })
            } # end if x
         })


         observe({
            x = input$genomeSample
            if(length(x)){
               output$selectedIndivSample <- renderText(x)
               ## Render CN plot using dlfoo2::segmentPlotter
               if(x %in% names(gdc_cn)){
                  output$genomePlot <- renderPlot(dlfoo2::segmentPlotter(gdc_cn[[x]], sample_name = x))
                  }else{
                  output$genomePlot <- renderPlot(dlfoo2::ggPlotEmpty("Copy number profile not available"))
                  } # end CN plot

               # Plot individual selected sample in all three data sets
               output$tsne_plot_indiv1_name <- renderText(names(tsne_list)[1])

               output$tsne_plot_indiv1 <- renderPlot(
                  tsnePlotWrapSampleHighlight(my_df = subset(tsne_list[[1]], perplexity==input$perplex & annotation==input$annot),
                     my_perplex = input$perplex, my_samples = x, pointsize = input$pointSize*0.75, Alpha = input$pointAlpha, showGuide=F))
               output$tsne_plot_tax1 <- renderPlot(
                  tsnePlotWrap(my_df = subset(tsne_list[[1]], perplexity==input$perplex & annotation=="taxonomy_dl"),
                     my_perplex = input$perplex, my_annot = "taxonomy_dl", pointsize = input$pointSize, Alpha = input$pointAlpha, showGuide = F))
               output$tsne_plot_indiv2_name <- renderText(names(tsne_list)[2])
               output$tsne_plot_tax2 <- renderPlot(
                  tsnePlotWrap(my_df = subset(tsne_list[[2]], perplexity==input$perplex2 & annotation=="taxonomy_dl"),
                     my_perplex = input$perplex2, my_annot = "taxonomy_dl", pointsize = input$pointSize, Alpha = input$pointAlpha, showGuide = F))
               output$tsne_plot_indiv2 <- renderPlot(
                  tsnePlotWrapSampleHighlight(my_df = subset(tsne_list[[2]], perplexity==input$perplex2 & annotation==input$annot),
                     my_perplex = input$perplex2, my_samples = x, pointsize = input$pointSize*0.75, Alpha = input$pointAlpha, showGuide=F))
               output$tsne_plot_indiv3_name <- renderText(names(tsne_list)[3])
               output$tsne_plot_tax3 <- renderPlot(
                  tsnePlotWrap(my_df = subset(tsne_list[[3]], perplexity==input$perplex3 & annotation=="taxonomy_dl"),
                     my_perplex = input$perplex3, my_annot = "taxonomy_dl", pointsize = input$pointSize, Alpha = input$pointAlpha, showGuide = F))
               output$tsne_plot_indiv3 <- renderPlot(
                  tsnePlotWrapSampleHighlight(my_df = subset(tsne_list[[3]], perplexity==input$perplex3 & annotation==input$annot),
                     my_perplex = input$perplex3, my_samples = x, pointsize = input$pointSize*0.75, Alpha = input$pointAlpha, showGuide=F))

               # Plot tsne with indicv samples


               }# end if x
         })

      ## Clinical data (SurvPlot for selected samples)
      # observe({
      #    xx=brushedPoints(subset(tsne_df,  perplexity==input$perplex & annotation==input$annot), input$plot_brush)$patient_barcode
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





#' Specialized function for plotting tsne table with pdata atacched.
#' @description  Should be sorted for perplexity
#' @family plot
#' @family tsne
#' @family shinyPlot
#' @param tsne_tab results list (tSNEres) from run.tsne-analysis
#' @param annotation.column should be a column in the pdata object to for annotation.
#' @param plot.title Title of plot
#' @param color_key what color key to use for annotation. Detaults to dlfoo2::color_subtypes (but could be set to e.g. dlfoo2::color_annotations)
#' @param shape.key what shape key to use. dafaults to dlfoo2::shape_subtypes
#' @param do.gradient set to FALSE if to manually override numeric check for gradient or manual fil
#' @param color.key.gradient what gradient.
#' @param gradient.nomalize TRUE or FALSE if to z normalize annotation value (x-mean(x))/sd(x)
#' @param my.size point size
#' @param my.alpha point alpha
#' @param my.stroke stroke width
#' @param my.stroke.color stroke line colors
#' @param showGuide if to show guide
#' @param coord.fixed if coord_fixed ggplot param
#' @return a ggplot object
#' @export
tsnePlot_special <- function(tsne_tab, annotation.column=NULL, plot.title="", my.size=2, my.alpha=0.7, my.stroke=0.05, my.stroke.color="gray10",
do.gradient=T, palette.gradient=rev(dlfoo2::palette_gradients[["red2white2blue"]]), gadient.na.col=NA, gradient.normalize=F,
showGuide=T, coord.fixed=T, shape.key=dlfoo2::shape_subtypes
){
   require(Rtsne)
   require(dplyr)
   require(tidyverse)
   require(tibble)

   color.key.annotation <- c(dlfoo2::color_subtypes, dlfoo2::color_annotations)
   #shape.key <- dlfoo2::shape_subtypes

   if(!(annotation.column %in% colnames(tsne_tab))) stop("annotation.column not in selected pdata object")

   my_df <- tsne_tab
   my_df$annotation <- tsne_tab[,annotation.column]

   # if fill manual or fill gradient
   if(!is.numeric(my_df$annotation)){
      do.fill = "manual"
      if(!is.factor(my_df$annotation)){
         my_df <- my_df %>% mutate(annotation = factor(annotation)) %>% arrange(annotation)
         color.key <- color.key.annotation[sort(levels(my_df$annotation))]
         }
      if(is.factor(my_df$annotation)){
         my_df <- my_df %>% arrange(annotation)
         color.key <- color.key.annotation[sort(levels(my_df$annotation))]
         }
      }
   if(is.numeric(my_df$annotation) && do.gradient){
      do.fill = "gradient"
      if(gradient.normalize) my_df$annotation <- my_df$annotation - mean(my_df$annotation, na.rm=T) / sd(my_df$annotation, na.rm=T)
      #color.key.gradient = palette_gradientRamp(x = my_df$annotation, my.pal = palette.gradient, palette.sat = gradient.sat, na.col = gadient.na.col)
      color.key.gradient = palette.gradient
      }

   ## if do.shape
   do.shape = F
   if(any(levels(my_df$annotation) %in% names(shape.key))){
      do.shape = T
      # shape.key <- shape_subtypes[sort(levels(my_df$annotation))]
      shape.key <- shape_subtypes[levels(my_df$annotation)]
      }

## GGPLOT
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
            #scale_fill_gradientn(colours=color.key.gradient$col[order(color.key.gradient$val)], na.value=NA)
            scale_fill_gradientn(colours=color.key.gradient, na.value=NA) +
            guides(fill=guide_legend(title=annotation.column))
      }
   if(showGuide==F){
          g <- g + guides(fill=FALSE, shape=FALSE, color=FALSE)
   }
    if(coord.fixed==T) {g <- g + coord_fixed(ratio = 1)}
return(g)
} # end me tSNE






#' function pcaPlot. Plots annotations from pca
#' @family pca
#' @family PlotWrap
#' @param pca.res results object list from pca -analysis
#' @param pdata defaults to pdata_utc_nt
#' @param annotation.column should be a column in the pdata object to for annotation.
#' @param plot.dims what dimensions to plot. defaults to c(1,2)
#' @param color_key what color key to use for annotation. Detaults to dlfoo2::color_subtypes (but could be set to e.g. dlfoo2::color_annotations)
#' @param do.gradient set to FALSE if to manually override numeric check for gradient or manual fil
#' @param color.key.gradient what gradient.
#' @param gradient.nomalize TRUE or FALSE if to z normalize annotation value (x-mean(x)) sd(x)
#' @param gradient.squish.probs if to cap squish gradient to remove effects on color scale of outliers. vector of tow percent values. Uses quantile-probs. Defaults to c(0.05, 0.95),i.e. cap scale at 5% and 95% of values
#' @param my.size point size
#' @param my.alpha point alpha
#' @param my.stroke stroke width
#' @param my.stroke.color stroke line colors
#' @param showGuide if to show guide
#' @return a ggplot object
#' @export
pcaPlot = function(pca.res, pdata=NULL, annotation.column="taxonomy_published", plot.dims = c(1,2), my.size=2, my.alpha=0.7, my.stroke=0.05, my.stroke.color="gray50",
   color.key.annotation=dlfoo2::color_subtypes,
   do.gradient=T, palette.gradient=rev(dlfoo2::palette_gradients[["red2white2blue"]]), gadient.na.col=NA, gradient.normalize=F, gradient.squish.probs=c(0.05, 0.95),
   showGuide=T
   ){

   require(Rtsne)
   require(dplyr)
   require(tidyverse)
   require(tibble)
   require(reshape2)
   if(class(pca.res)!="prcomp") stop("pca.res must be prcomp object")
   if(is.null(pdata)) pdata <- dlfoo2data::pdata_utc_nt
   if(!(annotation.column %in% colnames(pdata))) stop("annotation.column not in selected pdata object")




   my_dims <- paste0("PC",plot.dims)
   if(!any(my_dims %in% colnames(pca.res$x))){stop("Dimensions not among available dims ")}
   dim_mat <- pca.res$x[,my_dims]
   colnames(dim_mat) <- c("dim1","dim2")
   pca_df <- cbind(data.frame(sample_id = rownames(pca.res$x), stringsAsFactors = F), dim_mat)
   my_df <- dplyr::left_join(pca_df, pdata)
   my_df$annotation <- my_df[, annotation.column]
   my_df = as_tibble(my_df)

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


   ## Plot standard taxonomy
   ## ----------------------
      g <- ggplot(my_df) +
            aes(x=dim1, y=dim2) +
            ggtitle(label = annotation.column, subtitle = my_dims) +
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

   return(g)
   } # end me tSNE



#' Specialized function for highlighting single samples in a tSNE plot
#' @description Should be sorted for perplexity
#' @family plot
#' @family tsne
#' @family shinyPlot
#' @param tsne_tab results list (tSNEres) from run.tsne-analysis
#' @param my.samples should be a column in the pdata object to for annotation.
#' @param plot.title plot title
#' @param shape.val vector of length 2 with values for samples and non-samples resp.
#' @param col.val vector of length 2 with values for samples and non-samples resp.
#' @param fill.val vector of length 2 with values for samples and non-samples resp.
#' @param stroke.val vector of length 2 with values for samples and non-samples resp.
#' @param alpha.val vector of length 2 with values for samples and non-samples resp.
#' @param coord.fixed if coord_fixed ggplot param
#' @param shape.key if to plot annotations in different shapes (connected to annotation)
#' @return a ggplot object
#' @export
tsnePlot_sample <- function(tsne_tab, my_samples, plot.title="",
   shape.val = c(4,21), col.val=c("tomato4","gray"), fill.val=c("tomato4","gray"),  size.val=c(2.5, 0.25), stroke.val=c(2.5,0.5), alpha.val=c(1,0.75),
   add_sample_text = T, coord.fixed=T, shape.key=dlfoo2::shape_subtypes
   ){
   require(ggrepel)

   my_df <- tsne_tab
   my_df$plot_val = NA
   my_df$plot_val[my_df$sample_id %in% my_samples] <- "qq"
   my_df <- my_df[(order(my_df$plot_val, na.last = F)),]
   my_df$plot_text <- NA
   my_df$plot_text[my_df$sample_id %in% my_samples] <- as.character(my_df$sample_id[my_df$sample_id %in% my_samples] )


   par_list <- list(
      shape_val = shape.val,
      col_val = col.val,
      fill_val = fill.val,
      size_val = size.val,
      stroke_val = stroke.val,
      alpha_val = alpha.val
   )
   if(!all(unlist(lapply(par_list, function(x) length(x)))==2)) stop("ERROR: all value pars must be of length 2")

   my_pars <- lapply(par_list, function(x, y=my_df$plot_val){
      yy <- factor(y, exclude=NULL)
      xx <- factor(x, exclude=NULL)
      levels(yy) <- xx
      return(as.character(yy))
      })

      my_df2 <-  cbind(my_df, as.data.frame(do.call("cbind", my_pars), stringsAsFactors=F))

     g <-  ggplot(my_df2,
               aes(x=dim1, y=dim2)) +
               labs(x="tSNE dim1", y="tSNE dim2", colour='') +
               ggtitle((plot.title), label = "selected samples") +
               theme(axis.title.y = element_text(size=24)) +
               theme(axis.title.x = element_text(size=24)) +
               theme(legend.text = element_text(size=18)) +
               theme(title = element_text(size=24)) +
               #guides(fill=FALSE, shape=FALSE, color=FALSE) +
               geom_point(
                  colour=as.character(my_df2$col_val),
                  shape=as.numeric(my_df2$shape_val),
                  fill=as.character(my_df2$fill_val),
                  stroke=as.numeric(my_df2$stroke_val),
                  size = as.numeric(my_df2$size_val),
                  alpha=as.numeric(my_df2$alpha_val))

     if(add_sample_text) {g <- g + geom_label_repel(aes(label = plot_text),
                           box.padding   = 0.35,
                           point.padding = 0.5,
                           segment.color = 'grey50')}
      if(coord.fixed==T) {g <- g + coord_fixed(ratio = 1)}

         return(g)
      }


         # tsne_gex_pancan <- readRDS("~/PROJECTS/RCC_ccpRCC_2019/tSNE/tsne_PanCanGDC_Center_sd1_6211genes.rds")
         # pdata <- dlfoo2data::pdata_tcgaPanCan
         # pdata <- pdata %>% filter(platform == "IlluminaHiSeq_RNASeqV2")
         # pdata <- pdata %>% filter(sample_id %in% unique(rownames(tsne_gex_pancan$perplex_10$Y)))
         # pdata <- pdata %>% dplyr::left_join(dlfoo2data::pdata_utc_nt %>% select(sample_id, taxonomy_published, uro_taxonomy))
         # pdata <- pdata %>% arrange(uro_taxonomy)

         # tSNEres <- tsne_gex_pancan
         # tSNEres <-  readRDS("~/PROJECTS/RCC_ccpRCC_2019/tSNE/tsne_UTCnt_Center_sd05.rds")

         # perplex = 25
         # annotation.column="taxonomy_published"

         #       source("~/PROJECTS/RCC_ccpRCC_2019//R_SCRIPTS/R_ccpRCC_2019_SOURCE.R")
         # require(gplots)
         # require(ggplot2)
         #
         # tsne_gex_utc <- readRDS("~/PROJECTS/RCC_ccpRCC_2019/tSNE/tsne_UTCnt_Center_sd05.rds")
         # tsne_me_utc <- readRDS("~/RESOURCES/Methylation_450/Analyses/VMP_CorrGex/utc_me_1126set_VMP_CorrGEX/RCC_1305set_ECR_tSNE/RCC_1305set_ECR_tSNEtSneResults.rds")
         # tsne_mir_utc <- readRDS("~/PROJECTS/RCC_ccpRCC_2019/miRNA/tSNE/tsne_mir_utc_nt_sd05/tsne_mir_utc_nt_sd05tSneResults.rds")
         # tsne_rppa_utc <- readRDS("~/PROJECTS/RCC_ccpRCC_2019/Protein_RPPA/tSNE/tsne_utcnt/tsne_utcnttSneResults.rds")
         # tsne_gex_pancan <- readRDS("~/PROJECTS/RCC_ccpRCC_2019/tSNE/tsne_PanCanGDC_Center_sd1_6211genes.rds")
         #
         # pdata <- dlfoo2data::pdata_utc_nt
         # str(sig_df)
         # pdata <- pdata %>% dplyr::left_join(sig_df)
         # pdata <- pdata %>%
         #    mutate(tax_bxp = tax_simp) %>%
         #    mutate(tax_bxp = factor(tax_bxp, levels=c("ccRCC", "pRCC_e1", "pRCC_e2",
         #       "pCIMP","dCIMP","ccpRCC","chRCC","chONC","pONC","mesRCC"))) %>%
         #    mutate(tax_tsne = tax_dl) %>%
         #    mutate(tax_tsne = if_else(as.character(tax_tsne) %in% c("Urobasal","GU","Basal/SCCL","SC/NE","Mes-Inf","Inf_other"), "BLCA", as.character(tax_tsne))) %>%
         #    mutate(tax_tsne = factor(tax_tsne, levels=rev(c("crtx","crtmed","med","infl1","infl2","ccRCC_e2","ccRCC_e1","ccRCC_e3", "pRCC_e1a","pRCC_e1b", "pRCC_e2",
         #       "pCIMP","dCIMP","ccpRCC","chRCC","chONC","pONC","mesRCC",
         #       "bladder_n","BLCA"))))
         #
         # table(pdata$tax_bxp)
         # # dlfoo2::palette_gradients[["red2white2blue"]]
         # # rev(RColorBrewer::brewer.pal(11, "Spectral"))[-c(1:2)],"black")
         # # my_grad <- c(rev(RColorBrewer::brewer.pal(9,"GnBu"))[-1], RColorBrewer::brewer.pal(9,"YlOrRd"), "black")
         #
         # my_pal <- dlfoo2::color_subtypes
         # my_pal["BLCA"] <- "burlywood1"
         # my_pal["bladder_n"] <- "burlywood3"
   # my_grad <- c(
   #       RColorBrewer::brewer.pal(11,"RdBu")[-c(1:6)],
   #       rev(RColorBrewer::brewer.pal(11,"Spectral")[-c(1:6,10:11)]),
   #       RColorBrewer::brewer.pal(9,"YlOrRd")[-(1:3)])

         # tSNEres=tsne_gex_utc
         # annotation.column = "HNF4_expr"
         # color.key.annotation = my_pal
         # perplex = 25
         # my.size=2.5
         # my.alpha = 0.8
         # gradient.normalize = T
         # my.stroke=0.05
         # my.stroke.color="gray10"
         # do.gradient=T
         # gadient.na.col=NA
         # showGuide=T
         # perplex =  25
         # palette.gradient=my_grad


#' retrieve beta levels for specific probe
#' @description GET EXPRESSION :: Retrieve expression (exprs) for defined probe(s) and eSets - ESETS NOT IN WORKSPACE Returns a list [[1]] : matrix with exopression values of probe(s) [[2]]
#' @description pData info (sample_type and cancer_type) probe.id=c('cg13723431','cg14765933') eset.file =
#' @family old functions
#' @family legacy
#' @family misc
#' @family methylation
#' @param x TCGA aliquot barcode id ('-' or '_' separated)
#' @return numeric vector with beta values
#' @export
dl_getMethylationBeta = function(eset_filename = "/Users/david/RESOURCES/GDAC/DataSets/rcc_me450_867set_gdac.norm.Rdata", my_probes) {
    my.es <- get(load(eset_filename))
    u <- match(as.character(my_probes), featureNames(my.es))
    u <- u[!is.na(u)]
    if(length(u) == 0){stop("no probe found")}
    return(my.es[u,])
   }





#' Wrapper for tSNE analysis
#'
#' \code{tsneWrapper} performs tSNE Aanalysis using Rtsne function in stats package
#' @description OBSOLETE function
#' @description Sorry that I completely forgot to send the scripts. Here is the TSNE script…what were the other ones you were wanting? Plotting functions for TSNE objects?
#' @description Also not that the perplexity values that are set here are between 5 and 50, however this is for running TSNE on patients t(exprs(eset)) and not genes (exprs(eset)). I normally set my perplexity values between 400 and 1000 with increments of 50. seq(400,1000,50).
#' @family OLD
#' @aliases run.rtsne
#' @param x matrix with expression data
#' @return tsne restults - tSNEres
#' @export
#'
tsneWrapper <- function(
    x=NULL, #matrix with expression data
    theta=0.5,
    max_iter=5000,
    outfile=NULL,
    floor_value=NULL,
    z_score=F,
    log=F,
    variance=NULL,
    perplex=c(5:50)
   ){

   require(Rtsne)
   if(log){
      x<-log(x,2)
      x[which(x<0)] <-0
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

    number<-dim(x)[1]
    #How many samples are in the matrix
    #min_perplex<-400 # As suggested by author, minimal perplexity=5
    #max_perplex<-401#round((number-1)/3,digits=0) #Set the max perplexity to (number-1)/3
    #if(max_perplex>50){max_perplex=50} #set a maximum perplexity to try

    ## res<-list(perplex)
    res <- vector(mode = "list", length=length(perplex))
    names(res) <- paste0("perplex_",perplex)

    cat("#T-SNE result file from Rtsne\n", file=outfile, append=FALSE)
    for(i in 1:length(perplex)) {
      set.seed(42) #set the seed for every iteration
      res[[i]] <- Rtsne::Rtsne(as.matrix(x), theta=theta, max_iter=max_iter, perplexity=perplex[i], check_duplicates = F) #Perform Rtsne function
      rownames(res[[i]]$Y)<-rownames(x) #Add the samplename names

      #Header and result, cat sets the first col
      cat("\n#H:perplexity_",perplex[i],"\t", file=outfile, append=TRUE)
      write.table(res[[i]]$Y, file=outfile, col.names = F, quote=F, sep="\t", append=T)
    }
    #rtsne.data<-list(min_perplex=min(perplex), max_perplex=max(perplex), res=res)
    #return(rtsne.data)
    return(res)
   }


#' Older function that uses vanilla lines to generate genome
#' @family PlotWrap
#' @family genomics
#' @family old functions
#' @aliases genomeHeatmap
#' @param x data frame with segments
#' @param sample_vec what samples to use
#' @param add.wg.space = 10e6,
#' @param plot.chrom = c(1:23),
#' @param my.pal = rev(brewer.pal(11,"PuOr")[c(3:9)])
#' @return plot to quartz
#' @export
#'
plotGenomeHeatmap <- function(x, sample_vec, add.wg.space = 10e6, plot.chrom = c(1:23), my.pal = rev(brewer.pal(11,"PuOr")[c(3:9)])){
      chrom <- toChr(plot.chrom)	## Uses toChr function (external) to translate numeric chr's to character
   	u <- which(x$chrom %in% chrom)
   	if(length(u)<1){stop("no segments matched specified plot.chrom")}
   	x <- x[u,]
   	rm(u)
   	sampleIDs <- unique(x$ID)

   	# For each sample, segments are plotted in order of increasing (absolute) log.ratio :: i.e. change order within the x data frame
   	# this to not 'overplot' small high-level alterations (cosmetics...)
   	# So order x according to i) sample_vec amd ii) seg.mean
		u <- order(abs(x$seg.mean), decreasing=F)
		x <- x[u,]
	  	rm(u)
	   u <- match(x$ID, sample_vec)
		uu <- order(u, decreasing=F)
		x <- x[uu,]
   	rm(u, uu)

   	# Define genome coordinates for all segments
   	# Coordinates based on full genome positions
   	wg.start <- genomicPosition(pos=x$loc.start, chrom=x$chrom, force.full.genome.positions=T, add.space=add.wg.space)
   	wg.end <- wg.start + (x$loc.end-x$loc.start)
   	sample_pos_vec <- match(x$ID, unique(x$ID))
   	color_vec <- heatcol_creator2(x$seg.mean, palette.sat = c(-0.5, 0.5), my.pal = my.pal)$col

      plot.new()
      plot.window(xlim=c(0, max(sample_pos_vec)), ylim=c(-1*genomicPosition(pos = 156040895, chrom = 23, add.space=add.wg.space, force.full.genome.positions = T) ,0))
      rect.args <- list(
                     xleft = sample_pos_vec-1,
                     ybottom = -1*wg.start,
                     xright = sample_pos_vec,
                     ytop = -1*wg.end,
                     col = as.character(color_vec),
                     border = NA, lwd=0.1)
         do.call("rect", args = rect.args)
    } # edn foo genome heatmap




#' OLD function tsnePlot
#' Plots RCC and BLCA taxonomy annotations
#' @family tsne
#' @family old
#' @family PlotWrap
#' @param tSNEres results list (tSNEres) from run.tsne-analysis
#' @param perplexity what perplexity to use
#' @param my.size point size
#' @param my.alpha point alpha
#' @return a ggplot object
#' @export
tsnePlotTaxOLD <- function(tSNEres, perplex = NULL, my.size=2, my.alpha=0.7){
   require(Rtsne)
   require(dplyr)
   require(tidyverse)

   data("pdata_utc_nt", package="dlfoo2")
   data("color_subtypes", package="dlfoo2")
   data("shape_subtypes", package="dlfoo2")

   pdata <- pdata_utc_nt[,c("sample_id", "taxonomy_simplified", "taxonomy_dl", "sample_type2")]

   tsne_df <- do.call("rbind", lapply(tSNEres, tsne2df))
   available.perplexes <- sort(unique(tsne_df$perplexity))
   if(!is.null(perplex)) if(!(perplex %in% available.perplexes)){stop("perplex not among available perplexes ",paste(available.perplexes,collapse = " "))}
   tsne_df <- dplyr::left_join(tsne_df, pdata)
   #tsne_df$taxonomy_dl <- addNA(tsne_df$taxonomy_dl)

   ## Plot standard taxonomy
   ## ----------------------
   message("... ... plotting")
   shape_key <- shape_subtypes[(sort(unique(tsne_df$taxonomy_dl)))]

   if(is.null(perplex)){
      message("... ... ... plotting all perplexes")
      g_tsne <- ggplot(tsne_df) +
            aes(x = dim1, y = dim2) +
            facet_wrap(~perplexity, scales = "free") +
            theme(panel.background = element_rect(fill = 'gray95', colour = 'antiquewhite3', size = 0.5)) +
            theme(panel.grid.major =  element_line(colour = 'lightsteelblue2'), panel.grid.minor =  element_line(colour = 'mistyrose3')) +
            labs(x="dim1", y="dim2" ) +
            theme(axis.title.y = element_text(size=12)) +
            theme(axis.title.x = element_text(size=12)) +
            theme(axis.text = element_text(size=9)) +
            theme(title = element_text(size=9)) +
            theme(strip.text = element_text(size=12)) +
            theme(title = element_text(size=12)) +
            ggtitle(label = "tSNE Taxonomy") +
            geom_point(aes(fill=taxonomy_dl, shape=taxonomy_dl), colour="gray20",stroke=0.25, size=my.size, alpha=my.alpha) +
            scale_fill_manual(values = color_subtypes, na.value = "black") +
            scale_shape_manual(values = as.numeric(shape_key), na.value = 4)
         } # end plot all perplxes

   if(is.numeric(perplex)){
      message("... ... ... plotting perplex ", perplex)
      g_tsne <- ggplot(subset(tsne_df, perplexity==perplex)) +
            aes(x = dim1, y = dim2) +
            theme(panel.background = element_rect(fill = 'gray95', colour = 'antiquewhite3', size = 0.5)) +
            theme(panel.grid.major =  element_line(colour = 'lightsteelblue2'), panel.grid.minor =  element_line(colour = 'mistyrose3')) +
            labs(x="dim1", y="dim2" ) +
            theme(axis.title.y = element_text(size=12)) +
            theme(axis.title.x = element_text(size=12)) +
            theme(axis.text = element_text(size=9)) +
            theme(title = element_text(size=9)) +
            theme(strip.text = element_text(size=12)) +
            theme(title = element_text(size=12)) +
            ggtitle(label = "tSNE Taxonomy") +
            geom_point(aes(fill=taxonomy_dl, shape=taxonomy_dl), colour="gray20",stroke=0.25, size=my.size, alpha=my.alpha) +
            scale_fill_manual(values = color_subtypes, na.value = "black") +
            scale_shape_manual(values = as.numeric(shape_key), na.value = 4)
         } # end plot all perplxe
   return(g_tsne)
   }

#' Generic plot for enrichment analyses v1.0
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
enrichmentPlot_bubble <- function(
      ea.tab=NULL, size, pval,
      order.by.clustering=T,
      pval.pal = c(rev(RColorBrewer::brewer.pal(9,"YlOrRd")[3:9]), RColorBrewer::brewer.pal(9,"PuBuGn")[3:9]),
      reverse.x=F, reverse.y=F
   ){
   require(tibble)
   require(tidyr)
   require(ggplot2)


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






#'
#' #' Normalize Illumina Methylation Ringner/Staaf
#' #' @description function meNormalizeBetaV2 (filemame instead of object)
#' #' @description Version 2.0 20190108
#' #' @family methylation
#' #' @family genoset
#' #' @family genomics
#' #' @family transformation
#' #' @param fileName object of class GenoSet
#' #' @param runFolder directory where processed genoset object should be saved to
#' #' @param projectLabel character vector; name of project
#' #' @param probe.annotations.file file name for probe annotations that inlcude 'INFINIUM_DESIGN_TYPE' column: type I / II ILMN probe sets
#' #' @return genoset
#' #' @export
#' meNormalizeBetaV2 <- function(fileName, runFolder, projectLabel=NULL, probe.annotations.file='~/RESOURCES/Annotation Files/probe.annotations.Methylation.Ringner.20031112.Rdata'){
#'    require(genoset)
#'    if(length(grep(".rds", fileName))!=1) stop("fileName must point to .rds object")
#'    x <- readRDS(fileName)
#'    if(class(x)!="GenoSet") stop("gset is not of class 'GenoSet'")
#'    runName <- paste0(gsub(".rds","",basename(fileName)),"_norm")
#'    message("... creating runName from input filename. Adding '_norm: '", runName)
#'
#'
#'    find.max<-function(x,span=1,method="strict"){
#'       # x = density vector or series vector with values
#'       # span = number of items on each side of the tested point. 1 = window of 3, i-1,i,i+1
#'       #function supposes that the data is correctly sorted,
#'       # and that the data is relatively smooth
#'       maxima<-c()
#'       for(t in 1:length(x)){
#'          span.low<-t-span
#'          span.high<-t+span
#'          if(method=="strict"){
#'             if( (span.low<0) | (span.high>length(x))){
#'                #not valid for a check as the window spans outside the vector indexes
#'             }else{
#'                if(which.max(x[span.low:span.high])==(span+1)) maxima<-c(maxima,t) #this is a local maxima
#'             }}}
#'       maxima.order<-order(x[maxima],decreasing=TRUE)
#'       return(maxima[maxima.order])
#'    } # end funtion find.max
#'
#'    ## select peaks function
#'    select.peaks<-function(m){
#'       peaks<-c(2,-1)
#'       for(t in 1:nrow(m)){
#'          if((m$x[t]<0.4) & (m$x[t]<peaks[1]) & (m$y[t]>(0.05*m$y[1]))) {peaks[1]=m$x[t]}
#'          if((m$x[t]>0.6) & (m$x[t]>peaks[2]) & (m$y[t]>(0.05*m$y[1]))) {peaks[2]=m$x[t]}
#'       } # inf for t
#'       if(peaks[1]==2) {peaks[1]=0}
#'       if(peaks[2]==-1) {peaks[2]=1}
#'       return(peaks)
#'    } # end function select peaks
#'
#'
#'    # myDir <- file.path(runFolder, runName)
#'    myDir <- runFolder
#'    if(!dir.exists(myDir)) dir.create(myDir)
#'
#'
#'
#'    # Log
#'    logFile <- file.path(myDir, "logfile_Normalization.txt")
#'    cat(file = logFile, "PARAMS\n")
#'    cat(file = logFile, append = T, "runName:", runName, "\n")
#'    cat(file = logFile, append = T, "runFolder:", runFolder,"\n")
#'    cat(file = logFile, append = T, "Processd data directory:", myDir,"\n")
#'    cat(file = logFile, append = T, "Probe info file:", probe.annotations.file,"\n")
#'
#'    cat(file = logFile, append = T, "\nDATA\n")
#'    cat(file = logFile, append = T, "dim input genoSet: ", dim(x),"\n")
#'    if(length(metadata(x))>0){
#'       lapply(names(metadata(x)), function(y){
#'          cat(file = logFile, append = T, y,": ",metadata(x)[[y]],"\n")
#'       })
#'    }
#'    cat(file = logFile, append=T, "\n------\n")
#'
#'
#'    message("... reading probe annotations. matching genoset to probes in file")
#'    probe.annotations <- get(load(probe.annotations.file))
#'    my_probes <- intersect(rownames(x), probe.annotations$ILMNID)
#'    cat(file = logFile, append = T, "n features intersect with probe info file: ", length(my_probes),"\n")
#'    x <- x[my_probes]
#'    probe.annotations <- probe.annotations[my_probes,]
#'    # str(probe.annotations)
#'    message("... Defining Type I/II probes")
#'    probes_I = probe.annotations[(which(probe.annotations$INFINIUM_DESIGN_TYPE=="I")),]$ILMNID
#'    probes_II = probe.annotations[(which(probe.annotations$INFINIUM_DESIGN_TYPE=="II")),]$ILMNID
#'
#'    data = x[,,"beta"]
#'    dim(data)
#'    data_I=data[which(!is.na(match(rownames(data),probes_I))),]
#'    data_II=data[which(!is.na(match(rownames(data),probes_II))),]
#'
#'    peaks <- data.frame(unmet_I=rep(NA,ncol(data)),met_I=rep(NA,ncol(data)),unmet_II=rep(NA,ncol(data)),met_II=rep(NA,ncol(data)))
#'    rownames(peaks)=colnames(data)
#'
#'    plotFile = file.path(myDir, 'peaks_preNorm.pdf')
#'    pdf(plotFile)
#'    message("... plotting peaks pre normalization")
#'    for(sample in 1:ncol(data)) {
#'       cat("loop", sample, "\n")
#'       d_I=density(data_I[,sample],kernel="epanechnikov",bw=0.02,from=0,to=1,na.rm=T)
#'       d_II=density(data_II[,sample],kernel="epanechnikov",bw=0.02,from=0,to=1,na.rm=T)
#'       #d_Ir=density(data_Ir[,sample],kernel="epanechnikov",bw=0.02,from=0,to=1,na.rm=T)
#'       #d_Ig=density(data_Ig[,sample],kernel="epanechnikov",bw=0.02,from=0,to=1,na.rm=T)
#'
#'       plot(d_I,col="blue",main="")
#'       lines(d_II,col="orange")
#'       # lines(d_Ir,col="red")
#'       #lines(d_Ig,col="green")
#'       title(colnames(data)[sample])
#'       max_I<-find.max(d_I$y,span=3);
#'       max_I=cbind(max_I,d_I$x[max_I])
#'       max_I=cbind(max_I,d_I$y[max_I])
#'       colnames(max_I)=c("index","x","y")
#'       max_II<-find.max(d_II$y,span=3);
#'       max_II=cbind(max_II,d_II$x[max_II])
#'       max_II=cbind(max_II,d_II$y[max_II])
#'       colnames(max_II)=c("index","x","y")
#'       peaks_I=select.peaks(as.data.frame(max_I))
#'       peaks_II=select.peaks(as.data.frame(max_II))
#'       abline(v=peaks_I,col="blue")
#'       abline(v=peaks_II,col="orange")
#'       peaks[colnames(data)[sample],]=append(peaks_I,peaks_II)
#'    } # end sample loop
#'    dev.off()
#'
#'    ## Normalize!
#'    #####################
#'    message("... normalizing")
#'    data.norm=as.matrix(data)
#'    rows_I=which(!is.na(match(rownames(data.norm),probes_I)))
#'    rows_II=which(!is.na(match(rownames(data.norm),probes_II)))
#'    tmp=as.matrix(data)
#'    for(sample in 1:ncol(data.norm)) {
#'       data.norm[rows_I,sample]=(tmp[rows_I,sample]-peaks[sample,1])/(peaks[sample,2]-peaks[sample,1])
#'       data.norm[rows_II,sample]=(tmp[rows_II,sample]-peaks[sample,3])/(peaks[sample,4]-peaks[sample,3])
#'       }
#'    rm(tmp)
#'    gc()
#'
#'    data.norm=as.data.frame(data.norm)
#'    data.norm[]=lapply(data.norm, function(x) ifelse (x<0,0,x))
#'    data.norm[]=lapply(data.norm, function(x) ifelse (x>1,1,x))
#'
#'    data.norm_I=data.norm[which(!is.na(match(rownames(data.norm),probes_I))),]
#'    data.norm_II=data.norm[which(!is.na(match(rownames(data.norm),probes_II))),]
#'
#'    plotFile = file.path(myDir, 'peaks_postNorm.pdf')
#'    #file.name = paste(eset.dir, eset.name,'_peaks2.pdf',sep='')
#'    pdf(plotFile)
#'    op<-par(mfrow=c(2,1))
#'    for(sample in 1:ncol(data)) {
#'
#'       d_I=density(data_I[,sample],kernel="epanechnikov",bw=0.02,from=0,to=1,na.rm=T)
#'       d_II=density(data_II[,sample],kernel="epanechnikov",bw=0.02,from=0,to=1,na.rm=T)
#'       d.norm_I=density(data.norm_I[,sample],kernel="epanechnikov",bw=0.02,from=0,to=1,na.rm=T)
#'       d.norm_II=density(data.norm_II[,sample],kernel="epanechnikov",bw=0.02,from=0,to=1,na.rm=T)
#'
#'       plot(d_I,col="blue",main="")
#'       lines(d_II,col="orange")
#'       abline(v=peaks[sample,1:2],col="blue")
#'       abline(v=peaks[sample,3:4],col="orange")
#'       title(colnames(data.norm)[sample])
#'
#'       plot(d.norm_I,col="blue",main="")
#'       lines(d.norm_II,col="orange")
#'
#'       }
#'    par(op)
#'    dev.off()
#'
#'    stopifnot(all(dim(data) == dim(data.norm)))
#'    message("... creating normalized genoset")
#'    gs <- x
#'    message("... ... rounding beta values to 4 digits")
#'    gs[,,"beta"] <- round(data.norm,4)
#'
#'
#'    message("... adding project metadata")
#'    if(!is.null(projectLabel)) metadata(gs)[["projectLabel"]] <- projectLabel
#'    metadata(gs)[["meNormalizeBeta"]] <- timestamp()
#'    cat(file = logFile, append=T, "OUTPUT\n")
#'    cat(file = logFile, append=T, "\n------\n")
#'    cat(file = logFile, append = T, "dim normalized genoSet:", dim(gs) ,"\n")
#'    lapply(names(metadata(gs)), function(y){
#'          cat(file = logFile, append = T, y,": ",metadata(x)[[y]],"\n")
#'       })
#'    print(gs)
#'
#'    message("... saving normalized genoSet to file")
#'    outFile <- paste0(myDir, "/", runName, ".rds")
#'    saveRDS(gs, outFile)
#'    cat(file = logFile, append = T, "normalized genoSet saved to:", myDir ,"\n")
#'    message("... DONE!")
#'
#' } # end function meNormalizeBeta
#'
