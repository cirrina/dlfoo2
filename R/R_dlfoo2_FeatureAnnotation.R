#
# sig_list <- list(
#    clark_2019_kidney = list(species="mm", file="~/RESOURCES/scRNAseq/Clark_2019_scRNA_mmKidney/Supplemental/TableS1_Clark2019_CellType_DataBase.txt",
#          cols=c("Marker.Gene.Symbol","Cell.Type"), type="gene_symbol"),
#    hochane_2019_fetal = list(species="hs", file="~/RESOURCES/scRNAseq/Hochane_2019_FetalKidney/Supplemental/Table_s3_top4_cell_cluster_markers_AllCandidates.txt",
#          cols=c("ensembl_gene_id","cluster.of.interest"), type="ensembl_gene_id"),
#    menon_2018_fetal = list(species="hs", file="~/RESOURCES/scRNAseq/Menon_2018_scRNA_Fetal_Kidney/Annotations/Manon_2018_Tab1_Cluster_Anchors.txt",
#          cols=c("SYMBOL","population"), type="gene_symbol"),
#    park_2018_kidney = list(species="mm", file="~/RESOURCES/scRNAseq/Park_Susztak_2018_scSEQ_mm_Kidney/Supplemental/Population_markers_s1s2.txt",
#          cols=c("SYMBOL.mm","population"), type="gene_symbol"),
#    wang_2018_fetal = list(species="hs", file="~/RESOURCES/scRNAseq/Wang_2018_scRNA_FetalKidney/Supplemental/TableS2_ClusterGenes_13clusters.txt",
#          cols=c("SYMBOL","cluster"), type="gene_symbol"),
#    wu_2018_allograft = list(species="hs", file="~/RESOURCES/scRNAseq/Wu_2018_scSeq_KidneyAllograft/Supplemental/Wu_2018_nSeq_KidneyAllograft_STable_CellTypeMarkers.txt",
#          cols=c("gene","cell.type"), type="gene_symbol"),
#    young_2018_litterature = list(species="hs", file="~/RESOURCES/scRNAseq/Young_2018_Science_scRNA_kidney/Annotations/Young_anchor_genes_litterature.txt",
#          cols=c("Gene","Marker_of"), type="gene_symbol"),
#    young_2018_anchors = list(species="hs", file="~/RESOURCES/scRNAseq/Young_2018_Science_scRNA_kidney/Annotations/Young_population_anchors_Figures.txt",
#       cols=c("SYMBOL","population"), type="gene_symbol")
#    )
#
#
# sigs <- lapply(names(sig_list), function(x) {
#    y <-read.delim(file = sig_list[[x]]$file, as.is=T)
#    y <- as_tibble(y) %>% dplyr::select(sig_list[[x]]$cols, everything())
#    colnames(y)[1:2] <- c(sig_list[[x]]$type, "cell.type")
#
#    symbol_headers <- c("hgnc_symbol","mgi_symbol","rgd_symbol")
#    names(symbol_headers) <- c("hs","mm","rn")
#    if(sig_list[[x]]$type=="gene_symbol"){
#       colnames(y)[1] <- symbol_headers[sig_list[[x]]$species]
#       }
#    yy <- dlfoo2::featureGetBM(query.id = y[,1], query.type = sig_list[[x]]$type, query.species = sig_list[[x]]$species)
   # if(sig_list[[x]]$species0=="mm") return(dplyr::left_join(y, yy$hs.chain, by="mgi_symbol"))
   # if(sig_list[[x]]$species=="hs") return(dplyr::left_join(y, as_tibble(yy[["hs"]]), by=sig_list[[x]]$type))
   # })


#' Uses Biomart to retrieve gene info.
#' @description Use one gene/feature id to create object with multiple feature annotations, e.g. symbol, entrez, ENSG, and for multiple species
#' @description  The featureIdObjectSignature is can be used in downstream analyses
#' @description  For homologs - enseble is always used as primary id type when collapsing duplicate rows
#' @description NOTE!: update 20190916 . entrezgene changed to entrezgene_id
#' @family feature annotation
#' @family biomart
#' @param query.id Input one or multiple gene identifiers as character vector (default is hgnc_gene symbol)
#' @param query.type What type of identifer: c("hgnc_symbol","entrezgene_id","ensembl_gene_id")
#' @param query.species What species (hs, mm, or rn)
#' @param hs.genome if to retrieve ENSG gene info from hg19 database instead of hg38. defauts to hg38
#' @param mm.genome only mm10 supported so far
#' @param rn.genome only rn6
#' @param collapse.to.unique if some input id's will generate nultiple output entries (e.g. one symbol may be associated with multiple ENSGs, e.e.  1:n). Should these be collapsed by piping or not?
#' @param keep.chromosomes Gebns on what chrosomosomes to keep. Defaults to c(1:22, "X","Y","MT"). set to NULL if to keep all genomic chromsomes, contigs and scaffolds
#' @param transcript.biotypes what transcript biotypes to keep. Defaults to 'protein_coding' only. To hkepp all, set to NULL.
#' @return object list with slots
#' @export
featureGetBM <- function(
   query.id = 'all_genes',
   query.type=c("ensembl_gene_id","gene_symbol","entrezgene_id"),
   query.species = c("hs","mm","rn"),
   hs.genome = c('hg38','hg19'),
   mm.genome = c('mm10'),
   rn.genome = c('rn6'),
   keep.chromosomes = c(1:22, "X","Y","MT"),
   transcript.biotypes = c("protein_coding"),
   collapse.to.unique=TRUE,
   get_orthologs=FALSE
      ){
   require(biomaRt)
   require(GenomicRanges)
   require(genoset)
   require(dplyr)
   require(tibble)
   require(tidyr)

   # check input args
   query.id <- as.character(query.id)
   query.type = match.arg(query.type, choices = c("gene_symbol","entrezgene_id","ensembl_gene_id"), several.ok = F)
   query.species = match.arg(query.species, choices =  c("hs","mm","rn"), several.ok = F)
   hs.genome = match.arg(hs.genome, choices = c('hg38','hg19'), several.ok = F)
   mm.genome = match.arg(mm.genome, choices= c("mm10"), several.ok = F)
   rn.genome = match.arg(rn.genome, choices= c("rn6"), several.ok = F)
   biotype_choices <- c("3prime_overlapping_ncRNA","antisense","bidirectional_promoter_lncRNA","IG_C_gene","IG_C_pseudogene","IG_D_gene","IG_J_gene","IG_J_pseudogene","IG_pseudogene","IG_V_gene","IG_V_pseudogene","lincrna","lincRNA","macro_lncRNA","miRNA","misc_RNA","Mt_rRNA","Mt_tRNA","non_coding","non_stop_decay","nonsense_mediated_decay","polymorphic_pseudogene","processed_pseudogene","processed_transcript","protein_coding","pseudogene","retained_intron","ribozyme","rRNA","rRNA_pseudogene","scaRNA","scRNA","sense_intronic","sense_overlapping","snoRNA","snRNA","sRNA","TEC","TR_C_gene","TR_D_gene","TR_J_gene","TR_J_pseudogene","TR_V_gene","TR_V_pseudogene","transcribed_processed_pseudogene","transcribed_unitary_pseudogene","transcribed_unprocessed_pseudogene","translated_processed_pseudogene","unitary_pseudogene","unprocessed_pseudogene","vaultRNA")
   # if(is.null(transcript.biotypes)) transcript.biotypes<-biotype_choices
   if(!is.null(transcript.biotypes))  transcript.biotypes <- match.arg(transcript.biotypes, choices=biotype_choices, several.ok = T)

   symbol_headers <- c("hgnc_symbol","mgi_symbol","rgd_symbol")
   names(symbol_headers) <- c("hs","mm","rn")
   if(query.type=="gene_symbol"){
      query.type <- symbol_headers[query.species]
      }
   if(hs.genome == "hg19" &&  query.type=="hgnc_symbol") message(" ... Beware that some genes have changed Ensembl ID between hg19 and hg38. \n ... Thus, some 100's genes (ENSGs) in TCGA/GDC data chromosome locations will not available in biomart. \n ... ... e.g. HNF1B and uses a different ENSG in hg19. These must be otained using entrezgene id or gene symbol instead")

   # add args
   species_vec <- c("hs","mm","rn")

   ## Create output list
   featureId <- list(
      query.id=query.id,
      query.type=query.type,
      query.species=query.species,
      hs=NULL, hs.chain=NULL,
      mm=NULL, mm.chain=NULL,
      rn=NULL, rn.chain=NULL,
      hs.genome=hs.genome, mm.genome=mm.genome, rn.genome=rn.genome)

   message("\n<<< START >>>\n --> Fetching gene identifiers for ",length(query.id)," ids using bomaRt\n")
   MartAttributesList <- list(
      hs = (c("hgnc_symbol","entrezgene_id","ensembl_gene_id","start_position", "end_position", "chromosome_name","strand", "transcript_biotype")),
      mm = (c("mgi_symbol","entrezgene_id","ensembl_gene_id","start_position", "end_position", "chromosome_name","strand", "transcript_biotype")),
      rn = (c("rgd_symbol","entrezgene_id","ensembl_gene_id","start_position", "end_position", "chromosome_name","strand", "transcript_biotype")))

   ## Defeine & query species-specific mart
   MartArgsList <- list(
      hg19 = list(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org",path="/biomart/martservice",dataset="hsapiens_gene_ensembl"),
      hg38 = list(biomart="ensembl", path="/biomart/martservice", dataset="hsapiens_gene_ensembl"),
      mm10 = list(biomart="ensembl", path="/biomart/martservice",dataset="mmusculus_gene_ensembl"),
      rn6 = list(biomart="ensembl", path="/biomart/martservice",dataset="rnorvegicus_gene_ensembl"))
   myMartArgs <- MartArgsList[[ featureId[[ paste0(featureId$query.species,".genome")]] ]]
   message("--> calling useMart")
   myMart <-  do.call("useMart", myMartArgs)

   message("--> calling getBM function")
   if(featureId$query.id[1]=="all_genes"){
      myBM <- getBM(attributes=MartAttributesList[[featureId$query.species]], mart=myMart)
      }else{
      myBM <- getBM(attributes=MartAttributesList[[featureId$query.species]], values=featureId$query.id, filter=featureId$query.type, mart=myMart)
      }
   if(nrow(myBM)<1){stop("Cannot find ",query.id," in Ensembl database")}
   message("   --> OK!")

   ## Polish table
   ## ............
   message("... ... retrieved data frame with ", nrow(myBM), " rows")
   myBM$entrezgene_id <- as.character(myBM$entrezgene_id)
   myBM$strand <- ifelse(myBM$strand=="-1", "-", myBM$strand)
   myBM$strand <- ifelse(myBM$strand=="1", "+", myBM$strand)
   myBM$strand <- ifelse(myBM$strand %in% c("+","-"), myBM$strand, "*")

   ## Remove any NA/blank in primary identifier
   ## ................
   myBM_orig <- myBM
   if(any(is.na(myBM[,featureId$query.type])) | any(myBM[,featureId$query.type]=="") ){
      message(" ... Removing NA/missing primary_ids: ", featureId$query.type)
      myBM <- myBM[!c(is.na(myBM[,featureId$query.type]) | myBM[,featureId$query.type]==""),]
      message("... ... data frame of ", nrow(myBM), " rows after filtering blank/NA")
      }

   ## Polish chromsome
   ## ................
   if(!is.null(keep.chromosomes)){
      message("--> Filter")
      message(" ... ... keeping only selected chromosomes: ", paste(keep.chromosomes, collapse = ", "))
      myBM <- myBM[myBM$chromosome_name %in% keep.chromosomes,]
      message("... ... data frame of ", nrow(myBM), " rows")
      }

   ## Polish transcript biotype
   ## ................
    if(!is.null(transcript.biotypes)){
      message("--> Filter")
      message(" ... ... keeping only selected transcript biotypes: ", paste(transcript.biotypes, collapse = ", "))
      myBM <- myBM[myBM$transcript_biotype %in% transcript.biotypes,]
      message("... ... data frame of ", nrow(myBM), " rows")
      }


   ## Is query table 1:n of primary query.type
   ## ..............
   if(length(unique(myBM[,featureId$query.type])) < nrow(myBM) ){
      if(collapse.to.unique){
         message("\n--> Collapsing rows/features in biomart tables")
         mult_list <- split(myBM[,], myBM[,featureId$query.type])
         u <- unlist(lapply(mult_list, nrow))
         message(" ... ... mutiple matches (1:n) found for: ", length(which(u>1)), " id's")
         message(" ... ... coercing these using pipes '|'")

         mult_list[u>1] <- lapply(mult_list[u>1], function(x){
            # 1: for genomic location - if chromosome doesnt match - set all to NA
               if(!all(x$chromosome_name==x$chromosome[1])){
                  x$chromosome_name=NA
                  x$start_position=NA
                  x$end_position=NA}
               if(!all(x$strand==x$strand[1])){
                  x$strand="*"}
               # 2: for genomic location: use the widest region
               x$start_position <- min(x$start_position)
               x$end_position <- max(x$end_position)
               # 3: Pipe remaining identifiers (remove NA/blank)
               apply(x, 2, function(y){
                  yy <- unique(y)
                  yy <- yy[!is.na(yy)]
                  yy <- yy[!yy==""]
                  paste(yy, collapse = "|")
                  })
               })

         myBM_merged <- do.call("rbind", mult_list)
         if(any(duplicated(myBM_merged[,featureId$query.type]))) stop("FAILED GENERATING UNIQUE PRIMARY IDs")
         message("   --> OK!")
         message("... ... data frame of ", nrow(myBM_merged), " rows after merging rows")

         if(any(myBM_merged$start_position=="")){
            message(" ... Removing primary_ids with no/multiple mapped location(s) after merging")
            myBM_merged <- myBM_merged[!(myBM_merged$start_position==""),]
            message("... ... data frame of ", nrow(myBM_merged))
         }
         myBM <- myBM_merged
         message("\n   --> DONE!")
         message("... ... done polishing primary biomart object\n ... ... ", nrow(myBM), " rows after merging rows")
      }}


      ## Create  genomic ranges object
      ## .......................
      message("\n   --> Generating GRanges object of biomart table")
      probesGR <- GRanges(seqnames = paste0("chr", myBM$chromosome_name), strand=as.factor(myBM$strand),
             ranges=IRanges(start=as.numeric(myBM$start_position), end=as.numeric(myBM$end_position)))
      elementMetadata(probesGR) <- myBM[,-match(c("start_position", "end_position", "chromosome_name","strand"), colnames(myBM))]
      genome(probesGR) <- featureId[[ paste0(featureId$query.species, ".genome")]]
      probesGR <- toGenomeOrder(probesGR)
      if(collapse.to.unique) names(probesGR) <- elementMetadata(probesGR)[,featureId$query.type]

      metadata(probesGR)[["query.time"]] <- timestamp()
      metadata(probesGR)[["query.id"]] <- query.id
      metadata(probesGR)[["query.type"]] <- featureId$query.type
      metadata(probesGR)[["biomart.table"]] <- myBM_orig

      featureId[[featureId$query.species]] = probesGR





      ## ENSEMBL ORTHOLOGS - Always performed through ENSEMBL id
      ## ...........
      if(!get_orthologs) cat("... get_orthologs set to FALSE (default). No species orthologs are obtained")
      if(get_orthologs){
      species_homologs <- species_vec[!species_vec %in% featureId$query.species]
      message("\n --> Get orthologs")
      message("... Retreaving ensembl homologs from ", paste(species_homologs, collapse = " and "))
      MartArgsHomologs <- list(
         hs = c("hsapiens_homolog_ensembl_gene"),
         mm = c("mmusculus_homolog_ensembl_gene"),
         rn = c("rnorvegicus_homolog_ensembl_gene"))

      #myMartHomolog <- lapply(species_homologs, function(x){
      for(x in species_homologs){
         message(" --> getBM for ", MartArgsHomologs[[x]])

         hom_attributes <- c("ensembl_gene_id", MartArgsHomologs[[x]])

         homologBM = getBM(
            attributes=c(hom_attributes),
            values=elementMetadata(featureId[[featureId$query.species]])[,featureId$query.type],
            filter=featureId$query.type,
            mart=myMart)

         if(nrow(homologBM)<1){
            message("... ... no ensembl homologs found for ", x)
            #getBM(attributes="mmusculus_homolog_ensembl_gene", values=featureId$hs$ensembl_gene_id, filter="ensembl_gene_id", mart=mart_hs)
            y <- as.data.frame(matrix(nrow=0, ncol=length(MartAttributesList[[x]])))
            colnames(y) <- MartAttributesList[[x]]
            featureId[[x]] <- y
            next(x)}

         ## Use ensembl homologs to build BM table
         message("   --> OK")
         yMart <- do.call("useMart", MartArgsList[[ featureId[[ paste0(x, ".genome")]] ]])
         message(" --> getBM of ",nrow(homologBM)," retrieved ensembl homologs ")
         y = getBM(attributes=MartAttributesList[[x]], values=homologBM[,MartArgsHomologs[[x]]], filter="ensembl_gene_id", mart=yMart)
         message("   --> OK!")


         ## Polish table
         ## ............
         message("... ... retrieved data frame with ", nrow(y), " rows")
         y$entrezgene_id <- as.character(y$entrezgene_id)
         y$strand <- ifelse(y$strand=="-1", "-", y$strand)
         y$strand <- ifelse(y$strand=="1", "+", y$strand)
         y$strand <- ifelse(y$strand %in% c("+","-"), y$strand, "*")

         y_orig <- y

         ## Polish chromsome
         ## ................
         if(!is.null(keep.chromosomes)){
            message("--> Filter")
            message(" ... ... keeping only selected chromosomes: ", paste(keep.chromosomes, collapse = ", "))
            y <- y[y$chromosome_name %in% keep.chromosomes,]
            message("... ... data frame of ", nrow(y), " rows")
            }

         ## Polish transcript biotype
         ## ................
          if(!is.null(transcript.biotypes)){
            message("--> Filter")
            message(" ... ... keeping only selected transcript biotypes: ", paste(transcript.biotypes, collapse = ", "))
            y <- y[y$transcript_biotype %in% transcript.biotypes,]
            message("... ... data frame of ", nrow(y), " rows")
            }


         ## Is query table 1:n of primary query.type
         ## ..............
         if(length(unique(y[,"ensembl_gene_id"])) < nrow(y) ){
            if(collapse.to.unique){
               message("\n--> Collapsing rows/features in biomart tables")
               mult_list <- split(y[,], y[,"ensembl_gene_id"])
               u <- unlist(lapply(mult_list, nrow))
               message(" ... ... mutiple matches (1:n) found for: ", length(which(u>1)), " id's")
               message(" ... ... coercing these using pipes '|'")

               mult_list[u>1] <- lapply(mult_list[u>1], function(x){
                  # 1: for genomic location - if chromosome doesnt match - set all to NA
                     if(!all(x$chromosome_name==x$chromosome[1])){
                        x$chromosome_name=NA
                        x$start_position=NA
                        x$end_position=NA}
                     if(!all(x$strand==x$strand[1])){
                        x$strand="*"}
                     # 2: for genomic location: use the widest region
                     x$start_position <- min(x$start_position)
                     x$end_position <- max(x$end_position)
                     # 3: Pipe remaining identifiers (remove NA/blank)
                     apply(x, 2, function(y){
                        yy <- unique(y)
                        yy <- yy[!is.na(yy)]
                        yy <- yy[!yy==""]
                        paste(yy, collapse = "|")
                        })
                     })

               y_merged <- do.call("rbind", mult_list)
               if(any(duplicated(y_merged[,"ensembl_gene_id"]))) stop("FAILED GENERATING UNIQUE PRIMARY IDs")
               message("   --> OK!")
               message("... ... data frame of ", nrow(y_merged), " rows after merging rows")

               if(any(y_merged$start_position=="")){
                  message(" ... Removing primary_ids with no/multiple mapped location(s) after merging")
                  y_merged <- y_merged[!(y_merged$start_position==""),]
                  message("... ... data frame of ", nrow(y_merged))
               }
               y <- y_merged
               message("\n   --> DONE!")
               message("... ... done polishing primary biomart object\n ... ... ", nrow(y), " rows after merging rows")
            }}




         ## Create  genomic ranges object
         ## .......................
         message("\n   --> Generating GRanges object of biomart table")
         probesGR <- GRanges(seqnames = paste0("chr", y$chromosome_name), strand=as.factor(y$strand),
                ranges=IRanges(start=as.numeric(y$start_position), end=as.numeric(y$end_position)))
         elementMetadata(probesGR) <- y[,-match(c("start_position", "end_position", "chromosome_name","strand"), colnames(y))]
         genome(probesGR) <- featureId[[ paste0(x, ".genome")]]
         probesGR <- toGenomeOrder(probesGR)
         if(collapse.to.unique) names(probesGR) <- elementMetadata(probesGR)[,"ensembl_gene_id"]

         metadata(probesGR)[["query.time"]] <- timestamp()
         metadata(probesGR)[["query.id"]] <- homologBM
         metadata(probesGR)[["query.type"]] <- "ensembl_gene_id"
         metadata(probesGR)[["biomart.table"]] <- y_orig

         featureId[[x]] = probesGR


         ## CHAIN TABLE - Merge with table for primary species - Chain
         names(y)[which(names(y)=="ensembl_gene_id")] <- MartArgsHomologs[[x]]
         y_chain <- dplyr::left_join(homologBM, y, by=MartArgsHomologs[[x]])
         if(featureId$query.type!="ensembl_gene_id"){
            y_chain <- dplyr::left_join(as.data.frame(featureId[[featureId$query.species]]), y_chain, by="ensembl_gene_id")
         }
         featureId[[paste0(x,".chain")]] = as_tibble(y_chain)

      } # end x in species_homologs
      } # end if orthologs
      return(featureId)
   }



#' Function to cerate Genomic Ranges fdata object of gene coordinates using Biomart
#' previously genes2genomicRanges and probes2genomicRanges
#' @family granges
#' @family feature annotation
#' @param ids Input identfier - A vector of >1 gene ensembl ids. if 'all' then all genes in mart are retrieved
#' @param id_type what identifier type for biomart lookup - one of c("ensembl_gene_id","hgnc_symbol", "entrezgene")
#' @param collapse.to.unique If output should contain uniqe ids (per row) or if multiple rows are allowed for out. The non-collapsed table will be stored in metadata slot
#' @param primary_id If collapse.to.unique = T. This id type will be used as primary (per row). Other Ids will be piped. If multiople regions are present per primary region - The first instance will be acquired.
#' @param genome.version what genome should be used - one of c('hg19','grch38')
#' @param keep.chromosomes what chrosomosomes to keep. set to NULL if to keep all chromosomes
#' @return GRanges object with genomic start/stop for supplied identifiers (or all if 'all' option used)
#' @export
#'
features2genomicRanges = function(
   ids = "all",
   id_type = c("ensembl_gene_id","hgnc_symbol", "entrezgene_id"),
   collapse.to.unique = T,
   primary_id = NULL,
   genome.version = c('hg19','hg38'),
   keep.chromosomes = c(1:22, "X","Y","MT")
   ){

   require(biomaRt)
   require(GenomicRanges)
   require(genoset)

   ids <- as.character(ids)
   id_type = match.arg(id_type, several.ok = F)
   genome.version = match.arg(genome.version, several.ok = F)
   if(is.null(primary_id)) primary_id=id_type
   primary_id = match.arg(arg = primary_id, choices = c("ensembl_gene_id","hgnc_symbol", "entrezgene_id"), several.ok=F)
   message("... primary_id used: ", primary_id)

   useMartArgs <- list(
      hg19 = list(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org",
                      path="/biomart/martservice",dataset="hsapiens_gene_ensembl"),
      grch38 = list(biomart="ensembl", path="/biomart/martservice",dataset="hsapiens_gene_ensembl")
         )

   if(genome.version == "hg19") message("  !!! note !!! \nFor some 100's TCGA/GDC gene (ensembl id's) chr locations are not available in biomart - this since some hg38 ENSGs are newly created in hg38 (eg. HNF1B) and uses a different ENSG in hg19, try obtain using entrezgene id")

   ## Call useMart
   message("..")
   message("... calling ensembl mart for genome ", genome.version)
   my_mart <-  do.call("useMart", useMartArgs[[genome.version]])

   ## Genomic positions for ENSGs
   my_attributes <- c("ensembl_gene_id" ,"start_position", "end_position", "chromosome_name","strand", "hgnc_symbol", "entrezgene_id","transcript_biotype")
   #my_attributes <- c("ensembl_gene_id" ,"start_position", "end_position", "chromosome_name","strand")
   if(length(ids)>=1 & ids[1]!="all"){
      probe_df <- getBM(attributes=my_attributes, filters=id_type,  values=ids, mart=my_mart)
      message("... ... ... gchr locations obtained for ", nrow((probe_df)), " ensembl_gene_ids from list of ", length(unique(ids)), " ids")
      if(nrow(probe_df)==0) stop("no matching identifiers in mart - correct id_type?")
      }
   if(ids[1]=='all'){
      message("getting all ENSG gene entries")
      probe_df <- getBM(attributes=my_attributes, mart=my_mart)
      }
   message("... ... data frame of ", nrow(probe_df), " rows")
   probe_df$entrezgene_id <- as.character(probe_df$entrezgene_id)
   probe_df$strand <- ifelse(probe_df$strand=="-1", "-", probe_df$strand)
   probe_df$strand <- ifelse(probe_df$strand=="1", "+", probe_df$strand)
   probe_df$strand <- ifelse(probe_df$strand %in% c("+","-"), probe_df$strand, "*")

   ## Polish chromsome
   ## ................
      # keep only
   probe_df_orig <- probe_df
   if(!is.null(keep.chromosomes)){
      message(" ... ... keeping chromosomes: ", paste(keep.chromosomes, collapse = ", "))
      probe_df <- probe_df[probe_df$chromosome_name %in% keep.chromosomes,]
      message("... ... data frame of ", nrow(probe_df), " rows")
      }

   ## Remove any NA/blank in primary identifier
   ## ................
   if(any(is.na(probe_df[,primary_id])) | any(probe_df[,primary_id]=="") ){
      message(" ... Removing NA/missing primary_ids")
      probe_df <- probe_df[!c(is.na(probe_df[,primary_id]) | probe_df[,primary_id]==""),]
      message("... ... data frame of ", nrow(probe_df), " rows after filtering")
      }


   ## Keep one row per primary (input) identifier
   ## ...................
      # if primary idntifier is not ENSG - then multiple regions may occur per enetry
      # Other identifiers are piped '|'
      # The primary identifier is kept but alternate regions
   if(collapse.to.unique){
      if(any(duplicated(probe_df[,primary_id]))){
         message("... Duplicated primary_ids exist. collapsing dataframe on these. If genomic regions differ - the widest region is given.")
         probe_df_split <- split(probe_df[ ,colnames(probe_df)], probe_df[ ,primary_id])
         u <- unlist(lapply(probe_df_split, nrow))
         # for all with duplicated entries run....
         probe_df_split[u>1] <- lapply(probe_df_split[u>1], function(x){
               # 1: for genomic location - if chromosome doesnt match - set all to NA
               if(!all(x$chromosome_name==x$chromosome[1])){
                  x$chromosome_name=NA
                  x$start_position=NA
                  x$end_position=NA}
               if(!all(x$strand==x$strand[1])){
                  x$strand="*"}
               # 2: for genomic location: use the widest region
               x$start_position <- min(x$start_position)
               x$end_position <- max(x$end_position)
               # 3: Pipe remaining identifiers (remove NA/blank)
               apply(x, 2, function(y){
                  yy <- unique(y)
                  yy <- yy[!is.na(yy)]
                  yy <- yy[!yy==""]
                  paste(yy, collapse = "|")
                  })
               })

         probe_df_merged <- do.call("rbind", probe_df_split)
         if(any(duplicated(probe_df_merged[,primary_id]))) stop("FAILED GENERATING UNIQUE PRIMARY IDs")
         message("... ... data frame of ", nrow(probe_df_merged), " rows after merging rows")

         if(any(probe_df_merged$start_position=="")){
            message(" ... Removing primary_ids with no/multiple mapped location(s) after merging")
            probe_df_merged <- probe_df_merged[!(probe_df_merged$start_position==""),]
            message("... ... data frame of ", nrow(probe_df_merged))
         }

         probe_df <- probe_df_merged
         } # end if any duplicated
      } # end if collapse


   ## Identifiers with no (multiple) genomic location after merge
   ## .......................


   ## Create  genomic ranges object
   ## .......................
   probesGR <- GRanges(seqnames = paste0("chr", probe_df$chromosome_name), strand=as.factor(probe_df$strand),
          ranges=IRanges(start=as.numeric(probe_df$start_position), end=as.numeric(probe_df$end_position)))
   elementMetadata(probesGR) <- probe_df[,-match(c("start_position", "end_position", "chromosome_name","strand"), colnames(probe_df))]
   genome(probesGR) <- genome.version
   probesGR <- toGenomeOrder(probesGR)
   if(collapse.to.unique) names(probesGR) <- elementMetadata(probesGR)[,primary_id]

   metadata(probesGR)[["biomart"]] <- timestamp()
   metadata(probesGR)[["ids"]] <- ids
   metadata(probesGR)[["primary_id"]] <- primary_id
   metadata(probesGR)[["biomart.original.table"]] <- probe_df_orig
   return(probesGR)
}





#' Get homo sapiens gene symbols from a vector of entrez gene ids using appropriate org.eg.db package.
#' @family feature annotation
#' @family old functions
#' @param entrez_id, a vector of entrez ids. May include NAs and duplicates.
#' @param keep.all.symbols, if to include entrez for which no sympol was found in output table.
#' @param species, define species (e.g. "Hs", "Mm", "Rn", "Dm"...). Defines what org.eg.db package that will be loaded.
#' @return a feature data table with valid entrez ids (ENTREZID) and their respective gene symbols (SYMBOL)
#' @export
entrez2symbol<-function(entrez.id, keep.all.symbols=F, species="Hs"){
   ## species defines organism, e.g. "Hs" do use the human Entrez db ("org.Hs.eg.db")
   requireNamespace("org.Hs.eg.db")
   requireNamespace("AnnotationDbi")

   species.package<-paste("org.",species,".eg.db",sep="")
  print(paste("Organism =",species,"- Requires", species.package))
  do.call(require, args=AnnotationDbi::as.list(species.package))
  species.sub<-paste("org.",species,".eg",sep="")

	# Check no duplicates
	u<-which(duplicated(entrez.id)==T)
	if(length(u)>0){
		entrez.id<-entrez.id[-u]
		cat("removing",length(u),"duplicate id","\n")
		}
	# Gene Symbols
	assign("x", eval(parse(text=paste(species.sub,"SYMBOL",sep=""))))
  #x <- org.Hs.egSYMBOL
	mapped_genes <- mappedkeys(x)
	xx <- AnnotationDbi::as.list(x[mapped_genes])
	y<-xx[as.character(entrez.id)]

	yy<-unlist(y, recursive=2)
	if(length(y)!=length(yy)){
      print(paste(length(y)-length(yy),"accession number(s) not found"))
	   }


	out.table<-data.frame(ENTREZID=names(yy), SYMBOL=yy, row.names=names(yy), stringsAsFactors=F)

  if(keep.all.symmbols){
    yy <- unlist(y)
    out.table<-data.frame(ENTREZID=names(yy), SYMBOL=yy, row.names=names(yy), stringsAsFactors=F)
    ot<-merge(data.frame(ENTREZID=yy, stringsAsFactors=F), out.table, by.x="ENTREZID", by.y="ENTREZID", all.x=T, sort=F)
    u<-match(entrez.id, ot$ENTREZID)
    ot<-ot[u,]
    out.table=ot
    }

   out.table$ENTREZID <- sapply(out.table$ENTREZID, as.character)
  return(out.table)
	} # end entrez2symbol





#' Get homo sapiens chromosome locations from a vector of entrez gene ids using org.Hs.eg.db pack
#' @family feature annotation
#' @family old functions
#' @param entrez_id, a vector of entrez ids. May include NAs and duplicates.
#' @return a feature data table with valid entrez ids and their respective chr
#' @export
entrez2chr <- function(entrez_id) {
    requireNamespace("org.Hs.eg.db")
   requireNamespace("AnnotationDbi")

    # Check no duplicates
    u <- which(duplicated(entrez_id) == T)
    if (length(u) > 0) {
        entrez_id <- entrez_id[-u]
        cat("removing", length(u), "duplicate id", "\n")
    }

    # CHR
    x <- org.Hs.eg.db::org.Hs.egCHR
    mapped_genes <- AnnotationDbi::mappedkeys(x)
    xx <- AnnotationDbi::as.list(x[mapped_genes])
    y <- xx[as.character(entrez_id)]

    yy <- unlist(y, recursive = 2)
    if (length(y) != length(yy)) {
        print(paste(length(y) - length(yy), "accession number(s) not found"))
    }

    if (length(which(duplicated(names(yy)) == T)) > 0) {
        cat("\n ... ", length(which(duplicated(names(yy)) == T)), "entries with multiple chrs found\n")
        yy <- yy[!duplicated(names(yy))]
    }

    out.table <- data.frame(ENTREZID = names(yy), CHR = yy, row.names = names(yy), stringsAsFactors = F)
    return(out.table)
}  # end entrez2chr




#' update Illumina v4 probe ids using the IlluminaHumanv4.db package
#' @family feature annotation
#' @param probe_ids, vector containing valid Illumina probe IDs
#' @return a data frame with feature data
#' @export

illuminav4_2fdata <- function(probe_ids){
   requireNamespace("illuminaHumanv4.db", "methods")
   requireNamespace("AnnotationDbi")
   requireNamespace("DBI")

   cat('\n... Loading package illuminaHumanv4.db ...\n')
   fdata <- illuminaHumanv4fullReannotation()
   fdata <- fdata[,c('IlluminaID','ProbeQuality','CodingZone','GenomicLocation','ProbeSequence','ReporterGroupName')]
   u <- match(probe_ids, fdata$IlluminaID)
   if(any(is.na(u))) stop('not all probe ids are present in illumina v4 package')

   # ENTREZID
   x<-illuminaHumanv4ENTREZREANNOTATED
   mapped_probes <- mappedkeys(x)
   xx <- AnnotationDbi::as.list(x[mapped_probes])
   if(length(xx) > 0) {
      y <- data.frame(IlluminaID=names(unlist(xx[probe_ids])),ENTREZID=unlist(xx[probe_ids]))
      fdata <- merge(y, fdata, by.x='IlluminaID', by.y='IlluminaID', all.y=T)
      }

   # SYMBOL
   x<- illuminaHumanv4SYMBOLREANNOTATED
   mapped_probes <- mappedkeys(x)
   xx <- AnnotationDbi::as.list(x[mapped_probes])
   if(length(xx) > 0) {
      y <- data.frame(IlluminaID=names(unlist(xx[probe_ids])),SYMBOL=unlist(xx[probe_ids]))
      fdata<-merge(y, fdata, by.x='IlluminaID', by.y='IlluminaID', all.y=T)
      }

   # ENSEMBL
   x<-illuminaHumanv4ENSEMBLREANNOTATED
   mapped_probes <- mappedkeys(x)
   xx <- AnnotationDbi::as.list(x[mapped_probes])
   if(length(xx) > 0) {
      y <- data.frame(IlluminaID=names(unlist(xx[probe_ids])),ENSEMBL=unlist(xx[probe_ids]))
      fdata=merge(y, fdata, by.x='IlluminaID', by.y='IlluminaID', all.y=T)
      }

   # Build table and return
   u <- match(probe_ids, fdata$IlluminaID)
   fdata.out <- fdata[u,]
   fdata.out <- data.frame(sapply(fdata.out, as.character), stringsAsFactors=F)
   return(fdata.out)
   } # end illuminav4_2fdata



#' use locally stored NCBI homologene data (homologene.data.build68.Rdata) to find inter-species homologs gene IDs
#' @family feature annotation
#' @family old functions
#' @param entrez.id, entrez_id
#' @param from.tax, from what species (Mm, Hs, Rn, Dm)
#' @param to.tax,  from what species (Mm, Hs, Rn, Dm)
#' @param keep.all.hids, logical if to keep
#' @return a data frame with from and to entrezid and homolgy ID - HID
#' @export
homoloGene<-function(entrez.id, from.tax="Mm", to.tax="Hs", keep.all.hids=T){

   tax.id.list<-list(
         Mm="10090",
         Hs="9606",
         Rn="10116",
         Dm="7227"
         )

    match.arg(arg = from.tax, choices = names(tax.id.list))
   match.arg(arg = to.tax, choices = names(tax.id.list))

  	print(paste("..finding homoloGenes in",to.tax,"for",from.tax,"ENTREZIDs"))


	load(file="/Users/david/R/R_FUNCTIONS_AND_REFERENCE_FILES/homologene.data.build68.Rdata")
	print(paste("...loading file 'homologene.data.build68.Rdata'"))

	entrez.id<-as.integer(entrez.id)
	if(length(which(is.na(entrez.id)>0))){
	print(paste("...removing",table(!is.na(entrez.id))[1],"NA-entries from input ENTREZIDs ::", table(!is.na(entrez.id))[2],"ENTREZIDs kept"))
	entrez.id<-entrez.id[!is.na(entrez.id)]}

	if(length(which(duplicated(entrez.id)==T))){
	print(paste("...removing",table(duplicated(entrez.id))[2],"duplicated entries from input ENTREZIDs ::", table(duplicated(entrez.id))[1],"ENTREZIDs kept"))
	entrez.id<-entrez.id[!duplicated(entrez.id)]}

	entrez.table<-data.frame(entrez.id, stringsAsFactors=F)

	from.hid<-homologene[homologene$Taxonomy_ID==tax.id.list[from.tax],c(1,3)]
	hid.table<-merge(entrez.table, from.hid, by.x="entrez.id", by.y="ENTREZID", all.x=T, sort=F)
	colnames(hid.table)[1]<-paste("ENTREZID",from.tax,sep=".")

	to.hid<-homologene[homologene$Taxonomy_ID==tax.id.list[to.tax],c(3,1)]
	out.table<-merge(hid.table, to.hid, by.x="HID",by.y="HID", sort=F)
	out.table[] <- lapply(out.table, as.character)
	print(paste("...returning list of",nrow(out.table),"matching ENTREZIDs"))

   cat('\n ... out of', length(entrez.id), 'unique Entrez IDs from tax', from.tax)
   cat('\n ...... ', length(unique(out.table$HID)), 'were matched to a HID')
   cat('\n ...... and these correspond to', length(unique(out.table$ENTREZID)), to.tax, 'Entrez IDs')

   dup.hid.mapping.ids = unique(out.table$HID[which(duplicated(out.table$HID)==T)])
   cat('\n...... ...',length(dup.hid.mapping.ids), 'have more than 1 HID::Entrez mapping')

   if(!keep.all.hids){
      cat('\n...... ... removing these')
      u = which(out.table$HID %in% dup.hid.mapping.ids)
      if(length(u)>0) out.table = out.table[-u,]}

   cat('\n...... Returning table with', nrow(out.table), 'Entrez :: HID :: Entrez connections')

   out.table[ ] <- lapply(out.table, as.character)
	return(out.table)
	} # end homologene
