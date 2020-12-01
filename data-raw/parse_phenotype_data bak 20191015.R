   require(reshape2)
   require(dplyr)
   options(dplyr.width = Inf)
   require(tibble)
   require(tidyverse)
   require(tidyr)


   ## ::::::::::::::::::::::::::::::::::
   ##       CURATE PDATA
   ## ::::::::::::::::::::::::::::::::::


   ##   STEP 1 & 2: PANCAN PDATA - including pancan survival
   ## ------------------------------------


   ## TCGA PanCan aliquot info - based on pancan 'quality' table
      tcga_pancan_samples <- readRDS("~/RESOURCES/TCGA_PANCAN/Processed_data/pantcga_sampleAnnotList.rds")
   # names(tcga_pancan_samples)
     tcga_pancan_samples$quality  # A tibble: 79,286 x 17

    clinical_pan <- readRDS(file="~/RESOURCES/TCGA_PANCAN/Processed_data/PanCan_Clinical_CDR_n11141.rds")
      table(clinical_pan$vital_status)
      table(clinical_pan$OS)

      clinical_pan <- clinical_pan %>%
         dplyr::rename(patient_barcode = bcr_patient_barcode)
         #filter(patient_barcode %in% my_df$patient_barcode) %>%
         # filter(vital_status != "[Discrepancy]")

      tcga_surv_tab <- clinical_pan


   ## Match up with PAN tcga qc annotations
   ## ::::::::::::::::::::::::
   pdata <- as_tibble(tcga_pancan_samples$quality) %>% # A tibble: 79,286 x 28
      dplyr::left_join(tcga_surv_tab, by="patient_barcode") %>%
      dplyr::select(sample_id, patient_barcode, cancer.type, aliquot_barcode, Do_not_use, birth_days_to, age_at_initial_pathologic_diagnosis, gender,
         race, ajcc_pathologic_tumor_stage, clinical_stage, histological_type, histological_grade, initial_pathologic_dx_year, menopause_status, OS, OS.time, DSS, DSS.time, DFI, DFI.time, PFI,
         PFI.time, vital_status,tumor_status,last_contact_days_to,death_days_to,cause_of_death, everything())
   table(pdata$cancer.type)
   pdata[pdata$cancer.type=="",] # 10 samples have missing cancer.type, remove these
   pdata <- pdata %>%
      #mutate(cancer.type = if_else(cancer.type=="", as.character(fu_group), cancer.type)) %>%
      dplyr::filter(cancer.type!="") %>%
      dplyr::arrange(cancer.type)
   table(is.na(pdata$cancer.type))






   ## iCluster
   ## ......
   str(tcga_pancan_samples)
   pdata <- pdata %>%  # A tibble: 79,276 x 51
      dplyr::left_join(tcga_pancan_samples$icluster) %>%
      mutate(iCluster = as.character(iCluster)) %>%
      mutate(iCluster = paste0("cluster_",iCluster)) %>%
      mutate(iCluster = if_else(grepl("_NT",sample_id), "", iCluster))
   table(pdata$iCluster)





   ## MASTER CALLS
   ##........................
      # (based on cnv sample- may not be matched aliquot)
      # >900 samples have been called in legacy data
      # These lack proper aliquot_barcodes -  instead like: KIRC-TCGA-BP-4347-Tumor
      # Can be both genome wide snp and dna seq analyses

      Hmisc::describe(tcga_pancan_samples$master_calls)
            #       Value          legacy_call legacy_maf_call        maf_call        snp_call
            # Frequency              558              35             146              56
            # Proportion           0.702           0.044           0.184           0.070

      # BUT!! no patient barcodes are duplicated - therefore merge as is on patient barcode
      master_calls <- tcga_pancan_samples$master_calls
      patient_barcodes <- unlist(lapply(master_calls$patient_barcode, function(x){
         paste(strsplit(x, split = "-")[[1]][1:3], collapse = "-")
         }))
      master_calls <- master_calls %>% mutate(patient_barcode = patient_barcodes)
      table(duplicated(master_calls$patient_barcode)) ##  62 duplicated patients
      master_calls <- master_calls[!duplicated(master_calls$patient_barcode),]

      pdata <- pdata %>%
         dplyr::left_join(
            master_calls %>% dplyr::select(-c(sample_id, aliquot_barcode, sample_barcode, vial_barcode, portion_barcode, analyte_barcode, platform, solution, cancer.type)))
      colnames(pdata)
      # Hmisc::describe(pdata)
      # A tibble: 79,276 x 58


   ## CIBERSORT
   ## -----------------
   Hmisc::describe(tcga_pancan_samples$cibersort)
   str(tcga_pancan_samples$cibersort)
      #    aliquot_barcode
      #     n  missing distinct
      # 11373        0    11312
      #
      #    sample_id
      #     n  missing distinct
      # 10469      904    10423
      #    aliquot_barcode
      #     n  missing distinct
      # 11373        0    11312
   tcga_pancan_samples$cibersort[is.na(tcga_pancan_samples$cibersort$sample_id), ]
   # Fix missing patient / sample ids
   cibersort <- tcga_pancan_samples$cibersort

   cibersort <- cibersort %>%
      mutate(sample_id_add = pdata$sample_id[match(aliquot_barcode, pdata$aliquot_barcode)]) %>%
      mutate(sample_id = if_else(is.na(sample_id) , sample_id_add , sample_id)) %>%
      dplyr::select(-sample_id_add)
   Hmisc::describe(cibersort) # A tibble: 11,373 x 34
   cibersort <- cibersort %>% filter(!is.na(sample_id))
   table(duplicated(cibersort$sample_id)) ## 64 duplicated sample ids - remove these
   cibersort <- cibersort[!duplicated(cibersort$sample_id), ]

   # merge on sample_id after removing na-patient ids
      pdata <- pdata %>%
         dplyr::left_join(cibersort %>%  dplyr::select(-c(patient_barcode, aliquot_barcode, sample_barcode, vial_barcode, portion_barcode, analyte_barcode, platform, cancer.type)))
      # A tibble: 79,276 x 83





   ## :::::::::::::::::
   ##    KEEP OBJECT
   ## ::::::::::::::::

   pdata_step1 <- pdata
   rm(tcga_pancan_samples, pdata,  clinical_pan)








   ## ::::::::::::::::::::::::::::::::::-
   ##   PDATA STEP 2 - specific taxonomies mutations etc
   ## ::::::::::::::::::::::::::::::::::


   ## get histoology and molecular grouping information from previous analyses including:
   # * Taxonomy paper (Chen 2016), Includes annotations on exclusion of samples for the three RCC TCGA marker papers
   # * RCC FGFR2 isoform switch paper ()
   # * ccRCC histology review (Mod Path Favazza 2017)
   # REMOVE OUTLIERS (Lindgren 2017)
   # ----------------
   pdata <- pdata_step1

   table(is.na(pdata$sample_id))

   pdata_chen <- read.delim("~/R_PKG/dlfoo2/data-raw/RCC_taxonomy_Chen2016_TableS2.txt", as.is=T)
   chen_cols <- c("sample_id","Path.review.notes","explicitly.removed.from.marker.paper",
            "VHL_silenced","TFE3.TFEB.fusion","TFE3.fusion","SWI.SNF.module","PI3K.module","TP53","VHL","SETD2","BAP1","PBRM1","KDM5C","NF2","PTEN","ARID1A","MICALCL",
            "STAG2", "SLC1A3", "CDKN1A",  "MTOR",   "MET", "SMARCB1", "TCEB1", "NFE2L2", "PIK3CA", "MLL3",
            "OX.PHOS","Cell.Cycle", "hypoxia","YAP1", "NRF2.ARE")
   pdata <- dplyr::left_join(pdata, pdata_chen[,chen_cols], by="sample_id") # A tibble: 10,223 x 95 ;; # A tibble: 79,276 x 115
   pdata

   ## Fix 0/1 as chacacter for binary annotations (to YES/NO)
   head(pdata_chen[,chen_cols])
   table(pdata_chen$VHL_silenced)
   chen_cols2 <- c("TP53","VHL","SETD2","BAP1","PBRM1","KDM5C","NF2","PTEN","ARID1A","MICALCL",
            "STAG2", "SLC1A3", "CDKN1A",  "MTOR",   "MET", "SMARCB1", "TCEB1", "NFE2L2", "PIK3CA", "MLL3")

   myfoo <- function(x){
      x[x=="1"] <- "YES"
      x[x=="0"] <- "NO"
      return(x)
   }

   pdata <- pdata %>% mutate_at(chen_cols2, myfoo)
   table(pdata$TP53)
   str(pdata)
   # Hmisc::describe(pdata)

   pdata_old <- read.delim("~/R_PKG/dlfoo2data/data-raw/pdata_utc_1435set_edit_v20170830.txt", as.is=T, na.strings = c("","NA"))
   str(pdata_old)
   pdata_old_cols <- c("sample_id","EOS_davis_2014","SARC_davis_2014","FGFR2_IIIb_Zhao","histology_davis_2014","variant_histology_comment","uro_taxonomy","chen_taxonomy","dl_n_taxonomy")
   table(duplicated(pdata_old[,pdata_old_cols]$sample_id))
   table(duplicated(pdata$sample_id))
   pdata <- dplyr::left_join(pdata, as_tibble(pdata_old[,pdata_old_cols]), by="sample_id")
   str(pdata) # str(pdata) # Classes ‘tbl_df’, ‘tbl’ and 'data.frame':	79276 obs. of  123 variables:

   table(pdata$cancer.type, useNA="always")
   temp_annot <- sub("(\\_.*?)\\_", "\\1-", pdata$sample_id)
   temp_annot <- sub("-.*","",temp_annot)
   temp_annot <- sub("_","-",temp_annot)

   pdata <- pdata %>%
      mutate(sample_type2 =   temp_annot) %>%
      mutate(uro_taxonomy = replace(uro_taxonomy, sample_type2=="BLCA-NT", "blca_normal")) %>%
      #dplyr::select(-c(dl_n_taxonomy2, taxonomy, taxonomy2)) %>%
      mutate(chen_taxonomy = as.character(factor(chen_taxonomy, exclude=c("","NA")))) %>%
      mutate(uro_taxonomy = as.character(factor(uro_taxonomy, exclude=c("","NA")))) %>%
      mutate(dl_n_taxonomy = as.character(factor(dl_n_taxonomy, exclude=c("","NA")))) %>%
      mutate(taxonomy_published = if_else(is.na(chen_taxonomy), uro_taxonomy, chen_taxonomy)) %>%
      mutate(taxonomy_published = if_else(is.na(taxonomy_published), dl_n_taxonomy, taxonomy_published)) %>%
      #dplyr::select(-c(chen_taxonomy, uro_taxonomy, dl_n_taxonomy)) %>%
      mutate(taxonomy_published = dplyr::recode(taxonomy_published, "NK-1"="crtx", "NK-2"="med", "NK-3"="crtmed","NK-4"="infl1","NK-5"="infl2"))
      # A tibble: 10,223 x 104

   table(pdata$taxonomy_published)
   table(pdata$sample_type2)

   pdata_step2 <- pdata







   ## :::::::::::::::::::::::::::::
   ##    ADD UPDATE DL TAXONOMY 2019
   ## :::::::::::::::::::::::::::::
   retax <- as_tibble(read.delim("~/PROJECTS/RCC_ccpRCC_2019/reTaxonomy/reTaxonomy_81samples_20190929.txt", as.is=T))
   retax[retax$sample_id=="KIRC_TP_A8OX",]
   str(retax)
   table(retax$reTax)
   sum(table(retax$reTax))

   pdata <- pdata %>%
      mutate(taxonomy_published = factor(taxonomy_published, levels=c(
         "crtx","crtmed","med","infl1","infl2",
         "CC-e.2","CC-e.1","CC-e.3",
         "P-e.1a","P-e.1b", "P-e.2","P.CIMP-e",
         "Ch-e","mixed",
         "blca_normal","Urobasal","GU","Basal/SCCL","SC/NE","Mes-Inf","Inf_other"))) %>%
      mutate(x = as.character(taxonomy_published)) %>%
      # mutate(x = if_else(is.na(taxonomy_chen), "NA", x)) %>%
      mutate(x = recode(x, "CC-e.1"="ccRCC_e1", "CC-e.2"="ccRCC_e2", "CC-e.3"="ccRCC_e3")) %>%
      mutate(x = recode(x, "Ch-e"="chRCC", "P-e.1a"="pRCC_e1a", "P-e.1b"="pRCC_e1b")) %>%
      mutate(x = recode(x, "P-e.2"="pRCC_e2", "P.CIMP-e"="pCIMP")) %>%
      mutate(x = if_else(as.character(sample_type2)=="BLCA-NT", "bladder_n", x)) %>%
      dplyr::left_join(retax %>% dplyr::select(-taxonomy_published, -iCluster), by="sample_id") %>%
      mutate(x = if_else(!is.na(reTax), reTax, x)) %>%
      mutate(y = x) %>%
      mutate(y = recode(y, "crtmed"="kidney_n", "crtx"="kidney_n", "med"="kidney_n", "infl1"="kidney_n", "infl2"="kidney_n")) %>%
      mutate(y = recode(y, "ccRCC_e2"="ccRCC", "ccRCC_e1"="ccRCC", "ccRCC_e3"="ccRCC")) %>%
      mutate(y = recode(y, "pRCC_e1b"="pRCC_e1", "pRCC_e1a"="pRCC_e1")) %>%
            mutate(y = gsub("outlier_.*","outlier",y)) %>%
      mutate(x = factor(x, levels=c("crtx","crtmed","med","infl1","infl2","ccRCC_e2","ccRCC_e1","ccRCC_e3", "pRCC_e1a","pRCC_e1b", "pRCC_e2",
         "pCIMP","dCIMP","ccpRCC","chRCC","chONC","pONC","mesRCC", "metanephric","mtscRCC",
         "outlier", "outlier_atypic","outlier_tfe","outlier_cystic","outlier_cimp","outlier_inf",
         "bladder_n","Urobasal","GU","Basal/SCCL","SC/NE","Mes-Inf","Inf_other"))) %>%
      mutate(y = factor(y, levels=c("kidney_n","ccRCC", "pRCC_e1", "pRCC_e2",
         "pCIMP","dCIMP","ccpRCC","chRCC","chONC","pONC","mesRCC","metanephric","mtscRCC","outlier",
         "bladder_n","Urobasal","GU","Basal/SCCL","SC/NE","Mes-Inf","Inf_other"))) %>%
      dplyr::rename(tax_dl = x) %>%
      dplyr::rename(tax_simp = y) %>%
      dplyr::select(-reTax)

   table(pdata$tax_simp[which(pdata$platform=="IlluminaHiSeq_RNASeqV2")])
   pdata[which(pdata$tax_simp=="outlier"),]
   data.frame(pdata[pdata$sample_id=="KIRC_TP_A8OX",])

   sum(table(pdata$taxonomy_published))
   sum(table(pdata$tax_dl))
   sum(table(pdata$tax_simp))



   ##  pdata ConsensuClsutering RCC k25
   ## ------------------------------
    my.k <- 25
      # my.cc <- dlfoo2::ConsensusCluster_process(cc.dir = "~/PROJECTS/RCC_ccpRCC_2019/ConsensusClust/ConsensusCluster_rcc_t_5000_880set_5000genes_reps.1000_maxK.25_innerLinkage.ward_sd.rank.5000/", k = my.k)
      # #stopifnot(identical(my.cc$sample.names, pdata$sample_id))
      # my_tab <- data.frame(sample_id=my.cc$sample.names, cc_rcc_t_k25=paste0("cc_", my.cc$cc), stringsAsFactors = F)
      # apa <- pdata
      #
      # pdata <- pdata %>%
      #    left_join(my_tab, by="sample_id") %>%
      #    mutate(cc_rcc_t_k25 = factor(cc_rcc_t_k25, levels=paste0("cc_",1:25)))
      # sum(table(pdata$taxonomy_published, pdata$cc_rcc_t_k25))
      # sum(table(pdata$cc_rcc_t_k25))
      #
      #
      # #pdata <- dlfoo2data::pdata_panCanGex
      #
      # cc_k25_grouping <- list(
      #    pRCC_a = paste0("cc_",c(23,25,16)),
      #    pRCC_b = paste0("cc_",c(20,24)),
      #    pRCC_c = paste0("cc_",c(9)),
      #    pRCC_cimp = paste0("cc_",c(22)),
      #    ccRCC_a = paste0("cc_",c(18,14,6,19,12)),
      #    ccRCC_b = paste0("cc_",c(10,17,8,21,15)),
      #    ccRCC_c = paste0("cc_",c(7,13,11,5)),
      #    chRCC = paste0("cc_",c(1))
      #    )
      # cc_k25_grouping <- unlist(cc_k25_grouping)
      # x <- names(cc_k25_grouping)
      # names(x) <- cc_k25_grouping

   cc_gex <- dlfoo2::ConsensusCluster_process("~/PROJECTS/RCC_ccpRCC_2019/ConsensusClust/ConsensusCluster_rcc_t_5000_880set_5000genes_reps.1000_maxK.25_innerLinkage.ward_sd.rank.5000/", k = 25)
    #stopifnot(identical(my.cc$sample.names, pdata$sample_id))

   cc_grouping_list <- list(
      cc_gex = list(
         cc = cc_gex,
         cc_ordered = unique(paste0("cc_",cc_gex$cc[cc_gex$i]))
      )
   )
   cc_grouping_list$cc_gex$cc_groups = list(
            chRCC = cc_grouping_list$cc_gex$cc_ordered[1],
            pRCC_a = cc_grouping_list$cc_gex$cc_ordered[2:4],
            ccRCC_a = cc_grouping_list$cc_gex$cc_ordered[5:9],
            pRCC_b = cc_grouping_list$cc_gex$cc_ordered[10:11],
            ccRCC_c = cc_grouping_list$cc_gex$cc_ordered[12:15],
            pRCC_c = cc_grouping_list$cc_gex$cc_ordered[16],
            ccRCC_b2 = cc_grouping_list$cc_gex$cc_ordered[17:20],
            oncRCC = cc_grouping_list$cc_gex$cc_ordered[21],
            ccRCC_b1 = cc_grouping_list$cc_gex$cc_ordered[22],
            ccpRCC = cc_grouping_list$cc_gex$cc_ordered[23],
            dCIMP = cc_grouping_list$cc_gex$cc_ordered[24],
            pCIMP = cc_grouping_list$cc_gex$cc_ordered[25]
         )
   x <- rep(names(cc_grouping_list$cc_gex$cc_groups), lapply(cc_grouping_list$cc_gex$cc_groups, length))
   names(x) <- unlist(cc_grouping_list$cc_gex$cc_ordered)


   my_tab <- data.frame(sample_id=cc_gex$sample.names, cc_rcc_t_k25=paste0("cc_", cc_gex$cc), stringsAsFactors = F)
   my_tab$cc_rcc_t_k25 <- factor(my_tab$cc_rcc_t_k25, levels=names(x))
   my_tab$cc_rcc_t_k25_group <-  x[as.character(my_tab$cc_rcc_t_k25)]


   pdata2 <- pdata %>%
         left_join(my_tab, by="sample_id") %>%
         mutate(cc_rcc_t_k25 = factor(cc_rcc_t_k25, levels=names(x)))

   table(pdata2$cc_rcc_t_k25)
   sum(table(pdata2$taxonomy_published, pdata2$cc_rcc_t_k25))
   sum(table(pdata2$cc_rcc_t_k25))
   apa <- pdata2 %>% dplyr::filter(platform== "IlluminaHiSeq_RNASeqV2")
   table(apa$cc_rcc_t_k25, apa$tax_simp)


   ##
   pdata3 <- pdata2 %>%
      mutate(tax_simp2 = as.character(cc_rcc_t_k25_group)) %>%
      #mutate(tax_simp2 = if_else(as.character(cc_rcc_t_k25) %in% names(x), as.character(cc_k25_grouping), tax_simp2)) %>%
      mutate(tax_simp2 = if_else(as.character(tax_simp) %in% c("chONC","pONC","metanephric","mtscRCC","mesRCC","outlier"), as.character(tax_simp), tax_simp2)) %>%
      mutate(tax_simp2 = if_else(is.na(tax_simp2), as.character(tax_simp), tax_simp2)) %>%
      mutate(tax_simp2 = factor(tax_simp2, levels=c(
            "kidney_n","ccRCC_a", "ccRCC_b1", "ccRCC_b2", "ccRCC_c",
            "pRCC_a","pRCC_b","pRCC_c","pRCC_cimp", "pCIMP","dCIMP","ccpRCC","chRCC","chONC","pONC","mesRCC","metanephric","mtscRCC","outlier",
            "bladder_n","Urobasal","GU","Basal/SCCL","SC/NE","Mes-Inf","Inf_other")))

   table(pdata3$tax_simp2, pdata3$tax_simp)

   # %>%
   #    mutate(tax_simp2 = if_else(tax_simp2 %in% cc_k25_grouping, x[tax_simp2], tax_simp2)) %>%
   #    mutate(tax_simp2 = if_else(tax_simp %in% c("dCIMP","ccpRCC","chONC","pONC","mesRCC","metanephric","mtscRCC","outlier"), as.character(tax_simp), tax_simp2)) %>%
   #    mutate(tax_simp2 = factor(tax_simp2, levels=c("kidney_n",paste0("ccRCC_a",1:5), paste0("ccRCC_b",1:5), paste0("ccRCC_c",1:4),
   #          paste0("pRCC_a",1:3), paste0("pRCC_b",1:2),
   #          "pRCC_c","pRCC_cimp", "pCIMP","dCIMP","ccpRCC","chRCC","chONC","pONC","mesRCC","metanephric","mtscRCC","outlier",
   #          "bladder_n","Urobasal","GU","Basal/SCCL","SC/NE","Mes-Inf","Inf_other")))
   # table(pdata3$tax_simp2)
   # #pdata <- apa2

   pdata_step3 <- pdata3
   pdata <- pdata3
   table(pdata$tax_simp2)

   ## :::::::::::::::::::::::::::::
   ## Add pdata from pan-tcga studies (based on dna barcodes ....)
   ## :::::::::::::::::::::::::::::
   annotations_pan <- readRDS("~/RESOURCES/TCGA_PANCAN/Processed_data/pantcga_sampleAnnotList.rds")$leuk_estimates
   table(duplicated(annotations_pan$sample_id))
   table(duplicated(pdata$portion_barcode))
   table(duplicated(annotations_pan$portion_barcode))
   table(is.na(pdata$portion_barcode))

   u <- match(pdata$portion_barcode, annotations_pan$portion_barcode)
   table(duplicated(u[!is.na(u)]))

   pdata <- pdata %>% mutate(leuk_estimate = annotations_pan$leuk_estimate[u])
   pdata # # A tibble: 79,276 x 132
   data.frame(pdata[pdata$sample_id=="KIRC_TP_A8OX",])

   # table(pdata$iCluster, pdata$sample_type2)








   ## :::::::::::::::::::::::::::::
   ##  Save data
   ## :::::::::::::::::::::::::::::
   pdata_panCanFull <- pdata
   table(pdata_panCanFull$tax_simp2)
   usethis::use_data(pdata_panCanFull, overwrite = T)







