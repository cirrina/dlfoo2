

   ## ::::::::::::::::::::::::::::::::::
   ##       CURATE PDATA
   ## ::::::::::::::::::::::::::::::::::


   ##   STEP 1 & 2: PANCAN PDATA - including pancan survival
   ## ------------------------------------


   ## TCGA PanCan aliquot info - based on pancan 'quality' table
      tcga_pancan_samples <- readRDS("~/RESOURCES/TCGA_PANCAN/Processed_data/pantcga_sampleAnnotList.rds")
   # names(tcga_pancan_samples)
     tcga_pancan_samples$quality  # A tibble: 79,286 x 17


   ## Match up w tcga pan clinical data :: Old variant
   ## NOTE!!! :: prior 20190920 clinical file  clinical_PANCAN_patient_with_followup.tsv was used (below)
   ## :::::::::::::::::
      ## This should be replaced with the TCGA-Clinical Data Resource (CDR) Outcome file.
            # clinical_pan <- readRDS(file = "/Users/david/RESOURCES/TCGA_PANCAN/Processed_data/pantcga_clinical.rds")
            # table(clinical_pan$vital_status)
            # clinical_pan <- clinical_pan %>%
            #    dplyr::rename(patient_barcode = bcr_patient_barcode) %>%
            #    #filter(patient_barcode %in% my_df$patient_barcode) %>%
            #    filter(vital_status != "[Discrepancy]")
            #
            # tcga_surv_tab <- clinical_pan %>% # A tibble: 10,901 x 749
            #             mutate(fu_group = factor(acronym)) %>%
            #             mutate(fu_status = as.numeric(recode(vital_status, "Alive"="0", "Dead"="1"))) %>%
            #             mutate(fu_time = as.numeric(if_else(fu_status==1, days_to_death, days_to_last_followup))) %>%
            #             filter(fu_time!="NA")
      clinical_pan <- readRDS(file="~/RESOURCES/TCGA_PANCAN/Processed_data/PanCan_Clinical_CDR_n11141.rds")
      table(clinical_pan$vital_status)
      table(clinical_pan$OS)

      clinical_pan <- clinical_pan %>%
         dplyr::rename(patient_barcode = bcr_patient_barcode)
         #filter(patient_barcode %in% my_df$patient_barcode) %>%
         # filter(vital_status != "[Discrepancy]")

      # tcga_surv_tab <- clinical_pan %>% # A tibble: 10,901 x 749
      #             mutate(fu_group = factor(acronym)) %>%
      #             mutate(fu_status = as.numeric(recode(vital_status, "Alive"="0", "Dead"="1"))) %>%
      #             mutate(fu_time = as.numeric(if_else(fu_status==1, days_to_death, days_to_last_followup))) %>%
      #             filter(fu_time!="NA")
      #       #
      tcga_surv_tab <- clinical_pan


   ## Match up with PAN tcga qc annotations
   ## ::::::::::::::::::::::::
   # pdata <- as_tibble(tcga_pancan_samples$quality) %>% # A tibble: 79,286 x 28
      # dplyr::left_join(tcga_surv_tab %>% dplyr::select(c(patient_barcode, bcr_patient_uuid, gender, fu_group, fu_status, fu_time, vital_status, days_to_birth,
      #    days_to_death, days_to_last_followup, age_at_initial_pathologic_diagnosis, race)),
      #    by="patient_barcode") %>%
      # dplyr::select(sample_id, patient_barcode, cancer.type, aliquot_barcode, Do_not_use, gender, fu_group, fu_status, fu_time, vital_status, days_to_birth,
      #    days_to_death, days_to_last_followup, age_at_initial_pathologic_diagnosis, race, bcr_patient_uuid, everything()) %>%
      # mutate(cancer.type = if_else(cancer.type=="", as.character(fu_group), cancer.type)) %>%
      # dplyr::arrange(cancer.type)

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

   pdata_tcgaPanCan <- pdata
   rm(tcga_pancan_samples, pdata,  clinical_pan)






   ## ------------------
   ##  UTC PDATA STEP 2
   ## ------------------

   # NOTE!!! Curated pan_gdc object includes some 772 samples that are run on GA platform ('IlluminaGA_RNASeqV2').
         # COAD-TP ESCA-NT LAML-TB   OV-TP   OV-TR READ-TP STAD-NT STAD-TP UCEC-NT UCEC-TP
         #  171       1       6     137       2      72       2       1      11     369
      ## Create separate pdata/dataset for GA-samples


   # NOTE 2!!! All OV-samples on GA platform are also listed as Do_not_use == True. (and listed as platform=="" in quality)
      # Redaction:Genotype mismatch:Pairwise comparisons of expression results between platforms HT_HG-U133A and AgilentG4502A_07 suggests that results from this aliquot were from a different…
      # Make sure that the OV samples are note used

   # NOTE 3!!! some reåplicates in names exists, i.e. _r1 - remove these samples !!! - ALWAYS USE ALIQUOT ID TO MATCH WHEN POSSIBLE




   ## 1. Remove _r replicates
   sampleNames(pan_gdc)[grep("_r", sampleNames(pan_gdc))]
   pan_gdc <- pan_gdc[,-grep("_r", sampleNames(pan_gdc))]

   ## 2. Create pdata_pan & pdata_pan_GA
   pdata_pan <- pData(pan_gdc)
   pdata_pan <- as_tibble(pdata_pan[,c("sample_id","aliquot_barcode", "sample_type","sample_type2")])
   pdata_pan_GA <- pdata_pan

   ## 3. Filter HiSeq samples pdata
   apa <- pdata_tcgaPanCan %>% filter(platform=="IlluminaHiSeq_RNASeqV2")
   my_samples <- intersect(apa$aliquot_barcode, pdata_pan$aliquot_barcode)
   str(my_samples) # [1:10223]
   sampleNames(pan_gdc)
   pdata_pan <- pdata_pan %>%  # A tibble: 10,223 x 31
      dplyr::filter(aliquot_barcode %in% my_samples) %>%
      dplyr::left_join(pdata_tcgaPanCan %>% dplyr::select(-sample_id) %>% dplyr::filter(platform=="IlluminaHiSeq_RNASeqV2") %>% dplyr::filter(aliquot_barcode %in% my_samples), by="aliquot_barcode")
   table(duplicated(pdata_pan$aliquot_barcode))
   table(duplicated(pdata_pan$sample_id))

   ## 4. Filter GASeq samples (here OV samples are removed sinve platform =="")
   apa <- pdata_tcgaPanCan %>% filter(platform=="IlluminaGA_RNASeqV2")
   my_samples <- intersect(apa$aliquot_barcode, pdata_pan_GA$aliquot_barcode)
   str(my_samples) #  chr [1:643]
   pdata_pan_GA <- pdata_pan_GA %>%  ## A tibble: 643 x 85
      dplyr::filter(aliquot_barcode %in% my_samples) %>%
      dplyr::left_join(pdata_tcgaPanCan %>% dplyr::select(-sample_id) %>% dplyr::filter(platform=="IlluminaGA_RNASeqV2") %>% dplyr::filter(aliquot_barcode %in% my_samples), by="aliquot_barcode")
   table(duplicated(pdata_pan_GA$aliquot_barcode))
   table(duplicated(pdata_pan_GA$sample_id))

   ## 4. Create filtered eSets. pan_gdc & pan_gdc_GA
   pan_gdc_GA <- pan_gdc[,pdata_pan_GA$sample_id]
   pan_gdc <- pan_gdc[, pdata_pan$sample_id]
   str(pdata_pan[pdata_pan$Do_not_use=="True", ])




   # write.table(pdata_pan[pdata_pan$Do_not_use=="True", ], file="~/PROJECTS/RCC_ccpRCC_2019/reTaxonomy/pan_gdc_102235set_PanCan_DoNotUse_n170.txt", sep="\t", row.names=F)
   #str(tcga_pancan_samples)




   ## STEP 3 - UTC PDATA (chen etc)
   ## -----------------------

   ## get histoology and molecular grouping information from previous analyses including:
   # * Taxonomy paper (Chen 2016), Includes annotations on exclusion of samples for the three RCC TCGA marker papers
   # * RCC FGFR2 isoform switch paper ()
   # * ccRCC histology review (Mod Path Favazza 2017)
   # REMOVE OUTLIERS (Lindgren 2017)
   # ----------------
   table(is.na(pdata_pan$sample_id))
   str(sampleNames(pan_gdc))
   #identical(sampleNames(pan_gdc), pdata$sample_id)

   pdata_chen <- read.delim("~/R_PKG/dlfoo2/data-raw/RCC_taxonomy_Chen2016_TableS2.txt", as.is=T)
   chen_cols <- c("sample_id","Path.review.notes","explicitly.removed.from.marker.paper",
            "VHL_silenced","TFE3.TFEB.fusion","TFE3.fusion","SWI.SNF.module","PI3K.module","TP53","VHL","SETD2","BAP1","PBRM1","KDM5C","NF2","PTEN","ARID1A","MICALCL",
            "STAG2", "SLC1A3", "CDKN1A",  "MTOR",   "MET", "SMARCB1", "TCEB1", "NFE2L2", "PIK3CA", "MLL3",
            "OX.PHOS","Cell.Cycle", "hypoxia","YAP1", "NRF2.ARE")
   pdata <- dplyr::left_join(pdata_pan, pdata_chen[,chen_cols], by="sample_id") # A tibble: 10,223 x 95

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
   # pdata %>% dplyr::mutate_at("VHL", list(~myfoo))
   # pdata <- pdata %>% mutate_at(chen_cols2, list(~myfoo))


   pdata <- pdata %>% mutate_at(chen_cols2, myfoo)
   table(pdata$TP53)
   str(pdata)
   Hmisc::describe(pdata)

   pdata_old <- read.delim("~/R_PKG/dlfoo2data/data-raw/pdata_utc_1435set_edit_v20170830.txt", as.is=T, na.strings = c("","NA"))
   str(pdata_old)
   pdata_old_cols <- c("sample_id","EOS_davis_2014","SARC_davis_2014","FGFR2_IIIb_Zhao","histology_davis_2014","variant_histology_comment","uro_taxonomy","chen_taxonomy","dl_n_taxonomy")
   table(duplicated(pdata_old[,pdata_old_cols]$sample_id))
   table(duplicated(pdata$sample_id))
   pdata <- dplyr::left_join(pdata, as_tibble(pdata_old[,pdata_old_cols]), by="sample_id")
   str(pdata) # str(pdata) # A tibble: 10,223 x 103


   pdata <- pdata %>%
      #mutate(sample_type2 = pData(utc_nt)$sample_type2) %>%
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




   ## :::::::::::::::::::::::::::::
   ##    ADD UPDATE DL TAXONOMY 2019
   ## :::::::::::::::::::::::::::::
   # pdata <- pdata_pan_df[sampleNames(utc_nt),]
   #stopifnot(identical(pdata$sample_id, sampleNames(utc_nt)))

   apa <- pdata
   #  retax <- as_tibble(read.delim("~/PROJECTS/RCC_ccpRCC_2019/reTaxonomy/reTaxonomy_79samples_20190613.txt", as.is=T))
   # redo update 5 september - outlier CIMP to pCIMP.
   retax <- as_tibble(read.delim("~/PROJECTS/RCC_ccpRCC_2019/reTaxonomy/reTaxonomy_81samples_20190929.txt", as.is=T))
   str(retax)
   table(retax$reTax)
   sum(table(retax$reTax))

   # retax[retax$sample_id=="KIRC_TP_A8OX",]
   # data.frame(pdata_pre[pdata_pre$sample_id=="KIRC_TP_A8OX",])
   # data.frame(pdata_tcgaPanCan[pdata_tcgaPanCan$sample_id=="KIRC_TP_A8OX",])
   #

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

   table(pdata$tax_simp)
   pdata[which(pdata$tax_simp=="outlier"),]
   #data.frame(pdata[pdata$sample_id=="KIRC_TP_A8OX",])


   sum(table(pdata$taxonomy_published))
   sum(table(pdata$tax_dl))
   sum(table(pdata$tax_simp))


   # ##  pdata ConsensuClsutering RCC/BLCA k25
   # ## ------------------------------
   # my.k <- 25
   # my.cc <- dlfoo2::ConsensusCluster_process(cc.dir = "~/PROJECTS/RCC_ccpRCC_2019/ConsensusClust/ConsensusCluster_utc_nt_5000_1437set_5000genes_reps.1000_maxK.25_innerLinkage.ward_sd.rank.5000/", k = my.k)
   # #stopifnot(identical(my.cc$sample.names, pdata$sample_id))
   # my_tab <- data.frame(sample_id=my.cc$sample.names, cc_1437set_k25=paste0("cc_", my.cc$cc), stringsAsFactors = F)
   #    # apa <- pdata
   # pdata <- pdata %>%
   #    left_join(my_tab, by="sample_id") %>%
   #    mutate(cc_1437set_k25 = factor(cc_1437set_k25, levels=paste0("cc_",1:25)))
   # sum(table(pdata$taxonomy_published, pdata$cc_1437set_k25))
   # sum(table(pdata$cc_1437set_k25))
   #
   #
   # my.k <- 5
   # my.cc <- dlfoo2::ConsensusCluster_process(cc.dir = "~/PROJECTS/RCC_ccpRCC_2019/ConsensusClust/ConsensusCluster_utc_nt_5000_1437set_5000genes_reps.1000_maxK.25_innerLinkage.ward_sd.rank.5000/", k = my.k)
   # #stopifnot(identical(my.cc$sample.names, pdata$sample_id))
   # my_tab <- data.frame(sample_id=my.cc$sample.names, cc_1437set_k5=paste0("cc_", my.cc$cc), stringsAsFactors = F)
   #    # apa <- pdata
   # pdata <- pdata %>%
   #    left_join(my_tab, by="sample_id") %>%
   #    mutate(cc_1437set_k5 = factor(cc_1437set_k5, levels=paste0("cc_",1:5)))
   # sum(table(pdata$taxonomy_published, pdata$cc_1437set_k5))
   # sum(table(pdata$cc_1437set_k5))
   # table(pdata$cc_1437set_k5, pdata$taxonomy_published)


   ##  pdata ConsensuClsutering RCC k25
   ## ------------------------------
   my.k <- 25
   my.cc <- dlfoo2::ConsensusCluster_process(cc.dir = "~/PROJECTS/RCC_ccpRCC_2019/ConsensusClust/ConsensusCluster_rcc_t_5000_880set_5000genes_reps.1000_maxK.25_innerLinkage.ward_sd.rank.5000/", k = my.k)
   #stopifnot(identical(my.cc$sample.names, pdata$sample_id))
   my_tab <- data.frame(sample_id=my.cc$sample.names, cc_rcc_t_k25=paste0("cc_", my.cc$cc), stringsAsFactors = F)
   apa <- pdata

   pdata <- pdata %>%
      left_join(my_tab, by="sample_id") %>%
      mutate(cc_rcc_t_k25 = factor(cc_rcc_t_k25, levels=paste0("cc_",1:25)))
   sum(table(pdata$taxonomy_published, pdata$cc_rcc_t_k25))
   sum(table(pdata$cc_rcc_t_k25))


   #pdata <- dlfoo2data::pdata_panCanGex

   cc_k25_grouping <- list(
      pRCC_a = paste0("cc_",c(23,25,16)),
      pRCC_b = paste0("cc_",c(20,24)),
      pRCC_c = paste0("cc_",c(9)),
      pRCC_cimp = paste0("cc_",c(22)),
      ccRCC_a = paste0("cc_",c(18,14,6,19,12)),
      ccRCC_b = paste0("cc_",c(10,17,8,21,15)),
      ccRCC_c = paste0("cc_",c(7,13,11,5)),
      chRCC = paste0("cc_",c(1))
      )
   cc_k25_grouping <- unlist(cc_k25_grouping)
   x <- names(cc_k25_grouping)
   names(x) <- cc_k25_grouping

   apa2 <- pdata %>%
      mutate(tax_simp2 = as.character(tax_simp)) %>%
      mutate(tax_simp2 = if_else(as.character(cc_rcc_t_k25) %in% cc_k25_grouping, as.character(cc_rcc_t_k25), tax_simp2)) %>%
      mutate(tax_simp2 = if_else(tax_simp2 %in% cc_k25_grouping, x[tax_simp2], tax_simp2)) %>%
      mutate(tax_simp2 = if_else(tax_simp %in% c("dCIMP","ccpRCC","chONC","pONC","mesRCC","metanephric","mtscRCC","outlier"), as.character(tax_simp), tax_simp2)) %>%
      mutate(tax_simp2 = factor(tax_simp2, levels=c("kidney_n",paste0("ccRCC_a",1:5), paste0("ccRCC_b",1:5), paste0("ccRCC_c",1:4),
            paste0("pRCC_a",1:3), paste0("pRCC_b",1:2),
            "pRCC_c","pRCC_cimp", "pCIMP","dCIMP","ccpRCC","chRCC","chONC","pONC","mesRCC","metanephric","mtscRCC","outlier",
            "bladder_n","Urobasal","GU","Basal/SCCL","SC/NE","Mes-Inf","Inf_other")))
   table(apa2$tax_simp2)
   pdata <- apa2




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
   pdata # A tibble: 10,223 x 105
   stopifnot(identical(pdata$sample_id, pdata_pan$sample_id))
   pdata_pan <- pdata
   table(pdata_pan$iCluster, pdata_pan$sample_type2)

   # also add for GenomeAnalyzer pdata_pan_GA
   u <- match(pdata_pan_GA$portion_barcode, annotations_pan$portion_barcode)
   table(duplicated(u[!is.na(u)]))
   pdata_GA <- pdata_pan_GA %>% mutate(leuk_estimate = annotations_pan$leuk_estimate[u]) # # A tibble: 643 x 86
   stopifnot(identical(pdata_GA$sample_id, pdata_pan_GA$sample_id))
   pdata_pan_GA <- pdata_GA # A tibble: 643 x 86




   ## ::::::::::::::::::::::
   ##    UPDATE/KEEP PAN GDC GEX OBJECT
   ## ::::::::::::::::::::::
   stopifnot(identical(sampleNames(pan_gdc), pdata_pan$sample_id))
   pdata_pan_df <- as.data.frame(pdata_pan)
   rownames(pdata_pan_df) <- pdata_pan_df$sample_id
   pData(pan_gdc) <- pdata_pan_df
   pan_gdc_counts.temp <- pan_gdc_counts
   pan_gdc_counts <- pan_gdc_counts[,pdata_pan_df$sample_id]
   pData(pan_gdc_counts) <- pdata_pan_df


   ## Update GA object
   stopifnot(identical(sampleNames(pan_gdc_GA), pdata_pan_GA$sample_id))
   pdata_pan_df <- as.data.frame(pdata_pan_GA)
   rownames(pdata_pan_df) <- pdata_pan_df$sample_id
   pData(pan_gdc_GA) <- pdata_pan_df
   pan_gdc_counts_GA <- pan_gdc_counts.temp[,pdata_pan_df$sample_id]
   pData(pan_gdc_counts_GA) <- pdata_pan_df
   # pData(pan_gdc_counts) <- pdata_pan_df

   rm(apa, my_samples)


   ## ::::::::::::::::::::::
   ##    UPDATE UTC GEX object
   ## ::::::::::::::::::::::

   my_index <- match(sampleNames(utc_nt), pdata_pan$sample_id)
   pdata_utc_nt <- pdata_pan[pdata_pan$sample_id %in% sampleNames(utc_nt) ,]
   stopifnot(identical(pdata_utc_nt$sample_id, sampleNames(utc_nt)))
   # rownames(pdata_utc_nt) <- pdata_utc_nt$sample_id

   pdata_utc_nt_df <- as.data.frame(pdata_utc_nt)
   rownames(pdata_utc_nt_df) <- pdata_utc_nt_df$sample_id
   pData(utc_nt) <- pdata_utc_nt_df
   pData(utc_nt_counts) <- pdata_utc_nt_df
   pData(rcc_nt) <- pdata_utc_nt_df[sampleNames(rcc_nt), ]
   pData(rcc_nt_counts) <- pdata_utc_nt_df[sampleNames(rcc_nt_counts), ]
   pData(rcc_t) <- pdata_utc_nt_df[sampleNames(rcc_t), ]
   pData(rcc_t_counts) <- pdata_utc_nt_df[sampleNames(rcc_t_counts), ]
   rm(pdata, pdata_pan_df)





   ## ::::::::::::::::::::::
   ##    SAVE ALL OBJECTS
   ## ::::::::::::::::::::::

   ## pdata_tcgaPanCan
   # pdata_panCanFull <- pdata_tcgaPanCan
   # usethis::use_data(pdata_panCanFull, overwrite = T)
   #
   # ## pdata_pan ::  named pdata_panGex (pan pdata but only Hiseq GEX matched wtih pan_gdc)
   # pdata_panCanGex <- pdata_pan
   # usethis::use_data(pdata_panCanGex, overwrite = T)
   #
   #
   # usethis::use_data(pdata_utc_nt, overwrite = T)
   # usethis::use_data(tcga_surv_tab, overwrite = T)
   #
   # usethis::use_data(pan_gdc, overwrite = T)
   # usethis::use_data(pan_gdc_GA, overwrite = T)
   #
   #
   # usethis::use_data(utc_nt, overwrite = T)
   # usethis::use_data(rcc_nt, overwrite = T)
   # usethis::use_data(rcc_t, overwrite = T)
   #
   #
   # usethis::use_data(pan_gdc_counts, overwrite = T)
   # usethis::use_data(pan_gdc_counts_GA, overwrite = T)
   #
   # usethis::use_data(utc_nt_counts, overwrite = T)
   # usethis::use_data(rcc_nt_counts, overwrite = T)
   # usethis::use_data(rcc_t_counts, overwrite = T)
   #
   #
   #







