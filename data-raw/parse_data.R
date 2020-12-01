setwd("~/R_PKG/dlfoo2/")

usethis::proj_set(rstudioapi::getActiveProject())




## -----------------------------------
##    Sample Color tables and annotations etc
## ------------------------------------

#usethis::proj_set(rstudioapi::getActiveProject())

   color_tab_subtypes <- read.delim(file="/Volumes/Crab2000/RESOURCES/SAMPLE_ANNOTATIONS/COLOR_KEYS/color_tab_subtypes_v1.1.txt")
      color_subtypes <- as.character(color_tab_subtypes$color)
      names(color_subtypes) <- as.character(color_tab_subtypes$sample)
      color_subtypes["BLCA"]
      # Fix na-blank
      #color_subtypes[names(color_subtypes)==""]

      shape_subtypes <- as.character(color_tab_subtypes$shape)
      names(shape_subtypes) <- as.character(color_tab_subtypes$sample)

   annot_tab <- read.delim("/Volumes/Crab2000/RESOURCES/SAMPLE_ANNOTATIONS/COLOR_KEYS/annotation_names.txt")
      color_annotations <- as.character(annot_tab$color)
      names(color_annotations) <- as.character(annot_tab$sample)

   usethis::use_data(color_subtypes, overwrite = T)
   usethis::use_data(shape_subtypes, overwrite = T)
   usethis::use_data(color_annotations, overwrite = T)


   color_celltypes <- read.delim(file="/Volumes/Crab2000/RESOURCES/SAMPLE_ANNOTATIONS/COLOR_KEYS/color_tab_cellTypes_v1.0.txt")
   usethis::use_data(color_celltypes, overwrite = T)




# es_ccle = readRDS(file="/Volumes/MacPro2TB/RESOURCES/CCLE/Curated/ccle_1156set_rpkm_DepMap_18q3_20180718_log2.rds")
# usethis::use_data(es_ccle, overwrite = T)
#
#
#
# clinical_pan <- readRDS(file = "/Users/david/RESOURCES/MSig_DL_Samples/pantcga_clinical.rds")
# usethis::use_data(clinical_pan, overwrite = T)



## Centromere data
load("~/RESOURCES/Annotation Files/centromere.limits.Rdata")
usethis::use_data(centromere.limit.hg17, overwrite = T)
usethis::use_data(centromere.limit.hg18, overwrite = T)
usethis::use_data(centromere.limit.hg19, overwrite = T)

















   # # Prapare data. Curated data in  "~/R/R_PUBLIC_DATA/R_scripts/R_Curate_Data R_objects_v3"
   # # ............................................
   #    require(Biobase)
   #    require(dplyr)
   #    options(dplyr.width = Inf)
   #    require(tibble)
   #    require(tidyverse)
   #    require(tidyr)
   #
   #    setwd("~/R_PKG/dlfoo2")
   #    df_hypox <- as.tibble(read.delim("data-raw/plot.table.hypoxia.public.txt", as.is=T))()
   #    df_renal <- as.tibble(read.delim("data-raw/plot.table.renal.public.txt", as.is=T))
   #
   # load_obj <- function(f){
   #     env <- new.env()
   #     nm <- load(f, env)[1]
   #     env[[nm]]}
   #
   # eset_list <- lapply(df_hypox$eset, function(x){
   #       load_obj(paste0("data-raw/",x,".rdata"))})
   # names(eset_list) <- df_hypox$eset
   #
   # lapply(names(eset_list), function(x){
   #    xx <- eset_list[[x]]
   #    cat("”n ...................")
   #    cat("\n ... ", x)
   #
   #    cat("\n ... ...  ", annotation(xx))
   #    cat("\n ... ...  ", head(featureNames(xx)))
   #    print(head(fData(xx)))
   #    print(head(pData(xx)[,1:7]))
   #    print(Hmisc::describe(exprs(xx)[,1]))
   #    })

