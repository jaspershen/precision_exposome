##avoid source
no_function()

sxtTools::setwd_project()
rm(list = ls())

setClass(
  Class = "databaseClass",
  representation(
    database.info = "list",
    spectra.info = "data.frame",
    spectra.data = "list"
  ),
  prototype = list(
    database.info = list(),
    spectra.info = data.frame(matrix(nrow = 0, ncol = 0), stringsAsFactors = FALSE),
    spectra.data = list()
  )
)


setwd("data_20200511/database/")

##------------------------------------------------------------------------------
#####list 1
library(tidyverse)
##list1
list1 <- readxl::read_xlsx("List 1.xlsx")

list1 <- list1[-1,]

list1 <- 
  list1 %>% 
  dplyr::rename("Lab.ID" = "Chemical Name",
                "Formula" = "Molecular Formula",
                "mz" = "Mono-isotopic Mass") %>% 
  dplyr::mutate(Compound.name = Lab.ID,
                RT = NA,
                CAS.ID = NA, 
                HMDB.ID = NA, 
                KEGG.ID = NA, 
                # Formula = NA, 
                mz.pos = NA,
                mz.neg = NA, 
                Submitter = NA) %>% 
  dplyr::select(Lab.ID, Compound.name, mz, RT, CAS.ID, 
                HMDB.ID, KEGG.ID, 
                Formula, mz.pos, 
                mz.neg, Submitter, everything())

database.info <- list(
  "Version" = "0.0.1",
  "Source" = "Xiaotao Shen",
  "Link" = "shenxt.me",
  "Creater" = "Xiaotao Shen",
  "Email" = "shenxt1990@163.com",
  "RT" = FALSE
)

spectra.info <- list1

Spectra <- list()

list1_ms1_database <- 
  new(
    Class = "databaseClass",
    database.info = database.info,
    spectra.info = spectra.info,
    spectra.data = Spectra
  )


save(list1_ms1_database, file = "list1_ms1_database")

sxtTools::setwd_project()
setwd("data_20200511/database/")
library(tidyverse)


##------------------------------------------------------------------------------
##list2
list2 <- readxl::read_xlsx("List 2.xlsx")


list2 <- 
  list2 %>% 
  dplyr::rename("Compound.name" = "Chemical Name",
                "Formula" = "MolecularFormula",
                "mz" = "Monoisotopic_Mass") %>% 
  dplyr::mutate(Lab.ID = Compound.name,
                RT = NA,
                CAS.ID = NA, 
                HMDB.ID = NA, 
                KEGG.ID = NA, 
                # Formula = NA, 
                mz.pos = NA,
                mz.neg = NA, 
                Submitter = NA) %>% 
  dplyr::select(Lab.ID, Compound.name, mz, RT, CAS.ID, 
                HMDB.ID, KEGG.ID, 
                Formula, mz.pos, 
                mz.neg, Submitter, everything())


database.info <- list(
  "Version" = "0.0.1",
  "Source" = "Xiaotao Shen",
  "Link" = "shenxt.me",
  "Creater" = "Xiaotao Shen",
  "Email" = "shenxt1990@163.com",
  "RT" = FALSE
)

spectra.info <- list2

Spectra <- list()

list2_ms1_database <- 
  new(
    Class = "databaseClass",
    database.info = database.info,
    spectra.info = spectra.info,
    spectra.data = Spectra
  )
  
save(list2_ms1_database, file = "list2_ms1_database")


##------------------------------------------------------------------------------
##nonspecific_biomarkers
nonspecific_biomarkers <- readxl::read_xlsx("Nonspecific biomarkers.xlsx")

nonspecific_biomarkers <- 
  nonspecific_biomarkers[c(1:7),]

nonspecific_biomarkers <- 
  nonspecific_biomarkers %>% 
  dplyr::rename("Compound.name" = "Chemical",
                "mz" = "Mass") %>% 
  dplyr::mutate(Lab.ID = Compound.name,
                RT = NA,
                CAS.ID = NA, 
                HMDB.ID = NA, 
                KEGG.ID = NA, 
                Formula = NA,
                mz.pos = NA,
                mz.neg = NA, 
                Submitter = NA) %>% 
  dplyr::select(Lab.ID, Compound.name, mz, RT, CAS.ID, 
                HMDB.ID, KEGG.ID, 
                Formula, mz.pos, 
                mz.neg, Submitter, everything())

database.info <- list(
  "Version" = "0.0.1",
  "Source" = "Xiaotao Shen",
  "Link" = "shenxt.me",
  "Creater" = "Xiaotao Shen",
  "Email" = "shenxt1990@163.com",
  "RT" = FALSE
)

spectra.info <- nonspecific_biomarkers

Spectra <- list()

nonspecific_biomarkers_ms1_database <- 
  new(
    Class = "databaseClass",
    database.info = database.info,
    spectra.info = spectra.info,
    spectra.data = Spectra
  )



save(nonspecific_biomarkers_ms1_database, file = "nonspecific_biomarkers_ms1_database")




##------------------------------------------------------------------------------
###############select_exposome
select_exposome <- readxl::read_xlsx("Select exposome.xlsx")

# cas_id <- lapply(select_exposome$MetabID, function(x){
#   metflow2::transID(query = x, from = "Chemical name", to = "CAS")
# })
# 
# lapply(cas_id, function(x){
#   nrow(x)
# }) %>% unlist()
# 
# cas_id <- 
# cas_id %>% 
#   do.call(rbind, .)
# 
# idx <- which(is.na(cas_id$CAS))
# cas_id$`Chemical name`[idx]


select_exposome <- 
  select_exposome %>% 
  dplyr::rename("Compound.name" = "MetabID",
                "mz" = "mz") %>% 
  dplyr::mutate(Lab.ID = Compound.name,
                RT = NA,
                CAS.ID = NA, 
                HMDB.ID = NA, 
                KEGG.ID = NA, 
                Formula = NA,
                mz.pos = NA,
                mz.neg = NA, 
                Submitter = NA) %>% 
  dplyr::select(Lab.ID, Compound.name, mz, RT, CAS.ID, 
                HMDB.ID, KEGG.ID, 
                Formula, mz.pos, 
                mz.neg, Submitter, everything())

save(select_exposome, file = "select_exposome")

database.info <- list(
  "Version" = "0.0.1",
  "Source" = "Xiaotao Shen",
  "Link" = "shenxt.me",
  "Creater" = "Xiaotao Shen",
  "Email" = "shenxt1990@163.com",
  "RT" = FALSE
)

spectra.info <- select_exposome

Spectra <- list()

select_exposome_ms1_database <- 
  new(
    Class = "databaseClass",
    database.info = database.info,
    spectra.info = spectra.info,
    spectra.data = Spectra
  )



save(select_exposome_ms1_database, file = "select_exposome_ms1_database")





##------------------------------------------------------------------------------
###############specific_biomarker
specific_biomarker <- readr::read_csv("Specific biomarkers.csv")

specific_biomarker <- 
  specific_biomarker %>% 
  dplyr::rename("Compound.name" = "Name",
                "mz" = "Mono. mass",
                "Lab.ID" = "ID") %>% 
  dplyr::mutate(
                RT = NA,
                CAS.ID = NA, 
                HMDB.ID = NA, 
                KEGG.ID = NA, 
                Formula = NA,
                mz.pos = NA,
                mz.neg = NA, 
                Submitter = NA) %>% 
  dplyr::select(Lab.ID, Compound.name, mz, RT, CAS.ID, 
                HMDB.ID, KEGG.ID, 
                Formula, mz.pos, 
                mz.neg, Submitter, everything())

database.info <- list(
  "Version" = "0.0.1",
  "Source" = "Xiaotao Shen",
  "Link" = "shenxt.me",
  "Creater" = "Xiaotao Shen",
  "Email" = "shenxt1990@163.com",
  "RT" = FALSE
)

spectra.info <- specific_biomarker

Spectra <- list()

specific_biomarker_ms1_database <- 
  new(
    Class = "databaseClass",
    database.info = database.info,
    spectra.info = spectra.info,
    spectra.data = Spectra
  )



save(specific_biomarker_ms1_database, file = "specific_biomarker_ms1_database")







##------------------------------------------------------------------------------
###############t3db
t3db <- readr::read_csv("T3DB.csv")

t3db <- 
  t3db %>% 
  dplyr::rename("Compound.name" = "Name",
                "mz" = "Monoisotopic Mass",
                "Lab.ID" = "T3DB ID") %>% 
  dplyr::mutate(
    RT = NA,
    CAS.ID = NA, 
    HMDB.ID = NA, 
    KEGG.ID = NA, 
    Formula = NA,
    mz.pos = NA,
    mz.neg = NA, 
    Submitter = NA) %>% 
  dplyr::select(Lab.ID, Compound.name, mz, RT, CAS.ID, 
                HMDB.ID, KEGG.ID, 
                Formula, mz.pos, 
                mz.neg, Submitter, everything())

t3db <- 
  t3db %>% 
  dplyr::filter(!stringr::str_detect(Status, "Inorganic compounds"))%>% 
  dplyr::filter(!stringr::str_detect(Toxicity, "Gas"))

database.info <- list(
  "Version" = "0.0.1",
  "Source" = "Xiaotao Shen",
  "Link" = "shenxt.me",
  "Creater" = "Xiaotao Shen",
  "Email" = "shenxt1990@163.com",
  "RT" = FALSE
)

spectra.info <- t3db

Spectra <- list()

t3db_ms1_database <- 
  new(
    Class = "databaseClass",
    database.info = database.info,
    spectra.info = spectra.info,
    spectra.data = Spectra
  )



save(t3db_ms1_database, file = "t3db_ms1_database")



##------------------------------------------------------------------------------
##blood exposome database
blood_exposome <- readr::read_csv("BloodExpsomeDatabase_version_1.0.csv")

blood_exposome <- 
  blood_exposome %>% 
  dplyr::rename("Compound.name" = "Compound Name",
                "mz" = "ExactMass",
                # "Lab.ID" = "T3DB ID",
                "CAS.ID" = "PubChem CID",
                "HMDB.ID" = "HMDB ID", 
                "KEGG.ID" = "KEGG ID", 
                "Formula" = "Molecular Formula",
                ) %>% 
  dplyr::mutate(
    Lab.ID = paste("BE", 1:nrow(blood_exposome), sep = ""),
    RT = NA,
    mz.pos = NA,
    mz.neg = NA, 
    Submitter = NA) %>% 
  dplyr::select(Lab.ID, Compound.name, mz, RT, CAS.ID, 
                HMDB.ID, KEGG.ID, 
                Formula, mz.pos, 
                mz.neg, Submitter, everything())

# blood_exposome <- 
#   blood_exposome %>% 
#   dplyr::filter(!stringr::str_detect(Status, "Inorganic compounds"))%>% 
#   dplyr::filter(!stringr::str_detect(Toxicity, "Gas"))

database.info <- list(
  "Version" = "0.0.1",
  "Source" = "Xiaotao Shen",
  "Link" = "shenxt.me",
  "Creater" = "Xiaotao Shen",
  "Email" = "shenxt1990@163.com",
  "RT" = FALSE
)

spectra.info <- blood_exposome

Spectra <- list()

blood_exposome_ms1_database <- 
  new(
    Class = "databaseClass",
    database.info = database.info,
    spectra.info = spectra.info,
    spectra.data = Spectra
  )



save(blood_exposome_ms1_database, file = "blood_exposome_ms1_database")





















