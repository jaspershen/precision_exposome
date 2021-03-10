sxtTools::setwd_project()
library(tidyverse)
setwd("data_20200511/metabolome/")
rm(list = ls())

load("clinic_data")
load("met_data")
load("met_tag")

met_data <-
  met_data %>%
  dplyr::filter(SubjectID == "69-001") %>%
  dplyr::filter(CollectionDate >= "2016-01-12",
                CollectionDate <= "2016-03-03")

variable_info <- met_tag

expression_data <- met_data

sample_info <- 
  expression_data %>% 
  dplyr::select(sample_id = SampleID, 
                subject_id = SubjectID, 
                CollectionDate,
                CL1, CL2, CL3, CL4)

expression_data <- 
  expression_data %>% 
  dplyr::select(-c(SampleID, SubjectID, 
                CollectionDate,
                CL1, CL2, CL3, CL4))

expression_data <-
  t(expression_data) %>% 
  as.data.frame()
  

colnames(expression_data) <- sample_info$sample_id

rownames(expression_data)

variable_info <- 
  variable_info %>% 
  dplyr::distinct(Compounds_ID, .keep_all = TRUE)

variable_info <- 
  variable_info %>% 
  dplyr::filter(Compounds_ID %in% rownames(expression_data))

variable_info <- 
  variable_info[match(rownames(expression_data), variable_info$Compounds_ID),]

variable_info$Compounds_ID == rownames(expression_data)

variable_info <- 
  variable_info %>% 
  dplyr::select(Compounds_ID, everything()) %>% 
  dplyr::rename(peak_name = Compounds_ID)

save(variable_info, file = "variable_info")
save(sample_info, file = "sample_info")
save(expression_data, file = "expression_data")


library(openxlsx)
wb = createWorkbook()
modifyBaseFont(wb, fontSize = 12, fontName = "Arial Narrow")
addWorksheet(wb, sheetName = "Sample information", gridLines = TRUE)
addWorksheet(wb, sheetName = "Variable information", gridLines = TRUE)
addWorksheet(wb, sheetName = "Expression data", gridLines = TRUE)
freezePane(wb, sheet = 1, firstRow = TRUE, firstCol = TRUE) 
freezePane(wb, sheet = 2, firstRow = TRUE, firstCol = TRUE) 
freezePane(wb, sheet = 3, firstRow = TRUE, firstCol = FALSE) 
writeDataTable(wb, sheet = 1, x = sample_info,
               colNames = TRUE, rowNames = FALSE)
writeDataTable(wb, sheet = 2, x = variable_info,
               colNames = TRUE, rowNames = FALSE)
writeDataTable(wb, sheet = 3, x = expression_data,
               colNames = TRUE, rowNames = FALSE)

saveWorkbook(wb, "metabolome_data.xlsx", overwrite = TRUE)




# clinic_data <-
#   clinic_data %>%
#   dplyr::filter(SubjectID == "69-001") %>%
#   dplyr::rename(subject_id = SubjectID,
#                 sample_id = SampleID)
# save(clinic_data, file = "clinic_data")
load("clinic_data")

###match using database from Peng
# variable_info_pos <- 
#   variable_info %>% 
#   dplyr::filter(stringr::str_detect(peak_name, "p")) %>% 
#   dplyr::select(peak_name, Mass) %>% 
#   dplyr::mutate(rt = stringr::str_split(peak_name, "\\_") %>% 
#                   lapply(function(x)x[3]) %>% 
#                   unlist()
#                   ) %>% 
#   dplyr::mutate(rt = as.numeric(rt) * 60) %>% 
#   dplyr::rename(name = peak_name, mz = Mass)
# 
# variable_info_pos$mz <- 
#   stringr::str_split(variable_info_pos$mz, pattern = "\\_") %>% 
#   lapply(function(x){x[1]}) %>% unlist() %>% as.numeric()
# 
# variable_info_neg <- 
#   variable_info %>% 
#   dplyr::filter(stringr::str_detect(peak_name, "n")) %>% 
#   dplyr::select(peak_name, Mass) %>% 
#   dplyr::mutate(rt = 
#                   stringr::str_extract(peak_name, 
#                                        "[0-9]{1}\\.[0-9]{1}")) %>% 
#   dplyr::mutate(rt = as.numeric(rt) * 10) %>% 
#   dplyr::rename(name = peak_name, mz = Mass)
# 
# variable_info_neg$mz <- 
#   stringr::str_split(variable_info_neg$mz, pattern = "\\_") %>% 
#   lapply(function(x){x[1]}) %>% unlist() %>% as.numeric()
# 
# variable_info_pos_hilic <- 
#   variable_info_pos %>% 
#   dplyr::filter(stringr::str_detect(name, "HILIC"))
# 
# variable_info_pos_rplc <- 
#   variable_info_pos %>% 
#   dplyr::filter(stringr::str_detect(name, "RPLC"))
# 
# variable_info_neg_hilic <- 
#   variable_info_neg %>% 
#   dplyr::filter(stringr::str_detect(name, "HILIC"))
# 
# variable_info_neg_rplc <- 
#   variable_info_neg %>% 
#   dplyr::filter(stringr::str_detect(name, "RPLC"))
# 
# write.csv(variable_info_pos, "variable_info_pos.csv", row.names = FALSE)
# write.csv(variable_info_pos_hilic, "variable_info_pos_hilic.csv", row.names = FALSE)
# write.csv(variable_info_pos_rplc, "variable_info_pos_rplc.csv", row.names = FALSE)
# 
# write.csv(variable_info_neg, "variable_info_neg.csv", row.names = FALSE)
# write.csv(variable_info_neg_hilic, "variable_info_neg_hilic.csv", row.names = FALSE)
# write.csv(variable_info_neg_rplc, "variable_info_neg_rplc.csv", row.names = FALSE)
# 
# 
# library(metID)
# 
# result1_pos <- identify_metabolites(ms1.data = "variable_info_pos.csv", 
#                                 polarity = "positive", ce = 'all',
#                                 column = "rp",
#                                 total.score.tol = 0.5,
#                                 candidate.num = 3,
#                                 threads = 3, 
#                                 database = "list1_ms1_database")
# 
# 
# result2_pos <- identify_metabolites(ms1.data = "variable_info_pos.csv", 
#                                 polarity = "positive", ce = 'all',
#                                 column = "rp",
#                                 total.score.tol = 0.5,
#                                 candidate.num = 3,
#                                 threads = 3, 
#                                 database = "list2_ms1_database")
# 
# result3_pos <- identify_metabolites(ms1.data = "variable_info_pos.csv", 
#                                 polarity = "positive",
#                                 ce = 'all',
#                                 column = "rp",
#                                 total.score.tol = 0.5,
#                                 candidate.num = 3,
#                                 threads = 3, 
#                                 database = "nonspecific_biomarkers_ms1_database")
# 
# result4_pos <- identify_metabolites(ms1.data = "variable_info_pos.csv",
#                                     polarity = "positive",
#                                     ce = 'all',
#                                     column = "rp",
#                                     total.score.tol = 0.5,
#                                     candidate.num = 3,
#                                     threads = 3,
#                                     database = "select_exposome_ms1_database")
# # load("select_exposome")
# # result4_pos <- mz_match(ms1.table = variable_info_pos, 
# #                         database = select_exposome, 
# #                         mz.error.tol = 25)
# 
# result5_pos <- identify_metabolites(ms1.data = "variable_info_pos.csv", 
#                                     polarity = "positive",
#                                     ce = 'all',
#                                     column = "rp",
#                                     total.score.tol = 0.5,
#                                     candidate.num = 3,
#                                     threads = 3, 
#                                     database = "t3db_ms1_database")
# 
# result6_pos <- identify_metabolites(ms1.data = "variable_info_pos.csv", 
#                                     polarity = "positive",
#                                     ce = 'all',
#                                     column = "rp",
#                                     total.score.tol = 0.5,
#                                     candidate.num = 3,
#                                     threads = 3, 
#                                     database = "specific_biomarker_ms1_database")
# 
# result7_pos <- identify_metabolites(ms1.data = "variable_info_pos.csv", 
#                                     polarity = "positive",
#                                     ce = 'all',
#                                     column = "rp",
#                                     total.score.tol = 0.5,
#                                     candidate.num = 3,
#                                     threads = 3, 
#                                     database = "hmdbMS1Database0.0.1")
# 
# 
# result8_pos <- identify_metabolites(ms1.data = "variable_info_pos_hilic.csv", 
#                                     polarity = "positive",
#                                     ce = 'all',
#                                     rt.match.tol = 30,
#                                     column = "hilic",
#                                     total.score.tol = 0.5,
#                                     candidate.num = 3,
#                                     threads = 3, 
#                                     database = "msDatabase_hilic0.0.2")
# 
# result9_pos <- identify_metabolites(ms1.data = "variable_info_pos_rplc.csv", 
#                                     polarity = "positive",
#                                     ce = 'all',
#                                     rt.match.tol = 30,
#                                     column = "rp",
#                                     total.score.tol = 0.5,
#                                     candidate.num = 3,
#                                     threads = 3, 
#                                     database = "msDatabase_rplc0.0.2")
# 
# 
# annotation_table_pos1 <- 
#   get_identification_table(result1_pos,
#                            result2_pos,
#                            result4_pos,
#                            result5_pos,
#                            result6_pos,
#                            result7_pos,
#                            type = "old", 
#                            candidate.num = 1)
# 
# 
# annotation_table_pos1$Identification <-
#   annotation_table_pos1$Identification %>% 
#     lapply(function(x){
#       if(is.na(x)){
#         return(NA)
#       }else{
#        x <- stringr::str_split(x, "\\{\\}")[[1]] 
#        x <- grep("\\(M\\+H\\)|\\(2M|H2O", x, value = TRUE)
#        if(length(x) == 0){
#          return(NA)
#        }else{
#          paste(x, collapse = "{}")
#        }
#       }
#     }) %>% 
#     unlist()
# 
# annotation_table_pos1 <- 
# metID::trans2newStyle(identification.table = annotation_table_pos1)
# 
# annotation_table_pos1 <- 
#   annotation_table_pos1 %>% 
#   dplyr::filter(!is.na(Compound.name))
# 
# annotation_table_pos2 <- 
#   get_identification_table(result8_pos,
#                            type = "old", 
#                            candidate.num = 1)
# 
# annotation_table_pos2$Identification <-
#   annotation_table_pos2$Identification %>% 
#   lapply(function(x){
#     if(is.na(x)){
#       return(NA)
#     }else{
#       x <- stringr::str_split(x, "\\{\\}")[[1]] 
#       x <- grep("\\(M\\+H\\)|\\(2M|H2O", x, value = TRUE)
#       if(length(x) == 0){
#         return(NA)
#       }else{
#         paste(x, collapse = "{}")
#       }
#     }
#   }) %>% 
#   unlist()
# 
# annotation_table_pos2 <- 
#   metID::trans2newStyle(identification.table = annotation_table_pos2)
# 
# annotation_table_pos2 <- 
#   annotation_table_pos2 %>% 
#   dplyr::filter(!is.na(Compound.name))
# 
# 
# annotation_table_pos3 <- 
#   get_identification_table(result9_pos,
#                            type = "old", 
#                            candidate.num = 1)
# 
# 
# annotation_table_pos3$Identification <-
#   annotation_table_pos3$Identification %>% 
#   lapply(function(x){
#     if(is.na(x)){
#       return(NA)
#     }else{
#       x <- stringr::str_split(x, "\\{\\}")[[1]] 
#       x <- grep("\\(M\\+H\\)|\\(2M|H2O", x, value = TRUE)
#       if(length(x) == 0){
#         return(NA)
#       }else{
#         paste(x, collapse = "{}")
#       }
#     }
#   }) %>% 
#   unlist()
# 
# annotation_table_pos3 <- 
#   metID::trans2newStyle(identification.table = annotation_table_pos3)
# 
# 
# annotation_table_pos3 <- 
#   annotation_table_pos3 %>% 
#   dplyr::filter(!is.na(Compound.name))
# 
# 
# # result4_pos <- 
# #   result4_pos %>% 
# #   dplyr::filter(!is.na(Compound.name))
# # 
# # annotation_table_pos1 <- 
# #   annotation_table_pos1 %>% 
# #   dplyr::filter(!name %in% result4_pos$name) 
# # 
# # 
# # annotation_table_pos1 <- 
# #   rbind(annotation_table_pos1, result4_pos)
# 
# annotation_table_pos1 <- 
#   annotation_table_pos1 %>% 
#   dplyr::filter(!name %in% annotation_table_pos2$name) %>% 
#   dplyr::filter(!name %in% annotation_table_pos3$name)
# 
# 
# 
# annotation_table_pos <- 
#   rbind(annotation_table_pos1[,-4],
#         annotation_table_pos2, 
#         annotation_table_pos3)
# 
# write.csv(annotation_table_pos,
#           file = "annotation_table_pos.csv", 
#           row.names = FALSE)
# 
# 
# 
# #####negative
# result1_neg <- identify_metabolites(ms1.data = "variable_info_neg.csv", 
#                                     polarity = "negative", ce = 'all',
#                                     column = "rp",
#                                     total.score.tol = 0.5,
#                                     candidate.num = 3,
#                                     threads = 3, 
#                                     database = "list1_ms1_database")
# 
# 
# result2_neg <- identify_metabolites(ms1.data = "variable_info_neg.csv", 
#                                     polarity = "negative", ce = 'all',
#                                     column = "rp",
#                                     total.score.tol = 0.5,
#                                     candidate.num = 3,
#                                     threads = 3, 
#                                     database = "list2_ms1_database")
# 
# result3_neg <- identify_metabolites(ms1.data = "variable_info_neg.csv", 
#                                     polarity = "negative",
#                                     ce = 'all',
#                                     column = "rp",
#                                     total.score.tol = 0.5,
#                                     candidate.num = 3,
#                                     threads = 3, 
#                                     database = "nonspecific_biomarkers_ms1_database")
# 
# result4_neg <- identify_metabolites(ms1.data = "variable_info_neg.csv",
#                                     polarity = "negative",
#                                     ce = 'all',
#                                     column = "rp",
#                                     total.score.tol = 0.5,
#                                     candidate.num = 3,
#                                     threads = 3,
#                                     database = "select_exposome_ms1_database")
# # load("select_exposome")
# # result4_neg <- mz_match(ms1.table = variable_info_neg, 
# #                         database = select_exposome, 
# #                         mz.error.tol = 25)
# 
# result5_neg <- identify_metabolites(ms1.data = "variable_info_neg.csv", 
#                                     polarity = "negative",
#                                     ce = 'all',
#                                     column = "rp",
#                                     total.score.tol = 0.5,
#                                     candidate.num = 3,
#                                     threads = 3, 
#                                     database = "t3db_ms1_database")
# 
# result6_neg <- identify_metabolites(ms1.data = "variable_info_neg.csv", 
#                                     polarity = "negative",
#                                     ce = 'all',
#                                     column = "rp",
#                                     total.score.tol = 0.5,
#                                     candidate.num = 3,
#                                     threads = 3, 
#                                     database = "specific_biomarker_ms1_database")
# 
# result7_neg <- identify_metabolites(ms1.data = "variable_info_neg.csv", 
#                                     polarity = "negative",
#                                     ce = 'all',
#                                     column = "rp",
#                                     total.score.tol = 0.5,
#                                     candidate.num = 3,
#                                     threads = 3, 
#                                     database = "hmdbMS1Database0.0.1")
# 
# 
# result8_neg <- identify_metabolites(ms1.data = "variable_info_neg_hilic.csv", 
#                                     polarity = "negative",
#                                     ce = 'all',
#                                     rt.match.tol = 30,
#                                     column = "hilic",
#                                     total.score.tol = 0.5,
#                                     candidate.num = 3,
#                                     threads = 3, 
#                                     database = "msDatabase_hilic0.0.2")
# 
# result9_neg <- identify_metabolites(ms1.data = "variable_info_neg_rplc.csv", 
#                                     polarity = "negative",
#                                     ce = 'all',
#                                     rt.match.tol = 30,
#                                     column = "rp",
#                                     total.score.tol = 0.5,
#                                     candidate.num = 3,
#                                     threads = 3, 
#                                     database = "msDatabase_rplc0.0.2")
# 
# 
# annotation_table_neg1 <- 
#   get_identification_table(result1_neg,
#                            result2_neg,
#                            result4_neg,
#                            result5_neg,
#                            result6_neg,
#                            result7_neg,
#                            type = "old", 
#                            candidate.num = 1)
# 
# 
# annotation_table_neg1$Identification <-
#   annotation_table_neg1$Identification %>% 
#   lapply(function(x){
#     if(is.na(x)){
#       return(NA)
#     }else{
#       x <- stringr::str_split(x, "\\{\\}")[[1]] 
#       x <- grep("\\(M\\-H\\)|\\(2M|H2O", x, value = TRUE)
#       if(length(x) == 0){
#         return(NA)
#       }else{
#         paste(x, collapse = "{}")
#       }
#     }
#   }) %>% 
#   unlist()
# 
# annotation_table_neg1 <- 
#   metID::trans2newStyle(identification.table = annotation_table_neg1)
# 
# annotation_table_neg1 <- 
#   annotation_table_neg1 %>% 
#   dplyr::filter(!is.na(Compound.name))
# 
# annotation_table_neg2 <- 
#   get_identification_table(result8_neg,
#                            type = "old", 
#                            candidate.num = 1)
# 
# annotation_table_neg2$Identification <-
#   annotation_table_neg2$Identification %>% 
#   lapply(function(x){
#     if(is.na(x)){
#       return(NA)
#     }else{
#       x <- stringr::str_split(x, "\\{\\}")[[1]] 
#       x <- grep("\\(M\\-H\\)|\\(2M|H2O", x, value = TRUE)
#       if(length(x) == 0){
#         return(NA)
#       }else{
#         paste(x, collapse = "{}")
#       }
#     }
#   }) %>% 
#   unlist()
# 
# annotation_table_neg2 <- 
#   metID::trans2newStyle(identification.table = annotation_table_neg2)
# 
# annotation_table_neg2 <- 
#   annotation_table_neg2 %>% 
#   dplyr::filter(!is.na(Compound.name))
# 
# 
# annotation_table_neg3 <- 
#   get_identification_table(result9_neg,
#                            type = "old", 
#                            candidate.num = 1)
# 
# 
# annotation_table_neg3$Identification <-
#   annotation_table_neg3$Identification %>% 
#   lapply(function(x){
#     if(is.na(x)){
#       return(NA)
#     }else{
#       x <- stringr::str_split(x, "\\{\\}")[[1]] 
#       x <- grep("\\(M\\-H\\)|\\(2M|H2O", x, value = TRUE)
#       if(length(x) == 0){
#         return(NA)
#       }else{
#         paste(x, collapse = "{}")
#       }
#     }
#   }) %>% 
#   unlist()
# 
# annotation_table_neg3 <- 
#   metID::trans2newStyle(identification.table = annotation_table_neg3)
# 
# 
# annotation_table_neg3 <- 
#   annotation_table_neg3 %>% 
#   dplyr::filter(!is.na(Compound.name))
# 
# 
# # result4_neg <- 
# #   result4_neg %>% 
# #   dplyr::filter(!is.na(Compound.name))
# # 
# # annotation_table_neg1 <- 
# #   annotation_table_neg1 %>% 
# #   dplyr::filter(!name %in% result4_neg$name) 
# # 
# # 
# # annotation_table_neg1 <- 
# #   rbind(annotation_table_neg1, result4_neg)
# 
# annotation_table_neg1 <- 
#   annotation_table_neg1 %>% 
#   dplyr::filter(!name %in% annotation_table_neg2$name) %>% 
#   dplyr::filter(!name %in% annotation_table_neg3$name)
# 
# 
# 
# annotation_table_neg <- 
#   rbind(annotation_table_neg1[,-4],
#         annotation_table_neg2, 
#         annotation_table_neg3)
# 
# write.csv(annotation_table_neg,
#           file = "annotation_table_neg.csv", 
#           row.names = FALSE)


###internal exposome is from peng
internal_exposome <- readxl::read_xlsx("Internal exposome (1).xlsx", sheet = 2)






