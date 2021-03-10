##avoid source
no_function()

sxtTools::setwd_project()
library(tidyverse)
setwd("data_20200511/nasal_microbiome/")
rm(list = ls())

load("ns_data")

nasal_microbiome_data <-
  ns_data %>%
  dplyr::filter(SubjectID == "69-001") %>%
  dplyr::filter(CollectionDate > "2015-12-01",
                CollectionDate < "2016-06-01")

colnames(ns_data)

variable_info <- 
  colnames(nasal_microbiome_data)[-c(1:8)]

variable_info <- data.frame(variable_id = variable_info,
                            stringsAsFactors = FALSE)


variable_info$short_name <-
  stringr::str_split(variable_info$variable_id, "_", n = 2) %>%
  purrr::map(
    .f = function(x) {
      x[2]
    }
  ) %>%
  unlist()

variable_info$level <-
  stringr::str_split(variable_info$variable_id, "_", n = 2) %>%
  purrr::map(
    .f = function(x) {
      x[1]
    }
  ) %>%
  unlist()
  
expression_data <- nasal_microbiome_data[,-c(1:8)]

sample_info <- nasal_microbiome_data[,c(1:8)]

sample_info <- 
  sample_info %>% 
  dplyr::select(sample_id = SampleID, 
                subject_id = SubjectID, 
                CollectionDate,
                CL1, CL2, CL3, CL4)

expression_data <-
  t(expression_data) %>% 
  as.data.frame()

colnames(expression_data) <- sample_info$sample_id

rownames(expression_data)

variable_info <- 
  variable_info %>% 
  dplyr::filter(variable_id %in% rownames(expression_data))

variable_info$variable_id == rownames(expression_data)

save(variable_info, file = "variable_info")
save(sample_info, file = "sample_info")
save(expression_data, file = "expression_data")






