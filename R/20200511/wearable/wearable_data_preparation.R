##avoid source
no_function()

sxtTools::setwd_project()
library(tidyverse)
setwd("data_20200511/wearable/")
rm(list = ls())

anatomy_history <- readr::read_csv("anatomy-history.csv")


anatomy_history <- 
anatomy_history %>% 
  dplyr::filter(`Timestamp (local)` > "2015-12-1" & `Timestamp (local)` < "2016-5-30")

data_metrics <- readr::read_csv("day-metrics.csv")

data_metrics <- 
  data_metrics %>% 
  dplyr::filter(Date > "2015-12-1" & Date < "2016-5-30")


dim(anatomy_history)
dim(data_metrics)


colnames(anatomy_history)
colnames(data_metrics)


anatomy_history %>% 
  ggplot(aes(`Timestamp (UTC)`, Weight)) +
  geom_point()

anatomy_history <- 
  anatomy_history %>% 
  dplyr::mutate(Date = as.Date(`Timestamp (UTC)`)) %>% 
  dplyr::select(Date, everything())

expression_data <- 
  anatomy_history %>% 
  dplyr::left_join(data_metrics, by = c("Date"))

colnames(expression_data)

expression_data <-
  expression_data %>% 
  dplyr::distinct(Date, .keep_all = TRUE)

expression_data <- 
  expression_data[,c(1,4,5,8:23)]

variable_info <-
  data.frame(variable_id = colnames(expression_data)[-1],
             stringsAsFactors = FALSE)

sample_info <- 
  expression_data %>% 
  dplyr::select(Date) %>% 
  dplyr::mutate(sample_id = paste("sample", 1:nrow(expression_data), sep = "_"))

rownames(sample_info) <- sample_info$sample_id

expression_data <- 
expression_data %>% 
  dplyr::select(-Date)

expression_data <- 
  expression_data %>% 
  t() %>% 
  as.data.frame()

rownames(expression_data) == variable_info$variable_id

colnames(expression_data) <- sample_info$sample_id

save(expression_data, file = "expression_data")
save(variable_info, file = "variable_info")
save(sample_info, file = "sample_info")














