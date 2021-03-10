##avoid source
no_function()

##
sxtTools::setwd_project()
library(tidyverse)
load("data_20200511/microbiome/dna_sample_info")
dna_sample_info <- dna_sample_info[,1:16]

sxtTools::setwd_project()
setwd("data_20200511/lab_test/")

sample_info <- readxl::read_xlsx("Exposome_select_peng_metadata-DNA.xlsx")

expression_data <-
  sample_info[, c(1, 14:57)] %>%
  as.data.frame() %>% 
  dplyr::filter(!is.na(STARTING_DATE))

rownames(expression_data) <- as.character(expression_data$STARTING_DATE)
expression_data <- 
  expression_data %>% 
  dplyr::select(-STARTING_DATE)

variable_info <-
  sample_info[1:2, -1] %>%
  t() %>%
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "variable_id") %>% 
  dplyr::filter(variable_id %in% colnames(expression_data))

colnames(variable_info)[2:3] = c("true_name", "unit")

sample_info <- 
  sample_info %>% 
  dplyr::select(STARTING_DATE) %>% 
  dplyr::filter(!is.na(STARTING_DATE)) %>% 
  dplyr::mutate(sample_id = as.character(STARTING_DATE)) %>% 
  dplyr::select(sample_id, everything())

expression_data <- 
  t(expression_data) %>% 
  as.data.frame()

expression_data <-
  expression_data %>% 
  apply(2, function(x){
    x <- as.numeric(x)
  }) %>% 
  as.data.frame()

rownames(expression_data) <- variable_info$variable_id

sum(is.na(expression_data))

expression_data

# expression_data <-
#   impute::impute.knn(data = as.matrix(expression_data))$data %>%
#   as.data.frame()

sample_info$sample_id == dna_sample_info$date.start

sample_info <- 
  cbind(sample_info,
        dna_sample_info[,-1])

sample_info$location[sample_info$location == "Mike_background"] <- "Campus"

remove_idx <- 
expression_data %>% 
  purrr::map(function(x){
    all(is.na(x))
  }) %>% 
  unlist() %>% 
  which()

if(length(remove_idx) > 0){
  expression_data <-
    expression_data[,-remove_idx]
  sample_info <-
    sample_info[-remove_idx,]
}

sum(is.na(expression_data))

save(expression_data, file = "expression_data")
save(sample_info, file = "sample_info")
save(variable_info, file = "variable_info")


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
saveWorkbook(wb, "blood_test_data.xlsx", overwrite = TRUE)
