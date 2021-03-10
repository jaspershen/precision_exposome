##avoid source
no_function()

sxtTools::setwd_project()
library(tidyverse)
setwd("data_20200511/proteome/")
rm(list = ls())

load("proteome_data")

proteome_data <-
  proteome_data %>%
  dplyr::filter(SubjectID == "69-001") %>%
  dplyr::filter(CollectionDate > "2015-12-01",
                CollectionDate < "2016-06-01")

variable_info <- 
  colnames(proteome_data)[-c(1:7)]

variable_info <- data.frame(protein_id = variable_info, 
                            stringsAsFactors = FALSE)

expression_data <- proteome_data[,-c(1:7)]

sample_info <- proteome_data[,c(1:7)]

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
  dplyr::filter(protein_id %in% rownames(expression_data))

variable_info$protein_id == rownames(expression_data)

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
saveWorkbook(wb, "proteome_data.xlsx", overwrite = TRUE)



