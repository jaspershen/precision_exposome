sxtTools::setwd_project()
setwd("data_20200302/exposome_data/")
library(tidyverse)
metabolite_table <-
  readxl::read_xlsx("SmallMoleculeProfiling_20161129Summary_neg2.xlsx", sheet = 1)
peak_table <-
  readxl::read_xlsx("SmallMoleculeProfiling_20161129Summary_neg2.xlsx", sheet = 2)
phenotype_table <-
  readxl::read_xlsx("SmallMoleculeProfiling_20161129Summary_neg2.xlsx", sheet = 3)

dim(metabolite_table)

colnames(metabolite_table)[1] <- "peak_ID"

idx <- apply(metabolite_table, 1, function(x) {
  all(is.na(x))
}) %>%
  which()

idx2 <- c(1, idx, idx + 1, nrow(metabolite_table)) %>%
  sort()

idx2 <- matrix(idx2, ncol = 2, byrow = TRUE)

metabolite_table2 <-
  apply(idx2, 1, function(x) {
    x <- metabolite_table[x[1]:x[2], ]
    x <-
      x %>%
      dplyr::filter(!is.na(peak_ID))
    x
    x <-
      x %>%
      mutate(Note = x$peak_ID[1]) %>%
      dplyr::filter(!is.na(MetabID)) %>%
      dplyr::select(peak_ID, Note, everything())
  }) %>%
  bind_rows()


metabolite_tags <-
  metabolite_table2 %>%
  select(peak_ID:Fold_Enrich)

metabolite_table <-
  metabolite_table2 %>%
  select(-c(peak_ID:Fold_Enrich))

save(metabolite_tags, file = "metabolite_tags")
save(metabolite_table, file = "metabolite_table")

phenotype_table

colnames(phenotype_table)[1] <- "sample_ID"

##add the latitude and longitude for different locations
phenotype_table
colnames(phenotype_table)[1] <- "sample_id"
colnames(phenotype_table)[2] <- "start_date"
colnames(phenotype_table)[3] <- "end_date"
colnames(phenotype_table)[4] <- "location"

latitude <- c(
  37.408899,
  37.408899,
  37.408899,
  37.408899,
  37.408899,
  38.106940,
  38.106940,
  38.106940,
  38.106940,
  38.106940,
  46.601860,
  29.764424,
  37.408899,
  37.408899,
  37.408899,
  42.362238,
  37.408899,
  41.876097,
  38.552765,
  38.984273,
  42.362238
)

longitude <- c(
  -122.149814,
  -122.149814,
  -122.149814,
  -122.149814,
  -122.149814,
  -122.564569,
  -122.564569,
  -122.564569,
  -122.564569,
  -122.564569,
  -112.038093,
  -95.368194,
  -122.149814,
  -122.149814,
  -122.149814,
  -71.056961,
  -122.149814,
  -87.627601,
  -121.451996,
  -77.094888,
  -71.056961
)


phenotype_table <-
  data.frame(phenotype_table,
             longitude,
             latitude,
             stringsAsFactors = FALSE)


save(phenotype_table, file = "phenotype_table")




