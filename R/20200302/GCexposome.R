sxtTools::setwd_project()
setwd("data_20200302/GCexposome/")
dir()

gcexposome_data <-
  readxl::read_xlsx("GC_Orbitrap_Exposomics_Chemical_List_27Sept2018.xlsx")

load("HMDB.metabolite.data")

colnames(gcexposome_data)[4] <- "Compound.name"
colnames(gcexposome_data)[7] <- "mz"
colnames(gcexposome_data)[3] <- "PubChem.ID"
colnames(gcexposome_data)[5] <- "Formula"

gcexposome_data <- 
  gcexposome_data %>% 
  dplyr::select(Compound.name, mz, everything()) %>% 
  mutate(KEGG.ID = NA, HMDB.ID = NA, CAS.ID = NA)


gcExposome.data <- gcexposome_data %>% 
  as.data.frame()

save(gcExposome.data, file = "gcExposome.data", compress = "xz")


sxtTools::setwd_project()
setwd("data_20200302/GCexposome/")
peak_table <- readxl::read_xlsx("SmallMoleculeProfiling_20161129Summary_neg2.xlsx", sheet = 2)
peak_table <- peak_table[-1,]
peak_table$MetID[1:5]

peak_table <- 
  peak_table %>% 
  dplyr::select(name = MetID, mz, rt)

peak_table_pos <- 
  peak_table %>% 
  dplyr::filter(str_detect(name, "PM"))

peak_table_neg <- 
  peak_table %>% 
  dplyr::filter(str_detect(name, "NM"))

write.csv(peak_table_pos, "peak_table_pos.csv", row.names = FALSE)
write.csv(peak_table_neg, "peak_table_neg.csv", row.names = FALSE)

load("gcExposome.data")


tibble::as_tibble(gcExposome.data)

result_pos <- metID::mzIdentify(ms1.data = "peak_table_pos.csv", 
                                ms1.match.ppm = 15, 
                                polarity = "positive",
                                column = "rp", 
                                candidate.num = 3, 
                                database = "gcExposome.data")

result_neg <- metID::mzIdentify(ms1.data = "peak_table_neg.csv", 
                                ms1.match.ppm = 15, 
                                polarity = "negative",
                                column = "rp", 
                                candidate.num = 3, 
                                database = "gcExposome.data")

annotation_table_pos <- 
  metID::getIdentificationTable2(object = result_pos, candidate.num = 1, type = "new")


annotation_table_neg <- 
  metID::getIdentificationTable2(object = result_neg, candidate.num = 1, type = "new")



##only remain M+H, M-H, 2M+H and 2M-H

unique(annotation_table_pos$Adduct)

name_pos <- 
  annotation_table_pos %>% 
    dplyr::filter(!is.na(Adduct)) %>% 
    dplyr::filter(Adduct == "(M+H)+" | Adduct == "(2M+H)+") %>% 
    pull(name)
  

unique(annotation_table_neg$Adduct)
name_neg <-  
annotation_table_neg %>% 
    dplyr::filter(!is.na(Adduct)) %>% 
    dplyr::filter(Adduct == "(M-H)-" | Adduct == "(2M-H)-") %>% 
    pull(name)
  

annotation_table_pos[-match(name_pos, annotation_table_pos$name),-c(1,2,3)] <- NA

annotation_table_neg[-match(name_neg, annotation_table_neg$name),-c(1,2,3)] <- NA


write.csv(annotation_table_pos, "annotation_table_pos_gcExposome.csv", row.names = FALSE)
write.csv(annotation_table_neg, "annotation_table_neg_gcExposome.csv", row.names = FALSE)
