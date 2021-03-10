sxtTools::setwd_project()
setwd("data_20200302/metabolomics_data/")
load("clinic_data")
load("met_data")
load("met_tag")


met_data <-
  met_data %>% 
  dplyr::filter(SubjectID == "69-001")

save(met_data, file = "met_data")
save(met_tag, file = "met_tag")

clinic_data <- 
  clinic_data %>% 
  dplyr::filter(SubjectID == "69-001")

save(clinic_data, file = "clinic_data")

met_data

met_tag$Metabolite







