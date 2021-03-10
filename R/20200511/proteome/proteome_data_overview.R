sxtTools::setwd_project()
rm(list=ls())

load("data_20200511/proteome/sample_info")
load("data_20200511/proteome/expression_data")
load("data_20200511/proteome/variable_info")
setwd("data_20200511/proteome/")

dim(sample_info)
dim(variable_info)

variable_info$protein_id

###Nuclear factorerythroid 2-related factor 2 (NRF2)
grep("NFE2L2", variable_info$protein_id)
grep("NRF2", variable_info$protein_id)

##Kelch-like ECH-associated protein 1 (KEAP1)
grep("KEAP1", variable_info$protein_id)
grep("INRF2", variable_info$protein_id)
grep("KIAA0132", variable_info$protein_id)
grep("KLHL19", variable_info$protein_id)

##Nuclear factor kB (NF-kB)
grep("NF-kB", variable_info$protein_id)
grep("NF", variable_info$protein_id, value = TRUE)

library(ggplot2)

plot <- 
sample_info %>% 
  ggplot(aes(as.Date(CollectionDate), y = subject_id)) +
  geom_point(size = 5, color = "#008280FF") +
  theme_bw() +
  scale_x_continuous(trans = "date",
                     breaks = c(as.Date(sample_info$CollectionDate)),
                     labels = as.character(sample_info$CollectionDate)
  ) +
  labs(x = "", y = "") +
  theme(axis.title = element_text(size = 10),
        axis.text.x = element_text(size = 10, 
                                   angle = 45,
                                   vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 10),
        # panel.grid = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "pt")
  )


plot

# ggsave(plot, file = "proteome_sample_collection.pdf", 
#        width = 9, height = 7)




