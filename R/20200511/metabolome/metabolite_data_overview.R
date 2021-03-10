##avoid source
no_function()

sxtTools::setwd_project()
rm(list=ls())

load("data_20200511/metabolome/sample_info")
load("data_20200511/metabolome/expression_data")
load("data_20200511/metabolome/variable_info")
setwd("data_20200511/metabolome/")

# openxlsx::write.xlsx(sample_info, file = "metabolome_sample_info.xlsx")
# 
# openxlsx::write.xlsx(expression_data, file = "metabolome_expression_data.xlsx")
# 
# openxlsx::write.xlsx(variable_info, file = "metabolome_variable_info.xlsx")


dim(sample_info)
dim(variable_info)

library(ggplot2)

plot <- 
sample_info %>% 
  ggplot(aes(as.Date(CollectionDate), y = subject_id)) +
  geom_point(size = 5, color = "#80B1D3") +
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

ggsave(plot, file = "metabolome_sample_collection.pdf", 
       width = 9, height = 7)


ggsave(plot, file = "metabolome_sample_collection.png", 
       width = 9, height = 7)

