##avoid source
no_exist_function()

sxtTools::setwd_project()
rm(list = ls())
library(tidyverse)

# load data
###fiber information
fiber_sample_info <-
  readxl::read_xlsx("data_20200511/fiber/External-Internal Samples and Dates 19-11-06.xlsx", 
                    col_types = c("text", "date", "date", "text", "text", "text", "text", "text",
                                  "text", "text", "numeric", "numeric", "numeric",
                                  "numeric", "numeric", "numeric", "numeric", "numeric"))

fiber_sample_info = 
fiber_sample_info[2:15,] %>% 
  dplyr::select(start_date = `filter start`, 
                end_date = `filter collect`, 
                location = comments,
                class = Class) %>% 
  dplyr::mutate(start_date = as.Date(start_date),
                end_date = as.Date(end_date))

##exposome chemical
load("data_20200511/exposome/sample_info")
load("data_20200511/exposome/expression_data")

exposomeChemical_sample_info <- sample_info
exposomeChemical_expression_data <- expression_data

exposomeChemical_expression_data <-
  exposomeChemical_expression_data %>%
  t() %>% 
  as.data.frame() %>%
  rownames_to_column(var = "Sample_ID") %>%
  mutate(Sample_ID = stringr::str_replace(Sample_ID, "_[0-9]{1,2}", "")) %>%
  plyr::dlply(.variables = "Sample_ID") %>%
  lapply(function(x) {
    apply(x[, -1], 2, mean)
  }) %>%
  do.call(rbind, .) %>% 
  t() %>% 
  as.data.frame()

##metabolome
load("data_20200511/metabolome/sample_info")
load("data_20200511/metabolome/expression_data")
metabolome_sample_info <- sample_info
metabolome_expression_data <- expression_data

##proteome
load("data_20200511/proteome/sample_info")
load("data_20200511/proteome/expression_data")
proteome_sample_info <- sample_info
proteome_expression_data <- expression_data

##microbiome
load("data_20200511/microbiome/dna_sample_info")
load("data_20200511/microbiome/dna_expression_data")
exposomeBiological_sample_info <- dna_sample_info
exposomeBiological_expression_data <- dna_expression_data
exposomeBiological_expression_data <- 
  exposomeBiological_expression_data[grep("phylum", rownames(exposomeBiological_expression_data)),]
  
##cytokine
load("data_20200511/cytokine/sample_info")
load("data_20200511/cytokine/expression_data")
cytokine_sample_info <- sample_info
cytokine_expression_data <- expression_data

##gut microbiome
load("data_20200511/gut_microbiome/sample_info")
load("data_20200511/gut_microbiome/expression_data")
gutmicrobiome_sample_info <- sample_info
gutmicrobiome_expression_data <- expression_data

gutmicrobiome_sample_info <-
  gutmicrobiome_sample_info %>% 
  dplyr::distinct(CollectionDate, .keep_all = TRUE)

gutmicrobiome_expression_data <-
  gutmicrobiome_expression_data[,gutmicrobiome_sample_info$sample_id]

##environment
load("data_20200511/environment/sample_info")
load("data_20200511/environment/expression_data")
environment_sample_info <- sample_info
environment_expression_data <- expression_data

##lab test
load("data_20200511/lab_test/sample_info")
load("data_20200511/lab_test/expression_data")
lab_sample_info <- sample_info
lab_expression_data <- expression_data

setwd("data_analysis/study_information")

###fiber
fiber_sample_info1 =
  fiber_sample_info %>% 
  dplyr::filter(start_date >= "2016-01-12" & start_date <= "2016-03-03") %>% 
  dplyr::distinct(end_date, .keep_all = TRUE)

##exposome chemical
exposomeChemical_sample_info1 <- 
  exposomeChemical_sample_info %>% 
  dplyr::select(start_date = start_date, end_date = end_date, location = comments) %>% 
  dplyr::mutate(start_date = as.Date(start_date),
                end_date = as.Date(end_date),
                class = "exp") %>% 
  dplyr::filter(start_date >= "2016-01-12" & start_date <= "2016-03-03") %>% 
  dplyr::distinct(end_date, .keep_all = TRUE)

####exposome toxins and carcinogens
toxin_sample_info1 <- 
  metabolome_sample_info %>% 
  dplyr::select(start_date = CollectionDate,
                end_date = CollectionDate) %>% 
  dplyr::mutate(start_date = as.Date(start_date),
                end_date = as.Date(end_date),
                location = NA,
                class = "toxin") %>% 
  dplyr::filter(start_date >= "2016-01-12" & start_date <= "2016-03-03") %>% 
  dplyr::distinct(end_date, .keep_all = TRUE)

##metabolome
metabolome_sample_info1 <- 
  metabolome_sample_info %>% 
  dplyr::select(start_date = CollectionDate,
                end_date = CollectionDate) %>% 
  dplyr::mutate(start_date = as.Date(start_date),
                end_date = as.Date(end_date),
                location = NA,
                class = "met") %>% 
  dplyr::filter(start_date >= "2016-01-12" & start_date <= "2016-03-03") %>% 
  dplyr::distinct(end_date, .keep_all = TRUE)

##proteome
proteome_sample_info1 <- 
  proteome_sample_info %>% 
  dplyr::select(start_date = CollectionDate,
                end_date = CollectionDate) %>% 
  dplyr::mutate(start_date = as.Date(start_date),
                end_date = as.Date(end_date),
                location = NA,
                class = "pro") %>% 
  dplyr::filter(start_date >= "2016-01-12" & start_date <= "2016-03-03") %>% 
  dplyr::distinct(end_date, .keep_all = TRUE)

##microbiome
exposomeBiological_sample_info1 <- 
  exposomeBiological_sample_info %>% 
  dplyr::select(start_date = date.start,
                end_date = date.end,
                location) %>% 
  dplyr::mutate(start_date = as.Date(start_date),
                end_date = as.Date(end_date),
                class = "micro") %>% 
  dplyr::filter(start_date >= "2016-01-12" & start_date <= "2016-03-03") %>% 
  dplyr::distinct(end_date, .keep_all = TRUE)

##cytokine
cytokine_sample_info1 <- 
  cytokine_sample_info %>% 
  dplyr::select(start_date = CollectionDate,
                end_date = CollectionDate) %>% 
  dplyr::mutate(start_date = as.Date(start_date),
                end_date = as.Date(end_date),
                location = NA,
                class = "cytokine") %>% 
  dplyr::filter(start_date >= "2016-01-12" & start_date <= "2016-03-03") %>% 
  dplyr::distinct(end_date, .keep_all = TRUE)

##gut microbiome
gutmicrobiome_sample_info1 <- 
  gutmicrobiome_sample_info %>% 
  dplyr::select(start_date = CollectionDate,
                end_date = CollectionDate) %>% 
  dplyr::mutate(start_date = as.Date(start_date),
                end_date = as.Date(end_date),
                location = NA,
                class = "gutmicrobiome") %>% 
  dplyr::filter(start_date >= "2016-01-12" & start_date <= "2016-03-03") %>% 
  dplyr::distinct(end_date, .keep_all = TRUE)

##environment
environment_sample_info1 <- 
  environment_sample_info %>% 
  dplyr::select(start_date = date.start,
                end_date = date.end,
                location) %>% 
  dplyr::mutate(start_date = as.Date(start_date),
                end_date = as.Date(end_date),
                class = "environment") %>% 
  dplyr::filter(start_date >= "2016-01-12" & start_date <= "2016-03-03") %>% 
  dplyr::distinct(end_date, .keep_all = TRUE)

##lab test
lab_sample_info1 <- 
  lab_sample_info %>% 
  dplyr::select(start_date = date.start,
                end_date = date.end,
                location
                ) %>% 
  dplyr::mutate(start_date = as.Date(start_date),
                end_date = as.Date(end_date),
                class = "lab") %>% 
  dplyr::filter(start_date >= "2016-01-12" & start_date <= "2016-03-03") %>% 
  dplyr::distinct(end_date, .keep_all = TRUE)



class_color <- 
  c("7Environment" = ggsci::pal_d3()(10)[1],
    "9Exposome (chemical)" = ggsci::pal_d3()(10)[2],
    "6Metabolome" = ggsci::pal_d3()(10)[3],
    "4Proteome" = ggsci::pal_d3()(10)[4],
    "8Exposome (biological)" = ggsci::pal_d3()(10)[5],
    "1Gut microbiome" = ggsci::pal_d3()(10)[6],
    "3Blood test" = ggsci::pal_d3()(10)[7],
    "2Cytokine" = ggsci::pal_d3()(10)[8],
    "5Toxins and carcinogens" = ggsci::pal_d3()(10)[9]
  )



all_date = 
  unique(c(exposomeChemical_sample_info1$start_date,
           exposomeChemical_sample_info1$end_date,
           exposomeBiological_sample_info1$start_date,
           exposomeBiological_sample_info1$end_date,
           environment_sample_info1$start_date,
           environment_sample_info1$end_date,
           metabolome_sample_info1$start_date,
           proteome_sample_info1$start_date,
           gutmicrobiome_sample_info1$start_date,
           cytokine_sample_info1$start_date,
           lab_sample_info1$start_date
           )) %>% 
  sort()


temp_data1 = 
  rbind(exposomeChemical_sample_info1 %>% 
          dplyr::mutate(image  = NA),
        exposomeBiological_sample_info1 %>% 
          dplyr::mutate(image = NA),
        environment_sample_info1 %>% 
          dplyr::mutate(image = NA),
        metabolome_sample_info1 %>% 
          dplyr::mutate(image = "metabolome.png"),
        toxin_sample_info1 %>% 
          dplyr::mutate(image = "toxin.png"),
        proteome_sample_info1 %>% 
          dplyr::mutate(image = "proteome.png"),
        lab_sample_info1 %>% 
          dplyr::mutate(image = "lab.png"),
        cytokine_sample_info1 %>% 
          dplyr::mutate(image = "cytokine.png"),
        gutmicrobiome_sample_info1 %>% 
          dplyr::mutate(image = "microbiome.png")
        ) %>% dplyr::mutate(
    class = case_when(
      class == "exp" ~ "9Exposome (chemical)",
      class == "micro" ~ "8Exposome (biological)",
      class == "environment" ~ "7Environment",
      class == "met" ~ "6Metabolome",
      class == "toxin" ~ "5Toxins and carcinogens",
      class == "pro" ~ "4Proteome",
      class == "lab" ~ "3Blood test",
      class == "cytokine" ~ "2Cytokine",
      class == "gutmicrobiome" ~ "1Gut microbiome"
    ) 
  ) 

temp_data2 = 
  rbind(
    fiber_sample_info1 %>% 
      dplyr::select(start_date, end_date, location) %>% 
      dplyr::mutate(class = location,
                    image = NA,
                    ymin = 0, 
                    ymax = 1),
    fiber_sample_info1 %>% 
      dplyr::select(start_date, end_date, location = class) %>% 
      dplyr::mutate(class = location,
                    image = NA,
                    ymin = 1, 
                    ymax = 2)
  )

library(ggimage)


plot = 
ggplot(temp_data1) +
  geom_segment(
    aes(
      x = start_date,
      xend = end_date,
      y = class,
      yend = class,
      color = class
    ),
    arrow = arrow(
      angle = 90,
      type = "open",
      ends = "both",
      length = unit(x = 0.2, units = "cm")
    ),
    show.legend = FALSE,
    data = temp_data1 %>%
      dplyr::filter(
        class %in% c(
          "9Exposome (chemical)",
          "8Exposome (biological)",
          "7Environment"
        )
      )
  ) +
  geom_image(aes(x = start_date,
                 y = class,
                 image = image),
             data = temp_data1 %>%
               dplyr::filter(
                 class %in% c(
                   "6Metabolome",
                   "5Toxins and carcinogens",
                   "4Proteome",
                   "3Blood test",
                   "2Cytokine",
                   "1Gut microbiome"
                 )
               )) +
  geom_rect(
    aes(
      xmin = start_date,
      xmax = end_date,
      ymin = ymin - 2,
      ymax = ymax - 2,
      fill = class
    ),
    color = "black",
    data = temp_data2,
    show.legend = TRUE
  ) +
  scale_x_date(breaks = all_date,
               labels = all_date,
               expand = expansion(mult = c(0.05, 0))) +
  scale_color_manual(values = class_color) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  ) +
  scale_fill_manual(values = c(
    "NO" = "grey", 
    "arabinoxylan" = "skyblue", 
    "guar gum" = "skyblue",
    "Campus" = ggsci::pal_aaas()(n=10)[1],
    "SF and Novato, CA" = ggsci::pal_aaas()(n=10)[2],
    "Montana" = ggsci::pal_aaas()(n=10)[3],
    "Houston" = ggsci::pal_aaas()(n=10)[4]
  )) + 
  labs(x = "", y = "")


ggsave(plot = plot, filename = "detailed_info_plot.pdf", width = 10, height = 7)

