##avoid source
no_exist_function()

sxtTools::setwd_project()
rm(list = ls())
library(tidyverse)

# load data
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


###begin circlize
library(circlize)

##exposome chemical
exposomeChemical_sample_info1 <- 
  exposomeChemical_sample_info %>% 
  dplyr::select(date = start_date) %>% 
  dplyr::mutate(date = as.Date(date),
                class = "exp")

####exposome toxins and carcinogens
toxin_sample_info1 <- 
  metabolome_sample_info %>% 
  dplyr::select(date = CollectionDate) %>% 
  dplyr::mutate(date = as.Date(date),
                class = "toxin")

##metabolome
metabolome_sample_info1 <- 
  metabolome_sample_info %>% 
  dplyr::select(date = CollectionDate) %>% 
  dplyr::mutate(date = as.Date(date),
                class = "met")

##proteome
proteome_sample_info1 <- 
  proteome_sample_info %>% 
  dplyr::select(date = CollectionDate) %>% 
  dplyr::mutate(date = as.Date(date),
                class = "pro")

##microbiome
exposomeBiological_sample_info1 <- 
  exposomeBiological_sample_info %>% 
  dplyr::select(date = date.start) %>% 
  dplyr::mutate(date = as.Date(date),
                class = "micro")

##cytokine
cytokine_sample_info1 <- 
  cytokine_sample_info %>% 
  dplyr::select(date = CollectionDate) %>% 
  dplyr::mutate(date = as.Date(date),
                class = "cytokine")

##gut microbiome
gutmicrobiome_sample_info1 <- 
  gutmicrobiome_sample_info %>% 
  dplyr::select(date = CollectionDate) %>% 
  dplyr::mutate(date = as.Date(date),
                class = "gutmicrobiome")

##environment
environment_sample_info1 <- 
  environment_sample_info %>% 
  dplyr::select(date = date.start) %>% 
  dplyr::mutate(date = as.Date(date),
                class = "environment")

##lab test
lab_sample_info1 <- 
  lab_sample_info %>% 
  dplyr::select(date = date.start) %>% 
  dplyr::mutate(date = as.Date(date),
                class = "lab")

all_date <- 
  data.frame(date = unique(c(
    exposomeChemical_sample_info1$date,
    metabolome_sample_info1$date,
    toxin_sample_info1$date,
    proteome_sample_info1$date,
    exposomeBiological_sample_info1$date,
    cytokine_sample_info1$date,
    gutmicrobiome_sample_info1$date,
    environment_sample_info1$date,
    lab_sample_info1$date
  ))) %>% 
  dplyr::filter(date >= "2016-01-12" & date <= "2016-03-03")
  
exposomeChemical_sample_info1 <-
  all_date %>%
  dplyr::left_join(exposomeChemical_sample_info1, by = c("date")) %>%
  dplyr::mutate(value =
                  case_when(!is.na(class) ~ "Yes",
                            TRUE ~ "No")) %>%
  dplyr::mutate(class = class[!is.na(class)][1])

toxin_sample_info1 <-
  all_date %>% 
  dplyr::left_join(toxin_sample_info1, by = c("date")) %>% 
  dplyr::mutate(value = 
                  case_when(
                    !is.na(class) ~ "Yes",
                    TRUE ~ "No"
                  )) %>% 
  dplyr::mutate(class = class[!is.na(class)][1])

metabolome_sample_info1 <-
  all_date %>% 
  dplyr::left_join(metabolome_sample_info1, by = c("date")) %>% 
  dplyr::mutate(value = 
                  case_when(
                    !is.na(class) ~ "Yes",
                    TRUE ~ "No"
                  )) %>% 
  dplyr::mutate(class = class[!is.na(class)][1])

proteome_sample_info1 <-
  all_date %>% 
  dplyr::left_join(proteome_sample_info1, by = c("date")) %>% 
  dplyr::mutate(value = 
                  case_when(
                    !is.na(class) ~ "Yes",
                    TRUE ~ "No"
                  )) %>% 
  dplyr::mutate(class = class[!is.na(class)][1])

exposomeBiological_sample_info1 <-
  all_date %>% 
  dplyr::left_join(exposomeBiological_sample_info1, by = c("date")) %>% 
  dplyr::mutate(value = 
                  case_when(
                    !is.na(class) ~ "Yes",
                    TRUE ~ "No"
                  )) %>% 
  dplyr::mutate(class = class[!is.na(class)][1])

cytokine_sample_info1 <-
  all_date %>% 
  dplyr::left_join(cytokine_sample_info1, by = c("date")) %>% 
  dplyr::mutate(value = 
                  case_when(
                    !is.na(class) ~ "Yes",
                    TRUE ~ "No"
                  )) %>% 
  dplyr::mutate(class = class[!is.na(class)][1])

gutmicrobiome_sample_info1 <-
  all_date %>% 
  dplyr::left_join(gutmicrobiome_sample_info1, by = c("date")) %>% 
  dplyr::mutate(value = 
                  case_when(
                    !is.na(class) ~ "Yes",
                    TRUE ~ "No"
                  )) %>% 
  dplyr::mutate(class = class[!is.na(class)][1])

environment_sample_info1 <-
  all_date %>% 
  dplyr::left_join(environment_sample_info1, by = c("date")) %>% 
  dplyr::mutate(value = 
                  case_when(
                    !is.na(class) ~ "Yes",
                    TRUE ~ "No"
                  )) %>% 
  dplyr::mutate(class = class[!is.na(class)][1])

lab_sample_info1 <-
  all_date %>% 
  dplyr::left_join(lab_sample_info1, by = c("date")) %>% 
  dplyr::mutate(value = 
                  case_when(
                    !is.na(class) ~ "Yes",
                    TRUE ~ "No"
                  )) %>% 
  dplyr::mutate(class = class[!is.na(class)][1])

total_sample_info <-
  rbind(
    exposomeChemical_sample_info1,
    toxin_sample_info1,
    metabolome_sample_info1,
    proteome_sample_info1,
    exposomeBiological_sample_info1,
    cytokine_sample_info1,
    gutmicrobiome_sample_info1,
    environment_sample_info1,
    lab_sample_info1
  )

library(ComplexHeatmap)

value <- 
  c("Environment" = ggsci::pal_d3()(10)[1],
    "Exposome (chemical)" = ggsci::pal_d3()(10)[2],
    "Metabolome" = ggsci::pal_d3()(10)[3],
    "Proteome" = ggsci::pal_d3()(10)[4],
    "Exposome (biological)" = ggsci::pal_d3()(10)[5],
    "Gut microbiome" = ggsci::pal_d3()(10)[6],
    "Blood test" = ggsci::pal_d3()(10)[7],
    "Cytokine" = ggsci::pal_d3()(10)[8],
    "Toxins and carcinogens" = ggsci::pal_d3()(10)[9]
    )

total_sample_info <- 
  total_sample_info %>% 
  dplyr::mutate(date = as.character(date)) %>%
  dplyr::mutate(
    class = case_when(
      class == "exp" ~ "Exposome (chemical)",
      class == "micro" ~ "Exposome (biological)",
      class == "met" ~ "Metabolome",
      class == "pro" ~ "Proteome",
      class == "lab" ~ "Blood test",
      class == "cytokine" ~ "Cytokine",
      class == "gutmicrobiome" ~ "Gut microbiome",
      class == "environment" ~ "Environment",
      class == "toxin" ~ "Toxins and carcinogens"
    )
  ) %>% 
  dplyr::mutate(class = factor(class, levels = rev(
    c(
      "Exposome (chemical)",
      "Exposome (biological)",
      "Environment",
      "Metabolome",
      "Toxins and carcinogens",
      "Proteome",
      "Blood test",
      "Cytokine",
      "Gut microbiome"
    )
  )))
  

plot1 <- 
total_sample_info %>%
  ggplot(aes(date, class)) +
  geom_tile(aes(x = date, y = class, fill = value), color = "black", show.legend = FALSE) +
  scale_y_discrete(expand = expansion(mult = c(0, 0))) +
  scale_x_discrete(expand = expansion(mult = c(0, 0))) +
  scale_fill_manual(values = c("Yes" = "red", "No" = "white")) +
  theme_bw() +
  labs(x = "", y = "") +
  theme(axis.text = element_text(
    angle = 45,
    hjust = 1,
    vjust = 1
  ),
  axis.ticks.y = element_blank(),
  plot.margin = unit(c(0, 0, 0, 0), "pt"),
  panel.grid.minor = element_blank())

plot1

plot2 <-
  total_sample_info %>% 
    dplyr::filter(value == "Yes") %>% 
  dplyr::group_by(class) %>% 
  dplyr::summarise(n = n()) %>% 
    dplyr::ungroup() %>% 
  ggplot(aes(date, class)) +
    geom_bar(aes(x = n, y = class), stat = "identity", color = "black") +
    scale_x_continuous(expand = expansion(mult = c(0, 0))) +
    scale_y_discrete(expand = expansion(mult = c(0, 0))) +
    scale_fill_manual(values = c("Yes" = "red", "No" = "white")) +
    theme_classic() +
    labs(x = "Sample number", y = "") +
    theme(axis.text.x = element_text(
      size = 12
    ),
    axis.text.y = element_blank(),
    axis.title = element_text(size = 13),
    axis.ticks.y = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "pt"),
    panel.grid.minor = element_blank())

plot2

plot3 <-
  total_sample_info  %>% 
  dplyr::filter(value == "Yes") %>% 
  dplyr::group_by(date, class) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::ungroup() %>% 
  ggplot(aes(date, n)) +
  geom_bar(aes(x = date, y = n, fill = class), 
           stat = "identity", color = "black", show.legend = FALSE) +
  guides(fill = guide_legend(nrow = 1)) +
  scale_x_discrete(expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  scale_fill_manual(values = value) +
  guides(fill = guide_legend(nrow = 2)) +
  theme_classic() +
  labs(x = "", y = "Sample number") +
  theme(axis.text.y = element_text(
    size = 12
  ),
  axis.text.x = element_blank(),
  axis.title = element_text(size = 13),
  axis.ticks.x = element_blank(),
  plot.margin = unit(c(0, 0, 0, 0), "pt"),
  panel.grid.minor = element_blank(),
  legend.position = "top")

plot3  
    
library(patchwork)

plot <-
  {
    plot3 + plot2 + plot_layout(ncol = 2, widths = c(5, 1))
  } -
  {
    plot1 + plot2 + plot_layout(ncol = 2, widths = c(5, 1))
  } +
  plot_layout(ncol = 1, heights = c(1, 5))

plot

# ggsave(plot, filename = "sample_collection_all.pdf", width = 7, height = 5)



##circos plot
library(circlize)
dim(exposomeChemical_expression_data)
dim(exposomeBiological_expression_data)
dim(environment_expression_data)
dim(metabolome_expression_data)
dim(proteome_expression_data)
dim(lab_expression_data)
dim(cytokine_expression_data)
dim(gutmicrobiome_expression_data)

colnames(exposomeChemical_expression_data)
colnames(exposomeBiological_expression_data)

toxin_expression_data <- 
  metabolome_expression_data[1:30,]

colnames(exposomeChemical_expression_data) <-
  as.character(as.Date(exposomeChemical_sample_info$start_date))

colnames(exposomeBiological_expression_data) <-
  as.character(exposomeBiological_sample_info$date.start)

colnames(environment_expression_data) <-
  as.character(environment_sample_info$date.start)

colnames(metabolome_expression_data) <- 
  colnames(toxin_expression_data) <- 
  as.character(metabolome_sample_info$CollectionDate)

colnames(proteome_expression_data) <- 
  as.character(proteome_sample_info$CollectionDate)

colnames(lab_expression_data) <-
  as.character(lab_sample_info$date.start)

colnames(cytokine_expression_data) <- 
  as.character(cytokine_sample_info$CollectionDate)

colnames(gutmicrobiome_expression_data) <-
  as.character(gutmicrobiome_sample_info$CollectionDate)

all_date <-
  Reduce(union, list(colnames(exposomeChemical_expression_data),
                     colnames(exposomeBiological_expression_data),
                     colnames(environment_expression_data),
                     colnames(metabolome_expression_data),
                     colnames(toxin_expression_data),
                     colnames(proteome_expression_data),
                     colnames(lab_expression_data),
                     colnames(cytokine_expression_data),
                     colnames(gutmicrobiome_expression_data)))


all_date <-
  all_date[all_date >= "2016-01-12" & all_date <= "2016-03-03"] %>% 
  stringr::str_sort(numeric = TRUE)

template <-
  matrix(1, nrow = 1, ncol = length(all_date)) %>% 
  as.data.frame()

colnames(template) <- all_date

exposomeChemical_expression_data1 <-
  exposomeChemical_expression_data %>%
  dplyr::select(one_of(all_date)) %>%
  dplyr::left_join(template, 
                   by = intersect(colnames(template), 
                                  colnames(exposomeChemical_expression_data))) %>% 
  `[`(,all_date,) 

exposomeBiological_expression_data1 <-
  exposomeBiological_expression_data %>%
  dplyr::select(one_of(all_date)) %>%
  dplyr::left_join(template, 
                   by = intersect(colnames(template), 
                                  colnames(exposomeBiological_expression_data))) %>% 
  `[`(,all_date,)


environment_expression_data1 <-
  environment_expression_data %>%
  dplyr::select(one_of(all_date)) %>%
  dplyr::left_join(template, 
                   by = intersect(colnames(template), 
                                  colnames(environment_expression_data))) %>% 
  `[`(,all_date,)

metabolome_expression_data1 <-
  metabolome_expression_data %>%
  dplyr::select(one_of(all_date)) %>%
  dplyr::left_join(template, 
                   by = intersect(colnames(template), 
                                  colnames(metabolome_expression_data))) %>% 
  `[`(,all_date,)

toxin_expression_data1 <-
  toxin_expression_data %>%
  dplyr::select(one_of(all_date)) %>%
  dplyr::left_join(template, 
                   by = intersect(colnames(template), 
                                  colnames(toxin_expression_data))) %>% 
  `[`(,all_date,)

proteome_expression_data1 <-
  proteome_expression_data %>%
  dplyr::select(one_of(all_date)) %>%
  dplyr::left_join(template, 
                   by = intersect(colnames(template), 
                                  colnames(proteome_expression_data))) %>% 
  `[`(,all_date,)

lab_expression_data1 <-
  lab_expression_data %>%
  dplyr::select(one_of(all_date)) %>%
  dplyr::left_join(template, 
                   by = intersect(colnames(template), 
                                  colnames(lab_expression_data))) %>% 
  `[`(,all_date,)

cytokine_expression_data1 <-
  cytokine_expression_data %>%
  dplyr::select(one_of(all_date)) %>%
  dplyr::left_join(template, 
                   by = intersect(colnames(template), 
                                  colnames(cytokine_expression_data))) %>% 
  `[`(,all_date,)

gutmicrobiome_expression_data1 <-
  gutmicrobiome_expression_data %>%
  dplyr::select(one_of(all_date)) %>%
  dplyr::left_join(template, 
                   by = intersect(colnames(template), 
                                  colnames(gutmicrobiome_expression_data))) %>% 
  `[`(,all_date,)


par(mar = c(2, 2, 2, 2))
par(xpd = TRUE)
library(circlize)
set.seed(999)

## 18 time points
df <- data.frame(
  factors = factor(
    colnames(exposomeChemical_expression_data1),
    levels = colnames(exposomeChemical_expression_data1)
  ),
  x = 1,
  y = 1
)

circos.par(
  "track.height" = 0.07,
  start.degree = 90,
  clock.wise = TRUE,
  gap.after = c(rep(0.5, 17), 90),
  cell.padding = c(0, 0, 0, 0
  )
)

circos.initialize(factors = df$factors,
                  x = df$x,
                  xlim = c(0.5,1.5))

# exposome chemical data
col_fun <- colorRamp2(c(-2, 0, 2),
                      c(ggsci::pal_aaas()(n = 10)[1], 
                        "white", 
                        ggsci::pal_aaas()(n =10)[2]))

temp_data <-
  exposomeChemical_expression_data1[1:50,] %>%
  # exposomeChemical_expression_data1 %>% 
  apply(1, function(x){
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  }) %>% 
  t() %>% 
  as.data.frame() %>% 
  as.matrix()

circos.track(
  factors = df$factors,
  # x = df$x,
  y = df$y,
  ylim = c(0, 1),
  bg.border = ggsci::pal_d3()(10)[2],
  # bg.col = NA,
  track.height = 0.1,
  panel.fun = function(x, y) {
    name = get.cell.meta.data("sector.index")
    i = get.cell.meta.data("sector.numeric.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    
    #text direction (dd) and adjusmtents (aa)
    theta = circlize(mean(xlim), 1.3)[1, 1] %% 360
    dd <- ifelse(theta < 90 || theta > 270, "clockwise", "reverse.clockwise")
    aa = c(1, 0.5)
    if(theta < 90 || theta > 270){
      aa = c(0, 0.5)
    }  
    #plot country labels
    circos.text(
      x = mean(xlim),
      y = 1.2,
      labels = name,
      facing = dd,
      cex = 0.6,
      adj = aa
    )
    
    nr <- nrow(temp_data)
    nc <- ncol(temp_data)
    
    col_mat <- col_fun(temp_data)
    col_mat[is.na(col_mat)] <- "white"
    
    circos.rect(xleft = xlim[1], 
                ybottom = (seq(1, nr, 1) - 1)/nr,
                xright = xlim[2],
                ytop = seq(1, nr, 1)/nr,
                col = col_mat[,i], 
                border = NA)
  }
)

# exposome biological
col_fun <- colorRamp2(c(-2, 0, 2),
                      c(ggsci::pal_aaas()(n = 10)[3], 
                        "white", 
                        ggsci::pal_aaas()(n =10)[6]))

dim(exposomeBiological_expression_data1)

temp_data <-
  exposomeBiological_expression_data1 %>%
  # exposomeBiological_expression_data1[1:10,] %>% 
  apply(1, function(x){
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  }) %>% 
  t() %>% 
  as.data.frame() %>% 
  as.matrix()

circos.track(
  factors = df$factors,
  # x = df$x,
  y = df$y,
  ylim = c(0, 1),
  bg.border = ggsci::pal_d3()(10)[5],
  # bg.col = NA,
  track.height = 0.05,
  panel.fun = function(x, y) {
    name = get.cell.meta.data("sector.index")
    i = get.cell.meta.data("sector.numeric.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    nr <- nrow(temp_data)
    nc <- ncol(temp_data)
    
    col_mat <- col_fun(temp_data)
    col_mat[is.na(col_mat)] <- "white"
    
    circos.rect(xleft = xlim[1], 
                ybottom = (seq(1, nr, 1) - 1)/nr,
                xright = xlim[2],
                ytop = seq(1, nr, 1)/nr,
                col = col_mat[,i], 
                border = NA)
  }
)

# environment
col_fun <- colorRamp2(c(-2, 0, 2),
                      c(ggsci::pal_d3()(n = 10)[1], 
                        "white", 
                        ggsci::pal_d3()(n =10)[2]))

dim(environment_expression_data1)

temp_data <-
  environment_expression_data1 %>%
  # exposomeBiological_expression_data1[1:10,] %>% 
  apply(1, function(x){
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  }) %>% 
  t() %>% 
  as.data.frame() %>% 
  as.matrix()

circos.track(
  factors = df$factors,
  # x = df$x,
  y = df$y,
  ylim = c(0, 1),
  bg.border = ggsci::pal_d3()(10)[1],
  # bg.col = NA,
  track.height = 0.05,
  panel.fun = function(x, y) {
    name = get.cell.meta.data("sector.index")
    i = get.cell.meta.data("sector.numeric.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    nr <- nrow(temp_data)
    nc <- ncol(temp_data)
    
    col_mat <- col_fun(temp_data)
    col_mat[is.na(col_mat)] <- "white"
    
    circos.rect(xleft = xlim[1], 
                ybottom = (seq(1, nr, 1) - 1)/nr,
                xright = xlim[2],
                ytop = seq(1, nr, 1)/nr,
                col = col_mat[,i], 
                border = NA)
  }
)

# metabolome
col_fun <- colorRamp2(c(-2, 0, 2),
                      c(ggsci::pal_d3()(n = 10)[10], 
                        "white", 
                        ggsci::pal_d3()(n =10)[9]))

dim(metabolome_expression_data1)

temp_data <-
  # metabolome_expression_data1 %>%
  metabolome_expression_data1[1:100,] %>%
  apply(1, function(x){
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  }) %>% 
  t() %>% 
  as.data.frame() %>% 
  as.matrix()

circos.track(
  factors = df$factors,
  # x = df$x,
  y = df$y,
  ylim = c(0, 1),
  bg.border = ggsci::pal_d3()(10)[3],
  # bg.col = NA,
  track.height = 0.1,
  panel.fun = function(x, y) {
    name = get.cell.meta.data("sector.index")
    i = get.cell.meta.data("sector.numeric.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    nr <- nrow(temp_data)
    nc <- ncol(temp_data)
    
    col_mat <- col_fun(temp_data)
    col_mat[is.na(col_mat)] <- "white"
    
    circos.rect(xleft = xlim[1], 
                ybottom = (seq(1, nr, 1) - 1)/nr,
                xright = xlim[2],
                ytop = seq(1, nr, 1)/nr,
                col = col_mat[,i], 
                border = NA)
  }
)

# toxin
col_fun <- colorRamp2(c(-2, 0, 2),
                      c(ggsci::pal_d3()(n = 10)[10], 
                        "white", 
                        ggsci::pal_d3()(n =10)[9]))

dim(toxin_expression_data1)

temp_data <-
  toxin_expression_data1 %>%
  # exposomeBiological_expression_data1[1:10,] %>% 
  apply(1, function(x){
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  }) %>% 
  t() %>% 
  as.data.frame() %>% 
  as.matrix()

circos.track(
  factors = df$factors,
  # x = df$x,
  y = df$y,
  ylim = c(0, 1),
  bg.border = ggsci::pal_d3()(10)[9],
  # bg.col = NA,
  track.height = 0.05,
  panel.fun = function(x, y) {
    name = get.cell.meta.data("sector.index")
    i = get.cell.meta.data("sector.numeric.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    nr <- nrow(temp_data)
    nc <- ncol(temp_data)
    
    col_mat <- col_fun(temp_data)
    col_mat[is.na(col_mat)] <- "white"
    
    circos.rect(xleft = xlim[1], 
                ybottom = (seq(1, nr, 1) - 1)/nr,
                xright = xlim[2],
                ytop = seq(1, nr, 1)/nr,
                col = col_mat[,i], 
                border = NA)
  }
)

# proteome
col_fun <- colorRamp2(c(-2, 0, 2),
                      c(ggsci::pal_d3()(n = 10)[1],
                        "white",
                        ggsci::pal_d3()(n = 10)[2]))

dim(proteome_expression_data1)

temp_data <-
  # proteome_expression_data1 %>%
  proteome_expression_data1[1:100,] %>%
  apply(1, function(x) {
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  }) %>%
  t() %>%
  as.data.frame() %>%
  as.matrix()

circos.track(
  factors = df$factors,
  # x = df$x,
  y = df$y,
  ylim = c(0, 1),
  bg.border = ggsci::pal_d3()(10)[4],
  # bg.col = NA,
  track.height = 0.1,
  panel.fun = function(x, y) {
    name = get.cell.meta.data("sector.index")
    i = get.cell.meta.data("sector.numeric.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    nr <- nrow(temp_data)
    nc <- ncol(temp_data)
    
    col_mat <- col_fun(temp_data)
    col_mat[is.na(col_mat)] <- "white"
    
    circos.rect(
      xleft = xlim[1],
      ybottom = (seq(1, nr, 1) - 1) / nr,
      xright = xlim[2],
      ytop = seq(1, nr, 1) / nr,
      col = col_mat[, i],
      border = NA
    )
  }
)


# blood test
col_fun <- colorRamp2(c(-2, 0, 2),
                      c(ggsci::pal_futurama()(n = 10)[9],
                        "white",
                        ggsci::pal_futurama()(n = 10)[2]))

dim(lab_expression_data1)

temp_data <-
  lab_expression_data1 %>%
  # exposomeBiological_expression_data1[1:10,] %>%
  apply(1, function(x) {
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  }) %>%
  t() %>%
  as.data.frame() %>%
  as.matrix()

circos.track(
  factors = df$factors,
  # x = df$x,
  y = df$y,
  ylim = c(0, 1),
  bg.border = ggsci::pal_d3()(10)[7],
  # bg.col = NA,
  track.height = 0.05,
  panel.fun = function(x, y) {
    name = get.cell.meta.data("sector.index")
    i = get.cell.meta.data("sector.numeric.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    nr <- nrow(temp_data)
    nc <- ncol(temp_data)
    
    col_mat <- col_fun(temp_data)
    col_mat[is.na(col_mat)] <- "white"
    
    circos.rect(
      xleft = xlim[1],
      ybottom = (seq(1, nr, 1) - 1) / nr,
      xright = xlim[2],
      ytop = seq(1, nr, 1) / nr,
      col = col_mat[, i],
      border = NA
    )
  }
)

# cytokine
col_fun <- colorRamp2(c(-2, 0, 2),
                      c(ggsci::pal_d3()(n = 10)[1],
                        "white",
                        ggsci::pal_d3()(n = 10)[2]))

dim(cytokine_expression_data1)

temp_data <-
  cytokine_expression_data1 %>%
  # exposomeBiological_expression_data1[1:10,] %>%
  apply(1, function(x) {
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  }) %>%
  t() %>%
  as.data.frame() %>%
  as.matrix()

circos.track(
  factors = df$factors,
  # x = df$x,
  y = df$y,
  ylim = c(0, 1),
  bg.border = ggsci::pal_d3()(10)[8],
  # bg.col = NA,
  track.height = 0.05,
  panel.fun = function(x, y) {
    name = get.cell.meta.data("sector.index")
    i = get.cell.meta.data("sector.numeric.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    nr <- nrow(temp_data)
    nc <- ncol(temp_data)
    
    col_mat <- col_fun(temp_data)
    col_mat[is.na(col_mat)] <- "white"
    
    circos.rect(
      xleft = xlim[1],
      ybottom = (seq(1, nr, 1) - 1) / nr,
      xright = xlim[2],
      ytop = seq(1, nr, 1) / nr,
      col = col_mat[, i],
      border = NA
    )
  }
)

# gutmicrobiome
col_fun <- colorRamp2(c(-2, 0, 2),
                      c(ggsci::pal_futurama()(n = 10)[3],
                        "white",
                        ggsci::pal_futurama()(n = 10)[2]))

dim(gutmicrobiome_expression_data1)

temp_data <-
  gutmicrobiome_expression_data1 %>%
  # exposomeBiological_expression_data1[1:10,] %>%
  apply(1, function(x) {
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  }) %>%
  t() %>%
  as.data.frame() %>%
  as.matrix()

circos.track(
  factors = df$factors,
  # x = df$x,
  y = df$y,
  ylim = c(0, 1),
  bg.border = ggsci::pal_d3()(10)[6],
  # bg.col = NA,
  track.height = 0.05,
  panel.fun = function(x, y) {
    name = get.cell.meta.data("sector.index")
    i = get.cell.meta.data("sector.numeric.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    nr <- nrow(temp_data)
    nc <- ncol(temp_data)
    
    col_mat <- col_fun(temp_data)
    col_mat[is.na(col_mat)] <- "white"
    
    circos.rect(
      xleft = xlim[1],
      ybottom = (seq(1, nr, 1) - 1) / nr,
      xright = xlim[2],
      ytop = seq(1, nr, 1) / nr,
      col = col_mat[, i],
      border = NA
    )
  }
)


