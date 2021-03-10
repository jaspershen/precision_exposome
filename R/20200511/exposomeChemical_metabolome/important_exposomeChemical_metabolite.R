##avoid source
no_exist()

sxtTools::setwd_project()
source("R/20200511/tools.R")
library(tidyverse)
rm(list = ls())

####load data
load("data_20200511/exposome/variable_info")
exp_variable_info <- variable_info

load("data_20200511/metabolome/variable_info")
met_variable_info <- variable_info

met_annotation_info <- readr::read_csv("data_20200511/metabolome/annotation_table.csv")

setwd("data_analysis/expsomeChemical_metabolome/")

load("met_expression_data")
load("exp_expression_data")
load("met_sample_info")
load("exp_sample_info")

met_sample_info$sample_id == colnames(met_expression_data)
colnames(met_expression_data) <- as.character(met_sample_info$CollectionDate)
colnames(exp_expression_data) <- colnames(met_expression_data)

#######correlation data
exp_expression_data <- log(exp_expression_data + 1, 2)
met_expression_data <- log(met_expression_data + 1, 2)

load('cor_value')

cor_value1 <- 
  cor_value %>% 
  dplyr::filter(abs(cor) > 0.9 & p_value < 0.05)

dim(cor_value)
dim(cor_value1)

cor_value1 <- 
  cor_value1 %>% 
  dplyr::left_join(exp_variable_info[,c(1:2)], by = c("from" = "peak_ID")) %>% 
  dplyr::rename(exp_id = MetabID) %>% 
  dplyr::left_join(met_annotation_info[,c(1,5)], by = c("to" = "name")) %>% 
  dplyr::rename(met_id = Compound.name) %>% 
  dplyr::left_join(met_variable_info[,c(1,3)], by = c("to" = "peak_name")) %>% 
  dplyr::rename(met_id2 = Metabolite) %>% 
  dplyr::filter(!is.na(met_id))


dir.create("important_exposomeChemical_metabolite", showWarnings = FALSE)
setwd("important_exposomeChemical_metabolite/")

###how many exp for each met
met_num <- 
  cor_value1 %>% 
  dplyr::group_by(to) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::arrange(dplyr::desc(n)) %>% 
  dplyr::rename(exp_number_for_each_met = n)

##how many met for each exp
exp_num <- 
  cor_value1 %>% 
  dplyr::group_by(from) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::arrange(dplyr::desc(n)) %>% 
  dplyr::rename(met_number_for_each_exp = n)

cor_value1 <-
  cor_value1 %>% 
  dplyr::left_join(met_num, by = c("to")) %>% 
  dplyr::left_join(exp_num, by = c("from")) 

sxtTools::setwd_project()
setwd("data_analysis/expsomeChemical_metabolome/important_exposomeChemical_metabolite/")

exp_met_cor <- 
  cor_value1 %>% 
  dplyr::arrange(from)

write.csv(exp_met_cor, "met_exp_cor.csv", row.names = FALSE)

##
import_met <- readxl::read_xlsx("../cor_value_peng.xlsx", sheet = 2)

library(igraph)
library(ggraph)
library(tidygraph)

##pheatmap
library(pheatmap)

temp_data <- 
  met_expression_data[import_met$met_name,]

colnames(temp_data) <- c(1:5)

rownames(temp_data) <-
  import_met$met_id[match(rownames(temp_data),
                          import_met$met_name)]


temp_data <- 
  exp_expression_data[unique(cor_value1$exp_name[which(cor_value1$met_name %in% import_met$met_name)]),]

colnames(temp_data) <- c(1:5)

rownames(temp_data) <-
  exp_variable_info$MetabID[match(rownames(temp_data),
                                  exp_variable_info$peak_ID)]
 
annotation_row <-
  data.frame(name = rownames(temp_data), stringsAsFactors = FALSE) %>%
  dplyr::left_join(exp_variable_info[, c(2, 3)], by = c("name" = "MetabID")) %>% 
  tibble::column_to_rownames(var = "name")


##remove unknown
remove_idx <- which(annotation_row$Chemical_class == "Unknown")
temp_data <- temp_data[-remove_idx,,drop = FALSE]
annotation_row <- annotation_row[-remove_idx,,drop = FALSE]

chemical_class <- 
  annotation_row$Chemical_class %>% unique() %>% 
  sort()

library(ggsci)

col <-
  c(pal_futurama()(12), pal_jama()(3))

names(col) <- chemical_class

annotation_colors <- 
  list(Chemical_class = col)

plot <- 
  pheatmap(temp_data, scale = "row", 
           show_rownames = FALSE, 
           border_color = "white",
           annotation_row = annotation_row, 
           annotation_colors = annotation_colors,
           clustering_method = "ward.D")

ggsave(plot, filename = "important_exposome_heatmap.pdf", width = 10, height = 7)




###heatmap for all exposome
temp_data <- 
  exp_expression_data

colnames(temp_data) <- c(1:5)

rownames(temp_data) <-
  exp_variable_info$MetabID

annotation_row <-
  data.frame(name = rownames(temp_data), stringsAsFactors = FALSE) %>%
  dplyr::left_join(exp_variable_info[, c(2, 3)], by = c("name" = "MetabID")) %>% 
  tibble::column_to_rownames(var = "name")


##remove unknown
remove_idx <- which(annotation_row$Chemical_class == "Unknown")
temp_data <- temp_data[-remove_idx,,drop = FALSE]
annotation_row <- annotation_row[-remove_idx,,drop = FALSE]

chemical_class <- 
  annotation_row$Chemical_class %>% unique() %>% 
  sort()

library(ggsci)

col <-
  c(pal_futurama()(12), pal_jama()(3))

names(col) <- chemical_class

annotation_colors <- 
  list(Chemical_class = col)

plot <- 
  pheatmap(temp_data, scale = "row", 
           show_rownames = FALSE, 
           border_color = "white",
           annotation_row = annotation_row, 
           annotation_colors = annotation_colors,
           clustering_method = "ward.D")

ggsave(plot, filename = "all_exposome_heatmap.png", width = 10, height = 7)

ggsave(plot, filename = "all_important_exposome_heatmap.pdf", width = 10, height = 7)


###network for important metabolites
idx <- 20

edge_data <-  
  cor_value1 %>% 
  dplyr::filter(met_name %in% import_met$met_name[idx]) %>%
  dplyr::rename(from = exp_name, 
                to = met_id, 
                Correlation = estimate) %>% 
  dplyr::mutate(p.value = -log(p.value, 10))

node_data <- 
  cor_value1 %>% 
  dplyr::filter(met_name %in% import_met$met_name[idx]) %>%
  dplyr::rename(from = exp_name, to = met_id) %>% 
  dplyr::select(from, to) %>% 
  tidyr::pivot_longer(cols = c(from, to), 
                      names_to = "class", values_to = "node") %>% 
  dplyr::mutate(class1 = case_when(
    stringr::str_detect(class, "from") ~ "Exposome",
    TRUE ~ "Metabolome"
  )) %>% 
  dplyr::select(node, class1) %>% 
  dplyr::rename(Class = class1) %>%
  dplyr::distinct(node, .keep_all = TRUE)

temp_data <- 
  tidygraph::tbl_graph(nodes = node_data, 
                       edges = edge_data,
                       directed = TRUE) %>% 
  dplyr::mutate(Degree = centrality_degree(mode = 'all'))

plot <- 
  ggraph(temp_data,
         layout = 'linear',
         circular = TRUE) +
  geom_edge_arc(
    aes(width = p.value, 
        color = Correlation),
    show.legend = TRUE,
    arrow = arrow(length = unit(2, 'mm'))
  ) +
  geom_node_text(aes(label = ifelse(Class == "Metabolome", node, NA))) +
  geom_node_point(aes(color = Class, size = Degree), show.legend = TRUE) +
  ggraph::scale_edge_color_gradientn(colours = c(alpha("#8DD3C7", 1),
                                                 "white",
                                                 alpha("#FB8072", 1))) +
  ggraph::scale_edge_width(range = c(0.2,2)) +
  scale_size_continuous(range = c(1, 4)) +
  ggsci::scale_color_aaas() +
  theme_void() +
  theme(plot.background = element_rect(fill = "transparent", color = NA), 
        panel.background = element_rect(fill = "transparent", color = NA))

plot

ggsave(plot, filename = paste("important_met_network",idx,".pdf",sep=""), 
       width = 8, height = 7, bg = "transparent")

ggsave(plot, filename = paste("important_met_network",idx,".png",sep=""),
       width = 8, height = 7, , bg = "transparent")

import_met <- 
  apply(import_met, 1, function(x){
    x <- as.character(x)
    cor <- dplyr::filter(cor_value1, met_name == x[1]) %>% 
      dplyr::pull(estimate)
    c(sum(cor > 0), sum(cor < 0)) * 100/length(cor)
  }) %>% 
  t() %>% 
  cbind(import_met, .)


colnames(import_met)[c(12,13)] <- c("pos", "neg")

import_met <- 
import_met[,c(1,2,3,12,13)] %>% 
  dplyr::mutate(class = case_when(
    pos > 60 ~ "positive dominant",
    TRUE ~ "negative dominant"
  ))

write.csv(import_met, "import_met.csv")


important_exp <- 
  cor_value1 %>% 
  dplyr::group_by(exp_name) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::arrange(desc(n)) %>% 
  dplyr::filter(n >= 20)


important_exp <- 
  apply(important_exp, 1, function(x){
    x <- as.character(x)
    cor <- dplyr::filter(cor_value1, exp_name == x[1]) %>% 
      dplyr::pull(estimate)
    c(sum(cor > 0), sum(cor < 0)) * 100/length(cor)
  }) %>% 
  t() %>% 
  cbind(important_exp, .)


colnames(important_exp)[c(3,4)] <- c("positive", "negative")

write.csv(important_exp, "important_exp.csv", row.names = FALSE)

###correlation network for important exposome
idx <- 5

edge_data <-  
  cor_value1 %>% 
  dplyr::filter(exp_name %in% important_exp$exp_name[idx]) %>%
  dplyr::rename(from = exp_id, 
                to = met_id, 
                Correlation = estimate) %>% 
  dplyr::mutate(p.value = -log(p.value, 10))

node_data <- 
  cor_value1 %>% 
  dplyr::filter(exp_name %in% important_exp$exp_name[idx]) %>%
  dplyr::rename(from = exp_id, to = met_id) %>% 
  dplyr::select(from, to) %>% 
  tidyr::pivot_longer(cols = c(from, to), 
                      names_to = "class", values_to = "node") %>% 
  dplyr::mutate(class1 = case_when(
    stringr::str_detect(class, "from") ~ "Exposome",
    TRUE ~ "Metabolome"
  )) %>% 
  dplyr::select(node, class1) %>% 
  dplyr::rename(Class = class1) %>%
  dplyr::distinct(node, .keep_all = TRUE)

temp_data <- 
  tidygraph::tbl_graph(nodes = node_data, 
                       edges = edge_data,
                       directed = TRUE) %>% 
  dplyr::mutate(Degree = centrality_degree(mode = 'all'))

plot <- 
  ggraph(temp_data,
         layout = 'linear',
         circular = TRUE) +
  geom_edge_arc(
    aes(width = p.value, 
        color = Correlation),
    show.legend = TRUE,
    arrow = arrow(length = unit(2, 'mm'))
  ) +
  # geom_node_text(aes(label = ifelse(Class == "Exposome", node, NA))) +
  geom_node_text(aes(label = node)) +
  geom_node_point(aes(color = Class, size = Degree), show.legend = TRUE) +
  ggraph::scale_edge_color_gradientn(colours = c(alpha("#8DD3C7", 1),
                                                 "white",
                                                 alpha("#FB8072", 1))) +
  ggraph::scale_edge_width(range = c(0.2,2)) +
  scale_size_continuous(range = c(1, 4)) +
  # scale_color_manual(values = c("Exposome" = "#3B4992FF",
  #                               "Proteome" = "#008280FF")) +
  ggsci::scale_color_aaas() +
  theme_void() +
  theme(plot.background = element_rect(fill = "transparent", color = NA), 
        panel.background = element_rect(fill = "transparent", color = NA))

plot


ggsave(plot, filename = paste("important_exp_network",idx,".pdf",sep=""), 
       width = 8, height = 7, bg = "transparent")

ggsave(plot, filename = paste("important_exp_network",idx,".png",sep=""),
       width = 8, height = 7, , bg = "transparent")




###get the network for importance metabolites and its exposome
edge_data <-  
  cor_value1 %>% 
  dplyr::filter(met_name %in% import_met$met_name) %>%
  dplyr::rename(from = exp_id ,
                to = met_id, 
                Correlation = estimate) %>% 
  dplyr::mutate(p.value = -log(p.value, 10))

node_data <- 
  cor_value1 %>% 
  dplyr::filter(met_name %in% import_met$met_name) %>%
  dplyr::rename(from = exp_id, to = met_id) %>% 
  dplyr::select(from, to) %>% 
  tidyr::pivot_longer(cols = c(from, to), 
                      names_to = "class", values_to = "node") %>% 
  dplyr::mutate(class1 = case_when(
    stringr::str_detect(class, "from") ~ "Exposome",
    TRUE ~ "Metabolome"
  )) %>% 
  dplyr::select(node, class1) %>% 
  dplyr::rename(Class = class1) %>%
  dplyr::distinct(node, .keep_all = TRUE) %>% 
  dplyr::left_join(cor_value1[,c(2,7)] %>% dplyr::distinct(met_name, .keep_all = TRUE),
                   by = c("node" = "met_id")) %>% 
  dplyr::left_join(import_met[,c(1,6)], by = c("met_name")) %>% 
  dplyr::select(-met_name)

node_data$class[is.na(node_data$class)] <- "Exposome"

temp_data <- 
  tidygraph::tbl_graph(nodes = node_data, 
                       edges = edge_data,
                       directed = TRUE) %>% 
  dplyr::mutate(Degree = centrality_degree(mode = 'all'))

plot <- 
  ggraph(temp_data,
         layout = 'linear',
         circular = TRUE) +
  geom_edge_arc(
    aes(width = p.value, 
        color = Correlation),
    show.legend = TRUE,
    arrow = arrow(length = unit(2, 'mm'))
  ) +
  geom_node_text(aes(label = ifelse(Class == "Metabolome", node, NA)),
                 check_overlap = FALSE) +
  geom_node_point(aes(color = class, size = Degree), show.legend = TRUE) +
  ggraph::scale_edge_color_gradientn(colours = c(alpha("#8DD3C7", 1),
                                                 "white",
                                                 alpha("#FB8072", 1))) +
  ggraph::scale_edge_width(range = c(0.2,2)) +
  scale_size_continuous(range = c(1, 4)) +
  scale_color_manual(values = c(
    "Exposome" = ggsci::pal_aaas()(1),
    "positive dominant" = ggsci::pal_aaas(alpha = 0.7)(12)[6],
    "negative dominant" = ggsci::pal_aaas(alpha = 0.7)(12)[5]
  )) +
  # ggsci::scale_color_aaas() +
  theme_void() +
  theme(plot.background = element_rect(fill = "transparent", color = NA), 
        panel.background = element_rect(fill = "transparent", color = NA))

plot


ggsave(plot, filename = "important_metabolite_exposome.pdf", 
       width = 8, height = 7, bg = "transparent")

ggsave(plot, filename ="important_metabolite_exposome.png", 
       width = 8, height = 7, , bg = "transparent")
