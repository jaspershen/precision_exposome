###avoid source
no_function()


##load data
sxtTools::setwd_project()
library(tidyverse)
rm(list = ls())

source("R/20200511/piumet_code.R")

load("data_20200511/exposome/variable_info")
exp_variable_info <- variable_info

load("data_20200511/metabolome/variable_info")
met_variable_info <- variable_info

met_annotation_info <- readr::read_csv("data_20200511/metabolome/annotation_table.csv")

setwd("data_analysis/met_exp/")

load("met_expression_data")
load("exp_expression_data")
load("met_sample_info")
load("exp_sample_info")

met_sample_info$sample_id == colnames(met_expression_data)
colnames(met_expression_data) <- as.character(met_sample_info$CollectionDate)
colnames(exp_expression_data) <- colnames(met_expression_data)


#-----------------------------------------------------------------------------------
#####module analysis
##read results from PIUMet
##get the cluster 1 peak table
cluster <- readr::read_csv("exposome_cluster/cluster.csv")
load("cor_value")

cluster <- 
cluster %>% 
  dplyr::left_join(cor_value, by = c("peak" = "from")) %>% 
  dplyr::select(peak = to, cluster) %>% 
  dplyr::left_join(met_variable_info, by = c("peak" = "peak_name")) %>% 
  dplyr::select(peak, cluster, Mass) %>% 
  dplyr::mutate(
    polarity = case_when(
      stringr::str_detect(peak, "p") ~ "positive",
      stringr::str_detect(peak, "n") ~ "negative"
    )
  )

colnames(cluster)[1] <- "name"
colnames(cluster)[3] <- "mz"

cluster$mz <- 
  cluster$mz %>% 
  stringr::str_replace("\\_[0-9\\.]{1,6}", "") %>% 
  as.numeric()

cluster1 <- 
  cluster %>% 
  dplyr::filter(cluster == 1)

cluster2 <- 
  cluster %>% 
  dplyr::filter(cluster == 2)

plot1 <- 
readPIUMet(path = "exposome_cluster/piumet_output_cluster1", 
           marker_table = cluster1,
           size_range = c(0.5, 3),
           width_range = c(0.2, 0.5),
           text = FALSE)

plot1

plot2 <- 
  readPIUMet(path = "exposome_cluster/piumet_output_cluster2", 
             marker_table = cluster2,
             size_range = c(0.5, 3),
             width_range = c(0.2, 0.5),
             text = FALSE)

plot2

##pathway analysis for two cluster
load("hsa_pathway")
load("exposome_cluster/piumet_output_cluster1/Result/node_data")

hmdb_id <- node_data$HMDB_ID
hmdb_id <- hmdb_id[hmdb_id != " " & hmdb_id != "NULL"]


all_module_pathway <- 
  enrich_pathway(id = module_node$node,
                 pathway_database = hsa_pathway)

save(all_module_pathway, file = "all_module_pathway")


all_module_pathway$Pathway.name <- 
  all_module_pathway$Pathway.name %>% 
  stringr::str_replace_all(pattern = " - Homo sapiens \\(human\\)", replacement = "")

plot <- 
  all_module_pathway %>% 
  `[`(1:20, ,) %>% 
  dplyr::arrange(desc(p.value.fdr)) %>% 
  dplyr::mutate(Pathway.name = factor(Pathway.name, levels = Pathway.name)) %>% 
  # dplyr::filter(Overlap.frac > 0.25) %>% 
  dplyr::mutate(class = case_when(
    p.value.fdr < 0.01 & Overlap >= 10 ~ "Yes",
    TRUE ~ "No"
  )) %>% 
  ggplot(aes(x = Pathway.name, -log(p.value.fdr, 10))) +
  geom_segment(aes(x = Pathway.name, y = 0, 
                   xend = Pathway.name, 
                   yend = -log(p.value.fdr, 10), color = class), show.legend = FALSE) +
  geom_point(aes(size = Overlap, color = class), show.legend = TRUE) +
  scale_color_manual(values = c("Yes" = "#FB8072", "No" = "#D9D9D9")) +
  scale_size_continuous(range = c(3,8)) +
  theme_classic() +
  labs(y = "-log10(P-value, FDR)", x = "") +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12), 
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.position = c(1,0), legend.justification = c(1,0),
        legend.background = element_blank(),
        strip.background = element_rect(fill = "#0099B47F"),
        strip.text = element_text(color = "white", size = 13)) +
  ggplot2::coord_flip()

plot
ggsave(plot, filename = "all_module_pathway.png",
       height = 7, width = 8, bg = "transparent")

ggsave(plot, filename = "all_module_pathway.pdf",
       height = 7, width = 8, bg = "transparent")



###annotation
annotation_from_module <- 
  rbind(
    module_edge %>% 
      dplyr::filter(class == "p-c")
    # module2_edge %>% 
    #   dplyr::filter(class == "p-c")
    # module3_edge %>% 
    #   dplyr::filter(class == "p-c")
  ) %>% 
  dplyr::distinct(from, to, .keep_all = TRUE) %>% 
  dplyr::arrange(from)


load("kegg.compound.rda")


# marker_info %>% 
#   dplyr::filter(!is.na(Compound.name)) %>% 
#   dplyr::select(name, Compound.name, HMDB.ID, KEGG.ID) %>% 
#   dplyr::left_join(annotation_from_module, by = c("name" = "from"))

library(plyr)
annotation_from_module <- 
  annotation_from_module %>% 
  plyr::dlply(.(from))


annotation_from_module <- 
  annotation_from_module %>% 
  lapply(function(x) {
    data.frame(
      name = x[1,1],
      annotation = apply(x, 1, function(y){paste(y[-c(1,5)], collapse = ";")}) 
      %>% paste(collapse = "{}"),
      stringsAsFactors = FALSE
    ) 
  })

annotation_from_module <- 
  do.call(rbind, annotation_from_module)

annotation_from_module <- 
  as.data.frame(annotation_from_module)








##get the edge of module 2
idx <- 13
v <- 
  c(module[[idx]], 
    data.frame(ID = module[[idx]]) %>% 
      dplyr::left_join(annotation_result[,1:2], by = c("ID")) %>% 
      pull(peak) %>% 
      unique())

v <- v[!is.na(v)]
module13_network <- 
  igraph::subgraph(graph = module_network, 
                   v = v)

save(module13_network, file = 'module13_network')

plot <- 
  ggraph::ggraph(graph = module13_network, layout = "kk") +
  ggraph::geom_edge_link(aes(color = class)) +
  ggraph::scale_edge_color_manual(values = c("c-c" = "black", "p-c" = "grey")) +
  ggraph::geom_node_point(aes(color = class, shape = class, size = importance)) +
  scale_color_manual(values = c("Detected metabolite" = "#FB8072",
                                "Hidden metabolite" = "#D9D9D9",
                                "Peak" = "#80B1D3")) +
  scale_size_continuous(range = c(3,8)) +
  theme_void()  +
  theme(plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA)
  )

plot

ggsave(plot, filename = "module13_netwprk.png",
       width = 9, height = 7, 
       bg = "transparent")

ggsave(plot, filename = "module13_netwprk.pdf",
       width = 9, height = 7, 
       bg = "transparent")




##pathway analysis for module3
load("hsa_pathway")

module3_pathway <- 
  enrich_pathway(id = grep("C[0-9]{1,10}",
                           names(igraph::V(module3_network)), value = TRUE),
                 pathway_database = hsa_pathway)

save(module3_pathway, file = "module3_pathway")


module3_pathway$Pathway.name <- 
  module3_pathway$Pathway.name %>% 
  stringr::str_replace_all(pattern = " - Homo sapiens \\(human\\)", replacement = "")

plot <- 
  module3_pathway %>% 
  `[`(1:15, ,) %>% 
  dplyr::arrange(desc(p.value.fdr)) %>% 
  dplyr::mutate(Pathway.name = factor(Pathway.name, levels = Pathway.name)) %>% 
  # dplyr::filter(Overlap.frac > 0.25) %>% 
  dplyr::mutate(class = case_when(
    p.value.fdr < 0.01 & Overlap >= 10 ~ "Yes",
    TRUE ~ "No"
  )) %>% 
  ggplot(aes(x = Pathway.name, -log(p.value.fdr, 10))) +
  geom_segment(aes(x = Pathway.name, y = 0, 
                   xend = Pathway.name, 
                   yend = -log(p.value.fdr, 10), color = class), show.legend = FALSE) +
  geom_point(aes(size = Overlap, color = class), show.legend = TRUE) +
  scale_color_manual(values = c("Yes" = "#FB8072", "No" = "#D9D9D9")) +
  scale_size_continuous(range = c(3,8)) +
  theme_classic() +
  labs(y = "-log10(P-value, FDR)", x = "") +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12), 
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.position = c(1,0), legend.justification = c(1,0),
        legend.background = element_blank(),
        strip.background = element_rect(fill = "#0099B47F"),
        strip.text = element_text(color = "white", size = 13)) +
  ggplot2::coord_flip()

plot
ggsave(plot, filename = "module3_pathway.png",
       height = 7, width = 8, bg = "transparent")

ggsave(plot, filename = "module3_pathway.pdf",
       height = 7, width = 8, bg = "transparent")




####annotation check
# check_table <- 
# marker_info_new %>% 
#   dplyr::filter(!is.na(Compound.name), !is.na(annotation))
# 
# write.csv(check_table, file = "check_table.csv", row.names = FALSE)
