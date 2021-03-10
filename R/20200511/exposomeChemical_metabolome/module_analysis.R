##load data
sxtTools::setwd_project()
library(tidyverse)
rm(list = ls())

source("R/20200511/exp_met/module_tools.R")

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
### module analyis
load("kegg.adduct.pos.rda")
load("kegg.adduct.neg.rda")
load("kegg_network")

##remove the unnormal adduct of metabolites
kegg.adduct.pos <- 
  kegg.adduct.pos %>% 
  dplyr::filter(Column.rp == TRUE) %>% 
  dplyr::filter(Adduct %in% c("M+H", "M+H-H2O", "M+H-2H2O", "M+NH4", "M+Na", "2M+H"))

kegg.adduct.neg <- 
  kegg.adduct.neg %>% 
  dplyr::filter(Column.rp == TRUE) %>% 
  dplyr::filter(Adduct %in% c("M-H", "M+Cl", "M-H-H2O", "M+CH3COO"))


###########################################################################
##identify initial modules
##prepare data for module analysis
load("cor_value")

cor_value1 <- 
  cor_value %>% 
  dplyr::filter(abs(cor) > 0.9 & fdr < 0.05)


# cor_value1 <- 
#   cor_value1 %>% 
#   dplyr::left_join(exp_variable_info[,c(1:2)], by = c("from" = "peak_ID")) %>% 
#   dplyr::rename(exp_id = MetabID) %>% 
#   dplyr::left_join(met_annotation_info[,c(1,5)], by = c("to" = "name")) %>% 
#   dplyr::rename(met_id = Compound.name) %>% 
#   dplyr::left_join(met_variable_info[,c(1,3)], by = c("to" = "peak_name")) %>% 
#   dplyr::rename(met_id2 = Metabolite) %>% 
#   dplyr::filter(!is.na(met_id))

marker_info <- 
  cor_value1 %>% 
  # dplyr::filter(exp_name %in% cluster2) %>%
  dplyr::select(to) %>% 
  dplyr::distinct(to) %>% 
  dplyr::left_join(met_variable_info, by = c("to" = "peak_name")) %>% 
  dplyr::select(to, Mass) %>% 
  dplyr::rename(name = to) %>% 
  dplyr::mutate(mz = stringr::str_split(Mass, "_") %>% 
                  lapply(function(x) x[1]) %>% 
                  unlist() %>% 
                  as.numeric(),
                rt = stringr::str_split(Mass, "_") %>% 
                  lapply(function(x) x[2]) %>% 
                  unlist() %>% 
                  as.numeric()) %>% 
  dplyr::select(name, mz, rt) %>% 
  dplyr::mutate(rt = rt * 60)

mz.pos <-
  marker_info %>% dplyr::filter(stringr::str_detect(name, "p")) %>% pull(mz)

rt.pos <-
  marker_info %>% dplyr::filter(stringr::str_detect(name, "p")) %>% pull(rt)

mz.neg <-
  marker_info %>% dplyr::filter(stringr::str_detect(name, "n")) %>% pull(mz)

rt.neg <-
  marker_info %>% dplyr::filter(stringr::str_detect(name, "n")) %>% pull(mz)

annotation_result_pos <-
  keggMatch(
    mz = mz.pos,
    rt = rt.pos,
    metabolite = kegg.adduct.pos,
    mz.tol = 25,
    rt.tol = 1000000,
    polarity = "positive"
  )

dir.create("community_analysis")

save(annotation_result_pos, file = "community_analysis/annotation_result_pos")

load("community_analysis/annotation_result_pos")

annotation_result_neg <-
  keggMatch(
    mz = mz.neg,
    rt = rt.neg,
    metabolite = kegg.adduct.neg,
    mz.tol = 25,
    rt.tol = 1000000,
    polarity = "negative"
  )

save(annotation_result_neg, file = "community_analysis/annotation_result_neg")

load("community_analysis/annotation_result_neg")

metabolite.id <-
  c(
    lapply(annotation_result_pos, function(x){
      if(is.null(x)) return(NULL)
      x$ID
    }) %>%
      unlist() %>%
      unique(),
    lapply(annotation_result_neg, function(x){
      if(is.null(x)) return(NULL)
      x$ID
    }) %>%
      unlist() %>%
      unique()
  ) %>%
  unique()
# 
# #####get the initial modules
library(igraph)
initial_modules <-
  getModule(node.id = metabolite.id,
            metabolic.network = kegg_network,
            threads = 3,
            max.reaction = 3,
            progressbar = TRUE,
            calculate.impact = TRUE)


save(initial_modules, file = "community_analysis/initial_modules")
load("community_analysis/initial_modules") 

# module <- initial_modules$module[module.p < 0.05]
module <- initial_modules$module

save(module, file = "community_analysis/module")
load("community_analysis/module")

##nodule is the significant modules
library(ggraph)

##module 1
module1 <-
  data.frame(id = module[[1]], stringsAsFactors = FALSE) %>%
  mutate(detected = case_when(id %in% initial_modules$hidden.id ~ "No",
                              TRUE ~ 'Yes',))

module2 <-
  data.frame(id = module[[2]], stringsAsFactors = FALSE) %>%
  mutate(detected = case_when(id %in% initial_modules$hidden.id ~ "No",
                              TRUE ~ 'Yes',))

module3 <-
  data.frame(id = module[[3]], stringsAsFactors = FALSE) %>%
  mutate(detected = case_when(id %in% initial_modules$hidden.id ~ "No",
                              TRUE ~ 'Yes',))

module4 <-
  data.frame(id = module[[4]], stringsAsFactors = FALSE) %>%
  mutate(detected = case_when(id %in% initial_modules$hidden.id ~ "No",
                              TRUE ~ 'Yes',))

all_module <-
  data.frame(id = unlist(module) %>% unique(), 
             stringsAsFactors = FALSE) %>%
  mutate(detected = case_when(id %in% initial_modules$hidden.id ~ "No",
                              TRUE ~ 'Yes',))

##construct network according to modules
names(annotation_result_pos) <- 
  marker_info %>% 
  dplyr::filter(stringr::str_detect(name, "p")) %>% 
  pull(name)

names(annotation_result_neg) <- 
  marker_info %>% 
  dplyr::filter(stringr::str_detect(name, "n")) %>% 
  pull(name)

annotation_result_pos <- 
  purrr::map2(
    .x = annotation_result_pos,
    .y = names(annotation_result_pos),
    .f = function(x, y) {
      if(is.null(x)) return(NULL)
      data.frame(peak = y, x[,c(1,6,7,12)])
    }
  )

annotation_result_pos <- 
  annotation_result_pos %>% 
  dplyr::bind_rows() %>% 
  distinct()


annotation_result_pos <- 
  annotation_result_pos %>% 
  dplyr::filter(ID %in% unique(unlist(module)))


annotation_result_neg <- 
  purrr::map2(
    .x = annotation_result_neg,
    .y = names(annotation_result_neg),
    .f = function(x, y) {
      if(is.null(x)) return(NULL)
      data.frame(peak = y, x[,c(1,6,7,12)])
    }
  )

annotation_result_neg <- 
  annotation_result_neg %>% 
  dplyr::bind_rows() %>% 
  distinct()


annotation_result_neg <- 
  annotation_result_neg %>% 
  dplyr::filter(ID %in% unique(unlist(module)))


annotation_result <- rbind(annotation_result_pos, annotation_result_neg)





##get the edge of all module
module_edge <- 
  igraph::subgraph(graph = kegg_network, 
                   v = all_module$id) %>% 
  igraph::as_data_frame()

module_edge <- 
  data.frame(module_edge, mz.error = 0, Adduct = NA, class = "c-c") %>% 
  rbind(.,
        annotation_result %>% 
          dplyr::filter(ID %in% all_module$id) %>% 
          dplyr::select(from = peak, to = ID, 
                        mz.error = mz.error,
                        Adduct) %>% 
          data.frame(., class = "p-c")
  ) 

module_node <- 
  unique(c(module_edge$from, module_edge$to)) %>% 
  data.frame(node = ., stringsAsFactors = FALSE) %>% 
  dplyr::left_join(all_module, by = c("node" = "id")) %>% 
  dplyr::mutate(class = case_when(
    is.na(detected) ~ "Peak",
    detected == "Yes" ~ "Detected metabolite",
    detected == "No" ~ "Hidden metabolite"
    # TRUE ~ NA
  ))


module_network <- 
  igraph::graph_from_data_frame(d = module_edge, 
                                directed = FALSE, 
                                vertices = module_node)

degree <- 
  igraph::degree(graph = module_network, mode = "total", normalized = TRUE)

betweenness <- igraph::betweenness(graph = module_network, 
                                   directed = FALSE, 
                                   # weights = FALSE, 
                                   normalized = TRUE)

closeness <- igraph::closeness(graph = module_network, normalized = TRUE)

names(degree)
names(betweenness)
names(closeness)

importance <- 
  data.frame(degree, betweenness, closeness, stringsAsFactors = FALSE) %>% 
  apply(2, function(x){
    x/sd(x)
  }) %>% 
  as.data.frame() %>% 
  data.frame(node = names(degree), .) %>% 
  dplyr::mutate(importance = (degree + betweenness + closeness) /3 )

module_node <- 
  module_node %>% 
  dplyr::left_join(importance, by = "node")

module_network <- 
  igraph::graph_from_data_frame(d = module_edge, 
                                directed = FALSE, 
                                vertices = module_node)

save(module_network, file = 'module_network')
load("module_network")
plot <- 
  ggraph::ggraph(graph = module_network, layout = "auto") +
  ggraph::geom_edge_link(aes(color = class)) +
  ggraph::scale_edge_color_manual(values = c("c-c" = "black", "p-c" = "grey")) +
  ggraph::geom_node_point(aes(color = class, shape = class, size = importance)) +
  scale_color_manual(values = c("Detected metabolite" = "#FB8072",
                                "Hidden metabolite" = "#D9D9D9",
                                "Peak" = "#80B1D3")) +
  scale_size_continuous(range = c(1,5)) +
  theme_void() +
  theme(plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA)
  )
plot

save(module_network, file = "module_network")

ggsave(plot, filename = "all_module_netwprk.png",
       width = 9, height = 7, 
       bg = "transparent")

ggsave(plot, filename = "all_module_netwprk.pdf",
       width = 9, height = 7, 
       bg = "transparent")


##pathway analysis for module
load("hsa_pathway")

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
