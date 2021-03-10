##avoid source
no_function()

##load data
sxtTools::setwd_project()
library(tidyverse)
rm(list = ls())

source("R/20200511/exposomeChemical_metabolome/module_tools.R")
load("data_20200511/exposome/variable_info")
exp_variable_info <- variable_info

load("data_20200511/metabolome/variable_info")
met_variable_info <- variable_info

met_annotation_info <-
  readr::read_csv("data_20200511/metabolome/annotation_table.csv")

setwd("data_analysis/exposomeChemical_metabolome/")

load("met_expression_data")
load("exp_expression_data")
load("met_sample_info")
load("exp_sample_info")

met_sample_info$sample_id == colnames(met_expression_data)
colnames(met_expression_data) <-
  as.character(met_sample_info$CollectionDate)
colnames(exp_expression_data) <- colnames(met_expression_data)

#######correlation analysis
exp_expression_data <- log(exp_expression_data + 1, 2)
met_expression_data <- log(met_expression_data + 1, 2)

dim(exp_expression_data)
dim(met_expression_data)

###correct fiber for metabolomics
met_expression_data1 <-
  purrr::map(
    as.data.frame(t(met_expression_data)),
    .f = function(x) {
      temp_data <-
        data.frame(fiber = c(0,
                             0,
                             0,
                             1,
                             1),
                   x,
                   stringsAsFactors = FALSE)
      lm_result <- lm(formula = x ~ fiber, data = temp_data)
      lm_result$residuals
    }
  ) %>%
  do.call(rbind, .) %>%
  as.data.frame()

temp_data <-
  apply(met_expression_data1, 1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t()

colnames(temp_data) <-
  colnames(met_expression_data1) <-
  colnames(met_expression_data)

load("kegg.adduct.pos.rda")
load("kegg.adduct.neg.rda")
load("kegg_network")
load("kegg.compound.rda")

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

met_variable_info <-
  met_variable_info %>% 
  dplyr::mutate(mz = stringr::str_split(Mass, "_") %>% 
                  lapply(function(x) x[1]) %>% 
                  unlist() %>% 
                  as.numeric(),
                rt = stringr::str_split(Mass, "_") %>% 
                  lapply(function(x) x[2]) %>% 
                  unlist() %>% 
                  as.numeric())
  
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

# mz.pos <-
#   marker_info %>% dplyr::filter(stringr::str_detect(name, "p")) %>% pull(mz)
# 
# rt.pos <-
#   marker_info %>% dplyr::filter(stringr::str_detect(name, "p")) %>% pull(rt)
# 
# mz.neg <-
#   marker_info %>% dplyr::filter(stringr::str_detect(name, "n")) %>% pull(mz)
# 
# rt.neg <-
#   marker_info %>% dplyr::filter(stringr::str_detect(name, "n")) %>% pull(mz)
# 
# annotation_result_pos <-
#   keggMatch(
#     mz = mz.pos,
#     rt = rt.pos,
#     metabolite = kegg.adduct.pos,
#     mz.tol = 25,
#     rt.tol = 1000000,
#     polarity = "positive"
#   )
# 
# dir.create("community_analysis")
# 
# save(annotation_result_pos, file = "community_analysis/annotation_result_pos")

load("community_analysis/annotation_result_pos")

# annotation_result_neg <-
#   keggMatch(
#     mz = mz.neg,
#     rt = rt.neg,
#     metabolite = kegg.adduct.neg,
#     mz.tol = 25,
#     rt.tol = 1000000,
#     polarity = "negative"
#   )
# 
# save(annotation_result_neg, file = "community_analysis/annotation_result_neg")

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

# initial_modules <-
#   getModule(node.id = metabolite.id,
#             metabolic.network = kegg_network,
#             threads = 3,
#             max.reaction = 3,
#             progressbar = TRUE,
#             calculate.impact = TRUE)
# 
# save(initial_modules, file = "community_analysis/initial_modules")
load("community_analysis/initial_modules") 

# ref.activity.score <- getNullDistribution(
#   all.mz.pos = met_variable_info %>%
#     dplyr::filter(stringr::str_detect(peak_name, "p")) %>%
#     pull(mz),
#   all.mz.neg = met_variable_info %>%
#     dplyr::filter(stringr::str_detect(peak_name, "n")) %>%
#     pull(mz),
#   mz.tol = 25,
#   metabolic.network = kegg_network,
#   size.pos = length(mz.pos),
#   size.neg = length(mz.neg),
#   times = 100,
#   threads = 3,
#   progressbar = FALSE,
#   kegg.adduct.pos = kegg.adduct.pos,
#   kegg.adduct.neg = kegg.adduct.neg
# )
# 
# 
# save(ref.activity.score, file = "community_analysis/ref.activity.score")
load("community_analysis/ref.activity.score")

ref.as <- unlist(ref.activity.score) + 0.00000001

activity.score <- initial_modules$activity.score + 0.00000001
###
para <- MASS::fitdistr(x = ref.as, densfun = "gamma")[[1]]
# standard.distribution <- rgamma(n = length(ref.as), shape = para[1], rate = para[2])
standard.distribution <- rgamma(n = 100000, shape = para[1], rate = para[2])
# ks.test(x = ref.as, y = standard.distribution)

cdf <- ecdf(x = standard.distribution)
module.p <- 1 - cdf(activity.score)

which(module.p < 0.05)

plot = 
data.frame(x = standard.distribution) %>% 
  ggplot() +
  geom_density(aes(x = x)) +
  theme_bw() +
  ggrepel::geom_label_repel(aes(x = score, y = 0, label = module),
                            data = data.frame(module = paste("Module", 1:length(activity.score)),
                                              score = activity.score)[3,,drop = FALSE]) +
  geom_vline(xintercept = activity.score[3], color = "red") +
  labs(x = "Activity score", y = "Density") +
  theme(panel.grid.minor = element_blank(),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12))
plot
# ggsave(plot, filename = "module.p.pdf", width = 7, height = 7)


# 
# module <- initial_modules$module[module.p < 0.05]
# module
# 
# save(module, file = "community_analysis/module")
load("community_analysis/module")

##nodule is the significant modules
library(ggraph)

module <- module[1]

##module 1
module1 <-
  data.frame(id = module[[1]], stringsAsFactors = FALSE) %>%
  mutate(detected = case_when(id %in% initial_modules$hidden.id ~ "No",
                              TRUE ~ 'Yes',))

# module2 <-
#   data.frame(id = module[[2]], stringsAsFactors = FALSE) %>%
#   mutate(detected = case_when(id %in% initial_modules$hidden.id ~ "No",
#                               TRUE ~ 'Yes',))

all_module <-
  data.frame(id = unlist(module) %>% unique(), 
             stringsAsFactors = FALSE) %>%
  mutate(detected = case_when(id %in% initial_modules$hidden.id ~ "No",
                              TRUE ~ 'Yes',))

##construct network according to modules
# names(annotation_result_pos) <- 
#   marker_info %>% 
#   dplyr::filter(stringr::str_detect(name, "p")) %>% 
#   pull(name)
# 
# names(annotation_result_neg) <- 
#   marker_info %>% 
#   dplyr::filter(stringr::str_detect(name, "n")) %>% 
#   pull(name)
# 
# annotation_result_pos <- 
#   purrr::map2(
#     .x = annotation_result_pos,
#     .y = names(annotation_result_pos),
#     .f = function(x, y) {
#       if(is.null(x)) return(NULL)
#       data.frame(peak = y, x[,c(1,6,7,12)])
#     }
#   )
# 
# annotation_result_pos <- 
#   annotation_result_pos %>% 
#   dplyr::bind_rows() %>% 
#   distinct()
# 
# annotation_result_pos <- 
#   annotation_result_pos %>% 
#   dplyr::filter(ID %in% unique(unlist(module)))
# 
# annotation_result_neg <- 
#   purrr::map2(
#     .x = annotation_result_neg,
#     .y = names(annotation_result_neg),
#     .f = function(x, y) {
#       if(is.null(x)) return(NULL)
#       data.frame(peak = y, x[,c(1,6,7,12)])
#     }
#   )
# 
# annotation_result_neg <- 
#   annotation_result_neg %>% 
#   dplyr::bind_rows() %>% 
#   distinct()
# 
# 
# annotation_result_neg <- 
#   annotation_result_neg %>% 
#   dplyr::filter(ID %in% unique(unlist(module)))
# 
# 
# annotation_result <- rbind(annotation_result_pos, annotation_result_neg)

# save(annotation_result, file = "annotation_result")
load("annotation_result")

##remove redundant annotations
library(plyr)
annotation_result <- 
annotation_result %>%
  dplyr::distinct(ID, Adduct, Formula, mz.error, .keep_all = TRUE) %>%
  plyr::dlply(.variables = "ID") %>%
  purrr::map(
    .f = function(x) {
      if (nrow(x) == 1) {
        return(x)
      } else{
        if (any(x$Adduct == "M+H") | any(x$Adduct == "M-H")) {
          x <-
            x %>%
            dplyr::filter(Adduct == "M+H" | Adduct == "M-H")
        } else{
          x <-
            x %>%
            dplyr::filter(mz.error == min(mz.error))
        }
      }
      return(x)
    }
  ) %>% 
  do.call(rbind, .) %>% 
  as.data.frame()
  
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
  igraph::degree(graph = module_network, mode = "total", normalized = FALSE)

betweenness <- igraph::betweenness(graph = module_network,
                                   directed = FALSE,
                                   # weights = FALSE,
                                   normalized = TRUE)

closeness <- igraph::closeness(graph = module_network, normalized = TRUE)

names(degree)
# names(betweenness)
# names(closeness)

importance <- 
  data.frame(degree, betweenness, closeness, stringsAsFactors = FALSE) %>% 
  data.frame(node = names(degree), .) %>% 
  dplyr::mutate(importance = degree )

module_node <-
  module_node %>%
  dplyr::left_join(importance, by = "node") %>%
  dplyr::left_join(kegg.compound[,c("ID", "Name")], by = c("node" = "ID")) %>%
  dplyr::mutate(true_name = stringr::str_split(Name, ";") %>% 
                  lapply(function(x)x[1]) 
                %>% unlist())

module_network <- 
  igraph::graph_from_data_frame(d = module_edge, 
                                directed = FALSE, 
                                vertices = module_node)

# save(module_network, file = 'module_network')
load("module_network")

plot <-
  ggraph::ggraph(graph = module_network, layout = "kk") +
  ggraph::geom_edge_link(aes(color = class)) +
  ggraph::scale_edge_color_manual(values = c("c-c" = "black", "p-c" = "grey")) +
  ggraph::geom_node_point(aes(
    color = class,
    shape = class,
    size = importance
  )) +
  geom_node_text(aes(label = ifelse(
    class == "Detected metabolite",
    true_name, NA
  )),
  repel = TRUE) +
  scale_color_manual(
    values = c(
      "Detected metabolite" = ggsci::pal_d3()(n = 10)[4],
      "Hidden metabolite" = ggsci::pal_d3()(n = 10)[6],
      "Peak" = ggsci::pal_d3()(n = 10)[8]
    )
  ) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  scale_size_continuous(range = c(1, 10)) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  )

plot

# ggsave(plot, filename = "all_module_netwprk.pdf",
#        width = 9, height = 7, 
#        bg = "transparent")

##pathway analysis for module
load("hsa_pathway")

# all_module_pathway <-
#   enrich_pathway(id = module_node %>% 
#                    dplyr::filter(class == "Detected metabolite") %>% 
#                    dplyr::pull(node),
#                  pathway_database = hsa_pathway)
# 
# save(all_module_pathway, file = "all_module_pathway")
load("all_module_pathway")

all_module_pathway$Pathway.name <- 
  all_module_pathway$Pathway.name %>% 
  stringr::str_replace_all(pattern = " - Homo sapiens \\(human\\)", replacement = "")

# library(openxlsx)
# wb = createWorkbook()
# modifyBaseFont(wb, fontSize = 12, fontName = "Arial Narrow")
# addWorksheet(wb, sheetName = "KEGG pathway", gridLines = TRUE)
# addWorksheet(wb, sheetName = "Module node", gridLines = TRUE)
# addWorksheet(wb, sheetName = "Module edge", gridLines = TRUE)
# freezePane(wb, sheet = 1, firstRow = TRUE, firstCol = TRUE) 
# freezePane(wb, sheet = 2, firstRow = TRUE, firstCol = TRUE) 
# freezePane(wb, sheet = 3, firstRow = TRUE, firstCol = TRUE) 
# writeDataTable(wb, sheet = 1, x = all_module_pathway,
#                colNames = TRUE, rowNames = FALSE)
# writeDataTable(wb, sheet = 2, x = module_node,
#                colNames = TRUE, rowNames = FALSE)
# writeDataTable(wb, sheet = 3, x = module_edge,
#                colNames = TRUE, rowNames = FALSE)
# saveWorkbook(wb, "Exposome_chemical_metabolome_enriched_pahtway.xlsx", overwrite = TRUE)

plot <- 
  all_module_pathway %>% 
  dplyr::filter(p.value.fdr < 0.01) %>% 
  dplyr::filter(Overlap > 3) %>% 
  # `[`(1:30, ,) %>% 
  dplyr::arrange(p.value.fdr) %>% 
  dplyr::mutate(Pathway.name = factor(Pathway.name, levels = Pathway.name)) %>% 
  dplyr::mutate(class = case_when(
    p.value.fdr < 0.01 & Overlap >= 5 ~ "Yes",
    TRUE ~ "No"
  )) %>% 
  ggplot(aes(y = Pathway.name, x = -log(p.value.fdr, 10))) +
  geom_segment(aes(y = Pathway.name, x = 0, 
                   yend = Pathway.name, 
                   xend = -log(p.value.fdr, 10), color = class), 
               show.legend = FALSE) +
  geom_point(aes(size = Overlap, 
                 fill = class), 
             shape = 21,
             show.legend = TRUE) +
  scale_fill_manual(values = c("Yes" = ggsci::pal_aaas()(n=10)[2], 
                                "No" = "black")) +
  scale_color_manual(values = c("Yes" = ggsci::pal_aaas()(n=10)[2], 
                               "No" = "black")) +
  scale_size_continuous(range = c(5,15)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme_classic() +
  labs(x = "-log10(P-value, FDR)", y = "") +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12), 
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.position = c(1,0), legend.justification = c(1,0),
        legend.background = element_blank(),
        strip.background = element_rect(fill = "#0099B47F"),
        strip.text = element_text(color = "white", size = 13),
        axis.text.x = element_text(
          angle = 45, hjust = 1, vjust = 1
        )) +
  ggplot2::coord_flip()

plot

# ggsave(plot, filename = "all_module_pathway.pdf",
#        height = 7, width = 13, bg = "transparent")




##for each pathway, it is positive or negative correlation with exposome chemical
all_id <-
module_edge %>% 
  dplyr::filter(class == "p-c")

metabolite_id <- 
all_module_pathway %>% 
  dplyr::filter(p.value.fdr < 0.05 & Overlap > 3) %>% 
  dplyr::pull(Pathway.ID) %>% 
  purrr::map(function(x){
    temp_id <-
      hsa_pathway[[grep(x, names(hsa_pathway))]][hsa_pathway[[grep(x, names(hsa_pathway))]] %in% all_id$to]
    temp_variable <-
      all_id %>% 
      dplyr::filter(to %in% temp_id)
    x <- 
    temp_variable %>% 
      dplyr::select(from) %>% 
      dplyr::left_join(cor_value1[,c(1,2,3)] %>% 
                         dplyr::rename(chemical_id = from, 
                                       metabolite_id = to, 
                                       correlation = cor), 
                       by = c("from" = "metabolite_id"))
    colnames(x) <- c("metabolite_id", "chemical_id", "correlation")
    x
  })

metabolite_id <-
  purrr::map2(
    .x = metabolite_id,
    .y = all_module_pathway %>% 
      dplyr::filter(p.value.fdr < 0.05 & Overlap > 3) %>% 
      dplyr::pull(Pathway.name),
    .f = function(x, y) {
      data.frame(x, pathway = y) %>%
        dplyr::select(chemical_id, metabolite_id, pathway, correlation)
    }
  )

metabolite_id <-
  metabolite_id %>%
  do.call(rbind, .) %>%
  dplyr::distinct(.keep_all = TRUE) %>%
  dplyr::left_join(exp_variable_info[,c(1,2)],
                   by = c("chemical_id" = "peak_ID")) %>% 
    dplyr::select(MetabID, everything()) %>% 
    dplyr::select(-chemical_id) %>% 
    dplyr::rename(chemical_id = MetabID)


metabolite_id <- 
  data.frame(metabolite_id, pathway_class = "KEGG")

metabolite_id_data <-
  rbind(metabolite_id)

annotation_result_in_module =
  module_edge %>%
  dplyr::filter(class == "p-c") %>% 
  dplyr::left_join(kegg.compound[,c(1,2)], by = c("to" = "ID"))

annotation_result_in_module = 
annotation_result_in_module %>% 
  dplyr::left_join(met_annotation_info, by = c("from" = "name"))
  
# library(openxlsx)
# wb = createWorkbook()
# modifyBaseFont(wb, fontSize = 12, fontName = "Arial Narrow")
# addWorksheet(wb, sheetName = "Annotation result in module", gridLines = TRUE)
# freezePane(wb, sheet = 1, firstRow = TRUE, firstCol = TRUE) 
# writeDataTable(wb, sheet = 1, x = annotation_result_in_module,
#                colNames = TRUE, rowNames = FALSE)
# saveWorkbook(wb, "Annotation_result_in_module.xlsx", overwrite = TRUE)


###only remain immune system
metabolite_id_data2 <-
  metabolite_id_data
  # dplyr::filter(
  #   pathway_class == "GO" &
  #     pathway %in% c(
  #       "acute inflammatory response",
  #       "humoral immune response",
  #       "complement activation",
  #       "regulation of humoral immune response"
  #     ) |
  #     pathway_class == "KEGG" &
  #     pathway %in% c("Complement and coagulation cascades") |
  #     pathway_class == "Reactome" &
  #     pathway %in% c(
  #       "Complement cascade",
  #       "Regulation of Complement cascade",
  #       "Platelet activation, signaling and aggregation"
  #     )
  # )

temp_data1 <- 
  metabolite_id_data2 %>%
  dplyr::select(pathway, chemical_id, correlation) %>%
  dplyr::mutate(affect = case_when(correlation > 0 ~ "pos",
                                   correlation < 0 ~ "neg")) %>%
  dplyr::select(-c(correlation)) %>%
  dplyr::distinct() %>%
  dplyr::select(chemical_id, affect) %>%
  dplyr::mutate(chemical_id = factor(chemical_id, levels = names(sort(
    table(chemical_id),
    decreasing =
      FALSE
  )))) 

pal <-
  wesanderson::wes_palette(name = "Zissou1", n = 100, type = "continuous")

plot1 <- 
  temp_data1 %>% 
  ggplot(aes(y = chemical_id)) +
  geom_bar(aes(fill = affect), show.legend = FALSE,
           color = "black") +
  scale_x_reverse(expand = expansion(mult = c(0.1,0))) +
  scale_fill_manual(values = c("pos" = pal[100],
                               "neg" = pal[1])) +
  theme_bw() +
  labs(y = "", x= "Pathway number") +
  theme(panel.grid = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        axis.ticks.y = element_blank())

plot1

# ggsave(plot1, filename = "plot1.pdf", width = 7, height = 7)

temp_data2 <- 
  metabolite_id_data2 %>%
  dplyr::select(pathway, chemical_id, correlation) %>%
  dplyr::mutate(affect = case_when(correlation > 0 ~ "pos",
                                   correlation < 0 ~ "neg")) %>%
  dplyr::select(-c(correlation)) %>%
  dplyr::distinct() %>%
  dplyr::select(pathway, affect) %>%
  dplyr::mutate(pathway = factor(pathway, levels = names(sort(
    table(pathway),
    decreasing =
      FALSE
  )))) 

plot2 <- 
  temp_data2 %>% 
  ggplot(aes(y = pathway)) +
  geom_bar(aes(fill = affect), show.legend = FALSE,
           color = "black") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_fill_manual(values = c("pos" = pal[100],
                               "neg" = pal[1])) +
  theme_bw() +
  labs(y = "", x= "Compound number") +
  theme(panel.grid = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        axis.ticks.y = element_blank())

plot2

# ggsave(plot2, filename = "plot2.pdf", width = 7, height = 7)

edge_data <-
  rbind(
    metabolite_id_data2[, c("chemical_id", "metabolite_id", "correlation")] %>%
      dplyr::distinct(chemical_id, metabolite_id, .keep_all = TRUE) %>%
      dplyr::rename(from = chemical_id,
                    to = metabolite_id,
                    cor = correlation),
    metabolite_id_data2[, c("metabolite_id", "pathway", "correlation")] %>%
      dplyr::distinct(pathway, metabolite_id, .keep_all = TRUE) %>%
      dplyr::rename(from = metabolite_id,
                    to = pathway,
                    cor = correlation) %>%
      dplyr::mutate(cor = NA)
  )

node_data <-
  rbind(
    metabolite_id_data2[, "chemical_id", drop = FALSE] %>%
      dplyr::distinct() %>%
      dplyr::mutate(class = "chemical") %>%
      dplyr::rename(node = chemical_id),
    metabolite_id_data2[, "metabolite_id", drop = FALSE] %>%
      dplyr::distinct() %>%
      dplyr::mutate(class = "metabolite") %>%
      dplyr::rename(node = metabolite_id),
    metabolite_id_data2[, "pathway", drop = FALSE] %>%
      dplyr::distinct() %>%
      dplyr::mutate(class = "pathway") %>%
      dplyr::rename(node = pathway)
  )

node_data[node_data$class == "chemical", ] <-
  node_data[node_data$class == "chemical", ][match(levels(temp_data1$chemical_id),
                                                   node_data[node_data$class == "chemical", ]$node), ]

node_data[node_data$class == "pathway", ] <-
  node_data[node_data$class == "pathway", ][match(levels(temp_data2$pathway),
                                                  node_data[node_data$class == "pathway", ]$node), ]

node_data <- 
  node_data %>% 
  dplyr::left_join(annotation_result[,1:2] %>% 
                     dplyr::distinct(peak ,.keep_all = TRUE), 
                   by = c("node" = "peak")) %>% 
  dplyr::rename(true_name = ID) %>% 
  dplyr::mutate(true_name = case_when(
    !is.na(true_name) ~ true_name, 
    TRUE ~ node
  ))

node_data <- 
  node_data %>% 
  dplyr::mutate(node = factor(node, levels = node))

node_data$node

library(tidygraph)

temp_data <- 
  tidygraph::tbl_graph(nodes = node_data, 
                       edges = edge_data,
                       directed = TRUE) %>% 
  dplyr::mutate(Degree = centrality_degree(mode = 'all'))

pal <-
  wesanderson::wes_palette(name = "Zissou1", n = 100, type = "continuous")


my_layout <- create_layout(temp_data, 
                           layout = 'linear', )

my_layout$y[my_layout$class == "metabolite"] <- 5

my_layout$y[my_layout$class == "pathway"] <- 10

my_layout1 <-
  my_layout

my_layout1$x <-
  my_layout$y

my_layout1$y <-
  my_layout$x

# my_layout1[my_layout1$class == "chemical",] <-
#   dplyr::left_join(data.frame(node = levels(temp_data1$chemical_id)),
#                    my_layout1[my_layout1$class == "chemical", ], by = "node") %>%
#   dplyr::select(x, y, node, everything())
# 
# my_layout1[my_layout1$class == "pathway", ] <-
#   dplyr::left_join(data.frame(node = levels(temp_data2$pathway)),
#                    my_layout1[my_layout1$class == "pathway",], by = "node") %>%
#   dplyr::select(x, y, node, everything())

my_layout1$y[my_layout1$class == "chemical"] <-
  seq(from = 1, to = 100, length.out = sum(my_layout1$class == "chemical"))

my_layout1$y[my_layout1$class == "metabolite"] <-
  my_layout1$y[my_layout1$class == "metabolite"] <-
  seq(from = 1, to = 100, length.out = sum(my_layout1$class == "metabolite"))

my_layout1$y[my_layout1$class == "pathway"] <-
  my_layout1$y[my_layout1$class == "pathway"] <-
  seq(from = 1, to = 100, length.out = sum(my_layout1$class == "pathway"))

plot <-
  ggraph(my_layout1) +
  geom_edge_link(aes(color = cor),
                 show.legend = FALSE) +
  geom_node_point(aes(fill = class,
                      size = Degree,
                      shape = class),
                  show.legend = FALSE) +
  scale_shape_manual(values = c(
    "chemical" = 21,
    "metabolite" = 21,
    "pathway" = 22
  )) +
  scale_fill_manual(
    values = c(
      "chemical" = ggsci::pal_d3()(10)[2],
      "metabolite" = ggsci::pal_d3()(10)[3],
      "pathway" = "black"
    )
  ) +
  scale_color_manual(
    values = c(
      "chemical" = ggsci::pal_d3()(10)[2],
      "metabolite" = ggsci::pal_d3()(10)[3],
      "pathway" = "black"
    )
  ) +
  geom_node_text(
    aes(
      x = x * 1.05,
      y = y * 1,
      label = true_name,
      hjust = ifelse(class == "chemical", 1, 0),
      # angle = -((-node_angle(x, y) + 90) %% 180) + 90,
      size = 3,
      colour = class
    ),
    size = 3,
    alpha = 1, 
    show.legend = FALSE
  ) +
  guides(edge_color = ggraph::guide_edge_colorbar(title = "Spearman correlation"),
         fill = guide_legend(title = "Class", 
                             override.aes = list(size = 4, linetype = "blank")),
         size = guide_legend(title = "Degree", override.aes = list(linetype = 0))) +
  ggraph::scale_edge_color_gradientn(colours = pal) +
  ggraph::scale_edge_width(range = c(0.2, 2)) +
  scale_size_continuous(range = c(3, 10)) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
    # legend.position = c(1,0), legend.justification = c(1,0)
  )

plot

# ggsave(plot, filename = "chemical_metabolite_pathway.pdf", width = 15, height = 7)

