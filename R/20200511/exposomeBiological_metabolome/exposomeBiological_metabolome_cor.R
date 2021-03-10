##avoid source
no_function()

##load data
sxtTools::setwd_project()
library(tidyverse)
rm(list = ls())

load("data_20200511/metabolome/variable_info")
met_variable_info <- variable_info

met_annotation_info <- readr::read_csv("data_20200511/metabolome/annotation_table.csv")

setwd("data_analysis/exposomeBiological_metabolome/")

load("met_expression_data")
load("microbiome_expression_data")
load("met_sample_info")
load("microbiome_sample_info")

load("microbiome_variable_info")

met_sample_info$sample_id == colnames(met_expression_data)
colnames(met_expression_data) <- as.character(met_sample_info$CollectionDate)
colnames(microbiome_expression_data) <- colnames(met_expression_data)

###only remain species
microbiome_variable_info$variable_id == rownames(microbiome_expression_data)

microbiome_variable_info <- 
  microbiome_variable_info %>% 
  dplyr::filter(level == "genus")

microbiome_expression_data <- 
microbiome_expression_data[match(microbiome_variable_info$variable_id, 
                                 rownames(microbiome_expression_data)),]

dim(microbiome_expression_data)

#######correlation analysis
microbiome_expression_data <- log(microbiome_expression_data + 1, 2)
met_expression_data <- log(met_expression_data + 1, 2)

dim(microbiome_expression_data)
dim(met_expression_data)

###correct fiber for metabolomics
met_expression_data1 <- 
  purrr::map(as.data.frame(t(met_expression_data)), .f = function(x){
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
  }) %>% 
  do.call(rbind, .) %>% 
  as.data.frame()

temp_data <- 
  apply(met_expression_data1, 1, function(x){
    (x - mean(x))/sd(x)
  }) %>% 
  t() 

colnames(temp_data) <- 
colnames(met_expression_data1) <-
  colnames(met_expression_data)

##heatmap of metabolome
library(circlize)

col_fun = colorRamp2(breaks = seq(min(temp_data), max(temp_data), length.out = 90),
                     colors = 
                       viridis::magma(n = 100)[-c(1:10)], 
                     transparency = 0
)

plot <- 
temp_data %>% 
  ComplexHeatmap::Heatmap(cluster_columns = FALSE,
                          show_column_names = TRUE,
                          show_row_names = FALSE,
                          clustering_method_rows = "ward.D",
                          clustering_method_columns = "ward.D",
                          clustering_distance_columns = "euclidean",
                          clustering_distance_rows = "euclidean",
                          col = col_fun,
                          km = 2, border = TRUE, 
                          row_dend_reorder = TRUE, 
                          column_dend_reorder = TRUE,
                          column_names_rot = 45, 
                          name = "Z-score")

plot <- ggplotify::as.ggplot(plot)

plot

# ggsave(plot, filename = "metabolome_plot/metabolome_heatmap.pdf", width = 7, height = 7)

###remove some microbiome variables
remain_idx <-
  apply(microbiome_expression_data, 1, function(x) {
    length(unique(x)) != 1
  }) %>%
  which()

microbiome_expression_data <- 
  microbiome_expression_data[remain_idx,]

microbiome_variable_info <- 
  microbiome_variable_info[remain_idx,]

###heatmap of microbiome
plot <- 
  microbiome_expression_data %>% 
  apply(1, function(x){
    (x - mean(x))/sd(x)
  }) %>% 
  t() %>% 
  ComplexHeatmap::Heatmap(cluster_columns = FALSE,
                          show_column_names = TRUE,
                          show_row_names = FALSE,
                          clustering_method_rows = "ward.D",
                          clustering_method_columns = "ward.D",
                          clustering_distance_columns = "euclidean",
                          clustering_distance_rows = "euclidean",
                          col = col_fun,
                          km = 2, border = TRUE, 
                          row_dend_reorder = TRUE, 
                          column_dend_reorder = TRUE,
                          column_names_rot = 45, 
                          name = "Z-score")

plot <- ggplotify::as.ggplot(plot)

plot

# ggsave(plot, filename = "microbiome_plot/exposome_heatmap.pdf", width = 7, height = 7)


# ####calculate correlation between metabolome and exposome
# cor_value <-
#   cor(x = t(as.matrix(microbiome_expression_data)),
#       y = t(as.matrix(met_expression_data1)),
#       method = "spearman")
# 
# cor_value <-
#   cor_value %>%
#   as.data.frame() %>%
#   tibble::rownames_to_column(var = "from") %>%
#   tidyr::pivot_longer(-from, names_to = "to", values_to = "cor")
# 
# library(plyr)
# 
# p_value <-
#   purrr::map(
#     as.data.frame(t(cor_value)),
#     .f = function(x) {
#       value1 <- as.numeric(microbiome_expression_data[x[1], ])
#       value2 <- as.numeric(met_expression_data1[x[2], ])
#       cor.test(value1, value2, method = "spearman")$p.value
#     }
#   ) %>%
#   unlist()
# 
# cor_value <-
#   data.frame(cor_value, p_value, stringsAsFactors = FALSE)
# 
# plot(density(cor_value$p_value))
# 
# library(plyr)
# cor_value <-
#   cor_value %>%
#   plyr::dlply(.variables = .(from)) %>%
#   purrr::map(
#     .f = function(x) {
#       x <- x %>%
#         dplyr::filter(abs(cor) > 0.9)
#       fdr <- p.adjust(x$p_value, method = "fdr")
#       x <-
#         data.frame(x, fdr, stringsAsFactors = FALSE)
#       x
#     }
#   )
# 
# cor_value <-
#   cor_value %>%
#   do.call(rbind, .) %>%
#   as.data.frame()
# 
# cor_value <-
#   cor_value %>%
#   dplyr::filter(abs(cor) > 0.9 & cor != 1)
# 
# dim(cor_value)
# 
# save(cor_value, file = "cor_value")
load('cor_value')

cor_value1 <-
cor_value %>%
  dplyr::filter(abs(cor) > 0.9 & fdr < 0.05)

dim(cor_value1)
dim(cor_value)

sxtTools::setwd_project()
setwd("data_analysis/exposomeChemical_metabolome/exposomeChemical_metabolome_plot")
# for(idx in 11:20){
#   cat(idx, " ")
#   path1 <- file.path(cor_value1$from[idx])
#   dir.create(path1, showWarnings = FALSE)
#   temp_data <-
#     data.frame(date = as.character(met_sample_info$CollectionDate),
#                exp = as.numeric(microbiome_expression_data[cor_value1$from[idx],]),
#                met = as.numeric(met_expression_data1[cor_value1$to[idx],]),
#                stringsAsFactors = FALSE)
#   plot <-
#   temp_data %>%
#     ggplot(aes(exp, met)) +
#     geom_point() +
#     geom_smooth(method = "lm", color = "skyblue") +
#     ggrepel::geom_label_repel(aes(x = exp, met, label = date)) +
#     labs(x = paste("Exposome (Chemical): ", cor_value1$microbiome_id[idx], sep = ""),
#          y = paste("Metabolome: " ,cor_value1$met_id[idx]), sep = "") +
#     theme_bw() +
#     theme(axis.title = element_text(size = 13),
#           axis.text = element_text(size = 12),
#           plot.background = element_rect(fill = "transparent", color = NA),
#           panel.background = element_rect(fill = "transparent", color = NA)) +
#     annotate(geom = "text",
#              x = -Inf, y = Inf,
#              label = paste("Correlation: ",round(cor_value1$cor[idx], 2),
#                            "\nFDR adjusted P value: ",
#                            round(cor_value1$fdr[idx], 3), sep = ""),
#              vjust = 2, hjust = -1)
# 
#   name <- paste(cor_value1$from[idx], "_",
#                 cor_value1$to[idx], ".pdf", sep = "")
# 
#   ggsave(plot, filename = file.path(path1, name),
#          width = 7, height = 7, bg = "transparent")
# 
# }

sxtTools::setwd_project()  
setwd("data_analysis/exposomeBiological_metabolome")

cor_value1$from %>% unique()

##cluster exposome
###k-means
## fuzzy c-means clustring
###fuzzy k-means clustering
library(Mfuzz)

#first get the time point data together:
# bind that to the data frame
temp_data <- microbiome_expression_data

##scale
temp_data <- 
  temp_data %>% 
  apply(1, function(x) (x - mean(x))/sd(x)) %>% 
  t() %>% 
  as.data.frame()

time <- c(1:5)

temp_data <- rbind(time, temp_data)

row.names(temp_data)[1]<-"time"

#save it to a temp file so ti doesnt clutter up my blog directory
# tmp <- tempfile(tmpdir = "exposomeBiological_k_mean_clustering")
# 
# write.table(
#   temp_data,
#   file = tmp,
#   sep = '\t',
#   quote = F,
#   col.names = NA
# )
# 
# #read it back in as an expression set
# data <- table2eset(file=tmp)
# 
# # data.s <- standardise(data)
# data.s <- data
# 
# m1 <- mestimate(data.s)
# m1
# 
# Dmin(data.s, m=m1, crange=seq(2,22,1), repeats=3, visu=TRUE)
# 
# clust = 8
# 
# c <- mfuzz(data.s,c=clust,m=m1)
# save(c, file = file.path("exposomeBiological_k_mean_clustering","c"))
# load("exposomeBiological_k_mean_clustering/c")
# 
# mfuzz.plot(eset = data.s, 
#            # min.mem = 0.6,
#            cl = c, 
#            mfrow=c(2,4),
#            time.labels = time,
#            new.window = FALSE)
# 
# centers <- c$centers
# names(c$cluster) == rownames(c$membership)
# 
# cluster_info <-
#   data.frame(
#     variable_id = names(c$cluster),
#     c$membership,
#     cluster = c$cluster,
#     stringsAsFactors = FALSE
#   ) %>%
#   arrange(cluster)
# 
# openxlsx::write.xlsx(x = cluster_info, 
#                      file = file.path("exposomeBiological_k_mean_clustering/","cluster_info.xlsx"),
#                      asTable = TRUE)



# ####plot for each cluster
# for (cluster in 1:8) {
#   cat(cluster, " ")
#   cluster_data <-
#     cluster_info[cluster_info$cluster == cluster,] %>%
#     dplyr::select(1, 1 + cluster)
#   
#   colnames(cluster_data) <- c("variable_id", "membership")
#   
#   cluster_data <-
#     cluster_data %>%
#     dplyr::filter(membership > 0.3)
#   
#   openxlsx::write.xlsx(
#     cluster_data,
#     file =  file.path(
#       "exposomeBiological_k_mean_clustering/",
#       paste("cluster", cluster, ".xlsx", sep = "")
#     ),
#     asTable = TRUE
#   )
#   
#   temp_center <-
#     centers[cluster, , drop = TRUE] %>%
#     data.frame(ga = names(.),
#                value = .,
#                stringsAsFactors = FALSE) %>%
#     dplyr::mutate(ga = factor(ga, levels = ga))
#   
#   plot <-
#     temp_data[cluster_data$variable_id, ] %>%
#     apply(1, function(x) {
#       (x - mean(x)) / sd(x)
#     }) %>%
#     t() %>%
#     as.data.frame() %>%
#     data.frame(
#       membership = cluster_data$membership,
#       .,
#       stringsAsFactors = FALSE,
#       check.names = FALSE
#     ) %>%
#     tibble::rownames_to_column(var = "variable_id") %>%
#     tidyr::pivot_longer(
#       cols = -c(variable_id, membership),
#       names_to = "ga",
#       values_to = "value"
#     ) %>%
#     dplyr::mutate(ga = factor(ga, levels = unique(ga))) %>%
#     ggplot(aes(ga, value, group = variable_id)) +
#     geom_line(aes(color = membership)) +
#     scale_color_gradientn(colours = c(RColorBrewer::brewer.pal(n = 10, name = "OrRd")[c(1:7)])) +
#     labs(x = "", y = "Scaled intensity") +
#     theme_bw() +
#     theme(
#       legend.position = c(0, 1),
#       legend.justification = c(0, 1),
#       panel.grid.minor = element_blank(),
#       axis.title = element_text(size = 13),
#       axis.text = element_text(size = 12),
#       axis.text.x = element_text(
#         angle = 45,
#         hjust = 1,
#         vjust = 1,
#         size = 12
#       ),
#       panel.background = element_rect(fill = "transparent", color = NA),
#       plot.background = element_rect(fill = "transparent", color = NA),
#       legend.background = element_rect(fill = "transparent", color = NA)
#     ) +
#     geom_line(
#       mapping = aes(ga, value, group = 1),
#       data = temp_center,
#       size = 2
#     )
#   
#   plot
#   
#   ggsave(
#     plot,
#     filename =  file.path(
#       "exposomeBiological_k_mean_clustering/",
#       paste("cluster", cluster, ".pdf", sep = "")
#     ),
#     width = 7,
#     height = 7
#   )
# }



###correlation network for exposome biological and metabolome
cor_value1$from %>% unique() %>% length()
cor_value1$to %>% unique() %>% length()

library(igraph)
library(ggraph)
library(tidygraph)


###network for all the exposome biological and metabolome
edge_data <-  
  cor_value1 %>% 
  # dplyr::filter(from %in% cluster1) %>%
  dplyr::rename(from = from, 
                to = to, 
                Correlation = cor) %>% 
  dplyr::mutate(fdr = -log(fdr, 10))

node_data <- 
  cor_value1 %>% 
  # dplyr::filter(from %in% cluster1) %>%
  dplyr::rename(from = from, to = to) %>% 
  dplyr::select(from, to) %>% 
  tidyr::pivot_longer(cols = c(from, to), 
                      names_to = "class", values_to = "node") %>% 
  dplyr::mutate(class1 = case_when(
    stringr::str_detect(class, "from") ~ "Exposome biological",
    TRUE ~ "Metabolome"
  )) %>% 
  dplyr::select(node, class1) %>% 
  dplyr::rename(Class = class1) %>%
  dplyr::distinct(node, .keep_all = TRUE)

node_data <- 
  node_data %>% 
  dplyr::arrange(Class)

node_data <-
  node_data %>%
  dplyr::left_join(met_annotation_info[, c("name", "Compound.name")], by = c("node" = "name")) %>%
  dplyr::mutate(compound.name = case_when(!is.na(Compound.name) ~ Compound.name,
                                          is.na(Compound.name) ~ node)) %>%
  dplyr::select(node, Class, compound.name)

temp_data <- 
  tidygraph::tbl_graph(nodes = node_data, 
                       edges = edge_data,
                       directed = TRUE) %>% 
  dplyr::mutate(Degree = centrality_degree(mode = 'all'))

pal <-
  wesanderson::wes_palette(name = "Zissou1", n = 100, type = "continuous")

plot1 <-
  ggraph(temp_data,
         layout = 'linear',
         circular = TRUE) +
  geom_edge_arc(aes(color = Correlation
                    # width = fdr
                    ),
                show.legend = TRUE) +
  geom_node_point(aes(fill = Class,
                      size = Degree),
                  shape = 21,
                  show.legend = TRUE) +
  scale_fill_manual(values = c(
    "Exposome biological" = ggsci::pal_d3()(10)[5],
    "Metabolome" = ggsci::pal_d3()(10)[3]
  )) +
  # geom_node_text(
  #   aes(
  #     x = x * 1.05,
  #     y = y * 1.05,
  #     label = compound.name,
  #     hjust = 'outward',
  #     angle = -((-node_angle(x, y) + 90) %% 180) + 90,
  #     size = 3,
  #     colour = Class
  #   ),
  #   size = 3,
  #   alpha = 1, 
  #   show.legend = FALSE
  # ) +
  guides(edge_width = guide_legend(title = "-log10(FDR adjusted P value)", 
                                   override.aes = list(shape = NA)),
         edge_color = ggraph::guide_edge_colorbar(title = "Spearman correlation"),
         fill = guide_legend(title = "Class", 
                             override.aes = list(size = 4, linetype = "blank")),
         size = guide_legend(title = "Degree", override.aes = list(linetype = 0))) +
  ggraph::scale_edge_color_gradientn(colours = pal) +
  ggraph::scale_edge_width(range = c(0.1, 1)) +
  scale_size_continuous(range = c(0.3, 3)) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
  )

plot1

# ggsave(
#   plot1,
#   filename = "exposomeBiological_metabolome_correlation_network.pdf",
#   width = 8.5,
#   height = 7,
#   bg = "transparent"
# )




###only remain the top 25% correlation
###network for all the exposome biological and metabolome
edge_data <-  
  cor_value1 %>% 
  dplyr::rename(from = from, 
                to = to, 
                Correlation = cor) %>% 
  dplyr::mutate(fdr = -log(fdr, 10))

node_data <- 
  cor_value1 %>% 
  # dplyr::filter(from %in% cluster1) %>%
  dplyr::rename(from = from, to = to) %>% 
  dplyr::select(from, to) %>% 
  tidyr::pivot_longer(cols = c(from, to), 
                      names_to = "class", values_to = "node") %>% 
  dplyr::mutate(class1 = case_when(
    stringr::str_detect(class, "from") ~ "Exposome biological",
    TRUE ~ "Metabolome"
  )) %>% 
  dplyr::select(node, class1) %>% 
  dplyr::rename(Class = class1) %>%
  dplyr::distinct(node, .keep_all = TRUE)

node_data <- 
  node_data %>% 
  dplyr::arrange(Class)

node_data <-
  node_data %>%
  dplyr::left_join(met_annotation_info[, c("name", "Compound.name")], by = c("node" = "name")) %>%
  dplyr::mutate(compound.name = case_when(!is.na(Compound.name) ~ Compound.name,
                                          is.na(Compound.name) ~ node)) %>%
  dplyr::select(node, Class, compound.name)

top_node <- 
table(c(edge_data$from, edge_data$to)) %>% 
  `>`(15) %>% 
  which() %>% 
  names()

node_data <-
  node_data %>% 
  dplyr::filter(node %in% top_node)

edge_data <-
  edge_data %>%
  dplyr::filter(from %in% top_node & to %in% top_node)

node_data <-
 node_data %>% 
  dplyr::filter(node %in% unique(c(edge_data$from, edge_data$to)))


temp_data <- 
  tidygraph::tbl_graph(nodes = node_data, 
                       edges = edge_data,
                       directed = TRUE) %>% 
  dplyr::mutate(Degree = centrality_degree(mode = 'all'))

pal <-
  wesanderson::wes_palette(name = "Zissou1", n = 100, type = "continuous")

plot1 <-
  ggraph(temp_data,
         layout = 'linear',
         circular = TRUE) +
  geom_edge_arc(aes(color = Correlation
                    # width = fdr
  ),
  show.legend = TRUE) +
  geom_node_point(aes(fill = Class,
                      size = Degree),
                  shape = 21,
                  show.legend = TRUE) +
  scale_fill_manual(values = c(
    "Exposome biological" = ggsci::pal_d3()(10)[5],
    "Metabolome" = ggsci::pal_d3()(10)[3]
  )) +
  # geom_node_text(
  #   aes(
  #     x = x * 1.05,
  #     y = y * 1.05,
  #     label = compound.name,
  #     hjust = 'outward',
  #     angle = -((-node_angle(x, y) + 90) %% 180) + 90,
  #     size = 3,
  #     colour = Class
  #   ),
  #   size = 3,
#   alpha = 1, 
#   show.legend = FALSE
# ) +
guides(edge_width = guide_legend(title = "-log10(FDR adjusted P value)", 
                                 override.aes = list(shape = NA)),
       edge_color = ggraph::guide_edge_colorbar(title = "Spearman correlation"),
       fill = guide_legend(title = "Class", 
                           override.aes = list(size = 4, linetype = "blank")),
       size = guide_legend(title = "Degree", override.aes = list(linetype = 0))) +
  ggraph::scale_edge_color_gradientn(colours = pal) +
  ggraph::scale_edge_width(range = c(0.1, 1)) +
  scale_size_continuous(range = c(1, 5)) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
  )

plot1

# ggsave(
#   plot1,
#   filename = "exposomeBiological_metabolome_correlation_network_example.pdf",
#   width = 8.5,
#   height = 7,
#   bg = "transparent"
# )