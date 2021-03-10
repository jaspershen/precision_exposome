##avoid source
no_function()

sxtTools::setwd_project()
library(tidyverse)
rm(list = ls())

##load data
###exposome
load("data_20200511/exposome/variable_info")
exp_variable_info <- variable_info

###lab_test
load("data_20200511/lab_test/variable_info")
lab_test_variable_info <- variable_info

setwd("data_analysis/exposomeChemical_blood_test/")

load("lab_test_expression_data")
load("exp_expression_data")
load("lab_test_sample_info")
load("exp_sample_info")

lab_test_sample_info$sample_id == colnames(lab_test_expression_data)
colnames(lab_test_expression_data) <-
  as.character(lab_test_sample_info$date.start)
colnames(exp_expression_data) <- colnames(lab_test_expression_data)

#######correlation analysis
exp_expression_data <- log(exp_expression_data + 1, 2)
lab_test_expression_data <- log(lab_test_expression_data + 1, 2)

dim(exp_expression_data)
dim(lab_test_expression_data)

###correct fiber for lab_test
lab_test_expression_data1 <-
  purrr::map(
    as.data.frame(t(lab_test_expression_data)),
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


###output the plot between raw data correct data
dir.create("raw_adjust_plot")

for(i in 1:nrow(lab_test_expression_data)){
  cat(i, " ")
  temp_data <-
    data.frame(raw = as.numeric(lab_test_expression_data[i,]),
               adjust = as.numeric(lab_test_expression_data1[i,]))
  correlation <- cor(temp_data$raw, temp_data$adjust)
  plot <- 
    ggplot(temp_data, aes(raw, adjust)) +
    geom_point() +
    geom_smooth(method = "lm") +
    labs(title = rownames(lab_test_expression_data)[i]) +
    annotate(geom = "text",
             x = -Inf,
             y = Inf,
             label = correlation, hjust = -1, vjust = 1)
  
  ggsave(plot, filename = file.path("raw_adjust_plot", 
                                    paste(rownames(lab_test_expression_data)[i], ".pdf", sep = "")),
         width = 7, height = 7)
  
  
}


temp_data <-
  lab_test_expression_data1

colnames(temp_data) <-
  colnames(lab_test_expression_data1) <-
  colnames(lab_test_expression_data)

##heatmap of lab_test
library(circlize)

col_fun = colorRamp2(
  breaks = seq(min(temp_data), max(temp_data), length.out = 90),
  colors =
    viridis::magma(n = 100)[-c(1:10)],
  transparency = 0
)

plot <-
  temp_data %>%
  ComplexHeatmap::Heatmap(
    cluster_columns = FALSE,
    show_column_names = TRUE,
    show_row_names = FALSE,
    clustering_method_rows = "ward.D",
    clustering_method_columns = "ward.D",
    clustering_distance_columns = "euclidean",
    clustering_distance_rows = "euclidean",
    col = col_fun,
    km = 2,
    border = TRUE,
    row_dend_reorder = TRUE,
    column_dend_reorder = TRUE,
    column_names_rot = 45,
    name = "Z-score"
  )

plot <- ggplotify::as.ggplot(plot)

plot

ggsave(plot,
       filename = "lab_test_plot/lab_test_heatmap.pdf",
       width = 7,
       height = 7)

###heatmap of exposome
temp_data <-
  exp_expression_data %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t()

col_fun = colorRamp2(
  breaks = seq(min(temp_data), max(temp_data), length.out = 90),
  colors =
    viridis::magma(n = 100)[-c(1:10)],
  transparency = 0
)

plot <-
  temp_data %>%
  ComplexHeatmap::Heatmap(
    cluster_columns = FALSE,
    show_column_names = TRUE,
    show_row_names = FALSE,
    clustering_method_rows = "ward.D",
    clustering_method_columns = "ward.D",
    clustering_distance_columns = "euclidean",
    clustering_distance_rows = "euclidean",
    col = col_fun,
    km = 2,
    border = TRUE,
    row_dend_reorder = TRUE,
    column_dend_reorder = TRUE,
    column_names_rot = 45,
    name = "Z-score"
  )

plot <- ggplotify::as.ggplot(plot)

plot

ggsave(plot,
       filename = "exposome_plot/exposome_heatmap.pdf",
       width = 7,
       height = 7)

# ###calculate correlation between lab_test and exposome
# cor_value <-
#   cor(x = t(as.matrix(exp_expression_data)),
#       y = t(as.matrix(lab_test_expression_data1)),
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
#   purrr::map(as.data.frame(t(cor_value)), .f = function(x){
#     value1 <- as.numeric(exp_expression_data[x[1],])
#     value2 <- as.numeric(lab_test_expression_data1[x[2],])
#     cor.test(value1, value2, method = "spearman")$p.value
#   }) %>%
#   unlist()
# 
# cor_value <-
#   data.frame(cor_value, p_value, stringsAsFactors = FALSE)
# 
# plot(density(cor_value$p_value))
# 
# library(plyr)
# cor_value <-
# cor_value %>%
#   plyr::dlply(.variables = .(from)) %>%
#   purrr::map(.f = function(x){
#     x <- x %>%
#       dplyr::filter(abs(cor) > 0.9)
#     fdr <- p.adjust(x$p_value, method = "fdr")
#     x <-
#       data.frame(x, fdr, stringsAsFactors = FALSE)
#     x
#   })
# 
# cor_value <-
# cor_value %>%
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

cor_value1 <-
  cor_value1 %>%
  dplyr::left_join(exp_variable_info[, c(1:2)], by = c("from" = "peak_ID")) %>%
  dplyr::rename(exp_id = MetabID) %>%
  dplyr::mutate(lab_test_id = to) %>%
  dplyr::filter(!is.na(lab_test_id))

sxtTools::setwd_project()
setwd("data_analysis/exposomeChemical_blood_test/exposomeChemical_lab_test_plot")

# for(idx in 1:nrow(cor_value1)){
#   cat(idx, " ")
#   path1 <- file.path(cor_value1$from[idx])
#   dir.create(path1, showWarnings = FALSE)
#   temp_data <-
#     data.frame(date = as.character(lab_test_sample_info$date.start),
#                exp = as.numeric(exp_expression_data[cor_value1$from[idx],]),
#                pro = as.numeric(lab_test_expression_data1[cor_value1$to[idx],]),
#                stringsAsFactors = FALSE)
#   plot <-
#   temp_data %>%
#     ggplot(aes(exp, pro)) +
#     geom_point() +
#     geom_smooth(method = "lm", color = "skyblue") +
#     ggrepel::geom_label_repel(aes(x = exp, pro, label = date)) +
#     labs(x = paste("Exposome (chemical): ", cor_value1$exp_id[idx], sep = ""),
#          y = paste("Lab test: " ,cor_value1$lab_test_id[idx]), sep = "") +
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
setwd("data_analysis/exposomeChemical_blood_test/")

cor_value1$from %>% unique()

###correlation network for pro and exp
cor_value1$from %>% unique() %>% length()
cor_value1$to %>% unique() %>% length()

library(igraph)
library(ggraph)
library(tidygraph)

###network for all the exposome and lab_test
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
                      names_to = "class",
                      values_to = "node") %>%
  dplyr::mutate(class1 = case_when(
    stringr::str_detect(class, "from") ~ "Exposome",
    TRUE ~ "lab_test"
  )) %>%
  dplyr::select(node, class1) %>%
  dplyr::rename(Class = class1) %>%
  dplyr::distinct(node, .keep_all = TRUE)

node_data <-
  node_data %>%
  dplyr::arrange(Class)

node_data <-
  node_data %>%
  dplyr::left_join(exp_variable_info[, c("peak_ID", "MetabID")], by = c("node" = "peak_ID")) %>%
  dplyr::mutate(compound.name = case_when(!is.na(MetabID) ~ MetabID,
                                          TRUE ~ node)) %>%
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
  geom_edge_arc(aes(color = Correlation),
                show.legend = TRUE) +
  geom_node_point(aes(fill = Class,
                      size = Degree),
                  shape = 21,
                  show.legend = TRUE) +
  scale_fill_manual(values = c(
    "Exposome" = ggsci::pal_d3()(10)[2],
    "lab_test" = ggsci::pal_d3()(10)[8]
  )) +
  geom_node_text(
    aes(
      x = x * 1.05,
      y = y * 1.05,
      label = compound.name,
      hjust = 'outward',
      angle = -((-node_angle(x, y) + 90) %% 180) + 90,
      size = 3,
      colour = Class
    ),
    size = 3,
    alpha = 1,
    show.legend = FALSE
  ) +
  guides(
    edge_width = guide_legend(title = "-log10(FDR adjusted P value)",
                              override.aes = list(shape = NA)),
    edge_color = ggraph::guide_edge_colorbar(title = "Spearman correlation"),
    fill = guide_legend(
      title = "Class",
      override.aes = list(size = 4, linetype = "blank")
    ),
    size = guide_legend(title = "Degree", override.aes = list(linetype = 0))
  ) +
  ggraph::scale_edge_color_gradientn(colours = pal) +
  ggraph::scale_edge_width(range = c(0.2, 2)) +
  scale_size_continuous(range = c(5, 10)) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
  )

plot1

ggsave(
  plot1,
  filename = "exposome_lab_test_correlation_network.pdf",
  width = 8.5,
  height = 7,
  bg = "transparent"
)

ggsave(
  plot1,
  filename = "exposome_lab_test_correlation_network.png",
  width = 8.5,
  height = 7,
  bg = "transparent"
)
