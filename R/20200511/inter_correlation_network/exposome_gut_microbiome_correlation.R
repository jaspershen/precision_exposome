##avoid source
no_function()

##load data
sxtTools::setwd_project()
library(tidyverse)
rm(list = ls())

#####exposomeChemical
load("data_analysis/exposomeChemical_gut_microbiome/cor_value")
exposomeChemical_gut_microbiome_cor <-
  cor_value %>%
  dplyr::filter(abs(cor) > 0.9 & p_value < 0.05)

dim(exposomeChemical_gut_microbiome_cor)

load("data_analysis/exposomeChemical_gut_microbiome/exp_expression_data1")
exposomeChemical_expression_data <- exp_expression_data1
load("data_analysis/exposomeChemical_gut_microbiome/exp_sample_info")
exposomeChemical_sample_info <- exp_sample_info

load("data_analysis/exposomeChemical_gut_microbiome/gut_microbiome_expression_data1")
gut_microbiome_expression_data1 <- gut_microbiome_expression_data1
load("data_analysis/exposomeChemical_gut_microbiome/gut_microbiome_sample_info")
gut_microbiome_sample_info1 <- gut_microbiome_sample_info

#####exposomeBiological
load("data_analysis/exposomeBiological_gutmicrobiome/cor_value")
exposomeBiological_gut_microbiome_cor <-
  cor_value %>%
  dplyr::filter(abs(cor) > 0.9 & p_value < 0.05) %>%
  dplyr::mutate(from = paste("exposome", from, sep = "_"))

load("data_analysis/exposomeBiological_gutmicrobiome/exposomeBiological_expression_data1")
exposomeBiological_expression_data <- exposomeBiological_expression_data1
load("data_analysis/exposomeBiological_gutmicrobiome/exposomeBiological_sample_info")
exposomeBiological_sample_info <- exposomeBiological_sample_info

load("data_analysis/exposomeBiological_gutmicrobiome/gutmicrobiome_expression_data1")
gut_microbiome_expression_data2 <- gutmicrobiome_expression_data1
load("data_analysis/exposomeBiological_gutmicrobiome/gutmicrobiome_sample_info")
gut_microbiome_sample_info2 <- gutmicrobiome_sample_info

exposomeChemical_sample_info$start_date
gut_microbiome_sample_info1$CollectionDate

exposomeBiological_sample_info$date.start
gut_microbiome_sample_info2$CollectionDate

cor_value <-
  rbind(
    exposomeChemical_gut_microbiome_cor,
    exposomeBiological_gut_microbiome_cor
  )

dim(cor_value)

###load variable_info
load("data_20200511/exposome/variable_info")
exposomeChemical_variable_info <-
  variable_info

load("data_20200511/microbiome/dna_variable_info")
exposomeBiological_variable_info <-
  dna_variable_info %>%
  dplyr::mutate(variable_id = paste("exposome", variable_id, sep = "_"))

load("data_20200511/gut_microbiome/variable_info")
gut_microbiome_variable_info <-
  variable_info

sxtTools::setwd_project()
setwd("data_analysis/exposome_gut_microbiome_correlation_network/")

library(igraph)
library(ggraph)
library(tidygraph)

edge_data <-
  cor_value %>%
  dplyr::rename(from = from,
                to = to,
                Correlation = cor) %>%
  dplyr::mutate(fdr = -log(fdr, 10))

node_data <-
  cor_value %>%
  dplyr::select(from, to) %>%
  tidyr::pivot_longer(cols = c(from, to),
                      names_to = "class",
                      values_to = "node") %>%
  dplyr::mutate(
    class1 = case_when(
      node %in% exposomeChemical_variable_info$peak_ID ~ "Exposome (chemical)",
      node %in% exposomeBiological_variable_info$variable_id ~ "Exposome (biological)",
      node %in% gut_microbiome_variable_info$variable_id ~ "Gut microbiome"
    )
  ) %>%
  dplyr::select(node, class1) %>%
  dplyr::rename(Class = class1) %>%
  dplyr::distinct(node, .keep_all = TRUE) %>%
  dplyr::arrange(Class) %>%
  dplyr::left_join(exposomeChemical_variable_info[, c("peak_ID", "MetabID")],
                   by = c("node" = "peak_ID")) %>%
  dplyr::mutate(true_name = case_when(!is.na(MetabID) ~ MetabID,
                                      TRUE ~ node)) %>%
  dplyr::select(node, Class, true_name) %>%
  dplyr::mutate(Class = factor(
    Class,
    levels = c(
      "Exposome (chemical)",
      "Exposome (biological)",
      "Gut microbiome"
    )
  ))

node_data$true_name <-
  node_data$true_name %>% 
  stringr::str_replace("exposome_", "")

value <-
  c(
    "Environment" = ggsci::pal_d3()(10)[1],
    "Exposome (chemical)" = ggsci::pal_d3()(10)[2],
    "Metabolome" = ggsci::pal_d3()(10)[3],
    "Proteome" = ggsci::pal_d3()(10)[4],
    "Exposome (biological)" = ggsci::pal_d3()(10)[5],
    "Gut microbiome" = ggsci::pal_d3()(10)[6],
    "Blood test" = ggsci::pal_d3()(10)[7],
    "Cytokine" = ggsci::pal_d3()(10)[8],
    "Toxins and carcinogens" = ggsci::pal_d3()(10)[9]
  )

alpha_value <-
  c(
    "Environment" = 1,
    "Exposome (chemical)" = 0.5,
    "Metabolome" = 0.5,
    "Proteome" = 0.5,
    "Exposome (biological)" = 0.5,
    "Gut microbiome" = 0.5,
    "Blood test" = 0.5,
    "Cytokine" = 0.5,
    "Toxins and carcinogens" = 0.5
  )

# node_data$true_name <-
#   node_data$true_name %>%
#   stringr::str_replace("genus_", "")

edge_info <- edge_data

edge_info <- 
edge_info %>% 
  dplyr::left_join(node_data[,c(1,3)], by = c("from" = "node")) %>% 
  dplyr::rename(from_name = true_name) %>% 
  dplyr::left_join(node_data[,c(1,3)], by = c("to" = "node")) %>% 
  dplyr::rename(to_name = true_name)

edge_info$from_number <-
  purrr::map(edge_info$from, function(y){
    sum(y == edge_info$from)
  }) %>% 
  unlist()

edge_info$to_number <-
  purrr::map(edge_info$to, function(y){
    sum(y == edge_info$to)
  }) %>% 
  unlist()

# write.csv(edge_info, "edge_info.csv", row.names = FALSE)

total_graph <-
  tidygraph::tbl_graph(nodes = node_data,
                       edges = edge_data,
                       directed = TRUE) %>%
  dplyr::mutate(Degree = centrality_degree(mode = 'all'))

pal <-
  wesanderson::wes_palette(name = "Zissou1", n = 100, type = "continuous")

library(igraph)

g <- total_graph

V(g)$type <- bipartite_mapping(g)$type

coords <- 
  create_layout(g, layout = "bipartite") %>%
  dplyr::select(x, y) 

# coords$index <- 1:nrow(coords)
# coords <-
#   coords %>% 
#   dplyr::arrange(x)
# 
# coords$x[coords$y == 0] <- 
#   seq(0, 332, length.out = sum(coords$y == 0))
# 
# coords$x[coords$y == 1] <- 
#   seq(0, 332, length.out = sum(coords$y == 1))

coords$y[coords$y == 0] <- 0.3

coords <- 
coords %>% 
  dplyr::select(x,y) %>% 
  dplyr::mutate(theta = x / (max(x) + 1) * 2 * pi,
         r = y + 1, 
         x = r * cos(theta), 
         y = r * sin(theta))

my_graph <-
  create_layout(graph = g,
                layout = "manual", 
                x = coords$x,
                y = coords$y
                # node.position = coords
                )

plot <-
  ggraph(my_graph,
         layout = 'bipartite') +
  geom_edge_link(aes(color = Correlation),
                 show.legend = TRUE) +
  geom_node_point(
    aes(fill = Class,
        size = Degree,
        alpha = Class),
    shape = 21,
    # alpha = 0.6,
    show.legend = TRUE
  ) +
  scale_fill_manual(values = value) +
  scale_color_manual(values = value) +
  scale_alpha_manual(values = alpha_value) +
  geom_node_text(
    aes(
      x = x * 1.03,
      y = y * 1.03,
      label = ifelse(Class == "Gut microbiome", true_name, true_name),
      hjust = ifelse(Class == "Gut microbiome", "inward",'outward'),
      angle = -((-node_angle(x, y) + 90) %% 180) + 90,
      size = ifelse(Class == "Gut microbiome", 2, 4),
      # size = 1,
      colour = Class
    ),repel = FALSE,
    # size = 3,
    alpha = 1,
    show.legend = FALSE
  ) +
  guides(
    edge_width = guide_legend(title = "-log10(FDR adjusted P value)",
                              override.aes = list(shape = NA)),
    edge_color = ggraph::guide_edge_colorbar(title = "Spearman correlation"),
    fill = guide_legend(
      title = "Class",
      override.aes = list(size = 7, linetype = "blank")
    ),
    size = guide_legend(title = "Degree", override.aes = list(linetype = 0))
  ) +
  ggraph::scale_edge_color_gradientn(colours = pal) +
  ggraph::scale_edge_width(range = c(0.2, 2)) +
  scale_size_continuous(range = c(0.5, 5)) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
  )

plot

# write.csv(edge_data, "edge_data.csv", row.names = FALSE)
# save(total_graph, file = "total_graph")
load("total_graph")

dim(edge_data)
dim(node_data)

# ggsave(
#   plot,
#   filename = "inter_omics_correlation_network.pdf",
#   width = 8.5,
#   height = 7,
#   bg = "transparent"
# )

sum(edge_data$from %in% node_data$node[node_data$Class == "Exposome (chemical)"])
sum(edge_data$from %in% node_data$node[node_data$Class == "Exposome (biological)"])
sum(edge_data$from %in% node_data$node[node_data$Class == "Environment"])

sum(edge_data$to %in% node_data$node[node_data$Class == "Blood test"])
sum(edge_data$to %in% node_data$node[node_data$Class == "Cytokine"])
sum(edge_data$to %in% node_data$node[node_data$Class == "Gut microbiome"])
sum(edge_data$to %in% node_data$node[node_data$Class == "Metabolome"])
sum(edge_data$to %in% node_data$node[node_data$Class == "Proteome"])

# write.csv(edge_data, "edge_data.csv", row.names = FALSE)
# write.csv(node_data, "node_data.csv", row.names = FALSE)

dim(edge_data)
dim(node_data)

table(node_data$Class)

###edge distributation
temp_data <-
  cor_value %>%
  dplyr::select(from, to) %>%
  dplyr::mutate(
    var1_class =
      case_when(
        from %in% exposomeChemical_variable_info$peak_ID ~ "Exposome (chemical)",
        from %in% exposomeBiological_variable_info$variable_id ~ "Exposome (biological)",
        # from %in% environment_variable_info$variable_id ~ "Environment",
        # from %in% blood_test_variable_info$variable_id ~ "Blood test",
        # from %in% cytokine_variable_info$variable_id ~ "Cytokine",
        # from %in% gut_microbiome_variable_info$variable_id ~ "Gut microbiome",
        # from %in% metabolome_variable_info$peak_name ~ "Metabolome",
        from %in% gut_microbiome_variable_info$variable_id ~ "Gut microbiome"
      )
  ) %>%
  dplyr::mutate(
    var2_class =
      case_when(
        to %in% exposomeChemical_variable_info$peak_ID ~ "Exposome (chemical)",
        to %in% exposomeBiological_variable_info$variable_id ~ "Exposome (biological)",
        # to %in% environment_variable_info$variable_id ~ "Environment",
        # to %in% blood_test_variable_info$variable_id ~ "Blood test",
        # to %in% cytokine_variable_info$variable_id ~ "Cytokine",
        # to %in% gut_microbiome_variable_info$variable_id ~ "Gut microbiome",
        # to %in% metabolome_variable_info$peak_name ~ "Metabolome",
        to %in% gut_microbiome_variable_info$variable_id ~ "Gut microbiome"
      )
  )

library(ggalluvial)
temp_data1 <- temp_data %>%
  dplyr::mutate(id = paste(from, to, var1_class, var2_class, sep = "_"))

# temp_data1 <-
#   temp_data %>%
#   t() %>%
#   as.data.frame() %>%
#   purrr::map(
#     .f = function(x) {
#       if (x[4] == "Environment") {
#         x <- x[c(2, 1, 4, 3)]
#       }
#       as.character(x)
#     }
#   ) %>%
#   do.call(rbind, .) %>%
#   as.data.frame() %>%
#   dplyr::mutate(id = paste(V1, V2, V3, V4, sep = "_"))
# 
# colnames(temp_data1)[1:4] <-
#   colnames(temp_data)
# 
# rownames(temp_data1) <- NULL

temp_data2 <-
  temp_data1 %>%
  dplyr::select(id, var1_class, var2_class) %>%
  tidyr::pivot_longer(cols = -id,
                      names_to = "class",
                      values_to = "data")  %>%
  dplyr::mutate(freq = 1) %>%
  dplyr::mutate(data = factor(
    data,
    levels = c(
      "Exposome (chemical)",
      "Exposome (biological)",
      "Environment",
      "Blood test",
      "Cytokine",
      "Gut microbiome",
      "Metabolome",
      "Proteome"
    )
  ))

value <-
  c(
    "Environment" = ggsci::pal_d3()(10)[1],
    "Exposome (chemical)" = ggsci::pal_d3()(10)[2],
    "Metabolome" = ggsci::pal_d3()(10)[3],
    "Proteome" = ggsci::pal_d3()(10)[4],
    "Exposome (biological)" = ggsci::pal_d3()(10)[5],
    "Gut microbiome" = ggsci::pal_d3()(10)[6],
    "Blood test" = ggsci::pal_d3()(10)[7],
    "Cytokine" = ggsci::pal_d3()(10)[8],
    "Toxins and carcinogens" = ggsci::pal_d3()(10)[9]
  )

plot <-
  ggplot(temp_data2,
         aes(
           x = class,
           y = freq,
           stratum = data,
           alluvium = id,
           fill = data,
           label = data
         )) +
  scale_x_discrete(expand = c(.1, .1)) +
  ggalluvial::geom_flow(show.legend = FALSE) +
  labs(x = "", y = "") +
  scale_fill_manual(values = value) +
  ggalluvial::geom_stratum(alpha = 1,
                           color = "black",
                           show.legend = FALSE) +
  # geom_text(stat = "stratum", size = 3) +
  theme_bw() +
  theme(
    legend.position = "top",
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  )

plot

# ggsave(
#   plot,
#   file = file.path("edge_information.pdf"),
#   width = 5,
#   height = 7,
#   bg = "transparent"
# )


####node distributation
plot <-
  node_data %>%
  dplyr::group_by(Class) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Class = factor(
    Class,
    levels = c(
      "Exposome (chemical)",
      "Exposome (biological)",
      "Environment",
      "Blood test",
      "Cytokine",
      "Gut microbiome",
      "Metabolome",
      "Proteome"
    ) %>%
      rev()
  )) %>%
  ggplot(aes(x = 2, y = Class)) +
  geom_point(shape = 21,
             aes(fill = Class, size = n),
             show.legend = FALSE) +
  theme_void() +
  scale_size_continuous(range = c(5, 30)) +
  scale_fill_manual(values = value) +
  geom_text(aes(x = 2, y = Class, label = n))

plot

# ggsave(
#   plot,
#   file = file.path("node_information.pdf"),
#   width = 5,
#   height = 7,
#   bg = "transparent"
# )






####example graph in figure
load("total_graph")
edge_data1 = 
  igraph::as_data_frame(total_graph)

node_data1 = 
  igraph::vertex_attr(total_graph) %>% 
  purrr::map(function(x){
    as.character(x)
  }) %>% 
  do.call(cbind, .) %>% 
  as.data.frame() %>% 
  dplyr::mutate(Degree = as.numeric(Degree))

edge_data1$from = 
  node_data1$node[edge_data1$from]

edge_data1$to = 
  node_data1$node[edge_data1$to]

node_data1 = 
  node_data1 %>% 
  dplyr::filter((Class != "Gut microbiome") | (
    Class == "Gut microbiome" & Degree > 20
  ))

edge_data1 = 
edge_data1 %>%
  dplyr::filter((from %in% node_data1$node[node_data1$Class == "Gut microbiome"]) |
                  (to %in% node_data1$node[node_data1$Class == "Gut microbiome"]))

node_data1 =
  node_data1 %>% 
  dplyr::filter(node %in% unique(c(edge_data1$from, edge_data1$to)))


total_graph1 <-
  tidygraph::tbl_graph(nodes = node_data1,
                       edges = edge_data1,
                       directed = TRUE)

pal <-
  wesanderson::wes_palette(name = "Zissou1", n = 100, type = "continuous")

library(igraph)

g <- total_graph1

V(g)$type <- bipartite_mapping(g)$type

coords <- 
  create_layout(g, layout = "bipartite") 

coords$index = 1:nrow(coords)

coords$y[coords$y == 0] <- 0.5

# coords1 = coords[coords$y == -0.1,] %>% 
#   dplyr::arrange(x)
# coords2 = coords[coords$y != -0.1,] %>% 
#   dplyr::arrange(x)
# 
# coords1$x = 
#   seq(
#     from = min(coords1$x),
#     to = max(coords1$x),
#     length.out = nrow(coords1)
#   )
# 
# coords = 
#   rbind(coords1,
#         coords2) %>% 
#   as.data.frame() %>% 
#   dplyr::arrange(index)


coords <- 
  coords %>% 
  dplyr::select(x,y) %>% 
  dplyr::mutate(theta = x / (max(x) + 1) * 2 * pi,
                r = y + 1, 
                x = r * cos(theta), 
                y = r * sin(theta))

my_graph <-
  create_layout(graph = g,
                layout = "manual", 
                x = coords$x,
                y = coords$y
                # node.position = coords
  )

plot <-
  ggraph(my_graph,
         layout = 'bipartite') +
  geom_edge_diagonal(aes(color = Correlation),
                 show.legend = TRUE) +
  geom_node_point(
    aes(fill = Class,
        size = Degree,
        alpha = Class),
    shape = 21,
    alpha = 1,
    show.legend = TRUE
  ) +
  scale_fill_manual(values = value) +
  scale_color_manual(values = value) +
  scale_alpha_manual(values = alpha_value) +
  geom_node_text(
    aes(
      x = x * 1.03,
      y = y * 1.03,
      label = ifelse(Class == "Gut microbiome", true_name, true_name),
      hjust = ifelse(Class == "Gut microbiome", "inward",'outward'),
      angle = -((-node_angle(x, y) + 90) %% 180) + 90,
      # size = ifelse(Class == "Gut microbiome", 4, 4),
      size = 10,
      colour = Class
    ),repel = FALSE,
    # size = 3,
    alpha = 1,
    show.legend = FALSE
  ) +
  guides(
    edge_width = guide_legend(title = "-log10(FDR adjusted P value)",
                              override.aes = list(shape = NA)),
    edge_color = ggraph::guide_edge_colorbar(title = "Spearman correlation"),
    fill = guide_legend(
      title = "Class",
      override.aes = list(size = 7, linetype = "blank")
    ),
    size = guide_legend(title = "Degree", override.aes = list(linetype = 0))
  ) +
  ggraph::scale_edge_color_gradientn(colours = pal) +
  ggraph::scale_edge_width(range = c(0.2, 2)) +
  scale_size_continuous(range = c(0.5, 5)) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
  )

plot

# ggsave(
#   plot,
#   filename = "inter_omics_correlation_network_example.pdf",
#   width = 8.5,
#   height = 7,
#   bg = "transparent"
# )



temp_name <- "Roseburia"
grep(temp_name,edge_data$from_name, value = TRUE) %>% length()

idx1 <- 
grep(temp_name,edge_data$to, value = FALSE)
idx1

edge_info$from[idx1] %>% unique()
edge_info$to[idx1] %>% unique()

temp_data1 <-
  gut_microbiome_expression_data2[match(edge_info$to[idx1] %>% unique(),
                                       rownames(gut_microbiome_expression_data2)), ]

temp_data2 <-
  exposomeChemical_expression_data[match(edge_info$from[idx1] %>% unique(),
                                         rownames(exposomeChemical_expression_data)), ] %>% 
  tibble::rownames_to_column(var = "node") %>% 
  dplyr::filter(!is.na(`2016-01-12`)) %>% 
  dplyr::left_join(node_data[,c(1,3)], by = "node") %>% 
  dplyr::select(-node) %>% 
  tibble::column_to_rownames(var = "true_name")

temp_data3 <-
  exposomeBiological_expression_data[match(
    edge_info$from[idx1] %>% unique() %>%
      stringr::str_replace_all("exposome_", ""),
    rownames(exposomeBiological_expression_data)
  ),] %>%
  tibble::rownames_to_column(var = "node") %>%
  dplyr::filter(!is.na(`2016-01-12`)) %>%
  tibble::column_to_rownames(var = "node")
  
temp_data <- 
  rbind(
    temp_data1,
    temp_data2,
    temp_data3
  )

library(ComplexHeatmap)

temp_data <- 
  apply(temp_data, 1, function(x){
    (x - mean(x)) / sd(x)
  }) %>% 
  t() %>% 
  as.data.frame()
library(circlize)
col_fun = colorRamp2(
  breaks = seq(min(temp_data), max(temp_data), length.out = 90),
  colors =
    viridis::magma(n = 100)[-c(1:10)],
  transparency = 0
)

plot1 <-
  temp_data %>%
  ComplexHeatmap::Heatmap(
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    show_column_names = TRUE,
    show_row_names = TRUE,
    clustering_method_rows = "ward.D",
    clustering_method_columns = "ward.D",
    clustering_distance_columns = "euclidean",
    clustering_distance_rows = "euclidean",
    col = col_fun,
    rect_gp = gpar(col = "white"),
    border = TRUE,
    row_dend_reorder = TRUE,
    column_dend_reorder = TRUE,
    column_names_rot = 45,
    name = "Z-score"
  )

plot1 <- ggplotify::as.ggplot(plot1)
plot1

temp_edge_data <- 
edge_info[idx1,] %>% 
  dplyr::select(from, to, Correlation) 

temp_node_data <-
  node_data[match(rev(rownames(temp_data)), node_data$true_name), ] %>%
  dplyr::mutate(node = factor(node, levels = node))

value <-
  c(
    "Environment" = ggsci::pal_d3()(10)[1],
    "Exposome (chemical)" = ggsci::pal_d3()(10)[2],
    "Metabolome" = ggsci::pal_d3()(10)[3],
    "Proteome" = ggsci::pal_d3()(10)[4],
    "Exposome (biological)" = ggsci::pal_d3()(10)[5],
    "Gut microbiome" = ggsci::pal_d3()(10)[6],
    "Blood test" = ggsci::pal_d3()(10)[7],
    "Cytokine" = ggsci::pal_d3()(10)[8],
    "Toxins and carcinogens" = ggsci::pal_d3()(10)[9]
  )

alpha_value <-
  c(
    "Environment" = 1,
    "Exposome (chemical)" = 0.5,
    "Metabolome" = 0.5,
    "Proteome" = 0.5,
    "Exposome (biological)" = 0.5,
    "Gut microbiome" = 0.5,
    "Blood test" = 0.5,
    "Cytokine" = 0.5,
    "Toxins and carcinogens" = 0.5
  )

temp_total_graph <-
  tidygraph::tbl_graph(nodes = temp_node_data,
                       edges = temp_edge_data,
                       directed = TRUE) %>%
  dplyr::mutate(Degree = centrality_degree(mode = 'all'))

pal <-
  wesanderson::wes_palette(name = "Zissou1", n = 100, type = "continuous")

library(igraph)

plot2 <-
  ggraph(temp_total_graph,
         layout = 'linear', circular = FALSE) +
  geom_edge_arc(
    aes(color = Correlation),
    show.legend = TRUE,
    arrow = grid::arrow(),
    strength = -0.2
  ) +
  geom_node_point(
    aes(fill = Class),
    shape = 21,
    size = 4,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = value) +
  scale_color_manual(values = value) +
  # geom_node_text(
  #   aes(
  #     x = x * 1.03,
  #     y = y * 1.03,
  #     label = true_name,
  #     hjust = "outward",
  #     angle = 0,
  #     size = 3,
  #     colour = Class
  #   ),repel = FALSE,
  #   # size = 3,
  #   alpha = 1,
  #   show.legend = FALSE
  # ) +
  guides(
    # edge_width = guide_legend(title = "-log10(FDR adjusted P value)",
    #                           override.aes = list(shape = NA)),
    edge_color = ggraph::guide_edge_colorbar(title = "Spearman correlation")
    # fill = guide_legend(title = "Class",
    #                     override.aes = list(size = 7, linetype = "blank"))
  ) +
  ggraph::scale_edge_color_gradientn(colours = pal) +
  ggraph::scale_edge_width(range = c(0.2, 2)) +
  scale_size_continuous(range = c(0.5, 5)) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "left"
  ) +
  coord_flip()

plot2

library(patchwork)

plot <- 
plot2 + plot1 +
  patchwork::plot_layout(ncol = 2, widths = c(1,5))

plot

# ggsave(plot, filename = paste(temp_name, "heatmap.pdf"), width = 9, height = 6)

