##avoid source
no_function()

##load data
sxtTools::setwd_project()
library(tidyverse)
rm(list = ls())

load("data_analysis/exposomeChemical_environment/temp")
load("data_analysis/exposomeChemical_environment/diff")
temp1 = temp
diff1 = diff

load("data_analysis/exposomeBiological_environment/temp")
load("data_analysis/exposomeBiological_environment/diff")
temp2 = temp
diff2 = diff

load("data_analysis/exposomeChemical_exposomeBiological/temp")
load("data_analysis/exposomeChemical_exposomeBiological/diff")
temp3 = temp
diff3 = diff

load("data_analysis/exposomeChemical_environment/cor_value")
exposomeChemical_environment_cor <-
  cor_value %>%
  dplyr::filter(abs(cor) > 0.9 & p_value < 0.05)

load("data_analysis/exposomeBiological_environment/cor_value")
exposomeBiological_environment_cor <-
  cor_value %>%
  dplyr::filter(abs(cor) > 0.9 & p_value < 0.05)

load("data_analysis/exposomeChemical_exposomeBiological/cor_value")
exposomeChemical_exposomeBiological_cor <-
  cor_value %>%
  dplyr::filter(abs(cor) > 0.9 & p_value < 0.05)

head(exposomeBiological_environment_cor)
head(exposomeChemical_environment_cor)
head(exposomeChemical_exposomeBiological_cor)

cor_value <-
  rbind(
    exposomeBiological_environment_cor,
    exposomeChemical_environment_cor,
    exposomeChemical_exposomeBiological_cor
  )

###load variable_info
load("data_20200511/exposome/variable_info")
exposomeChemical_variable_info <-
  variable_info

load("data_20200511/microbiome/dna_variable_info")
exposomeBiological_variable_info <-
  dna_variable_info

load("data_20200511/environment/variable_info")
environment_variable_info <-
  variable_info

sxtTools::setwd_project()
setwd("data_analysis/exposome_cloud/")

library(igraph)
library(ggraph)
library(tidygraph)


###plot show sample matching
temp_edge1 = diff1 %>%
  dplyr::mutate(
    from = paste("exposomeChemical", exposomeChemical.date, sep = "_"),
    to = paste("environment", environment.date, sep = "_")
  ) %>%
  dplyr::select(from, to)

temp_node1 =
  data.frame(
    node = c(
      paste(diff1$environment.class, diff1$environment.date, sep = "_"),
      paste(
        "exposomeChemical",
        diff1$exposomeChemical.date,
        sep = "_"
      )
    ),
    true_name = c(diff1$environment.date,
                  diff1$exposomeChemical.date),
    Class = c(diff1$environment.class, diff1$exposomeChemical.class)
  )

temp_edge2 = diff2 %>%
  dplyr::mutate(
    from = paste("exposomeBiological", exposomeBiological.date, sep = "_"),
    to = paste("environment", environment.date, sep = "_")
  ) %>% 
  dplyr::select(from, to)

temp_node2 =
  data.frame(
    node = c(
      paste("environment", diff1$environment.date, sep = "_"),
      paste(
        "exposomeBiological",
        diff2$exposomeBiological.date,
        sep = "_"
      )
    ),
    true_name = c(diff2$environment.date,
                  diff2$exposomeBiological.date),
    Class = c(diff2$environment.class, diff2$exposomeBiological.class)
  )


temp_edge3 = diff3 %>%
  dplyr::mutate(
    from = paste("exposomeBiological", microbiome.date, sep = "_"),
    to = paste("exposomeChemical", exp.date, sep = "_")
  ) %>% 
  dplyr::select(from, to)

temp_node3 =
  data.frame(
    node = c(
      paste("exposomeChemical", diff3$exp.date, sep = "_"),
      paste(
        "exposomeBiological",
        diff3$microbiome.date,
        sep = "_"
      )
    ),
    true_name = c(diff3$exp.date,
                  diff3$microbiome.date),
    Class = c(diff3$exp.class, diff3$microbiome.class)
  )


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


temp_edge =
  rbind(temp_edge1,
        temp_edge2,
        temp_edge3) %>% 
  dplyr::distinct(.keep_all = TRUE)

temp_node =
  rbind(temp_node1,
        temp_node2,
        temp_node3) %>% 
  dplyr::distinct(.keep_all = TRUE)

temp_node$Class[temp_node$Class == "exp"] = "exposomeChemical"
temp_node$Class[temp_node$Class == "microbiome"] = "exposomeBiological"

temp_node$Class[temp_node$Class == "exposomeChemical"] = "Exposome (chemical)"
temp_node$Class[temp_node$Class == "exposomeBiological"] = "Exposome (biological)"
temp_node$Class[temp_node$Class == "environment"] = "Environment"

temp_node = 
  temp_node %>% 
  dplyr::distinct(.keep_all = TRUE)
  
match_graph <-
  tidygraph::tbl_graph(nodes = temp_node,
                       edges = temp_edge,
                       directed = FALSE) %>%
  dplyr::mutate(Degree = centrality_degree(mode = 'all'))

plot = 
ggraph(match_graph,
       layout = 'linear',
       circular = TRUE) +
  geom_edge_diagonal(show.legend = TRUE) +
  geom_node_point(aes(fill = Class, size = Degree),
                  shape = 21,
                  show.legend = TRUE) +
  scale_fill_manual(values = value) +
  scale_color_manual(values = value) +
  geom_node_text(
    aes(
      x = x * 1.05,
      y = y * 1.05,
      label = true_name,
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
    fill = guide_legend(
      title = "Class",
      override.aes = list(size = 7, linetype = "blank")
    ),
    size = guide_legend(title = "Degree", override.aes = list(linetype = 0))
  ) +
  scale_size_continuous(range = c(6, 10)) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
  )
plot
# ggsave(plot, filename = "match_graph.pdf", width = 19, height = 7)



######exposome cloud
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
      node %in% environment_variable_info$variable_id ~ "Environment"
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
  dplyr::select(node, Class, true_name)


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

node_data$true_name <-
  node_data$true_name %>%
  stringr::str_replace("genus_", "")

###output node data and edge data
library(openxlsx)
wb = createWorkbook()
modifyBaseFont(wb, fontSize = 12, fontName = "Arial Narrow")
addWorksheet(wb, sheetName = "Node information", gridLines = TRUE)
addWorksheet(wb, sheetName = "Edge information", gridLines = TRUE)
freezePane(wb, sheet = 1, firstRow = TRUE, firstCol = TRUE) 
freezePane(wb, sheet = 2, firstRow = TRUE, firstCol = TRUE) 
writeDataTable(wb, sheet = 1, x = node_data,
               colNames = TRUE, rowNames = FALSE)
writeDataTable(wb, sheet = 2, x = edge_data %>% dplyr::select(from, to, everything()),
               colNames = TRUE, rowNames = FALSE)

saveWorkbook(wb, "exposome_cloud.xlsx", overwrite = TRUE)

exposome_cloud <-
  tidygraph::tbl_graph(nodes = node_data,
                       edges = edge_data,
                       directed = TRUE) %>%
  dplyr::mutate(Degree = centrality_degree(mode = 'all'))

pal <-
  wesanderson::wes_palette(name = "Zissou1", n = 100, type = "continuous")

plot <-
  ggraph(exposome_cloud,
         layout = 'linear',
         circular = TRUE) +
  geom_edge_arc(aes(color = Correlation),
                show.legend = TRUE) +
  geom_node_point(aes(fill = Class,
                      size = Degree),
                  shape = 21,
                  show.legend = TRUE) +
  scale_fill_manual(values = value) +
  scale_color_manual(values = value) +
  geom_node_text(
    aes(
      x = x * 1.05,
      y = y * 1.05,
      label = true_name,
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
      override.aes = list(size = 7, linetype = "blank")
    ),
    size = guide_legend(title = "Degree", override.aes = list(linetype = 0))
  ) +
  ggraph::scale_edge_color_gradientn(colours = pal) +
  ggraph::scale_edge_width(range = c(0.2, 2)) +
  scale_size_continuous(range = c(3, 15)) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
  )

plot

# ggsave(
#   plot,
#   filename = "exposome_cloud_network.pdf",
#   width = 8.5,
#   height = 7,
#   bg = "transparent"
# )

edge_data2 <- 
edge_data %>% 
  dplyr::left_join(exposomeChemical_variable_info[,1:2], by = c("from" = "peak_ID")) %>% 
  dplyr::mutate(true_name = 
                  case_when(
                    is.na(MetabID) ~ from,
                    TRUE ~ MetabID
                  )) %>% 
  dplyr::select(true_name, everything()) %>% 
  dplyr::select(-MetabID)

# write.csv(edge_data2, "edge_data.csv", row.names = FALSE)
# write.csv(node_data, "node_data.csv", row.names = FALSE)

# save(exposome_cloud, file = "exposome_cloud")

dim(edge_data)
dim(node_data)

table(node_data$Class)

sxtTools::setwd_project()
load("data_20200511/environment/expression_data")
environment_expression_data <- expression_data

load("data_20200511/environment/sample_info")
environment_sample_info <- sample_info

load("data_20200511/environment/variable_info")
environment_variable_info <- variable_info

head(environment_variable_info)


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

sxtTools::setwd_project()
setwd("data_analysis/exposome_cloud/")

plot <- 
igraph::vertex_attr(graph = exposome_cloud) %>% 
  do.call(cbind, .) %>% 
  as.data.frame() %>% 
  dplyr::filter(Class != "Environment") %>% 
  dplyr::mutate(Degree = as.numeric(Degree)) %>% 
  dplyr::filter(Degree > 1) %>% 
  dplyr::arrange(Class, Degree) %>% 
  dplyr::mutate(true_name = factor(true_name, levels = true_name)) %>% 
  ggplot(aes(y = true_name, x = Degree)) +
  geom_segment(aes(y = true_name, yend = true_name, x = 0, xend = Degree, 
                   color = Class), show.legend = FALSE) +
  geom_point(aes(fill = Class), shape = 21, 
             size = 8,
             show.legend = FALSE) +
  scale_fill_manual(values = value) +
  scale_color_manual(values = value) +
  labs(y = "", x = "Degree") +
  geom_text(aes(y = true_name, x = Degree + 1, 
                label = true_name) 
            # angle = 90, hjust = 1
            ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text = element_text(size = 12),
    axis.text.y = element_blank(),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA)
  ) 
plot

# ggsave(plot, filename = "degree_distributation.pdf", width = 7, height = 7)

##get the subnetwork with degree >= 4
idx <- 
  match(c("Tricholoma", "Cylindrobasidium", "Piriformospora",
          "Diisononyl phthalate", "Butylated triphenyl phosphate"),
        igraph::vertex_attr(graph = exposome_cloud, name = "true_name"))
name <- 
igraph::vertex_attr(graph = exposome_cloud, name = "node")[idx]

name <- 
edge_data %>%
  dplyr::filter(from %in% name | to %in% name) %>% 
  dplyr::select(from, to) 

name <- unique(c(name$from, name$to))

library(igraph)

idx <- match(name, igraph::vertex_attr(graph = exposome_cloud, name = "node"))

subnetwork <-
  igraph::induced_subgraph(graph = exposome_cloud, v = idx)

plot <-
  ggraph(subnetwork,
         layout = 'kk',
         circular = FALSE) +
  geom_edge_link(aes(color = Correlation),
                show.legend = FALSE) +
  geom_node_point(aes(fill = Class,
                      size = Degree),
                  shape = 21,
                  show.legend = FALSE) +
  scale_fill_manual(values = value) +
  scale_color_manual(values = value) +
  geom_node_text(
    aes(
      x = x * 1.05,
      y = y * 1.05,
      label = true_name,
      # hjust = 'outward',
      # angle = -((-node_angle(x, y) + 90) %% 180) + 90,
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
      override.aes = list(size = 7, linetype = "blank")
    ),
    size = guide_legend(title = "Degree", override.aes = list(linetype = 0))
  ) +
  ggraph::scale_edge_color_gradient2(low = ggsci::pal_aaas()(n=10)[1],
                                     mid = "white",
                                     high = ggsci::pal_aaas()(n=10)[2]) +
  ggraph::scale_edge_width(range = c(0.2, 2)) +
  scale_size_continuous(range = c(3, 15)) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
  )

plot

# ggsave(
#   plot,
#   filename = "exposome_cloud_subnetwork.pdf",
#   width = 10,
#   height = 7,
#   bg = "transparent"
# )







