##avoid source
no_function()
rm(list = ls())
sxtTools::setwd_project()
##https://blog.csdn.net/fjsd155/article/details/84726785

library(pls)
data(yarn)
data(mtcars)

# data(gasoline, package = "pls")
# library(pls)
# library(plsVarSel)
# pls <- plsr(octane ~ NIR, ncomp = 10, validation = "LOO", data = gasoline)
# comp <- which.min(pls$validation$PRESS)
# X <- unclass(gasoline$NIR)
# vip <- VIP(pls, comp)
# sr <- SR (pls, comp, X)
# smc <- sMC(pls, comp, X)
# lw <- LW (pls, comp)
# rc <- RC (pls, comp)
# urc <- URC(pls, comp)
# frc <- FRC(pls, comp)
# mrm <- mRMR(pls, 401, X)$score
# matplot(scale(cbind(vip, sr, smc, lw, rc, urc, frc, mrm)), type = 'l')


##load data
#####exposomeChemical
load("data_20200511/exposome/expression_data")
exposomeChemical_expression_data <- expression_data
load("data_20200511/exposome/sample_info")
exposomeChemical_sample_info <- sample_info
load("data_20200511/exposome/variable_info")
exposomeChemical_variable_info <- variable_info

exposomeChemical_expression_data <-
  exposomeChemical_expression_data %>%
  t() %>% 
  as.data.frame() %>%
  tibble::rownames_to_column(var = "Sample_ID") %>%
  dplyr::mutate(Sample_ID = stringr::str_replace(Sample_ID, "_[0-9]{1,2}", "")) %>%
  plyr::dlply(.variables = "Sample_ID") %>%
  lapply(function(x){
    apply(x[,-1], 2, mean)
  }) %>%
  do.call(rbind, .) %>% 
  t() %>% 
  as.data.frame()

##exposomeBiological
###exposome biological
load("data_20200511/microbiome/dna_expression_data")
exposomeBiological_expression_data <- dna_expression_data
load("data_20200511/microbiome/dna_sample_info")
exposomeBiological_sample_info <- dna_sample_info
load("data_20200511/microbiome/dna_variable_info")
exposomeBiological_variable_info <- dna_variable_info

###environment
load("data_20200511/environment/expression_data")
environment_expression_data <- expression_data
load("data_20200511/environment/sample_info")
environment_sample_info <- sample_info
load("data_20200511/environment/variable_info")
environment_variable_info <- variable_info

####load blood_test data
load("data_20200511/lab_test/expression_data")
blood_test_expression_data <- expression_data
load("data_20200511/lab_test/sample_info")
blood_test_sample_info <- sample_info
load("data_20200511/lab_test/variable_info")
blood_test_variable_info <- variable_info

##cytokine
load("data_20200511/cytokine/expression_data")
cytokine_expression_data <- expression_data
load("data_20200511/cytokine/sample_info")
cytokine_sample_info <- sample_info
load("data_20200511/cytokine/variable_info")
cytokine_variable_info <- variable_info

###
setwd("data_analysis/exposome_blood_test_cytokine")

###data preparation
exposomeBiological_sample_info <- 
  exposomeBiological_sample_info %>% 
  dplyr::rename(start_date = date.start)

exposomeChemical_sample_info <- 
  exposomeChemical_sample_info %>% 
  dplyr::rename(start_date = start_date)

environment_sample_info <- 
  environment_sample_info %>% 
  dplyr::rename(start_date = STARTING_DATE)

blood_test_sample_info <- 
  blood_test_sample_info %>% 
  dplyr::rename(start_date = STARTING_DATE)

cytokine_sample_info <- 
  cytokine_sample_info %>% 
  dplyr::rename(start_date = CollectionDate)

exposomeChemical_sample_info$start_date
exposomeBiological_sample_info$start_date
environment_sample_info$start_date

blood_test_sample_info$start_date
cytokine_sample_info$start_date


####exposome vs blood_test
intersect_sample_date <- 
  Reduce(f = intersect, x = list(
    as.character(as.Date(exposomeChemical_sample_info$start_date)),
    as.character(as.Date(exposomeBiological_sample_info$start_date)),
    as.character(as.Date(environment_sample_info$start_date)),
    as.character(as.Date(blood_test_sample_info$start_date))
  ))

exposomeChemical_sample_info <- 
  exposomeChemical_sample_info %>% 
  dplyr::filter(as.character(start_date) %in% intersect_sample_date)

exposomeChemical_expression_data <- 
  exposomeChemical_expression_data %>% 
  dplyr::select(one_of(exposomeChemical_sample_info$sample_id))

exposomeBiological_sample_info <- 
  exposomeBiological_sample_info %>% 
  dplyr::filter(as.character(start_date) %in% intersect_sample_date)

exposomeBiological_expression_data <- 
  exposomeBiological_expression_data %>% 
  dplyr::select(one_of(exposomeBiological_sample_info$sample_id))


environment_sample_info <- 
  environment_sample_info %>% 
  dplyr::filter(as.character(start_date) %in% intersect_sample_date)

environment_expression_data <- 
  environment_expression_data %>% 
  dplyr::select(one_of(environment_sample_info$sample_id))




load("data_analysis/exposomeChemical_cytokine/cor_value")
exposomeChemical_cytokine_cor <-
  cor_value %>%
  dplyr::filter(abs(cor) > 0.9 & p_value < 0.05)
dim(exposomeChemical_cytokine_cor)

load("data_analysis/exposomeChemical_gut_microbiome/cor_value")
exposomeChemical_gut_microbiome_cor <-
  cor_value %>%
  dplyr::filter(abs(cor) > 0.9 & p_value < 0.05)

dim(exposomeChemical_gut_microbiome_cor)

load("data_analysis/exposomeChemical_metabolome/cor_value")
exposomeChemical_metabolome_cor <-
  cor_value %>%
  dplyr::filter(abs(cor) > 0.9 & p_value < 0.05)

dim(exposomeChemical_metabolome_cor)

load("data_analysis/exposomeChemical_proteome/cor_value")
exposomeChemical_proteome_cor <-
  cor_value %>%
  dplyr::filter(abs(cor) > 0.9 & p_value < 0.05)

dim(exposomeChemical_proteome_cor)

#####exposomeBiological
load("data_analysis/exposomeBiological_blood_test/cor_value")
exposomeBiological_blood_test_cor <-
  cor_value %>%
  dplyr::filter(abs(cor) > 0.9 & p_value < 0.05) %>%
  dplyr::mutate(from = paste("exposome", from, sep = "_"))

load("data_analysis/exposomeBiological_cytokine/cor_value")
exposomeBiological_cytokine_cor <-
  cor_value %>%
  dplyr::filter(abs(cor) > 0.9 & p_value < 0.05) %>%
  dplyr::mutate(from = paste("exposome", from, sep = "_"))

load("data_analysis/exposomeBiological_gutmicrobiome/cor_value")
exposomeBiological_gut_microbiome_cor <-
  cor_value %>%
  dplyr::filter(abs(cor) > 0.9 & p_value < 0.05) %>%
  dplyr::mutate(from = paste("exposome", from, sep = "_"))

load("data_analysis/exposomeBiological_metabolome/cor_value")
exposomeBiological_metabolome_cor <-
  cor_value %>%
  dplyr::filter(abs(cor) > 0.9 & p_value < 0.05) %>%
  dplyr::mutate(from = paste("exposome", from, sep = "_"))

load("data_analysis/exposomeBiological_proteome/cor_value")
exposomeBiological_proteome_cor <-
  cor_value %>%
  dplyr::filter(abs(cor) > 0.9 & p_value < 0.05) %>%
  dplyr::mutate(from = paste("exposome", from, sep = "_"))


#####environment
load("data_analysis/environment_blood_test/cor_value")
environment_blood_test_cor <-
  cor_value %>%
  dplyr::filter(abs(cor) > 0.9 & p_value < 0.05)

load("data_analysis/environment_cytokine/cor_value")
environment_cytokine_cor <-
  cor_value %>%
  dplyr::filter(abs(cor) > 0.9 & p_value < 0.05)

load("data_analysis/environment_gutmicrobiome/cor_value")
environment_gut_microbiome_cor <-
  cor_value %>%
  dplyr::filter(abs(cor) > 0.9 & p_value < 0.05)

load("data_analysis/environment_metabolome/cor_value")
environment_metabolome_cor <-
  cor_value %>%
  dplyr::filter(abs(cor) > 0.9 & p_value < 0.05) %>% 
  dplyr::mutate(from1 = to, to1 = from) %>% 
  dplyr::select(-c(from, to)) %>% 
  dplyr::select(from1, to1, everything()) %>% 
  dplyr::rename(from = from1, to = to1)

load("data_analysis/environment_proteome/cor_value")
environment_proteome_cor <-
  cor_value %>%
  dplyr::filter(abs(cor) > 0.9 & p_value < 0.05) %>% 
  dplyr::mutate(from1 = to, to1 = from) %>% 
  dplyr::select(-c(from, to)) %>% 
  dplyr::select(from1, to1, everything()) %>% 
  dplyr::rename(from = from1, to = to1)

cor_value <-
  rbind(
    exposomeChemical_blood_test_cor,
    exposomeChemical_cytokine_cor,
    exposomeChemical_gut_microbiome_cor,
    exposomeChemical_metabolome_cor,
    exposomeChemical_proteome_cor,
    
    exposomeBiological_blood_test_cor,
    exposomeBiological_cytokine_cor,
    exposomeBiological_gut_microbiome_cor,
    exposomeBiological_metabolome_cor,
    exposomeBiological_proteome_cor,
    
    environment_blood_test_cor,
    environment_cytokine_cor,
    environment_gut_microbiome_cor,
    environment_metabolome_cor,
    environment_proteome_cor
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

load("data_20200511/environment/variable_info")
environment_variable_info <-
  variable_info

load("data_20200511/metabolome/variable_info")
metabolome_variable_info <-
  variable_info

load("data_20200511/cytokine/variable_info")
cytokine_variable_info <-
  variable_info

load("data_20200511/gut_microbiome/variable_info")
gut_microbiome_variable_info <-
  variable_info

load("data_20200511/proteome/variable_info")
proteome_variable_info <-
  variable_info

load("data_20200511/lab_test/variable_info")
blood_test_variable_info <-
  variable_info

sxtTools::setwd_project()
setwd("data_analysis/inter_omics_correlation_network/")

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
      node %in% environment_variable_info$variable_id ~ "Environment",
      node %in% blood_test_variable_info$variable_id ~ "Blood test",
      node %in% cytokine_variable_info$variable_id ~ "Cytokine",
      node %in% gut_microbiome_variable_info$variable_id ~ "Gut microbiome",
      node %in% metabolome_variable_info$peak_name ~ "Metabolome",
      node %in% proteome_variable_info$protein_id ~ "Proteome"
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
  dplyr::left_join(metabolome_variable_info[, c("peak_name", "Metabolite")],
                   by = c("node" = "peak_name")) %>%
  dplyr::mutate(true_name = case_when(!is.na(Metabolite) ~ Metabolite,
                                      TRUE ~ true_name)) %>%
  dplyr::select(node, Class, true_name) %>%
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
    "Exposome (chemical)" = 1,
    "Metabolome" = 0.5,
    "Proteome" = 0.5,
    "Exposome (biological)" = 1,
    "Gut microbiome" = 0.5,
    "Blood test" = 0.5,
    "Cytokine" = 0.5,
    "Toxins and carcinogens" = 0.5
  )

node_data$true_name <-
  node_data$true_name %>%
  stringr::str_replace("genus_", "")

total_graph <-
  tidygraph::tbl_graph(nodes = node_data,
                       edges = edge_data,
                       directed = TRUE) %>%
  dplyr::mutate(Degree = centrality_degree(mode = 'all'))

pal <-
  wesanderson::wes_palette(name = "Zissou1", n = 100, type = "continuous")

plot <-
  ggraph(total_graph,
         layout = 'stress') +
  geom_edge_link(aes(color = Correlation),
                 show.legend = TRUE,
                 alpha = 0.5) +
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
  # geom_node_text(
  #   aes(
  #     x = x * 1.05,
  #     y = y * 1.05,
  #     label = true_name,
  #     hjust = 'outward',
  #     angle = -((-node_angle(x, y) + 90) %% 180) + 90,
  #     size = 3,
  #     colour = Class
  #   ),
  #   size = 3,
#   alpha = 1,
#   show.legend = FALSE
# ) +
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

save(total_graph, file = "total_graph")
load("total_graph")

dim(edge_data)
dim(node_data)

ggsave(
  plot,
  filename = "inter_omics_correlation_network.pdf",
  width = 10,
  height = 7,
  bg = "transparent"
)

sum(edge_data$from %in% node_data$node[node_data$Class == "Exposome (chemical)"])
sum(edge_data$from %in% node_data$node[node_data$Class == "Exposome (biological)"])
sum(edge_data$from %in% node_data$node[node_data$Class == "Environment"])

sum(edge_data$to %in% node_data$node[node_data$Class == "Blood test"])
sum(edge_data$to %in% node_data$node[node_data$Class == "Cytokine"])
sum(edge_data$to %in% node_data$node[node_data$Class == "Gut microbiome"])
sum(edge_data$to %in% node_data$node[node_data$Class == "Metabolome"])
sum(edge_data$to %in% node_data$node[node_data$Class == "Proteome"])

write.csv(edge_data, "edge_data.csv", row.names = FALSE)
write.csv(node_data, "node_data.csv", row.names = FALSE)

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
        from %in% environment_variable_info$variable_id ~ "Environment",
        from %in% blood_test_variable_info$variable_id ~ "Blood test",
        from %in% cytokine_variable_info$variable_id ~ "Cytokine",
        from %in% gut_microbiome_variable_info$variable_id ~ "Gut microbiome",
        from %in% metabolome_variable_info$peak_name ~ "Metabolome",
        from %in% proteome_variable_info$protein_id ~ "Proteome"
      )
  ) %>%
  dplyr::mutate(
    var2_class =
      case_when(
        to %in% exposomeChemical_variable_info$peak_ID ~ "Exposome (chemical)",
        to %in% exposomeBiological_variable_info$variable_id ~ "Exposome (biological)",
        to %in% environment_variable_info$variable_id ~ "Environment",
        to %in% blood_test_variable_info$variable_id ~ "Blood test",
        to %in% cytokine_variable_info$variable_id ~ "Cytokine",
        to %in% gut_microbiome_variable_info$variable_id ~ "Gut microbiome",
        to %in% metabolome_variable_info$peak_name ~ "Metabolome",
        to %in% proteome_variable_info$protein_id ~ "Proteome"
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

ggsave(
  plot,
  file = file.path("edge_information.pdf"),
  width = 5,
  height = 7,
  bg = "transparent"
)


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

ggsave(
  plot,
  file = file.path("node_information.pdf"),
  width = 5,
  height = 7,
  bg = "transparent"
)

total_graph


###community analysis
library(igraph)
all_subnetworks <-
  igraph::cluster_walktrap(graph = total_graph, steps = 4)

save(all_subnetworks, file = "all_subnetworks")

load("all_subnetworks")

plot <- 
  ggplot(
    data.frame(index = 1:length(all_subnetworks$modularity),
               modu = all_subnetworks$modularity, stringsAsFactors = FALSE),
    aes(index, modu) 
  ) +
  geom_vline(xintercept = which.max(all_subnetworks$modularity), 
             linetype = 2, colour = "#800000B2") + 
  labs(x = "Community analysis iteration", y = "Modularity") +
  geom_line(colour = "black") +
  # geom_point() +
  theme_bw() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12))

plot <-
  plot + 
  ggplot2::annotate(geom = "point", 
                    x = which.max(all_subnetworks$modularity),
                    y = max(all_subnetworks$modularity), 
                    size = 3, 
                    colour = "#FFA319FF") +
  annotate(geom = "text", 
           x = which.max(all_subnetworks$modularity),
           y = max(all_subnetworks$modularity), 
           label = paste("(",  which.max(all_subnetworks$modularity),
                         ",", 
                         max(all_subnetworks$modularity) %>% round(3),
                         ")"),
           size = 5,
           colour = "#FFA319FF"
  )

plot

ggsave(plot, filename = "all_modularity.pdf", width = 7, height = 7)
ggsave(plot, filename = "all_modularity.png", width = 7, height = 7)

all_subnetworks$membership

total_node_info <- 
igraph::vertex_attr(total_graph) %>% 
  purrr::map(function(x){
    as.character(x)
  }) %>% 
  do.call(cbind, .) %>% 
  as.data.frame() %>% 
  data.frame(membership = all_subnetworks$membership) %>% 
  dplyr::mutate(Degree = as.numeric(Degree))

table(total_node_info$membership)

which(table(total_node_info$membership) >= 5)

idx <- 2

for(idx in unique(total_node_info$membership)) {
  cat(idx, " ")
  
  subnetwork <-
    igraph::induced_subgraph(graph = total_graph,
                             v = which(total_node_info$membership == idx))
  
  path <- paste("subnetwork", idx, sep = "_")
  dir.create(path)
  
  pal <-
    wesanderson::wes_palette(name = "Zissou1", n = 100, type = "continuous")
  
  plot <-
    ggraph(subnetwork,
           layout = 'linear',
           circular = TRUE) +
    geom_edge_arc(
      strength = 1,
      aes(color = Correlation),
      alpha = 1,
      show.legend = FALSE,
      edge_width = 0.5
      # arrow = arrow(length = unit(2, 'mm'))
    ) +
    geom_node_point(
      aes(fill = Class,
          size = Degree),
      alpha = 1,
      show.legend = FALSE,
      shape = 21
    ) +
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
    guides(color = guide_legend(override.aes = list(size = 3))) +
    ggraph::scale_edge_color_gradient2(low = pal[1],
                                       mid = "white",
                                       high = pal[100]) +
    scale_size_continuous(range = c(3, 10)) +
    scale_color_manual(values = value) +
    scale_fill_manual(values = value) +
    ggraph::theme_graph() +
    theme(
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA),
      legend.position = "right",
      legend.background = element_rect(fill = "transparent", color = NA)
    )
  
  dir.create(paste("subnetwork", idx, sep = "_"))
  
  ggsave(
    plot,
    filename = file.path(path, paste("subnetwork", idx, "_circle.pdf", sep = "")),
    width = 7,
    height = 7,
    bg = "transparent"
  )
}



###output each subnetwork
idx <- 2
idx <- which(table(membership(communities = all_subnetworks)) >= 5)
subnetwork_info <- vector(mode = "list", length = max(idx))

for(i in idx){
  cat(i, " ")
  
  subnetwork <-
    subgraph(graph = all_omics_graph,
             v = which(membership(all_subnetworks) == i))
  
  temp_node <- 
    igraph::vertex_attr(subnetwork)
  
  temp_node <- lapply(temp_node, as.character) %>% 
    do.call(cbind, .) %>% 
    as.data.frame()
  temp_edge <- igraph::as_data_frame(subnetwork)
  
  temp_edge$from <- temp_node$node[temp_edge$from]
  temp_edge$to <- temp_node$node[temp_edge$to]
  
  dir.create(paste("subnetwork", i, sep = "_"))
  write.csv(temp_node, 
            file = file.path(paste("subnetwork", i, sep = "_"), 
                             paste("subnetwork", i, "node.csv", sep = "_")), 
            row.names = FALSE)
  write.csv(temp_edge,
            file = file.path(
              paste("subnetwork", i, sep = "_"),
              paste("subnetwork", i, "edge.csv", sep = "_")
            ),
            row.names = FALSE)
  
  node_num <- dim(temp_node)[1]
  edge_num <- dim(temp_edge)[1]
  node_num2 <- table(temp_node$class)
  num <- c(node_num, edge_num, node_num2)
  names(num)[1:2] <- c("node", "edge")
  num2 <- rep(0, 9)
  names(num2) <- c("node", "edge", "RNA", "Protein",
                   "Metabolite", "Cytokine", "Virginal microbiome", 
                   "Gut microbiome", 
                   "Clinical")
  
  num2[match(names(num), names(num2))] <- num
  
  subnetwork_info[[i]] <- c(subnetwork = i, num2)
}


subnetwork_info <- 
  subnetwork_info %>% 
  do.call(rbind, .) %>% 
  as.data.frame()

write.csv(subnetwork_info, "subnetwork_info.csv", row.names = FALSE)


plot <-
  subnetwork_info %>%
  ggplot(aes(node, edge)) +
  geom_point(
    fill = "black",
    shape = 21,
    aes(size = node),
    show.legend = FALSE,
    alpha = 0.8
  ) +
  ggrepel::geom_text_repel(aes(node, edge,
                               label = ifelse(node > 100, subnetwork, NA)),
                           hjust = 2.5) +
  ggforce::facet_zoom(
    xlim = c(0, 70),
    ylim = c(0, 80),
    horizontal = FALSE,
    zoom.size = 1
  ) +
  theme_bw() +
  labs(x = "Node number", y = "Edge number") +
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right",
    legend.background = element_rect(fill = "transparent", color = NA)
  ) +
  # cowplot::theme_cowplot() +
  theme(legend.position = "top")

plot

ggsave(
  plot,
  filename = "subnetwork_plot.pdf",
  width = 7,
  height = 7,
  bg = "transparent"
)

ggsave(
  plot,
  filename = "subnetwork_plot.png",
  width = 7,
  height = 7,
  bg = "transparent"
)

