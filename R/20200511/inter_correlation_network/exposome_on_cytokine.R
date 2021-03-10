##avoid source
no_function()
rm(list = ls())
sxtTools::setwd_project()
##https://blog.csdn.net/fjsd155/article/details/84726785

library(pls)
data(yarn)
data(mtcars)

#####exposomeChemical
load("data_analysis/exposomeChemical_cytokine/cor_value")
exposomeChemical_cytokine_cor <-
  cor_value %>%
  dplyr::filter(abs(cor) > 0.9 & p_value < 0.05)

dim(exposomeChemical_cytokine_cor)

#####exposomeBiological
load("data_analysis/exposomeBiological_cytokine/cor_value")
exposomeBiological_cytokine_cor <-
  cor_value %>%
  dplyr::filter(abs(cor) > 0.9 & p_value < 0.05) %>%
  dplyr::mutate(from = paste("exposome", from, sep = "_"))

#####environment
load("data_analysis/environment_cytokine/cor_value")
environment_cytokine_cor <-
  cor_value %>%
  dplyr::filter(abs(cor) > 0.9 & p_value < 0.05) %>%
  dplyr::mutate(from = paste("exposome", from, sep = "_"))

cor_value <-
  rbind(
    exposomeChemical_cytokine_cor,
    exposomeBiological_cytokine_cor,
    environment_cytokine_cor
  )

dim(cor_value)

##load data
#####exposomeChemical
load("data_20200511/exposome/expression_data")
exposomeChemical_expression_data <- expression_data
load("data_20200511/exposome/sample_info")
exposomeChemical_sample_info <- sample_info
load("data_20200511/exposome/variable_info")
exposomeChemical_variable_info <- variable_info

###cytokine class
cytokine_class <- readxl::read_xlsx("data_analysis/cytokine_analysis/cy_cla.xlsx")

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
exposomeChemical_variable_info <- 
  exposomeChemical_variable_info %>% 
  dplyr::filter(peak_ID %in% unique(c(cor_value$from,cor_value$to)))

exposomeChemical_expression_data <-
  exposomeChemical_expression_data[match(
    exposomeChemical_variable_info$peak_ID,
    rownames(exposomeChemical_expression_data)
  ),]

##exposomeBiological
###exposome biological
load("data_20200511/microbiome/dna_expression_data")
exposomeBiological_expression_data <- dna_expression_data
load("data_20200511/microbiome/dna_sample_info")
exposomeBiological_sample_info <- dna_sample_info
load("data_20200511/microbiome/dna_variable_info")
exposomeBiological_variable_info <- dna_variable_info

exposomeBiological_variable_info <- 
exposomeBiological_variable_info %>% 
  dplyr::filter(level == "genus")

exposomeBiological_expression_data <- 
  exposomeBiological_expression_data[match(exposomeBiological_variable_info$variable_id, 
                                           rownames(exposomeBiological_expression_data)),]

exposomeBiological_variable_info <-
  exposomeBiological_variable_info %>%
  dplyr::filter(variable_id %in% stringr::str_replace(unique(c(cor_value$from, cor_value$to)),"exposome_", ""))

exposomeBiological_expression_data <-
  exposomeBiological_expression_data[match(
    exposomeBiological_variable_info$variable_id,
    rownames(exposomeBiological_expression_data)
  ),]

###environment
load("data_20200511/environment/expression_data")
environment_expression_data <- expression_data
load("data_20200511/environment/sample_info")
environment_sample_info <- sample_info
load("data_20200511/environment/variable_info")
environment_variable_info <- variable_info

environment_expression_data <- 
  apply(environment_expression_data, 1, function(x){
    x[is.na(x)] <- min(x, na.rm = TRUE)
    x
  }) %>% 
  t() %>% 
  as.data.frame()

environment_variable_info <-
  environment_variable_info %>%
  dplyr::filter(variable_id %in% stringr::str_replace(unique(c(cor_value$from, cor_value$to)),"exposome_", ""))

environment_expression_data <-
  environment_expression_data[match(
    environment_variable_info$variable_id,
    rownames(environment_expression_data)
  ),]

####load cytokine data
load("data_20200511/cytokine/expression_data")
cytokine_expression_data <- expression_data
load("data_20200511/cytokine/sample_info")
cytokine_sample_info <- sample_info
load("data_20200511/cytokine/variable_info")
cytokine_variable_info <- variable_info

cytokine_variable_info <- 
cytokine_variable_info %>% 
  dplyr::left_join(cytokine_class, by = "variable_id") %>% 
  dplyr::filter(classification != "Control")

cytokine_expression_data <- 
  cytokine_expression_data[match(cytokine_variable_info$variable_id, rownames(cytokine_expression_data)),] 

cytokine_variable_info <-
  cytokine_variable_info %>%
  dplyr::filter(variable_id %in% stringr::str_replace(unique(c(cor_value$from, cor_value$to)),"exposome_", ""))

cytokine_expression_data <-
  cytokine_expression_data[match(
    cytokine_variable_info$variable_id,
    rownames(cytokine_expression_data)
  ),]

###
setwd("data_analysis/exposome_cytokine")

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

cytokine_sample_info <- 
  cytokine_sample_info %>% 
  dplyr::rename(start_date = CollectionDate)

exposomeChemical_sample_info$start_date
exposomeBiological_sample_info$start_date
environment_sample_info$start_date

match_date <- 
  data.frame(exposome = c("2016-01-15",
                          "2016-01-19",
                          "2016-01-25",
                          "2016-02-25",
                          "2016-03-03"), 
             cytokine = c("2016-01-15",
                          "2016-01-19",
                          "2016-01-26",
                          "2016-02-24",
                          "2016-03-03")
             )
####exposome vs cytokine
exposomeChemical_sample_info <- 
  exposomeChemical_sample_info %>% 
  dplyr::filter(as.character(start_date) %in% match_date$exposome)

exposomeChemical_expression_data <- 
  exposomeChemical_expression_data %>% 
  dplyr::select(one_of(exposomeChemical_sample_info$sample_id))

exposomeBiological_sample_info <- 
  exposomeBiological_sample_info %>% 
  dplyr::filter(as.character(start_date) %in% match_date$exposome)

exposomeBiological_expression_data <- 
  exposomeBiological_expression_data %>% 
  dplyr::select(one_of(exposomeBiological_sample_info$sample_id))

environment_sample_info <- 
  environment_sample_info %>% 
  dplyr::filter(as.character(start_date) %in% match_date$exposome)

environment_expression_data <- 
  environment_expression_data %>% 
  dplyr::select(one_of(environment_sample_info$sample_id))

cytokine_sample_info <- 
  cytokine_sample_info %>% 
  dplyr::filter(as.character(start_date) %in% match_date$cytokine)

cytokine_expression_data <- 
  cytokine_expression_data %>% 
  dplyr::select(one_of(cytokine_sample_info$sample_id))

colnames(exposomeChemical_expression_data) =
  colnames(exposomeBiological_expression_data) =
  colnames(environment_expression_data) =
  colnames(cytokine_expression_data) =
  as.character(as.Date(exposomeChemical_sample_info$start_date))

####PCA analysis for each exposome
##exposomeChemical
temp_exposomeChemical_expression_data <-
  log(exposomeChemical_expression_data + 1, 2) %>% 
  apply(1, function(x){
    (x - mean(x))/sd(x)
  }) %>% 
  t() %>% 
  as.data.frame()

# exposomeChemical_pca <-
#   prcomp(
#     x = t(temp_exposomeChemical_expression_data),
#     center = FALSE,
#     scale. = FALSE
#   )
# 
# save(exposomeChemical_pca, file = "exposomeChemical_pca")
load("exposomeChemical_pca")

idx <-
  which(summary(exposomeChemical_pca)$importance[3, ] > 0.8)[1]

exposomeChemical_pc <- 
exposomeChemical_pca$x[,1:1, drop = FALSE] %>% 
  t() %>% 
  as.data.frame()

####exposomeBiological
temp_exposomeBiological_expression_data <-
  log(exposomeBiological_expression_data + 1, 2) %>% 
  apply(1, function(x){
    (x - mean(x))/sd(x)
  }) %>% 
  t() %>% 
  as.data.frame()

sum(is.na(temp_exposomeBiological_expression_data))
which(is.na(temp_exposomeBiological_expression_data), arr.ind = TRUE)
temp_exposomeBiological_expression_data[36,]

##remove the feature which are all is NA
remove_idx <-
  apply(temp_exposomeBiological_expression_data, 1, function(x) {
    all(is.na(x))
  }) %>% which()

temp_exposomeBiological_expression_data[remove_idx,]
if(length(remove_idx) > 0){
  temp_exposomeBiological_expression_data <- 
  temp_exposomeBiological_expression_data[-remove_idx,]
}

# exposomeBiological_pca <-
#   prcomp(
#     x = t(temp_exposomeBiological_expression_data),
#     center = FALSE,
#     scale. = FALSE
#   )
# 
# save(exposomeBiological_pca, file = "exposomeBiological_pca")
load("exposomeBiological_pca")

idx <-
  which(summary(exposomeBiological_pca)$importance[3, ] > 0.8)[1]

exposomeBiological_pc <- 
  exposomeBiological_pca$x[,1, drop = FALSE] %>% 
  t() %>% 
  as.data.frame()

####environment
temp_environment_expression_data <-
  environment_expression_data %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

sum(is.na(temp_environment_expression_data))
which(is.na(temp_environment_expression_data), arr.ind = TRUE)
temp_environment_expression_data[6,]
temp_environment_expression_data[8,]

# environment_pca <-
#   prcomp(
#     x = t(temp_environment_expression_data),
#     center = FALSE,
#     scale. = FALSE
#   )
# 
# save(environment_pca, file = "environment_pca")
load("environment_pca")

idx <-
  which(summary(environment_pca)$importance[3, ] > 0.8)[1]

environment_pc <- 
  environment_pca$x[,1,drop = FALSE] %>% 
  t() %>% 
  as.data.frame()

###cytokine
##adjust fiber
temp_cytokine_expression_data <- 
  purrr::map(as.data.frame(t(cytokine_expression_data)), .f = function(x){
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

temp_cytokine_expression_data <-
  temp_cytokine_expression_data %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

sum(is.na(temp_cytokine_expression_data))

# cytokine_pca <-
#   prcomp(
#     x = t(temp_cytokine_expression_data),
#     center = FALSE,
#     scale. = FALSE
#   )
# 
# save(cytokine_pca, file = "cytokine_pca")
load("cytokine_pca")

idx <-
  which(summary(cytokine_pca)$importance[3, ] > 0.8)[1]

cytokine_pc <- 
  cytokine_pca$x[,1,drop = FALSE] %>% 
  t() %>% 
  as.data.frame()

#####multiple linear regression
rownames(exposomeChemical_pc) <- 
  paste("exposomeChemical", rownames(exposomeChemical_pc), sep = "_")

rownames(exposomeBiological_pc) <- 
  paste("exposomeBiological", rownames(exposomeBiological_pc), sep = "_")

rownames(environment_pc) <- 
  paste("environment", rownames(environment_pc), sep = "_")

exposome_pc <- rbind(exposomeChemical_pc,
                     exposomeBiological_pc,
                     environment_pc)

# total_r2 <- 
# temp_cytokine_expression_data %>% 
#   t() %>% 
#   as.data.frame() %>% 
#   purrr::map(.f = function(x){
#     temp_data <-
#       rbind(y = x,
#             exposome_pc) %>%
#       t() %>%
#       as.data.frame()
#     
#     lm_object <-
#       lm(
#         formula = y ~ .,
#         data = temp_data
#       ) 
#     summary(lm_object)$r.squared
#   }) %>% 
#   unlist()
# save(total_r2, file = "total_r2")
load("total_r2")

###PLSR
library(pls)
library(plsVarSel)

# total_vip <- 
#   temp_cytokine_expression_data %>% 
#   t() %>% 
#   as.data.frame() %>% 
#   purrr::map(.f = function(x){
#     temp_data <-
#       rbind(y = x,
#             exposome_pc) %>%
#       t() %>%
#       as.data.frame()
#     
#     plsr_object <-
#       plsr(y ~ ., ncomp = 3, validation = "LOO", data = temp_data)
#     vip <- VIP(plsr_object, which.min(plsr_object$validation$PRESS))
#     vip/sum(vip)
#   })
# 
# save(total_vip, file = "total_vip")
load("total_vip")

# total_explain <- 
# purrr::map2(
#   .x = total_r2,
#   .y = total_vip,
#   .f = function(x, y) {
#     z <- x * y
#     z <- c(z, 1 - sum(z))
#     names(z)[4] <- "Other"
#     names(z) <- stringr::str_replace_all(names(z), "_PC1", "")
#     z
#   }
# )
# 
# total_explain <- 
# purrr::map2(
#   .x = names(total_explain),
#   .y = total_explain,
#   .f = function(x, y) {
#     cbind(variable_id = x, data.frame(percent = y)) %>%
#       tibble::rownames_to_column(var = "class") %>%
#       dplyr::select(variable_id, everything())
#   }
# ) %>% 
#   do.call(rbind, .) %>% 
#   as.data.frame()
# 
# save(total_explain, file = "total_explain")
load("total_explain")

value1 <- 
  c("environment" = ggsci::pal_d3()(10)[1],
    "exposomeChemical" = ggsci::pal_d3()(10)[2],
    "Metabolome" = ggsci::pal_d3()(10)[3],
    "Proteome" = ggsci::pal_d3()(10)[4],
    "exposomeBiological" = ggsci::pal_d3()(10)[5],
    "Gut microbiome" = ggsci::pal_d3()(10)[6],
    "Blood test" = ggsci::pal_d3()(10)[7],
    "Cytokine" = ggsci::pal_d3()(10)[8],
    "Toxins and carcinogens" = ggsci::pal_d3()(10)[9],
    "Other" = "black"
  )

other_percent <- 
  total_explain %>% 
  dplyr::filter(class == "Other") %>% 
  dplyr::arrange(percent)

# write.csv(other_percent, "other_percent_new.csv", row.names = FALSE)

total_explain_output = 
  total_explain %>%
  dplyr::mutate(percent = percent * 100) %>% 
  tidyr::pivot_wider(names_from = class,values_from = percent)

total_explain_output <-
  total_explain_output %>%
  dplyr::left_join(cytokine_variable_info[, 1:2],
                   by = "variable_id") 

# library(openxlsx)
# wb = createWorkbook()
# modifyBaseFont(wb, fontSize = 12, fontName = "Times New Roma")
# addWorksheet(wb, sheetName = "Cytokine", gridLines = TRUE)
# freezePane(wb,
#            sheet = 1,
#            firstRow = TRUE,
#            firstCol = TRUE)
# writeDataTable(
#   wb,
#   sheet = 1,
#   x = total_explain_output,
#   colNames = TRUE,
#   rowNames = FALSE
# )
# saveWorkbook(wb, "total_explain.xlsx", overwrite = TRUE)

total_explain <- 
total_explain %>%
  dplyr::mutate(percent = percent * 100) %>%
  dplyr::mutate(class = factor(class,
                               levels = rev(
                                 c(
                                   "exposomeChemical",
                                   "exposomeBiological",
                                   "environment",
                                   "Other"
                                 )
                               ))) %>%
  dplyr::mutate(variable_id = factor(variable_id,
                                     levels = other_percent$variable_id))
text_color <- 
total_explain$variable_id %>% levels() %>% 
  data.frame(variable_id = .) %>% 
  dplyr::left_join(cytokine_variable_info, by = "variable_id") %>% 
  dplyr::pull(classification)

text_color[text_color == "Anti-inflammatory"] <- "red"
text_color[text_color == "Proinflammatory"] <- "blue"
text_color[text_color == "Proinflammatory/Anti-inflammatory"] <- "purple"

plot1 <-
  total_explain %>%
  ggplot(aes(variable_id, y = percent)) +
  geom_bar(aes(fill = class), stat = "identity") +
  scale_fill_manual(values = value1) +
  scale_x_discrete(expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  ggplot2::theme_bw() +
  labs(x = "", y = "Percent (%)") +
  theme(axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    vjust = 1,
    color = text_color
  ),
  axis.ticks.x = element_blank(),
  legend.position = "bottom")

plot1

###get the network of all the exposome with cytokine
library(igraph)
library(ggraph)
library(tidygraph)

edge_data <-
  cor_value %>%
  dplyr::rename(from = from,
                to = to,
                Correlation = cor) %>%
  dplyr::mutate(fdr = -log(fdr, 10)) %>% 
  dplyr::mutate(from = stringr::str_replace(from, "exposome_", "")) %>% 
  dplyr::mutate(to = stringr::str_replace(to, "exposome_", "")) 

node_data <-
  cor_value %>%
  dplyr::select(from, to) %>%
  tidyr::pivot_longer(cols = c(from, to),
                      names_to = "class",
                      values_to = "node") %>%
  dplyr::mutate(node = stringr::str_replace(node, "exposome_", "")) %>% 
  dplyr::mutate(
    class1 = case_when(
      node %in% exposomeChemical_variable_info$peak_ID ~ "Exposome (chemical)",
      node %in% exposomeBiological_variable_info$variable_id ~ "Exposome (biological)",
      node %in% environment_variable_info$variable_id ~ "Environment",
      node %in% cytokine_variable_info$variable_id ~ "Cytokine"
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
      "Environment",
      "Cytokine"
    )
  ))

node_data$true_name <-
  node_data$true_name %>% 
  stringr::str_replace("exposome_", "")

##remove control cytokine
node_data <-
  node_data %>%
  dplyr::filter(!node %in% cytokine_class$variable_id[cytokine_class$classification == "Control"])

edge_data <- 
  edge_data %>% 
  dplyr::filter(!from %in% cytokine_class$variable_id[cytokine_class$classification == "Control"],
                !to %in% cytokine_class$variable_id[cytokine_class$classification == "Control"])

library(plyr)
node_data <- 
node_data %>%
  plyr::dlply(.variables = .(Class)) %>% 
  purrr::map(.f = function(x){
    if(unique(x$Class) == 'Cytokine'){
      x <- 
        x %>% 
        dplyr::left_join(data.frame(node = levels(total_explain$variable_id)),
                         .,
                         by = "node")
    }else{
      x <- x
    }
    x
  }) %>% 
  do.call(rbind, .) %>% 
  as.data.frame()

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

library(igraph)

g <- total_graph

V(g)$type <- bipartite_mapping(g)$type

coords <-
  create_layout(g, layout = "bipartite") %>%
  dplyr::select(x, y, node, Class) %>%
  tibble::rowid_to_column("index") %>% 
  dplyr::arrange(x) %>%
  dplyr::mutate(
    y = case_when(
      Class == "Exposome (chemical)" ~ 1,
      Class == "Exposome (biological)" ~ 1,
      Class == "Environment" ~ 1,
      Class == "Cytokine" ~ -8
    )
  ) %>% 
  plyr::dlply(.variables = .(Class)) %>% 
  purrr::map(.f = function(z){
      z$x <- 
        seq(0, 130, length.out = sum(nrow(z)))
      z
  }) %>% 
  do.call(rbind, .) %>% 
  as.data.frame() %>% 
  dplyr::arrange(index) 

library(plyr)
coords <-
  coords %>%
  plyr::dlply(.variables = .(Class)) %>%
  purrr::map(
    .f = function(z) {
      if (unique(z$Class) == "Exposome (chemical)") {
        z <-
          z %>%
          dplyr::mutate(
            theta = x / (max(x) + 1) * 2 * pi,
            r = y + 2,
            x = r * cos(theta),
            y = r * sin(theta)
          ) %>%
          dplyr::mutate(x = x/5) 
      }
      
      if (unique(z$Class) == "Exposome (biological)") {
        z <-
          z %>%
          dplyr::mutate(
            theta = x / (max(x) + 1) * 2 * pi,
            r = y + 4,
            x = r * cos(theta),
            y = r * sin(theta)
          ) %>%
          dplyr::mutate(x = x/5)
      }
      
      if (unique(z$Class) == "Environment") {
        z <-
          z %>%
          dplyr::mutate(
            theta = x / (max(x) / 2 + 1) * 2 * pi,
            r = y,
            x = r * cos(theta),
            y = r * sin(theta)
          ) %>%
          dplyr::mutate(x = x) 
        
        if (nrow(z) <= 3) {
          z$x <- mean(z$x)
          z$y <- seq(min(coords$y[coords$Class != "Cytokine"]), 
                     max(coords$y[coords$Class != "Cytokine"]), 
                     length.out = nrow(z))
        }
        
      }
      
      if (unique(z$Class) == "Cytokine") {
        z <-
          z %>%
          dplyr::mutate(theta = 0, r = 0)
      }
      return(z)
    }
  ) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  dplyr::arrange(index)

coords$y[coords$Class == "Environment"] <- 
seq(from = min(coords$y[coords$Class != "Cytokine"]),
    to = max(coords$y[coords$Class != "Cytokine"]),
    length.out = sum(coords$Class == "Environment")
)

angle <- 
  coords %>% 
  dplyr::mutate(angle = -((-node_angle(x, y) + 90) %% 180) + 90) %>% 
  dplyr::pull(angle)

coords2 <-
  coords %>%
  plyr::dlply(.variables = .(Class)) %>%
  purrr::map(
    .f = function(z) {
      if (unique(z$Class) == "Exposome (chemical)") {
        z <-
          z %>%
          dplyr::mutate(x = x - 4)
      }
      
      if (unique(z$Class) == "Exposome (biological)") {
        z <-
          z %>%
          dplyr::mutate(x = x)
      }
      
      if (unique(z$Class) == "Environment") {
        z <-
          z %>%
          dplyr::mutate(x = x  + 3)
      }
      
      if (unique(z$Class) == "Cytokine") {
        z <- z
      }
      return(z)
    }
  ) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  dplyr::arrange(index)


range(coords2$x[coords2$Class != "Cytokine"])
coords2$x[coords2$Class == "Cytokine"] <-
  seq(min(coords2$x[coords2$Class != "Cytokine"]),
      max(coords2$x[coords2$Class != "Cytokine"]), 
      length.out = sum(coords2$Class == "Cytokine"))

my_graph <-
  create_layout(
    graph = g,
    layout = "manual",
    x = coords2$x,
    y = coords2$y
    # node.position = coords
  )

plot2 <-
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
      label = ifelse(Class == "Cytokine", NA, true_name),
      hjust = "outward",
      angle = angle,
      size = 3,
      colour = Class
    ),
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
  labs(x = "", y = "") +
  theme(
    legend.position = "top",
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) 
  # coord_fixed(ratio = 0.5)
# plot(coords$x, coords$y)

plot2

library(patchwork)

plot <- 
plot2 + 
  plot1 +
  patchwork::plot_layout(ncol = 1, heights = c(1,1))

plot

# ggsave(plot,
#        filename = "exposome_on_cytokine.pdf",
#        width = 10,
#        height = 7.5)

edge_info <- edge_data
edge_info <- 
  edge_info %>% 
  dplyr::left_join(node_data, by = c("from" = "node")) %>% 
  dplyr::rename(from_class = Class, from_name = true_name) %>% 
  dplyr::left_join(node_data, by = c("to" = "node")) %>% 
  dplyr::rename(to_class = Class, to_name = true_name) 

from_num <-
  edge_info %>% 
  dplyr::group_by(from) %>% 
  dplyr::summarise(from_num = n())

to_num <-
  edge_info %>% 
  dplyr::group_by(to) %>% 
  dplyr::summarise(to_num = n())

edge_info <- 
  edge_info %>% 
  dplyr::left_join(from_num, by = "from") %>% 
  dplyr::left_join(to_num, by = "to") 

write.csv(edge_info, file = "exposome vs cytokine_edge_info.csv", row.names = FALSE) 

