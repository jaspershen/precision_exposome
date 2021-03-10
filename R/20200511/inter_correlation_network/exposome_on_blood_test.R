##avoid source
no_function()
rm(list = ls())
sxtTools::setwd_project()
##https://blog.csdn.net/fjsd155/article/details/84726785

library(pls)
library(tidyverse)
data(yarn)
data(mtcars)

#####exposomeChemical
load("data_analysis/exposomeChemical_blood_test/cor_value")
exposomeChemical_blood_test_cor <-
  cor_value %>%
  dplyr::filter(abs(cor) > 0.9 & p_value < 0.05)

dim(exposomeChemical_blood_test_cor)

#####exposomeBiological
load("data_analysis/exposomeBiological_blood_test/cor_value")
exposomeBiological_blood_test_cor <-
  cor_value %>%
  dplyr::filter(abs(cor) > 0.9 & p_value < 0.05) %>%
  dplyr::mutate(from = paste("exposome", from, sep = "_"))

#####environment
load("data_analysis/environment_blood_test/cor_value")
environment_blood_test_cor <-
  cor_value %>%
  dplyr::filter(abs(cor) > 0.9 & p_value < 0.05) %>%
  dplyr::mutate(from = paste("exposome", from, sep = "_"))

cor_value <-
  rbind(
    exposomeChemical_blood_test_cor,
    exposomeBiological_blood_test_cor,
    environment_blood_test_cor
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

environment_variable_info <-
  environment_variable_info %>%
  dplyr::filter(variable_id %in% stringr::str_replace(unique(c(cor_value$from, cor_value$to)),"exposome_", ""))

environment_expression_data <-
  environment_expression_data[match(
    environment_variable_info$variable_id,
    rownames(environment_expression_data)
  ),]

####load blood_test data
load("data_20200511/lab_test/expression_data")
blood_test_expression_data <- expression_data
load("data_20200511/lab_test/sample_info")
blood_test_sample_info <- sample_info
load("data_20200511/lab_test/variable_info")
blood_test_variable_info <- variable_info

blood_test_expression_data <- 
  blood_test_expression_data[match(blood_test_variable_info$variable_id, rownames(blood_test_expression_data)),] 

blood_test_variable_info <-
  blood_test_variable_info %>%
  dplyr::filter(variable_id %in% stringr::str_replace(unique(c(cor_value$from, cor_value$to)),"exposome_", ""))

blood_test_expression_data <-
  blood_test_expression_data[match(
    blood_test_variable_info$variable_id,
    rownames(blood_test_expression_data)
  ),]

###
setwd("data_analysis/exposome_blood_test")

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

exposomeChemical_sample_info$start_date
exposomeBiological_sample_info$start_date
environment_sample_info$start_date

blood_test_sample_info$start_date

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

colnames(exposomeChemical_expression_data) =
  colnames(exposomeBiological_expression_data) =
  colnames(environment_expression_data) =
  colnames(blood_test_expression_data) =
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

exposomeChemical_pca <-
  prcomp(
    x = t(temp_exposomeChemical_expression_data),
    center = FALSE,
    scale. = FALSE
  )

idx <-
  which(summary(exposomeChemical_pca)$importance[3, ] > 0.8)[1]

exposomeChemical_pc <-
  exposomeChemical_pca$x[, 1:idx] %>%
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
  exposomeBiological_pca$x[,1:idx] %>% 
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
  environment_pca$x[,1:idx] %>% 
  t() %>% 
  as.data.frame()

rownames(environment_pc) <- paste("PC", 1:nrow(environment_pc), sep = "")

###blood_test
##adjust fiber
temp_blood_test_expression_data <- 
  purrr::map(as.data.frame(t(blood_test_expression_data)), .f = function(x){
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

temp_blood_test_expression_data <-
  temp_blood_test_expression_data %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

sum(is.na(temp_blood_test_expression_data))

# blood_test_pca <-
#   prcomp(
#     x = t(temp_blood_test_expression_data),
#     center = FALSE,
#     scale. = FALSE
#   )
# 
# save(blood_test_pca, file = "blood_test_pca")
load("blood_test_pca")
idx <-
  which(summary(blood_test_pca)$importance[3, ] > 0.8)[1]

blood_test_pc <- 
  blood_test_pca$x[,1:idx] %>% 
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
# temp_blood_test_expression_data %>%
#   t() %>%
#   as.data.frame() %>%
#   purrr::map(.f = function(x){
#     temp_data <-
#       rbind(y = x,
#             exposome_pc[c(1,3,6),]) %>%
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
#
# save(total_r2, file = "total_r2")

load("total_r2")

###PLSR
library(pls)
library(plsVarSel)

# total_vip <- 
#   temp_blood_test_expression_data %>% 
#   t() %>% 
#   as.data.frame() %>% 
#   purrr::map(.f = function(x){
#     temp_data <-
#       rbind(y = x,
#             exposome_pc[c(1,3,6),]) %>%
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
  dplyr::left_join(blood_test_variable_info[, 1:2],
                   by = "variable_id") %>%
  dplyr::arrange(percent)

# write.csv(other_percent, "other_percent.csv", row.names = FALSE)
library(tidyr)

# total_explain_output = 
#   total_explain %>%
#   dplyr::mutate(percent = percent * 100) %>% 
#   tidyr::pivot_wider(names_from = class,values_from = percent)
# 
# total_explain_output <-
#   total_explain_output %>%
#   dplyr::left_join(blood_test_variable_info[, 1:2],
#                    by = "variable_id") 
# 
# library(openxlsx)
# wb = createWorkbook()
# modifyBaseFont(wb, fontSize = 12, fontName = "Times New Roma")
# addWorksheet(wb, sheetName = "Blood test", gridLines = TRUE)
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
  dplyr::left_join(blood_test_variable_info[, 1:2],
                   by = "variable_id") %>% 
  dplyr::mutate(true_name = factor(true_name,
                                     levels = other_percent$true_name))

plot1 <-
  total_explain %>%
  ggplot(aes(true_name, y = percent)) +
  geom_bar(aes(fill = class), stat = "identity") +
  scale_fill_manual(values = value1) +
  scale_x_discrete(expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  ggplot2::theme_bw() +
  labs(x = "", y = "Percent (%)") +
  theme(axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    vjust = 1
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
      node %in% blood_test_variable_info$variable_id ~ "Blood test"
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
  dplyr::left_join(blood_test_variable_info[, c("variable_id", "true_name")] %>% 
                     dplyr::rename(blood_name = true_name,
                                   node = variable_id),
                   by = c("node")) %>%
  dplyr::mutate(true_name = case_when(!is.na(blood_name) ~ blood_name,
                                      TRUE ~ true_name)) %>%
  dplyr::select(node, Class, true_name) %>%
  dplyr::mutate(Class = factor(
    Class,
    levels = c(
      "Exposome (chemical)",
      "Exposome (biological)",
      "Environment",
      "Blood test"
    )
  ))

node_data$true_name <-
  node_data$true_name %>% 
  stringr::str_replace("exposome_", "")

library(plyr)
node_data <- 
  node_data %>%
  plyr::dlply(.variables = .(Class)) %>% 
  purrr::map(.f = function(x){
    if(unique(x$Class) == 'Blood test'){
      x <- 
        x %>%
        dplyr::left_join(data.frame(true_name = levels(total_explain$true_name)),
                         .,
                         by = "true_name") %>% 
        dplyr::select(node, Class, true_name)
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
      Class == "Blood test" ~ -10
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
            r = y + 4,
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
            r = y + 6,
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
          z$y <- seq(min(coords$y[coords$Class != "Blood test"]), 
                     max(coords$y[coords$Class != "Blood test"]), 
                     length.out = nrow(z))
        }
        
      }
      
      if (unique(z$Class) == "Blood test") {
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
  seq(from = min(coords$y[coords$Class != "Blood test"]),
      to = max(coords$y[coords$Class != "Blood test"]),
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
      
      if (unique(z$Class) == "Blood test") {
        z <- z
      }
      return(z)
    }
  ) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  dplyr::arrange(index)


range(coords2$x[coords2$Class != "Blood test"])
coords2$x[coords2$Class == "Blood test"] <-
  seq(min(coords2$x[coords2$Class != "Blood test"]),
      max(coords2$x[coords2$Class != "Blood test"]), 
      length.out = sum(coords2$Class == "Blood test"))

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
      label = ifelse(Class == "Blood test", NA, true_name),
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

# write.csv(edge_info, file = "exposome vs blood test_edge_info.csv", row.names = FALSE) 

# ggsave(plot2, filename = "exposome_on_blood_test2.pdf", width = 6.67, height = 4)

plot <- 
  plot2 + 
  plot1 +
  patchwork::plot_layout(ncol = 1, heights = c(1,1))

plot

# ggsave(plot, filename = "exposome_on_blood_test.pdf", width = 6.67, height = 8)

####read the range of blood test
sxtTools::setwd_project()

blood_test_range <-
  readxl::read_xlsx("data_analysis/lab_test_analysis/data_overview/Blood range.xlsx")

blood_test_range <- 
  blood_test_range %>% 
  dplyr::rename(true_name = `Full name`)

setwd("data_analysis/exposome_blood_test")

blood_test_variable_info <-
  blood_test_variable_info %>%
  dplyr::left_join(blood_test_range[, -2], by = c("variable_id")) 

temp_data <-
  blood_test_expression_data %>%
  tibble::rownames_to_column(var = "variable_id") %>%
  tidyr::pivot_longer(cols = -variable_id,
                      names_to = "date",
                      values_to = "value") %>%
  dplyr::left_join(blood_test_variable_info, by = "variable_id")


##remove the blood test which are all normal
library(plyr)
temp_data <-
  temp_data %>%
  plyr::dlply(.variables = .(variable_id)) %>%
  purrr::map(
    .f = function(x) {
      if (any(x$value < x$Minimum) | any(x$value > x$Maximum)) {
        return(x)
      }
    }
  ) %>%
  do.call(rbind, .) %>%
  as.data.frame()

##remove the blood test which changed after adjusted
temp_data <- 
  temp_data %>% 
  dplyr::filter(!variable_id %in% c("Cholesterol._Total", 
                                    "eGFR", "LDL_.Calculated.",
                                    "LDL.HDL_Ratio",
                                    "MCHC", "MCV"))

unique(temp_data$variable_id)

temp_data <-
  temp_data %>%
  dplyr::mutate(variable_id =
                  factor(variable_id,
                         levels = unique(total_explain$variable_id)[unique(total_explain$variable_id) %in% temp_data$variable_id]))

plot1 <- 
temp_data %>%
  ggplot(aes(date, value, group = true_name)) +
  geom_point() +
  geom_line(aes(),
            show.legend = FALSE) +
  labs(x = "", y = "Value") +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      size = 10,
      hjust = 1,
      vjust = 1
    ),
    axis.title = element_text(size = 10),
    panel.grid = element_blank()
  ) +
  geom_rect(
    aes(
      xmin = -Inf,
      xmax = Inf,
      ymin = Minimum,
      ymax = Maximum,
      group = true_name
    ),
    fill = ggsci::pal_aaas()(n = 10)[5],
    alpha = 0.5,
    data = blood_test_range %>%
      dplyr::filter(variable_id %in% temp_data$variable_id),
    inherit.aes = FALSE
  ) +
  facet_wrap(facets = vars(true_name),
             scales = "free_y")

plot1

# ggsave(plot1, filename = "blood_test_range_plot.pdf", width = 10, height = 7)


###change the range of blood test
temp_data$Maximum[temp_data$true_name == "Hemoglobin"] = 13.5
temp_data$Minimum[temp_data$true_name == "Urea nitrogen"] = 15

##get the relationship with each blood test
dir.create("final_blood_test_plot")

# for (idx in 1:length(as.character(unique(temp_data$variable_id)))) {
#   cat(idx, " ")
#   i <- as.character(unique(temp_data$variable_id))[idx]
#   temp1 <-
#     edge_data %>%
#     dplyr::filter(from %in% i |
#                     to %in% i)
#   
#   pal <-
#     wesanderson::wes_palette(name = "Zissou1", n = 100, type = "continuous")
#   
#   plot3 <- 
#   temp1 %>% 
#     dplyr::mutate(from = factor(from, levels = rev(from))) %>%
#     ggplot(aes(to, from)) +
#     geom_tile(aes(fill = Correlation), color = "white") +
#     labs(x = "", y = "") +
#     scale_fill_gradientn(colours = pal) +
#     theme_bw() +
#     scale_x_discrete(expand = expansion(mult = c(0,0))) +
#     scale_y_discrete(expand = expansion(mult = c(0,0))) +
#     theme(
#       panel.grid = element_blank(),
#       axis.text = element_blank(),
#       axis.ticks = element_blank(),
#       legend.position = "bottom"
#     )
#   
#   
#   # plot3 <-
#   #   temp1 %>%
#   #   dplyr::mutate(from = factor(from, levels = rev(from))) %>%
#   #   ggplot(aes(to, from)) +
#   #   geom_point(aes(fill = Correlation), shape = 21, size = 10) +
#   #   labs(x = "", y = "") +
#   #   scale_fill_gradientn(colours = pal) +
#   #   theme_bw() +
#   #   theme(
#   #     panel.grid = element_blank(),
#   #     axis.text = element_blank(),
#   #     axis.ticks = element_blank(),
#   #     legend.position = "bottom"
#   #   )
#   
#   temp1_1 <- 
#   exposomeChemical_expression_data[match(temp1$from, rownames(exposomeChemical_expression_data)),] %>%
#     tibble::rownames_to_column(var = "variable_id") %>%
#     dplyr::filter(!stringr::str_detect(variable_id, "NA"))
#   
#   temp1_1$variable_id <- 
#     exposomeChemical_variable_info$MetabID[match(temp1_1$variable_id, exposomeChemical_variable_info$peak_ID)]
#   
#   temp1_2 <- 
#   exposomeBiological_expression_data[match(temp1$from, rownames(exposomeBiological_expression_data)),] %>%
#     tibble::rownames_to_column(var = "variable_id") %>%
#     dplyr::filter(!stringr::str_detect(variable_id, "NA"))
#   
#   temp1_3 <- 
#   environment_expression_data[match(temp1$from, rownames(environment_expression_data)),] %>%
#     tibble::rownames_to_column(var = "variable_id") %>%
#     dplyr::filter(!stringr::str_detect(variable_id, "NA"))
#   
#   temp1_1 <-
#     temp1_1 %>% 
#     tibble::column_to_rownames(var = "variable_id") %>% 
#     `+`(1) %>% 
#     log(2) %>% 
#     apply(1, function(x){
#       (x - mean(x)) / sd(x)
#     }) %>% 
#     t() %>% 
#     as.data.frame()
#   
#   temp1_2 <-
#     temp1_2 %>% 
#     tibble::column_to_rownames(var = "variable_id") %>% 
#     `+`(1) %>% 
#     log(2) %>% 
#     apply(1, function(x){
#       (x - mean(x)) / sd(x)
#     }) %>% 
#     t() %>% 
#     as.data.frame()
#   
#   temp1_3 <-
#     temp1_3 %>% 
#     tibble::column_to_rownames(var = "variable_id") %>% 
#     `+`(1) %>% 
#     log(2) %>% 
#     apply(1, function(x){
#       (x - mean(x)) / sd(x)
#     }) %>% 
#     t() %>% 
#     as.data.frame()
#   
#   temp1 <- 
#     rbind(temp1_1,
#           temp1_2,
#           temp1_3)
#   
#   col_fun = colorRamp2(
#     breaks = seq(min(temp1), max(temp1), length.out = 90),
#     colors =
#       viridis::magma(n = 100)[-c(1:10)],
#     transparency = 0
#   )
#   
#   temp_idx <- which.max(nchar(rownames(temp1)))
#   if(nchar(rownames(temp1)[temp_idx]) < 35){
#     rownames(temp1)[temp_idx] <- 
#       paste(c(rownames(temp1)[temp_idx],
#             rep("_", 35 - nchar(rownames(temp1))[temp_idx])), collapse = "")
#   }
#   
#   library(ComplexHeatmap)
#   plot1 <-
#     temp1 %>%
#     ComplexHeatmap::Heatmap(
#       cluster_columns = FALSE,
#       cluster_rows = FALSE,
#       show_column_names = FALSE,
#       show_row_names = TRUE,
#       clustering_method_rows = "ward.D",
#       clustering_method_columns = "ward.D",
#       clustering_distance_columns = "euclidean",
#       clustering_distance_rows = "euclidean",
#       col = col_fun,
#       rect_gp = gpar(col = "white"),
#       border = TRUE,
#       row_dend_reorder = TRUE,
#       column_dend_reorder = TRUE,
#       column_names_rot = 45,
#       name = "Z-score",
#       cell_fun = function(j, i, x, y, width, height, fill) {
#         grid.text(sprintf("%.1f", temp1[i, j]), x, y, 
#                   gp = gpar(fontsize = 15, col = "white"))
#       }
#     )
#   
#   plot1 <- ggplotify::as.ggplot(plot1)
#   
#   temp2 <-
#     temp_data %>%
#     dplyr::filter(variable_id %in% i)
#   
#   plot2 <- 
#   temp2 %>%
#     ggplot(aes(date, value, group = true_name)) +
#     geom_rect(
#       aes(
#         xmin = -Inf,
#         xmax = Inf,
#         ymin = Minimum,
#         ymax = Maximum,
#       ),
#       fill = ggsci::pal_aaas()(n = 10)[5],
#       alpha = 0.5,
#       data = temp2 %>% dplyr::distinct(variable_id, .keep_all = TRUE)
#     ) +
#     geom_point(size = 4, 
#                shpae = 21,
#                fill = "red") +
#     geom_line(aes(),
#               show.legend = FALSE) +
#     labs(
#       x = "",
#       y = paste("Value (", temp2$unit[1], ")", sep = ""),
#       title = temp2$true_name[1]
#     ) +
#     theme_bw() +
#     theme(
#       axis.text.x = element_blank(),
#       axis.ticks.x = element_blank(),
#       axis.text.y = element_text(size = 10),
#       axis.title = element_text(size = 10),
#       panel.grid = element_blank()
#     )
#   
#   library(patchwork)
# 
# #   plot <-
# #     plot2 + plot1 + patchwork::plot_layout(ncol = 1)
#  
#   plot <-
#   (
#     patchwork::plot_spacer() + plot2 + patchwork::plot_spacer() + patchwork::plot_layout(widths = c(1.5, 14, 7.3))
#   ) / (plot3 + plot1 + patchwork::plot_layout(widths = c(1, 8)))
# 
#   plot  
#   ggsave(plot, filename = file.path("final_blood_test_plot", paste(temp2$true_name[1], ".pdf", sep = "")),
#          width = 10.58, height = 8)
# }


total_explain
