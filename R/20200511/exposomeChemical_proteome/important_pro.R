##load data
sxtTools::setwd_project()
library(tidyverse)
rm(list = ls())
####load proteome data
sxtTools::setwd_project()
load("data_20200511/proteome/expression_data")
pro_expression_data <- expression_data
load("data_20200511/proteome/sample_info")
pro_sample_info <- sample_info
load("data_20200511/proteome/variable_info")
pro_variable_info <- variable_info

load("data_20200511/exposome/expression_data")
exp_expression_data <- expression_data
load("data_20200511/exposome/sample_info")
exp_sample_info <- sample_info
load("data_20200511/exposome/variable_info")
exp_variable_info <- variable_info

setwd("data_analysis/pro_exp/")

dim(exp_variable_info)

dim(pro_variable_info)

temp1 <- 
  exp_sample_info$start_date %>% 
  as.Date() %>% 
  data.frame(date = ., 
             exp = 1,
             stringsAsFactors = FALSE)

temp2 <- 
  pro_sample_info$CollectionDate %>% 
  as.Date() %>% 
  data.frame(date = ., 
             pro = 1,
             stringsAsFactors = FALSE)
temp <- 
  temp1 %>% 
  dplyr::full_join(temp2, by = "date") %>% 
  arrange(date)

temp <- 
  temp %>% 
  tidyr::pivot_longer(cols = -date, names_to = "class", values_to = "value")

temp <- temp %>% 
  dplyr::filter(!is.na(value))


diff <- 
  lapply(pro_sample_info$CollectionDate %>% as.Date(), function(x){
    x <- x - as.Date(exp_sample_info$start_date) 
    as.Date(exp_sample_info$start_date) [which(x >=0 & x <= 7)]
  }) 

names(diff) <- pro_sample_info$CollectionDate %>% as.character()

diff <- 
  purrr::map2(diff, names(diff), .f = function(x,y){
    if(length(x) == 0){
      return(NULL)
    }
    z <- data.frame(y, x, stringsAsFactors = FALSE)
    colnames(z) <- c("pro", "exp")
    z
  }) 

diff <- diff[lapply(diff, is.null) %>% unlist() %>% `!`]

diff[[2]] <- diff[[2]][-c(1),]

diff[[3]] <- diff[[3]][-c(1),]

diff <- purrr::map2(.x = diff, .y = 1:length(diff), .f = function(x,y){
  data.frame(x, group = y, stringsAsFactors = FALSE)
}) %>% 
  do.call(rbind, .)

rownames(diff) <- NULL

diff <- 
  diff %>% 
  mutate(exp = as.character(exp), 
         pro = as.character(pro)) %>% 
  tidyr::pivot_longer(cols = -group,
                      names_to = "class", 
                      values_to = "date")

temp <-
  temp %>% 
  mutate(date = as.character(date)) %>% 
  left_join(diff, by = c("date","class"))

temp$group[is.na(temp$group)] <- "No"


diff

library(plyr)

diff <-
  diff %>% 
  plyr::dlply(.variables = .(group)) %>%
  lapply(function(x){
    x %>% 
      plyr::dlply(.variables = .(class)) %>% 
      do.call(cbind, .)
  }) %>% 
  do.call(rbind, .)


diff$exp.date <- as.Date(diff$exp.date)

diff$pro.date <- as.Date(diff$pro.date)

diff$exp.group <- as.character(diff$exp.group)


plot <- 
  temp %>% 
  mutate(date = as.Date(date)) %>% 
  ggplot(aes(date, class)) +
  geom_point(aes(color = group),
             alpha = 0.6, size = 4,
             shape = 16, show.legend = FALSE) +
  geom_segment(aes(x = exp.date, 
                   xend = pro.date, 
                   y = exp.class, 
                   yend = pro.class,
                   color = exp.group), data = diff, show.legend = FALSE) +
  scale_color_manual(values = c("No" = "grey", 
                                "1" = ggsci::pal_aaas()(8)[1],
                                "2" = ggsci::pal_aaas()(8)[2],
                                "3" = ggsci::pal_aaas()(8)[3],
                                "4" = ggsci::pal_aaas()(8)[4],
                                "5" = ggsci::pal_aaas()(8)[5])) +
  scale_x_continuous(trans = "date",
                     breaks = c(as.Date(temp$date)),
                     labels = as.character(temp$date)
  ) +
  # scale_y_continuous(breaks = c("pro", "exp"), labels = c("pro", "Exp")) +
  scale_y_discrete(breaks = c("pro", "exp"), labels = c("proteome", "Exposome")) +
  labs(x = "", y = "") +
  # ggrepel::geom_label_repel(aes(label = as.character(date))) +
  theme_bw() +
  theme(axis.title = element_text(size = 13),
        axis.text.x = element_text(size = 10, 
                                   angle = 45,
                                   vjust = 1, hjust = 1),
        axis.text = element_text(size = 10),
        panel.grid.minor = element_blank())


plot

#-------------------------------------------------------------------------------
###prepare data
dim(pro_expression_data)
dim(pro_sample_info)

pro_sample_info <-
  pro_sample_info %>% 
  dplyr::filter(CollectionDate %in% unique(diff$pro.date))

pro_expression_data <- 
  pro_expression_data %>% 
  dplyr::select(one_of(pro_sample_info$sample_id))

pro_sample_info <- 
  pro_sample_info %>% 
  dplyr::left_join(diff[,c("pro.date", "pro.group")] %>% dplyr::distinct(pro.date, pro.group), 
                   by = c("CollectionDate" = "pro.date")) %>% 
  dplyr::rename(group = pro.group)


exp_sample_info <-
  exp_sample_info %>% 
  dplyr::filter(as.Date(start_date) %in% diff$exp.date)


exp_expression_data <-
  exp_expression_data %>%
  dplyr::select(-contains("Blank")) %>% 
  t()

exp_expression_data <-
  exp_expression_data %>%
  as.data.frame() %>%
  rownames_to_column(var = "Sample_ID") %>%
  mutate(Sample_ID = stringr::str_replace(Sample_ID, "_[0-9]{1,2}", "")) %>%
  plyr::dlply(.variables = "Sample_ID") %>%
  lapply(function(x){
    apply(x[,-1], 2, mean)
  }) %>%
  do.call(rbind, .) %>% 
  t() %>% 
  as.data.frame()

exp_expression_data <- 
  exp_expression_data %>% 
  dplyr::select(one_of(exp_sample_info$sample_id))

exp_sample_info <- 
  exp_sample_info %>% 
  mutate(start_date = as.Date(start_date)) %>% 
  dplyr::left_join(diff[,c("exp.date", "exp.group")], by = c("start_date" = "exp.date")) %>% 
  dplyr::rename(group = exp.group)


dim(pro_sample_info)
dim(exp_sample_info)


##combine data
colnames(exp_expression_data) == exp_sample_info$sample_id

data <-
  t(exp_expression_data) %>% 
  data.frame(., group = exp_sample_info$group, stringsAsFactors = FALSE) %>% 
  plyr::dlply(.variables = .(group)) %>% 
  lapply(function(x){
    x <- 
      x %>% dplyr::select(-group)
    purrr::map(x, .f = function(x){mean(x)}) %>% unlist()
  }) %>% 
  do.call(rbind, .)

rownames(data) == exp_sample_info$group

exp_expression_data <- t(data) %>% as.data.frame()

exp_sample_info <-
  exp_sample_info %>% 
  plyr::dlply(.variables = .(group)) %>% 
  lapply(function(x){
    x[which.max(x$start_date),]
  }) %>% 
  do.call(rbind, .)


dim(exp_expression_data)

dim(pro_expression_data)

save(exp_expression_data, file = "exp_expression_data")
save(pro_expression_data, file = "pro_expression_data")

#######correlation analysis
exp_expression_data <- log(exp_expression_data + 1, 2) %>% 
  apply(1, function(x){
    (x - mean(x))/ sd(x)
  }) %>% 
  t() %>% 
  as.data.frame()

# pro_expression_data <- log(pro_expression_data + 1, 2)
load('cor_value')

cor_value1 <- 
  cor_value %>% 
  dplyr::filter(abs(estimate) > 0.8 & p.value < 0.05)

cor_value1 <- 
  cor_value1 %>% 
  dplyr::left_join(exp_variable_info[,c(1:2)], by = c("exp_name" = "peak_ID")) %>% 
  dplyr::rename(exp_id = MetabID)



pro_num <- 
  cor_value1 %>% 
  dplyr::group_by(pro_name) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::arrange(dplyr::desc(n))


cor_value1 <-
  cor_value1 %>% 
  dplyr::left_join(pro_num, by = c("pro_name")) %>% 
  dplyr::arrange(desc(n), desc(pro_name), desc(estimate))


dir.create("important_pro")
setwd("important_pro/")
####important proteome
important_pro <- 
  cor_value1 %>% 
  dplyr::group_by(pro_name) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::arrange(desc(n)) %>%  
  dplyr::filter(n >= 50)

# important_exp <- 
#   cor_value1 %>% 
#   dplyr::group_by(exp_name) %>% 
#   dplyr::summarise(n = n()) %>% 
#   dplyr::arrange(desc(n)) %>% 
#   dplyr::filter(n >= 20)



library(igraph)
library(ggraph)
library(tidygraph)

##pheatmap
library(pheatmap)

temp_data <- 
  pro_expression_data[important_pro$pro_name,]

colnames(temp_data) <- c(1:5)

# rownames(temp_data) <-
#   exp_variable_info$MetabID[match(rownames(temp_data),
#                                   exp_variable_info$peak_ID)]

plot <- 
  pheatmap(temp_data, scale = "row", 
           show_rownames = TRUE, border_color = "white")


ggsave(plot, filename = "important_proteome_heatmap.png", width = 10, height = 7)

ggsave(plot, filename = "important_proteome_heatmap.pdf", width = 10, height = 7)


temp_data <- 
  exp_expression_data[unique(cor_value1$exp_name[which(cor_value1$pro_name %in% important_pro$pro_name)]),]

colnames(temp_data) <- c(1:5)

annotation_row <-
  data.frame(name = rownames(temp_data), stringsAsFactors = FALSE) %>%
  dplyr::left_join(exp_variable_info[, c(1, 3)], by = c("name" = "peak_ID")) %>% 
  tibble::column_to_rownames(var = "name")

##remove unknown
remove_idx <- which(annotation_row$Chemical_class == "Unknown")
temp_data <- temp_data[-remove_idx,,drop = FALSE]
annotation_row <- annotation_row[-remove_idx,,drop = FALSE]


chemical_class <- 
  annotation_row$Chemical_class %>% unique() %>% 
  sort()

library(ggsci)
col <-
  c(pal_futurama()(12), pal_jama()(3))

names(col) <- chemical_class

annotation_colors <- 
  list(Chemical_class = col)

plot <- 
  pheatmap(temp_data, scale = "row", 
           show_rownames = FALSE, 
           border_color = "white",
           annotation_row = annotation_row, 
           annotation_colors = annotation_colors,
           clustering_method = "ward.D")


ggsave(plot, filename = "important_exposome_heatmap.png", width = 10, height = 7)

ggsave(plot, filename = "important_exposome_heatmap.pdf", width = 10, height = 7)


###correlation network for importtant proteome
idx <- 1

edge_data <-  
  cor_value1 %>% 
  dplyr::filter(pro_name %in% important_pro$pro_name[idx]) %>%
  dplyr::rename(from = exp_name, 
                to = pro_name, 
                Correlation = estimate) %>% 
  dplyr::mutate(p.value = -log(p.value, 10))

node_data <- 
  cor_value1 %>% 
  dplyr::filter(pro_name %in% important_pro$pro_name[idx]) %>%
  dplyr::rename(from = exp_name, to = pro_name) %>% 
  dplyr::select(from, to) %>% 
  tidyr::pivot_longer(cols = c(from, to), 
                      names_to = "class", values_to = "node") %>% 
  dplyr::mutate(class1 = case_when(
    stringr::str_detect(class, "from") ~ "Exposome",
    TRUE ~ "Proteome"
  )) %>% 
  dplyr::select(node, class1) %>% 
  dplyr::rename(Class = class1) %>%
  dplyr::distinct(node, .keep_all = TRUE)

temp_data <- 
  tidygraph::tbl_graph(nodes = node_data, 
                       edges = edge_data,
                       directed = TRUE) %>% 
  dplyr::mutate(Degree = centrality_degree(mode = 'all'))

plot <- 
  ggraph(temp_data,
         layout = 'linear',
         circular = TRUE) +
  geom_edge_arc(
    aes(width = p.value, 
        color = Correlation),
    show.legend = TRUE,
    arrow = arrow(length = unit(2, 'mm'))
  ) +
  geom_node_text(aes(label = ifelse(Class == "Proteome", node, NA))) +
  geom_node_point(aes(color = Class, size = Degree), show.legend = TRUE) +
  ggraph::scale_edge_color_gradientn(colours = c(alpha("#8DD3C7", 1),
                                                 "white",
                                                 alpha("#FB8072", 1))) +
  ggraph::scale_edge_width(range = c(0.2,2)) +
  scale_size_continuous(range = c(1, 4)) +
  scale_color_manual(values = c("Exposome" = "#3B4992FF",
                                "Proteome" = "#008280FF")) +
  theme_void() +
  theme(plot.background = element_rect(fill = "transparent", color = NA), 
        panel.background = element_rect(fill = "transparent", color = NA))

plot


ggsave(plot, filename = paste("important_pro_network",idx,".pdf",sep=""), 
       width = 8, height = 7, bg = "transparent")

ggsave(plot, filename = paste("important_pro_network",idx,".png",sep=""),
       width = 8, height = 7, , bg = "transparent")

important_pro <- 
  apply(important_pro, 1, function(x){
    x <- as.character(x)
    cor <- dplyr::filter(cor_value1, pro_name == x[1]) %>% 
      dplyr::pull(estimate)
    c(sum(cor > 0), sum(cor < 0)) * 100/length(cor)
  }) %>% 
  t() %>% 
  cbind(important_pro, .)

colnames(important_pro)[c(3,4)] <- c("pos", "neg")

important_pro <- 
  important_pro %>% 
  dplyr::mutate(class = case_when(
    pos > 50 ~ "positive dominant",
    TRUE ~ "negative dominant"
  ))

write.csv(important_pro, "important_pro.csv")



###correlation network for important exposome
idx <- 5

edge_data <-  
  cor_value1 %>% 
  dplyr::filter(exp_name %in% important_exp$exp_name[idx]) %>%
  dplyr::rename(from = exp_name, 
                to = pro_name, 
                Correlation = estimate) %>% 
  dplyr::mutate(p.value = -log(p.value, 10))

node_data <- 
  cor_value1 %>% 
  dplyr::filter(exp_name %in% important_exp$exp_name[idx]) %>%
  dplyr::rename(from = exp_name, to = pro_name) %>% 
  dplyr::select(from, to) %>% 
  tidyr::pivot_longer(cols = c(from, to), 
                      names_to = "class", values_to = "node") %>% 
  dplyr::mutate(class1 = case_when(
    stringr::str_detect(class, "from") ~ "Exposome",
    TRUE ~ "Proteome"
  )) %>% 
  dplyr::select(node, class1) %>% 
  dplyr::rename(Class = class1) %>%
  dplyr::distinct(node, .keep_all = TRUE)

temp_data <- 
  tidygraph::tbl_graph(nodes = node_data, 
                       edges = edge_data,
                       directed = TRUE) %>% 
  dplyr::mutate(Degree = centrality_degree(mode = 'all'))

plot <- 
  ggraph(temp_data,
         layout = 'linear',
         circular = TRUE) +
  geom_edge_arc(
    aes(width = p.value, 
        color = Correlation),
    show.legend = TRUE,
    arrow = arrow(length = unit(2, 'mm'))
  ) +
  geom_node_text(aes(label = ifelse(Class == "Proteome", node, NA))) +
  geom_node_text(aes(label = node)) +
  geom_node_point(aes(color = Class, size = Degree), show.legend = TRUE) +
  ggraph::scale_edge_color_gradientn(colours = c(alpha("#8DD3C7", 1),
                                                 "white",
                                                 alpha("#FB8072", 1))) +
  ggraph::scale_edge_width(range = c(0.2,2)) +
  scale_size_continuous(range = c(1, 4)) +
  scale_color_manual(values = c("Exposome" = "#3B4992FF",
                                "Proteome" = "#008280FF")) +
  theme_void() +
  theme(plot.background = element_rect(fill = "transparent", color = NA), 
        panel.background = element_rect(fill = "transparent", color = NA))

plot


ggsave(plot, filename = paste("important_exp_network",idx,".pdf",sep=""), 
       width = 8, height = 7, bg = "transparent")

ggsave(plot, filename = paste("important_exp_network",idx,".png",sep=""),
       width = 8, height = 7, , bg = "transparent")


important_pro <- 
  apply(important_pro, 1, function(x){
    x <- as.character(x)
    cor <- dplyr::filter(cor_value1, pro_name == x[1]) %>% 
      dplyr::pull(estimate)
    c(sum(cor > 0), sum(cor < 0)) * 100/length(cor)
  }) %>% 
  t() %>% 
  cbind(important_pro, .)

colnames(important_pro)[c(3,4)] <- c("pos", "neg")


write.csv(important_pro, "important_pro.csv")


important_exp <- 
  apply(important_exp, 1, function(x){
    x <- as.character(x)
    cor <- dplyr::filter(cor_value1, exp_name == x[1]) %>% 
      dplyr::pull(estimate)
    c(sum(cor > 0), sum(cor < 0)) * 100/length(cor)
  }) %>% 
  t() %>% 
  cbind(important_exp, .)

colnames(important_exp)[c(3,4)] <- c("pos", "neg")


write.csv(important_exp, "important_exp.csv")




#####importance protein and it's exposome
edge_data <-  
  cor_value1 %>% 
  dplyr::filter(pro_name %in% important_pro$pro_name) %>%
  dplyr::rename(from = exp_id ,
                to = pro_name, 
                Correlation = estimate) %>% 
  dplyr::mutate(p.value = -log(p.value, 10))

node_data <- 
  cor_value1 %>% 
  dplyr::filter(pro_name %in% important_pro$pro_name) %>%
  dplyr::rename(from = exp_id, to = pro_name) %>% 
  dplyr::select(from, to) %>% 
  tidyr::pivot_longer(cols = c(from, to), 
                      names_to = "class", values_to = "node") %>% 
  dplyr::mutate(class1 = case_when(
    stringr::str_detect(class, "from") ~ "Exposome",
    TRUE ~ "Proteome"
  )) %>% 
  dplyr::select(node, class1) %>% 
  dplyr::rename(Class = class1) %>%
  dplyr::distinct(node, .keep_all = TRUE) %>% 
  dplyr::left_join(important_pro[,c(1,5)], by = c("node" = "pro_name"))

node_data$class[is.na(node_data$class)] <- "Exposome"

temp_data <- 
  tidygraph::tbl_graph(nodes = node_data, 
                       edges = edge_data,
                       directed = TRUE) %>% 
  dplyr::mutate(Degree = centrality_degree(mode = 'all'))

plot <- 
  ggraph(temp_data,
         layout = 'linear',
         circular = TRUE) +
  geom_edge_arc(
    aes(width = p.value, 
        color = Correlation),
    show.legend = TRUE,
    arrow = arrow(length = unit(2, 'mm'))
  ) +
  geom_node_text(aes(label = ifelse(Class == "Proteome", node, NA))) +
  geom_node_point(aes(color = class, size = Degree), show.legend = TRUE) +
  ggraph::scale_edge_color_gradientn(colours = c(alpha("#8DD3C7", 1),
                                                 "white",
                                                 alpha("#FB8072", 1))) +
  ggraph::scale_edge_width(range = c(0.2,2)) +
  scale_size_continuous(range = c(1, 4)) +
  scale_color_manual(values = c(
    "Exposome" = ggsci::pal_aaas()(1),
    "positive dominant" = ggsci::pal_aaas(alpha = 0.7)(12)[6],
    "negative dominant" = ggsci::pal_aaas(alpha = 0.7)(12)[5]
  )) +
  # ggsci::scale_color_aaas() +
  theme_void() +
  theme(plot.background = element_rect(fill = "transparent", color = NA), 
        panel.background = element_rect(fill = "transparent", color = NA))

plot


ggsave(plot, filename = "important_protein_exposome.pdf", 
       width = 8, height = 7, bg = "transparent")

ggsave(plot, filename ="important_protein_exposome.png", 
       width = 8, height = 7, , bg = "transparent")


important_pro_expos <- 
cor_value1 %>% 
  dplyr::filter(pro_name %in% important_pro$pro_name) %>% 
  dplyr::arrange(pro_name)


write.csv(important_pro_expos, file = "pro_exp_cor.csv")
