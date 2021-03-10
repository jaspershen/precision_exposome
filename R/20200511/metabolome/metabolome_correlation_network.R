##load data
sxtTools::setwd_project()
library(tidyverse)
rm(list = ls())
####load metabolomics data
sxtTools::setwd_project()
load("data_20200511/metabolome/expression_data")
met_expression_data <- expression_data
load("data_20200511/metabolome/sample_info")
met_sample_info <- sample_info
load("data_20200511/metabolome/variable_info")
met_variable_info <- variable_info

met_annotation_info <- readr::read_csv("data_20200511/metabolome/annotation_table.csv")

met_annotation_info %>% dim()

variable_info <-
  variable_info %>% 
  left_join(met_annotation_info, by = c("peak_name" = "name"))

load("data_20200511/exposome/expression_data")
exp_expression_data <- expression_data
load("data_20200511/exposome/sample_info")
exp_sample_info <- sample_info
load("data_20200511/exposome/variable_info")
exp_variable_info <- variable_info

setwd("data_analysis/met_exp/")

dim(exp_variable_info)

dim(met_variable_info)

temp1 <- 
  exp_sample_info$start_date %>% 
  as.Date() %>% 
  data.frame(date = ., 
             exp = 1,
             stringsAsFactors = FALSE)

temp2 <- 
  met_sample_info$CollectionDate %>% 
  as.Date() %>% 
  data.frame(date = ., 
             met = 1,
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
  lapply(met_sample_info$CollectionDate %>% as.Date(), function(x){
    x <- x - as.Date(exp_sample_info$start_date) 
    as.Date(exp_sample_info$start_date) [which(x >=0 & x <= 7)]
  }) 

names(diff) <- met_sample_info$CollectionDate %>% as.character()

diff <- 
  purrr::map2(diff, names(diff), .f = function(x,y){
    if(length(x) == 0){
      return(NULL)
    }
    z <- data.frame(y, x, stringsAsFactors = FALSE)
    colnames(z) <- c("met", "exp")
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
         met = as.character(met)) %>% 
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

diff$met.date <- as.Date(diff$met.date)

diff$exp.group <- as.character(diff$exp.group)


plot <- 
  temp %>% 
  mutate(date = as.Date(date)) %>% 
  dplyr::mutate(class = factor(class, levels = c("met", "exp"))) %>% 
  ggplot(aes(date, class)) +
  geom_point(aes(color = group),
             alpha = 0.6, size = 4,
             shape = 16, show.legend = FALSE) +
  geom_segment(aes(x = exp.date, 
                   xend = met.date, 
                   y = exp.class, 
                   yend = met.class,
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
  # scale_y_continuous(breaks = c("met", "exp"), labels = c("Met", "Exp")) +
  scale_y_discrete(breaks = c("met", "exp"), labels = c("Metabolome", "Exposome")) +
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
dim(met_expression_data)
dim(met_sample_info)

met_sample_info <-
  met_sample_info %>% 
  dplyr::filter(CollectionDate %in% unique(diff$met.date))

met_expression_data <- 
  met_expression_data %>% 
  dplyr::select(one_of(met_sample_info$sample_id))

met_sample_info <- 
  met_sample_info %>% 
  dplyr::left_join(diff[,c("met.date", "met.group")] %>% dplyr::distinct(met.date, met.group), 
                   by = c("CollectionDate" = "met.date")) %>% 
  dplyr::rename(group = met.group)


exp_sample_info <-
  exp_sample_info %>% 
  dplyr::filter(as.Date(start_date) %in% diff$exp.date)


exp_expression_data <-
  exp_expression_data %>%
  select(-contains("Blank")) %>% 
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


dim(met_sample_info)
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

dim(met_expression_data)


#######correlation analysis
exp_expression_data <- log(exp_expression_data + 1, 2)
met_expression_data <- log(met_expression_data + 1, 2)


dir.create("met_correlation_network", showWarnings = FALSE)
setwd("met_correlation_network")

internal_exposome <- readxl::read_xlsx("Internal exposome.xlsx", sheet = 2)



met_variable_info <- 
  met_variable_info %>% 
  dplyr::left_join(internal_exposome[,c(1:2)], by = c("peak_name"))

met_variable_info$Chemical_class[is.na(met_variable_info$Chemical_class)] <- "Other"

table(met_variable_info$Chemical_class)



####only for internal exposome
temp_data <-  met_expression_data %>% 
  tibble::rownames_to_column(var = "name") %>% 
  dplyr::filter(name %in% met_variable_info$peak_name[met_variable_info$Chemical_class != "Other"]) %>%
  tibble::column_to_rownames(var = "name")

cor_value <- 
  cor(as.matrix(t(temp_data)), method = "spearman")

cor.test(as.numeric(temp_data[1,]), as.numeric(temp_data[20,]))

# cor_value <- 
#   Hmisc::rcorr(as.matrix(t(temp_data)), type = "spearman")
  # cor(as.matrix(t(temp_data)), method = "spearman")

cor_value[lower.tri(cor_value)] <- 0

cor_value <- 
  cor_value %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "from") %>% 
  tidyr::pivot_longer(-from, names_to = "to", values_to = "cor") %>% 
  dplyr::filter(abs(cor) > 0.3) %>% 
  dplyr::filter(cor != 1)

p_value <- 
  purrr::map(as.data.frame(t(cor_value)), .f = function(x){
    x <- as.character(x)
    cor.test(as.numeric(met_expression_data[x[1],]),
             as.numeric(met_expression_data[x[2],]), 
             method = "spearman")$p.value
  }) %>% 
  unlist()

p_adjust <- p.adjust(p_value, method = "fdr")


cor_value <- data.frame(cor_value, p_value, p_adjust, stringsAsFactors = FALSE)

# cor_value <- 
  # cor_value %>% 
  # dplyr::filter(cor > 0.8 & p_value < 0.05)
# dplyr::filter(cor > 0.3 & p_value < 0.1)

cor_value <- 
  cor_value %>% 
  dplyr::left_join(met_variable_info[,c(1,3)], by = c("from" = "peak_name")) %>% 
  dplyr::rename(from_id = Metabolite) %>% 
  dplyr::left_join(met_variable_info[,c(1,3)], by = c("to" = "peak_name")) %>% 
  dplyr::rename(to_id = Metabolite) 

library(igraph)
library(ggraph)
library(tidygraph)

edge_data <-  
  cor_value %>% 
  dplyr::select(-c(from, to)) %>% 
  dplyr::rename(from = from_id, 
                to = to_id, 
                Correlation = cor) %>% 
  dplyr::mutate(p.value = -log(p_value, 10))

node_data <- 
  cor_value %>% 
  dplyr::select(-c(from, to)) %>% 
  dplyr::rename(from = from_id, to = to_id) %>% 
  dplyr::select(from, to) %>% 
  tidyr::pivot_longer(cols = c(from, to), 
                      names_to = "class", values_to = "node") %>% 
  dplyr::select(node) %>% 
  dplyr::left_join(met_variable_info, by = c("node" = "Metabolite")) %>% 
  dplyr::distinct(node, .keep_all = TRUE)

# node_data <- 
#   node_data[node_data$Chemical_class != "Unknown",]

# edge_data <- 
#   edge_data %>% 
#   dplyr::filter((from %in% node_data$node) & (to %in% node_data$node))

temp_data <- 
  tidygraph::tbl_graph(nodes = node_data, 
                       edges = edge_data,
                       directed = FALSE) %>% 
  dplyr::mutate(Degree = centrality_degree(mode = 'all'))


chemical_class <- 
  node_data$Chemical_class %>% unique() %>% 
  sort()

library(ggsci)

col <-
  c(pal_futurama()(12), pal_jama()(3))

col <-
  pal_npg()(6)

names(col) <- chemical_class

plot <- 
  ggraph(temp_data,
         layout = 'kk',
         circular = FALSE) +
  geom_edge_link(
    aes(width = p.value, 
        color = Correlation),
    show.legend = TRUE,
    arrow = arrow(length = unit(2, 'mm'))
  ) +
  geom_node_point(aes(color = Chemical_class, size = Degree), show.legend = TRUE) +
  ggraph::scale_edge_color_gradientn(colours = c(alpha("#8DD3C7", 1),
                                                 "white",
                                                 alpha("#FB8072", 1))) +
  ggraph::scale_edge_width(range = c(0.2,1.5)) +
  scale_size_continuous(range = c(1, 8)) +
  scale_color_manual(values = col) +
  # ggsci::scale_color_aaas() +
  theme_void() +
  theme(plot.background = element_rect(fill = "transparent", color = NA), 
        panel.background = element_rect(fill = "transparent", color = NA))

plot

ggsave(plot, filename = "met_correlation_network.png", 
       width = 10, height = 7, bg = "transparent")

ggsave(plot, filename = "met_correlation_network.pdf", 
       width = 10, height = 7, bg = "transparent")


save(cor_value, file = "cor_value")
write.csv(cor_value, "cor_value.csv", row.names = FALSE)
