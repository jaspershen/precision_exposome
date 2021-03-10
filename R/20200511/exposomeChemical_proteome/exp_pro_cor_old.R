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
  dplyr::mutate(class = factor(class, levels = c("pro", "exp"))) %>% 
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


ggsave(plot, filename = "exp_pro_match.pdf", width = 10, height = 7)
ggsave(plot, filename = "exp_pro_match.png", width = 10, height = 7)
ggsave(plot, filename = "exp_pro_match2.png", width = 14, height = 7)



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

# cor_value <-
# purrr::map(.x = as.data.frame(t(exp_expression_data)), .f = function(y){
#   purrr::map(.x = as.data.frame(t(pro_expression_data)), .f = function(x){
#     test <- cor.test(x, y, method = "pearson")
#     broom::tidy(test)
#   })
# })
# 
# cor_value <-
#   cor_value %>%
#   lapply(function(x){
#     x <- do.call(rbind, x)
#     x$fdr <- p.adjust(x$p.value, method = "fdr")
#     x
#   })
# 
# cor_value <-
#   purrr::map2(cor_value, names(cor_value), .f = function(x, y){
#     x <- data.frame(exp_name = y,
#                     pro_name = pro_variable_info$protein_id,
#                     x[,c(1, 3, 9)], stringsAsFactors = FALSE)
#     x <-
#     x %>%
#       dplyr::filter(abs(estimate) > 0.3)
#     x
#   })
# 
# 
# cor_value <-
#   cor_value %>%
#   do.call(rbind, .)
# 
# cor_value %>%
#   dplyr::filter(estimate > 0.9) %>%
#   dplyr::sample_n(1)
# 
# 
# save(cor_value, file = "cor_value")

load('cor_value')

plot1 <-
  cor_value %>%
  dplyr::filter(exp_name %in% exp_variable_info$peak_ID[1:30]) %>%
  # dplyr::filter(abs(estimate) > 0.7) %>%
  dplyr::group_by(exp_name) %>%
  dplyr::summarise(n = sum(abs(estimate) > 0.5 & p.value < 0.05)) %>%
  ggplot(aes(exp_name, n)) +
  scale_y_continuous(limits = c(-100,50)) +
  geom_bar(stat = "identity", aes(fill = n), show.legend = FALSE) +
  scale_fill_gradientn(colours = c(
    alpha("#8DD3C7", 0.6),
    alpha("#8DD3C7", 1)
    # alpha("#FB8072", 0.4),
    # alpha("#FB8072", 1)
  )) +
  theme_void() +
  theme(
        plot.margin = unit(c(0, 0, 0, 0), "pt"),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA)) +
  coord_polar()

plot1


cor_value1 <-
  cor_value %>%
  dplyr::mutate(p.value = -log(p.value, 10))

plot2 <-
  cor_value1 %>%
  dplyr::filter(exp_name %in% exp_variable_info$peak_ID[1:30]) %>%
  ggplot() +
  geom_hline(yintercept = -log(0.05,10), color = "red") +
  geom_bar(aes(x = x, y = y),
           stat = "identity",
           fill = "transparent",
           color = "#8DD3C7",
           data = data.frame(x = exp_variable_info$peak_ID[1:30],
                             y = 3, stringsAsFactors = FALSE)) +
  geom_jitter(aes(x = exp_name, y = p.value,
                  color = estimate),
              shape = 16,
              show.legend = FALSE, size = 0.5) +
  # scale_y_continuous(limits = c(-2,5)) +
  scale_y_reverse(lim = c(3, 0)) +
  # scale_size_continuous(range = c(0.05,1)) +
  scale_colour_gradientn(colours = c(
    alpha("#8DD3C7", 1),
    # alpha("#8DD3C7", 0.1),
    # # "grey",
    # alpha("#FB8072", 0.1),
    # "red"
    alpha("#FB8072", 1)
  )) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
  ) +
  coord_polar()


plot2


ggsave(plot1, filename = "exp_pro_cor_plot1.pdf",
       width = 7, height = 7, bg = "transparent")
ggsave(plot1, filename = "exp_pro_cor_plot1.png",
       width = 7, height = 7,bg = "transparent")

ggsave(plot2, filename = "exp_pro_cor_plot2.pdf",
       width = 7, height = 7,
       bg = "transparent")
ggsave(plot2, filename = "exp_pro_cor_plot2.png",
       width = 7, height = 7,
       bg = "transparent")




cor_value1 <- 
cor_value %>% 
  dplyr::filter(abs(estimate) > 0.8 & p.value < 0.05)

cor_value1 <- 
  cor_value1 %>% 
  dplyr::left_join(exp_variable_info[,c(1:2)], by = c("exp_name" = "peak_ID")) %>% 
  dplyr::rename(exp_id = MetabID)
  
dir.create("correlation_pro_exp", showWarnings = FALSE)

# for(idx in 1:nrow(cor_value1)){
#   cat(idx, " ")
#   path1 <- file.path("correlation_pro_exp",
#                      cor_value1$exp_name[idx])
#   dir.create(path1, showWarnings = FALSE)
#   temp_data <-
#     data.frame(date = as.character(pro_sample_info$CollectionDate),
#                exp = as.numeric(exp_expression_data[cor_value1$exp_name[idx],]),
#                pro = as.numeric(pro_expression_data[cor_value1$pro_name[idx],]),
#                stringsAsFactors = FALSE)
#   plot <-
#   temp_data %>%
#     ggplot(aes(exp, pro)) +
#     geom_point() +
#     geom_smooth(method = "lm", color = "skyblue") +
#     ggrepel::geom_label_repel(aes(x = exp, pro, label = date)) +
#     labs(x = paste("Exposome: ", cor_value1$exp_id[idx], sep = ""),
#          y = paste("Proteome: " ,cor_value1$pro_name[idx]), sep = "") +
#     theme_bw() +
#     theme(axis.title = element_text(size = 13),
#           axis.text = element_text(size = 12),
#           plot.background = element_rect(fill = "transparent", color = NA),
#           panel.background = element_rect(fill = "transparent", color = NA)) +
#     annotate(geom = "text",
#              x = -Inf, y = Inf,
#              label = paste("Correlation: ",round(cor_value1$estimate[idx], 2),
#                            "\nP value: ", cor_value1$p.value[idx], sep = ""),
#              vjust = 2, hjust = -1)
# 
#   name <- paste(cor_value1$exp_name[idx], "_",
#                 cor_value1$pro_name[idx], ".png", sep = "")
# 
#   ggsave(plot, filename = file.path(path1, name),
#          width = 7, height = 7, bg = "transparent")
# 
# }

library(pheatmap)

cor_value1$exp_name %>% unique()

# for(i in 1:length(unique(cor_value1$exp_name))){
#   cat(i, " ")
#   temp_data <-
#     cor_value1 %>%
#     dplyr::filter(exp_name == unique(cor_value1$exp_name)[i])
# 
#   temp_data1 <-
#     exp_expression_data[temp_data$exp_name[1],]
# 
#   temp_data2 <-
#     pro_expression_data[temp_data$pro_name,]
# 
#   colnames(temp_data2) <- colnames(temp_data1)
# 
#   plot <-
#   pheatmap(rbind(temp_data1, temp_data2),
#            scale = "row", border_color = "white",
#            cluster_cols = FALSE)
# 
#   name <- file.path("correlation_pro_exp", unique(cor_value1$exp_name)[i],
#                     paste(unique(cor_value1$exp_name)[i], "_heatmap.png", sep = ""))
# 
#   ggsave(plot, filename = name,
#          width = 7, height = 7, bg = "transparent")
# 
# }


pro_num <- 
cor_value1 %>% 
  dplyr::group_by(pro_name) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::arrange(dplyr::desc(n))


cor_value1 <-
  cor_value1 %>% 
  dplyr::left_join(pro_num, by = c("pro_name")) %>% 
  dplyr::arrange(desc(n), desc(pro_name), desc(estimate))


write.csv(cor_value1[,-5], "cor_value.csv", row.names = FALSE)


##cluster exposome first
###k-means
## fuzzy c-means clustring
###fuzzy k-means clustering
library(Mfuzz)

#first get the time point data together:
# bind that to the data frame
temp_data <- exp_expression_data

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
tmp <- tempfile()

write.table(temp_data,file=tmp, sep='\t', quote = F, col.names=NA)

#read it back in as an expression set
data <- table2eset(file=tmp)

data.s <- standardise(data)


m1 <- mestimate(data.s)
m1


Dmin(data.s, m=m1, crange=seq(2,22,1), repeats=3, visu=TRUE)

clust = 2

c <- mfuzz(data.s,c=clust,m=m1)


mfuzz.plot(eset = data.s, 
           # min.mem = 0.6,
           cl = c, 
           # mfrow=c(1,1),
           time.labels = time,
           new.window = FALSE)

cluster <- 
  data.frame(peak = names(c$cluster), 
             cluster = c$cluster, 
             stringsAsFactors = FALSE) %>% 
  arrange(cluster)


write.csv(cluster, "exposome_cluster/cluster.csv")

cluster <- readr::read_csv("exposome_cluster/cluster.csv")

###pathway analysis for each exposome chemical PIUpro
cluster1 <- 
  cluster %>% 
  dplyr::filter(cluster == 1) %>% 
  dplyr::pull(peak)

cluster2 <- 
  cluster %>% 
  dplyr::filter(cluster == 2) %>% 
  dplyr::pull(peak)



###correlation network for pro and exp
cor_value1$exp_name %>% unique() %>% length()
cor_value1$pro_name %>% unique() %>% length()


library(igraph)
library(ggraph)
library(tidygraph)

###network for cluster 1
edge_data <-  
  cor_value1 %>% 
  dplyr::filter(exp_name %in% cluster1) %>%
  dplyr::rename(from = exp_name, 
                to = pro_name, 
                Correlation = estimate) %>% 
  dplyr::mutate(p.value = -log(p.value, 10))

node_data <- 
  cor_value1 %>% 
  dplyr::filter(exp_name %in% cluster1) %>%
  dplyr::rename(from = exp_name, to = pro_name) %>% 
  dplyr::select(from, to) %>% 
  tidyr::pivot_longer(cols = c(from, to), 
                      names_to = "class", values_to = "node") %>% 
  dplyr::mutate(class1 = case_when(
    stringr::str_detect(class, "from") ~ "Exposome",
    TRUE ~ "proteome"
  )) %>% 
  dplyr::select(node, class1) %>% 
  dplyr::rename(Class = class1) %>%
  dplyr::distinct(node, .keep_all = TRUE)

temp_data <- 
  tidygraph::tbl_graph(nodes = node_data, 
                       edges = edge_data,
                       directed = TRUE) %>% 
  dplyr::mutate(Degree = centrality_degree(mode = 'all'))

plot1 <- 
ggraph(temp_data,
       layout = 'linear',
       circular = TRUE) +
  geom_edge_arc(
    aes(width = p.value, 
        color = Correlation),
    show.legend = TRUE,
    arrow = arrow(length = unit(2, 'mm'))
  ) +
  geom_node_point(aes(color = Class, size = Degree), show.legend = TRUE) +
  ggraph::scale_edge_color_gradientn(colours = c(alpha("#8DD3C7", 1),
                                                 "white",
                                                 alpha("#FB8072", 1))) +
  ggraph::scale_edge_width(range = c(0.2,2)) +
  scale_size_continuous(range = c(1, 4)) +
  # ggsci::scale_color_aaas() +
  scale_color_manual(values = c("Exposome" = "#3B4992FF",
                                "proteome" = "#008280FF")) +
  theme_void() +
  theme(plot.background = element_rect(fill = "transparent", color = NA), 
        panel.background = element_rect(fill = "transparent", color = NA))

plot1


ggsave(plot1, filename = "cluster1_network.pdf",
       width = 8, height = 7, bg = "transparent")

ggsave(plot1, filename = "cluster1_network.png",
       width = 8, height = 7, , bg = "transparent")

###network for cluster 2
edge_data <-  
  cor_value1 %>% 
  dplyr::filter(exp_name %in% cluster2) %>%
  dplyr::rename(from = exp_name, 
                to = pro_name, 
                Correlation = estimate) %>% 
  dplyr::mutate(p.value = -log(p.value, 10))

node_data <- 
  cor_value1 %>% 
  dplyr::filter(exp_name %in% cluster2) %>%
  dplyr::rename(from = exp_name, to = pro_name) %>% 
  dplyr::select(from, to) %>% 
  tidyr::pivot_longer(cols = c(from, to), 
                      names_to = "class", values_to = "node") %>% 
  dplyr::mutate(class1 = case_when(
    stringr::str_detect(class, "from") ~ "Exposome",
    TRUE ~ "proteome"
  )) %>% 
  dplyr::select(node, class1) %>% 
  dplyr::rename(Class = class1) %>%
  dplyr::distinct(node, .keep_all = TRUE)

temp_data <- 
  tidygraph::tbl_graph(nodes = node_data, 
                       edges = edge_data,
                       directed = TRUE) %>% 
  dplyr::mutate(Degree = centrality_degree(mode = 'all'))

plot2 <- 
  ggraph(temp_data,
         layout = 'linear',
         circular = TRUE) +
  geom_edge_arc(
    aes(width = p.value, 
        color = Correlation),
    show.legend = TRUE,
    arrow = arrow(length = unit(2, 'mm'))
  ) +
  geom_node_point(aes(color = Class, size = Degree), show.legend = TRUE) +
  ggraph::scale_edge_color_gradientn(colours = c(alpha("#8DD3C7", 1),
                                                 "white",
                                                 alpha("#FB8072", 1))) +
  ggraph::scale_edge_width(range = c(0.2,2)) +
  scale_size_continuous(range = c(1, 4)) +
  scale_color_manual(values = c("Exposome" = "#3B4992FF",
                                "proteome" = "#008280FF")) +
  theme_void() +
  theme(plot.background = element_rect(fill = "transparent", color = NA), 
        panel.background = element_rect(fill = "transparent", color = NA))

plot2


ggsave(plot2, filename = "cluster2_network.pdf",
       width = 8, height = 7, bg = "transparent")

ggsave(plot2, filename = "cluster2_network.png",
       width = 8, height = 7, , bg = "transparent")



###each exposome chemical for pathway
cluster1_pro <- 
cor_value1 %>% 
  dplyr::filter(exp_name %in% cluster1) %>%
  dplyr::select(pro_name) %>% 
  dplyr::distinct(pro_name, .keep_all = TRUE)


cluster2_pro <- 
  cor_value1 %>% 
  dplyr::filter(exp_name %in% cluster2) %>%
  dplyr::select(pro_name) %>% 
  dplyr::distinct(pro_name, .keep_all = TRUE) 


###GSEA enrichment
library(clusterProfiler)
temp_data <- 
  pro_expression_data[unique(cor_value1$pro_name),]

correlation <- 
temp_data %>% 
  t() %>% 
  as.data.frame() %>% 
  purrr::map(function(x){
    cor(x, 1:5)
  }) %>% 
  unlist() 

correlation <- sort(correlation, decreasing = TRUE)

protein_list <-
  data.frame(protein_id = names(correlation),
             correlation,
             stringsAsFactors = FALSE)

protein_list_entrezid <-
  clusterProfiler::bitr(geneID = protein_list$protein_id, 
                        toType = "ENTREZID",
                        OrgDb = "org.Hs.eg.db",
                        fromType = "SYMBOL", 
                        drop = TRUE) %>% 
  dplyr::left_join(protein_list, by = c("SYMBOL" = "protein_id")) %>% 
  dplyr::distinct(ENTREZID, .keep_all = TRUE) 


protein_list_entrezid_name <- protein_list_entrezid$ENTREZID
protein_list_entrezid <- protein_list_entrezid$correlation

names(protein_list_entrezid) <- protein_list_entrezid_name



protein_list_uniprot <-
  clusterProfiler::bitr(geneID = protein_list$protein_id, 
                        toType = "UNIPROT",
                        OrgDb = "org.Hs.eg.db",
                        fromType = "SYMBOL", 
                        drop = TRUE) %>% 
  dplyr::left_join(protein_list, by = c("SYMBOL" = "protein_id")) %>% 
  dplyr::distinct(UNIPROT, .keep_all = TRUE)

protein_list_uniprot_name <- protein_list_uniprot$UNIPROT
protein_list_uniprot <- protein_list_uniprot$correlation

names(protein_list_uniprot) <- protein_list_uniprot_name


##GSEA GO
gse_kegg <-
  clusterProfiler::gseKEGG(
    organism = "hsa",
    verbose = TRUE,
    use_internal_data = FALSE,
    geneList = protein_list_uniprot,
    keyType = "uniprot",
    pvalueCutoff = 0.2,
    pAdjustMethod = "fdr",
    seed = TRUE
  )


save(gse_kegg, file = "gse_kegg")
load("gse_kegg")
gse_kegg@result$Description

idx <- 7
plot <- 
clusterProfiler::gseaplot(x = gse_kegg, 
                          geneSetID = idx, 
                          color.line = "#BB0021FF",
                          title = gse_kegg@result$Description[idx])

ggsave(plot, filename = paste(gse_kegg@result$ID[idx], ".png", sep = ""),
       width = 8, height = 7)

ggsave(plot, filename = paste(gse_kegg@result$ID[idx], ".pdf", sep = ""),
       width = 8, height = 7)


plot <- 
clusterProfiler::dotplot(gse_kegg)
ggsave(plot, filename = "gse_kegg_dotplot.pdf", width = 9, height = 7)
ggsave(plot, filename = "gse_kegg_dotplot.png", width = 9, height = 7)

plot <- 
clusterProfiler::cnetplot(gse_kegg)
plot
ggsave(plot, filename = "gse_kegg_cnetplot.pdf", width = 9, height = 7)
ggsave(plot, filename = "gse_kegg_cnetplot.png", width = 9, height = 7)

plot <- 
clusterProfiler::emapplot(gse_kegg)
plot


clusterProfiler::heatplot(gse_kegg)

clusterProfiler::ridgeplot(gse_kegg)

##GSEA GO
gse_go <-
  clusterProfiler::gseGO(
    geneList     = protein_list_entrezid,
    OrgDb        = org.Hs.eg.db,
    ont          = "ALL",
    keyType = "ENTREZID",
    # nPerm        = 1000,
    # minGSSize    = 100,
    # maxGSSize    = 500,
    pvalueCutoff = 0.5,
    pAdjustMethod = "fdr",
    verbose      = TRUE
  )


save(gse_go, file = "gse_go")
load("gse_go")

gse_go@result$Description

idx <- 13
plot <- 
  clusterProfiler::gseaplot(x = gse_go, 
                            geneSetID = idx, 
                            color.line = "#BB0021FF",
                            title = gse_go@result$Description[idx])
plot

ggsave(plot, filename = paste(gse_go@result$ID[idx], ".png", sep = ""),
       width = 8, height = 7)

ggsave(plot, filename = paste(gse_go@result$ID[idx], ".pdf", sep = ""),
       width = 8, height = 7)


plot <- 
  clusterProfiler::dotplot(gse_go)
plot
ggsave(plot, filename = "gse_go_dotplot.pdf", width = 9, height = 7)
ggsave(plot, filename = "gse_go_dotplot.png", width = 9, height = 7)

plot <- 
  clusterProfiler::cnetplot(gse_go)
plot
ggsave(plot, filename = "gse_go_cnetplot.pdf", width = 9, height = 7)
ggsave(plot, filename = "gse_go_cnetplot.png", width = 9, height = 7)


clusterProfiler::dotplot(gse_go)

clusterProfiler::cnetplot(gse_go)

clusterProfiler::emapplot(gse_go)

clusterProfiler::heatplot(gse_go)

clusterProfiler::ridgeplot(gse_go)




####important proteome
important_pro <- 
cor_value1 %>% 
  dplyr::group_by(pro_name) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::arrange(desc(n)) %>% 
  dplyr::filter(n >= 50)

important_exp <- 
  cor_value1 %>% 
  dplyr::group_by(exp_name) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::arrange(desc(n)) %>% 
  dplyr::filter(n >= 20)



library(igraph)
library(ggraph)
library(tidygraph)

##pheatmap
library(pheatmap)

temp_data <- 
  pro_expression_data[important_pro$pro_name,]

colnames(temp_data) <- c(1:5)

plot <- 
  pheatmap(temp_data, scale = "row", 
           show_rownames = TRUE, border_color = "white")

ggsave(plot, filename = "important_proteome_heatmap.png", width = 7, height = 7)

ggsave(plot, filename = "important_proteome_heatmap.pdf", width = 7, height = 7)


temp_data <- 
  exp_expression_data[unique(cor_value1$exp_name[which(cor_value1$pro_name %in% important_pro$pro_name)]),]

colnames(temp_data) <- c(1:5)

rownames(temp_data) <-
  exp_variable_info$MetabID[match(rownames(temp_data),
                                  exp_variable_info$peak_ID)]

plot <- 
  pheatmap(temp_data, scale = "row", 
           show_rownames = FALSE, border_color = NA)

ggsave(plot, filename = "important_exposome_heatmap.png", width = 7, height = 7)

ggsave(plot, filename = "important_exposome_heatmap.pdf", width = 7, height = 7)


###correlation network for importtant proteome
idx <- 11

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




