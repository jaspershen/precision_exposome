##avoid source
no_function()

##load data
sxtTools::setwd_project()
library(tidyverse)
rm(list = ls())

####load blood_test data
sxtTools::setwd_project()
load("data_20200511/lab_test/expression_data")
blood_test_expression_data <- expression_data
load("data_20200511/lab_test/sample_info")
blood_test_sample_info <- sample_info
load("data_20200511/lab_test/variable_info")
blood_test_variable_info <- variable_info

###exposome biological
load("data_20200511/microbiome/dna_expression_data")
exposomeBiological_expression_data <- dna_expression_data
load("data_20200511/microbiome/dna_sample_info")
exposomeBiological_sample_info <- dna_sample_info
load("data_20200511/microbiome/dna_variable_info")
exposomeBiological_variable_info <- dna_variable_info

setwd("data_analysis/exposomeBiological_blood_test/")

dim(exposomeBiological_variable_info)

dim(blood_test_variable_info)

temp1 <- 
  exposomeBiological_sample_info$date.start %>% 
  as.Date() %>% 
  data.frame(date = ., 
             exposomeBiological = 1,
             stringsAsFactors = FALSE)

temp2 <- 
  blood_test_sample_info$date.start %>% 
  as.Date() %>% 
  data.frame(date = ., 
             blood_test  = 1,
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
  lapply(blood_test_sample_info$date.start %>% as.Date(), function(x){
    x <- x - as.Date(exposomeBiological_sample_info$date.start) 
    as.Date(exposomeBiological_sample_info$date.start) [which(x >=0 & x <= 2)]
  }) 

names(diff) <- blood_test_sample_info$date.start %>% as.character()

diff <- 
  purrr::map2(diff, names(diff), .f = function(x,y){
    if(length(x) == 0){
      return(NULL)
    }
    z <- data.frame(y, x, stringsAsFactors = FALSE)
    colnames(z) <- c("blood_test", "exposomeBiological")
    z
  }) 

diff <- diff[lapply(diff, is.null) %>% unlist() %>% `!`]

diff <- purrr::map2(.x = diff, .y = 1:length(diff), .f = function(x,y){
  data.frame(x, group = y, stringsAsFactors = FALSE)
}) %>% 
  do.call(rbind, .)

rownames(diff) <- NULL

diff <- 
  diff %>% 
  mutate(exposomeBiological = as.character(exposomeBiological), 
         blood_test  = as.character(blood_test )) %>% 
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

diff$exposomeBiological.date <- as.Date(diff$exposomeBiological.date)

diff$blood_test.date <- as.Date(diff$blood_test.date)

diff$exposomeBiological.group <- as.character(diff$exposomeBiological.group)

plot <- 
  temp %>% 
  mutate(date = as.Date(date)) %>% 
  dplyr::mutate(class = factor(class, levels = c("blood_test", "exposomeBiological"))) %>% 
  ggplot(aes(date, class)) +
  geom_point(aes(fill = group),
             color = "black",
             alpha = 0.6, size = 4,
             shape = 21, show.legend = FALSE) +
  geom_segment(aes(x = exposomeBiological.date, 
                   xend = blood_test.date, 
                   y = exposomeBiological.class, 
                   yend = blood_test.class,
                   color = exposomeBiological.group), data = diff, show.legend = FALSE) +
  scale_fill_manual(values = c("No" = "grey", 
                               "1" = ggsci::pal_aaas()(8)[1],
                               "2" = ggsci::pal_aaas()(8)[2],
                               "3" = ggsci::pal_aaas()(8)[3],
                               "4" = ggsci::pal_aaas()(8)[4],
                               "5" = ggsci::pal_aaas()(8)[5])) +
  scale_x_continuous(trans = "date",
                     breaks = c(as.Date(temp$date)),
                     labels = as.character(temp$date)
  ) +
  # scale_y_continuous(breaks = c("blood_test ", "exposomeBiological"), labels = c("blood_test ", "Exp")) +
  scale_y_discrete(breaks = c("blood_test ", "exposomeBiological"), labels = c("Blood test", "Exposome biological")) +
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

ggsave(plot, filename = "exposomeBiological_blood_test_match.pdf", width = 10, height = 7)
ggsave(plot, filename = "exposomeBiological_blood_test_match2.pdf", width = 14, height = 7)


#-------------------------------------------------------------------------------
###prepare data
dim(blood_test_expression_data)
dim(blood_test_sample_info)

blood_test_sample_info <-
  blood_test_sample_info %>% 
  dplyr::filter(date.start %in% unique(diff$blood_test.date))

blood_test_expression_data <- 
  blood_test_expression_data %>% 
  dplyr::select(one_of(blood_test_sample_info$sample_id))

blood_test_sample_info <- 
  blood_test_sample_info %>% 
  dplyr::left_join(diff[,c("blood_test.date", "blood_test.group")] %>% 
                     dplyr::distinct(blood_test.date, blood_test.group), 
                   by = c("date.start" = "blood_test.date")) %>% 
  dplyr::rename(group = blood_test.group)

exposomeBiological_sample_info <-
  exposomeBiological_sample_info %>% 
  dplyr::filter(as.Date(date.start) %in% diff$exposomeBiological.date)

exposomeBiological_expression_data <-
  exposomeBiological_expression_data %>%
  select(-contains("Blank")) %>% 
  t()

exposomeBiological_expression_data <-
  exposomeBiological_expression_data %>%
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

exposomeBiological_expression_data <- 
  exposomeBiological_expression_data %>% 
  dplyr::select(one_of(exposomeBiological_sample_info$sample_id))

exposomeBiological_sample_info <- 
  exposomeBiological_sample_info %>% 
  mutate(date.start = as.Date(date.start)) %>% 
  dplyr::left_join(diff[,c("exposomeBiological.date", "exposomeBiological.group")], by = c("date.start" = "exposomeBiological.date")) %>% 
  dplyr::rename(group = exposomeBiological.group)

dim(blood_test_sample_info)
dim(exposomeBiological_sample_info)

##combine data
colnames(exposomeBiological_expression_data) == exposomeBiological_sample_info$sample_id

data <-
  t(exposomeBiological_expression_data) %>% 
  data.frame(., group = exposomeBiological_sample_info$group, stringsAsFactors = FALSE) %>% 
  plyr::dlply(.variables = .(group)) %>% 
  lapply(function(x){
    x <- 
      x %>% dplyr::select(-group)
    purrr::map(x, .f = function(x){mean(x)}) %>% unlist()
  }) %>% 
  do.call(rbind, .)

rownames(data) == exposomeBiological_sample_info$group

exposomeBiological_expression_data <- t(data) %>% as.data.frame()

exposomeBiological_sample_info <-
  exposomeBiological_sample_info %>% 
  plyr::dlply(.variables = .(group)) %>% 
  lapply(function(x){
    x[which.max(x$date.start),]
  }) %>% 
  do.call(rbind, .)

dim(exposomeBiological_expression_data)

dim(blood_test_expression_data)

save(exposomeBiological_expression_data, file = "exposomeBiological_expression_data")
save(blood_test_expression_data, file = "blood_test_expression_data")

save(exposomeBiological_sample_info, file = "exposomeBiological_sample_info")
save(blood_test_sample_info, file = "blood_test_sample_info")
save(exposomeBiological_variable_info, file = "exposomeBiological_variable_info")

