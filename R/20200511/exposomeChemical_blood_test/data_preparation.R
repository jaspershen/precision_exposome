##avoid source
no_function()

##load data
sxtTools::setwd_project()
library(tidyverse)
rm(list = ls())

####load lab_test
sxtTools::setwd_project()
load("data_20200511/lab_test/expression_data")
lab_test_expression_data <- expression_data
load("data_20200511/lab_test/sample_info")
lab_test_sample_info <- sample_info
load("data_20200511/lab_test/variable_info")
lab_test_variable_info <- variable_info

###exposome
load("data_20200511/exposome/expression_data")
exp_expression_data <- expression_data
load("data_20200511/exposome/sample_info")
exp_sample_info <- sample_info
load("data_20200511/exposome/variable_info")
exp_variable_info <- variable_info

setwd("data_analysis/exposomeChemical_blood_test/")

dim(exp_variable_info)

dim(lab_test_variable_info)

temp1 <- 
exp_sample_info$start_date %>% 
  as.Date() %>% 
  data.frame(date = ., 
             exp = 1,
             stringsAsFactors = FALSE)

temp2 <- 
lab_test_sample_info$date.start %>% 
  as.Date() %>% 
  data.frame(date = ., 
             lab_test = 1,
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
lapply(lab_test_sample_info$date.start %>% as.Date(), function(x){
  x <- abs(x - as.Date(exp_sample_info$start_date))
  as.Date(exp_sample_info$start_date)[which(x >=0 & x <= 0)]
}) 

names(diff) <- lab_test_sample_info$date.start %>% as.character()

####remove some duplicated samples
# diff$`2016-01-12` <- diff$`2016-01-12`[1]
# diff$`2016-01-13` <- diff$`2016-01-13`[2]
# diff$`2016-01-15` <- diff$`2016-01-15`[3]
# diff$`2016-01-16` <- diff$`2016-01-16`[4]
# diff$`2016-01-17` <- diff$`2016-01-17`[4]
# diff$`2016-02-24` <- diff$`2016-02-24`[3]
# 
# diff <- diff[-14]

diff <- 
purrr::map2(diff, names(diff), .f = function(x,y){
  if(length(x) == 0){
    return(NULL)
  }
  z <- data.frame(y, x, stringsAsFactors = FALSE)
  colnames(z) <- c("lab_test", "exp")
  z
}) 

diff <- diff[lapply(diff, is.null) %>% unlist() %>% `!`]

# diff[[2]] <- diff[[2]][-c(1),]
# 
# diff[[3]] <- diff[[3]][-c(1),]

diff <- purrr::map2(.x = diff, .y = 1:length(diff), .f = function(x,y){
  data.frame(x, group = y, stringsAsFactors = FALSE)
}) %>% 
  do.call(rbind, .)

rownames(diff) <- NULL

diff <- 
diff %>% 
  mutate(exp = as.character(exp), 
         lab_test = as.character(lab_test)) %>% 
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

diff$lab_test.date <- as.Date(diff$lab_test.date)

diff$exp.group <- as.character(diff$exp.group)

plot <- 
temp %>% 
  mutate(date = as.Date(date)) %>% 
  dplyr::mutate(class = factor(class, levels = c("lab_test", "exp"))) %>% 
  ggplot(aes(date, class)) +
  geom_point(aes(fill = group),
             color = "black",
             alpha = 0.6, size = 4,
             shape = 21, show.legend = FALSE) +
  geom_segment(aes(x = exp.date, 
                   xend = lab_test.date, 
                   y = exp.class, 
                   yend = lab_test.class,
                   color = exp.group), data = diff, show.legend = FALSE) +
  scale_fill_manual(values = c("No" = "grey", 
                                "1" = ggsci::pal_aaas()(11)[1],
                                "2" = ggsci::pal_aaas()(11)[2],
                                "3" = ggsci::pal_aaas()(11)[3],
                                "4" = ggsci::pal_aaas()(11)[4],
                                "5" = ggsci::pal_aaas()(11)[5],
                               "6" = ggsci::pal_aaas()(11)[6],
                               "7" = ggsci::pal_aaas()(11)[7],
                               "8" = ggsci::pal_aaas()(11)[8],
                               "9" = ggsci::pal_aaas()(11)[9],
                               "10" = ggsci::pal_aaas()(11)[10])) +
  scale_x_continuous(trans = "date",
                     breaks = c(as.Date(temp$date)),
                     labels = as.character(temp$date)
  ) +
  # scale_y_continuous(breaks = c("lab_test", "exp"), labels = c("lab_test", "Exp")) +
  scale_y_discrete(breaks = c("lab_test", "exp"), labels = c("lab_test", "Exposome")) +
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

# ggsave(plot, filename = "exp_lab_test_match.pdf", width = 10, height = 7)
# ggsave(plot, filename = "exp_lab_test_match2.pdf", width = 14, height = 7)

#-------------------------------------------------------------------------------
###prepare data
dim(lab_test_expression_data)
dim(lab_test_sample_info)

lab_test_sample_info <-
lab_test_sample_info %>% 
  dplyr::filter(date.start %in% unique(diff$lab_test.date))

lab_test_sample_info

lab_test_expression_data <- 
  lab_test_expression_data %>% 
  dplyr::select(one_of(lab_test_sample_info$sample_id))

lab_test_sample_info <- 
lab_test_sample_info %>% 
  dplyr::left_join(diff[,c("lab_test.date", "lab_test.group")] %>% dplyr::distinct(lab_test.date, lab_test.group), 
                   by = c("date.start" = "lab_test.date")) %>% 
  dplyr::rename(group = lab_test.group)

exp_sample_info <-
  exp_sample_info %>% 
  dplyr::filter(as.Date(start_date) %in% diff$exp.date)

# exp_expression_data <-
#   exp_expression_data %>%
#   select(-contains("Blank")) %>% 
#   t()

exp_expression_data <-
  exp_expression_data %>%
  t() %>% 
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

###browser
exp_expression_data <- 
  exp_expression_data %>% 
  dplyr::select(one_of(exp_sample_info$sample_id))

exp_sample_info <- 
  exp_sample_info %>% 
  mutate(start_date = as.Date(start_date)) %>% 
  dplyr::left_join(diff[,c("exp.date", "exp.group")], by = c("start_date" = "exp.date")) %>% 
  dplyr::rename(group = exp.group)

dim(lab_test_sample_info)
dim(exp_sample_info)

# ##combine data
# colnames(exp_expression_data) == exp_sample_info$sample_id
# 
# data <-
#   t(exp_expression_data) %>% 
#   data.frame(., group = exp_sample_info$group, stringsAsFactors = FALSE) %>% 
#   plyr::dlply(.variables = .(group)) %>% 
#   lapply(function(x){
#     x <- 
#       x %>% dplyr::select(-group)
#     purrr::map(x, .f = function(x){mean(x)}) %>% unlist()
#   }) %>% 
#   do.call(rbind, .)
# 
# rownames(data) == exp_sample_info$group
# 
# exp_expression_data <- t(data) %>% as.data.frame()
# 
# exp_sample_info <-
#   exp_sample_info %>% 
#   plyr::dlply(.variables = .(group)) %>% 
#     lapply(function(x){
#       x[which.max(x$start_date),]
#     }) %>% 
#     do.call(rbind, .)

dim(exp_expression_data)

dim(lab_test_expression_data)

save(exp_expression_data, file = "exp_expression_data")
save(lab_test_expression_data, file = "lab_test_expression_data")

save(exp_sample_info, file = "exp_sample_info")
save(lab_test_sample_info, file = "lab_test_sample_info")


