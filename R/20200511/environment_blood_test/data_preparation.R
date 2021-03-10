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

###environment
load("data_20200511/environment/expression_data")
environment_expression_data <- expression_data
load("data_20200511/environment/sample_info")
environment_sample_info <- sample_info
load("data_20200511/environment/variable_info")
environment_variable_info <- variable_info

setwd("data_analysis/environment_blood_test/")

dim(environment_variable_info)

dim(blood_test_variable_info)

temp1 <- 
  environment_sample_info$STARTING_DATE %>% 
  as.Date() %>% 
  data.frame(date = ., 
             environment = 1,
             stringsAsFactors = FALSE)

temp2 <- 
  blood_test_sample_info$STARTING_DATE %>% 
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
  lapply(blood_test_sample_info$STARTING_DATE %>% as.Date(), function(x){
    x <- x - as.Date(environment_sample_info$STARTING_DATE) 
    as.Date(environment_sample_info$STARTING_DATE) [which(x >=0 & x <= 1)]
  }) 

names(diff) <- blood_test_sample_info$STARTING_DATE %>% as.character()

diff <- 
  purrr::map2(diff, names(diff), .f = function(x,y){
    if(length(x) == 0){
      return(NULL)
    }
    z <- data.frame(y, x, stringsAsFactors = FALSE)
    colnames(z) <- c("blood_test", "environment")
    z
  }) 

diff <- diff[lapply(diff, is.null) %>% unlist() %>% `!`]

# diff[[2]] <- diff[[2]][-c(1),,drop = FALSE]

diff <- purrr::map2(.x = diff, .y = 1:length(diff), .f = function(x,y){
  data.frame(x, group = y, stringsAsFactors = FALSE)
}) %>% 
  do.call(rbind, .)

rownames(diff) <- NULL

diff <- 
  diff %>% 
  mutate(environment = as.character(environment), 
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


diff$environment.date <- as.Date(diff$environment.date)

diff$blood_test.date <- as.Date(diff$blood_test.date)

diff$environment.group <- as.character(diff$environment.group)

plot <- 
  temp %>% 
  mutate(date = as.Date(date)) %>% 
  dplyr::mutate(class = factor(class, levels = c("blood_test", "environment"))) %>% 
  ggplot(aes(date, class)) +
  geom_point(aes(fill = group),
             color = "black",
             alpha = 0.6, size = 4,
             shape = 21, show.legend = FALSE) +
  geom_segment(aes(x = environment.date, 
                   xend = blood_test.date, 
                   y = environment.class, 
                   yend = blood_test.class,
                   color = environment.group), data = diff, show.legend = FALSE) +
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
  # scale_y_continuous(breaks = c("blood_test ", "environment"), labels = c("blood_test ", "Exp")) +
  scale_y_discrete(breaks = c("blood_test ", "environment"), labels = c("Blood test", "Exposome biological")) +
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

ggsave(plot, filename = "environment_blood_test_match.pdf", width = 10, height = 7)
ggsave(plot, filename = "environment_blood_test_match2.pdf", width = 14, height = 7)


#-------------------------------------------------------------------------------
###prepare data
dim(blood_test_expression_data)
dim(blood_test_sample_info)

blood_test_sample_info <-
  blood_test_sample_info %>% 
  dplyr::filter(as.Date(STARTING_DATE) %in% unique(diff$blood_test.date))

blood_test_expression_data <- 
  blood_test_expression_data %>% 
  dplyr::select(one_of(blood_test_sample_info$sample_id))

blood_test_sample_info <- 
  blood_test_sample_info %>% 
  dplyr::mutate(STARTING_DATE = as.Date(STARTING_DATE)) %>% 
  dplyr::left_join(diff[,c("blood_test.date", "blood_test.group")] %>% 
                     dplyr::distinct(blood_test.date, blood_test.group), 
                   by = c("STARTING_DATE" = "blood_test.date")) %>% 
  dplyr::rename(group = blood_test.group)

environment_sample_info <-
  environment_sample_info %>% 
  dplyr::filter(as.Date(STARTING_DATE) %in% diff$environment.date)

environment_expression_data <-
  environment_expression_data %>%
  select(-contains("Blank")) %>% 
  t()

environment_expression_data <-
  environment_expression_data %>%
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

environment_expression_data <- 
  environment_expression_data %>% 
  dplyr::select(one_of(environment_sample_info$sample_id))

environment_sample_info <- 
  environment_sample_info %>% 
  mutate(STARTING_DATE = as.Date(STARTING_DATE)) %>% 
  dplyr::left_join(diff[,c("environment.date", "environment.group")], by = c("STARTING_DATE" = "environment.date")) %>% 
  dplyr::rename(group = environment.group)

dim(blood_test_sample_info)
dim(environment_sample_info)

##combine data
colnames(environment_expression_data) == environment_sample_info$sample_id

data <-
  t(environment_expression_data) %>% 
  data.frame(., group = environment_sample_info$group, stringsAsFactors = FALSE) %>% 
  plyr::dlply(.variables = .(group)) %>% 
  lapply(function(x){
    x <- 
      x %>% dplyr::select(-group)
    purrr::map(x, .f = function(x){mean(x)}) %>% unlist()
  }) %>% 
  do.call(rbind, .)

rownames(data) == environment_sample_info$group

environment_expression_data <- t(data) %>% as.data.frame()

environment_sample_info <-
  environment_sample_info %>% 
  plyr::dlply(.variables = .(group)) %>% 
  lapply(function(x){
    x[which.max(x$STARTING_DATE),]
  }) %>% 
  do.call(rbind, .)

dim(environment_expression_data)

dim(blood_test_expression_data)

save(environment_expression_data, file = "environment_expression_data")
save(blood_test_expression_data, file = "blood_test_expression_data")

save(environment_sample_info, file = "environment_sample_info")
save(blood_test_sample_info, file = "blood_test_sample_info")
save(environment_variable_info, file = "environment_variable_info")

