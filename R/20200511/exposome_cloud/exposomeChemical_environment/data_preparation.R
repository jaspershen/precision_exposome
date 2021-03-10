##avoid source
no_function()

##load data
sxtTools::setwd_project()
library(tidyverse)
rm(list = ls())

####load environment data
sxtTools::setwd_project()
load("data_20200511/environment/expression_data")
environment_expression_data <- expression_data
load("data_20200511/environment/sample_info")
environment_sample_info <- sample_info
load("data_20200511/environment/variable_info")
environment_variable_info <- variable_info

###exposome
load("data_20200511/exposome/expression_data")
exposomeChemical_expression_data <- expression_data
load("data_20200511/exposome/sample_info")
exposomeChemical_sample_info <- sample_info
load("data_20200511/exposome/variable_info")
exposomeChemical_variable_info <- variable_info

setwd("data_analysis/exposomeChemical_environment/")

dim(exposomeChemical_variable_info)

dim(environment_variable_info)

temp1 <- 
exposomeChemical_sample_info$start_date %>% 
  as.Date() %>% 
  data.frame(date = ., 
             exposomeChemical = 1,
             stringsAsFactors = FALSE)

temp2 <- 
environment_sample_info$date.start %>% 
  as.Date() %>% 
  data.frame(date = ., 
             environment = 1,
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
lapply(environment_sample_info$date.start %>% as.Date(), function(x){
  x <- x - as.Date(exposomeChemical_sample_info$start_date) 
  as.Date(exposomeChemical_sample_info$start_date) [which(x >=0 & x <= 0)]
}) 

names(diff) <- environment_sample_info$date.start %>% as.character()

diff <- 
purrr::map2(diff, names(diff), .f = function(x,y){
  if(length(x) == 0){
    return(NULL)
  }
  z <- data.frame(y, x, stringsAsFactors = FALSE)
  colnames(z) <- c("environment", "exposomeChemical")
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
  mutate(exposomeChemical = as.character(exposomeChemical), 
         environment = as.character(environment)) %>% 
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
  

diff$exposomeChemical.date <- as.Date(diff$exposomeChemical.date)

diff$environment.date <- as.Date(diff$environment.date)

diff$exposomeChemical.group <- as.character(diff$exposomeChemical.group)

plot <- 
temp %>% 
  mutate(date = as.Date(date)) %>% 
  dplyr::mutate(class = factor(class, levels = c("environment", "exposomeChemical"))) %>% 
  ggplot(aes(date, class)) +
  geom_point(aes(fill = group),
             color = "black",
             alpha = 0.6, size = 4,
             shape = 21, show.legend = FALSE) +
  geom_segment(aes(x = exposomeChemical.date, 
                   xend = environment.date, 
                   y = exposomeChemical.class, 
                   yend = environment.class,
                   color = exposomeChemical.group), data = diff, show.legend = FALSE) +
  scale_fill_manual(values = c("No" = "grey", 
                                "1" = ggsci::pal_aaas()(10)[1],
                                "2" = ggsci::pal_aaas()(10)[2],
                                "3" = ggsci::pal_aaas()(10)[3],
                                "4" = ggsci::pal_aaas()(10)[4],
                                "5" = ggsci::pal_aaas()(10)[5],
                               "6" = ggsci::pal_aaas()(10)[6],
                               "7" = ggsci::pal_aaas()(10)[7],
                               "8" = ggsci::pal_aaas()(10)[8],
                               "9" = ggsci::pal_aaas()(10)[9],
                               "10" = ggsci::pal_aaas()(10)[10]
                               )) +
  scale_x_continuous(trans = "date",
                     breaks = c(as.Date(temp$date)),
                     labels = as.character(temp$date)
  ) +
  # scale_y_continuous(breaks = c("environment", "exposomeChemical"), labels = c("environment", "Exp")) +
  scale_y_discrete(breaks = c("environment", "exposomeChemical"), 
                   labels = c("Environment", "Exposome chemical")) +
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

save(temp, file = "temp")
save(diff, file = "diff")
# ggsave(plot, filename = "exposomeChemical_environment_match.pdf", width = 10, height = 7)
# ggsave(plot, filename = "exposomeChemical_environment_match2.pdf", width = 14, height = 7)


#-------------------------------------------------------------------------------
###prepare data
dim(environment_expression_data)
dim(environment_sample_info)

environment_sample_info <-
  environment_sample_info %>% 
  dplyr::filter(date.start %in% unique(diff$environment.date))

environment_sample_info

environment_expression_data <- 
  environment_expression_data %>% 
  dplyr::select(one_of(environment_sample_info$sample_id))

environment_sample_info <- 
  environment_sample_info %>% 
  dplyr::left_join(diff[,c("environment.date", "environment.group")] %>% 
                     dplyr::distinct(environment.date, environment.group), 
                   by = c("date.start" = "environment.date")) %>% 
  dplyr::rename(group = environment.group)

exposomeChemical_sample_info <-
  exposomeChemical_sample_info %>% 
  dplyr::filter(as.Date(start_date) %in% diff$exposomeChemical.date)

exposomeChemical_expression_data <-
  exposomeChemical_expression_data %>%
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
exposomeChemical_expression_data <- 
  exposomeChemical_expression_data %>% 
  dplyr::select(one_of(exposomeChemical_sample_info$sample_id))

exposomeChemical_sample_info <- 
  exposomeChemical_sample_info %>% 
  mutate(start_date = as.Date(start_date)) %>% 
  dplyr::left_join(diff[,c("exposomeChemical.date", "exposomeChemical.group")], by = c("start_date" = "exposomeChemical.date")) %>% 
  dplyr::rename(group = exposomeChemical.group)

dim(environment_sample_info)
dim(exposomeChemical_sample_info)

##combine data
colnames(exposomeChemical_expression_data) == exposomeChemical_sample_info$sample_id

data <-
  t(exposomeChemical_expression_data) %>%
  data.frame(., group = exposomeChemical_sample_info$group, stringsAsFactors = FALSE) %>%
  plyr::dlply(.variables = .(group)) %>%
  lapply(function(x){
    x <-
      x %>% dplyr::select(-group)
    purrr::map(x, .f = function(x){mean(x)}) %>% unlist()
  }) %>%
  do.call(rbind, .)

rownames(data) == exposomeChemical_sample_info$group

exposomeChemical_expression_data <- t(data) %>% as.data.frame()

exposomeChemical_sample_info <-
  exposomeChemical_sample_info %>%
  plyr::dlply(.variables = .(group)) %>%
  lapply(function(x){
    x[which.max(x$start_date),]
  }) %>%
  do.call(rbind, .)

dim(exposomeChemical_expression_data)

dim(environment_expression_data)

remove_idx <- 
apply(environment_expression_data, 1, function(x){
  length(unique(x)) == 1
}) %>% 
  which() %>% 
  unname()

if(length(remove_idx) > 0){
  environment_expression_data <- environment_expression_data[-remove_idx,]
  environment_variable_info <- environment_variable_info[-remove_idx, ,drop = FALSE]
}


which(is.na(environment_expression_data), arr.ind = TRUE)

environment_expression_data <-
  environment_expression_data %>% 
  apply(1, function(x){
    x[is.na(x)] <- min(x, na.rm = TRUE)
    x
  }) %>% 
  t() %>% 
  as.data.frame()

save(exposomeChemical_expression_data, file = "exposomeChemical_expression_data")
save(environment_expression_data, file = "environment_expression_data")

save(exposomeChemical_sample_info, file = "exposomeChemical_sample_info")
save(environment_sample_info, file = "environment_sample_info")
save(exposomeChemical_variable_info, file = "exposomeChemical_variable_info")
