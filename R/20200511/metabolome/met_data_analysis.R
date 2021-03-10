####load metabolomics data
sxtTools::setwd_project()
load("data_20200511/metabolome/expression_data")
met_expression_data <- expression_data
load("data_20200511/metabolome/sample_info")
met_sample_info <- sample_info
load("data_20200511/metabolome/variable_info")
met_variable_info <- variable_info


##load exposome data
load("data_20200511/exposome/expression_data")
load("data_20200511/exposome/sample_info")
load("data_20200511/exposome/variable_info")

exp_expression_data <- expression_data

exp_sample_info <- sample_info

exp_variable_info <- variable_info


setwd("data_analysis/met_exp/")

