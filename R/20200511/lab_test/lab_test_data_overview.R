##avoid source
no_function()

sxtTools::setwd_project()
library(tidyverse)
load("data_20200511/lab_test/expression_data")
load("data_20200511/lab_test/variable_info")
load("data_20200511/lab_test/sample_info")

setwd("data_analysis/lab_test_analysis/data_overview")

temp_data <-
  log(expression_data + 1, 2) %>% 
  apply(1, function(x){
    (x - mean(x)) / sd(x)
  }) %>% 
  t()

colnames(temp_data) <- as.character(sample_info$date.start)

##heatmap of lab_test
library(circlize)
library(ComplexHeatmap)

col_fun = colorRamp2(
  breaks = seq(min(temp_data), max(temp_data), length.out = 90),
  colors =
    viridis::magma(n = 100)[-c(1:10)],
  transparency = 0
)

plot <- 
  temp_data %>% 
  ComplexHeatmap::Heatmap(cluster_columns = FALSE,
                          show_column_names = TRUE,
                          show_row_names = TRUE,
                          clustering_method_rows = "ward.D",
                          clustering_method_columns = "ward.D",
                          clustering_distance_columns = "euclidean",
                          clustering_distance_rows = "euclidean",
                          col = col_fun,
                          km = 2, border = TRUE, 
                          row_dend_reorder = TRUE, 
                          column_dend_reorder = TRUE,
                          column_names_rot = 45, 
                          name = "Z-score")

plot <- ggplotify::as.ggplot(plot)

plot

ggsave(plot, filename = "lab_test_heatmap.pdf", width = 7, height = 7)


plot <- 
temp_data %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "variable_id") %>%
  tidyr::pivot_longer(cols = -variable_id,
                      names_to = "date",
                      values_to = "value") %>%
  ggplot(aes(date, value, group = variable_id)) +
  labs(x = "", y = "Z-score") +
  geom_line(aes(), 
            show.legend = FALSE) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      size = 12,
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 13),
    panel.grid.minor = element_blank()
  )
plot
ggsave(plot, filename = "line_plot.pdf", width = 7, height = 7)




plot <-
  expression_data %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "variable_id") %>%
  tidyr::pivot_longer(cols = -variable_id,
                      names_to = "date",
                      values_to = "value") %>%
  dplyr::mutate(date = as.character(date)) %>%
  ggplot(data = ., aes(date, value, group = variable_id)) +
  labs(x = "", y = "Z-score") +
  geom_line(aes(),
            show.legend = FALSE) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 10),
    panel.grid.minor = element_blank()
  ) +
  facet_wrap(facets = vars(variable_id),scales = "free_y")

plot
ggsave(plot,
       filename = "all_line_plot.pdf",
       width = 15,
       height = 12)


blood_range <- readxl::read_xlsx("Blood range.xlsx")

match(blood_range$variable_id, variable_info$variable_id)

blood_range


plot <-
  expression_data %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "variable_id") %>%
  tidyr::pivot_longer(cols = -variable_id,
                      names_to = "date",
                      values_to = "value") %>%
  dplyr::mutate(date = as.character(date)) %>%
  ggplot(data = ., aes(date, value, group = variable_id)) +
  labs(x = "", y = "Z-score") +
  geom_line(aes(),
            show.legend = FALSE) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 10),
    # axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 10),
    panel.grid.minor = element_blank()
  ) +
geom_rect(
  aes(
    xmin = -Inf,
    xmax = Inf,
    ymin = Minimum,
    ymax = Maximum
  ),
  fill = ggsci::pal_aaas()(n = 10)[5],
  alpha = 0.5,
  data = blood_range,
  inherit.aes = FALSE
) +
  facet_wrap(facets = vars(variable_id),scales = "free_y")

plot
ggsave(plot,
       filename = "all_line_plot.pdf",
       width = 15,
       height = 12)








