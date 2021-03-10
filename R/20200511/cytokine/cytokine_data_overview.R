##avoid source
no_function()

sxtTools::setwd_project()
library(tidyverse)
load("data_20200511/cytokine/expression_data")
load("data_20200511/cytokine/variable_info")
load("data_20200511/cytokine/sample_info")

setwd("data_analysis/cytokine_analysis/data_overview")

temp_data <-
  log(expression_data + 1, 2) %>% 
  apply(1, function(x){
    (x - mean(x)) / sd(x)
  }) %>% 
  t()

colnames(temp_data) <- as.character(sample_info$CollectionDate)

##heatmap of cytokine
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
                          show_row_names = FALSE,
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

ggsave(plot, filename = "cytokine_heatmap.pdf", width = 7, height = 7)


plot <- 
temp_data %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "variable_id") %>%
  tidyr::pivot_longer(cols = -variable_id,
                      names_to = "date",
                      values_to = "value") %>%
  dplyr::mutate(class = case_when(
    date %in% c(
      "2016-01-15",
      "2016-01-19",
      "2016-01-26",
      "2016-02-24",
      "2016-03-03"
    ) ~ "YES",
    TRUE ~ "NO"
  )) %>%
  # dplyr::filter(variable_id %in% "IL17F") %>%
  ggplot(aes(date, value, group = variable_id)) +
  labs(x = "", y = "Z-score") +
  geom_line(aes(), 
            show.legend = FALSE) +
  geom_rect(
    aes(
      xmin = 4.5,
      xmax = 9.5,
      ymin = -Inf,
      ymax = Inf
    ),
    fill = ggsci::pal_aaas()(n = 10)[5],
    alpha = 0.5,
    data = data.frame(),
    inherit.aes = FALSE
  ) +
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
  temp_data %>%
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
  # geom_rect(
  #   aes(
  #     xmin = 4.5,
  #     xmax = 9.5,
  #     ymin = -Inf,
  #     ymax = Inf
  #   ),
  #   fill = ggsci::pal_aaas()(n = 10)[5],
  #   alpha = 0.5,
  #   data = data.frame(),
  #   inherit.aes = FALSE
  # ) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 10),
    panel.grid.minor = element_blank()
  ) +
  annotate(geom = "rect",
           xmin = 4.5,
           xmax = 9.5,
           ymin = -Inf,
           ymax = Inf,
           fill = ggsci::pal_aaas()(n = 10)[5],
           alpha = 0.5) +
  facet_wrap(facets = vars(variable_id))

plot
ggsave(plot,
       filename = "all_line_plot.pdf",
       width = 15,
       height = 12)
