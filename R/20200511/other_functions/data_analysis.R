sxtTools::setwd_project()

setwd("data_analysis/")
##load exposome data
load("../data_20200511/exposome/expression_data")
load("../data_20200511/exposome/sample_info")
load("../data_20200511/exposome/variable_info")

exp_expression_data <- expression_data

exp_sample_info <- sample_info

exp_variable_info <- variable_info

dim(expression_data)


##PCA analysis
data <- 
  expression_data %>% 
  select(-contains("Blank"))


data <- log(data + 1, 10)

data <- 
  apply(data, 1, function(x){
    (x - mean(x))/sd(x)
  })

data <- 
  data %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "Sample_ID") %>% 
  mutate(Sample_ID = stringr::str_replace(Sample_ID, "_[0-9]{1,2}", "")) %>% 
  plyr::dlply(.variables = "Sample_ID") %>% 
  lapply(function(x){
    apply(x[,-1], 2, mean)
  }) %>% 
  do.call(rbind, .)
  
pca_object <- 
  prcomp(x = data)

x <- pca_object$x

x[,1:4]

sample_info$location[grep("Campus", sample_info$location)] <- "Campus"

as.Date_origin <- function(x){
  as.Date(x, origin = '1970-01-01')
}

plot <- 
x %>%
  as.data.frame() %>%
  rownames_to_column(var = "sample_id") %>%
  left_join(sample_info, by = "sample_id") %>%
  mutate(date = as.integer(as.Date(start_date))) %>%
  # mutate(date = as.character(start_date)) %>%
  ggplot(aes(PC1, PC2, colour = date)) +
  geom_vline(xintercept = 0, linetype = 2, color = "black") +
  geom_hline(yintercept = 0, linetype = 2, color = "black") +
  geom_point(shape = 16) +
  guides(colour = guide_colourbar(title = "Date")) +
  # scale_colour_continuous() +
  scale_colour_gradientn(colours = c(
    alpha("#8DD3C7", 1),
    alpha("#8DD3C7", 0.4),
    alpha("#FB8072", 0.4),
    alpha("#FB8072", 1)
  ),
  labels=as.Date_origin) +
  ggrepel::geom_label_repel(aes(PC1, PC2, label = location)) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 13)) +
  labs(x = paste("PC1 (", round(summary(pca_object)$importance[2,1], 4)*100, "%)", sep = ""),
       y = paste("PC2 (", round(summary(pca_object)$importance[2,2], 4)*100, "%)", sep = "")) 

plot

ggsave(plot, file = "pca1.pdf", width = 8, height = 7)

ggsave(plot, file = "pca1.png", width = 8, height = 7)

plot <- 
x %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "sample_id") %>% 
  left_join(sample_info, by = "sample_id") %>% 
  mutate(date = as.character(start_date)) %>% 
  mutate(location = stringr::str_replace(location, "\\(.{1,20}\\)", "") %>% stringr::str_trim()) %>% 
  mutate(location = stringr::str_replace(location, "weekdays", "weekday")) %>% 
  ggplot(aes(PC1, PC2, colour = location)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_point() +
  ggsci::scale_color_aaas() +
  ggrepel::geom_label_repel(aes(PC1, PC2, label = date)) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 13)) +
  labs(x = paste("PC1 (", round(summary(pca_object)$importance[2,1], 4)*100, "%)", sep = ""),
       y = paste("PC2 (", round(summary(pca_object)$importance[2,2], 4)*100, "%)", sep = "")) 

plot

ggsave(plot, file = "pca2.pdf", width = 8, height = 7)
ggsave(plot, file = "pca2.png", width = 8, height = 7)


####load metabolomics data
load("../data_20200511/metabolome/expression_data")
met_expression_data <- expression_data
load("../data_20200511/metabolome/sample_info")
met_sample_info <- sample_info
load("../data_20200511/metabolome/variable_info")
met_variable_info <- variable_info


###PCA analysis of metabolomics
rsd <- 
  apply(met_expression_data, 1, function(x){
    sd(x)/mean(x)
  })

idx <- 
  which(rsd > 0.3)

temp_data <- met_expression_data[idx,]

temp_data <- 
  log(temp_data + 1, 10)

temp_data <- 
  apply(temp_data, 1, function(x){
    (x - mean(x)) / sd(x)
  })


pca_object <- 
  prcomp(x = temp_data)

x <- pca_object$x

x <- 
  x[,1:2] %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "sample_id") %>% 
  dplyr::left_join(met_sample_info, by = "sample_id")

x$CollectionDate <- 
  as.integer(as.Date(x$CollectionDate))


plot <- 
  x %>%
  ggplot(aes(PC1, PC2, colour = CollectionDate)) +
  geom_vline(xintercept = 0, linetype = 2, color = "black") +
  geom_hline(yintercept = 0, linetype = 2, color = "black") +
  geom_point(shape = 16, size = 3) +
  guides(colour = guide_colourbar(title = "Date")) +
  # scale_colour_continuous() +
  scale_colour_gradientn(colours = c(
    alpha("#8DD3C7", 1),
    alpha("#8DD3C7", 0.4),
    alpha("#FB8072", 0.4),
    alpha("#FB8072", 1)
  ),
  labels=as.Date_origin) +
  ggrepel::geom_label_repel(aes(PC1, PC2, 
                                label = as.character(as.Date_origin(CollectionDate)))) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 13)) +
  labs(x = paste("PC1 (", round(summary(pca_object)$importance[2,1], 4)*100, "%)", sep = ""),
       y = paste("PC2 (", round(summary(pca_object)$importance[2,2], 4)*100, "%)", sep = "")) 

plot

ggsave(plot, file = "met_pca.pdf", width = 8, height = 7)

ggsave(plot, file = "met_pca.png", width = 8, height = 7)









china <- 
  ne_countries(scale = "medium", returnclass = "sf", country = "united states")


ggplot(data = china) +
  geom_sf() +
  xlab("Longitude") + ylab("Latitude") +
  ggtitle("World map", subtitle = paste0("(", length(unique(world$NAME)), " countries)"))



ggplot(data = world) +
  geom_sf(aes(fill = pop_est)) +
  scale_fill_viridis_c(option = "plasma", trans = "sqrt")



ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-102.15, -74.12), ylim = c(7.65, 33.97), expand = FALSE)



china <- ne_countries(scale = "medium", returnclass = "sf", country = "china")

get_googlemap("waco texas", zoom = 12) %>% ggmap()


us <- c(left = -115, bottom = 35, right = -112, top = 40)
get_stamenmap(us, zoom = 5, maptype = "toner-lite") %>% ggmap() 


library("forcats")

# define helper
`%notin%` <- function(lhs, rhs) !(lhs %in% rhs)

# reduce crime to violent crimes in downtown houston
violent_crimes <- 
  crime %>% 
  dplyr::filter(
    offense %notin% c("auto theft", "theft", "burglary"),
    -95.39681 <= lon & lon <= -95.34188,
    29.73631 <= lat & lat <=  29.78400
  ) %>% 
  mutate(
    offense = fct_drop(offense),
    offense = fct_relevel(offense, c("robbery", "aggravated assault", "rape", "murder"))
  )


qmplot(lon, lat, data = violent_crimes, maptype = "toner-lite", color = I("red"))
