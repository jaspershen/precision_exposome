sxtTools::setwd_project()

setwd("data_analysis/")
load("../data/metabolite_table")
load("../data/metabolite_tags")
load("../data/phenotype_table")



dim(metabolite_table)


plot(density(metabolite_table$B001_1))
plot(density(metabolite_table$B001_2))
plot(density(metabolite_table$B001_3))
plot(density(metabolite_table$B002_1))


##PCA analysis
data <- 
  metabolite_table %>% 
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


x %>%
  as.data.frame() %>%
  rownames_to_column(var = "sample_ID") %>%
  left_join(phenotype_table, by = "sample_ID") %>%
  mutate(date = as.character(`filter collect`)) %>%
  ggplot(aes(PC1, PC2, colour = date)) +
  geom_point() +
  ggrepel::geom_label_repel(aes(PC1, PC2, label = date)) +
  theme_bw()



x %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "sample_ID") %>% 
  left_join(phenotype_table, by = "sample_ID") %>% 
  mutate(date = as.character(`filter collect`)) %>% 
  mutate(location = stringr::str_replace(comments, "\\(.{1,20}\\)", "") %>% stringr::str_trim()) %>% 
  mutate(location = stringr::str_replace(location, "weekdays", "weekday")) %>% 
  ggplot(aes(PC1, PC2, colour = comments)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_point() +
  ggrepel::geom_label_repel(aes(PC1, PC2, label = date)) +
  theme_bw() +
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 15))




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
