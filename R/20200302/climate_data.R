sxtTools::setwd_project()
rm(list=ls())

load("data_20200302/exposome_data/phenotype_table")

setwd("data_20200302/climate_data/")
library(rnoaa)
library(ropenaq)

# cities_table <- ropenaq::aq_cities(country = "US")
# 
# cities_table$city

##get air quality data
locations_table <-
  ropenaq::aq_locations(country = "US",
               parameter = "pm25")

library(rworldmap)
library(maps)
us_map_data <- map_data(map = "state")

us_map <- ggplot() + geom_map(
  data = us_map_data,
  map = us_map_data,
  aes(map_id = region, x = long, y = lat),
  fill = "#D9D9D9",
  colour = "black"
) +
  theme_void()
us_map

cal_map_data <- map_data(map = "state", region = "cal")

cal_map <- ggplot() + geom_map(
  data = cal_map_data,
  map = cal_map_data,
  aes(map_id = region, x = long, y = lat),
  fill = "#D9D9D9",
  colour = "black"
) +
  theme_void()
cal_map

locations_table <-
  locations_table %>%
  dplyr::filter(longitude > -130 & longitude < -100) %>%
  dplyr::filter(latitude > 25 & latitude < 50)

us_map +
  geom_point(
    data = locations_table,
    aes(x = longitude, y = latitude),
    size = 1,
    col = "#8DD3C7"
  ) +
  theme(
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none",
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_blank()
  ) +
  ggtitle("OpenAQ data sources with geographical coordinates") +
  theme(plot.title = element_text(lineheight = 1, face = "bold"))

us_map +
  geom_point(
    data = phenotype_table,
    aes(x = longitude, y = latitude),
    size = 2,
    col = "#FB8072"
  ) +
  ggrepel::geom_label_repel(data = phenotype_table, 
                            mapping = aes(x = longitude, y = latitude, 
                               label = location)) +
  theme(
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none",
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_blank()
  ) +
  # ggtitle("OpenAQ data sources with geographical coordinates") +
  theme(plot.title = element_text(lineheight = 1, face = "bold"))



##only for bay area
# cal_map2 <- 
  ggplot() + geom_map(
  data = cal_map_data,
  map = cal_map_data,
  aes(map_id = region, x = long, y = lat),
  fill = "#D9D9D9",
  colour = "black"
) +
  scale_x_continuous(limits = c(-123,-121.5)) +
  scale_y_continuous(limits = c(37.5,38)) +
  theme_bw() +
    ggrepel::geom_label_repel(data = phenotype_table, 
                              mapping = aes(x = longitude, y = latitude, 
                                            label = location)) 
    


station_data <- rnoaa::ghcnd_stations()

lat_lon_df <-
  phenotype_table %>%
  dplyr::select(location,
                latitude,
                longitude) %>%
  unique() %>%
  ungroup() %>%
  dplyr::rename(id = location) %>%
  mutate(id = factor(id))

lat_lon_df <- 
  lat_lon_df %>% 
  distinct(latitude, longitude, .keep_all = TRUE)

stations <-
  meteo_nearby_stations(
    lat_lon_df = as.data.frame(lat_lon_df),
    station_data = station_data,
    radius = 15,
    year_min = "2015",
    var = "all"
  )

# stations <-
#   unique(bind_rows(stations) %>% 
#            dplyr::select(-distance))


stations2 <-
  mapply(function(x,y){
    data.frame(location = x, 
               y,
               stringsAsFactors = FALSE) %>% 
      list()
  },
  x = names(stations),
  y = stations) %>% 
  bind_rows() %>% 
  distinct(id, .keep_all = TRUE)


library("ggmap")
# map <- ggmap::get_googlemap(center  = "450 serra mall, stanford, ca 94305, usa", 
#                       zoom = 11, 
#                       maptype = "roadmap",
#                       color = "color"
#                       # style = c(feature = "administrative.country",
#                       #           element = "labels",visibility = "off")
#                       )



us_map <-
  get_stamenmap(
    bbox = c(-127, 23, -62, 53),
    zoom = 5,
    maptype = "toner-lite",
    color = "bw"
  )

ggmap(us_map) +
  geom_point(
    aes(x = longitude, y = latitude),
    data = lat_lon_df,
    col = "#FB8072",
    size = 4
  ) +
  labs(x = "Longitude", y = "Latitude") +
  ggrepel::geom_label_repel(mapping = aes(x = longitude,
                                          y = latitude,
                                          label = id),
                            data = lat_lon_df, color = "#FB8072"
  ) +
  # geom_point(
  #   aes(x = longitude, y = latitude),
  #   data = stations2,
  #   col = "green",
  #   size = 1
  # ) +
  theme_bw()


cal_map <- get_stamenmap(
  bbox = c(-122.7, 37.3, -120.5, 38.6),
  zoom = 10,
  maptype = "toner-lite",
  color = "bw"
)


ggmap(cal_map) +
  geom_point(
    aes(x = longitude, y = latitude, color = id),
    data = lat_lon_df,
    size = 4, show.legend = FALSE
  ) +
  labs(x = "Longitude", y = "Latitude") +
  ggrepel::geom_text_repel(mapping = aes(x = longitude,
                                          y = latitude,
                                          label = id,
                                          color = id),
                            data = lat_lon_df, show.legend = FALSE
  ) +
  geom_point(
    aes(x = longitude, y = latitude, color = location),
    data = stations2,
    size = 2, show.legend = FALSE, shape = 17
  ) +
  theme_bw()

  
#https://www1.ncdc.noaa.gov/pub/data/ghcn/daily/readme.txt

monitors <- stations2$id

all_monitors_clean <-
  rnoaa::meteo_pull_monitors(monitors = monitors,
                      date_min = "2015-12-31",
                      date_max = "2016-05-01", 
                      var = "all") %>%
  dplyr::rename(day = date, location = id)

all_monitors_clean %>% head() %>% knitr::kable()

head(all_monitors_clean)

climate_data <- all_monitors_clean %>% 
  left_join(stations2, by = c("location" = "id")) %>% 
  dplyr::select(station = location, location = location.y, name, latitude, longitude, day,
                everything())

test <-
climate_data %>% 
  dplyr::filter(location == "Campus-home, weekend (No001_beads)") %>% 
  plyr::dlply(.variables = .(station))


save(climate_data, file = "climate_data")




location <- sort(phenotype_table$location)
location[grep("Campus", location)] <- "Campus"
location[grep("Boston", location)] <- "Boston"
location[grep("Porter", location)] <- "Campus"

location[grep("SF and Novato, CA", location)] <- "Novato, CA"

location[grep("UCSF, easter, UC Davis", location)] <- "UC Davis, CA"

phenotype_table$location <- location


phenotype_table %>% 
  ggplot(aes(x = location)) +
  geom_bar(aes(fill = location), show.legend = FALSE) +
  labs(x = "Location", y = "Sample number") +
  scale_fill_aaas() +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 13))

phenotype_table %>% 
  mutate(start_date = as.Date(start_date)) %>% 
  ggplot(aes(x = location)) +
  geom_point(aes(x = start_date, 
                 y = location, color = location),
             show.legend = FALSE,
             size = 4) +
  labs(x = "Date", y = "Location") +
  scale_color_aaas() +
  scale_x_continuous(trans = "date",
                     breaks = c(as.Date(phenotype_table$start_date)),
                     labels = as.character(phenotype_table$start_date)
                     ) +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text.x = element_text(size = 10, 
                                   angle = 45,
                                   vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 13),
        )

climate_data2 <- 
climate_data %>% 
  dplyr::filter(location == "Campus-home, weekend (No001_beads)") %>% 
  dplyr::filter(station == "USR0000CALT") %>% 
  mutate(day = as.Date(day)) %>% 
  mutate(tavg = tavg/10)

climate_data2 %>% 
ggplot(aes(x = day, 
           y = tavg),) +
  geom_point(aes(x = day, 
                 y = tavg),
             show.legend = FALSE,
             size = 3, color = "#B3DE69") +
  labs(x = "Date", y = "Temperature (C, avg)") +
  geom_smooth(color = '#B3DE69') +
  scale_color_aaas() +
  # scale_x_continuous(trans = "date",
  #                    breaks = c(as.Date(climate_data2$day)),
  #                    labels = as.character(climate_data2$day)
  # ) +
  theme_bw() +
  labs(title = "Campus") +
  theme(axis.title = element_text(size = 15),
        axis.text.x = element_text(size = 13, 
                                   # angle = 45,
                                   # vjust = 1, hjust = 1
                                   ),
        axis.text.y = element_text(size = 13),
  )



