# script to calculate pairwise distances between two coordinates

# load libraries
library(dplyr)
library(ggplot2)
library(maps) #to plot on map
library(geosphere) #to compute distance between coords

# load data
data = readRDS("DATA/Clough/clough_data.rds")
data 

# select study of interest
clough_data = data %>% 
  filter(Study_id == "11_Clough") %>%
  mutate(Network_id = sub("_.*", "", Network_id))

# explore unique coordinates
coords = clough_data %>% 
  distinct(Network_id, Latitude, Longitude)

# map
world <- map_data("world")

# plot
ggplot() +
  geom_polygon(
    data = world,
    aes(x = long, y = lat, group = group),
    fill = "grey90",
    colour = "white"
  ) +
  geom_point(
    data = coords,
    aes(x = Longitude, y = Latitude),
    size = 1
  ) +
  coord_quickmap(
    xlim = range(coords$Longitude, na.rm = TRUE) + c(-2.5, 2.5),
    ylim = range(coords$Latitude, na.rm = TRUE) + c(-2.5, 2.55)
  ) +
  theme_classic()


# all pairwise distances
distances = coords %>%
  rename(
    Network_1 = Network_id,
    Latitude_1 = Latitude,
    Longitude_1 = Longitude
  ) %>%
  cross_join( #compare with every possible combination
    coords %>%
      rename(
        Network_2 = Network_id,
        Latitude_2 = Latitude,
        Longitude_2 = Longitude
      )
  ) %>%
  filter(Network_1 < Network_2) %>% #to avoid identical rows
  rowwise() %>%
  mutate(
    Distance_km = distVincentyEllipsoid( #
      c(Longitude_1, Latitude_1),
      c(Longitude_2, Latitude_2)
    ) / 1000
  ) %>%
  ungroup() %>%
  arrange(Distance_km)

distances

saveRDS(clough_data, "clough_data.rds")
saveRDS(distances, "clough_pairwise_distance.rds")
