######################
#### Transect map ####
######################

library(readr)
library(tidyr)
library(dplyr)
library(broom)
library(ggplot2)
library(ggthemes)
library(mapdata)
library(maps)
library(rgbif)
library(CoordinateCleaner)
library(ggrepel)
library(png)
library(gridExtra)
library(colourpicker)

sites <- read.csv("~/Dropbox/master_data/site/site_details.csv")

(sites.map <- ggplot(sites, aes(x = Mean.Long, y = Mean.Lat)) +
    borders("worldHires", ylim = c(55, 58), xlim=c(-4.5,-3.2), colour = "gray40", fill = "gray40", size = 0.3) +
    # get a high resolution map of the world
    #theme_map() +
    coord_fixed(1.4)+
    theme_bw()+
    geom_point(alpha = 0.7, size = 1, colour = "aquamarine3"))  # alpha controls the transparency, 1 is fully transparent


UK <- map_data("world") %>%
  filter(region == "UK")
(sites.map <- ggplot(UK, aes(x = long, y = lat)) +
    borders("worldHires",  colour = "gray40", fill = "gray40", size = 0.3) +
    # get a high resolution map of the world
    #theme_map() +
    coord_fixed(1.4)+
    geom_point(data=sites, aes(x = Mean.Long, y = Mean.Lat), alpha = 0.7, size = 1, colour = "aquamarine3"))

ylim = c(55, 58), xlim=c(-4.5,-3.2),
UK <- map_data("world") %>%
  filter(region == "UK")
ggplot(UK, aes(x = long, y = lat)) +
  geom_polygon() +
  coord_fixed(1.4) 
