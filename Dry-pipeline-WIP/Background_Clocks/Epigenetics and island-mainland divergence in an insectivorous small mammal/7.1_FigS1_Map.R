##R Script for supplemental Figure S1 (map)
#Main author: Marie-Laurence Cossette

#set working directory
setwd("~/Desktop/SCHOOL/THESIS/Chapter 2 - Methyl/PAPER/Gitlab uploads/")

#load packages
library(dplyr)
library(ggplot2)
library(sf)
library(ggmap)

#make map of Canada that focuses on eastern part (want Ontario and Nova Scotia)
#coordinates we want
Canada <- c(
  left = -83,
  bottom = 41,
  right = -59.5,
  top = 47
)
#get map based on above coordinates
#choose 'terrain' as map type, and have it in color
Canada_map <- get_stamenmap(
  bbox = Canada, zoom=8, maptype = 'terrain',source = 'google',color = 'color')
map <- ggmap(Canada_map)
map

#same as above but adjust coordinates to zoom in on Bon Portage Island, NS
BPI <- c(
  left = -65.9,
  bottom = 43.35,
  right = -65.5,
  top = 43.65
)
BPI_map <- get_stamenmap(
  bbox = BPI, zoom=13, maptype = 'terrain',source = 'google',color = 'color')
BPI_map <- ggmap(BPI_map)
BPI_map


#same as above but adjust coordinates to zoom in on Long Island, NS
Long_I <- c(
  left = -66.5,
  bottom = 44.15,
  right = -65.40,
  top = 44.75
)
Long_I_map <- get_stamenmap(
  bbox = Long_I, zoom=13, maptype = 'terrain',source = 'google',color = 'color')
Long_I <- ggmap(Long_I_map)
Long_I

#saved each map on my computer and created the final Figure S1 on inkscape manually