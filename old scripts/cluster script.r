library(tidyverse)
library(TSclust)
library(lubridate)
library(cluster)
library(factoextra)
library(runner)

#Get case data for 
raw_data <- read.csv("https://api.coronavirus.data.gov.uk/v2/data?areaType=ltla&metric=newCasesBySpecimenDate&format=csv")
data <- arrange((raw_data %>% filter(date >= "2020-05-01")), date)

b = pivot_wider(data, id_cols = date, names_from = areaName, values_from = newCasesBySpecimenDate)
b <- b %>% select_if(~!any(is.na(.))) %>% select(-date) 


distances <- diss(b, "COR")

km.res <- kmeans(distances, 7, nstart = 10)
fviz_cluster(km.res, distances)
