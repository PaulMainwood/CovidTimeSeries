#Creates a plot of proportions of various selected variants vs a baseline variant
#After http://sonorouschocolate.com/covid19/index.php/UK_variant_comparison
#by @alexselby

library(lubridate)
library(tidyverse)
library(magrittr)
library(dplyr)
library(data.table)
library(stringr)

first_date <- as.Date("2022-06-1")
last_date <- as.Date("2022-09-01")

#Get case data - this takes ages 
all_data <- fread("https://cog-uk.s3.climb.ac.uk/phylogenetics/latest/cog_metadata.csv.gz")
our_data <- select(filter(all_data, sample_date >= first_date), lineage, sample_date)
#Keep all lineage names and select only the pure ones.

grp = group_by(our_data, sample_date, lineage) %>% summarise(number = n())

#Function to take new and base variant names and produce dataframes of one as 
#fraction of the other, and regress the log of the proportions by date to give 
#a best-fit line.
proportions_new_variant <- function(new_variant_name, base_variant_name) {
  base = select(filter(grp,  str_detect(base_variant_name, lineage)), sample_date, number)
  new = select(filter(grp, str_detect(new_variant_name, lineage)), sample_date, number)
  
  #Collect counts of each lineage in weighting by total of BA sequences
  df <- purrr::reduce(list(base, new), dplyr::left_join, by = 'sample_date')
  colnames(df) = c("sample_date", "base", "new")

  df <- mutate(df, new_frac = new / base)
  wts <- coalesce(1/(1/df$base + 1/df$new),0)
  
  #Regression with log2 transformation
  reg <- lm(log2(new_frac) ~ date(sample_date), data = df, weights = wts)
  
  return(list(df = df, reg = reg))
  }


#Use the function to grab whichever variants and base you want and plot. Manual
#awkwardness here, since when automate, get weird visuals.
plot(date(proportions_new_variant("BA.4.6", "BA.5.1")$df$sample_date), log2(proportions_new_variant("BA.4.6", "BA.5.1")$df$new_frac), pch = 19, col = 1, xlim = c(first_date, last_date), ylim = c(-12, 3), yaxt = "n", xaxt = "n", ylab = "Ratio to BA.5, % - log2 scale", xlab = "Date", main = "BA.2, BA.4, BA.5.1 as a ratio to BA.5 in COG-UK sequences")

yticks = seq(-12, 3)
#axis(side = 2, at = yticks, labels = c("0.024%","0.049%", "0.098%", "0.195%", "0.391%", "0.781%", "1.563%", "3.125%", "6.25%", "12.5%", "25%", "50%", "1x", "2x", "4x"))
axis(side = 2, at = yticks, labels = c("0.024%","0.049%", "0.098%", "0.195%", "0.391%", "0.781%", "1.563%", "3.125%", "6.25%", "12.5%", "25%", "50%", "1x", "2x", "4x", "8x"))
xticks = seq(first_date, last_date, 7)
axis(1, at = xticks, labels = format(date(xticks), format = "%d %b"))

points(date(proportions_new_variant("BA.5.2", "BA.5.1")$df$sample_date), pch = 19, log2(proportions_new_variant("BA.5.2", "BA.5.1")$df$new_frac), col = 2)
points(date(proportions_new_variant("BF.5", "BA.5.1")$df$sample_date), pch = 19, log2(proportions_new_variant("BF.5", "BA.5.1")$df$new_frac), col = 3)


abline(proportions_new_variant("BA.4.6", "BA.5.1")$reg, col = 2)
abline(proportions_new_variant("BA.5.2", "BA.5.1")$reg, col = 3)
abline(proportions_new_variant("BF.5", "BA.5.1")$reg, col = 4)
abline(h = 0, col = "gray", lty= "dashed")
legend("bottomright", inset = 0.05, legend = c("BA.2", "BA.4", "BA.5.1"), col = c(2, 3, 4), lty = c(1, 1))
title(sub="Plot of COG-UK data, by @paulmainwood, based on work from @alexselby1770. Regression weighted by inverse probability of selection; last few data points based on very small values", adj=1, line=4, font=1, cex.sub = 0.6)

