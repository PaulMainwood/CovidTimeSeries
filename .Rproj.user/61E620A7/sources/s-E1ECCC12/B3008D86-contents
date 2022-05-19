library(lubridate)
library(tidyverse)
library(magrittr)
library(dplyr)

#Get case data for the three we care about
#Note because of the lack of permanent links you have to go to https://sars2.cvr.gla.ac.uk/cog-uk/, select each and paste your session ID into the link manually.

ba4 <- read.csv("https://sars2.cvr.gla.ac.uk/cog-uk/session/531aeb6ca3d9166bbe077c4f7fb795fa/download/downloadConcern?w=")
ba5 <- read.csv("https://sars2.cvr.gla.ac.uk/cog-uk/session/531aeb6ca3d9166bbe077c4f7fb795fa/download/downloadConcern?w=")
#This one takes ages
ba2 <- read.csv("https://sars2.cvr.gla.ac.uk/cog-uk/session/531aeb6ca3d9166bbe077c4f7fb795fa/download/downloadConcern?w=")

first_date <- as.Date("2022-04-10")
last_date <- as.Date("2022-06-05")
occurrences <- function(df, start_date = "2022-04-10"){
  a <- df %>% count(sample_date) %>% filter(sample_date >= as.Date(start_date))
  return(a)
}

a = occurrences(ba2)
b = occurrences(ba4)
c = occurrences(ba5)

df <- purrr::reduce(list(a,b,c), dplyr::left_join, by = 'sample_date')
colnames(df) = c("sample_date", "ba2", "ba4", "ba5")
df <- mutate(df, ba4_frac = ba4 / (ba2 + ba4 + ba5))
df <- mutate(df, ba5_frac = ba5 / (ba2 + ba4 + ba5))

regba4 <- lm(log2(df$ba4_frac) ~ date(df$sample_date))
#regba5 <- lm(log2(df$ba5_frac) ~ date(df$sample_date))

plot(date(df$sample_date), log2(df$ba4_frac), col = 2, ylim = c(-13, 0), xlim = c(first_date, last_date), yaxt = "n", xaxt = "n", ylab = "Share of total, % - log2 scale", xlab = "Date", main = "BA.4 as a proportion of COG-UK sequences")
#points(date(df$sample_date), log2(df$ba5_frac), col = 3)
yticks = seq(-12, 0)
axis(side = 2, at = yticks, labels = c("0.024%","0.049%", "0.098%", "0.195%", "0.391%", "0.781%", "1.563%", "3.125%", "6.25%", "12.5%", "25%", "50%", "100%"))
xticks = seq(first_date, last_date, 7)
axis(1, at = xticks, labels = format(date(xticks), format = "%d %b"))

abline(regba4, col = 2)
#abline(regba5, col = 3)
abline(h = 0, col = "gray", lty= "dotted")
abline(h = -1, col = "gray", lty= "dashed")
title(sub="Plot of COG-UK data, by @paulmainwood, inspired by previous work from @alexselby1770", adj=1, line=4, font=1, cex.sub = 0.6)
