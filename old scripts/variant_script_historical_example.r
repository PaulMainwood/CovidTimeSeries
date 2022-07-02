library(lubridate)
library(tidyverse)
library(magrittr)
library(dplyr)
library(data.table)

first_date <- as.Date("2021-12-15")
last_date <- as.Date("2022-04-1")

#Get case data - this takes ages 
all_data <- fread("https://cog-uk.s3.climb.ac.uk/phylogenetics/latest/cog_metadata.csv.gz")
our_data <- select(filter(all_data, sample_date >= first_date), lineage, sample_date)
#Due to a notational quirk, we can collect lineages by just taking the first four characters of the lineage name.
our_data <- mutate(our_data, lineage = substr(lineage, 1, 4))

grp = group_by(our_data, sample_date, lineage) %>% summarise(number = n())

a = select(filter(grp, lineage == "BA.1"), sample_date, number)
b = select(filter(grp, lineage == "BA.2"), sample_date, number)

#Collect counts of each lineage in weighting by total of BA sequences
df <- purrr::reduce(list(a,b), dplyr::left_join, by = 'sample_date')
colnames(df) = c("sample_date", "ba1", "ba2")
df <- mutate(df, ba2_frac = ba2 / ba1)
df$wts = coalesce(df$ba1,0) + coalesce(df$ba2,0)

#Regression with log2 transformation until 12 Jan

df_to_10Jan = filter(df, sample_date < as.Date("2022-01-10"))
regba2 <- lm(log2(ba2_frac) ~ date(sample_date), data = df_to_10Jan, weights = df_to_10Jan$wts)


plot(date(df$sample_date), log2(df$ba2_frac), col = 2, ylim = c(-13, 4), xlim = c(first_date, last_date), yaxt = "n", xaxt = "n", ylab = "Ratio of BA.2 to BA.1, % - log2 scale", xlab = "Date", main = "BA.2 as a ratio to BA.1 in COG-UK sequences")

yticks = seq(-13, 3)
axis(side = 2, at = yticks, labels = c("0.012%", "0.024%","0.049%", "0.098%", "0.195%", "0.391%", "0.781%", "1.563%", "3.125%", "6.25%", "12.5%", "25%", "50%", "100%", "2x", "4x", "8x"))
xticks = seq(first_date, last_date, 7)
axis(1, at = xticks, labels = format(date(xticks), format = "%d %b"))

abline(regba2, col = 2)
abline(h = 0, col = "gray", lty= "dotted")
abline(h = -1, col = "gray", lty= "dashed")
title(sub="Plot of COG-UK data, by @paulmainwood, inspired by previous work from @alexselby1770", adj=1, line=4, font=1, cex.sub = 0.6)

#Cases prediction - using ONS prevalence
a = data.frame(sample_date = seq(first_date, last_date+10, 1))
a$predictionba4 = 2^(predict(regba4,a))
a$predictionba5 = 2^(predict(regba5,a))
may7 = 78000 #estimate from ONS incidence at similar point to current prevalence
ba2growth <- may7 * (1 + -0.022)^(seq(0:39))
b <- filter(a, sample_date >= as.Date("2022-05-07"))
b$ba2abs = ba2growth
b$ba4abs = b$predictionba4 * b$ba2
b$ba5abs = b$predictionba5 * b$ba2
b$total = b$ba2abs + b$ba4abs + b$ba5abs
plot(b$sample_date, b$total)
