library(lubridate)
library(tidyverse)
library(magrittr)
library(dplyr)
library(data.table)

first_date <- as.Date("2022-04-10")
last_date <- as.Date("2022-06-30")

#Get case data - this takes ages 
all_data <- fread("https://cog-uk.s3.climb.ac.uk/phylogenetics/latest/cog_metadata.csv.gz")
our_data <- select(filter(all_data, sample_date >= first_date), lineage, sample_date)
#Keep all lineage names and select only the pure ones.

grp = group_by(our_data, sample_date, lineage) %>% summarise(number = n())

a = select(filter(grp, lineage == "BA.2"), sample_date, number)
b = select(filter(grp, lineage == "BA.4"), sample_date, number)
c = select(filter(grp, lineage == "BA.5"), sample_date, number)
d = select(filter(grp, lineage == "BA.2.12.1"), sample_date, number)
e = select(filter(grp, lineage == "BA.5.1"), sample_date, number)


#Collect counts of each lineage in weighting by total of BA sequences
df <- purrr::reduce(list(a,b,c,d,e), dplyr::left_join, by = 'sample_date')
colnames(df) = c("sample_date", "ba2", "ba4", "ba5", "ba2121", "ba51")
df <- mutate(df, ba4_frac = ba4 / ba2)
df <- mutate(df, ba5_frac = ba5 / ba2)
df <- mutate(df, ba2121_frac = ba2121 / ba2)
df <- mutate(df, ba51_frac = ba51 / ba2)
wtsba4 = coalesce(1/(1/df$ba2 + 1/df$ba4),0)
wtsba5 = coalesce(1/(1/df$ba2 + 1/df$ba5),0)
wtsba2121 = coalesce(1/(1/df$ba2121 + 1/df$ba5),0)
wtsba51 = coalesce(1/(1/df$ba51 + 1/df$ba51),0)

#Regression with log2 transformation
regba4 <- lm(log2(ba4_frac) ~ date(sample_date), data = df, weights = wtsba4)
regba5 <- lm(log2(ba5_frac) ~ date(sample_date), data = df, weights = wtsba5)
regba2121 <- lm(log2(ba2121_frac) ~ date(sample_date), data = df, weights = wtsba2121)
regba51 <- lm(log2(ba51_frac) ~ date(sample_date), data = df, weights = wtsba51)


plot(date(df$sample_date), log2(df$ba4_frac), pch = 19, col = 2, xlim = c(first_date, last_date), ylim = c(-12, 2), yaxt = "n", xaxt = "n", ylab = "Ratio to BA.2, % - log2 scale", xlab = "Date", main = "BA.4, BA.5, BA.5.1 & BA.2.12.1 as a ratio to BA.2 in COG-UK sequences")
points(date(df$sample_date), pch = 19, log2(df$ba2121_frac), col = 4)
points(date(df$sample_date), pch = 19, log2(df$ba5_frac), col = 3)
points(date(df$sample_date), pch = 19, log2(df$ba51_frac), col = 5)


yticks = seq(-12, 2)
axis(side = 2, at = yticks, labels = c("0.024%","0.049%", "0.098%", "0.195%", "0.391%", "0.781%", "1.563%", "3.125%", "6.25%", "12.5%", "25%", "50%", "1x", "2x", "4x"))
xticks = seq(first_date, last_date, 7)
axis(1, at = xticks, labels = format(date(xticks), format = "%d %b"))

abline(regba4, col = 2)
abline(regba5, col = 3)
abline(regba2121, col = 4)
abline(regba51, col = 5)
abline(h = 0, col = "gray", lty= "dotted")
abline(h = -1, col = "gray", lty= "dashed")
legend("bottomright", inset = 0.05, legend = c("BA.4", "BA.5", "BA.2.12.1", "BA.5.1"), col = c(2, 3, 4, 5), lty = c(1, 1))
title(sub="Plot of COG-UK data, by @paulmainwood, based on work from @alexselby1770. Regression weighted by inverse probability of selection; last few data points based on very small values", adj=1, line=4, font=1, cex.sub = 0.6)

