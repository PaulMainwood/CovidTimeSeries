library(lubridate)
library(tidyverse)
library(magrittr)
library(dplyr)
library(data.table)

first_date <- as.Date("2022-04-10")
last_date <- as.Date("2022-06-05")

#Get case data - this takes ages 
all_data <- fread("https://cog-uk.s3.climb.ac.uk/phylogenetics/latest/cog_metadata.csv.gz")
our_data <- select(filter(all_data, sample_date >= first_date), lineage, sample_date)
our_data <- mutate(our_data, lineage = substr(lineage, 1, 4))

grp = group_by(our_data, sample_date, lineage) %>% summarise(number = n())

a = select(filter(grp, lineage == "BA.2"), sample_date, number)
b = select(filter(grp, lineage == "BA.4"), sample_date, number)
c = select(filter(grp, lineage == "BA.5"), sample_date, number)

#Note, all unweighted (so a day with 1 sequence counts as much as one with 100)
df <- purrr::reduce(list(a,b,c), dplyr::left_join, by = 'sample_date')
colnames(df) = c("sample_date", "ba2", "ba4", "ba5")
df <- mutate(df, ba4_frac = ba4 / (ba2 + ba4 + ba5))
df <- mutate(df, ba5_frac = ba5 / (ba2 + ba4 + ba5))

regba4 <- lm(log2(ba4_frac) ~ date(sample_date), data = df)
#regba5 <- lm(log2(ba5_frac) ~ date(sample_date), data = df)

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

