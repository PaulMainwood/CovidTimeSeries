#Creates a plot of proportions of various selected variants vs a baseline variant
#After http://sonorouschocolate.com/covid19/index.php/UK_variant_comparison
#by @alexselby

library(tidyverse)
library(magrittr)
library(data.table)
library(lubridate)

first_date <- as.Date("2022-06-15")
last_date <- as.Date("2022-08-22")

#Get case data - this takes ages 
all_data <- fread("https://cog-uk.s3.climb.ac.uk/phylogenetics/latest/cog_metadata.csv.gz")
our_data <- select(filter(all_data, sample_date >= first_date), lineage, sample_date)
#Keep all lineage names and select only the pure ones.

grp = group_by(our_data, sample_date, lineage) %>% summarise(number = n())

#Function to take new and base variant names and produce dataframes of one as 
#fraction of the other, and regress the log of the proportions by date to give 
#a best-fit line.
proportions_new_variant <- function(new_variant_name, base_variant_name, grp) {
  
  base <- grp %>% filter(str_detect(lineage, base_variant_name)) %>% summarise(number = sum(number)) %>% select(sample_date, number)
  new <- grp %>% filter(str_detect(lineage, new_variant_name)) %>% summarise(number = sum(number)) %>% select(sample_date, number)
  
  #Collect counts of each lineage in weighting by total of BA sequences
  df <- purrr::reduce(list(base, new), dplyr::left_join, by = 'sample_date')
  colnames(df) = c("sample_date", "base", "new")

  df <- mutate(df, new_frac = new / base, new_frac_log2 = log2(new / base), new_variant_name = new_variant_name, base_variant_name = base_variant_name)
  wts <- coalesce(1/(1/df$base + 1/df$new), 0)
  df <- mutate(df, wts = wts)

  #Regression with log2 transformation
  reg <- lm(log2(new_frac) ~ date(sample_date), weights = wts, data = df)
  
  return(list(df = df, reg = reg))
  }


#Use the function to grab whichever variants and base you want and plot. Manual
#awkwardness here, since when automate, get weird visuals.

base_variant = "BA.5.2"
new_variants = c("BA.4.6", "BA.5.1", "BF.5", "BA.2.75")
dfs <- sapply(new_variants, function(x) proportions_new_variant(x, base_variant, grp))
df <- bind_rows(dfs[1, ])

intercepts <- sapply(seq(1:length(dfs[2,])), function(x) dfs[2,][[x]]$coefficients[[1]])
slopes <- sapply(seq(1:length(dfs[2,])), function(x) dfs[2,][[x]]$coefficients[[2]])
ablines <- tibble(intercept = intercepts, slope = slopes, cols = new_variants)


ggplot(df, aes(x = sample_date, y = new_frac_log2, color = new_variant_name)) + 
  geom_point(size = 3, aes(alpha = base)) + 
  scale_alpha(guide = "none") +
  scale_x_date(date_breaks = "1 week", date_labels = "%d %b", limits = c(first_date, last_date)) + 
  xlab("Sample date") +
  scale_y_continuous(limits = c(-10, 1), breaks = seq(-10, 1), labels = c("0.10%","0.20%", "0.39%", "0.78%", "1.56%", "3.12%", "6.25%", "12.5%", "25%", "50%", "1x", "2x")) +
  ylab(paste0("Fraction of ", base_variant, " & subvariants, log2 scale")) +
  theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
  geom_hline(linetype = "dashed", colour = "black", yintercept = 0) +
  geom_abline(data = ablines, aes(intercept = intercept, slope = slope, colour = cols)) +
  labs(title = paste0("Selected variant frequency in COG-UK data as share of ", base_variant, " and subvariants"), caption = "Plot of COG-UK data, by @paulmainwood, based on earlier work from @alexselby1770. Regression weighted by inverse probability of selection; last few data points based on very small values, shown by fading") +
  theme(plot.caption = element_text(size = 7), plot.title = element_text(size = 16, hjust = 0.5), legend.position = c(0.94, 0.1), legend.title=element_blank())
