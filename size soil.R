#created on 21/08/2024
#created by Emily Thrift


# Load necessary libraries
library(dplyr)
library(lme4) 
library(tidyverse)
library(ggplot2)
library(glmmTMB)

df_size <- read_csv("Soil size.csv")


# Remove spaces from column names
colnames(df_size) <- make.names(colnames(df_size))


model1 <- lmer(Size ~  (1 | Site),  data = df_size)
summary(model1)

model2 <- lmer(Size ~ Land.cover  + (1 | Site),  data = df_size)
summary(model2)

anova(model1,model2)
anova(model2)

#model 1 lower AIC  showing land cover not significant

# Load the DHARMa package
library(DHARMa)

# Simulate residuals
sim_res <- simulateResiduals(fittedModel = model1)

# Plot diagnostics
plot(sim_res)
#QQplot shows good fit

# Overdispersion test using DHARMa
testDispersion(sim_res)

# Create a box plot of Size by Taxonomic.group

library(ggplot2)

#have a function which means you can just call this when you want to calculate standard error 
std_err <- function(x) {
  SE <- sd(x)/sqrt(length(x))
  return(SE)
}

# Filter the dataset to exclude Size values over 2
df_size_filtered <- df_size %>%
  filter(Size <= 2)

# Calculate mean and standard deviation for the filtered dataset
df_size_summary <- df_size_filtered %>%
  group_by(Land.cover) %>%
  summarise(size_mean = mean(Size),
            size_SD = sd(Size))

# Create the violin plot with mean as a dot and standard deviation as error bars
ggplot(df_size_filtered, aes(x = Land.cover, y = Size)) +
  geom_violin(fill = alpha("grey80", 0.4), color = "black") +  # Violin plot with lighter fill
  geom_point(data = df_size_summary, aes(x = Land.cover, y = size_mean), 
             color = "blue", size = 3) +                        # Blue dot for the mean
  geom_errorbar(data = df_size_summary, aes(x = Land.cover, y = size_mean, 
                                            ymin = size_mean - size_SD, 
                                            ymax = size_mean + size_SD),
                width = 0.2, color = "blue") +                 # Blue error bars for SD
  theme_minimal() +
  labs(x = "Land cover", y = "Size (mm)") +
  theme(
    axis.text.x = element_text(size = 14),       # Larger text for x-axis
    axis.text.y = element_text(size = 14),       # Larger text for y-axis
    axis.title.x = element_text(size = 16),      # Larger title for x-axis
    axis.title.y = element_text(size = 16),      # Larger title for y-axis
    panel.grid = element_blank(),                # Removes default grid lines
    axis.line = element_line(color = "black")    # Adds black axis lines
  ) +
  coord_cartesian(ylim = c(0, 2))                # Flexible y-axis limits to include SD



