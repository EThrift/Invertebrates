# Load necessary libraries
library(tidyverse)
library(MASS)
library(multcomp)
library(binom)

# Read the data

df_taxonomic_shape <- read_csv("taxonomic_shape.csv")


# Fix column names
colnames(df_taxonomic_shape) <- make.names(colnames(df_taxonomic_shape))


# Remove rows with missing or zero values in n and Total
df_taxonomic_shape <- df_taxonomic_shape %>%
  filter(!is.na(n) & !is.na(Total) & n > 0 & Total > 0 & n <= Total)


# Calculate Wilson confidence intervals
wilson_ci2 <- binom.confint(x = df_taxonomic_shape$n, n = df_taxonomic_shape$Total, method = "wilson")

# Add confidence intervals to the data frame
df_taxonomic_shape <- df_taxonomic_shape %>%
  mutate(lower = wilson_ci2$lower,
         upper = wilson_ci2$upper)

# Calculate label positions just above error bars
df_taxonomic_shape <- df_taxonomic_shape %>%
  mutate(label_y = upper * 100 + 3.5) 
# Summarize the total for each Taxonomic.group and add it as a separate column
# Summarize totals for each Taxonomic group
totals <- df_taxonomic_shape %>%
  group_by(Taxonomic.group) %>%
  summarize(Total = unique(Total), .groups = 'drop')

# Join totals back to the original data frame
df_taxonomic_shape <- df_taxonomic_shape %>%
  left_join(totals, by = "Taxonomic.group")

# Create the plot with single labels for totals under each group
ggplot(df_taxonomic_shape, aes(x = Taxonomic.group, y = Percentage, fill = Shape)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black", width = 0.7) +  # Clustered bars with black borders
  geom_errorbar(aes(ymin = lower * 100, ymax = upper * 100), 
                width = 0.3, position = position_dodge(width = 0.8)) +  # Error bars for each bar
  geom_text(aes(y = label_y, label = n),  # Show "n" at the top of each bar
            size = 4, position = position_dodge(width = 0.8), vjust = -0.5, fontface = "plain") +  # Text labels above bars
  # Add total label below each group, setting inherit.aes = FALSE to avoid issues with missing Shape
  geom_text(data = totals, aes(x = Taxonomic.group, y = -5, label = paste("N =", Total)), 
            vjust = 1.5, size = 4, inherit.aes = FALSE) +  # Adjust position below x-axis
  scale_fill_manual(values = c("Fibre" = "grey80", "Fragment" = "white")) +  # Set colors for "Fibre" and "Fragment"
  labs(x = "Taxonomic group", y = "% Presence") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12),  # Adjust x-axis text size
    axis.text.y = element_text(size = 12),  # Adjust y-axis text size
    axis.title.x = element_text(size = 14),  # Adjust x-axis title size
    axis.title.y = element_text(size = 14),  # Adjust y-axis title size
    legend.position = "right"                # Position legend on the right
  )
df_feeding_shape <- read.csv("feeding_shape.csv")

# Fix column names
colnames(df_feeding_shape) <- make.names(colnames(df_feeding_shape))

# Check for missing values
sum(is.na(df_feeding_shape$n))
sum(is.na(df_feeding_shape$Total))

# Check for invalid values
any(df_feeding_shape$n < 0)
any(df_feeding_shape$Total <= 0)
any(df_feeding_shape$n > df_feeding_shape$Total)

# Print rows with invalid values (if any)
df_feeding_shape %>% filter(is.na(n) | is.na(Total) | n < 0 | Total <= 0 | n > Total)

# Remove rows where Feeding.group is Omnivore
df_feeding_shape <- df_feeding_shape %>%
  filter(Feeding.group != "Omnivore")



# Calculate Wilson confidence intervals on the cleaned data
wilson_ci3 <- binom.confint(x = df_feeding_shape$n, n = df_feeding_shape$Total, method = "wilson")

# Add confidence intervals to the cleaned data frame
df_feeding_shape <- df_feeding_shape %>%
  mutate(lower = wilson_ci3$lower,
         upper = wilson_ci3$upper)

# Calculate label positions just above error bars
df_feeding_shape$label_y <- df_feeding_shape$upper * 100 + 3.5

# Summarize totals for each Taxonomic group
totals <- df_feeding_shape %>%
  group_by(Feeding.group) %>%
  summarize(Total = unique(Total), .groups = 'drop')

# Join totals back to the original data frame
df_feeding_shape <- df_feeding_shape %>%
  left_join(totals, by = "Feeding.group")

# Create the plot with single labels for totals under each group
ggplot(df_feeding_shape, aes(x = Feeding.group, y = Percentage, fill = Shape)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black", width = 0.7) +  # Clustered bars with black borders
  geom_errorbar(aes(ymin = lower * 100, ymax = upper * 100), 
                width = 0.3, position = position_dodge(width = 0.8)) +  # Error bars for each bar
  geom_text(aes(y = label_y, label = n),  # Show "n" at the top of each bar
            size = 4, position = position_dodge(width = 0.8), vjust = -0.5, fontface = "plain") +  # Text labels above bars
  # Add total label below each group, setting inherit.aes = FALSE to avoid issues with missing Shape
  geom_text(data = totals, aes(x = Feeding.group, y = -5, label = paste("N =", Total)), 
            vjust = 1.5, size = 4, inherit.aes = FALSE) +  # Adjust position below x-axis
  scale_fill_manual(values = c("Fibre" = "grey80", "Fragment" = "white")) +  # Set colors for "Fibre" and "Fragment"
  labs(x = "Trophic level", y = "% Presence") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12),  # Adjust x-axis text size
    axis.text.y = element_text(size = 12),  # Adjust y-axis text size
    axis.title.x = element_text(size = 14),  # Adjust x-axis title size
    axis.title.y = element_text(size = 14),  # Adjust y-axis title size
    legend.position = "right"                # Position legend on the right
  )
