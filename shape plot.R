 # Load necessary libraries
library(tidyverse)
library(MASS)
library(multcomp)
library(binom)

# Read the data

df_brand <- read_csv("Shape_plot.csv")

# Fix column names
colnames(df_brand) <- make.names(colnames(df_brand))

df_animal <- read_csv("Shape_plot_animal.csv")

# Fix column names
colnames(df_animal) <- make.names(colnames(df_animal))

df_type <- read_csv("Shape_plot_type.csv")

# Fix column names
colnames(df_type) <- make.names(colnames(df_type))
# Remove rows with missing or zero values in n and Total
df_brand <- df_brand %>%
  filter(!is.na(n) & !is.na(Total) & n > 0 & Total > 0 & n <= Total)


# Calculate Wilson confidence intervals
wilson_ci2 <- binom.confint(x = df_brand$n, n = df_brand$Total, method = "wilson")

# Add confidence intervals to the data frame
df_brand <- df_brand %>%
  mutate(lower = wilson_ci2$lower,
         upper = wilson_ci2$upper)

# Calculate label positions just above error bars
df_brand <- df_brand %>%
  mutate(label_y = upper * 100 + 3.5) 
# Summarize the total for each Taxonomic.group and add it as a separate column
# Summarize totals for each Taxonomic group
totals <- df_brand %>%
  group_by(Brand) %>%
  summarize(Total = unique(Total), .groups = 'drop')

# Join totals back to the original data frame
df_brand <- df_brand %>%
  left_join(totals, by = "Brand")

# Create the plot with single labels for totals under each group
ggplot(df_brand, aes(x = Brand, y = Percentage, fill = Shape)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black", width = 0.7) +  # Clustered bars with black borders
  geom_errorbar(aes(ymin = lower * 100, ymax = upper * 100), 
                width = 0.3, position = position_dodge(width = 0.8)) +  # Error bars for each bar
  geom_text(aes(y = label_y, label = n),  # Show "n" at the top of each bar
            size = 4, position = position_dodge(width = 0.8), vjust = -0.5, fontface = "plain") +  # Text labels above bars
  # Add total label below each group, setting inherit.aes = FALSE to avoid issues with missing Shape
  geom_text(data = totals, aes(x = Brand, y = -5, label = paste("N =", Total)), 
            vjust = 1.5, size = 3.5, inherit.aes = FALSE) +  # Adjust position below x-axis
  scale_fill_manual(values = c("Fibre" = "grey80", "Fragment" = "white")) +  # Set colors for "Fibre" and "Fragment"
  labs(x = "Taxonomic group", y = "% Presence") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10, angle = 90, hjust = 1),  # Rotate x-axis text
    axis.text.y = element_text(size = 10),  # Adjust y-axis text size
    axis.title.x = element_text(size = 14),  # Adjust x-axis title size
    axis.title.y = element_text(size = 14),  # Adjust y-axis title size
    legend.position = "right",               # Position legend on the right
    panel.grid = element_blank(),            # Removes default grid lines
    axis.line = element_line(color = "black") # Adds black axis lines
  ) + 
  scale_y_continuous(limits = c(0, 100))


# Remove rows with missing or zero values in n and Total
df_animal <- df_animal %>%
  filter(!is.na(n) & !is.na(Total) & n > 0 & Total > 0 & n <= Total)




# Calculate Wilson confidence intervals on the cleaned data
wilson_ci3 <- binom.confint(x = df_animal$n, n = df_animal$Total, method = "wilson")

# Add confidence intervals to the cleaned data frame
df_animal <- df_animal %>%
  mutate(lower = wilson_ci3$lower,
         upper = wilson_ci3$upper)

# Calculate label positions just above error bars
df_animal$label_y <- df_animal$upper * 100 + 3.5

# Summarize totals for each Taxonomic group
totals <- df_animal %>%
  group_by(Target.animal) %>%
  summarize(Total = unique(Total), .groups = 'drop')

# Join totals back to the original data frame
df_animal <- df_animal %>%
  left_join(totals, by = "Target.animal")

# Create the plot with single labels for totals under each group
ggplot(df_animal, aes(x = Target.animal, y = Percentage, fill = Shape)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black", width = 0.7) +  # Clustered bars with black borders
  geom_errorbar(aes(ymin = lower * 100, ymax = upper * 100), 
                width = 0.3, position = position_dodge(width = 0.8)) +  # Error bars for each bar
  geom_text(aes(y = label_y, label = n),  # Show "n" at the top of each bar
            size = 4, position = position_dodge(width = 0.8), vjust = -0.5, fontface = "plain") +  # Text labels above bars
  # Add total label below each group, setting inherit.aes = FALSE to avoid issues with missing Shape
  geom_text(data = totals, aes(x = Target.animal, y = -5, label = paste("N =", Total)), 
            vjust = 1.5, size = 4, inherit.aes = FALSE) +  # Adjust position below x-axis
  scale_fill_manual(values = c("Fibre" = "grey80", "Fragment" = "white")) +  # Set colors for "Fibre" and "Fragment"
  labs(x = "Target animal", y = "% Presence") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12),  # Adjust x-axis text size
    axis.text.y = element_text(size = 12),  # Adjust y-axis text size
    axis.title.x = element_text(size = 14),  # Adjust x-axis title size
    axis.title.y = element_text(size = 14),  # Adjust y-axis title size
    legend.position = "right",# Position legend on the right Position legend on the right,     
    panel.grid = element_blank(),           # Removes default grid lines
    axis.line = element_line(color = "black") # Adds black axis lines
  ) + 
  scale_y_continuous(limits = c(0, 100))

# Print rows with invalid values (if any)
df_type %>% filter(is.na(n) | is.na(Total) | n < 0 | Total <= 0 | n > Total)

# Remove rows where Type is Omnivore
df_type <- df_type %>%
  filter(Type != "Omnivore")



# Calculate Wilson confidence intervals on the cleaned data
wilson_ci3 <- binom.confint(x = df_type$n, n = df_type$Total, method = "wilson")

# Add confidence intervals to the cleaned data frame
df_type <- df_type %>%
  mutate(lower = wilson_ci3$lower,
         upper = wilson_ci3$upper)

# Calculate label positions just above error bars
df_type$label_y <- df_type$upper * 100 + 3.5

# Summarize totals for each Taxonomic group
totals <- df_type %>%
  group_by(Type) %>%
  summarize(Total = unique(Total), .groups = 'drop')

# Join totals back to the original data frame
df_type <- df_type %>%
  left_join(totals, by = "Type")

# Create the plot with single labels for totals under each group
ggplot(df_type, aes(x = Type, y = Percentage, fill = Shape)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black", width = 0.7) +  # Clustered bars with black borders
  geom_errorbar(aes(ymin = lower * 100, ymax = upper * 100), 
                width = 0.3, position = position_dodge(width = 0.8)) +  # Error bars for each bar
  geom_text(aes(y = label_y, label = n),  # Show "n" at the top of each bar
            size = 4, position = position_dodge(width = 0.8), vjust = -0.5, fontface = "plain") +  # Text labels above bars
  # Add total label below each group, setting inherit.aes = FALSE to avoid issues with missing Shape
  geom_text(data = totals, aes(x = Type, y = -5, label = paste("N =", Total)), 
            vjust = 1.5, size = 4, inherit.aes = FALSE) +  # Adjust position below x-axis
  scale_fill_manual(values = c("Fibre" = "grey80", "Fragment" = "white")) +  # Set colors for "Fibre" and "Fragment"
  labs(x = "Type", y = "% Presence") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12),  # Adjust x-axis text size
    axis.text.y = element_text(size = 12),  # Adjust y-axis text size
    axis.title.x = element_text(size = 14),  # Adjust x-axis title size
    axis.title.y = element_text(size = 14),  # Adjust y-axis title size
    legend.position = "right",# Position legend on the right Position legend on the right,     
    panel.grid = element_blank(),           # Removes default grid lines
    axis.line = element_line(color = "black") # Adds black axis lines
  ) + 
  scale_y_continuous(limits = c(0, 100))
