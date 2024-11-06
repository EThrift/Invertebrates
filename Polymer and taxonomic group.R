#created on 22/08/24
# By Emily Thrift
# script for a glm to assess the levels of polymer levels in both the taxonomic groups and trophic levels of different invertebrate samples tested 

#load libraries
library(tidyverse)
library(ggplot2)
library(lme4)
library(glmmTMB)
library(emmeans)


#load dataframe 

df_polymer <- read.csv("Polymer.csv")


# Remove spaces from column names
colnames(df_polymer) <- make.names(colnames(df_polymer))




df_taxonomic <- read.csv ("Taxonomic_group.csv")
colnames(df_taxonomic) <- make.names(colnames(df_taxonomic))

df_taxonomic

library(binom)

wilson_ci2 <- binom.confint(x = df_taxonomic$n, n = df_taxonomic$Total, method = "wilson")

# Add confidence intervals to the data using mutate
df_taxonomic <- df_taxonomic %>%
  mutate(lower = wilson_ci2$lower,
         upper = wilson_ci2$upper)
#
# Calculate label positions just above error bars
df_taxonomic$label_y <- df_taxonomic$upper * 100 + 3.5

# Summarize the total for each taxonomic group and create the label with "N = "
totals2 <- df_taxonomic %>% group_by(Taxonomic.group) %>% summarize(Total = unique(Total))
totals2$TotalLabel <- paste("N =", totals2$Total)

totals2

library(ggplot2)
library(ggpattern)

# Convert Polymer to a factor
df_taxonomic$Polymer <- factor(df_taxonomic$Polymer)

# Create a summary of mean percentages for ordering
mean_order <- df_taxonomic %>%
  group_by(Polymer) %>%
  summarize(mean_percentage = mean(Percentage, na.rm = TRUE)) %>%
  arrange(desc(mean_percentage)) %>%
  pull(Polymer)

ggplot(df_taxonomic, aes(x = Polymer, y = Percentage)) +  # Reorder points by Percentage in descending order
  geom_point(position = position_dodge(width = 0.8), size = 3, color = "black") +  # Use geom_point for scatter plot
  geom_errorbar(aes(ymin = lower * 100, ymax = upper * 100), 
                position = position_dodge(width = 0.8), width = 0.25) +
  geom_text(aes(y = label_y + 5, label = n),  # Adjust the y position for data labels
            position = position_dodge(width = 0.8), size = 3.5, color = "black") +
  labs(x = "Polymer Type", y = "% of Polymer") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),  # Adjust angle and size for better visibility
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",  # Remove the legend
    plot.margin = margin(1, 1, 1, 1, "cm")  # Increase the bottom margin
  ) +
  facet_wrap(~ Taxonomic.group, scales = "free_x", ncol = 6) 

colnames(df_taxonomic)


# Reshape data into long format for easier analysis
df_long <- df_polymer %>%
  pivot_longer(cols = 13:24, names_to = "Polymer", values_to = "Presence") %>%
  filter(Presence == 1)  # Keep only presence data

# Step 1: Calculate Simpson Diversity Index by taxonomic group
diversity_indices_simpson <- df_long %>%
  group_by(Taxonomic.group, Polymer) %>%
  summarise(Presence = sum(Presence), .groups = 'drop') %>%
  group_by(Taxonomic.group) %>%
  summarise(
    Simpson = 1 - sum((Presence / sum(Presence))^2),
    Richness = n()
  )

kruskal_test <- kruskal.test(Simpson ~ Taxonomic.group, data = diversity_indices_simpson)

# Extract the test statistic, degrees of freedom, and p-value
test_statistic <- kruskal_test_result$statistic
degrees_of_freedom <- kruskal_test_result$parameter
p_value <- kruskal_test_result$p.value

# Print the values
cat("Kruskal-Wallis chi-squared =", test_statistic, ", df =", degrees_of_freedom, ", p-value =", p_value, "\n")

#result has a p value of 0.41

# Ensure Taxonomic.group is a factor
diversity_indices_simpson$Taxonomic.group <- as.factor(diversity_indices_simpson$Taxonomic.group)


#no significant differences seen in the pairwise analysis

# Plotting Simpson Diversity Index across taxonomic groups
ggplot(diversity_indices_simpson, aes(x = Taxonomic.group, y = Simpson, fill = Taxonomic.group)) +
  geom_bar(stat = "identity", color = "black", show.legend = FALSE) +  # Bar plot
  labs(title = "Simpson Diversity Index of Polymers by Taxonomic Group",
       x = "Taxonomic Group", 
       y = "Simpson Diversity Index") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
# Check if Taxonomic.group is a factor
diversity_indices_simpson$Taxonomic.group <- as.factor(diversity_indices_simpson$Taxonomic.group)

# Run the Kruskal-Wallis test
kruskal_test_result <- kruskal.test(Simpson ~ Taxonomic.group, data = diversity_indices_simpson)

# Print the Kruskal-Wallis test result
# Extract the test statistic, degrees of freedom, and p-value
test_statistic <- kruskal_test_result$statistic
degrees_of_freedom <- kruskal_test_result$parameter
p_value <- kruskal_test_result$p.value

# Print the values
cat("Kruskal-Wallis chi-squared =", test_statistic, ", df =", degrees_of_freedom, ", p-value =", p_value, "\n")


library(FSA)

# Now run Dunn's test again
dunn_test_results <- dunnTest(Simpson ~ Taxonomic.group, data = diversity_indices_simpson, method = "bonferroni")

# Print the results
print(dunn_test_results)


# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(FactoMineR)
library(factoextra)

df_polymer
# Summarize polymer data per taxonomic group (wide format)
df_pca <- df_polymer %>%
  group_by(Taxonomic.group) %>%
  summarise(across(12:24, sum))  # Sum the presence/absence for each polymer

df_pca

# Perform PCA on the polymer columns
pca_result <- PCA(df_pca[, -1], graph = FALSE)  # Exclude the first column (Taxonomic.group)

# Check if Taxonomic.group exists and is a factor or character
df_pca$Taxonomic.group <- as.factor(df_pca$Taxonomic.group)

# Run PCA excluding the first column (Taxonomic.group), assuming the remaining columns are polymer data
pca_result <- PCA(df_pca[, -1], graph = FALSE)

pca_result

library(ggrepel)

# Run PCA excluding the first column (Taxonomic.group), assuming the remaining columns are polymer data
pca_result <- PCA(df_pca[, -1], graph = FALSE)

pca_result
# Extract the PCA variable coordinates and labels
pca_coords <- as.data.frame(pca_result$var$coord)  # Variable coordinates
arrow_labels <- rownames(pca_coords)  

# Check the content of pca_coords to ensure it's created correctly
print(head(pca_coords))

# Calculate label positions near the end of the arrows
pca_coords$label_x <- pca_coords$Dim.1 * 4.2  # Position label at 90% of the arrow length on x-axis
pca_coords$label_y <- pca_coords$Dim.2 * 4.2  # Position label at 90% of the arrow length on y-axis

# Plot with labels positioned near the arrowheads
fviz_pca_biplot(
  pca_result,
  label = "none",                      # Do not show default variable labels
  habillage = df_pca$Taxonomic.group,    # Color points by Trophic level
  addEllipses = TRUE,                  # Add concentration ellipses
  palette = "Set1"                     # Set color palette
) + 
  geom_point(aes(color = df_pca$Taxonomic.group), size = 4, shape = 16) +  # Round points only
  scale_color_brewer(palette = "Set1") +  # Use a color palette
  theme_minimal() +                       # Use a minimal theme
  theme(legend.position = "right") +      # Position the legend
  geom_text(
    data = pca_coords,                    # Use updated variable coordinates for label placement
    aes(x = label_x,                      # Adjusted x coordinate near the arrowhead
        y = label_y,                      # Adjusted y coordinate near the arrowhead
        label = arrow_labels),            # Use variable names for labels
    size = 4,                             # Adjust text size for clarity
    hjust = 0.5,                          # Center text horizontally
    vjust = -0.5                          # Adjust vertical alignment to avoid overlap with arrowhead
  )
