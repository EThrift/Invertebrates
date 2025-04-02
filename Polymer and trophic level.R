#created on 22/08/24
# By Emily Thrift
# script for a glm to assess the levels of polymer levels in trophic levels of di fferent invertebrate samples tested 

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

df_trophic <- read.csv("feeding_group.csv")
colnames(df_trophic) <- make.names(colnames(df_trophic))

df_trophic_talk <- read.csv("feeding_group_talk.csv")
colnames(df_trophic_talk) <- make.names(colnames(df_trophic_talk))
library(binom)
# Filter out rows where Trophic.level is "Omnivore"
df_trophic_talk <- df_trophic_talk %>%
  filter(Trophic.level != "Omnivore")

wilson_ci2 <- binom.confint(x = df_trophic$n, n = df_trophic$Total, method = "wilson")

# Add confidence intervals to the data using mutate
df_trophic <- df_trophic %>%
  mutate(lower = wilson_ci2$lower,
         upper = wilson_ci2$upper)
#
# Calculate label positions just above error bars
df_trophic$label_y <- df_trophic$upper * 100 + 3.5

# Summarize the total for each trophic group and create the label with "N = "
totals2 <- df_trophic %>% group_by(Trophic.level) %>% summarize(Total = unique(Total))
totals2$TotalLabel <- paste("N =", totals2$Total)

library(ggplot2)
library(ggpattern)


# Convert Polymer to a factor
df_trophic$Polymer <- factor(df_trophic$Polymer)

# Create a summary of mean percentages for ordering
mean_order <- df_trophic %>%
  group_by(Polymer) %>%
  summarize(mean_percentage = mean(Percentage, na.rm = TRUE)) %>%
  arrange(desc(mean_percentage)) %>%
  pull(Polymer)

ggplot(df_trophic, aes(x = Polymer, y = Percentage)) +  # Reorder points by Percentage in descending order
  geom_point(position = position_dodge(width = 0.8), size = 3, color = "black") +  # Use geom_point for scatter plot
  geom_errorbar(aes(ymin = lower * 100, ymax = upper * 100), 
                position = position_dodge(width = 0.8), width = 0.25) +
  geom_text(aes(y = label_y + 5, label = n),  # Adjust the y position for data labels
            position = position_dodge(width = 0.8), size = 3.5, color = "black") +
  labs(x = "Polymer Type", y = "% of Polymer") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 12, angle = 90, hjust = 1),  # Adjust angle and size for better visibility
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",  # Remove the legend
    plot.margin = margin(1, 1, 1, 1, "cm")  # Increase the bottom margin
  ) +
  ylim(0, 100) +
  facet_wrap(~ Trophic.level, scales = "free_x", ncol = 6) 

# Filter out rows where Trophic.level is "Omnivore"
df_trophic_talk <- df_trophic_talk %>%
  filter(Trophic.level != "Omnivore")

wilson_ci2 <- binom.confint(x = df_trophic_talk$n, n = df_trophic_talk$Total, method = "wilson")

# Add confidence intervals to the data using mutate
df_trophic_talk <- df_trophic_talk %>%
  mutate(lower = wilson_ci2$lower,
         upper = wilson_ci2$upper)
#
# Calculate label positions just above error bars
df_trophic_talk$label_y <- df_trophic_talk$upper * 100 + 3.5

# Summarize the total for each trophic group and create the label with "N = "
totals2 <- df_trophic_talk %>% group_by(Trophic.level) %>% summarize(Total = unique(Total))
totals2$TotalLabel <- paste("N =", totals2$Total)

# Define the custom color-blind friendly palette with 9 colors
color_blind_friendly_palette <- c(
  "Other" = "#1f77b4",        # Blue
  "Polyester" = "#ff7f0e",    # Orange
  "Polyethylene" = "#2ca02c", # Green
  "Polypropylene" = "#17becf",# Cyan
  "Polystyrene" = "#9467bd",  # Purple
  "Cellophane" = "#e377c2",   # Pink
  "Epoxy resin" = "red",  # Gray
  "Nylon" = "yellow",        # Yellow
  "Aluminium silicate" = "#8d8d8d"  # Dark Gray
)

# Ensure Polymer is a factor with "Other" as the first level in df_trophic_talk
df_trophic_talk$Polymer <- factor(df_trophic_talk$Polymer, 
                                  levels = c("Other", "Polyester", "Polyethylene", "Polypropylene", 
                                             "Polystyrene", "Cellophane", "Epoxy resin", 
                                             "Nylon", "Aluminium silicate"))

# Ensure Trophic.level is a factor in df_trophic_talk
df_trophic_talk$Trophic.level <- factor(df_trophic_talk$Trophic.level)

# Plot the data with custom colors and facet by Trophic.level
ggplot(df_trophic_talk, aes(x = Polymer, y = Percentage, fill = Polymer)) +  
  geom_col(position = position_dodge(width = 0.8), color = "black") +  # Bar chart
  geom_errorbar(aes(ymin = lower * 100, ymax = upper * 100), 
                position = position_dodge(width = 0.8), width = 0.25) +  # Error bars
  geom_text(aes(y = label_y + 5, label = paste(n)),  # Add "n = " before the n value
            position = position_dodge(width = 0.8), size = 2.5, color = "black") +  # Label count above bars
  labs(x = "Polymer Type in Trophic Levels", y = "% of Polymer") +  # Axis labels
  scale_fill_manual(values = color_blind_friendly_palette) +  # Apply custom color palette
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 12, angle = 90, hjust = 1),  # Rotate x-axis labels for readability
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",  # Remove the legend
    plot.margin = margin(1, 1, 1, 1, "cm"),  # Increase plot margin for spacing
    strip.text = element_text(size = 12, color = "black")  # Change facet label color
  ) +
  facet_wrap(~ Trophic.level, scales = "free_x", ncol = 6) +  # Facet by Trophic.level with free x scales
  ylim(0, 100)  # Set y-axis limits to go from 0 to 100


# Reshape data into long format for easier analysis
df_long <- df_polymer_trophic %>%
  pivot_longer(cols = 13:24, names_to = "Polymer", values_to = "Presence") %>%
  filter(Presence == 1)  # Keep only presence data

# Step 1: Calculate Simpson Diversity Index by taxonomic group
diversity_indices_simpson <- df_long %>%
  group_by(Trophic.level, Polymer) %>%
  summarise(Presence = sum(Presence), .groups = 'drop') %>%
  group_by(Trophic.level) %>%
  summarise(
    Simpson = 1 - sum((Presence / sum(Presence))^2),
    Richness = n()
  )

kruskal_test <- kruskal.test(Simpson ~ Trophic.level, data = diversity_indices_simpson)
# Extract the test statistic, degrees of freedom, and p-value
test_statistic <- kruskal_test$statistic
degrees_of_freedom <- kruskal_test$parameter
p_value <- kruskal_test$p.value

# Print the values
cat("Kruskal-Wallis chi-squared =", test_statistic, ", df =", degrees_of_freedom, ", p-value =", p_value, "\n")

#result has a p value of 0.41

# Ensure Taxonomic.group is a factor
diversity_indices_simpson$Trophic.level <- as.factor(diversity_indices_simpson$Trophic.level)


#no significant differences seen in the pairwise analysis

# Plotting Simpson Diversity Index across taxonomic groups
ggplot(diversity_indices_simpson, aes(x = Trophic.level, y = Simpson, fill = Trophic.level)) +
  geom_bar(stat = "identity", color = "black", show.legend = FALSE) +  # Bar plot
  labs(x = "Trophic.level", 
       y = "Simpson Diversity Index") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
# Check if Taxonomic.group is a factor
diversity_indices_simpson$Trophic.level <- as.factor(diversity_indices_simpson$Trophic.level)


# Now run Dunn's test again
dunn_test_results <- dunnTest(Simpson ~ Trophic.level, data = diversity_indices_simpson, method = "bonferroni")

# Print the results
print(dunn_test_results)


# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(FactoMineR)
library(factoextra)

# Summarize polymer data per taxonomic group (wide format)
df_mca <- df_polymer %>%
  group_by(Trophic.level) %>%
  summarise(across(12:24, sum))  # Adjust column indices if necessary

# Ensure Trophic.level is a factor
df_mca$Trophic.level <- as.factor(df_mca$Trophic.level)

# Convert remaining columns (polymer data) to factors for MCA
df_mca[, -1] <- lapply(df_mca[, -1], function(x) {
  as.factor(ifelse(x > 0, "Present", "Absent"))
})

# Perform MCA on the summarized data (excluding Trophic.level)
mca_result <- MCA(df_mca[, -1], graph = FALSE)

# Extract variable coordinates for labeling
mca_coords <- data.frame(mca_result$var$coord)
mca_coords$polymer <- rownames(mca_coords)  # Retain polymer names for labels

# Adjust positions for label placement
mca_coords$label_x <- mca_coords$Dim.1 * 1.2
mca_coords$label_y <- mca_coords$Dim.2 * 1.2


# Filter MCA coordinates to include only "Present" polymers
present_coords <- mca_coords %>%
  filter(grepl("Present", polymer)) %>%     # Keep only "Present" labels
  mutate(polymer = gsub("_Present", "", polymer))  # Remove "_Present" suffix from labels

# Create MCA biplot without red triangles but keep trophic level points

# Extract percentages of variance explained by each dimension
percent_var <- mca_result$eig[, 2]  # Second column contains the percentages


library(ggrepel)

# Load the saved axis limits from the first script
axis_limits <- readRDS("axis_limits.rds")

# Extract x and y limits
x_limit <- axis_limits$x_limit
y_limit <- axis_limits$y_limit

# Create MCA biplot for Trophic Level using the loaded axis limits for trophic levels make sure to load in axis from the taxonomic groups script first as it is preloaded
fviz_mca_biplot(
  mca_result,
  label = "none",                      
  habillage = df_mca$Trophic.level,    
  addEllipses = TRUE,                  
  ellipse.level = 0.95,                
  palette = "Set1",                    
  geom.var = "none",                   
  geom.ind = "point"                   
) +
  theme_minimal() +                    
  theme(
    legend.position = "right",         
    axis.text.x = element_text(size = 14),  
    axis.text.y = element_text(size = 14),  
    axis.title.x = element_text(size = 16), 
    axis.title.y = element_text(size = 16)  
  ) +
  labs(
    x = paste0("Dim 1 (", round(percent_var[1], 1), "%)"),
    y = paste0("Dim 2 (", round(percent_var[2], 1), "%)")
  ) +
  geom_text_repel(
    data = present_coords,             
    aes(x = label_x,                   
        y = label_y,                   
        label = polymer),              
    size = 4,                          
    max.overlaps = 50,                  
    box.padding = 0.5,                  
    point.padding = 0.5,                
    min.segment.length = 0.5,           
    segment.color = NA                 
  ) +
  scale_x_continuous(
    limits = x_limit                   
  ) +
  scale_y_continuous(
    limits = y_limit                   
  )

# Create MCA biplot with adjusted label positions using ggrepel
fviz_mca_biplot(
  mca_result,
  label = "none",                      # Remove default labels
  habillage = df_mca$Trophic.level,    # Color by Trophic.level
  addEllipses = TRUE,                  # Add ellipses for grouping
  ellipse.level = 0.95,                # Confidence level for ellipses
  palette = "Set1",                    # Set color palette
  geom.var = "none",                   # Disable variable points (red triangles)
  geom.ind = "point"                   # Keep individual points (colored circles)
) +
  theme_minimal() +                    # Minimal theme for clean plot
  theme(
    legend.position = "right",         # Adjust legend position
    axis.text.x = element_text(size = 14),  # X-axis text size
    axis.text.y = element_text(size = 14),  # Y-axis text size
    axis.title.x = element_text(size = 16), # X-axis title size
    axis.title.y = element_text(size = 16)  # Y-axis title size
  ) +
  # Add axis labels with explained variance percentages
  labs(
    x = paste0("Dim 1 (", round(percent_var[1], 1), "%)"),
    y = paste0("Dim 2 (", round(percent_var[2], 1), "%)")
  ) +
  # Use ggrepel to automatically adjust label positions and avoid overlap
  geom_text_repel(
    data = present_coords,             # Use filtered coordinates
    aes(x = label_x,                   # Adjusted x position for labels
        y = label_y,                   # Adjusted y position for labels
        label = polymer),              # Use polymer names as labels
    size = 4,                          # Label text size
    max.overlaps = 50,                  # Allow for more overlap before repositioning
    box.padding = 0.5,                  # Space between labels and points
    point.padding = 0.5,                # Space between label and points
    min.segment.length = 0.5,           # Minimum segment length for label connection lines
    segment.color = NA                 # Remove the connecting lines (black lines)
  ) +
  # Set consistent axis limits
  scale_x_continuous(
    limits = c(floor(min(c(mca_result$ind$coord[, 1], mca_result$ind$coord[, 2]))) - 1, 
               ceiling(max(c(mca_result$ind$coord[, 1], mca_result$ind$coord[, 2]))) + 1)
  ) +
  scale_y_continuous(
    limits = c(floor(min(c(mca_result$ind$coord[, 1], mca_result$ind$coord[, 2]))) - 1, 
               ceiling(max(c(mca_result$ind$coord[, 1], mca_result$ind$coord[, 2]))) + 1)
  )

