#created on 22/08/24
# By Emily Thrift
# script for a glm to  assess the levels of polymer levels in both the taxonomic groups and trophic levels of different invertebrate samples tested 

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

df_taxonomic_talk <- read.csv ("Taxonomic_group_talk.csv")
colnames(df_taxonomic_talk) <- make.names(colnames(df_taxonomic_talk))
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
    axis.text.x = element_text(size = 12, angle = 90, hjust = 1),  # Adjust angle and size for better visibility
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",  # Remove the legend
    plot.margin = margin(1, 1, 1, 1, "cm")  # Increase the bottom margin
  ) +
  facet_wrap(~ Taxonomic.group, scales = "free_x", ncol = 6) 


#Calculate Wilson Confidence Intervals for each taxonomic group
wilson_ci2 <- binom.confint(x = df_taxonomic_talk$n, n = df_taxonomic_talk$Total, method = "wilson")

#Add the confidence intervals to the data
df_taxonomic_talk <- df_taxonomic_talk %>%
  mutate(lower = wilson_ci2$lower,
         upper = wilson_ci2$upper)

# Calculate label positions just above the error bars
df_taxonomic_talk$label_y <- df_taxonomic_talk$upper * 100 + 3.5

#Summarize the total for each taxonomic group and create a label with "N = "
totals2 <- df_taxonomic_talk %>% 
  group_by(Taxonomic.group) %>%
  summarize(Total = unique(Total))
totals2$TotalLabel <- paste("N =", totals2$Total)

# Load RColorBrewer for color-blind friendly palettes
library(RColorBrewer)

# Define a color-blind friendly palette with 8 colors (adjusting for the number of Polymer types)
# Define a named color palette for each polymer
color_blind_friendly_palette <- c(
  "Other" = "#1f77b4",        # Blue
  "Polyester" = "#ff7f0e",    # Orange
  "Polyethylene" = "#2ca02c", # Green
  "Polypropylene" = "#17becf",# Cyan
  "Polystyrene" = "#9467bd",  # Purple
  "Cellophane" = "#e377c2",   # Pink
  "Epoxy resin" = "#7f7f7f",  # Gray
  "Nylon" = "yellow",        # Yellow
  "Aluminium silicate" = "#8d8d8d"  # Dark Gray
)



# Ensure Polymer is a factor with "Other" as the first level in df_taxonomic_talk
df_taxonomic_talk$Polymer <- factor(df_taxonomic_talk$Polymer, 
                                    levels = c("Other", "Polyester", "Polyethylene", "Polypropylene", 
                                               "Polystyrene", "Cellophane", "Epoxy resin", 
                                               "Nylon", "Aluminium silicate"))

# Plot the data with custom colors
ggplot(df_taxonomic_talk, aes(x = Polymer, y = Percentage, fill = Polymer)) +  
  geom_col(position = position_dodge(width = 0.8), color = "black") +  # Bar chart
  geom_errorbar(aes(ymin = lower * 100, ymax = upper * 100), 
                position = position_dodge(width = 0.8), width = 0.25) +  # Error bars
  geom_text(aes(y = label_y + 5, label = paste(n)),  # Add "n = " before the n value
            position = position_dodge(width = 0.8), size = 2.5, color = "black") +  # Label count above bars
  labs(x = "Polymer Type in Taxonomic Groups", y = "% of Polymer") +  # Axis labels
  scale_fill_manual(values = setNames(color_blind_friendly_palette, levels(df_taxonomic_talk$Polymer))) +  # Assign custom colors
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
  facet_wrap(~ Taxonomic.group, scales = "free_x", ncol = 6)  # Facet by Taxonomic.group with free x scales




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
df_mca <- df_polymer %>%
  group_by(Taxonomic.group) %>%
  summarise(across(12:24, sum))  # Sum the presence/absence for each polymer

df_mca


# Ensure Taxonomic.group is a factor
df_mca$Taxonomic.group <- as.factor(df_mca$Taxonomic.group)

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

# Create MCA biplot for Taxonomic Group
x_limit <- c(floor(min(c(mca_result$ind$coord[, 1], mca_result$ind$coord[, 2]))) - 1, 
             ceiling(max(c(mca_result$ind$coord[, 1], mca_result$ind$coord[, 2]))) + 1)

y_limit <- x_limit  # You can use the same limits for y-axis

# Save limits to file
saveRDS(list(x_limit = x_limit, y_limit = y_limit), "axis_limits.rds")


# Continue with your plot generation code...
fviz_mca_biplot(
  mca_result,
  label = "none",                      
  habillage = df_mca$Taxonomic.group,    
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
