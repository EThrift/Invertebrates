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




df_land <- read.csv ("Land_cover_polymer.csv")
colnames(df_land) <- make.names(colnames(df_land))

df_land <- read.csv ("Land_cover_polymer_talk.csv")
colnames(df_land) <- make.names(colnames(df_land))

df_land

library(binom)

wilson_ci2 <- binom.confint(x = df_land$n, n = df_land$Total, method = "wilson")

# Add confidence intervals to the data using mutate
df_land <- df_land %>%
  mutate(lower = wilson_ci2$lower,
         upper = wilson_ci2$upper)
#
# Calculate label positions just above error bars
df_land$label_y <- df_land$upper * 100 + 3.5

# Summarize the total for each taxonomic group and create the label with "N = "
totals2 <- df_land %>% group_by(Land.cover) %>% summarize(Total = unique(Total))
totals2$TotalLabel <- paste("N =", totals2$Total)

totals2

# Convert Polymer to a factor
df_land$Polymer <- factor(df_land$Polymer)

# Create a summary of mean percentages for ordering
mean_order <- df_land %>%
  group_by(Polymer) %>%
  summarize(mean_percentage = mean(Percentage, na.rm = TRUE)) %>%
  arrange(desc(mean_percentage)) %>%
  pull(Polymer)

ggplot(df_land, aes(x = Polymer, y = Percentage)) +  # Reorder points by Percentage in descending order
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
  facet_wrap(~ Land.cover, scales = "free_x", ncol = 6) 

#Calculate Wilson Confidence Intervals for each taxonomic group
wilson_ci2 <- binom.confint(x = df_land$n, n = df_land$Total, method = "wilson")

#Add the confidence intervals to the data
df_land <- df_land %>%
  mutate(lower = wilson_ci2$lower,
         upper = wilson_ci2$upper)

# Calculate label positions just above the error bars
df_land$label_y <- df_land$upper * 100 + 3.5

#Summarize the total for each taxonomic group and create a label with "N = "
totals2 <- df_land %>% 
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



# Ensure Polymer is a factor with "Other" as the first level in df_land
df_land$Polymer <- factor(df_land$Polymer, 
                                    levels = c("Other", "Polyester", "Polyethylene", "Polypropylene", 
                                               "Polystyrene", "Cellophane", "Epoxy resin", 
                                               "Nylon", "Aluminium silicate"))

# Plot the data with custom colors
ggplot(df_land, aes(x = Polymer, y = Percentage, fill = Polymer)) +  
  geom_col(position = position_dodge(width = 0.8), color = "black") +  # Bar chart
  geom_errorbar(aes(ymin = lower * 100, ymax = upper * 100), 
                position = position_dodge(width = 0.8), width = 0.25) +  # Error bars
  geom_text(aes(y = label_y + 5, label = paste(n)),  # Add "n = " before the n value
            position = position_dodge(width = 0.8), size = 2.5, color = "black") +  # Label count above bars
  labs(x = "Polymer Type in land cover types", y = "% of Polymer") +  # Axis labels
  scale_fill_manual(values = setNames(color_blind_friendly_palette, levels(df_land$Polymer))) +  # Assign custom colors
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
  facet_wrap(~ Land.cover, scales = "free_x", ncol = 6)  # Facet by Taxonomic.group with free x scales



colnames(df_land)


# Reshape data into long format for easier analysis
df_long <- df_polymer %>%
  pivot_longer(cols = 13:24, names_to = "Polymer", values_to = "Presence") %>%
  filter(Presence == 1)  # Keep only presence data

# Step 1: Calculate Simpson Diversity Index by taxonomic group
diversity_indices_simpson <- df_long %>%
  group_by(Land.cover, Polymer) %>%
  summarise(Presence = sum(Presence), .groups = 'drop') %>%
  group_by(Land.cover) %>%
  summarise(
    Simpson = 1 - sum((Presence / sum(Presence))^2),
    Richness = n()
  )

kruskal_test <- kruskal.test(Simpson ~ Land.cover, data = diversity_indices_simpson)

# Extract the  test statistic, degrees of freedom, and p-value
test_statistic <- kruskal_test_result$statistic
degrees_of_freedom <- kruskal_test_result$parameter
p_value <- kruskal_test_result$p.value

# Print the values
cat("Kruskal-Wallis chi-squared =", test_statistic, ", df =", degrees_of_freedom, ", p-value =", p_value, "\n")

#result has a p value of 0.41

# Ensure Taxonomic.group is a factor
diversity_indices_simpson$Land.cover <- as.factor(diversity_indices_simpson$Land.cover)


#no significant differences seen in the pairwise analysis

# Plotting Simpson Diversity Index across taxonomic groups
ggplot(diversity_indices_simpson, aes(x = Land.cover, y = Simpson, fill = Land.cover)) +
  geom_bar(stat = "identity", color = "black", show.legend = FALSE) +  # Bar plot
  labs(title = "Simpson Diversity Index of Polymers by Taxonomic Group",
       x = "Taxonomic Group", 
       y = "Simpson Diversity Index") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
# Check if Taxonomic.group is a factor
diversity_indices_simpson$Land.cover <- as.factor(diversity_indices_simpson$Land.cover)

# Run the Kruskal-Wallis test
kruskal_test_result <- kruskal.test(Simpson ~ Land.cover, data = diversity_indices_simpson)

# Print the Kruskal-Wallis test result
# Extract the test statistic, degrees of freedom, and p-value
test_statistic <- kruskal_test_result$statistic
degrees_of_freedom <- kruskal_test_result$parameter
p_value <- kruskal_test_result$p.value

# Print the values
cat("Kruskal-Wallis chi-squared =", test_statistic, ", df =", degrees_of_freedom, ", p-value =", p_value, "\n")


library(FSA)

# Now run Dunn's test again
dunn_test_results <- dunnTest(Simpson ~ Land.cover, data = diversity_indices_simpson, method = "bonferroni")

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
  group_by(Land.cover) %>%
  summarise(across(12:24, sum))  # Sum the presence/absence for each polymer

df_mca


# Ensure Taxonomic.group is a factor
df_mca$Land.cover <- as.factor(df_mca$Land.cover)

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

# Adjust axes to range from 2, 1, 0, -1, -2 based on your request
x_limit <- c(-2, 2)  # Adjusting to the range you specified (from -2 to 2)
y_limit <- c(-2, 2)  # Similarly adjusting the y-axis range

# Continue with your plot generation code...
fviz_mca_biplot(
  mca_result,
  label = "none",                      
  habillage = df_mca$Land.cover,    
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
library(ggalluvial)
library(tidyverse)
library(igraph)
library(ggraph)

df_long <- read.csv("table comparisions polymer short.csv")
# Ensure column names are correct

# Create edges for the network (showing how Polyester moves through trophic levels)
edges <- df_long %>%
  select(Trophic_level, Target) %>%
  distinct()

# Create graph object
graph <- graph_from_data_frame(edges, directed = TRUE)

# Improved network plot
ggraph(graph, layout = "tree") +  # Try "tree" or "stress"
  geom_edge_link(arrow = arrow(length = unit(6, 'mm')), 
                 end_cap = circle(4, 'mm'), alpha = 0.6) +
  geom_node_point(size = 6, color = "steelblue") +
  geom_node_text(aes(label = name), repel = TRUE, size = 6) +
  theme_void() +
  ggtitle("Polyester Propagation Through Trophic Levels")

# Convert Sites to numeric
df_long <- df_long %>%
  mutate(Value = as.numeric(Value))


# Check if 'Soil' is being read as NA
df_long <- df_long %>%
  mutate(Trophic_level = as.factor(Trophic_level),
         Target = as.factor(Target))

# Define color mapping (fix NA issue)
color_mapping <- c("Soil" = "#E69F00", 
                   "Detritivore" = "#56B4E9", 
                   "Herbivore" = "#F0E442", 
                   "Carnivore" = "#00A087")

# Ensure Soil is correctly interpreted and avoid NA values
df_long <- df_long %>%
  mutate(Trophic_level = factor(Trophic_level, levels = names(color_mapping)),
         Target = factor(Target, levels = names(color_mapping)))

# Plot the Sankey Diagram with x-axis removed & NO legend
ggplot(df_long, aes(axis1 = Trophic_level, axis2 = Target, y = Value)) +
  geom_alluvium(aes(fill = Trophic_level), alpha = 0.8) +  
  geom_stratum(aes(fill = Trophic_level), color = "black") +  
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), color = "black") +  
  scale_fill_manual(values = color_mapping, na.value = "#E69F00") +  
  theme_minimal() +
  theme(legend.position = "none",  # Completely remove the legend
        axis.text.x = element_blank(),  
        axis.ticks.x = element_blank(),  
        axis.title.x = element_blank(),  
        axis.text.y = element_blank(),  
        axis.ticks.y = element_blank())
library(RColorBrewer)

# Load the dataset
df <- read.csv("table comparisions polymer short.csv")

colnames(df)

# Define a 12-color colorblind-friendly palette
colorblind_palette <- c("#E69F00", "#56B4E9", "#009E73")



# Assign colors dynamically based on unique polymers
unique_polymers <- unique(df_long$Polymer)
polymer_colors <- setNames(colorblind_palette[1:length(unique_polymers)], unique_polymers)

# Convert factors to ensure proper ordering
df_long <- df_long %>%
  mutate(Trophic_level = factor(Trophic_level, levels = c("Soil", "Detritivore", "Herbivore", "Carnivore")),
         Target = factor(Target, levels = c("Soil", "Detritivore", "Herbivore", "Carnivore")),
         Polymer = factor(Polymer, levels = unique_polymers))  # Ensure correct order

# Plot the Sankey Diagram
ggplot(df_long, aes(axis1 = Trophic_level, axis2 = Target, y = Value)) +
  geom_alluvium(aes(fill = Polymer), alpha = 0.8, width = 1/12) +  
  geom_stratum(aes(fill = Trophic_level), color = "black") +  
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), color = "black") +  
  scale_fill_manual(values = polymer_colors) +
  theme_minimal() +
  labs(title = "Polymer Distribution Across Trophic Levels",
       y = "Number of Sites",
       fill = "Legend") +  # Adding the legend title here
  theme(legend.position = "right",  # Positioning the legend to the right
        axis.text.x = element_blank(),  
        axis.ticks.x = element_blank(),  
        axis.title.x = element_blank(),  
        axis.text.y = element_blank(),  
        axis.ticks.y = element_blank())


library(circlize)
library(dplyr)
library(RColorBrewer)  # Load RColorBrewer package for color palettes

# Load the data
data <- read.csv("table comparisions polymer short.csv")

# Create a matrix for the chord diagram
edges <- data %>%
  group_by(Polymer, Trophic_level, Target) %>%
  summarise(Value = sum(Value), .groups = 'drop')

# Adjust the Value column to only contain 0, 1, or 2
edges$Value <- ifelse(edges$Value > 2, 2, edges$Value)  # Set anything over 2 to 2
edges$Value <- round(edges$Value, 0)  # Ensure values are rounded to integers

# Define sector groups (Polymers)
polymers <- unique(edges$Polymer)

# Set up the color palette for different polymers using a color-blind friendly palette
poly_colors <- setNames(brewer.pal(min(length(polymers), 12), "Set1"), polymers)

# Prepare the links with Polymer as a color indicator for the chords
edges$col <- poly_colors[edges$Polymer]

# Generate the chord diagram
chordDiagram(
  x = edges[, c("Trophic_level", "Target", "Value")], 
  grid.col = poly_colors,   # Color the sectors by polymer
  col = edges$col,          # Color the chords by polymer
  transparency = 0.5,
  annotationTrack = c("name", "grid")  # Add name to grid to improve clarity
)

# Add a title
title("Polymer Influence on Trophic Level Interactions")

# Add a legend for polymers
legend("topright", 
       legend = names(poly_colors), 
       fill = poly_colors, 
       title = "Polymers", 
       cex = 0.8, 
       bty = "n", 
       border = "black")


df_polymer <- read.csv("Polymer_all.csv")


# Remove spaces from column names
colnames(df_polymer) <- make.names(colnames(df_polymer))
# Reshape data into long format for easier analysis
df_long <- df_polymer %>%
  pivot_longer(cols = 13:27, names_to = "Polymer", values_to = "Presence") %>%
  filter(Presence == 1)  # Keep only presence data

# Step 1: Calculate Simpson Diversity Index by taxonomic group
diversity_indices_simpson <- df_long %>%
  group_by(Land.cover, Polymer) %>%
  summarise(Presence = sum(Presence), .groups = 'drop') %>%
  group_by(Land.cover) %>%
  summarise(
    Simpson = 1 - sum((Presence / sum(Presence))^2),
    Richness = n()
  )

kruskal_test <- kruskal.test(Simpson ~ Land.cover, data = diversity_indices_simpson)

# Extract the test statistic, degrees of freedom, and p-value
test_statistic <- kruskal_test_result$statistic
degrees_of_freedom <- kruskal_test_result$parameter
p_value <- kruskal_test_result$p.value

# Print the values
cat("Kruskal-Wallis chi-squared =", test_statistic, ", df =", degrees_of_freedom, ", p-value =", p_value, "\n")

#result has a p value of 0.41

# Ensure Taxonomic.group is a factor
diversity_indices_simpson$Land.cover <- as.factor(diversity_indices_simpson$Land.cover)


#no significant differences seen in the pairwise analysis

# Plotting Simpson Diversity Index across taxonomic groups
ggplot(diversity_indices_simpson, aes(x = Land.cover, y = Simpson, fill = Land.cover)) +
  geom_bar(stat = "identity", color = "black", show.legend = FALSE) +  # Bar plot
  labs(title = "Simpson Diversity Index of Polymers by Taxonomic Group",
       x = "Taxonomic Group", 
       y = "Simpson Diversity Index") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
# Check if Taxonomic.group is a factor

