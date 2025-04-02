 #created on 22/08/24
# By Emily Thrift
# script for a glm to assess the levels of polymer levels in both soil samples from different land cover types

#load libraries
library(tidyverse)
library(ggplot2)
library(lme4)
library(glmmTMB)
library(emmeans)


#load dataframe 

df_polymer <- read.csv("Soil polymer.csv")


# Remove spaces from column names
colnames(df_polymer) <- make.names(colnames(df_polymer))

unique(df_polymer$Land.cover)  # Check unique values in the 'Land.cover' column
unique(df_polymer$Nylon)       # Check unique values in a polymer column


df_plot <- read.csv ("Soil polymer plot_talk.csv")
colnames(df_plot) <- make.names(colnames(df_plot))


df_plot

library(binom)

wilson_ci2 <- binom.confint(x = df_plot$n, n = df_plot$Total, method = "wilson")

# Add confidence intervals to the data using mutate
df_plot <- df_plot %>%
  mutate(lower = wilson_ci2$lower,
         upper = wilson_ci2$upper)
#
# Calculate label positions just above error bars
df_plot$label_y <- df_plot$upper * 100 + 3.5

# Summarize the total for each taxonomic group and create the label with "N = "
totals2 <- df_plot %>% group_by(Land.cover) %>% summarize(Total = unique(Total))
totals2$TotalLabel <- paste("N =", totals2$Total)

totals2

# Convert Polymer to a factor
df_plot$Polymer <- factor(df_plot$Polymer)

color_blind_friendly_palette <- c(
  "Other" = "#1f77b4",        # Blue
  "Polyester" = "#ff7f0e",    # Orange
  "Polyethylene" = "#2ca02c", # Green
  "Cellophane" = "#e377c2",   # Pink
  "Nylon" = "yellow"          # Yellow
)

ggplot(df_plot, aes(x = Polymer, y = Percentage, fill = Polymer)) +  
  geom_col(position = position_dodge(width = 0.8), color = "black") +  # Bar chart
  geom_errorbar(aes(ymin = lower * 100, ymax = upper * 100), 
                position = position_dodge(width = 0.8), width = 0.25) +  # Error bars
  geom_text(aes(y = label_y + 5, label = paste("n =", n)),  
            position = position_dodge(width = 0.8), size = 3.5, color = "black") +  # Label count above bars
  labs(x = "Polymer Type in soil samples", y = "% of Polymer") +  
  scale_fill_manual(values = color_blind_friendly_palette) +  # Use custom color palette
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 12, angle = 90, hjust = 1),  
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",  # Remove the legend
    plot.margin = margin(1, 1, 1, 1, "cm"),  
    strip.text = element_text(size = 12, color = "black")  
  ) +
  facet_wrap(~ Land.cover, scales = "free_x", ncol = 6)  # Facet by Land.cover with free x scales


# Reshape data into long format for easier analysis
df_long <- df_polymer %>%
  pivot_longer(cols = 5:13, names_to = "Polymer", values_to = "Presence") %>%
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
test_statistic <- kruskal_test$statistic
degrees_of_freedom <- kruskal_test$parameter
p_value <- kruskal_test$p.value

# Print the values
cat("Kruskal-Wallis chi-squared =", test_statistic, ", df =", degrees_of_freedom, ", p-value =", p_value, "\n")

#result has a p value of 0.36

# Ensure Land cover is a factor
diversity_indices_simpson$Land.cover <- as.factor(diversity_indices_simpson$Land.cover)


#no significant differences seen in the pairwise analysis

# Plotting Simpson Diversity Index across taxonomic groups
ggplot(diversity_indices_simpson, aes(x = Land.cover, y = Simpson, fill = Land.cover)) +
  geom_bar(stat = "identity", color = "black", show.legend = FALSE) +  # Bar plot
  labs(title = "Simpson Diversity Index of Polymers by Land cover",
       x = "Land cover", 
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

str(df_polymer)

# Summarize polymer data per land cover (wide format)

unique(df_polymer$Land.cover)  # Check unique values in the 'Land.cover' column
unique(df_polymer$Nylon)       # Check unique values in a polymer column

summary(df_polymer$Polyester)
summary(df_polymer$Polyethylene)


# Convert polymer columns to numeric
df_polymer <- df_polymer %>%
  mutate(across(Polyester:Polyethylene, as.numeric, na.rm = TRUE))

df_mca <- df_polymer %>%
  group_by(Land.cover) %>%
  summarise(across(Polyester:Polyethylene, sum, na.rm = TRUE))

df_mca



# Ensure land cover is a factor
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

