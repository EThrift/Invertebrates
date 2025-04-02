#created on 21/08/2024
#created by Emily Thrift


# Load necessary libraries
library(dplyr)
library(lme4)
library(tidyverse)
library(ggplot2)
library(glmmTMB)

df_size <- read_csv("Size.csv")

# Remove spaces from column names
colnames(df_size) <- make.names(colnames(df_size))


#model relationship between taxonomic group and size

df_size$Biomass_log <- log(df_size$Biomass + 1)  # Offset

library(lmerTest)

model1 <- lmer(Size ~ Land.cover + (1 | Site) + offset(Biomass_log), data = df_size, )
model2 <- lmer(Size ~ Land.cover + Taxonomic.group + (1 | Site) + offset(Biomass_log), data = df_size)


# Summary for fixed effects with F-statistics
summary(model1)
summary(model2)

# Get F-statistics for fixed effects
anova(model1, model2) 
anova(model2)

# Simulate residuals
sim_res <- simulateResiduals(fittedModel = model2)

# Plot diagnostics
plot(sim_res)
#QQplot shows good fit

# Overdispersion test using DHARMa
testDispersion(sim_res)



# Calculate the estimated marginal means for Taxonomic.group
emmeans_tax_size <- emmeans(model2, ~ Taxonomic.group)

# Perform pairwise comparisons
pairwise_comparisons <- contrast(emmeans_tax_size, method = "pairwise")

# View the results of the pairwise comparisons
summary(pairwise_comparisons)

library(ggplot2)


# Filter the dataset to exclude Size values over 2
df_size_filtered <- df_size %>%
  filter(Size <= 2)

# Reorder the Trophic.level factor
df_size_filtered <- df_size %>%
  filter(Size <= 2) %>%
  mutate(Trophic.level = factor(Trophic.level, levels = c("Herbivore", "Omnivore", "Carnivore", "Detritivore")))


# Calculate mean and standard deviation for the filtered data
df_size_summary <- df_size_filtered %>%
  group_by(Trophic.level) %>%
  summarise(size_mean = mean(Size),
            size_SD = sd(Size))


# Create the violin plot with mean as a dot and standard deviation as error bars
ggplot(df_size_filtered, aes(x = Trophic.level, y = Size)) +
  geom_violin(fill = alpha("grey80", 0.4), color = "black") +  # Violin plot with lighter fill
  geom_point(data = df_size_summary, aes(x = Trophic.level, y = size_mean), 
             color = "blue", size = 3) +                        # Blue dot for the mean
  geom_errorbar(data = df_size_summary, aes(x = Trophic.level, y = size_mean, 
                                            ymin = size_mean - size_SD, 
                                            ymax = size_mean + size_SD),
                width = 0.2, color = "blue") +                 # Blue error bars for SD
  theme_minimal() +
  labs(x = "Trophic level", y = "Size (mm)") +
  theme(
    axis.text.x = element_text(size = 14),       # Larger text for x-axis
    axis.text.y = element_text(size = 14),       # Larger text for y-axis
    axis.title.x = element_text(size = 16),      # Larger title for x-axis
    axis.title.y = element_text(size = 16),      # Larger title for y-axis
    panel.grid = element_blank(),                # Removes default grid lines
    axis.line = element_line(color = "black")    # Adds black axis lines
  ) +
  coord_cartesian(ylim = c(0, 2))   


                           

#model reladf_size#model relationship between size and trophic level
model1 <- lmer(Size ~  Land.cover + (1 | Site) + offset(Biomass_log), data = df_size)
summary(model1)

model2 <- lmer(Size ~Trophic.level + Land.cover + (1 | Site) + offset(Biomass_log), data = df_size)
summary(model2)

anova(model1,model2)

anova(model2)
#model 2 signifncantly better model chi = 31.4p = <0.0001 trophic level significantly impacts size 


# Simulate residuals
sim_res <- simulateResiduals(fittedModel = model2)

# Plot diagnostics
plot(sim_res)

#shows good qq and residual plot

# Overdispersion test using DHARMa
testDispersion(sim_res)
#not overdispersed p = 0.8

#pairwise comparisons 
# Calculate the estimated marginal means for Taxonomic.group
emmeans_tro_size <- emmeans(model2, ~ Trophic.level)

# Perform pairwise comparisons
pairwise_comparisons <- contrast(emmeans_tro_size, method = "pairwise")

# View the results of the pairwise comparisons
summary(pairwise_comparisons)

#carnivore - detritivore = 0.0001
#carnivore - herbivore = <0.0001
#herbivore &  detritivore significantly smaller than carnivore



# Filter the dataset to exclude Size values over 2
df_size_filtered <- df_size %>%
  filter(Size <= 2)

# Create a plot Trophic level

# Calculate mean and standard deviation
df_size_summary <- df_size_filtered %>%
  group_by(Taxonomic.group) %>%
  summarise(size_mean = mean(Size),
            size_SD = sd(Size))

# Reorder the Taxonomic.group variable
df_size_filtered$Taxonomic.group <- factor(df_size_filtered$Taxonomic.group, 
                                  levels = c("Coleoptera", "Isopoda", "Stylommatophora", 
                                             "Opisthopora", "Hemiptera", "Lepidoptera"))


# Create the violin plot with mean as a dot and standard deviation as error bars

ggplot(df_size_filtered, aes(x = Taxonomic.group, y = Size)) +
  geom_violin(fill = alpha("grey80", 0.4), color = "black") +  # Violin plot with lighter fill
  geom_point(data = df_size_summary, aes(x = Taxonomic.group, y = size_mean), 
             color = "blue", size = 3) +                        # Blue dot for the mean
  geom_errorbar(data = df_size_summary, aes(x = Taxonomic.group, y = size_mean, 
                                            ymin = size_mean - size_SD, 
                                            ymax = size_mean + size_SD),
                width = 0.2, color = "blue") +                 # Blue error bars for SD
  theme_minimal() +
  labs(x = "Taxonomic group", y = "Size (mm)") +
  theme(
    axis.text.x = element_text(size = 14),       # Larger text for x-axis
    axis.text.y = element_text(size = 14),       # Larger text for y-axis
    axis.title.x = element_text(size = 16),      # Larger title for x-axis
    axis.title.y = element_text(size = 16),      # Larger title for y-axis
    panel.grid = element_blank(),                # Removes default grid lines
    axis.line = element_line(color = "black")    # Adds black axis lines
  ) +
  coord_cartesian(ylim = c(0, 2))                # Flexible y-axis limits to include SD

#for land cover

# Calculate mean and standard deviation
df_size_summary <- df_size_filtered %>%
  group_by(Land.cover) %>%
  summarise(size_mean = mean(Size),
            size_SD = sd(Size))

# Reorder the Taxonomic.group variable
df_size_filtered$Land.cover <- factor(df_size_filtered$Land.cover)


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
#model relationship between size and overall biomass
