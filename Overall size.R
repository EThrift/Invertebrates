#created on 21/08/2024
#created by Emily Thrift


# Load necessary libraries
library(dplyr)
library(lme4)
library(tidyverse)
library(ggplot2)
library(glmmTMB)

df_size <- read_csv("Size.csv")

df_biomass <- read_csv("biomass.csv")

df_indi <- read_csv("Plastics_per_individual.csv")

# Remove spaces from column names
colnames(df_size) <- make.names(colnames(df_size))

# Remove spaces from column names
colnames(df_biomass) <- make.names(colnames(df_biomass))

# Remove spaces from column names
colnames(df_indi) <- make.names(colnames(df_indi))
#model relationship between taxonomic group and size

model1 <- glmmTMB(Size ~  Land.cover + (1 | Site),  data = df_size)
summary(model1)

model2 <- glmmTMB(Size ~ Land.cover + Taxonomic.group + (1 | Site),  data = df_size)
summary(model2)

anova(model1,model2)

#model 2 signficantly better model p = <0.001

# Load the DHARMa package
library(DHARMa)

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

# Create a box plot of Size by Taxonomic.group

library(ggplot2)

#have a function which means you can just call this when you want to calculate standard error 
std_err <- function(x) {
  SE <- sd(x)/sqrt(length(x))
  return(SE)
}

# Calculate mean and standard deviation
df_size_summary <- df_size %>%
  group_by(Trophic.level) %>%
  summarise(size_mean = mean(Size),
            size_SD = sd(Size))


# Create the violin plot with mean as a dot and standard deviation as error bars
ggplot(df_size, aes(x = Trophic.level, y = Size)) +
  geom_violin(fill = alpha("grey80", 0.4), color = "black") +  # Violin plot with lighter fill
  geom_point(data = df_size_summary, aes(y = size_mean), 
             color = "blue", size = 3) +                        # Red dot for the mean
  geom_errorbar(data = df_size_summary, aes(y = size_mean, 
                                            ymin = size_mean - size_SD, 
                                            ymax = size_mean + size_SD),
                width = 0.2, color = "blue") +                 # Red error bars for SD
  theme_minimal() +
  labs(x = "Trophic.level", y = "Size (mm)") +
  theme(
    axis.text.x = element_text(size = 14),  # Larger text for x-axis
    axis.text.y = element_text(size = 14),  # Larger text for y-axis
    axis.title.x = element_text(size = 16),  # Larger title for x-axis
    axis.title.y = element_text(size = 16)   # Larger title for y-axis
  ) +
  ylim(0, 6)  # Extend y-axis to a maximum of 6  # Extend y-axis to a maximum of 6
           

#model relationship between size and trophic level
model1 <- glmmTMB(Size ~  Land.cover + (1 | Site), data = df_size)
summary(model1)

model2 <- glmmTMB(Size ~Trophic.level + Land.cover + (1 | Site), data = df_size)
summary(model2)

anova(model1,model2)

#model 2 signifncantly better model p = <0.0001


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


# Create a plot Trophic level

# Calculate mean and standard deviation
df_size_summary <- df_size %>%
  group_by(Taxonomic.group) %>%
  summarise(size_mean = mean(Size),
            size_SD = sd(Size))

# Reorder the Taxonomic.group variable
df_size$Taxonomic.group <- factor(df_size$Taxonomic.group, 
                                  levels = c("Coleoptera", "Isopoda", "Stylommatophora", 
                                             "Opisthopora", "Hemiptera", "Lepidoptera"))


# Create the violin plot with mean as a dot and standard deviation as error bars
ggplot(df_size, aes(x = Taxonomic.group, y = Size)) +
  geom_violin(fill = alpha("grey80", 0.4), color = "black") +  # Violin plot with lighter fill
  geom_point(data = df_size_summary, aes(y = size_mean), 
             color = "blue", size = 3) +                        # Red dot for the mean
  geom_errorbar(data = df_size_summary, aes(y = size_mean, 
                                            ymin = size_mean - size_SD, 
                                            ymax = size_mean + size_SD),
                width = 0.2, color = "blue") +                 # Red error bars for SD
  theme_minimal() +
  labs(x = "Taxonomic group", y = "Size (mm)") +
  theme(
    axis.text.x = element_text(size = 14),  # Larger text for x-axis
    axis.text.y = element_text(size = 14),  # Larger text for y-axis
    axis.title.x = element_text(size = 16),  # Larger title for x-axis
    axis.title.y = element_text(size = 16)   # Larger title for y-axis
  ) +
  ylim(0, 6)  # Extend y-axis to a maximum of 6  # Extend y-axis to a maximum of 6

#model relationship between size and overall biomass

# Examine distributions of the log-transformed number of plastic pieces
df_size %>% 
  ggplot(aes(x = Size)) +  # Close aes() before geom_histogram
  geom_histogram(bins = 15, fill = "blue", alpha = 0.7) +  # Add fill and transparency for aesthetics
  labs(title = "Histogram of Log10 of Number of Plastic Pieces", 
       x = "Log10(Number of Plastic Pieces)", 
       y = "Frequency") +
  theme_minimal()  # Use a clean theme for better visibility

# 1. Plotting the Histogram of Biomass
df_size %>% 
  ggplot(aes(x = Biomass)) + 
  geom_histogram(bins = 30, fill = "blue", alpha = 0.7) +
  labs(title = "Histogram of Biomass (g)", 
       x = "Biomass (g)", 
       y = "Frequency") +
  theme_minimal()


# 3. No transformation applied; keep original scale
df_size <- df_size %>%
  mutate(Biomass = Biomass,  # Retain original scale
         Size = Size) # Retain original scale

# 4. Fit models without log transformation
MOD.1 <- glmmTMB(Plastic.per.gram ~ Biomass, data = df_biomass)

summary(MOD.1)
MOD.2 <- glmmTMB(Plastic.per.gram ~ Biomass + Taxonomic.group, data = df_biomass)

summary(MOD.2)

# 5. Compare models
anova(MOD.1, MOD.2, test = "Chisq")

library(emmeans)
#pairwise comparisons 
# Calculate the estimated marginal means for Taxonomic.group
emmeans_tax_group <- emmeans(MOD.2, ~ Taxonomic.group)

# Perform pairwise comparisons
pairwise_comparisons <- contrast(emmeans_tax_group, method = "pairwise")

# View the results of the pairwise comparisons
summary(pairwise_comparisons)

# Calculate mean and standard deviation
df_size_summary <- df_biomass %>%
  group_by(Taxonomic.group) %>%
  summarise(size_mean = mean(Plastic.per.gram),
            size_SD = sd(Plastic.per.gram))


# Create the violin plot with mean as a dot and standard deviation as error bars
ggplot(df_biomass, aes(x = Taxonomic.group, y = Plastic.per.gram )) +
  geom_violin(fill = alpha("grey80", 0.4), color = "black") +  # Violin plot with lighter fill
  geom_point(data = df_size_summary, aes(y = size_mean), 
             color = "blue", size = 3) +                        # Red dot for the mean
  geom_errorbar(data = df_size_summary, aes(y = size_mean, 
                                            ymin = size_mean - size_SD, 
                                            ymax = size_mean + size_SD),
                width = 0.2, color = "blue") +                 # Red error bars for SD
  theme_minimal() +
  labs(x = "Taxonomic group", y = "Plastic per gram") +
  theme(
    axis.text.x = element_text(size = 14),  # Larger text for x-axis
    axis.text.y = element_text(size = 14),  # Larger text for y-axis
    axis.title.x = element_text(size = 16),  # Larger title for x-axis
    axis.title.y = element_text(size = 16)   # Larger title for y-axis
  ) +
  ylim(0, 200)  # Extend y-axis to a maximum of 6  # Extend y-axis to a maximum of 6



### Build and evaluate model ###
MOD.1 <- glmmTMB( Plastic.per.gram ~  Biomass ,
                  data = df_biomass)

summary(MOD.1)



### Build and evaluate model ###
MOD.2 <- glmmTMB(Plastic.per.gram ~ Trophic.level +  Biomass,
                 data = df_biomass)

summary(MOD.2)

anova(MOD.1,
      MOD.2,
      test="Chi")
#model 2 better p = 0.003

library(DHARMa)

# Simulate residuals
sim_res <- simulateResiduals(fittedModel = MOD.2)

# Plot diagnostics
plot(sim_res)

#shows good qq and residual plot

# Overdispersion test using DHARMa
testDispersion(sim_res)
#not overdispersed p = 0.8

#pairwise comparisons 
# Calculate the estimated marginal means for Taxonomic.group
emmeans_tro_size <- emmeans(MOD.2, ~ Trophic.level)

# Perform pairwise comparisons
pairwise_comparisons <- contrast(emmeans_tro_size, method = "pairwise")

# View the results of the pairwise comparisons
summary(pairwise_comparisons)


df_biomass_filtered <- df_biomass %>%
  filter(Trophic.level != "omnivore")
# Calculate mean and standard deviation
df_size_summary <- df_biomass_filtered %>%
  group_by(Trophic.level) %>%
  summarise(size_mean = mean(Plastic.per.gram),
            size_SD = sd(Plastic.per.gram))



# Create the violin plot with mean as a dot and standard deviation as error bars
ggplot(df_biomass_filtered, aes(x = Trophic.level, y = Plastic.per.gram )) +
  geom_violin(fill = alpha("grey80", 0.4), color = "black") +  # Violin plot with lighter fill
  geom_point(data = df_size_summary, aes(y = size_mean), 
             color = "blue", size = 3) +                        # Red dot for the mean
  geom_errorbar(data = df_size_summary, aes(y = size_mean, 
                                            ymin = size_mean - size_SD, 
                                            ymax = size_mean + size_SD),
                width = 0.2, color = "blue") +                 # Red error bars for SD
  theme_minimal() +
  labs(x = "Trophic level", y = "Plastic per gram") +
  theme(
    axis.text.x = element_text(size = 14),  # Larger text for x-axis
    axis.text.y = element_text(size = 14),  # Larger text for y-axis
    axis.title.x = element_text(size = 16),  # Larger title for x-axis
    axis.title.y = element_text(size = 16)   # Larger title for y-axis
  ) +
  ylim(0, 200)  

# Extract unique Taxonomic groups from df_biomass
unique_taxonomic_groups <- unique(df_biomass$Taxonomic.group)

# Step 1: Create a new column for the adjusted biomass (standardized to 1g)
df_biomass$Standardized.Biomass <- 1  # All set to 1g

unique(df_biomass$Trophic.level)


# Step 2: Adjust the number of pieces of plastic
# Calculate the adjustment factor for the number of pieces
df_biomass$Adjusted.Plastic <- df_biomass$Number.of.pieces.of.plastics * (df_biomass$Standardized.Biomass / df_biomass$Biomass)

df_size_summary <- df_biomass %>%
  group_by(Trophic.level) %>%
  summarise(size_mean = mean(Adjusted.Plastic),
            size_SE = std_err(Adjusted.Plastic))




# Create a jitter plot to visualize the relationship between taxonomic group and adjusted plastic count
ggplot(df_biomass, aes(x = Trophic.level, y = Adjusted.Plastic, color = Trophic.level)) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.5) +
  labs(x = "Taxonomic Group", y = "Adjusted Number of Plastic Pieces (Standardized to 1g Biomass)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(df_biomass, aes(x = Trophic.level, y = Adjusted.Plastic )) +
  geom_jitter(width = 0.15, height = 0.1, alpha = 0.4) +
  geom_pointrange(data = df_size_summary,
                  aes(y = size_mean,
                      ymin = size_mean - size_SE,
                      ymax = size_mean + size_SE),
                  colour = "blue") +
  theme_minimal() + 
  labs(x = "Taxonomic group", y = "Plastic per gram of biomass") +
  geom_point(data = df_size_summary,
             aes(y = size_mean),
             size = 3, colour = "blue")

# Similar plot with geom_violin
ggplot(df_biomass, aes(x = Trophic.level, y = Adjusted.Plastic)) +
  geom_violin(fill = alpha("black", 0.2)) +
  geom_jitter(width = 0.05, height = 0, alpha = 0.4) +
  theme_bw() +
  geom_point(data = df_size_summary, aes(y = size_mean), size = 3, colour = "blue") +
  labs(x = "Trophic level", y = "Plastics per gram of biomass") +
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16)
  ) +
  scale_y_continuous(limits = c(0, 170), breaks = seq(0, 170, by = 20))

### Build and evaluate model ###
MOD.1 <- glmmTMB( Plastic.per.gram ~  Standardized.Biomass  ,
                  data = df_biomass)

summary(MOD.1)

### Build and evaluate model ###
MOD.2 <- glmmTMB(Plastic.per.gram ~ Trophic.level +  Standardized.Biomass,
                 data = df_biomass)

summary(MOD.2)

anova(MOD.1,
      MOD.2,
      test="Chi")

#model 2 p = 0.027

library(DHARMa)

# Simulate residuals
sim_res <- simulateResiduals(fittedModel = MOD.2)

# Plot diagnostics
plot(sim_res)

#shows good qq and residual plot

# Overdispersion test using DHARMa
testDispersion(sim_res)
#not overdispersed p = 0.8

#pairwise comparisons 
# Calculate the estimated marginal means for Taxonomic.group
emmeans_tro_size <- emmeans(MOD.2, ~ Trophic.level)

# Perform pairwise comparisons
pairwise_comparisons <- contrast(emmeans_tro_size, method = "pairwise")

# View the results of the pairwise comparisons
summary(pairwise_comparisons)

#Carnivore and herbivore p = 0.016

####per individual

library(glmmTMB)
# 4. Fit models without log transformation
MOD.1 <- glmmTMB(Plastic.per.individual ~ Land.cover, data = df_indi)

summary(MOD.1)
MOD.2 <- glmmTMB(Plastic.per.individual ~ Taxonomic.group + Land.cover, data = df_indi)
summary(MOD.2)

# 5. Compare models
anova(MOD.1, MOD.2, test = "Chisq")

library(emmeans)
#pairwise comparisons 
# Calculate the estimated marginal means for Taxonomic.group
emmeans_tax_group <- emmeans(MOD.2, ~ Taxonomic.group)

# Perform pairwise comparisons
pairwise_comparisons <- contrast(emmeans_tax_group, method = "pairwise")

# View the results of the pairwise comparisons
summary(pairwise_comparisons)

# Calculate mean and standard deviation
df_size_summary <- df_indi %>%
  group_by(Taxonomic.group) %>%
  summarise(size_mean = mean(Plastic.per.individual),
            size_SD = sd(Plastic.per.individual))


# Create the violin plot with mean as a dot and standard deviation as error bars
ggplot(df_indi, aes(x = Taxonomic.group, y = Plastic.per.individual )) +
  geom_violin(fill = alpha("grey80", 0.4), color = "black") +  # Violin plot with lighter fill
  geom_point(data = df_size_summary, aes(y = size_mean), 
             color = "blue", size = 3) +                        # Red dot for the mean
  geom_errorbar(data = df_size_summary, aes(y = size_mean, 
                                            ymin = size_mean - size_SD, 
                                            ymax = size_mean + size_SD),
                width = 0.2, color = "blue") +                 # Red error bars for SD
  theme_minimal() +
  labs(x = "Taxonomic group", y = "Plastic per individual") +
  theme(
    axis.text.x = element_text(size = 14),  # Larger text for x-axis
    axis.text.y = element_text(size = 14),  # Larger text for y-axis
    axis.title.x = element_text(size = 16),  # Larger title for x-axis
    axis.title.y = element_text(size = 16)   # Larger title for y-axis
  ) +
  ylim(0, 5)  # Extend y-axis to a maximum of 6  # Extend y-axis to a maximum of 6


# 4. Fit models without log transformation
MOD.1 <- glmmTMB(Plastic.per.individual ~ Land.cover, data = df_indi)

summary(MOD.1)
MOD.2 <- glmmTMB(Plastic.per.individual ~ Trophic.level + Land.cover, data = df_indi)
summary(MOD.2)

# 5. Compare models
anova(MOD.1, MOD.2, test = "Chisq")

library(emmeans)
#pairwise comparisons 
# Calculate the estimated marginal means for Taxonomic.group
emmeans_tax_group <- emmeans(MOD.2, ~ Trophic.level)

# Perform pairwise comparisons
pairwise_comparisons <- contrast(emmeans_tax_group, method = "pairwise")

# View the results of the pairwise comparisons
summary(pairwise_comparisons)

# Calculate mean and standard deviation
df_size_summary <- df_indi %>%
  group_by(Trophic.level) %>%
  summarise(size_mean = mean(Plastic.per.individual),
            size_SD = sd(Plastic.per.individual))


# Create the violin plot with mean as a dot and standard deviation as error bars
ggplot(df_indi, aes(x = Trophic.level, y = Plastic.per.individual )) +
  geom_violin(fill = alpha("grey80", 0.4), color = "black") +  # Violin plot with lighter fill
  geom_point(data = df_size_summary, aes(y = size_mean), 
             color = "blue", size = 3) +                        # Red dot for the mean
  geom_errorbar(data = df_size_summary, aes(y = size_mean, 
                                            ymin = size_mean - size_SD, 
                                            ymax = size_mean + size_SD),
                width = 0.2, color = "blue") +                 # Red error bars for SD
  theme_minimal() +
  labs(x = "Trophic level", y = "Plastic per individual") +
  theme(
    axis.text.x = element_text(size = 14),  # Larger text for x-axis
    axis.text.y = element_text(size = 14),  # Larger text for y-axis
    axis.title.x = element_text(size = 16),  # Larger title for x-axis
    axis.title.y = element_text(size = 16)   # Larger title for y-axis
  ) +
  ylim(0, 5)  # Extend y-axis to a maximum of 6  # Extend y-axis to a maximum of 6
