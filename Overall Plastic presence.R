#created on 20/08/24
# By Emily Thrift
# script for a glm to assess the levels of plastic presence in both the taxonomic groups and trophic levels of different invertebrate samples tested 

#load libraries
library(tidyverse)
library(ggplot2)
library(lme4)
library(glmmTMB)
library(emmeans)

#add in working location

setwd('C:/Users/Emily Thrift/Desktop
/Inverts')
list.files()

# Load data
df_presence <- read.csv("Plastic_presence.csv")
df_presence_land <- read.csv("plastic_presence_soil.csv")
# Remove spaces from column names
colnames(df_presence) <- make.names(colnames(df_presence))
colnames(df_presence_land) <- make.names(colnames(df_presence_land))

df_presence <- df_presence[df_presence$Trophic.level != "Omnivore", ]

df_presence$Biomass_log <- log(df_presence$Biomass + 1)  # Offset

table(df_presence$Plastic.present, df_presence$Taxonomic.group)
table(df_presence$Plastic.present, df_presence$Trophic.level)

summary(df_presence$Biomass)
table(df_presence$Land.cover, df_presence$Plastic.present)
table(df_presence$Taxonomic.group, df_presence$Plastic.present)

model1 <- glmmTMB(
  Plastic.present ~ (1 | Site) + Land.cover +  offset(Biomass_log), 
  data = df_presence, 
  family = binomial (link = logit)
)
summary(model1)

model2 <- glmmTMB(
  Plastic.present ~  Taxonomic.group + Land.cover + (1 | Site) + offset(Biomass_log), 
  data = df_presence, 
  family = binomial (link = logit)
)
summary(model2)

anova(model1, model2)

#model 2 with offset chi = 18.49 p = 0.0023
check_collinearity(model2)
# Load the DHARMa package
library(DHARMa)

table(df_presence$Taxonomic.group, df_presence$Plastic.present)


# Simulate residuals
sim_res <- simulateResiduals(fittedModel = model2)

performance::check_collinearity(model2)
testZeroInflation(sim_res)
testDispersion(sim_res)
testOutliers(sim_res, type = "bootstrap")

# Plot diagnostics
plot(sim_res)
#QQplot shows good fit

# Overdispersion test using DHARMa
testDispersion(sim_res)
#no evidence of overdispersion p = 0.88


# Calculate the estimated marginal means for Taxonomic.group
emmeans_tax_size <- emmeans(model2, ~ Taxonomic.group)

# Perform pairwise comparisons
pairwise_comparisons <- contrast(emmeans_tax_size, method = "pairwise")

# View the results of the pairwise comparisons
summary(pairwise_comparisons)

#Isopoda - Opisthopora (p = 0.0282) → Isopoda has significantly more than worms

library(binom)


# Calculate Percentage and Total by Taxonomic Group
df_percentage_total_taxonomic <- df_presence %>%
  group_by(Taxonomic.group) %>%
  summarise(Count = sum(Plastic.present), Total = n(), .groups = "drop") %>%
  mutate(Percentage = (Count / Total) * 100)

wilson_ci2 <- binom.confint(x = df_percentage_total_taxonomic$Count, 
                            n = df_percentage_total_taxonomic$Total, 
                            method = "wilson")


# Add confidence intervals to the data frame
df_percentage_total_taxonomic <- df_percentage_total_taxonomic %>%
  mutate(lower = wilson_ci2$lower,
         upper = wilson_ci2$upper)

# Calculate label positions just above error bars
df_percentage_total_taxonomic <- df_percentage_total_taxonomic %>%
  mutate(label_y = upper * 100 + 3.5) 

# Base plot
library(ggplot2)

df_percentage_total_taxonomic$Taxonomic.group <- factor(df_percentage_total_taxonomic$Taxonomic.group, 
                                                        levels = c("Lepidoptera", "Isopoda", "Coleoptera", 
                                                                   "Stylommatophora", "Hemiptera", "Opisthopora"))


# Create the plot
ggplot(df_percentage_total_taxonomic, aes(x = Taxonomic.group, y = Percentage, fill = Taxonomic.group)) +
  geom_bar(stat = "identity", fill = "grey90") +
  geom_errorbar(aes(ymin = lower * 100, ymax = upper * 100), width = 0.4) +
  geom_text(aes(y = label_y, label = paste(round(Count), "/", Total)), 
            size = 4, position = position_dodge(width = 0.9), fontface = "plain") +
  labs(x = "Taxonomic group", y = "% Plastic Presence") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    panel.grid = element_blank(),           # Removes default grid lines
    axis.line = element_line(color = "black") # Adds black axis lines
  ) +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.1)))

df_presence <- df_presence[df_presence$Trophic.level != "Omnivore", ]

modeltro <- glmmTMB(
  Plastic.present ~  (1 | Site) + Land.cover +  offset(Biomass_log), 
  data = df_presence, 
  family = binomial (link = logit))

summary(modeltro)

modeltro2 <- glmmTMB(
  Plastic.present ~  Land.cover + Trophic.level +(1 | Site) + offset(Biomass_log), 
  data = df_presence, 
  family = binomial (link = logit)
)
summary(modeltro2)

anova(modeltro, modeltro2)
check_collinearity(modeltro)
#with offset chi = 3.75, p 0.15

# Simulate residuals
sim_res2 <- simulateResiduals(fittedModel = modeltro)

# Plot diagnostics
plot(sim_res2)
#QQplot shows good fit

# Overdispersion test using DHARMa
testDispersion(sim_res2)
testOutliers(sim_res2, type = "bootstrap")
#no issues

modeltro <- glmmTMB(
  Plastic.present ~  (1 | Site) + Trophic.level +  offset(Biomass_log), 
  data = df_presence, 
  family = binomial (link = logit))

summary(modeltro)

modeltro2 <- glmmTMB(
  Plastic.present ~  Land.cover + Trophic.level +(1 | Site) + offset(Biomass_log), 
  data = df_presence, 
  family = binomial (link = logit)
)
summary(modeltro2)

anova(modeltro, modeltro2)

#with offset chi = 4.77, p 0.18

# Simulate residuals
sim_res2 <- simulateResiduals(fittedModel = modeltro)

# Plot diagnostics
plot(sim_res2)
#QQplot shows good fit

# Overdispersion test using DHARMa
testDispersion(sim_res2)
testOutliers(sim_res2, type = "bootstrap")
#no issues
# Calculate the estimated marginal means for Taxonomic.group
emmeans_tax_size <- emmeans(model2, ~ Land.cover)

# Perform pairwise comparisons
pairwise_comparisons <- contrast(emmeans_tax_size, method = "pairwise")

# View the results of the pairwise comparisons
summary(pairwise_comparisons)

#Isopoda - Opisthopora (p = 0.0282) → Isopoda has significantly more than worms


#plot the data for 

# Calculate Percentage and Total by Taxonomic Group
df_percentage_total_trophic <- df_presence %>%
  group_by(Trophic.level) %>%
  summarise(Count = sum(Plastic.present), Total = n(), .groups = "drop") %>%
  mutate(Percentage = (Count / Total) * 100)

wilson_ci2 <- binom.confint(x = df_percentage_total_trophic$Count, 
                            n = df_percentage_total_trophic$Total, 
                            method = "wilson")


# Add confidence intervals to the data frame
df_percentage_total_trophic <- df_percentage_total_trophic %>%
  mutate(lower = wilson_ci2$lower,
         upper = wilson_ci2$upper)

#have the labels a standard distance from the error bars for all bars 
df_percentage_total_trophic <- df_percentage_total_trophic %>%
  mutate(label_y = upper * 100 + 3.5)

df_percentage_total_trophic$Trophic.level <- factor(df_percentage_total_trophic$Trophic.level, 
                                                        levels = c("Herbivore","Omnivore","Carnivore", "Detritivore"))

# Create the plot
ggplot(df_percentage_total_trophic, aes(x = Trophic.level, y = Percentage, fill = Trophic.level)) +
  geom_bar(stat = "identity", fill = "grey90") +
  geom_errorbar(aes(ymin = lower * 100, ymax = upper * 100), width = 0.4) +
  geom_text(aes(y = label_y, label = paste(round(Count), "/", Total)), 
            size = 4, position = position_dodge(width = 0.9), fontface = "plain") +
  labs(x = "Trophic level", y = "% Plastic Presence") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    panel.grid = element_blank(),           # Removes default grid lines
    axis.line = element_line(color = "black") # Adds black axis lines
  ) +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.1)))


# Calculate standard 95% confidence intervals for both taxonomic and trophic groups
df_95ci_taxonomic <- df_percentage_total_taxonomic %>%
  mutate(
    Std_Error = sqrt((Count / Total) * (1 - (Count / Total)) / Total), # Standard error for proportion
    CI_Lower = (Count / Total) - 1.96 * Std_Error, # Lower bound
    CI_Upper = (Count / Total) + 1.96 * Std_Error  # Upper bound
  ) %>%
  mutate(
    CI_Lower = pmax(0, CI_Lower),      # Prevent lower bound from being negative
    CI_Upper = pmin(1, CI_Upper)       # Prevent upper bound from exceeding 1
  ) %>%
  select(Taxonomic.group, Count, Total, Percentage, CI_Lower, CI_Upper)

# Calculate standard 95% confidence intervals
df_95ci_trophic <- df_percentage_total_trophic %>%
  mutate(
    Std_Error = sqrt((Count / Total) * (1 - (Count / Total)) / Total), # Standard error for proportion
    CI_Lower = (Count / Total) - 1.96 * Std_Error, # Lower bound
    CI_Upper = (Count / Total) + 1.96 * Std_Error  # Upper bound
  ) %>%
  mutate(
    CI_Lower = pmax(0, CI_Lower),      # Prevent lower bound from being negative
    CI_Upper = pmin(1, CI_Upper)       # Prevent upper bound from exceeding 1
  ) %>%
  select(Trophic.level, Count, Total, Percentage, CI_Lower, CI_Upper)

# Calculate Percentage and Total 
df_percentage_total_land <- df_presence_land %>%
  group_by(Land.cover) %>%
  summarise(Count = sum(Plastic.present), Total = n(), .groups = "drop") %>%
  mutate(Percentage = (Count / Total) * 100)

wilson_ci2 <- binom.confint(x = df_percentage_total_land$Count, 
                            n = df_percentage_total_land$Total, 
                            method = "wilson")


# Add confidence intervals to the data frame
df_percentage_total_land <- df_percentage_total_land %>%
  mutate(lower = wilson_ci2$lower,
         upper = wilson_ci2$upper)

#have the labels a standard distance from the error bars for all bars 
df_percentage_total_land <- df_percentage_total_land %>%
  mutate(label_y = upper * 100 + 3.5)

ggplot(df_percentage_total_land, aes(x = Land.cover, y = Percentage, fill = Land.cover)) +
  geom_bar(stat = "identity", fill = "grey90") +
  geom_errorbar(aes(ymin = lower * 100, ymax = upper * 100), width = 0.4) +
  geom_text(aes(y = label_y, label = paste(round(Count), "/", Total)), 
            size = 4, position = position_dodge(width = 0.9), fontface = "plain") +
  labs(x = "Land cover", y = "% Plastic Presence") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    panel.grid = element_blank(),           # Removes default grid lines
    axis.line = element_line(color = "black") # Adds black axis lines
  ) +
  scale_y_continuous(
    limits = c(0, 105),          # Extend upper limit slightly
    breaks = seq(0, 100, 25),    # Ensure only integer values appear
    expand = expansion(mult = c(0, 0.1)) # Adjust bar spacing
  )

#soil



modelsoil <- glmmTMB(
  Plastic.present ~  (1 | Site) , 
  data = df_presence_land, 
  family = binomial (link = logit))

summary(modelsoil)

modeltro2 <- glmmTMB(
  Plastic.present ~  Land.cover +(1 | Site), 
  data = df_presence_land, 
  family = binomial (link = logit)
)
summary(modeltro2)

anova(modelsoil, modeltro2)

#with offset chi = 9.31, p = 0.053

check_collinearity(modelsoil)

# Simulate residuals
sim_res2 <- simulateResiduals(fittedModel = modelsoil)

# Plot diagnostics
plot(sim_res2)
#QQplot shows good fit

# Overdispersion test using DHARMa
testDispersion(sim_res2)
testOutliers(sim_res2, type = "bootstrap")
#no issues


