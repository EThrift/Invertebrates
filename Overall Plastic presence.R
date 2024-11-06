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

#load dataframe 

df_presence <- read.csv("Plastic_presence.csv")


# Remove spaces from column names
colnames(df_presence) <- make.names(colnames(df_presence))

model1 <- glmmTMB(Plastic.present ~  Land.cover + (1 | Site), 
               data = df_presence, 
               family = binomial(link = "logit"))
summary(model1)


model2 <- glmmTMB(Plastic.present ~ Taxonomic.group + Land.cover + (1 | Site), 
                data = df_presence, 
                family = binomial(link = "logit"))
summary(model2)

anova(model1, model2)

#choose model 2 lower AIC values

#model 2 p value <0.001

# Load the DHARMa package
library(DHARMa)

# Simulate residuals
sim_res <- simulateResiduals(fittedModel = model2)

# Plot diagnostics
plot(sim_res)
#QQplot shows good fit

# Overdispersion test using DHARMa
testDispersion(sim_res)
#no evidence of overdispersion p = 0.88

library(emmeans)
#pairwise comparisons 
# Calculate the estimated marginal means for Taxonomic.group
emmeans_tax_group <- emmeans(model2, ~ Taxonomic.group)

# Perform pairwise comparisons
pairwise_comparisons <- contrast(emmeans_tax_group, method = "pairwise")

# View the results of the pairwise comparisons
summary(pairwise_comparisons)



library(binom)


# Calculate Percentage and Total by Taxonomic Group
df_percentage_total_taxonomic <- df_presence %>%
  group_by(Taxonomic.group) %>%
  summarise(Count = sum(Plastic.present), Total = n(), .groups = "drop") %>%
  mutate(Percentage = (Count / Total) * 100)

# Calculate wilson confidence intervals
df_percentage_total_taxonomic <- df_percentage_total_taxonomic %>%
  rowwise() %>%
  mutate(
    Wilson_CI = list(binom.confint(Count, Total, method = "wilson")),
    Lower = Wilson_CI$lower,
    Upper = Wilson_CI$upper
  ) %>%
  ungroup() %>%
  select(-Wilson_CI)

#have the labels a standard distance from the error bars for all bars 
df_percentage_total_taxonomic <- df_percentage_total_taxonomic %>%
  mutate(label_y = Upper * 100 + 3.5)

# Base plot
library(ggplot2)

# Convert 'Taxonomic.group' to a factor with the specified order
df_percentage_total_taxonomic$Taxonomic.group <- factor(df_percentage_total_taxonomic$Taxonomic.group, 
                                                        levels = c("Lepidoptera", "Isopoda", "Coleoptera", 
                                                                   "Stylommatophora", "Hemiptera", "Opisthopora"))

# Create the plot
ggplot(df_percentage_total_taxonomic, aes(x = Taxonomic.group, y = Percentage, fill = Taxonomic.group)) +
  geom_bar(stat = "identity", fill = "grey90") +
  geom_errorbar(aes(ymin = Lower * 100, ymax = Upper * 100), width = 0.4) +
  geom_text(aes(y = label_y, label = paste(round(Count), "/", Total)), 
            size = 5, position = position_dodge(width = 0.9), fontface = "plain") +
  labs(x = "Taxonomic group", y = "% Plastic Presence") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16)
  ) +
  scale_y_continuous(limits = c(0, 60), expand = expansion(mult = c(0, 0.1)))



modeltro <- glmmTMB(Plastic.present ~  Land.cover + (1 | Site), 
                data = df_presence, 
                family = binomial(link = "logit"))
summary(modeltro)

modeltro2 <- glmmTMB(Plastic.present ~ Trophic.level + Land.cover + (1 | Site), 
                  data = df_presence, 
                  family = binomial(link = "logit"))
summary(modeltro2)

anova(modeltro, modeltro2)

#modeltro2 is significantly better p = 0.047


# Simulate residuals
sim_res2 <- simulateResiduals(fittedModel = modeltro2)

# Plot diagnostics
plot(sim_res2)
#QQplot shows good fit

# Overdispersion test using DHARMa
testDispersion(sim_res2)
#no evidence of overdispersion p = 0.92


#pairwise comparisons 
# Calculate the estimated marginal means for Taxonomic.group
emmeans_trophic_group <- emmeans(modeltro2, ~ Trophic.level)

# Perform pairwise comparisons
pairwise_comparisons2 <- contrast(emmeans_trophic_group, method = "pairwise")

# View the results of the pairwise comparisons
summary(pairwise_comparisons2)
#No significant differences seen between trophic levels

#plot the figures showing percentage and wilson confidence intervals for taxonomic group and trophic level




#plot the data for 

# Calculate Percentage and Total by Taxonomic Group
df_percentage_total_trophic <- df_presence %>%
  group_by(Trophic.level) %>%
  summarise(Count = sum(Plastic.present), Total = n(), .groups = "drop") %>%
  mutate(Percentage = (Count / Total) * 100)

# Calculate wilson confidence intervals
df_percentage_total_trophic <- df_percentage_total_trophic %>%
  rowwise() %>%
  mutate(
    Wilson_CI = list(binom.confint(Count, Total, method = "wilson")),
    Lower = Wilson_CI$lower,
    Upper = Wilson_CI$upper
  ) %>%
  ungroup() %>%
  select(-Wilson_CI)

#have the labels a standard distance from the error bars for all bars 
df_percentage_total_trophic <- df_percentage_total_trophic %>%
  mutate(label_y = Upper * 100 + 3.5)

# Create the plot
 ggplot(df_percentage_total_trophic, aes(x = Trophic.level, y = Percentage, fill = Trophic.level)) +
  geom_bar(stat = "identity", fill = "grey90") +  # Set fill color to dark grey
  geom_errorbar(aes(ymin = Lower * 100, ymax = Upper * 100), width = 0.4) +  # Use the correct column names
  geom_text(aes(y = label_y, label = paste(round(Count), "/", Total)),
            size = 5, position = position_dodge(width = 0.9), fontface = "plain") +  # Add text labels in plain font
  labs(x = "Trophic level", y = "% Plastic Presence") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 14),  # Adjust x-axis tick labels text size
    axis.text.y = element_text(size = 14),  # Adjust y-axis tick labels text size
    axis.title.x = element_text(size = 16),  # Adjust x-axis title text size
    axis.title.y = element_text(size = 16)   # Adjust y-axis title text size
  )



