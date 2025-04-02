#created on 22/08/24
# By Emily Thrift
# script for a glm to assess the differences in the fibre of plastics between taxonomic groups and trophic levels

#load libraries
library(tidyverse)
library(ggplot2)
library(lme4)
library(glmmTMB)
library(emmeans)

#add in working location

setwd('C:/Users/Emily Thrift/Desktop/Inverts')
list.files()

#load dataframe 

df_shape <- read.csv("Fibre.csv")


# Remove spaces from column names
colnames(df_shape) <- make.names(colnames(df_shape))

df_shape %>%
  group_by(Land.cover, Fibre) %>%
  summarise(count = n()) %>%
  spread(key = Fibre, value = count, fill = 0)


df_shape$Biomass_log <- log(df_shape$Biomass + 1)  # Offset

table(df_shape$Fibre, df_shape$Land.cover)
table(df_shape$Fibre, df_shape$Trophic.level)

# Subset the data to remove 'Omnivore' from Trophic.level
df_shape <- df_shape[df_shape$Trophic.level != "Omnivore", ]

model1 <- glmmTMB(Fibre ~   (1 | Site) +Land.cover + offset(Biomass_log), data = df_shape, 
                  family = binomial(link = logit))
summary(model1)

model2 <- glmmTMB(Fibre ~ Land.cover + Trophic.level + (1 | Site) + offset(Biomass_log), 
                  data = df_shape, 
                  family =binomial(link = logit))
summary(model2)
anova(model1, model2)

#chi = 5.62, p = 0.06

# Load the DHARMa package
library(DHARMa)

# Simulate residuals
sim_res <- simulateResiduals(fittedModel = model1)


# Plot diagnostics
plot(sim_res)
#QQplot shows good fit

# Overdispersion test using DHARMa
testDispersion(sim_res)
#no evidence of overdispersion p = 0.76

testOutliers(sim_res, type = "bootstrap")
#no issues

testZeroInflation(sim_res)

df_shape %>%
  group_by(Taxonomic.group, Fibre) %>%
  summarise(count = n()) %>%
  spread(key = Fibre, value = count, fill = 0)


model1 <- glmmTMB(Fibre ~  (1 | Site) + Land.cover + offset(Biomass_log), data = df_shape, 
                  family =binomial(link = logit))
summary(model1)

model2 <- glmmTMB(Fibre ~ Taxonomic.group + Land.cover + offset(Biomass_log), 
                  data = df_shape, 
                  family = binomial(link = logit))
summary(model2)
anova(model1, model2)

# model 1 is better chi = 3.88, p = 0.42 

library(performance)

# Check for multicollinearity in the entire model
check_collinearity(model1)
# Load the DHARMa package
library(DHARMa)

# Simulate residuals
sim_res <- simulateResiduals(fittedModel = model1)


# Plot diagnostics
plot(sim_res)
#QQplot shows good fit

# Overdispersion test using DHARMa
testDispersion(sim_res)
#no evidence of overdispersion p = 0.76

testOutliers(sim_res, type = "bootstrap")
#no issues

testZeroInflation(sim_res)

#check for perfect separation 
table(df_shape$Trophic.level, df_shape$Fibre)

# Exclude rows where Trophic level is omnivore due to a lack of results 
df_shape_filtered <- df_shape[!(df_shape$Trophic.level %in% c("Omnivore")), ]

df_shape_filtered$Biomass_log <- log(df_shape_filtered$Biomass + 1)  # Offset


Fibre_tro <- glmmTMB(Fibre
                   ~  (1 | Site) + Land.cover + offset(Biomass_log), 
                   data = df_shape_filtered, 
                   family = binomial(link=logit))
summary(Fibre_tro)

Fibre_tro_simplified <- glmmTMB(Fibre 
                     ~  Trophic.level + Land.cover + (1 | Site) + offset(Biomass_log), 
                     data = df_shape_filtered, 
                     family = binomial(link=logit))
summary(Fibre_tro_simplified)

anova(Fibre_tro, Fibre_tro_simplified)
library(performance)

# Check for multicollinearity in the entire model
check_collinearity(Fibre_tro)
#trophic level did not significant effect chi = 5.62, p = 0.061
# Simulate residuals
sim_res <- simulateResiduals(fittedModel = Fibre_tro)

# Plot diagnostics
plot(sim_res)
#QQplot shows bad fit deviations detected

# Overdispersion test using DHARMa
testDispersion(sim_res)

