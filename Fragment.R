#created on 22/08/24
# By Emily Thrift
# script for a glm to assess the differences in the shape of plastics between taxonomic groups and trophic levels

#load libraries
library(tidyverse)
library(ggplot2)
library(lme4)
install.packages("glmmTMB")
library(glmmTMB)
library(emmeans)

#add in working location

setwd('C:/Users/Emily Thrift/Desktop/Inverts')
list.files()

#load dataframe 

df_shape <- read.csv("Fragment.csv")


# Remove spaces from column names
colnames(df_shape) <- make.names(colnames(df_shape))

df_shape$Biomass_log <- log(df_shape$Biomass)

# Zero-inflated model using glmmTMB
Fragment_tax1 <- glmmTMB(Fragment
                                      ~ (1 | Site) + Land.cover + offset(Biomass_log), 
                                      data = df_shape, 
                                      family = binomial(link = "logit"))  

# View the summary 
summary(Fragment_tax1)

# Another model with Taxonomic.group added
Fragment_tax <- glmmTMB(Fragment
                                   ~Land.cover + Taxonomic.group + (1 | Site) + offset(Biomass_log), 
                                   data = df_shape, 
                                   family = binomial(link = "logit"))

# View the summary of this model
summary(Fragment_tax_simplified)

# Anova to compare models
anova(Fragment_tax1, Fragment_tax_simplified)
check_collinearity(Fragment_tax1)

#taxonomic group did not have a significant impact on model chi = 0.29, p = 0.99

# Simulate residuals
sim_res <- simulateResiduals(fittedModel = Fragment_tax_zero_inflated)

# Plot diagnostics
plot(sim_res)


# Overdispersion test using DHARMa
testDispersion(sim_res)
#not overdispersed p = 0.62




#check for perfect separation 
table(df_shape$Trophic.level, df_shape$Fragment)


# Exclude rows where Trophic level is omnivore due to a lack of results 
df_shape_filtered <- df_shape[!(df_shape$Trophic.level %in% c("Omnivore")), ]

# Ensure Biomass is properly log-transformed and scaled
df_shape_filtered$Biomass_log <- scale(log(df_shape_filtered$Biomass))

# Now fit the model
Fragment_tro <- glmmTMB(Fragment
                        ~ Land.cover + (1 | Site) + offset(Biomass_log), 
                        data = df_shape_filtered, 
                        family = binomial(link = "logit"))

# View the summary of the model
summary(Fragment_tro)


Fragment_tro_simplified <- glmmTMB(Fragment
                                   ~  Trophic.level + Land.cover+(1 | Site) + offset(Biomass_log), 
                                   data = df_shape_filtered, 
                                   family = binomial(link = "logit"))

summary(Fragment_tro_simplified)

anova(Fragment_tro,Fragment_tro_simplified)

#trophic level did not have a significant impact chi = 0.24, p = 0.88
check_collinearity(Fragment_tro)
# Simulate residuals
sim_res <- simulateResiduals(fittedModel = Fragment_tro)

# Plot diagnostics
plot(sim_res)


# Overdispersion test using DHARMa
testDispersion(sim_res)



