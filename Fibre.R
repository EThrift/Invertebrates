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

Fibre_tax <- glmmTMB(Fibre
                                ~  Land.cover + (1 | Site), 
                                data = df_shape, 
                                family = binomial(link = "logit"))
summary(Fibre_tax)

Fibre_tax_simplified <- glmmTMB(Fibre
                     ~ Land.cover +Taxonomic.group + (1 | Site), 
                     data = df_shape, 
                     family = binomial(link = "logit"))
summary(Fibre_tax_simplified)

anova(Fibre_tax, Fibre_tax_simplified)

# taxonomic group did not show signifcant impact p = 0.56

# Simulate residuals
sim_res <- simulateResiduals(fittedModel = Fibre_tax_simplified)

# Plot diagnostics
plot(sim_res)
#QQplot shows good fit 

# Overdispersion test using DHARMa
testDispersion(sim_res)
#not overdispersed p = 0.88

#check for perfect separation 
table(df_shape$Trophic.level, df_shape$Fibre)


# Exclude rows where Trophic level is omnivore due to a lack of results 
df_shape_filtered <- df_shape[!(df_shape$Trophic.level %in% c("Omnivore")), ]

Fibre_tro <- glmmTMB(Fibre
                   ~ Land.cover + (1 | Site), 
                   data = df_shape_filtered, 
                   family = binomial(link = "logit"))
summary(Fibre_tro)

Fibre_tro_simplified <- glmmTMB(Fibre
                     ~ Land.cover + Trophic.level +  (1 | Site), 
                     data = df_shape_filtered, 
                     family = binomial(link = "logit"))
summary(Fibre_tro_simplified)

anova(Fibre_tro, Fibre_tro_simplified)

#trophic level did not have significant effect

# Simulate residuals
sim_res <- simulateResiduals(fittedModel = Fibre_tro_simplified)

# Plot diagnostics
plot(sim_res)
#QQplot shows bad fit deviations detected

# Overdispersion test using DHARMa
testDispersion(sim_res)
