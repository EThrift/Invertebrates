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
Fragment_tax <- glmmTMB(Fragment
                   ~ Land.cover + (1 | Site), 
                   data = df_shape, 
                   family = binomial(link = "logit"))
summary(Fragment_tax)


Fragment_tax_simplified <- glmmTMB(Fragment
                                 ~Land.cover +  Taxonomic.group + (1 | Site), 
                                 data = df_shape, 
                                 family = binomial(link = "logit"))

summary(Fragment_tax_simplified)

anova(Fragment_tax, Fragment_tax_simplified)

#taxonomic group did not have a significant impact on model p = 0.076

# Simulate residuals
sim_res <- simulateResiduals(fittedModel = Fragment_tax)

# Plot diagnostics
plot(sim_res)


# Overdispersion test using DHARMa
testDispersion(sim_res)
#not overdispersed p = 0.62




#check for perfect separation 
table(df_shape$Trophic.level, df_shape$Fragment)

# Exclude rows where Trophic level is omnivore due to a lack of results 
df_shape_filtered <- df_shape[!(df_shape$Trophic.level %in% c("Omnivore")), ]


Fragment_tro <- glmmTMB(Fragment
                        ~ Land.cover + (1 | Site), 
                        data = df_shape_filtered, 
                        family = binomial(link = "logit"))
summary(Fragment_tro)


Fragment_tro_simplified <- glmmTMB(Fragment
                                   ~Land.cover + Trophic.level + (1 | Site), 
                                   data = df_shape_filtered, 
                                   family = binomial(link = "logit"))

summary(Fragment_tro_simplified)

anova(Fragment_tro,Fragment_tro_simplified)

#trophic level did not have a significant impact p = 0.11

# Simulate residuals
sim_res <- simulateResiduals(fittedModel = Fragment_tro)

# Plot diagnostics
plot(sim_res)
#QQplot shows bad fit deviations detected

# Overdispersion test using DHARMa
testDispersion(sim_res)

#not overdispersed p = 0.76


#pairwise comparisons 
# Calculate the estimated marginal means for Trophic.level
emmeans_tro_group <- emmeans(Fragment_tro, ~ Trophic.level)




