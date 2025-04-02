#created on 20/08/24
# By Emily Thrift
# script for a glm to assess the ratio of biomass

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

df_ratio <- read.csv("Actual MP ratio.csv")
df_weight <- read.csv("Plastic_weight_all.csv")
df_fragment <- read.csv("Fragment.csv")
df_fibre <- read.csv ("Fibre.csv")
df_size <- read.csv("Size.csv")


# Remove spaces from column names
colnames(df_ratio) <- make.names(colnames(df_ratio))
colnames(df_weight) <- make.names(colnames(df_weight))
colnames(df_fragment) <- make.names(colnames(df_fragment))
colnames(df_fibre) <- make.names(colnames(df_fibre))
colnames(df_size) <- make.names(colnames(df_size))
# Subset the data to remove 'Omnivore' from Trophic.level
df_fibre <- df_fibre[df_fibre$Trophic.level != "Omnivore", ]

# Subset the data to remove 'Omnivore' from Trophic.level
df_weight <- df_weight[df_weight$Trophic.level != "Omnivore", ]

# Subset the data to remove 'Omnivore' from Trophic.level
df_ratio <- df_ratio[df_ratio$Trophic.level != "Omnivore", ]

# Subset the data to remove 'Omnivore' from Trophic.level
df_fragment <- df_fragment[df_fragment$Trophic.level != "Omnivore", ]

# Subset the data to remove 'Omnivore' from Trophic.level
df_size <- df_size[df_size$Trophic.level != "Omnivore", ]
# Check the levels of Trophic.level after filtering
table(df_fibre$Trophic.level)

df_fragment$Biomass_log <- log(df_fragment$Biomass + 1)  # Offset
df_fragment$Average.biomass.per.individual <- log(df_fragment$Average.biomass.per.individual + 1)

df_fibre$Biomass_log <- log(df_fibre$Biomass + 1)  # Offset
df_fibre$Average.biomass.per.individual <- log(df_fibre$Average.biomass.per.individual + 1)

df_ratio$Biomass_log <- log(df_ratio$Biomass + 1)  # Offset
df_ratio$Average.biomass.per.individual <- log(df_ratio$Average.biomass.per.individual + 1)

df_weight$Biomass_log <- log(df_weight$Biomass + 1)  # Offset
df_weight$Average.biomass.per.individual <- log(df_weight$Average.biomass.per.individual + 1)

df_size$Biomass_log <- log(df_size$Biomass + 1)  # Offset
df_size$Average.biomass <- log(df_size$Average.biomass + 1)

model1 <- lmer(Size ~  Average.biomass + Land.cover + (1 | Site) +offset(Biomass_log), 
                  data = df_size, REML = TRUE)

summary(model1)

model2 <- lmer(Size ~ Average.biomass  +  Land.cover + Trophic.level +  (1 | Site) +offset(Biomass_log),
                  data = df_size, REML = TRUE)
summary(model2)
anova(model1, model2)
anova(model2)
#model 2 better chi = 24.76 p = <0.0001. F score 6.33, p = <0.0001

# Load the DHARMa package
library(DHARMa)

# Simulate residuals
sim_res <- simulateResiduals(fittedModel = model2)


# Plot diagnostics
plot(sim_res)
#QQplot shows good fit

# Overdispersion test using DHARMa
testDispersion(sim_res)
#no evidence of overdispersion p = 0.96


anova(model2)

#F = 14.59, p = <0.0001

# Calculate the estimated marginal means for Trophic level
emmeans_tax_size <- emmeans(model2, ~ Taxonomic.group)

# Perform pairwise comparisons
pairwise_comparisons <- contrast(emmeans_tax_size, method = "pairwise")

# View the results of the pairwise comparisons
summary(pairwise_comparisons)

#Carnivore - Detritivore	0.0002 (significant)
#Carnivore - Herbivore	<0.0001 (highly significant)

model1 <- lmer(Size ~  Average.biomass + Land.cover+ (1 | Site) +offset(Biomass_log), 
                  data = df_size, REML = FALSE)

summary(model1)

model2 <- lmer(Size ~ Average.biomass  +  Taxonomic.group + Land.cover + (1 | Site) +offset(Biomass_log), 
                  data = df_size, REML = FALSE)
summary(model2)
anova(model1, model2)

#model 2 betterchi= 29.61, p = <0.0001

anova(model2)

#F value = 6.85, p <0.0001

# Load the DHARMa package
library(DHARMa)

# Simulate residuals
sim_res <- simulateResiduals(fittedModel = model2)


# Plot diagnostics
plot(sim_res)
#QQplot shows good fit

# Overdispersion test using DHARMa
testDispersion(sim_res)
#no evidence of overdispersion p = 0.95

# Calculate the estimated marginal means for Taxonomic.group
emmeans_tax_size <- emmeans(model2, ~ Taxonomic.group)

# Perform pairwise comparisons
pairwise_comparisons <- contrast(emmeans_tax_size, method = "pairwise")

# View the results of the pairwise comparisons
summary(pairwise_comparisons)

#Coleoptera - Opisthopora	0.0211	
#Coleoptera - Stylommatophora	<0.0001	

#Coleoptera significantly larger

str(df_size)  # Check structure of the dataset
head(df_size) # View first few rows

str(df_ratio)


model1 <- lmer(Average.MP.per.individual ~ Average.biomass.per.individual + Land.cover + (1 | Site) + offset(Biomass_log), 
               data = df_ratio, REML = FALSE)

summary(model1)

model2 <- lmer(Average.MP.per.individual ~ Average.biomass.per.individual + Land.cover + Taxonomic.group +  (1 | Site) + offset(Biomass_log), data = df_ratio, REML = FALSE)

summary(model2)
anova(model1, model2)

#model 2 tax not sigificant CHI = 10.19, P = 0.069

anova(model2)

#F value = 2.19, p = 0.064


# Load the DHARMa package
library(DHARMa)

# Simulate residuals
sim_res <- simulateResiduals(fittedModel = model1)


# Plot diagnostics
plot(sim_res)
#QQplot shows good fit

# Overdispersion test using DHARMa
testDispersion(sim_res)
#no evidence of overdispersion p = 0.95


model1 <- lmer(Average.MP.per.individual ~ Average.biomass.per.individual + Land.cover + (1 | Site) + offset(Biomass_log), 
               data = df_ratio, REML = FALSE)

summary(model1)

model2 <- lmer(Average.MP.per.individual ~ Average.biomass.per.individual + Trophic.level  + Land.cover+ (1 | Site)  + offset(Biomass_log), 
                  data = df_ratio)
summary(model2)
anova(model1, model2)
anova(model2)

#model 1 is better average weight does not impact number of pieces of plastic model chi = 4.68, p = 0.09

anova(model2)

#f stat 3.42, p = 0.10


# Simulate residuals
sim_res2 <- simulateResiduals(fittedModel = model1)

# Plot diagnostics
plot(sim_res2)
#QQplot shows good fit

# Overdispersion test using DHARMa
testDispersion(sim_res2)
testOutliers(sim_res2, type = "bootstrap")

#pairwise comparisons 
# Calculate the estimated marginal means for Taxonomic.group
emmeans_tro_size <- emmeans(model2, ~ Land.cover)

# Perform pairwise comparisons
pairwise_comparisons <- contrast(emmeans_tro_size, method = "pairwise")

# View the results of the pairwise comparisons
summary(pairwise_comparisons)

model1 <- lmer(Average.MP.per.individual ~ Average.biomass.per.individual + Taxonomic.group + (1 | Site) + offset(Biomass_log), 
               data = df_ratio, REML = FALSE)

model2 <- lmer(Average.MP.per.individual ~ Average.biomass.per.individual + Taxonomic.group + Land.cover + (1 | Site) + offset(Biomass_log), 
               data = df_ratio, REML = FALSE)
summary(model2)
anova(model1, model2)
#model 1 is better average weight does not impact number of pieces of plastic model when including offset chi = 8.35, p = 0.039
anova(model2)
#F value 2.99, p = 0.036


# Load the DHARMa package
library(DHARMa)

# Simulate residuals
sim_res <- simulateResiduals(fittedModel = model2)


# Plot diagnostics
plot(sim_res)
#QQplot shows good fit

# Overdispersion test using DHARMa
testDispersion(sim_res)
#no evidence of overdispersion p = 0.95


model1 <- lmer(MPs.per.gram ~  Land.cover + (1 | Site) +offset(Biomass_log),
                  data = df_weight, REML = FALSE)

model2 <- lmer(MPs.per.gram ~  Taxonomic.group +  Land.cover + (1 | Site) + offset(Biomass_log), 
                  data = df_weight, REML = FALSE)
summary(model2)

anova(model1, model2)
anova(model2)
#model 1 taxonomic group does not impact mps per gram chi = 4.13 p = 0.53, F score = 0.83, p = 0.76
# Simulate residuals
sim_res2 <- simulateResiduals(fittedModel = model1)

# Plot diagnostics
plot(sim_res2)
#QQplot shows good fit


# Overdispersion test using DHARMa
testDispersion(sim_res2)
testOutliers(sim_res2, type = "bootstrap")

model1 <- lmer(MPs.per.gram ~  Land.cover + (1 | Site) +offset(Biomass_log),
               data = df_weight, REML = FALSE)

model2 <- lmer(MPs.per.gram ~  Trophic.level +  Land.cover + (1 | Site) + offset(Biomass_log), 
               data = df_weight, REML = FALSE)
summary(model2)

anova(model1, model2)
library(car)
vif(model2)

#model 1 taxonomic group does not impact plastic presence chi = 0.56  p = 0.75,F value = 0.28, p = 0.70
# Simulate residuals
sim_res2 <- simulateResiduals(fittedModel = model1)

# Plot diagnostics
plot(sim_res2)
#QQplot shows good fit


# Overdispersion test using DHARMa
testDispersion(sim_res2)
testOutliers(sim_res2, type = "bootstrap")


model1 <- glmmTMB(Plastic.present ~   (1 | Site) + Trophic.level + offset(Biomass_log), 
                  family = binomial(link = "logit"), 
                  data = df_weight)
summary(model1)
model2 <- glmmTMB(Plastic.present ~ Average.biomass.per.individual+ (1 | Site) + Trophic.level + Land.cover + offset(Biomass_log), 
                  family = binomial(link = "logit"), 
                  data = df_weight)

summary(model2)
anova(model1, model2)

library(car)
vif(model2)

#model 1 average biomass does impact plastic presence chi = 0.46 p = 0.49

# Load the DHARMa package
library(DHARMa)

# Simulate residuals
sim_res <- simulateResiduals(fittedModel = model1)

# Plot diagnostics
plot(sim_res)
#QQplot shows good fit

# Overdispersion test using DHARMa
testDispersion(sim_res)
#no evidence of overdispersion p = 0.88

# Summary statistics for Fibre and Fragments
summary(df_fibre$Fibre)
summary(df_fragment$Fragment)  # Assuming df_fragments exists for fragments

cor(df_fibre$Average.biomass.per.individual, df_fibre$Biomass)

# Check if Taxonomic.group is a factor
df_fibre$Taxonomic.group <- as.factor(df_fibre$Taxonomic.group)

# Check the levels of the factor
levels(df_fibre$Taxonomic.group)


df_shape %>%
  group_by(Taxonomic.group, Fibre) %>%
  summarise(count = n()) %>%
  spread(key = Fibre, value = count, fill = 0)

#Scale the continuous variable (Average.biomass.per.individual)
df_fibre$Average.biomass.per.individual <- scale(df_fibre$Average.biomass.per.individual)



model1 <- glmmTMB(Fibre ~ Average.biomass.per.individual + Land.cover + (1 | Site) + offset(Biomass_log), data = df_fibre, 
                  family = binomial(link = logit))
summary(model1)

# Re-run the model with the d variable
model2 <- glmmTMB(Fibre ~ Taxonomic.group + Average.biomass.per.individual + Land.cover + (1 | Site) +offset(Biomass_log), 
                         data = df_fibre, 
                         family = binomial (link = logit))
summary(model2)

# Check the summary of the model

anova(model1, model2)
library(car)
vif(model2)
#model 1 is better average weight does not impact number of fragments  model chi 3.90, p = 0.56


#choose model 1 better without interaction 

# Load the DHARMa package
library(DHARMa)

# Simulate residuals
sim_res <- simulateResiduals(fittedModel = model1)


# Plot diagnostics
plot(sim_res)
#QQplot shows good fit

# Overdispersion test using DHARMa
testDispersion(sim_res)
#no evidence of overdispersion p = 0.96

#no evidence of overdispersion p = 0.96
model1 <- glmmTMB(Fragment ~ Average.biomass.per.individual  + (1 | Site) + Land.cover + offset(Biomass_log), 
                  family = binomial (link = "logit"), 
                  data = df_fragment)

model2 <- glmmTMB(Fragment ~ Average.biomass.per.individual + Taxonomic.group + Land.cover +(1 | Site)  + offset(Biomass_log), 
                  family = binomial (link = "logit"), 
                  data = df_fragment)

summary(model2)
anova(model1, model2)
library(car)
vif(model2)
#model 1 is better average weight does not impact number of fragments  model chi 3.16, p = 0.67

#choose model 1 better without interaction 

# Load the DHARMa package
library(DHARMa)

# Simulate residuals
sim_res <- simulateResiduals(fittedModel = model1)


# Plot diagnostics
plot(sim_res)
#QQplot shows good fit

# Overdispersion test using DHARMa
testDispersion(sim_res)
#no evidence of overdispersion p = 0.88

model1 <- glmmTMB(Fragment ~ Average.biomass.per.individual  + (1 | Site)  + Land.cover + offset(Biomass_log), 
                  family = binomial (link = "logit"), 
                  data = df_fragment)

model2 <- glmmTMB(Fragment ~ Average.biomass.per.individual + Taxonomic.group + (1 | Site)  + offset(Biomass_log), 
                  family = binomial (link = "logit"), 
                  data = df_fragment)

summary(model2)
anova(model1, model2)
library(car)
vif(model2)
#chi 2.71, p = 0.25

#choose model 1 better without interaction 

# Load the DHARMa package
library(DHARMa)

# Simulate residuals
sim_res <- simulateResiduals(fittedModel = model1)


# Plot diagnostics
plot(sim_res)
#QQplot shows good fit

# Overdispersion test using DHARMa
testDispersion(sim_res)
#no evidence of overdispersion p = 0.88

# Remove rows where Trophic.level is 'Omnivore'
df_ratio <- df_ratio[df_ratio$Trophic.level != "Omnivore", ]

# Now use this new variable in your models
model1 <- glmmTMB(Number.of.pieces.of.plastics ~ Average.biomass.per.individual + Land.cover+ (1 | Site)  +  offset(Biomass_log),
                  family = gaussian,
                  data = df_ratio)

model2 <- glmmTMB(Number.of.pieces.of.plastics ~ Average.biomass.per.individual + Trophic.level + Land.cover+ (1 | Site)   + offset(Biomass_log),
                  family = gaussian,
                  data = df_ratio)

anova(model1, model2)
library(car)
vif(model2)
#model 1 better average biomass does not impact number of mps chi = 0.78, p=0.37. Without offset chi = 0.16, p = 0.68

model3 <- glmmTMB(Number.of.pieces.of.plastics ~ Average.biomass.per.individual + Trophic.level +  Average.biomass.per.individual * Trophic.level + Land.cover + (1 | Site)  +  offset(Biomass_log), 
                  family = gaussian, 
                  data = df_ratio)
summary (model3)

anova(model1, model3)

#choose model 2 lower AIC chi square 0.94, p = 0.81

# Load the DHARMa package
library(DHARMa)

# Simulate residuals
sim_res <- simulateResiduals(fittedModel = model1)

# Plot diagnostics
plot(sim_res)
#QQplot shows good fit

# Overdispersion test using DHARMa
testDispersion(sim_res)
#no evidence of overdispersion p = 0.88

model1 <- glmmTMB(MPs.per.gram ~  Land.cover + (1 | Site)  + offset(Biomass_log)  , 
                  family = gaussian, 
                  data = df_weight)

model2 <- glmmTMB(MPs.per.gram ~  Trophic.level +  Land.cover + (1 | Site) + offset(Biomass_log) , 
                  family = gaussian, 
                  data = df_weight)
summary(model2)

anova(model1, model2) 
library(car)
vif(model1)
#model 2 trophic level does not impacts model chi square 0.50, p = 0.77, without offset chi= 2.53, p = 0.46


# Load the DHARMa package
library(DHARMa)

# Simulate residuals
sim_res <- simulateResiduals(fittedModel = model1)

# Plot diagnostics
plot(sim_res)
#QQplot shows good fit

# Overdispersion test using DHARMa
testDispersion(sim_res)
#no evidence of overdispersion p = 0.88


model1 <- glmmTMB(Plastic.present ~ Average.biomass.per.individual + Land.cover + (1 | Site) + offset(Biomass_log),
                  family = binomial(link = "logit"), 
                  data = df_weight)
summary(model1)

model2 <- glmmTMB(Plastic.present ~ Average.biomass.per.individual + Trophic.level + Land.cover + (1 | Site) + offset(Biomass_log),
                  family = binomial(link = "logit"), 
                  data = df_weight)
summary(model2)


anova(model1, model2)

#model2 chi square = 16.02, p = <0.001 average biomass significantly improves model fit negative relationship between size and plastic presence. Without offset chi = 0.08, p = 0.77


model3 <- glmmTMB(Plastic.present ~ Average.biomass.per.individual + Trophic.level +  Average.biomass.per.individual * Trophic.level + Land.cover + 
                           (1 | Site) + offset(Biomass_log), 
                         family = binomial(link = "logit"), 
                         data = df_weight)

summary(model3)

anova(model1, model3)
library(car)
vif(model3)
# Simulate residuals
sim_res <- simulateResiduals(fittedModel = model1)

# Plot diagnostics
plot(sim_res)
#QQplot shows good fit

# Overdispersion test using DHARMa
testDispersion(sim_res)
#no evidence of overdispersion p = 0.82

model1 <- glmmTMB(Plastic.present ~ Average.biomass.per.individual + (1 | Site) + offset(Biomass_log),
                  family = binomial(link = "logit"), 
                  data = df_weight)
summary(model1)

model2 <- glmmTMB(Plastic.present ~ Average.biomass.per.individual + Land.cover + (1 | Site) + offset(Biomass_log),
                  family = binomial(link = "logit"), 
                  data = df_weight)
summary(model2)


anova(model1, model2)

#model2 chi square = 16.02, p = <0.001 average biomass significantly improves model fit negative relationship between size and plastic presence. Without offset chi = 0.08, p = 0.77


model3 <- glmmTMB(Plastic.present ~ Average.biomass.per.individual  +  Average.biomass.per.individual * Land.cover + Land.cover + 
                    (1 | Site) + offset(Biomass_log), 
                  family = binomial(link = "logit"), 
                  data = df_weight)

summary(model3)

anova(model2, model3)

# Simulate residuals
sim_res <- simulateResiduals(fittedModel = model2)

# Plot diagnostics
plot(sim_res)
#QQplot shows good fit

# Overdispersion test using DHARMa
testDispersion(sim_res)
#no evidence of overdispersion p = 0.82

# Calculate the estimated marginal means for Taxonomic.group
emmeans_tax_size <- emmeans(model2, ~ Land.cover)

# Perform pairwise comparisons
pairwise_comparisons <- contrast(emmeans_tax_size, method = "pairwise")

# View the results of the pairwise comparisons
summary(pairwise_comparisons)


table(df_fibre$Fibre, df_fibre$Trophic.level)
table(df_fibre$Fibre, df_fibre$Land.cover)

# Scale the continuous variable (Average.biomass.per.individual)
df_fibre$Average.biomass.per.individual <- scale(df_fibre$Average.biomass.per.individual)


model1 <- glmmTMB(Fibre
       ~  Land.cover +  (1 | Site) + offset(log(Biomass)), 
       data = df_fibre, 
       family = binomial (link = logit))
summary(model1)

model2 <- glmmTMB(Fibre ~ Average.biomass.per.individual + Trophic.level +Land.cover + (1 | Site) + offset(Biomass_log), 
                  family = binomial (link = logit), 
                  data = df_fibre)

summary(model2)
anova(model1, model2)

#model 2 is better average weight does not impact number of fibres  model chi = 2.43, p = 0.11. Without offset chi = 0.01 , p = 0.90

model3 <- glmmTMB(Fibre ~ Average.biomass.per.individual + Trophic.level +  Average.biomass.per.individual * Trophic.level + (1 | Site)+ Land.cover+ offset(Biomass_log), 
                  family = binomial (link = logit),
                  data = df_fibre)
summary(model3)

anova(model2, model3)

#choose model 2 better without interaction chi = 2.36, p = 0.30
library(car)
vif(model2)
# Load the DHARMa package
library(DHARMa)

# Simulate residuals
sim_res <- simulateResiduals(fittedModel = model2)


# Plot diagnostics
plot(sim_res)
#QQplot shows good fit

# Overdispersion test using DHARMa
testDispersion(sim_res)
#no evidence of overdispersion p = 0.92

# Calculate the estimated marginal means for Taxonomic.group
emmeans_tax_size <- emmeans(model2, ~ Trophic.level)

# Perform pairwise comparisons
pairwise_comparisons <- contrast(emmeans_tax_size, method = "pairwise")

# View the results of the pairwise comparisons
summary(pairwise_comparisons)


# Create a contingency table between Fragment and Land.cover
table(df_fragment$Fragment, df_fragment$Land.cover)


model1 <- glmmTMB(Fragment ~ Trophic.level  + (1 | Site)  + + offset(Biomass_log), 
                  family = binomial(link = "logit"), 
                  data = df_fragment)

model2 <- glmmTMB(Fragment ~ Average.biomass.per.individual + Trophic.level +  (1 | Site) + offset(Biomass_log), 
                  family = binomial (link = "logit"), 
                  data = df_fragment)

summary(model2)
anova(model1, model2)

#model 1 is better average weight does not impact number of fragments  model chi = 4.0 p = 0.045
library(car)
vif(model2)
model3 <- glmmTMB(Fragment ~ Average.biomass.per.individual + Trophic.level +  Average.biomass.per.individual * Trophic.level + (1 | Site) + offset(Biomass_log), 
                  family = binomial (link = "logit"), 
                  data = df_fragment)
summary(model3)


anova(model2, model3)

#choose model 2 better without interaction 

# Load the DHARMa package
library(DHARMa)

# Simulate residuals
sim_res <- simulateResiduals(fittedModel = model1)


# Plot diagnostics
plot(sim_res)
#QQplot shows good fit

# Overdispersion test using DHARMa
testDispersion(sim_res)
#no evidence of overdispersion

df_fragment$Biomass_log <- log(df_fragment$Biomass + 1)  # Offset
df_fragment$Average.biomass.per.individual <- log(df_fragment$Average.biomass.per.individual + 1)

model1 <- glmmTMB(Fragment ~ Average.biomass.per.individual  + (1 | Site)  + offset(Biomass_log), 
                  family = binomial (link= logit), 
                  data = df_fragment)

model2 <- glmmTMB(Fragment ~ Average.biomass.per.individual + Land.cover +  (1 | Site)+ offset(Biomass_log), 
                  family = binomial (link = "logit"), 
                  data = df_fragment)

summary(model2)
anova(model1, model2)
library(car)
vif(model2)
#model 1 is better average weight does not impact number of fragments  model chi = 1.95 p = 0.16


# Load the DHARMa package
library(DHARMa)

# Simulate residuals
sim_res <- simulateResiduals(fittedModel = model1)


# Plot diagnostics
plot(sim_res)
#QQplot shows good fit

# Overdispersion test using DHARMa
testDispersion(sim_res)
#no evidence of overdispersion p = 0.08


model1 <- glmmTMB(Number.of.pieces.of.plastics ~ Average.MP.per.individual + (1 | Site) + offset(Biomass_log) , 
                  family = gaussian, 
                  data = df_ratio)

model2 <- glmmTMB(Number.of.pieces.of.plastics ~ Average.MP.per.individual +  Land.cover + (1 | Site) + offset(Biomass_log), 
                  family = gaussian, 
                  data = df_ratio)

anova(model1, model2)

#model 1 better trophic level does not impact number of mps chi = 2.44, 0.48 in model without offset chi = 4.04, p = 0.25
library(car)
vif(model2)
model3 <- glmmTMB(Number.of.pieces.of.plastics ~ Average.biomass.per.individual + Trophic.level +  Average.biomass.per.individual * Trophic.level + Land.cover + (1 | Site) + offset(Biomass_log) , 
                  family = gaussian, 
                  data = df_ratio)

anova(model1, model3)

#choose model 2 lower AIC chi square 0, p = 1

# Load the DHARMa package
library(DHARMa)

# Simulate residuals
sim_res <- simulateResiduals(fittedModel = model1)

# Plot diagnostics
plot(sim_res)
#QQplot shows good fit

# Overdispersion test using DHARMa
testDispersion(sim_res)
#no evidence of overdispersion p = 0.88

# Filter out the outliers where Number.of.pieces.of.plastics > 20
df_filtered <- df_ratio %>%
  filter(Number.of.pieces.of.plastics <= 20)

# Create the plot with the filtered data
ggplot(df_filtered, aes(x = Average.biomass.per.individual, y = Number.of.pieces.of.plastics)) +
  geom_point(alpha = 0.6, color = "black") +  # Scatter plot points
  geom_smooth(method = "lm", se = TRUE, color = "blue") +  # Linear regression trend line
  labs(
    x = "Mean Biomass per Individual (g)",
    y = "Number of MPs"
  ) +
  theme_minimal()

library(ggsci)

ggplot(df_filtered, aes(x = Average.biomass.per.individual, y = Number.of.pieces.of.plastics)) +
  geom_point(aes(color = Taxonomic.group), size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  scale_color_manual(values = c(
    "#E69F00", # orange
    "#56B4E9", # sky blue
    "#009E73", # bluish green
    "#F0E442", # yellow
    "#0072B2", # blue
    "#6A3D9A", # dark purple
    "#CC79A7"  # reddish purple
  )) +  
  labs(
    x = "Mean Biomass per Individual (g)",
    y = "Number of pieces of MPs",
    color = "Taxonomic Group"
  ) + 
  scale_y_continuous(
    breaks = seq(0, max(df_filtered$Number.of.pieces.of.plastics, na.rm = TRUE), by = 2)  # Set interval of 2 starting from 0
  ) +  
  theme_minimal() +
  theme(
    legend.position = "right",
    text = element_text(size = 16),  # Set global text size
    axis.title = element_text(size = 16),  # Axis titles
    axis.text = element_text(size = 14),  # Axis labels
    legend.text = element_text(size = 14),  # Legend text
    legend.title = element_text(size = 14),  # Legend title
    axis.line = element_line(color = "black")  # Black axis lines
  )


df_filtered <- df_filtered %>%
  filter(Trophic.level != "Omnivore")


ggplot(df_filtered, aes(x = Average.biomass.per.individual, y = Number.of.pieces.of.plastics)) +
  geom_point(aes(color = Trophic.level), size = 3, alpha = 0.7) +
  geom_smooth(
    method = "lm", 
    se = TRUE, 
    color = "black", 
    formula = y ~ x, 
    fullrange = TRUE
  ) +
  scale_color_manual(values = c(
    "#56B4E9", # sky blue
    "#6A3D9A", # dark purple
    "#E69F00"  # orange
  )) +  
  labs(
    x = "Mean Biomass per Individual (g)",
    y = "Number of pieces of MPs",
    color = "Trophic level"
  ) +
  scale_y_continuous(limits = c(0, 6), breaks = seq(0, 6, by = 2)) +  # Adjust breaks for y-axis
  theme_minimal() +
  theme(
    legend.position = "right",
    text = element_text(size = 16),  # Set global text size
    axis.title = element_text(size = 16),  # Axis titles
    axis.text = element_text(size = 14),  # Axis labels
    legend.text = element_text(size = 14),  # Legend text
    legend.title = element_text(size = 14),  # Legend title
    axis.line = element_line(color = "black")  # Black axis lines
  )


# Create the violin plot with mean as a dot and standard deviation as error bars
ggplot(df_weight, aes(x = Taxonomic.group, y = MPs.per.gram)) +
  geom_violin(fill = alpha("grey80", 0.4), color = "black") +  # Violin plot with lighter fill
  geom_point(data = df_mps_summary, aes(x = Taxonomic.group, y = mps_mean), 
             color = "blue", size = 3) +                        # Blue dot for the mean
  geom_errorbar(data = df_mps_summary, aes(x = Taxonomic.group, y = mps_mean, 
                                           ymin = mps_mean - mps_SD, 
                                           ymax = mps_mean + mps_SD),
                width = 0.2, color = "blue") +                 # Blue error bars for SD
  theme_minimal() +
  labs(x = "Taxonomic Group", y = "MPs per gram") +
  theme(
    axis.text.x = element_text(size = 12),  # Rotate x-axis labels for better readability
    axis.text.y = element_text(size = 14),  # Larger text for y-axis
    axis.title.x = element_text(size = 16),  # Larger title for x-axis
    axis.title.y = element_text(size = 16),  # Larger title for y-axis
    panel.grid = element_blank(),           # Removes default grid lines
    axis.line = element_line(color = "black")  # Adds black axis lines
  )


df_weight <- df_filtered %>%
  filter(Trophic.level != "Omnivore")

summary(df_weight$MPs.per.gram)

df_mps_summary <- df_weight %>%
  group_by(Trophic.level) %>%
  summarise(
    mps_mean = mean(MPs.per.gram, na.rm = TRUE),
    mps_SD = sd(MPs.per.gram, na.rm = TRUE)
  )

# Create the violin plot with mean as a dot and standard deviation as error bars
ggplot(df_weight, aes(x = Trophic.level, y = MPs.per.gram)) +
  geom_violin(fill = alpha("grey80", 0.4), color = "black") +  # Violin plot with lighter fill
  geom_point(data = df_mps_summary, aes(x = Trophic.level, y = mps_mean), 
             color = "blue", size = 3) +                        # Blue dot for the mean
  geom_errorbar(data = df_mps_summary, aes(x = Trophic.level, y = mps_mean, 
                                           ymin = mps_mean - mps_SD, 
                                           ymax = mps_mean + mps_SD),
                width = 0.2, color = "blue") +                 # Blue error bars for SD
  theme_minimal() +
  labs(x = "Taxonomic Group", y = "MPs per gram") +
  theme(
    axis.text.x = element_text(size = 12),  # Rotate x-axis labels for better readability
    axis.text.y = element_text(size = 14),  # Larger text for y-axis
    axis.title.x = element_text(size = 16),  # Larger title for x-axis
    axis.title.y = element_text(size = 16),  # Larger title for y-axis
    panel.grid = element_blank(),           # Removes default grid lines
    axis.line = element_line(color = "black")  # Adds black axis lines
  )
