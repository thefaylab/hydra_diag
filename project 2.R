source("R/read.report.R")
source("R/gettables.R")
library(tidyverse)

# mackerel, herring, spiny dogfish and cod

### READ DATA ###

load("test-data/hydraDataList.rda")
repfile <- c("test-data/hydra_sim.rep")
#output<-reptoRlist(repfile)

### CATCH AND SURVEY BIOMASS ###
obs_surveyB <- hydraDataList$observedBiomass
obs_catchB <- hydraDataList$observedCatch

biorows <- dim(obs_surveyB)[1]
catchrows <- dim(obs_catchB)[1]

indexfits <- gettables(repfile, biorows, catchrows)

biomass<-indexfits[[1]] %>%
          mutate(species = hydraDataList$speciesList[species])

catch<-indexfits[[2]]%>%
  mutate(species = hydraDataList$speciesList[species])


##### CATCH SIMULATED DATA ######

pred_catch <- catch$pred_catch
cv_catch<-catch$cv

# create simulation object
sim_data <- NULL
isim <- 1

# replace index with simulated data
sim_list <- NULL
sim_list$index <- obs_catchB %>%
  mutate(obs = rnorm(nrow(.), pred_catch, cv_catch))
#do similar for other data tables

# store simulated object
sim_data[[isim]] <- sim_list

##### SURVEY BIOMASS SIMULATED DATA ######

pred_survey <- biomass$pred_bio
cv_survey<-biomass$cv

sim_list$survey <- obs_surveyB %>%
  mutate(obs = (rnorm(nrow(.), pred_survey, cv_survey)))
#do similar for other data tables

# store simulated object
sim_data[[isim]] <- sim_list


### CATCH AND SURVEY SIZE COMPOSITION ###

surv_size <- hydraDataList$observedSurvSize %>% tibble()
surv_size <- surv_size %>% pivot_longer(cols=6:ncol(.), names_to = "lenbin") %>% #filter(value != -999)%>%

  mutate(lenbin = as.integer(str_remove(lenbin, "sizebin")),
         label = rep("survey",nrow(.)),
         species = hydraDataList$speciesList[species])
#surv_size<- surv_size %>% filter(value != -999)
#surv_size$value[which(surv_size$value == -999)] = 0.00001

pred_surv<-output$pred_survey_size
surv_size$pred_surv<-pred_surv
nll_survey<-output$nll_survey_size
surv_size$nll_surv<-nll_survey

colnames(surv_size) <- c('number','year','species','type','InpN','lenbin','obs_value','label',
                         'pred_value','nll')


catch_size <- hydraDataList$observedCatchSize %>% tibble()
catch_size<-catch_size %>% pivot_longer(cols=7:ncol(.), names_to = "lenbin") %>%
  mutate(lenbin = as.integer(str_remove(lenbin, "sizebin")),
         label = rep("catch",nrow(.)),
         species = hydraDataList$speciesList[species])# %>% filter(value != -999)
#catch_size<- catch_size %>% filter(value != -999)
#catch_size$value[which(catch_size$value == -999)] = 0.00001

pred_catch<-output$pred_catch_size
catch_size$pred_catch<-pred_catch
nll_catch<-output$nll_catch_size
catch_size$nll_catch<-nll_catch

colnames(catch_size) <- c('number','area','year','species','type','InpN','lenbin','obs_value','label',
                          'pred_value','nll')


### DIET COMPOSITION ###

diet_comp <- hydraDataList$observedSurvDiet %>% tibble()
diet_comp<-diet_comp %>% pivot_longer(cols=6:ncol(.), names_to = "prey") %>%
  mutate(#lenbin = as.integer(str_remove(lenbin, "V")),
    species = hydraDataList$speciesList[species],
    label = rep("diet",nrow(.)))
#diet_comp<- diet_comp  %>% filter(value != -999)
#%>%
#  I()
pred_diet<-output$pred_dietprop

if (length(pred_diet)!=nrow(diet_comp)) diet_comp <- diet_comp %>% filter(value != 0)

diet_comp$pred_diet<-pred_diet
nll_diet<-output$nll_dietprop
diet_comp$nll_diet<-nll_diet

colnames(diet_comp) <- c('number','year','species','lenbin','InpN','prey','obs_value','label', 'pred_value','nll')





























