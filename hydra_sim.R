source("R/read.report.R")
source("R/gettables.R")
library(tidyverse)

### READ DATA ###

load("test-data/hydraDataList.rda")
repfile <- c("test-data/hydra_sim.rep")
output<-reptoRlist(repfile)

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
sim_list$catch <- obs_catchB %>%
  mutate(obs = rnorm(nrow(.), pred_catch, cv_catch))

sim_data[[isim]] <- sim_list


##### SURVEY BIOMASS SIMULATED DATA ######

pred_survey <- biomass$pred_bio
cv_survey<-biomass$cv

sim_list$survey <- obs_surveyB %>%
  mutate(obs = (rnorm(nrow(.), pred_survey, cv_survey)))

# store simulated object
sim_data[[isim]] <- sim_list


#### CATCH SIZE COMPOSITION SIMULATED DATA ####

catch_size <- hydraDataList$observedCatchSize %>% tibble()
catch_size<-catch_size %>% pivot_longer(cols=7:ncol(.), names_to = "lenbin") %>%
  mutate(lenbin = as.integer(str_remove(lenbin, "sizebin")),
         label = rep("catch",nrow(.)),
         species = hydraDataList$speciesList[species])# %>% filter(value != -999)
#catch_size<- catch_size %>% filter(value != -999)
catch_size$value[which(catch_size$value == -999)] = 0.000001
catch_size <- select(catch_size, -label)

pred_catchsize<-output$pred_catch_size
nll_catch<-output$nll_catch_size

inpN<-unique(catch_size$inpN)

sim_list$catchsize <- catch_size %>%
  mutate(obs = (rmultinom(nrow(.), inpN, pred_catchsize)))

# store simulated object
sim_data[[isim]] <- sim_list

write.csv(sim_list[["catchsize"]], file = "catchsize.csv", row.names = T)







#### SURVEY SIZE COMPOSITION SIMULATED DATA ####

surv_size <- hydraDataList$observedSurvSize %>% tibble()
surv_size <- surv_size %>% pivot_longer(cols=6:ncol(.), names_to = "lenbin") %>% #filter(value != -999)%>%

  mutate(lenbin = as.integer(str_remove(lenbin, "sizebin")),
         label = rep("survey",nrow(.)),
         species = hydraDataList$speciesList[species])
#surv_size<- surv_size %>% filter(value != -999)
surv_size$value[which(surv_size$value == -999)] = 0.000001
surv_size <- select(surv_size, -label)

pred_survsize<-output$pred_survey_size
nll_survey<-output$nll_survey_size

inpN<-unique(surv_size$inpN)

sim_list$survsize <- surv_size %>%
  mutate(obs = (rmultinom(nrow(.), inpN, pred_survsize)))

# store simulated object
sim_data[[isim]] <- sim_list

write.csv(sim_list[["survsize"]], file = "survsize.csv", row.names = T)



#### DIET COMPOSITION SIMULATED DATA ####

diet_comp <- hydraDataList$observedSurvDiet %>% tibble()
diet_comp<-diet_comp %>% pivot_longer(cols=6:ncol(.), names_to = "prey") %>%
  mutate(#lenbin = as.integer(str_remove(lenbin, "V")),
    species = hydraDataList$speciesList[species],
    label = rep("diet",nrow(.)))
diet_comp$value[which(diet_comp$value == -999)] = 0.000001
diet_comp <- select(diet_comp, -label)

#diet_comp<- diet_comp  %>% filter(value != -999)
#%>%
#  I()

pred_diet<-output$pred_dietprop
if (length(pred_diet)!=nrow(diet_comp)) diet_comp <- diet_comp %>% filter(value != 0)
nll_diet<-output$nll_dietprop

inpN<-unique(diet_comp$inpN)

sim_list$dietcomp <- diet_comp %>%
  mutate(obs = (rmultinom(nrow(.), inpN, pred_diet)))

# store simulated object
sim_data[[isim]] <- sim_list

write.csv(sim_list[["dietcomp"]], file = "dietcomp.csv", row.names = T)




















