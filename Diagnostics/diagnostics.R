rm(list = ls())
setwd("C:/Users/macristina.perez/Documents/UMassD/2021/Tesis/hydra_sim/Diagnostics")
library(tidyverse)

source("read.report.R")
output<-reptoRlist("hydra_sim.rep")
names(output)
head(output)

#number of catch and survey observations
#Ncatch_obs 720
#Nsurvey_obs 1680

#read survey data obs_survey_biomass(i), pred_survey_index(i) 
survey_biom<-data.frame(output$survey_biomass)
colnames(survey_biom) <- c('survey','year','species','obs_survey_biomass','cv','pred_survey_index','resid_survey','nll_survey')
#dim(survey_biom) #1680

#read survey size data obs_survey_size(i) pred_survey_size(i)
survey_size<-data.frame(output$survey_size)
colnames(survey_size) <- c('survey','year','species','type','effN','x','1','2','3','4','pred_survey_size','nll_survey')
survey_size$label<-rep(("survey"),each=1680)
#dim(survey_size) #1680

#read catch data obs_catch_biomass(i) pred_catch_biomass(i)
catch_biom<-data.frame(output$catch_data)
colnames(catch_biom) <- c('fleet','area','year','species','obs_catch_biomass','cv','pred_catch_biomass','resid_catch','nll_catch')
#dim(catch_biom) #720

#read catch size data obs_catch_size(i) pred_catch_size(i)
catch_size<-data.frame(output$catch_size)
colnames(catch_size) <- c('fleet','area','year','species','type','effN','x','1','2','3','4','pred_catch_size','nll_catch')
catch_size$label<-rep(("catch"),each=476)
#dim(catch_size) #476 ???????? 

#read diet prop data obs_dietprop(i) pred_dietprop(i) 
diet_prop<-data.frame(output$diet_prop)
colnames(diet_prop) <- c('survey','year','species','size_bin','InpN','wt_prey_1',	'wt_prey_2','wt_prey_3','wt_prey_4','wt_prey_5','wt_prey_6','wt_prey_7','wt_prey_8','wt_prey_9','wt_prey_10','wt_prey_11','other','pred_dietprop','nll_dietprop')
diet_prop$label<-rep(("diet"),each=4721)
#dim(diet_prop) #4721

#create table with catch size variables
Table_size <- data.frame(cbind(catch_size$year,catch_size$label,catch_size$fleet,catch_size$species, catch_size$'1', catch_size$'2', catch_size$'3', catch_size$'4',catch_size$effN, catch_size$pred_catch_size))
colnames(Table_size) <- c('Year','Label','Fleet/Survey','Species','X1','X2','X3','X4','effN','pred_value')

#create table with survey size variables
Table_survey <- data.frame(cbind(survey_size$year,survey_size$label,survey_size$survey,survey_size$species, survey_size$'1', survey_size$'2', survey_size$'3', survey_size$'4',survey_size$effN, survey_size$pred_survey_size))
colnames(Table_survey) <- c('Year','Label','Fleet/Survey','Species','X1','X2','X3','X4','effN','pred_value')

#not working yet
#create table with diet prop variables
#Table_dietprop <- data.frame(cbind(diet_prop$year,diet_prop$label,diet_prop$survey,diet_prop$species, diet_prop$'1', diet_prop$'2', diet_prop$'3', diet_prop$'4',diet_prop$effN, diet_prop$pred_dietprop))
#colnames(Table_size) <- c('Year','Label','Fleet/Survey','Species','X1','X2','X3','X4','effN','pred_value')

#join tables with catch and survey size variables
Table<-full_join(Table_size, Table_survey, by = NULL,copy = FALSE)

#save table in csv
#write.csv(Table, file = "Table_diagnostics", row.names = T)

###############################################################################################################################################




















pred_sizecatch<-output$pred_catch_size
length(pred_sizecatch)# 1458
catch_size<-output$nll_catch_size
#length(catch_size) 1458
pred_sizesurv<-output$pred_survey_size
#length(pred_sizesurv) 7189
surv_size<-output$nll_survey_size
#length(surv_size) 7189
pred_diet_prop<-output$pred_dietprop
#length(pred_diet_prop)22810 
diet_prop<-output$nll_dietprop
#length(diet_prop) 22810


