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

#read observed survey size data obs_survey_biomass(i), pred_survey_index(i) 
obs_survey_size<-data.frame(output$obs_survey_size)
colnames(obs_survey_size) <- c('survey','year','species','type','InpN','x','obs_1','obs_2','obs_3','obs_4')
obs_survey_size$label<-rep(("survey"),each=1680)
#dim(survey_biom) #1680

#read predicted survey size data obs_survey_size(i) pred_survey_size(i)
pred_survey_size<-data.frame(output$survey_size)
colnames(pred_survey_size) <- c('survey','year','species','type','effN','x','pred_1','pred_2','pred_3','pred_4','pred_survey_size','nll_survey')
pred_survey_size$label<-rep(("survey"),each=1680)
#dim(pred_survey_size) #1680

#read observed catch size data obs_catch_size(i) 
obs_catch_size<-data.frame(output$obs_catch_size)
colnames(obs_catch_size) <- c('fleet','area','year','species','type','InpN','x','obs_1','obs_2','obs_3','obs_4')
obs_catch_size$label<-rep(("catch"),each=476)
#dim(catch_size) #476 ???????? 

#read predicted catch size data obs_catch_size(i) pred_catch_size(i)
pred_catch_size<-data.frame(output$catch_size)
colnames(pred_catch_size) <- c('fleet','area','year','species','type','effN','x','pred_1','pred_2','pred_3','pred_4','pred_catch_size','nll_catch')
pred_catch_size$label<-rep(("catch"),each=476)
#dim(catch_size) #476 ???????? 

#read diet prop data obs_dietprop(i) pred_dietprop(i) 
diet_prop<-data.frame(output$diet_prop)
colnames(diet_prop) <- c('survey','year','species','size_bin','InpN','wt_prey_1',	'wt_prey_2','wt_prey_3','wt_prey_4','wt_prey_5','wt_prey_6','wt_prey_7','wt_prey_8','wt_prey_9','wt_prey_10','wt_prey_11','other','pred_dietprop','nll_dietprop')
diet_prop$label<-rep(("diet"),each=4721)
#dim(diet_prop) #4721

#create table with catch size variables
Table_size <- data.frame(cbind(pred_catch_size$year,pred_catch_size$label,pred_catch_size$fleet,pred_catch_size$species, obs_catch_size$'obs_1', obs_catch_size$'obs_2', obs_catch_size$'obs_3', obs_catch_size$'obs_4',pred_catch_size$'pred_1', pred_catch_size$'pred_2', pred_catch_size$'pred_3', pred_catch_size$'pred_4',pred_catch_size$effN, pred_catch_size$pred_catch_size))
colnames(Table_size) <- c('Year','Label','Fleet/Survey','Species','obs_1','obs_2','obs_3','obs_4','pred_1','pred_2','pred_3','pred_4','effN','pred_value')

#create table with survey size variables
Table_survey <- data.frame(cbind(pred_survey_size$year,pred_survey_size$label,pred_survey_size$survey,pred_survey_size$species, obs_survey_size$'obs_1', obs_survey_size$'obs_2', obs_survey_size$'obs_3', obs_survey_size$'obs_4', pred_survey_size$'pred_1', pred_survey_size$'pred_2', pred_survey_size$'pred_3', pred_survey_size$'pred_4',pred_survey_size$effN, pred_survey_size$pred_survey_size))
colnames(Table_survey) <- c('Year','Label','Fleet/Survey','Species','obs_1','obs_2','obs_3','obs_4','pred_1','pred_2','pred_3','pred_4','effN','pred_value')

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


