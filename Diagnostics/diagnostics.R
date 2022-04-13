rm(list = ls())
setwd("C:/Users/macristina.perez/Documents/UMassD/2021/Tesis/hydra_sim/Diagnostics")
library(tidyverse)
source("read.report.R")
output<-reptoRlist("hydra_sim.rep")
names(output)
head(output)

#Ncatch_obs 720
#Ncatch_size_obs 608
#Nsurvey_size_obs 1680
#Ndietprop_obs 4721
#catch_size_pred 1855

######### READ OBSERVED AND ESTIMATED SURVEY VALUES ###################

obs_survey<-read.table("hydra_sim_NOBA-ts.dat", skip=1695, nrows=1680, header=F)
colnames(obs_survey) <- c('number','year','species','type','InpN','1','2','3','4','5')
obs_survey<-obs_survey %>% pivot_longer(cols=6:10, names_to = "lenbin") %>%
  filter(value != -999)%>%         
  mutate(lenbin = as.integer(str_remove(lenbin, "V"))) 
obs_survey$label<-rep(("survey"),each=7189)
#dim(obs_survey) #7189

pred_surv<-output$pred_survey_size
obs_survey$pred_surv<-pred_surv
nll_survey<-output$nll_survey_size
obs_survey$nll_surv<-nll_survey

colnames(obs_survey) <- c('number','year','species','type','InpN','lenbin','obs_value','label',
                          'pred_value','nll')


######### READ OBSERVED AND PREDICTED CATCH VALUES ###################

obs_catch<-read.table("hydra_sim_NOBA-ts.dat", skip=4109, nrows=608, header=F)
colnames(obs_catch) <- c('number','area','year','species','type','InpN','1','2','3','4','5')
obs_catch<-obs_catch %>% pivot_longer(cols=7:11, names_to = "lenbin") %>%
  mutate(lenbin = as.integer(str_remove(lenbin, "V"))) %>%
  filter(value != -999) 
obs_catch$label<-rep(("catch"),each=1855)
#dim(obs_catch) #1855

pred_catch<-output$pred_catch_size
obs_catch$pred_catch<-pred_catch
nll_catch<-output$nll_catch_size
obs_catch$nll_catch<-nll_catch

colnames(obs_catch) <- c('number','area','year','species','type','InpN','lenbin','obs_value','label',
                         'pred_value','nll')


mydata<-full_join(obs_catch, obs_survey, by = NULL,copy = FALSE)

#### OK 


######### READ OBSERVED AND PREDICTED DIET PROPORTION VALUES ###################

obs_diet<-read.table("hydra_sim_NOBA-ts.dat", skip=4724, nrows=4721, header=F)
obs_diet$label<-rep(("diet"),each=4721)
dim(obs_diet) #4721

#I need to move the sprecies prop (wt_prey_1	wt_prey_2.....) down (as rows) and not to the side (as columns)


####### pred values 
pred_diet<-output$pred_dietprop
obs_diet$pred_diet<-pred_diet
nll_diet<-output$nll_dietprop
obs_diet$nll_diet<-nll_diet


colnames(obs_catch) <- c('number','area','year','species','type','InpN','lenbin','obs_value','label',
                         'pred_value','nll')



################# CREATE THE TABLE ######################################\

obs_values<-full_join(obs_catch, obs_survey, by = NULL,copy = FALSE)
#write.csv(obs_values, file = "observed_values.csv", row.names = T)






