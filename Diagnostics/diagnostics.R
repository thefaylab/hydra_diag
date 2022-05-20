rm(list = ls())
setwd("C:/Users/macristina.perez/Documents/GitHub/hydra_diag/Diagnostics")
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
obs_survey$pearson<-((obs_survey$value-obs_survey$pred_surv)/sqrt(obs_survey$pred_surv))

colnames(obs_survey) <- c('number','year','species','type','InpN','lenbin','obs_value','label',
                          'pred_value','nll','pearson')


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
obs_catch$pearson<-((obs_catch$value-obs_catch$pred_catch)/sqrt(obs_catch$pred_catch))

colnames(obs_catch) <- c('number','area','year','species','type','InpN','lenbin','obs_value','label',
                         'pred_value','nll','pearson')


diet_catch<-full_join(obs_catch, obs_survey, by = NULL,copy = FALSE)

######### READ OBSERVED AND PREDICTED DIET PROPORTION VALUES ###################

obs_diet<-read.table("hydra_sim_NOBA-ts.dat", skip=4724, nrows=4721, header=F)
#obs_diet$label<-rep(("diet"),each=4721)
obs_diet<-obs_diet %>% pivot_longer(cols=6:17, names_to = "lenbin") %>%
  mutate(lenbin = as.integer(str_remove(lenbin, "V"))) %>%
  filter(value >0 ) 
obs_diet$label<-rep(("diet"),each=22810)
#dim(obs_diet) #22810

pred_diet<-output$pred_dietprop
obs_diet$pred_diet<-pred_diet
nll_diet<-output$nll_dietprop
obs_diet$nll_diet<-nll_diet
obs_diet$pearson<-((obs_diet$value-obs_diet$pred_diet)/sqrt(obs_diet$pred_diet))
colnames(obs_diet) <- c('number','year','species','lenbin','InpN','prey','obs_value','label', 'pred_value','nll','pearson')


mydata<-full_join(diet_catch, obs_diet, by = NULL,copy = FALSE)
mydata$prey[is.na(mydata$prey)] <- -99
mydata$area[is.na(mydata$area)] <- 1
mydata$type[is.na(mydata$type)] <- 0


#write.csv(mydata, file = "mydata.csv", row.names = T)


#########################################################################################
#remotes::install_github("r4ss/r4ss")
library(r4ss)
devtools::source_url("https://github.com/r4ss/r4ss/blob/main/R/SSplotComps.R?raw=TRUE")
#########################################################################################


rm(list = ls())
setwd("C:/Users/macristina.perez/Documents/UMassD/Classes/Tesis/hydra_sim/Diagnostics")
data<-read.csv("mydata.csv", header = T)
head(data)

# plots from 1 to 10 --> for different fleets ... change to species 
# and 21 to 24 aggregated by year

#diagnostic plots for the fits to the composition data (both size composition & diet composition).
#oldnames = c("fleet", "year", "seas", "gender", "morph", "label"),
#newnames = c("Fleet", "Yr", "Seas", "Sex", "Morph", "Label")

source("SSplotComp.R")
replist<- (data)
names(replist)

SSplotComps(replist)


library(ggplot2)

sort(unique(replist$species))  # 1  2  3  4  5  6  7  8  9 10 11

temp = replist[which(replist$label == "catch"),]

#temp[which(temp$species == 2 & temp$year == 80 & temp$lenbin == 3),]


ggplot(temp[which(temp$species == 1),], aes(x = lenbin, y = pred_value)) +
  geom_line()+
  geom_line(aes(x = lenbin, y = obs_value), color = "red")+
  facet_wrap(.~year)















