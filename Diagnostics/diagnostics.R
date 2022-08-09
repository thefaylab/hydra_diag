
#  Diagnostic plots for Hydra
#  Last update 03/08/2022
#  Maria Cristina Perez

rm(list = ls())

#change directory
setwd("C:/Users/macristina.perez/Documents/GitHub/hydra_diag/Diagnostics")
library(tidyverse)
source("read.report.R")
output<-reptoRlist("hydra_sim.rep")# read hydra rep file
names(output)
head(output)

#Ncatch_obs 720
#Ncatch_size_obs 608
#Nsurvey_size_obs 1680
#Ndietprop_obs 4721
#catch_size_pred 1855

#### READ OBSERVED (.dat) AND ESTIMATED (.rep) SURVEY VALUES ####
# careful with "skip" and "nrows" value!

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


#### READ OBSERVED (.dat) AND PREDICTED (.rep) CATCH VALUES ####

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

#### READ OBSERVED (.dat) AND PREDICTED (.rep) DIET PROPORTION VALUES ####

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
mydata$prey[is.na(mydata$prey)] <- -99 # remove 99?s
mydata$area[is.na(mydata$area)] <- 1
mydata$type[is.na(mydata$type)] <- 0

# save data frame with observed and predicted values
# write.csv(mydata, file = "mydata.csv", row.names = T)


#### PLOTS ####
# 2 options --> read the csv file saved in the previous steps
#           --> or continue with mydata data frame

rm(list = ls())
setwd("C:/Users/macristina.perez/Documents/GitHub/hydra_diag/Diagnostics")
data<-read.csv("Diagnostics/mydata.csv", header = T)
#data<- select(data, -X)
data$residual<-ifelse(data$pearson<0,"negative","positive")
data$res_abs<-abs(data$pearson)
data<-as.data.frame(data)
head(data)

library(tidyverse)
library(ggplot2)
names(data)
str(data)

#select type of data "label= catch, survey or diet"
temp.catch = data[which(data$label == "catch"),]
temp.surv = data[which(data$label == "survey"),]
temp.diet = data[which(data$label == "diet"),]


#### plot 1 length composition plots by species (catch) ####

long_data <- temp.catch %>%
  pivot_longer(cols = c("pred_value","obs_value"),
               names_to = c("kind","junk"),names_sep = "_") %>%
  select(-junk)

#write.csv(long_data, file = "long_data.csv", row.names = T)


sp<-1

especies<-unique(long_data$species)
for (sp in especies) {

  temp_size<-long_data%>% filter(species == sp) %>%
    group_by(year) %>%
    summarize(mu_ss=mean(InpN))

  plot_catch<- long_data %>% filter (long_data$species==sp) %>%
    ggplot() +
    #aes(x=lenbin, y = value, col = kind) +
    geom_line(aes(x=lenbin, y = value, col = kind)) +
    facet_wrap(~year, dir="v") +
    geom_text(data=temp_size, aes(x = 4.8, y = 0.5, label = mu_ss), size=3) +# glue::glue("n={N}")), size=3)
    theme(legend.position = "bottom") +
    labs(col="") +
    guides(col = guide_legend(nrow = 1))

  ggsave(paste0("complot_catch_",sp,".jpeg"), plot_catch, width = 10, height = 7, dpi = 300)#, width=3000, height=2500, res=250)

}








#######
sp<-1

especies<-unique(temp.catch$species)
  for (sp in especies) {

temp_size<-temp.catch%>% filter(species == sp) %>%
  group_by(year) %>%
  summarize(mu_ss=mean(InpN))

plot_catch<-  temp.catch%>% filter(species == sp) %>%
    ggplot(aes(x = lenbin, y = obs_value), ylim=c(0,0.8)) +
    geom_line() + theme(title = element_text(angle = 0, hjust = 0.5, size=15, colour="black")) +
    geom_point() + labs(x="length bin", y="proportion value", title="Catch length comp by year") +
    geom_line(aes(x = lenbin, y = pred_value), color = "green") +
    facet_wrap(.~year, dir="v") +
    geom_text(data=temp_size, aes(x = 4.8, y = 0.5, label = mu_ss), size=3)# glue::glue("n={N}")), size=3)

ggsave(paste0("complot_catch_",sp,".jpeg"), plot_catch)#, width=3000, height=2500, res=250)

}

#### plot 1 length composition plots by species (survey) ####

sp<-1

especies<-unique(temp.surv$species)
for (sp in especies) {

  temp_size<-temp.surv %>% filter(species == sp & number==1) %>%
    group_by(year) %>%
    summarize(mu_ss=mean(InpN))

  plot_surv<-  temp.surv%>% filter(species == sp & number==1) %>%
    ggplot(aes(x = lenbin, y = obs_value), ylim=c(0,0.8)) +
    geom_line() + theme(title = element_text(angle = 0, hjust = 0.5, size=15, colour="black")) +
    geom_point() + labs(x="length bin", y="proportion value", title="Survey length comp by year") +
    geom_line(aes(x = lenbin, y = pred_value), color = "green") +
    facet_wrap(.~year, dir="v") +
    geom_text(data=temp_size, aes(x = 4.8, y = 0.5, label = mu_ss), size=3)# glue::glue("n={N}")), size=3)

  ggsave(paste0("complot_surv_",sp,".jpeg"), plot_surv)#, width=3000, height=2500, res=250)

}

#### plot 2 length-comp plots by species aggregated by year (catch and survey) ####

### catch

temporal2 = numeric()
especie = numeric(); especie = sort(unique(temp.catch$species)) # especies
for(e in 1:length(especie)){
  pos1 = numeric(); pos1 = which(temp.catch$species == especie[e])

  temporal1 = numeric()
  year = numeric(); year = sort(unique(temp.catch$year[pos1])) # ano
  for(y in 1:length(year)){
    pos2 = numeric(); pos2 = which(temp.catch$year[pos1] == year[y])

    N = numeric(); N = unique(temp.catch$InpN[pos1][pos2]) # sample size para el estrato especie ano
    n_ejemplares_obs = numeric(); n_ejemplares_obs = round(N * temp.catch$obs_value[pos1][pos2] , 0) # crea vector de proporcion por el sample size
    n_ejemplares_est = numeric(); n_ejemplares_est = round(N * temp.catch$pred_value[pos1][pos2], 0)

    tt1 = numeric()
    tt1 = data.frame(YEAR = year[y], BIN = temp.catch$lenbin[pos1][pos2], N = N, N_EJEM_OBS = n_ejemplares_obs, N_EJEM_EST = n_ejemplares_est)

    temporal1 = rbind(temporal1, tt1) # guarda la proporcion * sample size para la especie seleccionada y ano
  }

  bin1_obs = numeric(); bin1_obs = sum(temporal1$N_EJEM_OBS[which(temporal1$BIN == 1)])
  bin2_obs = numeric(); bin2_obs = sum(temporal1$N_EJEM_OBS[which(temporal1$BIN == 2)])
  bin3_obs = numeric(); bin3_obs = sum(temporal1$N_EJEM_OBS[which(temporal1$BIN == 3)])
  bin4_obs = numeric(); bin4_obs = sum(temporal1$N_EJEM_OBS[which(temporal1$BIN == 4)])
  bin5_obs = numeric(); bin5_obs = sum(temporal1$N_EJEM_OBS[which(temporal1$BIN == 5)])

  total_obs = sum(bin1_obs, bin2_obs, bin3_obs, bin4_obs, bin5_obs)
  prop_new_obs = numeric(); prop_new_obs = c(bin1_obs, bin2_obs, bin3_obs, bin4_obs, bin5_obs) / total_obs

  bin1_pred = numeric(); bin1_pred = sum(temporal1$N_EJEM_EST[which(temporal1$BIN == 1)])
  bin2_pred = numeric(); bin2_pred = sum(temporal1$N_EJEM_EST[which(temporal1$BIN == 2)])
  bin3_pred = numeric(); bin3_pred = sum(temporal1$N_EJEM_EST[which(temporal1$BIN == 3)])
  bin4_pred = numeric(); bin4_pred = sum(temporal1$N_EJEM_EST[which(temporal1$BIN == 4)])
  bin5_pred = numeric(); bin5_pred = sum(temporal1$N_EJEM_EST[which(temporal1$BIN == 5)])

  total_pred = sum(bin1_pred, bin2_pred, bin3_pred, bin4_pred, bin5_pred)
  prop_new_est = numeric(); prop_new_est = c(bin1_pred, bin2_pred, bin3_pred, bin4_pred, bin5_pred) / total_pred

  tt2 = numeric()
  tt2 = data.frame(ESPECIE = especie[e], BIN = seq(1, 5, 1), PROP_NEW_OBS = prop_new_obs, PROP_NEW_EST = prop_new_est)

  temporal2 = rbind(temporal2, tt2) # almacena por cada especie
}

tiff("complot_year_catch.jpeg", width=3000, height=2500, res=250)
ggplot(temporal2, aes(x = BIN, y = PROP_NEW_OBS)) +
  geom_line() +
  geom_point() +
  geom_line(aes(x = BIN, y = PROP_NEW_EST), color = "green") +
  #  annotate("text",  x = 4.3, y = 0.9, label = "n=", size=4) +
  theme(title = element_text(angle = 0, hjust = 0.5, size=15, colour="black")) +
  labs(x="length bin", y="proportion value", title="Length composition aggregated by year catch data") +
  facet_wrap(.~ESPECIE)
dev.off()

### survey

temporal2 = numeric()
especie = numeric(); especie = sort(unique(temp.surv$species)) # especies
for(e in 1:length(especie)){
  pos1 = numeric(); pos1 = which(temp.surv$species == especie[e])

  temporal1 = numeric()
  year = numeric(); year = sort(unique(temp.surv$year[pos1])) # ano
  for(y in 1:length(year)){
    pos2 = numeric(); pos2 = which(temp.surv$year[pos1] == year[y])

    N = numeric(); N = unique(temp.surv$InpN[pos1][pos2]) # sample size para el estrato especie ano
    n_ejemplares_obs = numeric(); n_ejemplares_obs = round(N * temp.surv$obs_value[pos1][pos2] , 0) # crea vector de proporcion por el sample size
    n_ejemplares_est = numeric(); n_ejemplares_est = round(N * temp.surv$pred_value[pos1][pos2], 0)

    tt1 = numeric()
    tt1 = data.frame(YEAR = year[y], BIN = temp.surv$lenbin[pos1][pos2], N = N, N_EJEM_OBS = n_ejemplares_obs, N_EJEM_EST = n_ejemplares_est)

    temporal1 = rbind(temporal1, tt1) # guarda la proporcion * sample size para la especie seleccionada y ano
  }

  bin1_obs = numeric(); bin1_obs = sum(temporal1$N_EJEM_OBS[which(temporal1$BIN == 1)])
  bin2_obs = numeric(); bin2_obs = sum(temporal1$N_EJEM_OBS[which(temporal1$BIN == 2)])
  bin3_obs = numeric(); bin3_obs = sum(temporal1$N_EJEM_OBS[which(temporal1$BIN == 3)])
  bin4_obs = numeric(); bin4_obs = sum(temporal1$N_EJEM_OBS[which(temporal1$BIN == 4)])
  bin5_obs = numeric(); bin5_obs = sum(temporal1$N_EJEM_OBS[which(temporal1$BIN == 5)])

  total_obs = sum(bin1_obs, bin2_obs, bin3_obs, bin4_obs, bin5_obs)
  prop_new_obs = numeric(); prop_new_obs = c(bin1_obs, bin2_obs, bin3_obs, bin4_obs, bin5_obs) / total_obs

  bin1_pred = numeric(); bin1_pred = sum(temporal1$N_EJEM_EST[which(temporal1$BIN == 1)])
  bin2_pred = numeric(); bin2_pred = sum(temporal1$N_EJEM_EST[which(temporal1$BIN == 2)])
  bin3_pred = numeric(); bin3_pred = sum(temporal1$N_EJEM_EST[which(temporal1$BIN == 3)])
  bin4_pred = numeric(); bin4_pred = sum(temporal1$N_EJEM_EST[which(temporal1$BIN == 4)])
  bin5_pred = numeric(); bin5_pred = sum(temporal1$N_EJEM_EST[which(temporal1$BIN == 5)])

  total_pred = sum(bin1_pred, bin2_pred, bin3_pred, bin4_pred, bin5_pred)
  prop_new_est = numeric(); prop_new_est = c(bin1_pred, bin2_pred, bin3_pred, bin4_pred, bin5_pred) / total_pred

  tt2 = numeric()
  tt2 = data.frame(ESPECIE = especie[e], BIN = seq(1, 5, 1), PROP_NEW_OBS = prop_new_obs, PROP_NEW_EST = prop_new_est)

  temporal2 = rbind(temporal2, tt2) # almacena por cada especie
}

tiff("complot_year_survey.jpeg", width=3000, height=2500, res=250)
ggplot(temporal2, aes(x = BIN, y = PROP_NEW_OBS)) +
  geom_line() +
  geom_point() +
  geom_line(aes(x = BIN, y = PROP_NEW_EST), color = "green") +
  # annotate("text",  x = 4.3, y = 0.9, label = "n=", size=4) +
  theme(title = element_text(angle = 0, hjust = 0.5, size=15, colour="black")) +
  labs(x="length bin", y="proportion value", title="Length composition aggregated by year survey data") +
  facet_wrap(.~ESPECIE)
dev.off()


#### plot 3 pearson residuals bubble plot ####

### catch, species= 1 to 11

tiff("bubbleplot_catch.jpeg", width=3000, height=2500, res=250)
ggplot(temp.catch, aes(x=year, y=lenbin, size = res_abs, color=factor(residual))) +
  geom_point(alpha=0.7) + theme(title = element_text(angle = 0, hjust = 0.5, size=15, colour="black")) +
  facet_wrap(.~species) + labs(x="year", y="length bin", title="Pearson residuals")
dev.off()

### survey, species = 1 to 11

tiff("bubbleplot_survey.jpeg", width=3000, height=2500, res=250)
ggplot(temp.surv, aes(x=year, y=lenbin, size = res_abs, color=factor(residual))) +
  geom_point(alpha=0.7) + theme(title = element_text(angle = 0, hjust = 0.5, size=15, colour="black")) +
  facet_wrap(.~species) + labs(x="year", y="length bin", title="Pearson residuals")
dev.off()


#### plot 4 proportions of stomach weights for each prey given a predator size class ####

predicted<- as.data.frame(cbind(temp.diet$number, temp.diet$year, temp.diet$species,
                                temp.diet$lenbin, temp.diet$pred_value, temp.diet$prey))
colnames(predicted) <- c('number','year','species','lenbin', 'prop', 'prey')
predicted$type2<-rep(("e"),each=22810)
predicted$sizefit<- paste0(predicted$lenbin,".",predicted$type2)

observed<- as.data.frame(cbind(temp.diet$number, temp.diet$year, temp.diet$species,
                               temp.diet$lenbin, temp.diet$obs_value, temp.diet$prey))
colnames(observed) <- c('number','year','species','lenbin', 'prop', 'prey')
observed$type2<-rep(("o"),each=22810)
observed$sizefit<- paste0(observed$lenbin,".",observed$type2)

pred_obs <- bind_rows(predicted, observed)# %>%
  #mutate(sizefit = paste0(lenbin,".",type2))

library(ggforce)

### species = 1

sp<-1

especies<-unique(pred_obs$species)
for (sp in especies) {

  plot_diet<-  pred_obs%>% filter(species == sp & number==1) %>%
    ggplot(aes(x = sizefit, y = prop, group = type2, fill = factor(prey))) +
    geom_col(position = "fill") +
    scale_x_discrete(limits = c("1.o","1.e","",
                                "2.o","2.e","",
                                "3.o","3.e","",
                                "4.o","4.e","",
                                "5.o","5.e",""),
                     breaks = c("1.o","1.e",NA,
                                "2.o","2.e",NA,
                                "3.o","3.e",NA,
                                "4.o","4.e",NA,
                                "5.o","5.e",NA),
                     labels = c("1.o","1.e","",
                                "2.o","2.e","",
                                "3.o","3.e","",
                                "4.o","4.e","",
                                "5.o","5.e","")) +
    coord_flip() +
    facet_wrap(~year) +
    theme_bw() +
    labs(x = "size & source (o=observed, e=expected)",
         fill = "prey",
         y = "proportion in diet") +
    scale_fill_brewer(type = "qual", palette = 3)
  # facet_wrap_paginate(~ , ncol = 4, nrow = 5, page = 4)

  ggsave(paste0("diet_plot_",sp,".jpeg"), plot_diet)#, width=3000, height=2500, res=250)

}

#### plot 5 survey sample size plots ####

names(temp.surv)
temp.surv$InpN
temp.surv$nll
length(temp.surv$InpN)
length(temp.surv$nll)

plot(temp.surv$nll, temp.surv$InpN)


tiff("samplesize_plot_1.jpeg", width=3000, height=2500, res=250)
pred_obs[which(pred_obs$number == 1 & pred_obs$species== 1),] %>%
  ggplot() +
  aes(x = sizefit, y = prop, group = type2, fill = factor(prey)) +
  geom_col(position = "fill") +

  dev.off()



#### plot 5 catch sample size plots ####


names(temp.catch)
temp.catch$InpN
temp.catch$nll
length(temp.catch$InpN)
length(temp.catch$nll)

plot(temp.catch$nll, temp.catch$InpN)


tiff("samplesize_plot_1.jpeg", width=3000, height=2500, res=250)
pred_obs[which(pred_obs$number == 1 & pred_obs$species== 1),] %>%
  ggplot() +
  aes(x = sizefit, y = prop, group = type2, fill = factor(prey)) +
  geom_col(position = "fill") +

dev.off()




#### END ####






#### TO DO LIST ####
# make a folder for each type of plot
# make a loop for each plot/prey/predator
# add legend, sample size (n=x) to each plot with independet values





# trying to add each sample size per year ..... NOT WORKING
tapply(temp.catch$InpN[which(temp.catch$species == 1)], temp.catch$year[which(temp.catch$species == 1)], function(x) length(unique(x)))
N<-as.numeric(t(tapply(temp.catch$InpN[which(temp.catch$species == 1)], temp.catch$year[which(temp.catch$species == 1)], unique)))
a<- a + annotate("text",  x = 4.0, y = 0.6, label = paste("Sample size:", N, sep=""), size=3)

N<-as.numeric(t(tapply(temp.catch$InpN[which(temp.catch$species == 1)], temp.catch$year[which(temp.catch$species == 1)], unique)))

# adding some text
data_text = data.frame(x = 20, y = 3.5, species = unique(temp.catch$species), label = c("605", "", "792", "","820", "", "820", "", "825", "", "830", "", "821", "", "822", "", "815", "", "826", "", "849", "", "823", "", "831", "", "818", "", "838", "", "835", "", "810", "", "838", "", "818", "", "831", "", "825", "", "804", "", "822", "", "846", "", "829", "", "824", "", "850", "", "816", "", "800", "", "821", "", "810", "", "812", "", "806", "", "825", "", "786", "", "811", "", "809", "", "815", "", "799", "", "816", "", "813", "", "807", "", "811", "", "816", "", "817", "", "794", "", "831", "", "787", "", "815", "", "808", "", "815", "", "789", "", "803", "", "809", "", "830", "", "804", "", "807", "", "820", "", "830", "", "825", "", "823", "", "818", "", "815", "", "804", "", "815", "", "836"))


ggplot(temp.catch[which(temp.catch$species == 1),], aes(x = lenbin, y = obs_value), ylim=c(0,0.8)) +
  geom_line() + theme(title = element_text(angle = 0, hjust = 0.5, size=15, colour="black")) +
  geom_point() + labs(x="length bin", y="proportion value", title="Catch length comp by year") +
  geom_line(aes(x = lenbin, y = pred_value), color = "green") +
  facet_wrap(.~year, dir="v") +
  geom_text(data = data_text, aes(x = x, y = y, label = label))

