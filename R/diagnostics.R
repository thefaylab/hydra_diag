
#  Diagnostic plots for Hydra
#  Last update 2022-08-12
#  Maria Cristina Perez

#rm(list = ls())  #Restart R session instead to properly clean your workspace

#change directory  - use Project  (hydra_diag.Rproj) to set working directory rather than hard-coding
#setwd("C:/Users/macristina.perez/Documents/GitHub/hydra_diag/Diagnostics")
library(tidyverse)
library(ggforce)


source("R/read.report.R")

# input & output files to use (these should be arguments in the function that calls this script)
#hydra report file
repfile <- "test-data/hydra_sim.rep"
#repfile <- "../hydra_sim/hydra_sim.rep"
#time series data file
tsfile <- "test-data/hydra_sim_NOBA-ts.dat"

#datafile object
load("test-data/hydraDataList.rda")
names(hydraDataList)

# read in output files
output<-reptoRlist(repfile)# read hydra rep file

#### READ OBSERVED (.dat) AND ESTIMATED (.rep) SURVEY VALUES ####
#obs_survey<-read.table(tsfile, skip=1695, nrows=1680, header=F)
obs_survey <- hydraDataList$observedSurvSize %>% tibble()

#colnames(obs_survey) <- c('number','year','species','type','InpN','1','2','3','4','5')  #todo - modify so number of bins is based on size of dataframe
obs_survey <- obs_survey %>% pivot_longer(cols=6:ncol(.), names_to = "lenbin") %>% #filter(value != -999)%>%
    mutate(lenbin = as.integer(str_remove(lenbin, "sizebin")),
         label = rep("survey",nrow(.)),
         species = hydraDataList$speciesList[species])
obs_survey$value[which(obs_survey$value == -999)] = 0.00001

#obs_survey$label<-rep(("survey"),each=7189)
#dim(obs_survey) #7189

pred_surv<-output$pred_survey_size
obs_survey$pred_surv<-pred_surv
nll_survey<-output$nll_survey_size
obs_survey$nll_surv<-nll_survey
obs_survey$pearson<-((obs_survey$value-obs_survey$pred_surv)/sqrt(obs_survey$pred_surv))

colnames(obs_survey) <- c('number','year','species','type','InpN','lenbin','obs_value','label',
                          'pred_value','nll','pearson')


#### READ OBSERVED (.dat) AND PREDICTED (.rep) CATCH VALUES ####

#obs_catch<-read.table(tsfile, skip=4109, nrows=608, header=F)
obs_catch <- hydraDataList$observedCatchSize %>% tibble()
#colnames(obs_catch) <- c('number','area','year','species','type','InpN','1','2','3','4','5')
obs_catch<-obs_catch %>% pivot_longer(cols=7:ncol(.), names_to = "lenbin") %>%
  mutate(lenbin = as.integer(str_remove(lenbin, "sizebin")),
         label = rep("catch",nrow(.)),
         species = hydraDataList$speciesList[species])# %>%  filter(value != -999)
obs_catch$value[which(obs_catch$value == -999)] = 0.00001

#obs_catch$label<-rep(("catch"),each=1855)
#dim(obs_catch) #1855

pred_catch<-output$pred_catch_size
obs_catch$pred_catch<-pred_catch
nll_catch<-output$nll_catch_size
obs_catch$nll_catch<-nll_catch
obs_catch$pearson<-((obs_catch$value-obs_catch$pred_catch)/sqrt(obs_catch$pred_catch))

colnames(obs_catch) <- c('number','area','year','species','type','InpN','lenbin','obs_value','label',
                         'pred_value','nll','pearson')


#diet_catch<-full_join(obs_catch, obs_survey, by = NULL,copy = FALSE)
diet_catch <- bind_rows(obs_catch, obs_survey)

#### READ OBSERVED (.dat) AND PREDICTED (.rep) DIET PROPORTION VALUES ####
#obs_diet<-read.table(tsfile, skip=4724, nrows=4721, header=F)
obs_diet <- hydraDataList$observedSurvDiet %>% tibble()

#obs_diet$label<-rep(("diet"),each=4721)
obs_diet<-obs_diet %>% pivot_longer(cols=6:ncol(.), names_to = "prey") %>%
  mutate(#lenbin = as.integer(str_remove(lenbin, "V")),
         species = hydraDataList$speciesList[species],
         label = rep("diet",nrow(.))) #%>% filter(value >0 )
obs_diet$value[which(obs_diet$value == -999)] = 0.00001

#obs_diet$label<-rep(("diet"),each=22810)
#dim(obs_diet) #56652

# observed data in hydradatalist 56652, predicted data in output file 56412

pred_diet<-output$pred_dietprop
obs_diet$pred_diet<-pred_diet
nll_diet<-output$nll_dietprop
obs_diet$nll_diet<-nll_diet
obs_diet$pearson<-((obs_diet$value-obs_diet$pred_diet)/sqrt(obs_diet$pred_diet))
colnames(obs_diet) <- c('number','year','species','lenbin','InpN','prey','obs_value','label', 'pred_value','nll','pearson')
write.csv(pred_diet, file = "outputs/pred_diet.csv", row.names = T)


#mydata<-full_join(diet_catch, obs_diet, by = NULL,copy = FALSE)
#mydata <- bind_rows(diet_catch, obs_diet)
mydata <- diet_catch %>%
  mutate(#prey = replace_na(prey, -99),
         area = replace_na(area, 1),
         type = replace_na(type, 0))
# mydata$prey[is.na(mydata$prey)] <- -99 # remove 99?s
# mydata$area[is.na(mydata$area)] <- 1
# mydata$type[is.na(mydata$type)] <- 0

# save data frame with observed and predicted values
write.csv(mydata, file = "outputs/sizecomps.csv", row.names = T)
write.csv(obs_diet, file = "outputs/survdiet.csv", row.names = T)

#### PLOTS ####
data <- read.csv("outputs/sizecomps.csv", header = T)
data<- select(data, -X)
data$residual<-ifelse(data$pearson<0,"negative","positive")
data$res_abs<-abs(data$pearson)
data<-as.data.frame(data)

#select type of data "label= catch, survey or diet"
temp.catch = data[which(data$label == "catch"),]
temp.surv = data[which(data$label == "survey"),]
#temp.diet = data[which(data$label == "diet"),]

data<- read.csv("outputs/survdiet.csv", header = T)
data<- select(data, -X)
data$residual<-ifelse(data$pearson<0,"negative","positive")
data$res_abs<-abs(data$pearson)
temp.diet = data[which(data$label == "diet"),]


#### plot 1 length composition plots by species (catch) ####

long_catch <- temp.catch %>%
  pivot_longer(cols = c("pred_value","obs_value"),
               names_to = c("kind","junk"),names_sep = "_") %>%
  select(-junk)

#sp<-1
plotdir <- "outputs/figures/comp_plots/by_species/"

especies<-unique(long_catch$species)
for (sp in especies) {

  temp_size<-long_catch%>% filter(species == sp) %>%
    group_by(year) %>%
    summarize(mu_ss=mean(InpN))

  plot_catch<- long_catch %>% filter (species==sp) %>%
    ggplot() +
    aes(x=lenbin, y = value) +
    geom_line(aes(col = kind)) +
    facet_wrap(~year, dir="v") +
    geom_text(data=temp_size, aes(x = 4.8, y = 0.5, label = mu_ss), size=3) +
    theme(legend.position = "bottom") +
    labs(col="") +
    guides(col = guide_legend(nrow = 1))

    ggsave(paste0(plotdir,"complot_catch_",sp,".jpeg"), plot_catch, width = 10, height = 7, dpi = 300)#, width=3000, height=2500, res=250)

}

#### plot 1 length composition plots by species (survey) ####

long_surv <- temp.surv %>%
  pivot_longer(cols = c("pred_value","obs_value"),
               names_to = c("kind","junk"),names_sep = "_") %>%
  select(-junk)

sp<-1

especies<-unique(long_surv$species)
for (sp in especies) {

  temp_size<-long_surv %>% filter(species == sp & number==1) %>%
    group_by(year) %>%
    summarize(mu_ss=mean(InpN))

  plot_surv<- long_surv %>% filter (species==sp& number==1) %>%
    ggplot() +
    aes(x=lenbin, y = value) +
    geom_line(aes(col = kind)) +
    facet_wrap(~year, dir="v") +
    geom_text(data=temp_size, aes(x = 4.5, y = 0.5, label = mu_ss), size=3) +
    theme(legend.position = "bottom") +
    labs(col="") +
    guides(col = guide_legend(nrow = 1))

  ggsave(paste0(plotdir,"complot_surv_",sp,".jpeg"), plot_surv, width = 10, height = 7, dpi = 300)#, width=3000, height=2500, res=250)

}

#### plot 2 length-comp plots by species aggregated by year (catch and survey) ####

plotdir <- "outputs/figures/comp_plots/by_year/"

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
  tt2 = data.frame(ESPECIE = especie[e], BIN = seq(1, 5, 1), PROP_OBS = prop_new_obs, PROP_EST = prop_new_est)

  temporal2 = rbind(temporal2, tt2) # almacena por cada especie
}


long_temp2 <- temporal2 %>%
  pivot_longer(cols = c("PROP_OBS","PROP_EST"),
               names_to = c("kind","junk"),names_sep = " ") %>%
  select(-junk)

complot_year<- long_temp2 %>%
  ggplot() +
  aes(x = BIN, y = value) +
  geom_line(aes(col=kind)) +
  theme(title = element_text(angle = 0, hjust = 0.5, size=15, colour="black")) +
  labs(x="length bin", y="proportion value", title="Length composition aggregated by year catch data") +
  facet_wrap(.~ESPECIE) +
  theme(legend.position = "bottom") +
  labs(col="") +
  guides(col = guide_legend(nrow = 1))

ggsave(paste0(plotdir,"complot_year_catch.jpeg"), complot_year)#, width=3000, height=2500, res=250)


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
  tt2 = data.frame(ESPECIE = especie[e], BIN = seq(1, 5, 1), PROP_OBS = prop_new_obs, PROP_EST = prop_new_est)

  temporal2 = rbind(temporal2, tt2) # almacena por cada especie
}

long_temp2 <- temporal2 %>%
  pivot_longer(cols = c("PROP_OBS","PROP_EST"),
               names_to = c("kind","junk"),names_sep = " ") %>%
  select(-junk)

complot_year<- long_temp2 %>%
  ggplot() +
  aes(x = BIN, y = value) +
  geom_line(aes(col=kind)) +
  theme(title = element_text(angle = 0, hjust = 0.5, size=15, colour="black")) +
  labs(x="length bin", y="proportion value", title="Length composition aggregated by year survey data") +
  facet_wrap(.~ESPECIE) +
  theme(legend.position = "bottom") +
  labs(col="") +
  guides(col = guide_legend(nrow = 1))

ggsave(paste0(plotdir,"complot_year_survey.jpeg"), complot_year)#, width=3000, height=2500, res=250)


#### plot 3 pearson residuals bubble plot ####
plotdir <- "outputs/figures/residuals/"

### catch, species= 1 to 11

tiff(paste0(plotdir,"bubbleplot_catch.jpeg"), width=3000, height=2500, res=250)
ggplot(temp.catch, aes(x=year, y=lenbin, size = res_abs, color=factor(residual))) +
  geom_point(alpha=0.7) + theme(title = element_text(angle = 0, hjust = 0.5, size=15, colour="black")) +
  facet_wrap(.~species) + labs(x="year", y="length bin", title="Pearson residuals")
dev.off()

### survey, species = 1 to 11

tiff(paste0(plotdir,"bubbleplot_survey.jpeg"), width=3000, height=2500, res=250)
ggplot(temp.surv, aes(x=year, y=lenbin, size = res_abs, color=factor(residual))) +
  geom_point(alpha=0.7) + theme(title = element_text(angle = 0, hjust = 0.5, size=15, colour="black")) +
  facet_wrap(.~species) + labs(x="year", y="length bin", title="Pearson residuals")
dev.off()


#### plot 4 proportions of stomach weights for each prey given a predator size class ####
plotdir <- "outputs/figures/diet/"

predicted<- as.data.frame(cbind(temp.diet$number, temp.diet$year, temp.diet$species,
                                temp.diet$lenbin, temp.diet$pred_value, temp.diet$prey))
colnames(predicted) <- c('number','year','species','lenbin', 'prop', 'prey')
predicted$type2<-"e" #rep(("e"),each=22810)
predicted$sizefit<- paste0(predicted$lenbin,".",predicted$type2)

observed<- as.data.frame(cbind(temp.diet$number, temp.diet$year, temp.diet$species,
                               temp.diet$lenbin, temp.diet$obs_value, temp.diet$prey))
colnames(observed) <- c('number','year','species','lenbin', 'prop', 'prey')
observed$type2<-"o" #rep(("o"),each=22810)
observed$sizefit<- paste0(observed$lenbin,".",observed$type2)

pred_obs <- bind_rows(predicted, observed)# %>%
  #mutate(sizefit = paste0(lenbin,".",type2))


### species = 1

sp<-1
nsize <- hydraDataList$Nsizebins
stringbit <- paste0(rep(1:nsize, each=2),c(".o",".e"))
limits_use <- rep("",3*nsize)
breaks_use <- rep(NA,3*nsize)
for (i in 1:nsize) {
  lo <- i*3-2
  hi <- i*3-1
  limits_use[lo:hi] <- paste0(rep(i, 2),c(".o",".e"))
  breaks_use[lo:hi] <- limits_use[lo:hi]
}

 paste0(rep(1:nsize, each=2),c(".o",".e"))


especies<-unique(pred_obs$species)
for (sp in especies) {

  plot_diet<-  pred_obs%>% filter(species == sp & number==1) %>%
    ggplot(aes(x = sizefit, y = prop, group = type2, fill = factor(prey))) +
    geom_col(position = "fill") +
    scale_x_discrete(limits = limits_use,
                     # c("1.o","1.e","",
                     #            "2.o","2.e","",
                     #            "3.o","3.e","",
                     #            "4.o","4.e","",
                     #            "5.o","5.e",""),
                     breaks = breaks_use, #c("1.o","1.e",NA,
                                # "2.o","2.e",NA,
                                # "3.o","3.e",NA,
                                # "4.o","4.e",NA,
                                # "5.o","5.e",NA),
                     labels = limits_use) + #c("1.o","1.e","",
                                # "2.o","2.e","",
                                # "3.o","3.e","",
                                # "4.o","4.e","",
                                # "5.o","5.e","")) +
    coord_flip() +
    facet_wrap(~year) +
    theme_bw() +
    labs(x = "size & source (o=observed, e=expected)",
         fill = "prey",
         y = "proportion in diet") +
    scale_fill_brewer(type = "qual", palette = 3)
  # facet_wrap_paginate(~ , ncol = 4, nrow = 5, page = 4)

  ggsave(paste0(plotdir,"diet_plot_",sp,".jpeg"), plot_diet, width = 10, height = 9, dpi = 300)#, width=3000, height=2500, res=250)

}

#### plot 5 survey sample size plots ####
plotdir <- "outputs/figures/"

dat_surv <- read.csv("outputs/sizecomps.csv", header = T)
dat_surv<-dat_surv %>% filter (label=="survey")

temp_surv = numeric()
number = numeric(); number = sort(unique(dat_surv$number))
for (n in 1:length(number)) {
  pos0 = numeric(); pos0 = which(dat_surv$number == number[n])

  species = numeric(); species = sort(unique(dat_surv$species[pos0]))
  for (e in 1:length(species)) {
    pos1 = numeric(); pos1 = which(dat_surv$species[pos0] == species[e])

    year = numeric(); year = sort(unique(dat_surv$year[pos0][pos1])) #pos1 hace que se trabaje por ano, seleccion de valores que cumplen esa condicion que estan guardados en el pos1
    for (y in 1:length(year)) {
      pos2 = numeric(); pos2 = which(dat_surv$year[pos0][pos1] == year[y])

      total=numeric()
      sum1=sum(dat_surv$pred_value[pos0][pos1][pos2]*(1-dat_surv$pred_value[pos0][pos1][pos2]))
      sum2=sum((dat_surv$obs_value[pos0][pos1][pos2] - dat_surv$pred_value[pos0][pos1][pos2])^2)
      total= sum1/sum2

      tt1 = numeric()
      tt1 = data.frame(number=number[n],species=species[e],Year=year[y],InpN=unique(dat_surv$InpN[pos0][pos1][pos2]), EffN=total*unique(dat_surv$InpN[pos0][pos1][pos2]))
      temp_surv = rbind(temp_surv, tt1)
    }
  }
}

#write.csv(temp_surv, file = "outputs/ss.csv", row.names = T)

#sp<-1

for (isurvey in unique(temp_surv$number)) {

  plot_size <- temp_surv %>% filter(number==isurvey) %>%
    ggplot() +
    aes(x=EffN, y = InpN) +
    geom_point() +
    facet_wrap(~species)

  ggsave(paste0(plotdir,"samplesize_survey_",isurvey,".jpeg"), plot_size, width = 10, height = 7, dpi = 300)#, width=3000, height=2500, res=250)

}


#### plot 5 catch sample size plots ####

dat_catch <- read.csv("outputs/sizecomps.csv", header = T)
dat_catch<-dat_catch %>% filter (label=="catch")

temp_catch = numeric()
number = numeric(); number = sort(unique(dat_catch$number))
for (n in 1:length(number)) {
  pos0 = numeric(); pos0 = which(dat_catch$number == number[n])

  species = numeric(); species = sort(unique(dat_catch$species[pos0]))
  for (e in 1:length(species)) {
    pos1 = numeric(); pos1 = which(dat_catch$species[pos0] == species[e])

    year = numeric(); year = sort(unique(dat_catch$year[pos0][pos1])) #pos1 hace que se trabaje por ano, seleccion de valores que cumplen esa condicion que estan guardados en el pos1
    for (y in 1:length(year)) {
      pos2 = numeric(); pos2 = which(dat_catch$year[pos0][pos1] == year[y])

      total=numeric()
      sum1=sum(dat_catch$pred_value[pos0][pos1][pos2]*(1-dat_catch$pred_value[pos0][pos1][pos2]))
      sum2=sum((dat_catch$obs_value[pos0][pos1][pos2] - dat_catch$pred_value[pos0][pos1][pos2])^2)
      total= sum1/sum2

      tt1 = numeric()
      tt1 = data.frame(number=number[n],species=species[e],Year=year[y],InpN=unique(dat_catch$InpN[pos0][pos1][pos2]), EffN=total*unique(dat_catch$InpN[pos0][pos1][pos2]))
      temp_catch = rbind(temp_catch, tt1)
    }
  }
}

temp_catch
#write.csv(temp_catch, file = "outputs/cc.csv", row.names = T)

sp<-1

for (ifleet in unique(temp_catch[,1])) {

  plot_size<- temp_catch %>% filter(number==ifleet) %>%
    ggplot() +
    aes(x=EffN, y = InpN) +
    geom_point() +
    facet_wrap(~species)

  ggsave(paste0(plotdir,"samplesize_catch_",ifleet,".jpeg"), plot_size, width = 10, height = 7, dpi = 300)#, width=3000, height=2500, res=250)

}


#### END ####



