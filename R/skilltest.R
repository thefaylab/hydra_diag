# where are the results
d.name <- "~/Documents/0_Data/MSmods/hydra-est/hydra_sim/results/"

# where is the data input
data_object = "~/Documents/0_Data/MSmods/hydra-est/hydra_sim/inputRdatalists/hydraDataList_5bin_0comp.rds"
hydraDataList <- readRDS(data_object)

# 5 bin biomass output
list_05 <- c(
  "05bin",
  "NOBA05bin_5k",
  "NOBA05bin_5k2",
  "ignore/05bin_old",
  "ignore/NOBA05bin",
  "ignore/NOBA05bin_2",
  "ignore/NOBA05bin_5k"
)

# 10 bin biomass output
list_10 <- c(
  "10bin",
  "NOBA10bin_5k",
  "NOBA10bin_5k2",
  "ignore/10bin_old",
  "ignore/NOBA10bin",
  "ignore/NOBA10bin_2",
  "ignore/NOBA10bin_5k"
)

# list of repfiles
repfiles<-list.files(path=paste0(d.name, list_05), pattern="\\.rep$", full.names = TRUE)
repfiles<-list.files(path=paste0(d.name, list_10), pattern="\\.rep$", full.names = TRUE)

# read rep files and save est_bio
source("../R/read.report.R")

output<-purrr::map(repfiles, reptoRlist)
est_bio <- map(output, "EstBsize") %>%
  #walk(as.data.frame) %>%
  map_dfr(as.data.frame, .id = "model") %>%
  I()

nmodel <- length(repfiles)

stepperyr <- nrow(est_bio)/length(repfiles)/hydraDataList$Nyrs/length(hydraDataList$speciesList)

est_bio <- est_bio %>%  #output$EstBsize %>%
  #rowSums() %>%
  tibble() %>%
  mutate(bio = rowSums(across(where(is.numeric)))) %>%
  mutate(species = rep(rep(hydraDataList$speciesList, each = hydraDataList$Nyrs*stepperyr),nmodel),
         year  = rep(rep(1:(hydraDataList$Nyrs*stepperyr),hydraDataList$Nspecies),nmodel),
         year = (1-1/stepperyr) + year / stepperyr,  #5 time steps per year
         log_bio = ifelse(bio>0,log(bio),NA)) %>%
    #model = as.factor(model)) %>%
  I()

# true bio for each
truthfile = "~/Documents/0_Data/ms-keyrun/simulated-data/atlantisoutput/NOBA_sacc_38/nordic_runresults_01omlist_ss.rds"

# read true object
omlist_ss <- readRDS(truthfile)

#Number of years
nyears <- omlist_ss$runpar$nyears
total_sample <- omlist_ss$runpar$tstop/omlist_ss$runpar$outputstep

# hardcode this that only Poseidon knows
# or anyone who looks at mskeyrun/data-raw/build_simdata.R
fitstart = 40
fitend = 120

#survey species inherited from omlist_ss
atlsurvspp <- omlist_ss$species_ss
# survey season and other time dimensioning parameters
# generalized timesteps all models
atlnoutsteps <- omlist_ss$runpar$tstop/omlist_ss$runpar$outputstep
atltimeall <- c(0:atlnoutsteps)
atlstepperyr <- if(omlist_ss$runpar$outputstepunit=="days") 365/omlist_ss$runpar$toutinc


# user specified fit start and times if different from full run
fitstartyr <- ifelse(!is.null(fitstart), fitstart-1, 0)
fitendyr <- ifelse(!is.null(fitend), fitend, total_sample)

atlantis_full <- c(1:total_sample)
mod_burnin <- fitstartyr*atlstepperyr+1
fit_nyears <- fitendyr-fitstartyr
fit_ntimes <- fit_nyears*atlstepperyr
fittimes <- atlantis_full[mod_burnin:(mod_burnin+fit_ntimes-1)]
#fit_timesteps <- seq(fittimes[stepperyr], max(fittimes), by=stepperyr) #last timestep
#fit_years <- unique(floor(fittimes/stepperyr)) #from Christine's new sardine_config.R
fittimes.days <- if(omlist_ss$runpar$outputstepunit=="days") fittimes*omlist_ss$runpar$outputstep


truebio <- omlist_ss$truetotbio_ss %>%
  dplyr::filter(time %in% fittimes.days) %>%
  dplyr::mutate(year = (time/365),
                year = year-fitstartyr,
                biomass = atoutput,
                logbio = log(biomass)) %>%
  I()

# compare time series on one plot to truth (natural)
plotB <-ggplot() +
  geom_line(data=truebio, aes(x=year,y=biomass),
            alpha = 10/10) +
  geom_line(data=est_bio, aes(x=year,y=bio, colour=model),
             alpha = 5/10) +
  theme_minimal() +
  theme(legend.position = "top") +
  labs(x = "Year",
       y = "Biomass (t)",
       title = "Total biomass skill, 10 bin Hydra estimated (color) vs True Atlantis (black)",
       colour="")

# compare time series on one plot to truth (log)
plotB <-ggplot() +
  geom_line(data=truebio, aes(x=year,y=logbio),
            alpha = 10/10) +
  geom_line(data=est_bio, aes(x=year,y=log_bio, colour=model),
            alpha = 5/10) +
  theme_minimal() +
  theme(legend.position = "top") +
  labs(x = "Year",
       y = "Log Biomass (t)",
       title = "Total biomass skill, 5 bin Hydra estimated (color) vs True Atlantis (black)",
       colour="")


plotB +
  facet_wrap(~species, scales="free")

# calculate metrics for each
#skill assessment metrics to evaluate model fit
# following from http://heuristically.wordpress.com/2013/07/12/calculate-rmse-and-mae-in-r-and-sas/
# Function that returns Root Mean Squared Error
rmse <- function(error)
{
  sqrt(mean(error^2, na.rm=T))
}

# Function that returns Mean Absolute Error (SAME AS AAE FOR US)
aae <- function(error)
{
  mean(abs(error), na.rm=T)
}

# Function that returns Mean Absolute Error (SAME AS AE FOR US)
ae <- function(error)
{
  mean(error, na.rm=T)
}

# Function that returns Modeling Efficiency MEF
mef <- function(obs, error)
{
  obsavg <- mean(obs, na.rm=T)
  obserr <- obs - obsavg
  (sum(obserr^2)-sum(error^2, na.rm=T))/sum(obserr^2)
}

# for 5 timestep model only (5bin)
# need to match correctly for other outputs
# Atlantis happens to have 5 output steps per year
trueb <- as.tibble(truebio) %>%
  dplyr::mutate(year = as.character(year)) %>%
  dplyr::select(species, year, biomass, logbio)

est_bio_err <- est_bio %>%
  dplyr::mutate(year = as.character(year)) %>%
  dplyr::left_join(trueb) %>%
  dplyr::select(model, species, year, estbio = bio, estlogB = log_bio, truebio = biomass, truelogB = logbio) %>%
  dplyr::mutate(error = estbio - truebio,
                lerror = estlogB - truelogB)

skillsp <- est_bio_err %>%
  dplyr::group_by(model, species) %>%
  dplyr::summarise(RMSE = rmse(error),
                   AAE = aae(error),
                   AE = ae(error),
                   MEF = mef(estbio, error)) %>%
  tidyr::pivot_longer(-c(model, species), names_to = "metric", values_to = "value")

skillmod <- est_bio_err %>%
  dplyr::group_by(model) %>%
  dplyr::summarise(RMSE = rmse(error),
                   AAE = aae(error),
                   AE = ae(error),
                   MEF = mef(estbio, error))

plotskill <- ggplot(skillspl) +
  geom_boxplot(aes(x = species, y = value))+ #, colour=model))+
  theme_bw() +
  theme(axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      size = 12)) +
  facet_wrap(~metric, scales = "free_y")

skillspl <- est_bio_err %>%
  dplyr::group_by(model, species) %>%
  dplyr::summarise(RMSE = rmse(lerror),
                   AAE = aae(lerror),
                   AE = ae(lerror),
                   MEF = mef(estlogB, lerror)) %>%
  tidyr::pivot_longer(-c(model, species), names_to = "metric", values_to = "value")

plotmef <- ggplot(skillspl %>% filter(metric=="MEF")) +
  geom_hline(yintercept=0) +
  geom_jitter(aes(x = species, y = value, colour=model), width = 0.2)+
  theme_bw() +
  theme(axis.text.x = element_text(
    angle = 90,
    hjust = 1,
    size = 12))


skillmodl <- est_bio_err %>%
  dplyr::group_by(model) %>%
  dplyr::summarise(RMSE = rmse(error),
                   AAE = aae(error),
                   AE = ae(error),
                   MEF = mef(estlogB, error))
