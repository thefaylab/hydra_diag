# Function reptoRlist cannot read the tables of fits to survey and catch indices correctly.
# Eventually fix in `hydra_sim` code report section, but for now read from already generated rep files:
# reads in same repfile
# takes number of survey bio and catch index rows which can be read in from data object dimensions
# returns list with 2 elements:
#  1. survey index observed and predicted with residual and nll
#  2. fishery catch observed and predicted with residual and nll
#
# usage example
#
# obs_surveyB <- hydraDataList$observedBiomass
# obs_catchB <- hydraDataList$observedCatch
# biorows <- dim(obs_surveyB)[1]
# catchrows <- dim(obs_catchB)[1]
#
# indexfits <- gettables(repfile, biorows, catchrows)


gettables <- function(repfile, biorows, catchrows){
  # scan repfile to the line "table of fits to survey"
  # next line is "survey biomass data, predicted, residual, nll"
  # columns are the -ts file inputs with predicted B, residual, and nll appended
  # based on https://stackoverflow.com/questions/37663246/extract-data-between-a-pattern-from-a-text-file
  reptext <- readLines(repfile)

  tablelist <- split(reptext, cumsum(grepl(pattern = 'table of fits to', reptext)))

  survlist <- as.list(tablelist$`1`[1:biorows])


  sv <- lapply(survlist[-c(1:2)], function(x) read.table(textConnection(x), header = FALSE))
  survdf <- do.call(rbind, sv)
  names(survdf) <- c("survey","year", "species", "biomass", "cv",
                     "pred_bio", "residual", "nll")

  catchlist <- as.list(tablelist$`2`[1:catchrows])

  ct <- lapply(catchlist[-c(1:2)], function(x) read.table(textConnection(x), header = FALSE))
  catchdf <- do.call(rbind, ct)
  names(catchdf) <- c("fishery", "area" ,   "year" ,   "species", "catch",  "cv",
                      "pred_catch", "residual", "nll")

  tablist <- list(survdf, catchdf)

  return(tablist)

}



