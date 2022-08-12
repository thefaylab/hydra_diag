# Test script for running hydra_sim from R, capturing output, 
# fitting other mods, etc
# S Gaichas 1 August 2014
# UPDATED for MS model testing, July 2015

#cleanup workspace, if desired
rm(list=objects()) 


#function to read cout output from hydra_sim into an R list
reptoRlist = function(fn)
{
  ifile=scan(fn,what="character",flush=T,blank.lines.skip=F,quiet=T)
  idx=sapply(as.double(ifile),is.na)
  vnam=ifile[idx] #list names
  nv=length(vnam) #number of objects
  A=list()
  ir=0
  for(i in 1:nv)
  {
    ir=match(vnam[i],ifile)
    if(i!=nv) irr=match(vnam[i+1],ifile) else irr=length(ifile)+1
    dum<-NA
    if(irr-ir==2) dum<-as.double(scan(fn,skip=ir,nlines=1,quiet=T,what=""))
    if(irr-ir>2) dum<-as.matrix(read.table(fn,skip=ir,nrow=irr-ir-1,fill=T))
    if(is.numeric(dum))#Logical test to ensure dealing with numbers
    {
      A[[ vnam[i ]]]<-dum
    }
  }
  return(A)
}

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

SpNames<-c("dogfish", "skate", "herring", "cod", "haddock", "ytail_fl", "wint_fl", "mackerel", "silverhake", "goosefish")
#SpNames<-c("dogfish", "wskate", "Aherring", "Acod", "haddock", "ytailfl", "wintfl", "Amack", "silverhake", "goosefish")
SpNames2<-c("elasmo1", "elasmo2", "pelagic1", "groundfish1", "groundfish2", "flatfish1", "flatfish2", "pelagic2", "groundfish3", "groundfish4")
spcol<-c("gray47", "brown", "green", "red", "blue", "gold3", "darkorange2", "violet", "steelblue", "black")
#spbor<-rep("white",10)

############################################################

# run simulation model, save outputs, make Kraken & SAS inputs

#change error type
for(er in c("noerr", "q", "surverr", "q_surverr")){


#number of simulation model runs to do
nSims<-115


if(Sys.info()['sysname']=="Windows"){
  setwd("C:/0_Data/Size_Struct_Ecomod/MSmod_testing/")  
  dir.create(outdir<-paste("simbio_", er, sep=""))
  dir.create(file.path(outdir, "crashed"))
  
  #run hydra_sim n times seeded consecutively: windows 
  for(n in 1:nSims){
     runmod<-paste("..\\Code\\hydra_sim.exe -ainp hydra_sim_", er, ".pin -sim ",n, sep="")  #writes csv file to current working dir
     #runmod<-paste("..\\Code\\hydra_sim.exe -sim",n)  #writes csv file to current working dir
     outfile<-paste(outdir,"\\simtest_",n,".out", sep="")
     write(shell(runmod, intern=T), outfile) #saves full output as well as csv
     krakenin<-paste("simout_",n,".csv", sep="")
     file.rename("simout.csv", file.path(outdir, krakenin))
     nantest<-paste("findstr nan ", outfile)  #find crashes
     if(shell(nantest)==0){
       #copy outfile and krakenin to crashed directory
       file.copy(outfile, file.path(outdir, "crashed"))
       #delete outfile and krakenin from working directory
       file.remove(outfile)
       file.remove(file.path(outdir, krakenin))
     }
  }
  shell("DEL variance fmin.log eigv.rpt") #clean up
 
  #read in simulated datasets, capture and save key pieces
  f<-list.files(path=outdir, pattern="\\.out$")
  
  A <- lapply(f, function(x) reptoRlist(paste(outdir, x, sep='\\'))) #all in one huge list, no loop
  
  names(A)<-paste("seed_", lapply(A, '[[', "rseed"), sep="")  #assign seed to name outputs

  true_bio<-lapply(A, '[[', "avByr")
  sim_bio<-lapply(A, '[[', "est_survey_biomass")

  true_totcatch<-lapply(A, '[[', "catch_biomass")
  sim_totcatch<-lapply(A, '[[', "est_catch_biomass")
  
  saveRDS(true_bio, file.path(outdir, paste("trueB_", er,".rds", sep="")))
  saveRDS(sim_bio, file.path(outdir, paste("simB_", er,".rds", sep="")))
  saveRDS(true_totcatch, file.path(outdir, paste("trueC_", er,".rds", sep="")))
  saveRDS(sim_totcatch, file.path(outdir, paste("simC_", er,".rds", sep="")))
  

  #write files out in SAS input format (columns of catch then B) for each seed
  for (n in names(sim_bio)){
    simoutdat<-cbind(seq(1:dim(sim_bio[[n]])[2]), t(sim_totcatch[[n]]), t(sim_bio[[n]]))
    write.table(simoutdat, file.path(outdir, paste("sim2SAS_", n, ".csv", sep="")), sep=",", col.names=F, row.names=F)  
  }

}#end if windows loop

if(Sys.info()['sysname']=="Linux"){
  setwd("/home/sgaichas/Data/Projects/Size_Struct_Ecomod/MSmod_testing")
  dir.create(outdir<-paste("simbio_", er, sep=""))
  dir.create(file.path(outdir, "crashed"))
  #run hydra_sim n times seeded consecutively: linux
  for(n in 1:nSims){
    #runmod<-paste("../code_aug1/hydra_sim -sim",n)  #writes csv file to current working dir
    runmod<-paste("../code_aug1/hydra_sim -ainp hydra_sim_", er, ".pin -sim ",n, sep="")  #writes csv file to current working dir
    outfile<-paste(outdir,"/simtest_",n,".out", sep="")
    write(system(runmod, intern=T), outfile) #saves full output as well as csv
    krakenin<-paste("simout_",n,".csv", sep="")
    file.rename("simout.csv", file.path(outdir, krakenin))
    nantest<-paste("grep nan", outfile)  #find crashes
    if((system(nantest))==0){
      #copy outfile and krakenin to crashed directory
      file.copy(outfile, file.path(outdir, "crashed"))
      #delete outfile and krakenin from working directory
      file.remove(outfile)
      file.remove(file.path(outdir, krakenin))
    }
  }
  system("rm variance fmin.log eigv.rpt") #clean up
  
  #read in simulated datasets, capture and save key pieces
  f<-list.files(path=outdir, pattern="\\.out$")
  
  A <- lapply(f, function(x) reptoRlist(paste(outdir, x, sep='/'))) #all in one huge list, no loop
  
  names(A)<-paste("seed_", lapply(A, '[[', "rseed"), sep="")  #assign seed to name outputs
  
  true_bio<-lapply(A, '[[', "avByr")
  sim_bio<-lapply(A, '[[', "est_survey_biomass")
  
  true_totcatch<-lapply(A, '[[', "catch_biomass")
  sim_totcatch<-lapply(A, '[[', "est_catch_biomass")
  
  saveRDS(true_bio, file.path(outdir, paste("trueB_", er,".rds", sep="")))
  saveRDS(sim_bio, file.path(outdir, paste("simB_", er,".rds", sep="")))
  saveRDS(true_totcatch, file.path(outdir, paste("trueC_", er,".rds", sep="")))
  saveRDS(sim_totcatch, file.path(outdir, paste("simC_", er,".rds", sep="")))
  
  
  #write files out in SAS input format (columns of catch then B) for each seed
  for (n in names(sim_bio)){
    simoutdat<-cbind(seq(1:dim(sim_bio[[n]])[2]), t(sim_totcatch[[n]]), t(sim_bio[[n]]))
    write.table(simoutdat, file.path(outdir, paste("sim2SAS_", n, ".csv", sep="")), sep=",", col.names=F, row.names=F)  
  }
  
}#end if linux loop

}#end run sim model, save outputs

###############################################

#run kraken, save results back to R

#change error type
for(er in c("noerr", "q", "surverr", "q_surverr")){
  
  wdir<-paste("simbio_", er, sep="")
  sim_bio<-readRDS(file.path(wdir, paste("simB_", er, ".rds", sep="")))

  library(stringr)

  # prepare regular expression
  regexp <- "[[:digit:]]+"

  # process string
  nl<-str_extract(names(sim_bio), regexp)

  for(n in nl){  
    runKraken<-paste("..\\Kraken_Console\\Kraken_Console -p \"..\\Kraken_Console\\param_Test1_SS.csv\" -fit -o \"",wdir,"\\outputtest",n,".txt\" -d \"",wdir,"\\simout_",n,".csv\" -id ",n, sep="")
    shell(runKraken)
  }

  k<-list.files(path=wdir,pattern="\\.txt$")
  K<-lapply(k, function(x) reptoRlist(paste(wdir, x, sep='\\')))

  names(K)<-paste("seed_", lapply(K, '[[', "rSeed"), sep="")  #assign seed to name outputs

  kraken_bio<-lapply(K, '[[', "Kraken_Biomass")
  kraken_catch<-lapply(K, '[[', "Kraken_Catch")

  saveRDS(kraken_bio, file.path(wdir, paste("krakB_", er,".rds", sep="")))
  saveRDS(kraken_catch, file.path(wdir, paste("krakC_", er,".rds", sep="")))

}

###############################################

# Run SAS models, save output to R

#change error type
for(er in c("noerr", "q", "surverr", "q_surverr")){
  
  wdir<-paste("simbio_", er, sep="")
  sim_bio<-readRDS(file.path(wdir, paste("simB_", er, ".rds", sep="")))
  
  #run Production-LV model
  SASoutLV<-list()

  for (n in names(sim_bio)){
    #runSAS<-paste("\"C:\\Program Files\\SASHome\\x86\\SASFoundation\\9.3\\sas.exe\" -sysin ..\\SAS_MSmod\\GB_LV_TypeI_test.sas ICON -NOSPLASH -LOG ..\\SAS_MSmod\\test2.log -PRINT testout3.csv -SYSPARM  hydratestdat_1.csv")
    #runSAS<-paste("\"C:\\Program Files\\SASHome\\x86\\SASFoundation\\9.3\\sas.exe\" -sysin ..\\SAS_MSmod\\GB_LV_TypeI_test.sas ICON -NOSPLASH -LOG ..\\SAS_MSmod\\test2.log -PRINT testout3.csv -SYSPARM  sim2SAS_", n, ".csv", sep="")
    runSAS<-paste("\"C:\\Program Files\\SASHome\\x86\\SASFoundation\\9.3\\sas.exe\" -sysin ..\\SAS_MSmod\\LV_SAS.sas ICON -NOSPLASH -LOG ..\\SAS_MSmod\\LVmod.log -PRINT LVout.csv -SYSPARM  ",wdir,"\\sim2SAS_", n, ".csv", sep="")
    shell(runSAS)
    SASoutLV[[n]]<-read.csv("..\\SAS_MSmod\\LVout.csv")
  }

  #save output
  saveRDS(SASoutLV, file.path(wdir, paste("SASoutLV_", er,".rds", sep="")))
  
  #run Production-LV model with q
  SASadjLVq<-list()
  SASoutLVq<-list()
  SASparsLVq<-list()
  SAS_Badj<-list()
  
  for (n in names(sim_bio)){
    #runSAS<-paste("\"C:\\Program Files\\SASHome\\x86\\SASFoundation\\9.3\\sas.exe\" -sysin ..\\SAS_MSmod\\GB_LV_TypeI_test.sas ICON -NOSPLASH -LOG ..\\SAS_MSmod\\test2.log -PRINT testout3.csv -SYSPARM  hydratestdat_1.csv")
    #runSAS<-paste("\"C:\\Program Files\\SASHome\\x86\\SASFoundation\\9.3\\sas.exe\" -sysin ..\\SAS_MSmod\\GB_LV_TypeI_test.sas ICON -NOSPLASH -LOG ..\\SAS_MSmod\\test2.log -PRINT testout3.csv -SYSPARM  sim2SAS_", n, ".csv", sep="")
    runSAS<-paste("\"C:\\Program Files\\SASHome\\x86\\SASFoundation\\9.3\\sas.exe\" -sysin ..\\SAS_MSmod\\LVq_SAS.sas ICON -NOSPLASH -LOG ..\\SAS_MSmod\\LVqmod.log -PRINT LVqout.csv -SYSPARM  ",wdir,"\\sim2SAS_", n, ".csv", sep="")
    shell(runSAS)
    SASoutLVq[[n]]<-read.csv("..\\SAS_MSmod\\LVqout.csv")
    SASparsLVq[[n]]<-read.csv("..\\SAS_MSmod\\LVqpars.csv")
    SAS_Badj[[n]]<-mapply('/', SASoutLVq[[n]][,4:13],SASparsLVq[[n]][21:30])
    SASadjLVq[[n]]<-cbind(SASoutLVq[[n]][,1:3],SAS_Badj[[n]],SASoutLVq[[n]][,14:33])
  }
  
  #save output
  saveRDS(SASoutLVq, file.path(wdir, paste("SASoutLVq_", er,".rds", sep="")))
  saveRDS(SASadjLVq, file.path(wdir, paste("SASadjLVq_", er,".rds", sep="")))
  
  
  #run Delay Difference model
  SASoutDD<-list()

  for (n in names(sim_bio)){
    #runSAS<-paste("\"C:\\Program Files\\SASHome\\x86\\SASFoundation\\9.3\\sas.exe\" -sysin ..\\SAS_MSmod\\GB_LV_TypeI_test.sas ICON -NOSPLASH -LOG ..\\SAS_MSmod\\test2.log -PRINT testout3.csv -SYSPARM  hydratestdat_1.csv")
    #runSAS<-paste("\"C:\\Program Files\\SASHome\\x86\\SASFoundation\\9.3\\sas.exe\" -sysin ..\\SAS_MSmod\\GB_LV_TypeI_test.sas ICON -NOSPLASH -LOG ..\\SAS_MSmod\\test2.log -PRINT testout3.csv -SYSPARM  sim2SAS_", n, ".csv", sep="")
    runSAS<-paste("\"C:\\Program Files\\SASHome\\x86\\SASFoundation\\9.3\\sas.exe\" -sysin ..\\SAS_MSmod\\DD_SAS.sas ICON -NOSPLASH -LOG ..\\SAS_MSmod\\DDmod.log -PRINT DDout.csv -SYSPARM  ",wdir,"\\sim2SAS_", n, ".csv", sep="")
    shell(runSAS)
    SASoutDD[[n]]<-read.csv("..\\SAS_MSmod\\DDout.csv")
  }
  
  #save output
  saveRDS(SASoutDD, file.path(wdir, paste("SASoutDD_", er,".rds", sep="")))
  
}  

###############################################

#compare kraken and SAS results with true hydra results

#change error type
for(er in c("noerr", "q", "surverr", "q_surverr")){
  
  wdir<-paste("simbio_", er, sep="")
  
  sim_bio<-readRDS(file.path(wdir, paste("simB_", er, ".rds", sep="")))
  true_bio<-readRDS(file.path(wdir, paste("trueB_", er, ".rds", sep="")))
  sim_totcatch<-readRDS(file.path(wdir, paste("simC_", er, ".rds", sep="")))
  true_totcatch<-readRDS(file.path(wdir, paste("trueC_", er, ".rds", sep="")))

  kraken_bio<-readRDS(file.path(wdir, paste("krakB_", er, ".rds", sep="")))
  kraken_catch<-readRDS(file.path(wdir, paste("krakC_", er, ".rds", sep="")))
  
  SASoutLV<-readRDS(file.path(wdir, paste("SASoutLV_", er, ".rds", sep="")))
  SASoutLVq<-readRDS(file.path(wdir, paste("SASoutLVq_", er, ".rds", sep="")))
  SASadjLVq<-readRDS(file.path(wdir, paste("SASadjLVq_", er, ".rds", sep="")))
  SASoutDD<-readRDS(file.path(wdir, paste("SASoutDD_", er, ".rds", sep="")))

  # Plot true and sim bio each species across all seeds
  png(paste("True_B_sim_", er,".png", sep=""), width=1200, height=1200, units="px", pointsize=20)
  par(mfrow=c(4,3))
  par(mar=c(2,2,3,2)+0.1)
  par(oma=c(2,2,2,0))
  true_bio<-lapply(true_bio, function(x) replace(x, is.infinite(x),NA))
  sim_bio<-lapply(sim_bio, function(x) replace(x, is.infinite(x),NA))
  maxplott<-as.data.frame(lapply(true_bio,FUN=function(x)apply(x,MARGIN=1,FUN=max, na.rm=T)))
  maxplots<-as.data.frame(lapply(sim_bio,FUN=function(x)apply(x,MARGIN=1,FUN=max, na.rm=T)))
  #minplot<-as.data.frame(lapply(krak_relerr,FUN=function(x)apply(x,MARGIN=1,FUN=min, na.rm=T)))
  for(i in 1:dim(true_bio[[1]])[1]){
    plot(true_bio[[1]][i,], ylim=c(0,max(max(maxplott[i,]), max(maxplots[i,]))), main=SpNames[i], type="l", col = rgb(0, 0, 1, 0.3))
    for (j in 1:length(true_bio)){
      lines(true_bio[[j]][i,], col = rgb(0, 0, 1, 0.3))
      lines(sim_bio[[j]][i,], col = rgb(1, 0, 0, 0.3))
    }
    #abline(h=0.0, col="blue")
  }
  mtext(paste("True and Simulated Biomass ", length(true_bio), "runs", er), outer=T, side=3)
  mtext("Biomass", outer=T, side=2)
  mtext("year", outer=T, side=1)
  dev.off()
  
  singleseed<-0
  
if(singleseed==1){
  #each species time series for a single seed: Kraken
  #biomass
  for (j in 1:length(true_bio)){
    par(mfrow=c(4,3))
    par(mar=c(2,2,3,2)+0.1)
    par(oma=c(2,2,2,0))
    for(i in 1:dim(true_bio[[1]])[1]){
      maxplot<-max(true_bio[[j]][i,], sim_bio[[j]][i,], kraken_bio[[j]][i,])
      plot(seq(1:53), true_bio[[j]][i,], ylim=c(0,maxplot), main=SpNames[i], col="blue")
      points(seq(1:53), sim_bio[[j]][i,], col="red")
      lines(seq(1:53), kraken_bio[[j]][i,])
    }
    mtext(paste("Kraken Biomass ", names(true_bio[j]), er), outer=T, side=3)
  }

  #each species time series for a single seed
  #catch
  #for (j in 1:length(true_totcatch)){
  #  par(mfrow=c(4,3))
  #  par(mar=c(2,2,3,2)+0.1)
  #  par(oma=c(2,2,2,0))
  #  for(i in 1:dim(true_totcatch[[1]])[1]){
  #    maxplot<-max(true_totcatch[[j]][i,], sim_totcatch[[j]][i,], kraken_catch[[j]][i,])
  #    plot(seq(1:53), true_totcatch[[j]][i,], ylim=c(0,maxplot), main=SpNames[i])
  #    lines(seq(1:53), kraken_catch[[j]][i,])
  #    lines(seq(1:53), sim_totcatch[[j]][i,])
  #  }
  #  mtext(paste("Total catch ", names(true_totcatch[j]), er), outer=T, side=3)
  #}

  #each species time series for a single seed: SAS LV model
  #biomass
  for (j in 1:length(true_bio)){
    par(mfrow=c(4,3))
    par(mar=c(2,2,3,2)+0.1)
    par(oma=c(2,2,2,0))
    for(i in 1:dim(true_bio[[1]])[1]){
      maxplot<-max(true_bio[[j]][i,], sim_bio[[j]][i,], SASoutLV[[j]][,i+3], na.rm=T)
      plot(true_bio[[j]][i,-(1:4)], ylim=c(0,maxplot), main=SpNames[i], col="blue")
      points(sim_bio[[j]][i,-(1:4)], col="red")      
      lines(SASoutLV[[j]][,i+3])
    }
    mtext(paste("SAS LV Biomass ", names(true_bio[j]), er), outer=T, side=3)
  }

  #each species time series for a single seed: SAS LV q model
  #biomass
  for (j in 1:length(true_bio)){
    par(mfrow=c(4,3))
    par(mar=c(2,2,3,2)+0.1)
    par(oma=c(2,2,2,0))
    for(i in 1:dim(true_bio[[1]])[1]){
      maxplot<-max(true_bio[[j]][i,], sim_bio[[j]][i,], SASoutLVq[[j]][,i+3], na.rm=T)
      #maxplot<-max(true_bio[[j]][i,], sim_bio[[j]][i,], SASadjLVq[[j]][,i+3], na.rm=T)
      plot(true_bio[[j]][i,-(1:4)], ylim=c(0,maxplot), main=SpNames[i], col="blue")
      points(sim_bio[[j]][i,-(1:4)], col="red")      
      lines(SASoutLVq[[j]][,i+3])
      #lines(SASadjLVq[[j]][,i+3])
    }
    mtext(paste("SAS LV q Biomass raw", names(true_bio[j]), er), outer=T, side=3)
  }

  #each species time series for a single seed: SAS LV q model adj
  #biomass
  for (j in 1:length(true_bio)){
    par(mfrow=c(4,3))
    par(mar=c(2,2,3,2)+0.1)
    par(oma=c(2,2,2,0))
    for(i in 1:dim(true_bio[[1]])[1]){
      #maxplot<-max(true_bio[[j]][i,], sim_bio[[j]][i,], SASoutLVq[[j]][,i+3], na.rm=T)
      maxplot<-max(true_bio[[j]][i,], sim_bio[[j]][i,], SASadjLVq[[j]][,i+3], na.rm=T)
      plot(true_bio[[j]][i,-(1:4)], ylim=c(0,maxplot), main=SpNames[i], col="blue")
      points(sim_bio[[j]][i,-(1:4)], col="red")      
      #lines(SASoutLVq[[j]][,i+3])
      lines(SASadjLVq[[j]][,i+3])
    }
    mtext(paste("SAS LV q Biomass adj", names(true_bio[j]), er), outer=T, side=3)
  }
  
  #each species time series for a single seed: SAS DD model
  #biomass
  for (j in 1:length(true_bio)){
    par(mfrow=c(4,3))
    par(mar=c(2,2,3,2)+0.1)
    par(oma=c(2,2,2,0))
    for(i in 1:dim(true_bio[[1]])[1]){
      maxplot<-max(true_bio[[j]][i,], SASoutDD[[j]][,i+3], na.rm=T)
      plot(true_bio[[j]][i,-(1:4)], ylim=c(0,maxplot), main=SpNames[i], col="blue")
      points(sim_bio[[j]][i,-(1:4)], col="red")      
      lines(SASoutDD[[j]][,i+3])
    }
    mtext(paste("SAS DD Biomass ", names(true_bio[j]), er), outer=T, side=3)
  }
  
  
}#end singleseed 

  #calculate error for plotting and skill metrics

  krak_err<-mapply('-',kraken_bio,true_bio,SIMPLIFY=FALSE)
  krak_relerr<-mapply('/', krak_err, true_bio, SIMPLIFY=FALSE)
  
  # reshape SAS output bio to have same dimensions as true

  SASLV_bio<-lapply(SASoutLV, "[", c(4:13))
  SASLV_bio<-lapply(SASLV_bio, t)
  true_bioSAS<-lapply(true_bio, function(x)x[,-c(1:4)])
  SASLV_err<-mapply('-',SASLV_bio,true_bioSAS,SIMPLIFY=FALSE)
  SASLV_relerr<-mapply('/', SASLV_err,true_bioSAS, SIMPLIFY=FALSE)

  #SASLVq_bio<-lapply(SASoutLVq, "[", c(4:13))
  SASLVq_bio<-lapply(SASadjLVq, "[", c(4:13))
  SASLVq_bio<-lapply(SASLVq_bio, t)
  SASLVq_err<-mapply('-',SASLVq_bio,true_bioSAS,SIMPLIFY=FALSE)
  SASLVq_relerr<-mapply('/', SASLVq_err,true_bioSAS, SIMPLIFY=FALSE)
  
  SASDD_bio<-lapply(SASoutDD, "[", c(4:13))
  SASDD_bio<-lapply(SASDD_bio, t)
  SASDD_err<-mapply('-',SASDD_bio,true_bioSAS,SIMPLIFY=FALSE)
  SASDD_relerr<-mapply('/', SASDD_err,true_bioSAS, SIMPLIFY=FALSE)
  
  #ensemble models--equal weight LEAVE OUT SAS LV q
  ensSAS_bio<-mapply('+',SASLV_bio, SASDD_bio, SIMPLIFY=FALSE)
  ensAll_bio<-mapply('+',ensSAS_bio, lapply(kraken_bio, function(x)x[,-c(1:4)]), SIMPLIFY=FALSE)
  ensAll_bio<-lapply(ensAll_bio, function(x) x/3)
  ensAll_err<-mapply('-',ensAll_bio,true_bioSAS,SIMPLIFY=FALSE)
  ensAll_relerr<-mapply('/', ensAll_err,true_bioSAS, SIMPLIFY=FALSE)
  
if(singleseed==1){
    
  #each species time series for a single seed: all models
  for (j in 1:length(true_bio)){
    if (j %in% c(13, 43, 96)) png(paste("Allmod_singlefit_sim",j, "_", er, ".png", sep=""), width=1200, height=1200, units="px", pointsize=20)
    par(mfrow=c(4,3))
    par(mar=c(2,2,3,2)+0.1)
    par(oma=c(2,2,2,0))
    for(i in 1:dim(true_bio[[1]])[1]){
      maxplot<-max(true_bio[[j]][i,], sim_bio[[j]][i,], kraken_bio[[j]][i,], SASoutDD[[j]][,i+3], SASoutLV[[j]][,i+3], na.rm=T)
      plot(true_bio[[j]][i,-(1:4)], ylim=c(0,maxplot), main=SpNames[i], col="blue")
      if(er %in% c("q", "surverr", "q_surverr")) points(sim_bio[[j]][i,-(1:4)], col="red")    
      lines(kraken_bio[[j]][i,-(1:4)], col="darkgreen", lwd=2)
      lines(SASoutLV[[j]][,i+3], col="mediumvioletred", lwd=2)
      lines(SASoutDD[[j]][,i+3], col="mediumslateblue", lwd=2)
      lines(ensAll_bio[[j]][i,], col="black", lwd=2)
    }
    plot(1, type = "n", axes=FALSE, xlab="", ylab="")
    plot_colors <- c("blue","red", "darkgreen", "mediumvioletred", "mediumslateblue", "black")
    legend(x = "top",inset = 0,
           legend = c("True biomass", "Observed biomass", "Kraken estimate","SAS LV estimate", "SAS DD estimate", "Ensemble"), 
           col=plot_colors, pch=c(1,1,NA, NA, NA, NA), lwd=c(NA, NA, 3, 3, 3, 3), bty="n", cex=1.7)
    mtext(paste("All model biomass time series fits", names(true_bio[j]), er), outer=T, side=3)
    dev.off()
  }
  
}#end singleseed 

  # Kraken error for each species across all seeds
  png(paste("Kraken_Berr_", er,".png", sep=""), width=1200, height=1200, units="px", pointsize=20)
  par(mfrow=c(4,3))
  par(mar=c(2,2,3,2)+0.1)
  par(oma=c(2,2,2,0))
  krak_relerr<-lapply(krak_relerr, function(x) replace(x, is.infinite(x),NA))
  maxplot<-as.data.frame(lapply(krak_relerr,FUN=function(x)apply(x,MARGIN=1,FUN=max, na.rm=T)))
  minplot<-as.data.frame(lapply(krak_relerr,FUN=function(x)apply(x,MARGIN=1,FUN=min, na.rm=T)))
  for(i in 1:dim(krak_relerr[[1]])[1]){
#    plot(krak_relerr[[1]][i,], ylim=c(min(minplot[i,],-1.0),max(maxplot[i,],1.0)), main=SpNames[i], type="l")
    plot(krak_relerr[[1]][i,], ylim=c(-1,10), main=SpNames[i], type="l", col = rgb(0, 0, 0, 0.3))
    for (j in 1:length(krak_relerr)){
      lines(krak_relerr[[j]][i,], col = rgb(0, 0, 0, 0.3))
    }
    abline(h=0.0, col="blue", lwd=3)
  }
  mtext(paste("Kraken Biomass Error ", length(true_bio), "runs", er), outer=T, side=3)
  mtext("error", outer=T, side=2)
  mtext("year", outer=T, side=1)
  dev.off()

#  png(paste("Kraken_Berr_", er,".png", sep=""), width=1200, height=1200, units="px", pointsize=20)
#  par(mfrow=c(4,3))
#  par(mar=c(2,2,3,2)+0.1)
#  par(oma=c(2,2,2,0))
#  krak_relerr<-lapply(krak_relerr, function(x) replace(x, is.infinite(x),NA))
#  maxplot<-as.data.frame(lapply(krak_relerr,FUN=function(x)apply(x,MARGIN=1,FUN=quantile, probs=0.95, na.rm=T)))
#  minplot<-as.data.frame(lapply(krak_relerr,FUN=function(x)apply(x,MARGIN=1,FUN=quantile, probs=0.05,na.rm=T)))
#  for(i in 1:dim(krak_relerr[[1]])[1]){
#    #plot(krak_relerr[[1]][i,], ylim=c(min(minplot[i,],-1.0),max(maxplot[i,],1.0)), main=SpNames[i], type="l")
#    plot(krak_relerr[[1]][i,], ylim=c(-1,10), main=SpNames[i], type="l")
#    for (j in 1:length(krak_relerr)){
#    lines(krak_relerr[[j]][i,])
#    }
#      abline(h=0.0, col="blue")
#  }
#  mtext(paste("Kraken Biomass Error ", length(true_bio), "runs", er), outer=T, side=3)
#  mtext("error", outer=T, side=2)
#  mtext("year", outer=T, side=1)
#  dev.off()


  # SAS LV error for each species across all seeds
  png(paste("SASLV_Berr_", er,".png", sep=""), width=1200, height=1200, units="px", pointsize=20)
  par(mfrow=c(4,3))
  par(mar=c(2,2,3,2)+0.1)
  par(oma=c(2,2,2,0))
  SASLV_relerr<-lapply(SASLV_relerr, function(x) replace(x, is.infinite(x),NA))
  maxplot<-as.data.frame(lapply(SASLV_relerr,FUN=function(x)apply(x,MARGIN=1,FUN=max, na.rm=T)))
  minplot<-as.data.frame(lapply(SASLV_relerr,FUN=function(x)apply(x,MARGIN=1,FUN=min, na.rm=T)))
  for(i in 1:dim(SASLV_relerr[[1]])[1]){
#    plot(SASLV_relerr[[1]][i,], ylim=c(min(minplot[i,],-1.0),max(maxplot[i,],1.0)), main=SpNames[i], type="l")
    plot(SASLV_relerr[[1]][i,], ylim=c(-1,10), main=SpNames[i], type="l", col = rgb(0, 0, 0, 0.3))
    for (j in 1:length(SASLV_relerr)){
      lines(SASLV_relerr[[j]][i,], col = rgb(0, 0, 0, 0.3))
    }
    abline(h=0.0, col="blue", lwd=3)
  }
  mtext(paste("SAS LV Biomass Error ", length(true_bio), "runs", er), outer=T, side=3)
  mtext("error", outer=T, side=2)
  mtext("year", outer=T, side=1)
  dev.off()

  # SAS LV q error for each species across all seeds
#  png(paste("SASLVq_Berr_", er,".png", sep=""), width=1200, height=1200, units="px", pointsize=20)
#  par(mfrow=c(4,3))
#  par(mar=c(2,2,3,2)+0.1)
#  par(oma=c(2,2,2,0))
#  SASLVq_relerr<-lapply(SASLVq_relerr, function(x) replace(x, is.infinite(x),NA))
#  maxplot<-as.data.frame(lapply(SASLVq_relerr,FUN=function(x)apply(x,MARGIN=1,FUN=max, na.rm=T)))
#  minplot<-as.data.frame(lapply(SASLVq_relerr,FUN=function(x)apply(x,MARGIN=1,FUN=min, na.rm=T)))
#  for(i in 1:dim(SASLVq_relerr[[1]])[1]){
#    plot(SASLVq_relerr[[1]][i,], ylim=c(min(minplot[i,],-1.0),max(maxplot[i,],1.0)), main=SpNames[i], type="l")
#    for (j in 1:length(SASLVq_relerr)){
#      lines(SASLVq_relerr[[j]][i,])
#    }
#    abline(h=0.0, col="blue")
#  }
#  mtext(paste("SAS LV q Biomass Error ", length(true_bio), "runs", er), outer=T, side=3)
#  mtext("error", outer=T, side=2)
#  mtext("year", outer=T, side=1)
#  dev.off()

# SAS DD error for each species across all seeds
  png(paste("SASDD_Berr_", er,".png", sep=""), width=1200, height=1200, units="px", pointsize=20)
  par(mfrow=c(4,3))
  par(mar=c(2,2,3,2)+0.1)
  par(oma=c(2,2,2,0))
  SASDD_relerr<-lapply(SASDD_relerr, function(x) replace(x, is.infinite(x),NA))
  maxplot<-as.data.frame(lapply(SASDD_relerr,FUN=function(x)apply(x,MARGIN=1,FUN=max, na.rm=T)))
  minplot<-as.data.frame(lapply(SASDD_relerr,FUN=function(x)apply(x,MARGIN=1,FUN=min, na.rm=T)))
  for(i in 1:dim(SASDD_relerr[[1]])[1]){
#    plot(SASDD_relerr[[1]][i,], ylim=c(min(minplot[i,],-1.0),max(maxplot[i,],1.0)), main=SpNames[i], type="l")
    plot(SASDD_relerr[[1]][i,], ylim=c(-1,10), main=SpNames[i], type="l", col = rgb(0, 0, 0, 0.3))
    for (j in 1:length(SASDD_relerr)){
      lines(SASDD_relerr[[j]][i,], col = rgb(0, 0, 0, 0.3))
    }
    abline(h=0.0, col="blue", lwd=3)
  }
  mtext(paste("SAS DD Biomass Error ", length(true_bio), "runs", er), outer=T, side=3)
  mtext("error", outer=T, side=2)
  mtext("year", outer=T, side=1)
  dev.off()

  # simple ensemble error for each species across all seeds
  png(paste("ensAll_Berr_", er,".png", sep=""), width=1200, height=1200, units="px", pointsize=20)
  par(mfrow=c(4,3))
  par(mar=c(2,2,3,2)+0.1)
  par(oma=c(2,2,2,0))
  ensAll_relerr<-lapply(ensAll_relerr, function(x) replace(x, is.infinite(x),NA))
  maxplot<-as.data.frame(lapply(ensAll_relerr,FUN=function(x)apply(x,MARGIN=1,FUN=max, na.rm=T)))
  minplot<-as.data.frame(lapply(ensAll_relerr,FUN=function(x)apply(x,MARGIN=1,FUN=min, na.rm=T)))
  for(i in 1:dim(ensAll_relerr[[1]])[1]){
#    plot(ensAll_relerr[[1]][i,], ylim=c(min(minplot[i,],-1.0),max(maxplot[i,],1.0)), main=SpNames[i], type="l")
    plot(ensAll_relerr[[1]][i,], ylim=c(-1,10), main=SpNames[i], type="l", col = rgb(0, 0, 0, 0.3))
    for (j in 1:length(ensAll_relerr)){
      lines(ensAll_relerr[[j]][i,], col = rgb(0, 0, 0, 0.3))
    }
    abline(h=0.0, col="blue", lwd=3)
  }
  mtext(paste("Model Ensemble Biomass Error ", length(true_bio), "runs", er), outer=T, side=3)
  mtext("error", outer=T, side=2)
  mtext("year", outer=T, side=1)
  dev.off()
  
  
  #one page, four MEF boxplots
  png(paste(er,"_MEF.png", sep=""), width=1200, height=1200, units="px", pointsize=20)
  par(mfrow=c(2,2))
  par(mar=c(4,2,2,1)+0.1)
  par(oma=c(0,0,0,0))
  
  #loop over model errors for each skill metric
  #moderr<-list(krak_err, SASLV_err, SASLVq_err, SASDD_err)
  #names(moderr)<-c("krak", "SASLV", "SASLVq", "SASDD")
  moderr<-list(krak_err, SASLV_err, SASDD_err, ensAll_err)
  names(moderr)<-c("krak", "SASLV", "SASDD", "ensAll")
  for(e in 1:length(moderr)){
    
    skill<-list()

    #rmse by species for each seed
    sp.rmse<-t(as.data.frame(rapply(moderr[[e]], function(i) apply(i, 1, rmse), how="replace")))
    colnames(sp.rmse)<-SpNames
    skill$sp.rmse<-sp.rmse
    
    #rmse for all species combined by seed
    skill$tot.rmse<-unlist(lapply(moderr[[e]], rmse))
    #mapply(rmse, krak_err) #does the same thing

    #aae by species for each seed
    sp.aae<-t(as.data.frame(rapply(moderr[[e]], function(i) apply(i, 1, aae), how="replace")))
    colnames(sp.aae)<-SpNames
    skill$sp.aae<-sp.aae
    
    #aae for all species combined by seed
    skill$tot.aae<-unlist(lapply(moderr[[e]], aae))

    #ae by species for each seed
    sp.ae<-t(as.data.frame(rapply(moderr[[e]], function(i) apply(i, 1, ae), how="replace")))
    colnames(sp.ae)<-SpNames
    skill$sp.ae<-sp.ae
    
    #ae for all species combined by seed
    skill$tot.ae<-unlist(lapply(moderr[[e]], ae))

    #mef by species for each seed BUT rapply requires single argument function
    mef_out<-matrix(rep(0),length(true_bio),dim(true_bio[[1]])[1])
    for(i in 1:length(true_bio)){
      for(j in 1:dim(true_bio[[1]])[1]){
       mef_out[i,j]<-mef(true_bio[[i]][j,], moderr[[e]][[i]][j,])
      }
    }
    rownames(mef_out)<-names(true_bio)
    colnames(mef_out)<-SpNames
    skill$sp.mef<-mef_out

    #mef for all species combined by seed
    skill$tot.mef<-unlist(mapply(mef, true_bio, moderr[[e]]))
    
    saveRDS(skill, file.path(wdir, paste(names(moderr)[e],"_", er,".rds", sep="")))

    #boxplot(skill$sp.mef, las=3, main=paste(names(moderr)[e],"", er," MEF"), ylim=c(min(-1, min(skill$sp.mef, na.rm=T)), 1))
    #boxplot(cbind(skill$sp.mef,"all"=skill$tot.mef), col=c(spcol, "white"), notch=T, las=3, main=paste(names(moderr)[e],"", er," MEF"), ylim=c(-15, 1))
    boxplot(cbind(skill$sp.mef,"all"=skill$tot.mef), col=c(spcol, "white"), notch=T, las=3, ylim=c(-15, 1))
    abline(h=0, lty=3)
    abline(h=1, lty=3)
    #abline(v=10.5, lty=3)
    
  }
  dev.off()

}# end compare loop



#correlation for each species by seed NOT DONE FOR SAS MODS
cor_out<-matrix(rep(0),length(true_bio),dim(true_bio[[1]])[1])

for(i in 1:length(true_bio)){
  for(j in 1:dim(true_bio[[1]])[1]){
    cor_out[i,j]<-cor(true_bio[[i]][j,], kraken_bio[[i]][j,], method="spearman")
  }
}
rownames(cor_out)<-names(true_bio)
colnames(cor_out)<-SpNames

#correlation for all time series together by seed
mapply(cor.test,kraken_bio,true_bio,method="spearman",SIMPLIFY=FALSE)

}


###############################################


