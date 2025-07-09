
##### ---------------------------------------------------------------
##### Fit population-level model for `Known-aged Breeder` using simulated data
##### Author: C. George Glen
##### Date: 01/17/2023 
##### ---------------------------------------------------------------
rm(list=ls()) # remove everything currently held in the R memory
graphics.off() # close all open graphics window
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


wd.path.name  <- "/Users/constantingeorgeglen/Documents/University/UFL/RESEARCH/Projects_Writing/PhDChapters/ch2.senescenceCTC/submission/code-upload/reproductiveMod"
setwd(wd.path.name)


##### Load the data
source(paste0(wd.path.name,"/0.load.libs.R"))
source(paste0(wd.path.name,"/0.POPmodFun.v6.R"))
load(paste0(wd.path.name,"/test.pop.data.RData"))
#load(paste0(wd.path.name,"/test.cohort-unstruct.data.RData"))


##### Load all the model texts
model.path = paste0(wd.path.name)
modellist  = list.files(path = model.path, pattern = ".*.txt")
modeltexts = lapply(modellist, function(x) source(paste0(model.path,"/",x)))


##### Some preamble data prep
keep    <- c("GSSmle","GSSci")
models  <- c(fit=fGSSPois, pred=fGSSPoishofXgY)
data    <- list( Ymat=Ymat, agematrix=agematrix, 
                 turtlematrix=turtlematrix, trueagematrix=trueagematrix,
                 Xmat=Xmat, lens=lens, nind=nind)


##### Data cloning settings
K     <- 1



##### MCMC settings
MCMClist   <- list('n.iter'=10000, 'n.burn'=5000)
#MCMClist   <- lapply(MCMClist, function(v) v / 100)
MCMClist$n.thin         <- 10
MCMClist$n.adapt        <- 1000
MCMClist$n.chains.bayes <- 3
MCMClist$n.chains.dc    <- 3


MCMC.initials <- function(){
  list('Kparam'=9900,'theta'=0.09,'sigma'=1.1) ### sigma here is on the log scale
}

##### Run the models
GSSpopulationmod   <- GSSPopmodfun(MCMClist=MCMClist, 
                                   data=data, K=K, model=fGSSPois, 
                                   parallel = F, inits = MCMC.initials,
                                   priors = list(K.shape=110, K.rate=0.01, 
                                                 theta.a=1, theta.b=20))


GSSpopulationmod$GSSmle ## ML estimates
GSSpopulationmod$GSSsd  ## SE in ML estimates
GSSpopulationmod$GSSci  ## CI for ML estimates 


#### Plot the model
dat <- knownagetest
dat <- dat %>% rename( cumeggs=cumout ) %>%
  group_by(TAGID) %>% 
  mutate(ISRP=AGE-lag(AGE)) ## create a variable for the waiting time between breeding seasons
# plot(table(dat$ISRP))

# tau <- 1/median(dat$ISRP, na.rm=TRUE)
tau <- 0.9
data.list <- list(models       = models,
                  MCMClist     = MCMClist,
                  rawdata      = dat,
                  data         = data,
                  nsims        = 20000, 
                  add.ind.pred = F, add.traj = F,
                  timemeasure  = "Age", K=GSSpopulationmod$K,
                  fitted.model = GSSpopulationmod,
                  pred.future  = F, pNfuture=20,
                  tau          = tau, my.cols="grey40")

pop.predict <- do.call(pop.predict.trajectory, data.list)
abline( reg = lm(cumeggs ~ AGE, data=dat ), lty="dashed", lwd=4, col="red" )





