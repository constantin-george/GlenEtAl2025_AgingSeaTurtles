
##### ---------------------------------------------------------------
##### Fit population-level model for `Known-aged Breeder`  using simulated data
##### Author: C. George Glen
##### Date: 01/17/2023 
##### ---------------------------------------------------------------
wd.path.name  <- "/Users/constantingeorgeglen/Documents/University/UFL/RESEARCH/Projects_Writing/PhDChapters/ch2.senescenceCTC/submission/code-upload/reproductiveMod"
setwd(wd.path.name)


##### Load the data
source(paste0(wd.path.name,"/0.load.libs.R"))
source(paste0(wd.path.name,"/0.POPmodFun.v6.R"))
load(paste0(wd.path.name,"/test.pop.data.RData"))


##### Load all the model texts
model.path = paste0(wd.path.name)
modellist  = list.files(path = model.path, pattern = ".*.txt")
modeltexts = lapply(modellist, function(x) source(paste0(model.path,"/",x)))



##### Some preamble data prep
keep    <- c("GSSmle","GSSci")
models  <- c(fit=fGSSPois, pred=fGSSPoishofXgY)
data    <- list( Ymat=Ymat, agematrix=agematrix, 
                 turtlematrix=turtlematrix, 
                 Xmat=Xmat, lens=lens, nind=nind)



##### Data cloning settings
Kvec  <- 2^(seq(1,3,1))
### [1]    2    4    8   


##### MCMC settings
MCMClist   <- list('n.iter'=50000, 'n.burn'=500,
                   'n.thin'=5,'n.adapt'=1000,
                   'n.chains.bayes'=5,'n.chains.dc'=1)



MCMC.initials <- function(){
  list('Kparam'=14810,'theta'=0.06759,'sigma'=1.115) ### sigma here is on the log scale
}

GssPriors <- list(theta.a  = 3, theta.b  = 3, K.shape  = 115, K.rate   = 0.01)

##### Run the models
GSSpopulationident <- IdentifiableGSSPopmodfun(MCMClist=MCMClist, inits=MCMC.initials, data=data, 
                                               Kvec= Kvec, model=fGSSPois, 
                                               parallel.bayes  = T, parallel.dc = F,
                                               priors=GssPriors, n.cores  = 5)

##### plot for identifiability
plot.identifiable.pop(Kvec=Kvec, data = GSSpopulationident )

