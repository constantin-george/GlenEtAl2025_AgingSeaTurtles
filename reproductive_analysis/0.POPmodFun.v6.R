##### ---------------------------------------------------------------
##### Model function
##### Author: George Glen
##### Date: 01/17/2023 
##### ---------------------------------------------------------------

##### Currently the following function are written and work
##### 1: Run a DC model by 
#####   - Population
##### 2: Assess parameter identifiability by 
#####   - Population
##### 3: Fitted model/ future predictions by 
#####   - Population

##### ---------------------------------------------------------------
##### Identifiability function
##### ---------------------------------------------------------------

##### This function is designed to perform data cloning by population where you can set distribution
##### Current arguments are
##### 1: MCMC settings - niter, n.chains, n.adapt, n.update, thin
##### 2: param         - parameters to keep track of
##### 3: priors        - for k, sigma, and theta
##### 4: data 
#####     - Ymat is a matrix of cumulative reproductive output 
#####     - agematrix is a matrix of the differences in ages between breeding season i and i+1
#####     - turtlematrix is a matrix of the differences in cumout between breeding season i and i+1
##### rows are females and columns are ages


##### This function is designed to preform data cloning by population
##### Same arguments as above except femaleID is removed

#' Arguments
#' @param MCMClist a list contain the MCMC setting for JAGS
#' @param data list containing matrices of cumulative reproductive output per female (Ymat),
#'             annual reproductive output (turtlematrix), and 
#'             waiting times between reproductive seasons (agematrix)
#' @param K specifying the number of clones for data cloning
#' @param model JAGS model 
#' @param plot if TRUE returns MCMC diagnostic plots for the parameters (params) of interest
#' @param parallel if TRUE runs the MCMC chains in parallel
#' @param inits list of initial values for the parameters of interest
#' @param priors list of priors for theta and kappa (Kparm):
#'               theta is assumed to have a beta prior with shape parameters theta.a and theta.b
#'               kappa is assumed to have a gamma prior with shape parameters K.shape and K.rate
#' @param params parameters to monitor
#' @return list containing K clones, 
#'                         the MLE, sd for the MLE, and wald type-ci
#'                         the fitted DC object
#'                         the MCMC samples (all.nodes)

GSSPopmodfun <- function(MCMClist = MCMClist, 
                         data     = data, 
                         K        = K,
                         model    = model,
                         plot     = F,
                         parallel = F,
                         inits    = inits,
                         priors   = priors, n.cores = 10,
                         params   = c("Kparam","theta","sigma") ){
  
  #### Load in the data
  Ymat         <- data$Ymat
  agematrix    <- data$agematrix
  turtlematrix <- data$turtlematrix
  lens         <- if(!is.numeric(data$lens)){ as.numeric(data$lens) }else{data$lens}
  nind         <- data$nind
  
  #### array[,,,] = row, column, array number
  m              <- dim(turtlematrix)[1]
  n              <- dim(turtlematrix)[2]
  dcYmat         <- array(Ymat, c(m,n,K))
  dcturtlematrix <- array(turtlematrix, c(m,n,K))
  dcagematrix    <- array(agematrix, c(m,n,K))
  dclens         <- lens
  
  #### Set up the model
  ssm.list <- append(list(Y       = dcYmat,
                          lens    = lens, 
                          ind     = nind, 
                          Xtrans  = dcagematrix,
                          Ytrans  = dcturtlematrix,
                          K       = K), priors)
  
  #### Run the model
  #### Make sure that even if parallel=T but nchain=1, parallel=F
  if(MCMClist$n.chains.dc==1) parallel <- F
  
  if(parallel == F){
    GSSfit <- dclone::jags.fit(data     = ssm.list, 
                               params   = params, 
                               model    = model, 
                               inits    = inits,
                               n.chains = MCMClist$n.chains.dc,
                               n.adapt  = MCMClist$n.adapt,
                               n.update = MCMClist$n.burn,
                               n.iter   = MCMClist$n.iter,
                               thin     = MCMClist$n.thin)
  }else{
    
    ### Load in the necessary packages
    require(foreach)
    require(parallel)
    require(doParallel)
    require(snow)
    require(snowfall)
    
    #### Set up the cluster
    # n.cores <- max(1L, detectCores(), na.rm = TRUE)
    # n.cores <- min(parallelly::freeConnections(), n.cores)
    cl      <- parallel::makePSOCKcluster(n.cores)
    parallel::clusterEvalQ(cl, library(dclone))
    dclone::parLoadModule(cl,"lecuyer")
    
    #### Run the model
    GSSfit <- dclone::jags.parfit(cl       = cl,
                                  data     = ssm.list,
                                  params   = params,
                                  model    = model,
                                  inits    = inits,
                                  n.chains = MCMClist$n.chains.dc,
                                  n.adapt  = MCMClist$n.adapt,
                                  n.update = MCMClist$n.burn,
                                  n.iter   = MCMClist$n.iter,
                                  thin     = MCMClist$n.thin)
    parallel::stopCluster(cl)
  }
  
  ### Estimates
  GSSmle    <- summary(GSSfit)[[1]][,1]
  GSSsdmle  <- sqrt(K)*(summary(GSSfit)[[1]][,2])
  all.nodes <- do.call(rbind, GSSfit)
  
  #### Wald-type CI
  fish.inv         <- K*var(all.nodes)
  alpha            <- 0.05
  z.alpha.half     <- qnorm(p=(1-alpha/2))
  K.std.error      <- z.alpha.half*sqrt(fish.inv[1,1])
  sigma.std.error  <- z.alpha.half*sqrt(fish.inv[2,2])
  theta.std.error  <- z.alpha.half*sqrt(fish.inv[3,3])
  K.cis            <- c(GSSmle[1]-K.std.error,GSSmle[1]+K.std.error)
  sigma.cis        <- c(GSSmle[2]-sigma.std.error,GSSmle[2]+sigma.std.error)
  theta.cis        <- c(GSSmle[3]-theta.std.error,GSSmle[3]+theta.std.error)
  GSSci            <- cbind(K.cis,sigma.cis,theta.cis)
  
  #### Show a plot
  if(plot==T) pairs(GSSfit)
  
  #### Print out results of all relevant calcualtion
  print(GSSmle); print(GSSci); 
  
  
  invisible(list(K           = K,
                 GSSmle      = GSSmle, 
                 GSSsdmle    = GSSsdmle, 
                 GSSci       = GSSci, 
                 model.file  = GSSfit,
                 MCMCsamples = all.nodes))
}


##### ---------------------------------------------------------------
##### Identifiability function
##### ---------------------------------------------------------------

##### This function is designed to preform data cloning to assess for parameter identifiability for the population
#' Arguments
#' @param MCMClist a list contain the MCMC setting for JAGS
#' @param data list containing matrices of cumulative reproductive output per female (Ymat),
#'             annual reproductive output (turtlematrix), and 
#'             waiting times between reproductive seasons (agematrix)
#' @param Kvec specifying a vecctor of clones for data cloning to iterate over
#' @param model JAGS model 
#' @param parallel.bayes if TRUE runs the MCMC chains in parallel for the bayesian model
#' @param parallel.dc if TRUE runs the MCMC chains in parallel for the DC model
#' @param inits list of initial values for the parameters of interest
#' @param priors list of priors for theta and kappa (Kparm):
#'               theta is assumed to have a beta prior with shape parameters theta.a and theta.b
#'               kappa is assumed to have a gamma prior with shape parameters K.shape and K.rate
#' @param n.cores number of cores to assign when computing in parallel
#' @param params parameters to monitor
#' @return list containing the vector of K clones, 
#'                         the bayesian object (GSSbayesfit)
#'                         the DC object (GSSdcfit)
#'                         an array with the changes in the marginal variance (clone.array.parm)
#'                         an array with the changes in the dominant eigenvalue (eigen.mat.parm)

IdentifiableGSSPopmodfun <- function(MCMClist = MCMClist, 
                                     data     = data, 
                                     Kvec     = Kvec,
                                     model    = model,
                                     inits    = inits,
                                     priors   = priors,
                                     n.cores  = 4,
                                     params   = c("Kparam","theta","sigma"),
                                     parallel.bayes  = F, parallel.dc = F){
  
  ##### ---------------------------------------------------------------
  ##### Model setup 
  ##### ---------------------------------------------------------------
  
  #### Load in the data
  Ymat         <- data$Ymat
  agematrix    <- data$agematrix
  turtlematrix <- data$turtlematrix
  lens         <- if(!is.numeric(data$lens)){ as.numeric(data$lens) }else{data$lens}
  nind         <- data$nind
  
  
  #### array[,,,] = row, column, array number
  K              <- 1
  m              <- dim(turtlematrix)[1]
  n              <- dim(turtlematrix)[2]
  dcYmat         <- array(Ymat, c(m,n,K))
  dcturtlematrix <- array(turtlematrix, c(m,n,K))
  dcagematrix    <- array(agematrix, c(m,n,K))
  dclens         <- lens
  
  #### Set up the model
  ssm.bayes.list <- append(list(Y       = dcYmat,
                                lens    = lens, 
                                ind     = nind, 
                                Xtrans  = dcagematrix,
                                Ytrans  = dcturtlematrix,
                                K       = K), priors)
  
  
  ##### ---------------------------------------------------------------
  ##### Fit the models
  ##### ---------------------------------------------------------------
  
  
  ### Load in the necessary packages
  require(foreach)
  require(parallel)
  require(doParallel)
  require(snow)
  require(snowfall)
  
  #### Make sure that even if parallel=T but nchain=1, parallel=F
  if(MCMClist$n.chains.bayes==1) parallel.bayes <- F
  if(MCMClist$n.chains.dc==1)    parallel.dc <- F
  
  #### Run a Bayesian model
  if(parallel.bayes == F){
    GSSbayesfit <- jags.fit(data     = ssm.bayes.list, 
                            params   = params, 
                            model    = model, 
                            inits    = inits,
                            n.chains = MCMClist$n.chains.bayes, 
                            n.adapt  = MCMClist$n.adapt, 
                            n.update = MCMClist$n.burn, 
                            n.iter   = MCMClist$n.iter, 
                            thin     = MCMClist$n.thin)
  }else{
    #### Set up the cluster
    cl.bayes <- parallel::makePSOCKcluster(n.cores)
    parallel::clusterEvalQ(cl.bayes, library(dclone))
    dclone::parLoadModule(cl.bayes,"lecuyer")
    
    #### Run the model
    GSSbayesfit <- dclone::jags.parfit(cl       = cl.bayes,
                                       data     = ssm.bayes.list, 
                                       params   = params, 
                                       model    = model, 
                                       inits    = inits,
                                       n.chains = MCMClist$n.chains.bayes, 
                                       n.adapt  = MCMClist$n.adapt, 
                                       n.update = MCMClist$n.burn, 
                                       n.iter   = MCMClist$n.iter, 
                                       thin     = MCMClist$n.thin)
    parallel::stopCluster(cl.bayes)
  }
  GSSbayesparams   <- summary(GSSbayesfit)[[1]][,1]
  
  
  #### Set up empty matrix to read data into
  clone.array.parm <- array(NA,
                            dim = c(length(params), 4, length(Kvec)), 
                            dimnames = list(params,
                                            c("Mean","SD","Naive SE","Time-series SE"),
                                            as.character(Kvec)))
  eigen.mat.parm   <- vector(length=length(Kvec))
  
  
  #### Variance covariance matrix and max eigen values
  var.postr  <- var(GSSbayesfit[[1]])
  max.eigen  <- max(eigen(var.postr)$values)
  NoClone.NI <- as.matrix(GSSbayesfit[[1]])
  
  
  ### prepare the foreach
  if( parallel.dc==F ){
    cl.K <- parallel::makeCluster(n.cores, outfile="")
  }else{
    cl.K <- parallel::makeCluster(n.cores/2, outfile="")
  }
  doParallel::registerDoParallel(cl.K)
  
  
  res <- foreach(i = 1:length(Kvec), 
                 .combine=rbind, 
                 .export = ls(.GlobalEnv)) %dopar% {
                   
                   #### Create a dummy list
                   GSSdclist <- list()
                   
                   #### extract the kth clone
                   K   <- Kvec[i]
                   
                   
                   #### Set up the model
                   dcYmat         <- array(Ymat, c(m,n,K))
                   dcturtlematrix <- array(turtlematrix, c(m,n,K))
                   dcagematrix    <- array(agematrix, c(m,n,K))
                   dclens         <- lens
                   
                   
                   #### Set up the model
                   ssm.dc.list <- append(list(Y       = dcYmat,
                                              lens    = dclens, 
                                              ind     = nind, 
                                              Xtrans  = dcagematrix,
                                              Ytrans  = dcturtlematrix,
                                              K       = K), priors)
                   
                   
                   #### Run the model
                   if(parallel.dc == F){
                     GSSdcfit <- dclone::jags.fit(data     = ssm.dc.list, 
                                                  params   = params, 
                                                  model    = model, 
                                                  inits    = inits,
                                                  n.chains = MCMClist$n.chains.dc,
                                                  n.adapt  = MCMClist$n.adapt,
                                                  n.update = MCMClist$n.burn,
                                                  n.iter   = MCMClist$n.iter,
                                                  thin     = MCMClist$n.thin)
                   }else{
                     
                     
                     cl.dc <- parallel::makePSOCKcluster(n.cores/2)
                     parallel::clusterEvalQ(cl.dc, library(dclone))
                     dclone::parLoadModule(cl.dc,"lecuyer")
                     
                     #### Run the model
                     GSSdcfit <- dclone::jags.parfit(cl       = cl.dc,
                                                     data     = ssm.dc.list,
                                                     params   = params,
                                                     model    = model,
                                                     inits    = inits,
                                                     n.chains = MCMClist$n.chains.dc,
                                                     n.adapt  = MCMClist$n.adapt,
                                                     n.update = MCMClist$n.burn,
                                                     n.iter   = MCMClist$n.iter,
                                                     thin     = MCMClist$n.thin)
                     parallel::stopCluster(cl.dc)
                   }
                   GSSdclist <- append(GSSdclist, GSSdcfit)
                   return( GSSdclist )
                 }
  parallel::stopCluster(cl.K) # Stop the parallel backend
  
  
  ### compute the relavant statistics to test for identifiability
  
  # Variance covariance matrix and max eigen values
  for( i in seq_along(Kvec)){
    clone.array.parm[,,i]   <- summary(res[[i]])[[1]]
    var.kth.post            <- var(res[[i]])
    eigen.mat.parm[i]     <- (eigen(var.kth.post)$values)/max.eigen # note here max.eigen is for no cloning data
  }
  
  #### compute mle and sd for the last K
  Klast     <- Kvec[length(Kvec)]
  GSSdcfit  <- res[[length(Kvec)]]
  GSSmle    <- summary(GSSdcfit)[["statistics"]][,"Mean"]
  GSSsdmle  <- sqrt(Klast)*(summary(GSSdcfit)[["statistics"]][,"SD"])
  
  
  ### Combine all MCMC chains if a list
  if( is.list(GSSdcfit) ){
    all.nodes <- do.call(rbind, GSSdcfit)
  }else{
    all.nodes <- GSSdcfit
  }
  
  #### Wald-type CI
  fish.inv           <- Klast*var(all.nodes)
  alpha              <- 0.05
  z.alpha.half       <- qnorm(p=(1-alpha/2))
  diagnonal.elements <- diag(fish.inv)
  mle.names          <- names(GSSmle)
  
  #### For K
  K.std.error <- z.alpha.half*sqrt(diagnonal.elements["Kparam"])
  GSSKci      <- c(GSSmle[mle.names=="Kparam"]-K.std.error,GSSmle[mle.names=="Kparam"]+K.std.error)
  
  #### For theta
  theta.std.error <- z.alpha.half*sqrt(diagnonal.elements["theta"])
  GSSThetaci      <- c(GSSmle[mle.names=="theta"]-theta.std.error,GSSmle[mle.names=="theta"]+theta.std.error)
  names(GSSKci)   <- names(GSSThetaci) <- c("lwrCI","upperCI")
  
  
  #### For sigma
  sigma.std.error  <- z.alpha.half*sqrt(diagnonal.elements["sigma"])
  GSSSigmaci        <- c(GSSmle[mle.names=="sigma"]-sigma.std.error,GSSmle[mle.names=="sigma"]+sigma.std.error)
  names(GSSSigmaci) <- c("lwrCI","upperCI")
  GSSci <- list(GSSKci=GSSKci,GSSThetaci=GSSThetaci,GSSSigmaci=GSSSigmaci)
  
  ##### save the dc objects into a list
  final.dc.obj <- list(K = Klast, GSSmle = GSSmle, 
                       GSSsdmle = GSSsdmle, GSSci = GSSci, 
                       model.file = GSSdcfit, 
                       MCMCsamples = all.nodes)
  
  #### return a list
  return(invisible(list(fullidentobj=res, Kvec=Kvec,
                        GSSbayesfit=GSSbayesfit,
                        GSSdcfit=final.dc.obj,
                        clone.array.parm=clone.array.parm,
                        eigen.mat.parm=eigen.mat.parm)))
}




##### ---------------------------------------------------------------
##### Plots
##### ---------------------------------------------------------------

#' Arguments
#' @param data list containing the output from the IdentifiableGSSPopmodfun function
#' @param Kvec specifying a vecctor of clones for data cloning to iterate over
#' @param params parameters to plot

plot.identifiable.pop <- function(data=data, Kvec=Kvec, cex=15,
                                  params=c("Kparam","theta","sigma")){
  require(ggpubr)
  
  ##### ---------------------------------------------------------------
  ##### Prepare data to be plotted
  ##### ---------------------------------------------------------------
  GSSbayesfit      <- data$GSSbayesfit
  GSSdcfit         <- data$GSSdcfit
  clone.array.parm <- data$clone.array.parm
  eigen.mat.parm   <- data$eigen.mat.parm
  
  
  #### melting separate No clone
  NoClone.NI.parm       <- as.matrix((GSSbayesfit)[[1]])
  melted.nocloneNI.parm <- reshape2::melt(NoClone.NI.parm)
  vline.parm            <- summarise(group_by(melted.nocloneNI.parm,Var2), mean = mean(value))
  
  #### Last clone
  Clone.NI.parm       <- as.matrix(GSSdcfit$model.file) # This should have the information about the last clone 512
  melted.cloneNI.parm <- reshape2::melt(Clone.NI.parm)
  vline.ni.parm      <- summarise(group_by(melted.cloneNI.parm,Var2), mean = mean(value))
  
  ## Accessing the mean from clone.array for all clones
  mean.mat.parm<-matrix(NA,
                        nrow = length(Kvec),ncol = length(params),
                        dimnames = list(as.character(Kvec),params))
  for(i in 1:length(Kvec)){
    mean.mat.parm[i,]<-t(clone.array.parm[,1,i])
  }
  
  ## Accessing standard deviation from the clone.array for all clones
  sd.mat.parm<-matrix(NA,
                      nrow = length(Kvec), ncol = length(params),
                      dimnames = list(as.character(Kvec),params))
  for (i in 1:length(Kvec)){
    sd.mat.parm[i,]<-t(clone.array.parm[,2,i])
  }
  
  # Getting variances and plot
  var.mat.parm    <- sd.mat.parm**2 # square of sd 
  melted.var.parm <- reshape2::melt(var.mat.parm)
  dat.eigen.parm  <- data.frame(Kvec,eigen.mat.parm)
  colnames(melted.var.parm) <- c("Var1" , "Parameters" ,"value")
  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  # cbbPalette <- viridis::plasma(5, begin = 0, end = 0.8)
  favcol     <- c("skyblue4","aquamarine4")
  
  new.labs <- c( as.character(expression(kappa)), as.character(expression(sigma)), 
                 as.character(expression(theta)))
  levels(melted.var.parm$Parameters) <- levels(melted.cloneNI.parm$Var2) <- levels(melted.nocloneNI.parm$Var2) <- new.labs
  
  theme_cex <- theme(
    plot.title = element_text(size = 20),              # Main title
    axis.title.x = element_text(size = 20),            # X-axis label
    axis.title.y = element_text(size = 20),            # Y-axis label
    strip.text = element_text(size = 20)               # Facet labels
  )
  
  # Function to create exactly three pretty-spaced ticks
  breaks_exact <- function(limits, size=3) {
    pretty_breaks <- pretty( seq(limits[1], limits[2], length.out = size)  ) 
    return(round(pretty_breaks,3))
  }
  
  
  # Figure 1A no clone
  ni.dc1 <- ggplot(melted.nocloneNI.parm, aes(x=value)) + 
    geom_density(alpha=.75, fill="#0072B2") + xlab("") +
    facet_wrap(~ Var2, scales = "free", labeller = label_parsed) +
    ylab(" Density")+ ggtitle("Bayesian analysis (no clones)")+
    theme_pubr(base_size = cex) + 
    theme_cex + scale_x_continuous(breaks = breaks_exact)
  
  # Figure 1B clone 128 
  ni.dc2 <- ggplot(melted.cloneNI.parm, aes(x=value)) + 
    geom_density(alpha=.75,fill="#0072B2")+ xlab("Parameter value") +
    facet_wrap(~ Var2, scales = "free",labeller = label_parsed) +
    ggtitle(paste0( max(Kvec), " clones" ))+ylab("Density")+
    theme_pubr(base_size = cex) + 
    theme_cex + scale_x_continuous(breaks = breaks_exact)
  
  # Figure 1C Marginal variances
  math.labs <- c( expression(kappa), 
                  expression(sigma), 
                  expression(theta))
  ni.dc3<-ggplot(data=melted.var.parm, 
                 aes(x=Var1,y=value,group=Parameters,colour=Parameters,shape=Parameters)) +
    geom_line(size=1) + geom_point(size=2)+
    scale_shape_discrete(labels=math.labs) + 
    scale_colour_manual(values=cbbPalette, labels=math.labs) +
    ggtitle("Marginal variance")+ xlab("Number of clones")+ ylab("Marginal variance ")+
    theme_pubr(base_size = cex) + 
    theme(
      legend.text = element_text(size = 20),
      plot.title = element_text(size = 20),              # Main title
      axis.title.x = element_text(size = 20),            # X-axis label
      axis.title.y = element_text(size = 20),            # Y-axis label
      strip.text = element_text(size = 20)               # Facet labels
    )
  
  # Figure 2D Eigen values
  ni.dc4 <- ggplot(data=dat.eigen.parm, aes(x=Kvec,y=eigen.mat.parm)) +
    geom_line(color="aquamarine4",size=1) + geom_point(color="aquamarine4",size=3)+
    geom_line(aes(x=Kvec,y=1/Kvec), color="grey80",size=1)+
    ggtitle("Eigenvalues")+ xlab("Number of clones")+ ylab("Eigenvalues ")+ 
    theme_pubr(base_size = cex) # + scale_x_continuous(breaks = dat.eigen.parm$Kvec)
  
  figure1 <- ggarrange(ni.dc1,ni.dc4,ni.dc2,ni.dc3,
                       labels = c("(A)", "(B)", "(C)", "(D)"),
                       font.label = list(face = "plain"),
                       nrow = 2,ncol = 2,
                       heights = c(1, 1),  # Adjust height of each row; second row slightly smaller
                       widths = c(1, 0.5))    # Adjust width of each column; second column less wide)
  
  print(figure1) 
}







######## Mortality and survival functions
Gompertz.Mortality <- function(t=t, b=b, a=NULL, c=NULL, shape="simple" ){ 
  b0 <- b[1]; b1 <- b[2]
  if(shape=="simple")    haz <- exp( b0+b1*t ) ### a * exp( b * t ) 
  if(shape=="makeham")   haz <- c + exp( b0+b1*t ) 
  if(shape=="bathtub"){ 
    a0 <- a[1]; a1 <- a[2]   
    haz <- exp( a0-a1*t ) + c + exp( b0+b1*t ) 
  }
  return(haz)
}
Gompertz.Survival  <- function( t=t, b=b, a=NULL, c=NULL, shape="simple" ){ 
  b0 <- b[1]; b1 <- b[2]
  if(shape=="simple") { surv <- exp((exp(b0)/b1)*(1-exp(b1*t))) }
  if(shape=="makeham"){ surv <- exp(-c*t + (exp(b0)/b1)*(1-exp(b1*t))) }
  if(shape=="bathtub"){
    a0 <- a[1]; a1 <- a[2]
    surv <- exp(exp(a0)/a1*(exp(-a1*t)-1)-c*t+(exp(b0)/b1)*(1-exp(b1*t)))
  }
  return(surv)
}
SE <- function(x){ sd(x)/sqrt(length(x)) }



##### ---------------------------------------------------------------
##### Fitted model/ future predictions functions
##### ---------------------------------------------------------------

#' Arguments
#' @param K specifying the number of clones for data cloning
#' @param model JAGS model 
#' @param MCMClist a list contain the MCMC setting for JAGS
#' @param data list containing matrices of cumulative reproductive output per female (Ymat),
#'             annual reproductive output (turtlematrix), and 
#'             waiting times between reproductive seasons (agematrix)
#' @param rawdata data frame of the reproductive CTC data
#' @param timemeasure either BSEA (breeding season) or BAGE (breeding age)
#' @param my.cols colour of the individual tracectors per female if add.ind.pred = TRUE
#' @param nsims number of simulations to compute the average trajectory
#' @param bayes fits a bayesian model rather than a DC model if TRUE
#' @param add.ind.pred draws predicted trajectories for each female when TRUE
#' @param add.traj draws each iteration of nsims if TRUE
#' @param fitted.model object fitted with GSSPopmodfun()
#' @param tau the waiting time between reproductive season scaled to BAGE needed to plot global trend
#' @param plot if TRUE returns MCMC diagnostic plots for the parameters (params) of interest
#' @param parallel if TRUE runs the MCMC chains in parallel
#' @param inits list of initial values for the parameters of interest
#' @param priors list of priors for theta and kappa (Kparm):
#'               theta is assumed to have a beta prior with shape parameters theta.a and theta.b
#' @return plot of the population model

pop.predict.trajectory <- function(K = K, models = models, MCMClist = MCMClist,
                                   data = data, rawdata = rawdata,
                                   timemeasure = c("BSEA","BAGE"),
                                   my.cols="darkmagenta", # c("darkmagenta", "firebrick2", "darkorange")
                                   nsims = 5000, alpha = 0.005, plot.range=TRUE,
                                   bayes = T, add.ind.pred = F, add.traj = F, 
                                   parallel = T, include.legend=F, pNfuture=NULL,
                                   fitted.model = NULL, tau = NULL, pred.future =T,
                                   priors = priors, inits = NULL
){
  if(!any(.packages() %in% "scales")) require(scales)
  if(!any(.packages() %in% "dclone")) require(dclone)
  if(!any(.packages() %in% "magrittr")) require(magrittr)
  
  if(!is.numeric(data$lens)) data$lens <- as.numeric(data$lens)
  if(is.null(fitted.model)){
    #### Make sure that even if parallel=T but nchain=1, parallel=F
    if(MCMClist$n.chains.dc==1) parallel <- F
    
    fitted.model <- GSSPopmodfun(MCMClist = MCMClist, 
                                 data     = data, 
                                 K        = K, 
                                 inits    = inits,
                                 model    = models$fit, 
                                 parallel = parallel)
  }
  
  Kparam    <- fitted.model$GSSmle[['Kparam']]
  theta     <- fitted.model$GSSmle[['theta']]
  sigma     <- fitted.model$GSSmle[['sigma']]
  sigsq.inv <- 1/(sigma^2)
  
  ssm.list <- list(Y         = data$Ymat,
                   lens      = data$lens, 
                   ind       = data$nind, 
                   Xtrans    = data$agematrix,
                   Ytrans    = data$turtlematrix,
                   Kparam    = Kparam,
                   theta     = theta,
                   sigma     = sigma)
  
  kalman.pred <- dclone::jags.fit(data     = ssm.list, 
                                  params   = "phi", 
                                  model    = models$pred, 
                                  n.chains = MCMClist$n.chains.dc, 
                                  n.adapt  = MCMClist$n.adapt, 
                                  n.update = MCMClist$n.burn, 
                                  n.iter   = MCMClist$n.iter, 
                                  thin     = MCMClist$n.thin)
  
  
  #### ----------------------------------------------------------------------------------------
  #### Predictions
  #### ----------------------------------------------------------------------------------------
  
  ### time measures: BSEA, BAGE, or absolute age
  if(timemeasure=="BSEA"){ 
    Time    <- rawdata$BSEA
    tempmat <- data$Xmat
    tempmat[!is.na(tempmat)] <- 1
    Timemat  <- t(apply(tempmat, 1, function(x) cumsum(x) ))
  } else if(timemeasure=="BAGE"){ 
    if(any(colnames(rawdata) %in% "BAGE")) colnames(rawdata)[which(colnames(rawdata) %in% "BAGE")] <- "AGE"
    Time     <- rawdata$AGE
    Timemat  <- data$Xmat
  }else if(timemeasure=="Age"){
    if(any(colnames(rawdata) %in% "Age")) colnames(rawdata)[which(colnames(rawdata) %in% "Age")] <- "AGE"
    # Time     <- rawdata$AGE
    # Timemat  <- data$trueagematrix
    Timelist <- list()
    ifelse(min(rawdata$AGE)==0, rawdata$AGE <- rawdata$AGE+1, rawdata$AGE <- rawdata$AGE)
    Timelist[paste0("Time")]     <- list(rawdata$AGE)
    Timelist[paste0("Timemat")]  <- list(data[["trueagematrix"]])
    attach(Timelist)
  }else{ print("No suitable time measure specified") }
  
  
  #### Individual trajectories
  PAR        <- summary(kalman.pred)$statistics
  CI         <- summary(kalman.pred)$quantiles 
  Ypred      <- as.vector(PAR[,"Mean"])
  cumoutList <- c()
  
  for(i in 1:length(lens) ){
    init.out = data$Ymat[,1][i]
    mat.temp = data.frame(ID         =i,
                          Times      = na.omit(as.vector(Timemat[i,])), 
                          MeanCUMOUT = cumsum(c(init.out, PAR[gsub(".*?([0-9]+).*", "\\1", row.names(PAR)) == i, 1])),
                          LowerCI    = cumsum(c(init.out, CI[gsub(".*?([0-9]+).*", "\\1", row.names(CI)) == i, 1])),
                          UpperCI    = cumsum(c(init.out, CI[gsub(".*?([0-9]+).*", "\\1", row.names(CI)) == i, 5])))
    cumoutList <- c(cumoutList, list(mat.temp))
  }
  cumoutDataframe <- do.call(rbind, cumoutList) |> as.data.frame()
  
  #### Define the jumps between reproductive episodes
  if(is.null(tau)){
    Xtrans <- 1 # rep(1, (max(data$Xmat, na.rm = T)))
  }else{ 
    Xtrans <- tau #rep(tau, (max(data$Xmat, na.rm = T))) 
  }
  
  #### initial reproductive output
  initAGE   <- 1
  Y1        <- mean(data$Ymat[,'1']) #### use the average number of eggs for first years
  minAGE    <- min(cumoutDataframe$Times)
  maxAGE    <- max(cumoutDataframe$Times) # length(unique(cumoutDataframe$Times)) ### max number of unique "ages" = 36
  rangeAGE  <- (maxAGE-minAGE)+1 ## number of ages with reproductive data
  
  #### ----------------------------------------------------------------------
  #### Bootstrap CI
  pop.mat      <- c()
  det.jump.mat <- c()
  
  ## number of years to predict into the future: if NULL = 0
  pNfuture <- ifelse(is.null(pNfuture), 0, pNfuture) 
  pNfuture <- ifelse(pred.future==FALSE, 0, pNfuture) 
  
  for (j in 1:nsims){
    
    pop.mat.temp      <- c(Y1, rep(NA, rangeAGE-1) ) 
    det.jump.mat.temp <- c(Y1, rep(NA, rangeAGE-1) ) 
    temp.rangeAGE     <- rangeAGE+pNfuture # max time plus future years
    
    #### Inner loop is stepping across years
    for (t in (initAGE+1):temp.rangeAGE){ ### For i in 2:36 (BAGE) :19 (BSEA)
      
      #### Make sure the trend does not get greater than the carrying capacity
      if(pop.mat.temp[t-1] > Kparam) break
      
      #### Draw values of the stochastic distributions  
      #### Project population size one time step with those values of a
      
      #### difference between the two means
      det.jump             <- (Kparam * exp( log(pop.mat.temp[(t-1)]/Kparam) * 
                                               exp(-theta * Xtrans))) - pop.mat.temp[(t-1)] 
      mu                   <- log(det.jump) - (1/sigsq.inv * 0.5)
      phi                  <- rlnorm(n=1, meanlog=mu, sdlog=sigma)
      pop.mat.temp[t]      <- pop.mat.temp[t-1] + phi
      det.jump.mat.temp[t] <- det.jump
      ### Save in lists
    } # close t
    
    ## make sure that lengths of pop.mat.temp and det.jump.mat.temp are all the same
    if( length(pop.mat.temp) < (temp.rangeAGE) ){
      lenNAs       <- (temp.rangeAGE) - length(pop.mat.temp)
      pop.mat.temp      <- c( pop.mat.temp, rep(NA,lenNAs) )
      det.jump.mat.temp <- c( det.jump.mat.temp, rep(NA,lenNAs) )
    }
    
    pop.mat      <- cbind(pop.mat, pop.mat.temp)
    det.jump.mat <- cbind(det.jump.mat, det.jump.mat.temp)
    
  }# close j
  
  
  #### ----------------------------------------------------------------------
  NfullList  <- tibble::data_frame(iter   = rep(1:nsims, each = temp.rangeAGE),
                                   bsea   = rep(1:temp.rangeAGE, nsims),
                                   class  = rep( c(rep("Observed", (initAGE)),
                                                   rep("predicted", (temp.rangeAGE-1) )), nsims ),
                                   mat    = as.vector( unlist(qpcR:::cbind.na(pop.mat)) ) ) %>% as.data.frame()
  NfullList$mat <- ifelse(NfullList$mat>Kparam, Kparam, NfullList$mat)
  pop.mat       <- apply(pop.mat,2,FUN=function(x) ifelse(x > Kparam, Kparam, x))
  PRED          <- apply(pop.mat, 1,  
                         FUN=function(x){ quantile(x, prob=c(0.025, 0.5, 0.975), na.rm=T  ) })|> t()
  PRED          <- data.frame( type=NA, lci=PRED[,1], mu=PRED[,2], uci=PRED[,3])
  PRED$type[1]        <- "observed"
  PRED$type[2:rangeAGE] <- "predicted"
  if( pred.future==T & pNfuture>0 ) PRED$type[(rangeAGE+1):temp.rangeAGE] <- "future_trajectory"
  
  
  #### ----------------------------------------------------------------------
  X1         <- PRED['50%',]
  LCI        <- PRED['2.5%',]
  UCI        <- PRED['97.5%',]
  
  
  #### ----------------------------------------------------------------------
  #### Model plots
  #### ----------------------------------------------------------------------
  x          <- Time
  y          <- rawdata$cumeggs
  pch        <- rawdata$estCS
  
  #### ----------------------------------------------------------------------
  #### Reproductive model plot
  #### ----------------------------------------------------------------------
  par(mfrow=c(1,1),mar=c(5,6,2,4))
  
  ### x-axus labels for age
  if(timemeasure=="Age"){
    # minAGE      <- min(rawdata$AGE, na.rm = T)
    pred.xs     <- seq(1, temp.rangeAGE, by=1) ## 7 is the min age for a reproductive turtle
    plot.maxAGE <- temp.rangeAGE+(minAGE-1)
    xlab <- paste0( timemeasure, " (years)" )
  }else{
    plot.maxAGE <- temp.rangeAGE
    if( timemeasure=="BSEA" ) pred.xs <- ((initAGE):temp.rangeAGE)
    if( timemeasure=="BAGE" ) pred.xs <- ((initAGE):temp.rangeAGE)-1
    xlab <- timemeasure
  }
  
  plot(x, rawdata$cumeggs,
       xlim=c(min(Timemat, na.rm = T), 
              ifelse( ifelse(plot.range==TRUE & pNfuture==1, TRUE, FALSE),
                      max(Timemat, na.rm = T),
                      max( max(Timemat, na.rm = T), plot.maxAGE, na.rm = T)) ), ## 7 is the min age for a reproductive turtle
       ylim=c(0, max(data$Ymat, Kparam, na.rm = T)), 
       type="p",  col=alpha(my.cols, alpha=0.5), cex=2,
       cex.axis=1.5, cex.lab=1.5, cex.main=2,
       xlab=xlab, ylab="Cumulative reproductive output")
  for(i in 1:length(data$lens)){
    lines(as.vector(Timemat[i,]), 
          ts( as.vector( t(data$Ymat[i,]) ) ), 
          type="b", pch=21, bg=alpha(my.cols, alpha=0.9), cex=2.5)
  }
  
  if(add.ind.pred==T){
    for(i in 1:length( data$lens ) ){
      points( cumoutDataframe$Times[cumoutDataframe$ID==i], cumoutDataframe$MeanCUMOUT[cumoutDataframe$ID==i],  type="l", pch=16, col="darkred", cex=1.1)
      lines( cumoutDataframe$Times[cumoutDataframe$ID==i], cumoutDataframe$LowerCI[cumoutDataframe$ID==i],  type="l", lty="dashed", pch=16, col="red3", cex=1.1)
      lines( cumoutDataframe$Times[cumoutDataframe$ID==i], cumoutDataframe$UpperCI[cumoutDataframe$ID==i],  type="l", lty="dashed", pch=16, col="red3", cex=1.1)
    }
  }
  
  #### ----------------------------------------------------------------------
  ### mean and 95% quantiles of the trend
  maxAGEPlot <- max(cumoutDataframe$Times)
  if(timemeasure=="Age"){
    lines(pred.xs[1:maxAGEPlot]+(minAGE-1), PRED[1:maxAGEPlot,3], type="l", pch=21, lwd=4, cex=2, bg = my.cols)
    lines(pred.xs[1:maxAGEPlot]+(minAGE-1), PRED[1:maxAGEPlot,2], lty=2, lwd=2, bg = my.cols)
    lines(pred.xs[1:maxAGEPlot]+(minAGE-1), PRED[1:maxAGEPlot,4], lty=2, lwd=2, bg = my.cols)
    abline(h=Kparam, lty="dashed", lwd=2.5)
  }else{
    lines(pred.xs[1:maxAGEPlot], PRED[1:maxAGEPlot,3], type="l", pch=21, lwd=4, cex=2, bg = my.cols)
    lines(pred.xs[1:maxAGEPlot], PRED[1:maxAGEPlot,2], lty=2, lwd=2, bg = my.cols)
    lines(pred.xs[1:maxAGEPlot], PRED[1:maxAGEPlot,4], lty=2, lwd=2, bg = my.cols)
    abline(h=Kparam, lty="dashed", lwd=2.5)
  }
  
  ### Trajectory +N years
  if(pred.future==T & pNfuture>0){
    traj.xs <- pred.xs[rangeAGE:temp.rangeAGE] # prediction
    
    # add a vertical line to separate trajectories (past last age of observed reproduction)
    segments( x0=max(rangeAGE)+(minAGE-1), x1=max(rangeAGE)+(minAGE-1), 
              y0=PRED[max(rangeAGE), 2], y1=PRED[max(rangeAGE), 4], 
              ## -6 because line 765: +6 because 7 is the min age for a reproductive turtle
              lty=1, lwd=2, col = "#525064FF")
    
    x.ribbon <- c(traj.xs, rev(traj.xs))+(minAGE-1)
    y.ribbon <- c(PRED[(rangeAGE):temp.rangeAGE,2], 
                  rev(PRED[(rangeAGE):temp.rangeAGE,4]))
    polygon(x.ribbon, y.ribbon,
            col = adjustcolor("#525064FF", 0.3), border = NA)
    
    lines(traj.xs+(minAGE-1), PRED[(rangeAGE):temp.rangeAGE,3], type="l", pch=21, lwd=4, cex=2, lty=1) # prediction
    lines(traj.xs+(minAGE-1), PRED[(rangeAGE):temp.rangeAGE,2], type="l", pch=21, lwd=2, cex=2, lty=2) # lower CI
    lines(traj.xs+(minAGE-1), PRED[(rangeAGE):temp.rangeAGE,4], type="l", pch=21, lwd=2, cex=2, lty=2) # upper CI
    
  }
  
  if(add.traj==T){
    for(i in 1:nsims){ 
      plot.data <- NfullList[which(NfullList$iter == i),] 
      lines((((initAGE):temp.rangeAGE)+(minAGE-1)), 
            plot.data$mat[(initAGE):temp.rangeAGE], 
            col = scales::alpha("grey30", alpha=alpha), lwd = 2)
    }
  }
  if(include.legend==T){
    legend(locator(1), 
           legend   = c("Observed cumulative output",
                        "Predicted cumulative output"),
           col      = c(my.cols,"black"),
           pch      = c(16, NA), 
           lty      = c(NA, "solid"),
           bty      = "n", pt.cex = 2, cex = 0.9, 
           text.col = "black", horiz = F, inset = c(0.1, 0.1)) 
  }
  out <- list(cumoutDataframe=cumoutDataframe, 
              NfullList=NfullList, PRED=PRED)
  return( out )
}


