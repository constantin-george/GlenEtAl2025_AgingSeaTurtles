#### ---------------------------------------------------------------------------
#### phylogenetic generalized least squares: SFLC aging rate
#### Author: C. George Glen and Maria Torres Sanchez 
#### Date: 12/12/2023
#### ---------------------------------------------------------------------------

rm(list=ls()) # remove everything currently held in the R memory
graphics.off() # close all open graphics window
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


#### ---------------------------------------------------------------------------
#### Preamble
#### ---------------------------------------------------------------------------


#### Libraries
load.lib <- function(){
  library(phytools)
  library(ape)
  library(car)
  library(nlme)
  library(geiger)
  library(ggplot2)
  library(ggtree)
  library(readODS)
  library(phylobase)
  library(adephylo)
  library(phylolm)
  library(caper)
  library(bbmle)
  library(MuMIn)
  library(phylolm) 
  library(tidyverse)
  library(emmeans)
  library(multcomp)
  library(patchwork)
  library(rr2)
}; load.lib()



## load R data and tree
data.path <- "/Users/constantingeorgeglen/Documents/University/UFL/RESEARCH/Projects_Analysis/RSenescence_CTC/Main/2.CTCAnalyses/PhylogeneticAnalysis/scripts/20240516/PGLS/"
load( paste0( data.path, "/submissionclean/model.data.RData") )
model.data$lnAnnualRepOut <- log(model.data$AnnualRepOut)



#### ---------------------------------------------------------------------------
#### Slow-fast continuum hypothesis from Reinke et al 2022
#### Hypothesis: aging rate is inversely related to annual fecundity and 
####             age at first reproduction
#### ---------------------------------------------------------------------------


### model formula
model.forms <- list( 
  "mFULLfixed" = formula(Mean.AR.Female ~ (lnAnnualRepOut + lnMASS + lnafr) + status)
) 


### fit models
fit.NULL <- gls(model.forms$mFULLfixed, correlation = NULL, method="ML",
                control=glsControl(opt="nlminb"), data = model.data, na.action = "na.omit")
fit.OU <- phylolm(model.forms$mFULLfixed, data = model.data, upper.bound = 1e50,
                  model = "OUrandomRoot", phy=model.tree, boot =1000,
                  measurement_error = F)
fit.lam <- phylolm(model.forms$mFULLfixed, data = model.data, lower.bound = 1e-20,
                   model = "lambda", phy=model.tree, boot =1000)
fit.BM <- phylolm(model.forms$mFULLfixed, data = model.data,
                  model = "BM", phy=model.tree, boot =1000)
fit.BK <- gls(model.forms$mFULLfixed, data = model.data,
              correlation = corBlomberg(0.1, phy = model.tree, 
                                        form = ~species.num, fixed = T),
              method="ML", control=glsControl(opt="nlminb") )

# IC comparison
options(digits=5)
ICres <- ICtab( fit.OU, fit.lam, fit.BK, fit.BM, fit.NULL, type="BIC",
                weights = T, delta = T, base = T, logLik = T)
ICres
#          logLik BIC    dLogLik dBIC   df weight
# fit.NULL  117.0 -209.0  322.2     0.0 6  0.8   
# fit.OU    117.0 -204.9  322.2     4.2 7  0.1   
# fit.lam   117.0 -204.9  322.2     4.2 7  0.1   
# fit.BK   -173.8  372.5   31.4   581.5 6  <0.001
# fit.BM   -205.2  435.4    0.0   644.4 6  <0.001


### Calculate partial and total R2s
rr2::R2_lik(fit.OU); rr2::R2_lik(fit.lam); 
# [1] 0.096248
# [1] 0.096248


#### ---------------------------------------------------------------------------
### pick the best model
best.mod <- fit.NULL



#### ---------------------------------------------------------------------------------
#### inference best model
performance::check_model(fit.NULL)

summary(best.mod)
vif(best.mod)


### check cor and cov-var matrices
# best.mod$varBeta; # covariance
anova(best.mod)
multcomp::glht(best.mod) |> summary()

### confidence intervals
confint(best.mod)
# intervals(best.mod)
intervals(best.mod, level=0.95, which="coef") 


## create a table of model results
sjPlot::tab_model(best.mod, 
                  show.intercept = T,show.stat = T, show.se = T, 
                  show.dev = T, show.loglik = T, 
                  show.icc =T, show.ngroups = T, show.re.var = T, p.style = "numeric",
                  collapse.ci = T, vcov.type = "HC0", digits = 2, p.val = "wald")




# diagnostics
hist(best.mod$residuals, breaks = 20) #check for normality of the residuals
qqnorm(best.mod$residuals) #check for normality of the residuals
qqline(best.mod$residuals) #check for normality of the residuals




#### ---------------------------------------------------------------------------
#### Plots
#### ---------------------------------------------------------------------------
library(viridis)
library(colorspace)

# load Reinke 2022 predictions
load("../reinke2022ReptilesALLPredDAT.RData")
load("../reinke2022ReptilesWOTPredDAT.RData")

# my.col <- c("#377EB8","#FF7F00"); 
my.col.default <- viridis(2, begin=0.2, end=0.6, option = "magma")
my.col.class <- c( my.col.default, rep("darkgreen",2) )
# my.col <- c("royalblue3","forestgreen")

pch <- c(21,22,23,24); cex <- 8; base_size <- 25
my.col.lightened <- colorspace::lighten(my.col.class, 0.3)
my.col.scale <- scales::alpha( 
  #c("grey20","grey30"), 
  my.col.class,
  1
)

## function to log transform if needed
fn   <- function( x, logT=F ) return(if( logT==T){log(x)}else{x})
logT <- F

# lnAnnualRepOut holding everything else at their mean value
Niter    <- 1000
parMLE   <- coef(best.mod)
new.data <- by(model.data, model.data$status, function(x){
  x.data <- data.frame( lnAnnualRepOut=seq(min(x$lnAnnualRepOut), max(x$lnAnnualRepOut), length.out=Niter),
                        lnafr=mean(x$lnafr),
                        lnMASS=mean(x$lnMASS))
  x.data$status <- unique(x$status)
  return(x.data)
} )
new.data <- do.call(rbind, new.data)

fit.mat <- matrix(predict(best.mod, new.data, se=T)$fit, nrow=Niter) |> as.data.frame()
se.mat  <- matrix(predict(best.mod, new.data, se=T)$se.fit, nrow=Niter)|> as.data.frame()
new.data$fit    <- tidyr::gather(fit.mat, value = "fit")[,2]
new.data$fit.se <- tidyr::gather(se.mat, value = "fit.se")[,2]

# add a label for testudines in class
model.data$class <- "testudines"

# main plot 1
pg1 <- ggplot(model.data ) + 
  geom_point(aes( x=lnAnnualRepOut, y=fn(Mean.AR.Female, logT=logT), shape=status, 
                  group=status, fill=status), 
             size=cex, col="black") + 
  geom_ribbon( data=new.data, aes( x=lnAnnualRepOut, 
                                   ymin = fn(fit - fit.se, logT=logT), 
                                   ymax = fn(fit + fit.se, logT=logT), 
                                   group=status, col=status, fill=status), 
               alpha=0.5, linetype="longdash")+
  geom_line( data=new.data, aes( x=lnAnnualRepOut, y=fn(fit, logT=logT), group=status), 
             col="black", lwd=1.5 ) +
  
  # Reinke 2022 predictions
  geom_point( data=ectodata %>% filter(Group %in% "Reptiles" & 
                                         !(Type %in% "Turtles" | Type %in% "Squamates")),
              aes(x=lnAF, y=GompertzSlope, shape = Type), size=cex, fill="darkgreen") +
  geom_line( data=reinke2022reptilesALL_AF, aes(x=lnAF, y=fn(fit, logT=logT)), 
             lwd=1, linetype="dashed", col="darkgreen") +
  geom_line( data=reinke2022reptilesWOT_AF, aes(x=lnAF, y=fn(fit, logT=logT)), 
             lwd=1, linetype="solid", col="darkgreen") +
  
  
  # plot settings
  scale_y_continuous(breaks = fn(pretty( seq(min(model.data$Mean.AR.Female), 
                                             max(model.data$Mean.AR.Female, 
                                                 reinke2022reptilesWOT_AF$fit), 
                                             length.out=8)), logT=logT )) + 
  scale_shape_manual( values=pch, labels=c("Testudines: Captive",
                                           "Testudines: Wild",
                                           "Crocodilians","Tuataras") ) + 
  scale_color_manual( values = my.col.default) + 
  scale_fill_manual( values = my.col.lightened) +
  ggpubr::theme_pubr(base_size = base_size) + 
  theme(legend.text = element_text(size = 12)) + 
  labs(x="ln(Annual fecundity)", y="") +
  guides(shape=guide_legend(title = "", 
                            override.aes = list( fill = my.col.lightened)), 
         col="none", fill = "none") 
pg1


### lnafr
new.data <- by(model.data, model.data$status, function(x){
  x.data <- data.frame( lnafr=seq(min(x$lnafr), max(x$lnafr), length.out=Niter),
                        lnAnnualRepOut=mean(x$lnAnnualRepOut),
                        lnMASS=mean(x$lnMASS))
  x.data$status <- unique(x$status)
  return(x.data)
} )
new.data <- do.call(rbind, new.data)

fit.mat <- matrix(predict(best.mod, new.data, se=T)$fit, nrow=100) |> as.data.frame()
se.mat  <- matrix(predict(best.mod, new.data, se=T)$se.fit, nrow=100)|> as.data.frame()
new.data$fit    <- tidyr::gather(fit.mat, value = "fit")[,2]
new.data$fit.se <- tidyr::gather(se.mat, value = "fit.se")[,2]

pg2 <- ggplot(model.data ) + 
  geom_point(aes( x=lnafr, y=Mean.AR.Female, shape=status, 
                  group=status, fill=status), 
             size=cex, col="black") + 
  geom_ribbon( data=new.data, aes( x=lnafr, ymin = fit - fit.se, 
                                   ymax = fit + fit.se, 
                                   group=status, col=status, fill=status), 
               alpha=0.5, linetype="longdash")+
  geom_line( data=new.data, aes( x=lnafr, y=fit, 
                                 group=status), 
             col="black", lwd=1.5 ) +
  # Reinke 2022 predictions
  geom_point( data=ectodata %>% filter(Group %in% "Reptiles" & 
                                         !(Type %in% "Turtles" | Type %in% "Squamates")),
              aes(x=lnAFR, y=GompertzSlope, shape = Type), size=cex, fill="darkgreen") +
  geom_line( data=reinke2022reptilesALL_AFR, aes(x=lnAFR, y=fit), 
             lwd=1, linetype="dashed", col="darkgreen") +
  geom_line( data=reinke2022reptilesWOT_AFR, aes(x=lnAFR, y=fit), 
             lwd=1, linetype="solid", col="darkgreen") +
  
  # plot settings
  scale_y_continuous(breaks = pretty( seq(min(model.data$Mean.AR.Female), 
                                          max(model.data$Mean.AR.Female, 0.5), 
                                          length.out=8) ), 
                     limits = c(-0.04, 0.5) ) + 
  scale_shape_manual( values=pch, labels=c("Testudines: Captive",
                                           "Testudines: Wild",
                                           "Crocodilians","Tuataras") ) + 
  scale_color_manual( values = my.col.default) + 
  scale_fill_manual( values = my.col.lightened) +
  ggpubr::theme_pubr(base_size = base_size)+
  labs(x="ln(AFR)", y="Aging rate") + theme(legend.text = element_text(size = 12)) + 
  lims(x=range(model.data$lnafr))+
  guides(shape=guide_legend(title = "", 
                            override.aes = list( fill = my.col.lightened)), 
         col="none", fill = "none")
pg2


### body mass
new.data <- by(model.data, model.data$status, function(x){
  x.data <- data.frame( lnMASS=seq(min(x$lnMASS), max(x$lnMASS), length.out=Niter),
                        lnAnnualRepOut=mean(x$lnAnnualRepOut),
                        lnafr=mean(x$lnafr))
  x.data$status <- unique(x$status)
  return(x.data)
} )
new.data <- do.call(rbind, new.data)

fit.mat <- matrix(predict(best.mod, new.data, se=T)$fit, nrow=100) |> as.data.frame()
se.mat  <- matrix(predict(best.mod, new.data, se=T)$se.fit, nrow=100)|> as.data.frame()
new.data$fit    <- tidyr::gather(fit.mat, value = "fit")[,2]
new.data$fit.se <- tidyr::gather(se.mat, value = "fit.se")[,2]

pg3 <- ggplot(model.data ) + 
  geom_point(aes( x=lnMASS, y=Mean.AR.Female, shape=status, 
                  group=status, fill=status), 
             size=cex, col="black") + 
  geom_ribbon( data=new.data, aes( x=lnMASS, ymin = fit - fit.se, 
                                   ymax = fit + fit.se, 
                                   group=status, col=status, fill=status), 
               alpha=0.5, linetype="longdash")+
  geom_line( data=new.data, aes( x=lnMASS, y=fit, 
                                 group=status), 
             col="black", lwd=1.5 ) +
  # Reinke 2022 predictions
  geom_point( data=ectodata %>% filter(Group %in% "Reptiles" &
                                         !(Type %in% "Turtles" | Type %in% "Squamates")),
              aes(x=lnBM, y=GompertzSlope, shape = Type), size=cex, fill="darkgreen") +
  geom_line( data=reinke2022reptilesALL_MASS, aes(x=lnBM, y=fit),
             lwd=1, linetype="dashed", col="darkgreen") +
  geom_line( data=reinke2022reptilesWOT_MASS, aes(x=lnBM, y=fit),
             lwd=1, linetype="solid", col="darkgreen") +
  
  # plot settings
  scale_y_continuous(breaks = pretty( seq(min(model.data$Mean.AR.Female), 
                                          max(model.data$Mean.AR.Female, 0.5), 
                                          length.out=8) ), 
                     limits = c(-0.05, 0.5) ) +
  scale_shape_manual( values=pch, labels=c("Testudines: Captive",
                                           "Testudines: Wild",
                                           "Crocodilians","Tuataras") ) + 
  scale_color_manual( values = my.col.default) + 
  scale_fill_manual( values = my.col.lightened) +
  ggpubr::theme_pubr(base_size = base_size)+
  labs(x="ln(Mass)", y="") + theme(legend.text = element_text(size = 12)) + 
  lims(x=c( min(model.data$lnMASS), max(model.data$lnMASS, 5.5) ))+
  guides(shape=guide_legend(title = "", 
                            override.aes = list( fill = my.col.lightened)), 
         col="none", fill = "none")
pg3

(pgSFLC <- (pg1/pg2/pg3) +
    plot_layout(guides='collect') & theme(legend.position='top')) 
# ggsave(pgSFLC, width = 12, height = 12,filename = "PGLSPlotSFLCagingrateWILDCAPTIVEpoly.png")



