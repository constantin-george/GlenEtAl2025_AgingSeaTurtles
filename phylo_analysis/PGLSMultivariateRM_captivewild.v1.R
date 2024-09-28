#### ---------------------------------------------------------------------------
#### phylogenetic generalized least squares: captive and wild
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




#### ---------------------------------------------------------------------------
#### Phylogenetic generalized least squares for rep.allocation
#### ---------------------------------------------------------------------------

## define a response variable to swap out multiple
model.data$response <- model.data$lnrepMASS
response.lab        <-  "ln(reproductive mass)"


### model formula
model.forms <- list( 
  "mFULLfixed" = formula(response ~ (lnmeanEX + Mean.AR.Female + lnafr + lnMASS) + 
                           status)
)

### fit models
fit.NULL <- gls(model.forms$mFULLfixed, correlation = NULL, method="ML",
                control=glsControl(opt="nlminb"), data = model.data)
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
ICres <- bbmle::ICtab( fit.OU, fit.lam, fit.BK, fit.BM, fit.NULL, type="BIC",
                       weights = T, delta = T, base = T, logLik = T)
ICres
#          logLik BIC    dLogLik dBIC   df weight
# fit.NULL  -51.4  131.9  215.7     0.0 7  0.8   
# fit.OU    -51.4  136.1  215.7     4.2 8  0.1   
# fit.lam   -51.4  136.1  215.7     4.2 8  0.1   
# fit.BK   -235.8  500.7   31.3   368.8 7  <0.001
# fit.BM   -267.1  563.3    0.0   431.4 7  <0.001



### Calculate partial and total R2s
rr2::R2_lik(fit.NULL); rr2::R2_lik(fit.OU); rr2::R2_lik(fit.lam); 
# [1] 0.82034
# [1] 0.82034
# [1] 0.82034


#### ---------------------------------------------------------------------------
### pick the best model
best.mod <- fit.NULL


#### ----------------------------------------------------------------------------------------
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
best.mod$modelStruct


## create a table of model results
sjPlot::tab_model(best.mod, 
                  show.intercept = T,show.stat = T, show.se = T, 
                  show.dev = T, show.loglik = T, 
                  show.icc =T, show.ngroups = T, show.re.var = T, p.style = "numeric",
                  collapse.ci = T, vcov.type = "HC0", digits = 2, p.val = "wald")


### compare captive withg wild
emmeans(best.mod, ~ status, regrid = "response")
contrast(emmeans(best.mod, ~ status, regrid = "response"),list(c(-1,1)))
# contrast estimate    SE df t.ratio p.value
# c(-1, 1)   0.0219 0.175 58   0.125  0.9012


# diagnostics
hist(best.mod$residuals, breaks = 20) #check for normality of the residuals
qqnorm(best.mod$residuals) #check for normality of the residuals
qqline(best.mod$residuals) #check for normality of the residuals



#### ---------------------------------------------------------------------------
#### Plots
#### ---------------------------------------------------------------------------
library(viridis)
library(colorspace)

# my.col <- c("#377EB8","#FF7F00"); 
my.col <- viridis(2, begin=0.2, end=0.6, option = "magma")
# my.col <- c("royalblue3","forestgreen")

pch <- c(21,22); cex <- 8; base_size <- 25
my.col.lightened <- colorspace::lighten(my.col, 0.3)
my.col.scale <- scales::alpha( 
  #c("grey20","grey30"), 
  my.col,
  1
)



# Aging rate holding everything else at their mean value
Niter    <- 1000
parMLE   <- coef(best.mod)
new.data <- by(model.data, model.data$status, function(x){
  x.data <- data.frame( Mean.AR.Female=seq(min(x$Mean.AR.Female), max(x$Mean.AR.Female), length.out=Niter),
                        lnmeanEX=mean(x$lnmeanEX),
                        lnMASS=mean(x$lnMASS),
                        lnafr=mean(x$lnafr))
  x.data$status <- unique(x$status)
  return(x.data)
} )
new.data <- do.call(rbind, new.data)

fit.mat <- matrix(predict(best.mod, new.data, se=T)$fit, nrow=Niter) |> as.data.frame()
se.mat  <- matrix(predict(best.mod, new.data, se=T)$se.fit, nrow=Niter)|> as.data.frame()
new.data$fit    <- tidyr::gather(fit.mat, value = "fit")[,2]
new.data$fit.se <- tidyr::gather(se.mat, value = "fit.se")[,2]

pg1 <- ggplot(model.data ) + 
  geom_point(aes( x=Mean.AR.Female, y=response, shape=status, 
                  group=status, fill=status), 
             size=cex, col="black") + 
  geom_ribbon( data=new.data, aes( x=Mean.AR.Female, ymin = fit - fit.se, 
                                   ymax = fit + fit.se, 
                                   group=status, col=status, fill=status), 
               alpha=0.5, linetype="longdash")+
  geom_line( data=new.data, aes( x=Mean.AR.Female, y=fit, 
                                 group=status), 
             col="black", lwd=1.5 ) +
  # geom_smooth( aes( x=Mean.AR.Female, y=corr.rep.allocation, group=status), col="grey20", method="lm", se=F) +
  # geom_vline( xintercept = 0.11, lty="dashed", col="black", lwd=1.1 ) +
  scale_y_continuous(breaks = pretty( seq(min(model.data$response), 
                                          max(model.data$response), 
                                          length.out=5) ),
                     limits = c(NA, NA)) + 
  scale_shape_manual( values=pch, labels=c('Captive', "Wild") )+ 
  scale_color_manual( values = my.col) + 
  scale_fill_manual( values = my.col.lightened) +
  ggpubr::theme_pubr(base_size = base_size)+
  labs(x="Aging rate", y="")+
  guides(shape=guide_legend(title = "Type", 
                            override.aes = list( fill = my.col.lightened)), 
         col="none", fill = "none") 
pg1


### EX
new.data <- by(model.data, model.data$status, function(x){
  x.data <- data.frame( lnmeanEX=seq(min(x$lnmeanEX), max(x$lnmeanEX), length.out=1000),
                        Mean.AR.Female=mean(x$Mean.AR.Female),
                        lnMASS=mean(x$lnMASS),
                        lnafr=mean(x$lnafr))
  x.data$status <- unique(x$status)
  return(x.data)
} )
new.data <- do.call(rbind, new.data)

fit.mat <- matrix(predict(best.mod, new.data, se=T)$fit, nrow=100) |> as.data.frame()
se.mat  <- matrix(predict(best.mod, new.data, se=T)$se.fit, nrow=100)|> as.data.frame()
new.data$fit    <- tidyr::gather(fit.mat, value = "fit")[,2]
new.data$fit.se <- tidyr::gather(se.mat, value = "fit.se")[,2]

pg2 <- ggplot(model.data ) + 
  geom_point(aes( x=lnmeanEX, y=response, shape=status, 
                  group=status, fill=status), 
             size=cex, col="black") + 
  geom_ribbon( data=new.data, aes( x=lnmeanEX, ymin = fit - fit.se, 
                                   ymax = fit + fit.se, 
                                   group=status, col=status, fill=status), 
               alpha=0.5, linetype="longdash")+
  geom_line( data=new.data, aes( x=lnmeanEX, y=fit, 
                                 group=status), 
             col="black", lwd=1.5 ) +
  # geom_smooth( aes( x=Mean.AR.Female, y=corr.rep.allocation, group=status), col="grey20", method="lm", se=F) +
  # geom_vline( xintercept = 0.11, lty="dashed", col="black", lwd=1.1 ) +
  scale_y_continuous(breaks = pretty( seq(min(model.data$response), 
                                          max(model.data$response), length.out=5) ),
                     limits = c(NA, NA)) + 
  scale_shape_manual( values=pch, labels=c('Captive', "Wild") )+ 
  scale_color_manual( values = my.col) + 
  scale_fill_manual( values = my.col.lightened) +
  ggpubr::theme_pubr(base_size = base_size)+
  labs(x="ln(Life expectancy)", y=response.lab)+
  guides(shape=guide_legend(title = "Type", 
                            override.aes = list( fill = my.col.lightened)), 
         col="none", fill = "none") 
pg2



### for log of body mass
new.data <- by(model.data, model.data$status, function(x){
  x.data <- data.frame( lnMASS=seq(min(x$lnMASS), max(x$lnMASS), length.out=1000),
                        Mean.AR.Female=mean(x$Mean.AR.Female),
                        lnmeanEX=mean(x$lnmeanEX),
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
  geom_point(aes( x=lnMASS, y=response, shape=status,
                  group=status, fill=status),
             size=cex, col="black") +
  geom_ribbon( data=new.data, aes( x=lnMASS, ymin = fit - fit.se,
                                   ymax = fit + fit.se,
                                   group=status, col=status, fill=status),
               alpha=0.5, linetype="longdash")+
  geom_line( data=new.data, aes( x=lnMASS, y=fit,
                                 group=status),
             col="black", lwd=1.5 ) +
  # geom_smooth( aes( x=Mean.AR.Female, y=corr.rep.allocation, group=status), col="grey20", method="lm", se=F) +
  # geom_vline( xintercept = 0.11, lty="dashed", col="black", lwd=1.1 ) +
  scale_shape_manual( values=pch, labels=c('Captive', "Wild") )+
  scale_color_manual( values = my.col) +
  scale_fill_manual( values = my.col.lightened) +
  ggpubr::theme_pubr(base_size = base_size)+
  labs(x="ln(Mass)", y="")+
  guides(shape=guide_legend(title = "Type",
                            override.aes = list( fill = my.col.lightened)),
         col="none", fill = "none")
pg3


#### AFR
new.data <- by(model.data, model.data$status, function(x){
  x.data <- data.frame( lnafr=seq(min(x$lnafr), max(x$lnafr), length.out=1000),
                        Mean.AR.Female=mean(x$Mean.AR.Female),
                        lnmeanEX=mean(x$lnmeanEX),
                        lnMASS=mean(x$lnMASS))
  x.data$status <- unique(x$status)
  return(x.data)
} )
new.data <- do.call(rbind, new.data)

fit.mat <- matrix(predict(best.mod, new.data, se=T)$fit, nrow=100) |> as.data.frame()
se.mat  <- matrix(predict(best.mod, new.data, se=T)$se.fit, nrow=100)|> as.data.frame()
new.data$fit    <- tidyr::gather(fit.mat, value = "fit")[,2]
new.data$fit.se <- tidyr::gather(se.mat, value = "fit.se")[,2]

pg4 <- ggplot(model.data ) + 
  geom_point(aes( x=lnafr, y=response, shape=status, 
                  group=status, fill=status), 
             size=cex, col="black") + 
  geom_ribbon( data=new.data, aes( x=lnafr, ymin = fit - fit.se, 
                                   ymax = fit + fit.se, 
                                   group=status, col=status, fill=status), 
               alpha=0.5, linetype="longdash")+
  geom_line( data=new.data, aes( x=lnafr, y=fit, 
                                 group=status), 
             col="black", lwd=1.5 ) +
  # geom_smooth( aes( x=lnafr, y=corr.rep.allocation, group=status), col="grey20", method="lm", se=F) +
  # geom_vline( xintercept = 0.11, lty="dashed", col="black", lwd=1.1 ) +
  scale_y_continuous(breaks = pretty( seq(min(model.data$response), 
                                          max(model.data$response), length.out=5) ),
                     limits = c(NA, NA)) + 
  scale_shape_manual( values=pch, labels=c('Captive', "Wild") )+ 
  scale_color_manual( values = my.col) + 
  scale_fill_manual( values = my.col.lightened) +
  ggpubr::theme_pubr(base_size = base_size)+
  labs(x="ln(AFR)", y="")+
  guides(shape=guide_legend(title = "Type", 
                            override.aes = list( fill = my.col.lightened)), 
         col="none", fill = "none") 
pg4

## reproductive mass
(pgFULL <- ((pg1/pg2/pg4)|pg3) + plot_layout(guides='collect') & theme(legend.position='top'))
# ggsave(pgFULL, width = 18, height = 13,filename = "PGLSPlotRM_WILDCAPTIVEpoly.png")




library(RColorBrewer); library(ggsci)
# cols <- viridis(10, option = "plasma")
cols <- pal_jco()(10)
pd <- position_dodge(0.3)
(gp1 <- ggplot(model.data, aes(group=species)) + 
    geom_point(aes( status, lnrepMASS, fill=family, shape=status ), 
               col="black", size=5, position=pd) + 
    # geom_line(aes(status, lnrepMASS, col=family), 
    #           position=pd) +
    geom_boxplot( aes(status, lnrepMASS), alpha=0, inherit.aes=F) + 
    scale_shape_manual( values=c(21,22), labels=c("Captive", "Wild") )+
    scale_x_discrete(labels=c("Captive", "Wild")) +
    scale_fill_manual(values = cols) + scale_color_manual(values = cols) +
    guides( fill=guide_legend(title = "Family", override.aes = list(shape=21)),
            shape="none", col="none") + 
    ggpubr::geom_signif(aes(status, lnrepMASS), 
                        test = "t.test",
                        comparisons = list(c("captive", "wild")), 
                        map_signif_level = F, y_position = 4)+
    labs(x="", y="ln(Reproductive Mass)") + ggpubr::theme_pubr(base_size = 25))

(gp2 <- ggplot(model.data, aes(group=species)) + 
    geom_point(aes( status, lnmeanEX, fill=family, shape=status ), 
               col="black", size=5, position=pd) + 
    # geom_line(aes(status, lnmeanEX, col=family), 
    #           position=pd) +
    geom_boxplot( aes(status, lnmeanEX), alpha=0, inherit.aes=F) + 
    scale_shape_manual( values=c(21,22), labels=c("Captive", "Wild") )+
    scale_x_discrete(labels=c("Captive", "Wild")) +
    scale_fill_manual(values = cols) + scale_color_manual(values = cols) +
    guides( fill=guide_legend(title = "Family", override.aes = list(shape=21)),
            shape="none", col="none") + 
    ggpubr::geom_signif(aes(status, lnmeanEX), test = "t.test",
                        comparisons = list(c("captive", "wild")), 
                        map_signif_level = F, y_position = 4.1)+
    labs(x="", y="ln(Life expectancy)") + ggpubr::theme_pubr(base_size = 25))

(gp3 <- ggplot(model.data, aes(group=species)) + 
    geom_point(aes( status, Mean.AR.Female, fill=family, shape=status ), 
               col="black", size=5, position=pd) + 
    # geom_line(aes(status, Mean.AR.Female, col=family), 
    #           position=pd) +
    geom_boxplot( aes(status, Mean.AR.Female), alpha=0, inherit.aes=F) + 
    scale_shape_manual( values=c(21,22), labels=c("Captive", "Wild") )+
    scale_x_discrete(labels=c("Captive", "Wild")) +
    scale_fill_manual(values = cols) + scale_color_manual(values = cols) +
    guides( fill=guide_legend(title = "Family", override.aes = list(shape=21)),
            shape="none", col="none") + 
    ggpubr::geom_signif(aes(status, Mean.AR.Female), test = "t.test",
                        comparisons = list(c("captive", "wild")), 
                        map_signif_level = F, y_position = 0.25)+
    labs(x="", y="Aging rate") + ggpubr::theme_pubr(base_size = 25))

(gp4 <- ggplot(model.data, aes(group=species)) + 
    geom_point(aes( status, lnMASS, fill=family, shape=status ), 
               col="black", size=5, position=pd) + 
    # geom_line(aes(status, lnMASS, col=family), 
    #           position=pd) +
    geom_boxplot( aes(status, lnMASS), alpha=0, inherit.aes=F) + 
    scale_shape_manual( values=c(21,22), labels=c("Captive", "Wild") )+
    scale_x_discrete(labels=c("Captive", "Wild")) +
    scale_fill_manual(values = cols) + scale_color_manual(values = cols) +
    guides( fill=guide_legend(title = "Family", override.aes = list(shape=21)),
            shape="none", col="none") + 
    ggpubr::geom_signif(aes(status, lnMASS), test = "t.test",
                        comparisons = list(c("captive", "wild")), 
                        map_signif_level = F, y_position = 6.2)+
    labs(x="", y="ln(Mass)") + ggpubr::theme_pubr(base_size = 25))

(gp5 <- ggplot(model.data, aes(group=species)) + 
    geom_point(aes( status, lnafr, fill=family, shape=status ), 
               col="black", size=5, position=pd) + 
    # geom_line(aes(status, lnafr, col=family), 
    #           position=pd) +
    geom_boxplot( aes(status, lnafr), alpha=0, inherit.aes=F) + 
    scale_shape_manual( values=c(21,22), labels=c("Captive", "Wild") )+
    scale_x_discrete(labels=c("Captive", "Wild")) +
    scale_fill_manual(values = cols) + scale_color_manual(values = cols) +
    guides( fill=guide_legend(title = "Family", override.aes = list(shape=21)),
            shape="none", col="none") + 
    ggpubr::geom_signif(aes(status, lnafr), test = "t.test",
                        comparisons = list(c("captive", "wild")), 
                        map_signif_level = F, y_position = 3.5)+
    labs(x="", y="ln(AFR)") + ggpubr::theme_pubr(base_size = 25))

(gp <- (gp1|gp3|gp5|gp2|gp4) + plot_layout(guides='collect') & theme(legend.position='top'))
# ggsave(gp, width = 20, height = 8,filename = "PGLSPlotALLcovariates.png")
 

















