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
data.path <- "/Users/constantingeorgeglen/Documents/University/UFL/RESEARCH/Projects_Analysis/Analysis_CTC/Main/2.CTCAnalyses/PhylogeneticAnalysis/scripts/20240516/PGLS/"
load( paste0( data.path, "/submissionclean/model_data2025.RData") )
model.data$lnAnnualRepOut <- log(model.data$AnnualRepOut)

## load excel data
# model.data2 <- xlsx::read.xlsx("/Users/constantingeorgeglen/Documents/University/UFL/RESEARCH/Projects_Writing/PhDChapters/ch2.senescenceCTC/submission/Nature/supplementaryData/DataS4_testudinetraitdata.xlsx", sheetIndex = 1)
# model.data2 <- model.data2 %>% filter(!(location %in% c("CTC - OTN", "CTC - MSN")))
# model.data2$lnAnnualRepOut <- log(model.data2$AnnualRepOut)
# model.data2$lnafr    <- log(model.data2$afr)
# model.data2$lnMASS   <- log(model.data2$femaleMassKG)
# model.data2$status   <- as.factor(model.data2$status)
# model.data2 <- dplyr::select(model.data2, c(lnAnnualRepOut, Mean.AR.Female, lnafr,
#                                             lnMASS, status))


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
# fit.NULL  120.5 -215.8  331.9     0.0 6  0.802 
# fit.OU    120.5 -211.6  331.9     4.2 7  0.099 
# fit.lam   120.5 -211.6  331.9     4.2 7  0.099 
# fit.BK   -179.6  384.3   31.8   600.1 6  <0.001
# fit.BM   -211.4  448.0    0.0   663.8 6  <0.001


### Calculate partial and total R2s
rr2::R2_lik(fit.OU); rr2::R2_lik(fit.lam); 
# [1] 0.080175
# [1] 0.080175


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


### compare captive withg wild
emmeans(best.mod, ~ status, regrid = "response")
contrast(emmeans(best.mod, ~ status, regrid = "response"),list(c(-1,1)))
# contrast estimate    SE df t.ratio p.value
# c(-1, 1)   0.0151 0.0123 61   1.231  0.2232



# Type I ANOVA table
lm.NULL <- lm(model.forms$mFULLfixed,data = model.data, 
              na.action = "na.omit")
anova_type1 <- anova(lm.NULL) # tests factors sequentially
anova_type1

# Total Sum of Squares (SST)
sst_multiple <- sum(anova_type1$"Sum Sq") # Sum of all SS including residuals
prop_var_type1 <- anova_type1$"Sum Sq" / sst_multiple
names(prop_var_type1) <- rownames(anova_type1)
print(round(prop_var_type1,2))
sum(prop_var_type1[c(1:4)])
# [1] 0.080175


#### ---------------------------------------------------------------------------
#### Plots
#### ---------------------------------------------------------------------------
library(viridis)
library(colorspace)
library(rphylopic)

# load Reinke 2022 predictions
load("reinke2022ReptilesALLPredDAT.RData")
load("reinke2022ReptilesWOTPredDAT.RData")


# my.col <- c("#377EB8","#FF7F00"); 
my.col.default <- viridis(2, begin=0.2, end=0.6, option = "magma")
my.col.class <- c( my.col.default, 
                   rep( viridis(1, begin=0.6, end=0.8, option = "mako"), 2 ) )

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

## get turtle phylopic
turtle <- rphylopic::get_uuid("Chelonia mydas")

# captive green turtle data
captivegreenturtle <- model.data[model.data$species=="Chelonia_mydas",] %>%
  mutate(name = "Chelonia mydas")
captivegreenturtle$image <- turtle
phylo_pic_cex <- 0.36


# main plot 1
plotdata <- model.data %>% filter(species!="Chelonia_mydas")
pg1 <- ggplot(plotdata ) + 
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
  # geom_point( data=ectodata %>% filter(Group %in% "Reptiles" &
  #                                        !(Type %in% "Turtles" | Type %in% "Squamates")),
  #             aes(x=lnAF, y=GompertzSlope, shape = Type), size=cex, fill=my.col.class[3]) +
  geom_line( data=reinke2022reptilesALL_AF, aes(x=lnAF, y=fn(fit, logT=logT)),
             lwd=1, linetype="dashed", col=my.col.class[3]) +
  geom_line( data=reinke2022reptilesWOT_AF, aes(x=lnAF, y=fn(fit, logT=logT)),
             lwd=1, linetype="solid", col=my.col.class[3]) +
  
  ## add a turtle for a phylopic
  rphylopic::geom_phylopic(data=captivegreenturtle, aes(x = lnAnnualRepOut,
                                                        y = fn(Mean.AR.Female, logT=logT),
                                                        name = name),
                           inherit.aes = FALSE, fill = my.col.default[[1]],
                           width = phylo_pic_cex+0.08,
                           angle=-18) +
  
  # plot settings
  scale_y_continuous(breaks = fn(pretty( seq(min(model.data$Mean.AR.Female), 
                                             max(model.data$Mean.AR.Female, 
                                                 reinke2022reptilesWOT_AF$fit, 0.6), 
                                             length.out=8)), logT=logT )) + 
  scale_shape_manual( values=pch, labels=c("Testudines: Captive", 
                                           "Testudines: Wild") ) + 
  scale_color_manual( values = my.col.default) + 
  scale_fill_manual( values = my.col.lightened) +
  ggpubr::theme_pubr(base_size = base_size) + 
  theme(legend.text = element_text(size = 12)) + 
  labs(x="ln(Annual fecundity)", y="Aging rate") +
  guides(shape=guide_legend(title = "", 
                            override.aes = list( fill = my.col.lightened[1:2])), 
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

pg2 <- ggplot(plotdata ) + 
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
  # geom_point( data=ectodata %>% filter(Group %in% "Reptiles" & 
  #                                        !(Type %in% "Turtles" | Type %in% "Squamates")),
  #             aes(x=lnAFR, y=GompertzSlope, shape = Type), size=cex, fill=my.col.class[3]) +
  geom_line( data=reinke2022reptilesALL_AFR, aes(x=lnAFR, y=fit),
             lwd=1, linetype="dashed", col=my.col.class[3]) +
  geom_line( data=reinke2022reptilesWOT_AFR, aes(x=lnAFR, y=fit),
             lwd=1, linetype="solid", col=my.col.class[3]) +
  
  ## add a turtle for a phylopic
  rphylopic::geom_phylopic(data=captivegreenturtle, aes(x = lnafr, y = Mean.AR.Female,
                                                        name = name), 
                           inherit.aes = FALSE, fill = my.col.default[[1]], 
                           width = phylo_pic_cex/2.15, 
                           angle=-18) +
  
  # plot settings
  scale_y_continuous(breaks = pretty( seq(min(model.data$Mean.AR.Female), 
                                          max(model.data$Mean.AR.Female, 0.6), 
                                          length.out=8) ), 
                     limits = c(-0.05, 0.55) ) + 
  scale_shape_manual( values=pch, labels=c("Testudines: Captive",
                                           "Testudines: Wild") ) + 
  scale_color_manual( values = my.col.default) + 
  scale_fill_manual( values = my.col.lightened) +
  ggpubr::theme_pubr(base_size = base_size)+
  labs(x="ln(AFR)", y="Aging rate") + theme(legend.text = element_text(size = 12)) + 
  lims(x=range(model.data$lnafr))+
  guides(shape=guide_legend(title = "", 
                            override.aes = list( fill = my.col.lightened[1:2])), 
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

pg3 <- ggplot( plotdata ) + 
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
  # geom_point( data=ectodata %>% filter(Group %in% "Reptiles" &
  #                                        !(Type %in% "Turtles" | Type %in% "Squamates")),
  #             aes(x=lnBM, y=GompertzSlope, shape = Type), size=cex, fill=my.col.class[3]) +
  geom_line( data=reinke2022reptilesALL_MASS, aes(x=lnBM, y=fit),
             lwd=1, linetype="dashed", col=my.col.class[3]) +
  geom_line( data=reinke2022reptilesWOT_MASS, aes(x=lnBM, y=fit),
             lwd=1, linetype="solid", col=my.col.class[3]) +
  
  ## add a turtle for a phylopic
  rphylopic::geom_phylopic(data=captivegreenturtle, aes(x = lnMASS, y = Mean.AR.Female,
                                                        name = name), 
                           inherit.aes = FALSE, fill = my.col.default[[1]], 
                           width = phylo_pic_cex+0.22, 
                           angle=-18) +
  
  # plot settings
  scale_y_continuous(breaks = pretty( seq(min(model.data$Mean.AR.Female), 
                                          max(model.data$Mean.AR.Female, 0.5), 
                                          length.out=8) ), 
                     limits = c(-0.05, 0.55) ) +
  scale_shape_manual( values=pch, labels=c("Testudines: Captive","Testudines: Wild"
                                           #"Crocodilians","Tuataras"
  ) ) + 
  scale_color_manual( values = my.col.default) + 
  scale_fill_manual( values = my.col.lightened) +
  ggpubr::theme_pubr(base_size = base_size)+
  labs(x="ln(Mass)", y="Aging rate") + theme(legend.text = element_text(size = 12)) + 
  lims(x=c( min(model.data$lnMASS), max(model.data$lnMASS, 5.5) ))+
  guides(shape=guide_legend(title = "", 
                            override.aes = list( fill = my.col.lightened[1:2])), 
         col="none", fill = "none")
pg3

(pgSFLC <- (pg1+labs(tag="(A)")+theme(plot.tag = element_text(size = 20))|
              pg2+labs(y="",tag="(B)")+theme(plot.tag = element_text(size = 20))|
              pg3+labs(y="",tag="(C)")+theme(plot.tag = element_text(size = 20))) + 
    plot_layout(guides='collect') & theme(legend.position='top', legend.text = element_text(size = 20))) 
ggsave(pgSFLC, width = 20, height = 10, filename = "PGLSPlotSFLCagingrateWILDCAPTIVEpoly.png")



