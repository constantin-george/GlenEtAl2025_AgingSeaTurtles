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
data.path <- "/Users/constantingeorgeglen/Documents/University/UFL/RESEARCH/Projects_Analysis/Analysis_CTC/Main/2.CTCAnalyses/PhylogeneticAnalysis/scripts/20240516/PGLS/"
load( paste0( data.path, "/submissionclean/model_data2025.RData") )


## load excel data
# model.data2 <- xlsx::read.xlsx("/Users/constantingeorgeglen/Documents/University/UFL/RESEARCH/Projects_Writing/PhDChapters/ch2.senescenceCTC/submission/Nature/supplementaryData/DataS4_testudinetraitdata.xlsx", sheetIndex = 1)
# model.data2 <- model.data2 %>% filter(!(location %in% c("CTC - OTN", "CTC - MSN")))
# model.data2$response <- log(model.data2$clutchMassKG)
# model.data2$lnmeanEX <- log(model.data2$Mean.EX.Female)
# model.data2$lnafr    <- log(model.data2$afr)
# model.data2$lnMASS   <- log(model.data2$femaleMassKG)
# model.data2$status   <- as.factor(model.data2$status)
# model.data2 <- dplyr::select(model.data2, c(response, lnmeanEX, Mean.AR.Female, lnafr,
#                                             lnMASS, status))

#### ---------------------------------------------------------------------------
#### Phylogenetic generalized least squares for reproductive mass
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
# fit.NULL  -52.5  134.2  222.8     0.0 7  0.802 
# fit.OU    -52.5  138.4  222.8     4.2 8  0.099 
# fit.lam   -52.5  138.4  222.8     4.2 8  0.099 
# fit.BK   -243.6  516.5   31.7   382.2 7  <0.001
# fit.BM   -275.3  579.8    0.0   445.6 7  <0.001


### Calculate partial and total R2s
rr2::R2_lik(fit.NULL); rr2::R2_lik(fit.OU); rr2::R2_lik(fit.lam); 
# [1] 0.81823
# [1] 0.81823
# [1] 0.81823


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


# diagnostics
hist(best.mod$residuals, breaks = 20) #check for normality of the residuals
qqnorm(best.mod$residuals) #check for normality of the residuals
qqline(best.mod$residuals) #check for n ormality of the residuals



### compare captive versus wild
emmeans(best.mod, ~ status, regrid = "response")
contrast(emmeans(best.mod, ~ status, regrid = "response"),list(c(-1,1)))
# contrast estimate    SE df t.ratio p.value
# c(-1, 1)   0.0485 0.167 60   0.290  0.7727



# Type I ANOVA table
lm.NULL <- lm(response ~ lnMASS + (lnmeanEX + Mean.AR.Female + lnafr) + status, 
              data = model.data, 
              na.action = "na.omit")
# anova_type2 <- car::Anova(lm.NULL, type = 2) # tests factors sequentially
anova_type1 <- anova(lm.NULL) # tests factors sequentially
anova_type1

# Total Sum of Squares (SST)
sst_multiple <- sum(anova_type1$"Sum Sq") # Sum of all SS including residuals
prop_var_type1 <- anova_type1$"Sum Sq" / sst_multiple
names(prop_var_type1) <- rownames(anova_type1)
print(round(prop_var_type1,2))
# lnMASS       lnmeanEX Mean.AR.Female          lnafr         status      Residuals 
#   0.81           0.01           0.00           0.00           0.00           0.18 

## proportion of variance explained by factors other than body mass
sum(prop_var_type1[c(2:5)])
# [1] 0.0055058

## proportion of variation explained by body mass
prop_var_type1[1]
# lnMASS 
# 0.81272 

## proportion of unexplained variation
prop_var_type1[6]
# Residuals 
#   0.18177 

#### ---------------------------------------------------------------------------
#### Plots
#### ---------------------------------------------------------------------------
library(viridis)
library(colorspace)
my.col <- viridis(2, begin=0.2, end=0.6, option = "magma")
pch <- c(21,22); cex <- 8; base_size <- 25
my.col.lightened <- colorspace::lighten(my.col, 0.3)
my.col.scale <- scales::alpha( 
  #c("grey20","grey30"), 
  my.col,
  1
)

# captive green turtle data
captivegreenturtle <- model.data[model.data$species=="Chelonia_mydas",] %>%
  mutate(name = "Chelonia mydas")
phylo_pic_cex <- 0.018
plotdata <- model.data %>% filter(species!="Chelonia_mydas")


# Aging rate holding everything else at their mean value
Niter    <- 1000
parMLE   <- coef(best.mod)
new.data <- by(model.data, model.data$status, function(x){
  x.data <- data.frame( Mean.AR.Female=seq(min(x$Mean.AR.Female), max(x$Mean.AR.Female), length.out=Niter),
                        lnmeanEX=mean(x$lnmeanEX, na.rm = TRUE),
                        lnMASS=mean(x$lnMASS, na.rm = TRUE),
                        lnafr=mean(x$lnafr, na.rm = TRUE))
  x.data$status <- factor(unique(x$status), levels = levels(x$status))
  return(x.data)
} )
new.data <- do.call(rbind, new.data)

fit.mat <- matrix(predict(best.mod, new.data, se=T)$fit, nrow=Niter) |> as.data.frame()
se.mat  <- matrix(predict(best.mod, new.data, se=T)$se.fit, nrow=Niter)|> as.data.frame()
new.data$fit    <- tidyr::gather(fit.mat, value = "fit")[,2]
new.data$fit.se <- tidyr::gather(se.mat, value = "fit.se")[,2]

pg1 <- ggplot(plotdata ) + 
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
  ## add a turtle for a phylopic
  rphylopic::geom_phylopic(data=captivegreenturtle, aes(x = Mean.AR.Female,
                                                        y = response,
                                                        name = name),
                           inherit.aes = FALSE, fill = my.col[[1]],
                           width = phylo_pic_cex-0.0008,
                           angle=-18) +
  scale_y_continuous(breaks = pretty( seq(min(model.data$response), 
                                          max(model.data$response), 
                                          length.out=5) ),
                     limits = c(NA, NA)) + 
  scale_shape_manual( values=pch, labels=c("Testudines: Captive","Testudines: Wild") )+ 
  scale_color_manual( values = my.col) + 
  scale_fill_manual( values = my.col.lightened) +
  ggpubr::theme_pubr(base_size = base_size)+
  labs(x="Aging rate", y="")+
  guides(shape=guide_legend(title = "", 
                            override.aes = list( fill = my.col.lightened)), 
         col="none", fill = "none") 
pg1


### EX
new.data <- by(model.data, model.data$status, function(x){
  x.data <- data.frame( lnmeanEX=seq(min(x$lnmeanEX), max(x$lnmeanEX), length.out=1000),
                        Mean.AR.Female=mean(x$Mean.AR.Female, na.rm = TRUE),
                        lnMASS=mean(x$lnMASS, na.rm = TRUE),
                        lnafr=mean(x$lnafr, na.rm = TRUE))
  x.data$status <- factor(unique(x$status), levels = levels(x$status))
  return(x.data)
} )
new.data <- do.call(rbind, new.data)

fit.mat <- matrix(predict(best.mod, new.data, se=T)$fit, nrow=100) |> as.data.frame()
se.mat  <- matrix(predict(best.mod, new.data, se=T)$se.fit, nrow=100)|> as.data.frame()
new.data$fit    <- tidyr::gather(fit.mat, value = "fit")[,2]
new.data$fit.se <- tidyr::gather(se.mat, value = "fit.se")[,2]

pg2 <- ggplot(plotdata ) + 
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
  ## add a turtle for a phylopic
  rphylopic::geom_phylopic(data=captivegreenturtle, aes(x = lnmeanEX,
                                                        y = response,
                                                        name = name),
                           inherit.aes = FALSE, fill = my.col[[1]],
                           width = phylo_pic_cex+0.18,
                           angle=-18) +
  scale_y_continuous(breaks = pretty( seq(min(model.data$response), 
                                          max(model.data$response), length.out=5) ),
                     limits = c(NA, NA)) + 
  scale_shape_manual( values=pch, labels=c("Testudines: Captive","Testudines: Wild") )+ 
  scale_color_manual( values = my.col) + 
  scale_fill_manual( values = my.col.lightened) +
  ggpubr::theme_pubr(base_size = base_size)+
  labs(x="ln(Life expectancy)", y="ln(RM)")+
  guides(shape=guide_legend(title = "", 
                            override.aes = list( fill = my.col.lightened)), 
         col="none", fill = "none") 
pg2



### for log of body mass
new.data <- by(model.data, model.data$status, function(x){
  x.data <- data.frame( lnMASS=seq(min(x$lnMASS), max(x$lnMASS), length.out=1000),
                        Mean.AR.Female=mean(x$Mean.AR.Female, na.rm = TRUE),
                        lnmeanEX=mean(x$lnmeanEX, na.rm = TRUE),
                        lnafr=mean(x$lnafr, na.rm = TRUE))
  x.data$status <- factor(unique(x$status), levels = levels(x$status))
  return(x.data)
} )
new.data <- do.call(rbind, new.data)

fit.mat <- matrix(predict(best.mod, new.data, se=T)$fit, nrow=100) |> as.data.frame()
se.mat  <- matrix(predict(best.mod, new.data, se=T)$se.fit, nrow=100)|> as.data.frame()
new.data$fit    <- tidyr::gather(fit.mat, value = "fit")[,2]
new.data$fit.se <- tidyr::gather(se.mat, value = "fit.se")[,2]

pg3 <- ggplot(plotdata ) +
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
  ## add a turtle for a phylopic
  rphylopic::geom_phylopic(data=captivegreenturtle, aes(x = lnMASS,
                                                        y = response,
                                                        name = name),
                           inherit.aes = FALSE, fill = my.col[[1]],
                           width = phylo_pic_cex+0.4,
                           angle=-18) +
  scale_shape_manual( values=pch, labels=c("Testudines: Captive","Testudines: Wild") )+
  scale_color_manual( values = my.col) +
  scale_fill_manual( values = my.col.lightened) +
  ggpubr::theme_pubr(base_size = base_size)+
  labs(x="ln(Mass)", y="")+
  guides(shape=guide_legend(title = "",
                            override.aes = list( fill = my.col.lightened)),
         col="none", fill = "none")
pg3


#### AFR
new.data <- by(model.data, model.data$status, function(x){
  x.data <- data.frame( lnafr=seq(min(x$lnafr), max(x$lnafr), length.out=1000),
                        Mean.AR.Female=mean(x$Mean.AR.Female, na.rm = TRUE),
                        lnmeanEX=mean(x$lnmeanEX, na.rm = TRUE),
                        lnMASS=mean(x$lnMASS, na.rm = TRUE))
  x.data$status <- factor(unique(x$status), levels = levels(x$status))
  return(x.data)
} )
new.data <- do.call(rbind, new.data)

fit.mat <- matrix(predict(best.mod, new.data, se=T)$fit, nrow=100) |> as.data.frame()
se.mat  <- matrix(predict(best.mod, new.data, se=T)$se.fit, nrow=100)|> as.data.frame()
new.data$fit    <- tidyr::gather(fit.mat, value = "fit")[,2]
new.data$fit.se <- tidyr::gather(se.mat, value = "fit.se")[,2]

pg4 <- ggplot(plotdata ) + 
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
  ## add a turtle for a phylopic
  rphylopic::geom_phylopic(data=captivegreenturtle, aes(x = lnafr,
                                                        y = response,
                                                        name = name),
                           inherit.aes = FALSE, fill = my.col[[1]],
                           width = phylo_pic_cex+0.11,
                           angle=-18) +
  scale_y_continuous(breaks = pretty( seq(min(model.data$response), 
                                          max(model.data$response), length.out=5) ),
                     limits = c(NA, NA)) + 
  scale_shape_manual( values=pch, labels=c("Testudines: Captive","Testudines: Wild") )+ 
  scale_color_manual( values = my.col) + 
  scale_fill_manual( values = my.col.lightened) +
  ggpubr::theme_pubr(base_size = base_size)+
  labs(x="ln(AFR)", y="")+
  guides(shape=guide_legend(title = "", 
                            override.aes = list( fill = my.col.lightened)), 
         col="none", fill = "none") 
pg4


## reproductive mass
(pgFULL <- ((pg1+labs(tag="(A)")+theme(plot.tag = element_text(size = 20)))/
              (pg2+labs(tag="(B)")+theme(plot.tag = element_text(size = 20)))/
              (pg4+labs(tag="(C)")+theme(plot.tag = element_text(size = 20)))|
              (pg3+labs(tag="(D)")+theme(plot.tag = element_text(size = 20)))) + 
    plot_layout(guides='collect') & theme(legend.position='top'))
ggsave(pgFULL, width = 18, height = 13,filename = "PGLSPlotRM_WILDCAPTIVEpoly.png")
