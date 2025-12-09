#--------------------------------------------------------------
#Ben Neely
#11/18/2025
#Model factors associated with consumption of key prey items
#--------------------------------------------------------------

## Clear R
cat("\014")  
rm(list=ls())

## Install and load packages
## Checks if package is installed, installs if not, activates for current session
if("tidyverse" %in% rownames(installed.packages()) == FALSE) {install.packages("tidyverse")}
library(tidyverse)

if("lme4" %in% rownames(installed.packages()) == FALSE) {install.packages("lme4")}
library(lme4)

if("emmeans" %in% rownames(installed.packages()) == FALSE) {install.packages("emmeans")}
library(emmeans)

if("performance" %in% rownames(installed.packages()) == FALSE) {install.packages("performance")}
library(performance)

if("broom.mixed" %in% rownames(installed.packages()) == FALSE) {install.packages("broom.mixed")}
library(broom.mixed)

if("patchwork" %in% rownames(installed.packages()) == FALSE) {install.packages("patchwork")}
library(patchwork)

## Set ggplot theme
pubtheme=theme_classic()+
  theme(panel.grid=element_blank(), 
        panel.background=element_blank(),
        plot.background=element_blank(),
        panel.border=element_rect(fill="transparent"),
        axis.title=element_text(size=22,color="black",face="bold"),
        axis.text=element_text(size=18,color="black"),
        legend.position="inside",
        legend.position.inside=c(0.01,0.99),
        legend.justification=c("left","top"),
        legend.title=element_text(size=16),
        legend.text=element_text(size=14),
        legend.background=element_rect(fill=NA),
        strip.text=element_text(size=16,face="bold"),
        panel.spacing=unit(0.75,"cm"))
options(scipen=999)

## Read in raw data
dat=read_csv("lmbdiets.csv")

## Filter data to retain bass and pertinent columns
## We'll also get rid of seven rows with no fish length
## There are some instances of the same id having different tl so fix that
## The mollusc line just fixes a typo in the database
lmb=dat%>%
  filter(Species=="MICSAL",
         `Fish TL (mm)`>0)%>%
  mutate(sn1=paste(StomachName,`Fish TL (mm)`,sep="_"))%>%
  select(id=sn1,
         tl=`Fish TL (mm)`,
         fineitem=`Subcategory`,
         coarseitem=`Overarching Category`)%>%
  mutate(coarseitem=case_when(coarseitem=="Mollusc" ~ "Mollusca",TRUE ~ coarseitem),
         size_class=case_when(tl<150 ~ "small",
                              tl>=150 & tl<250 ~ "medium",
                              tl>=250 ~ "large"))

## Extract impd, month, year, and species from id
## Regex change:
# We added '_.*$' to the end.
# This tells R: "Find the 6-character species code that is followed by an underscore, 
# and then ignore everything else (the length) until the end of the string."
lmb1=lmb%>%
  extract(id,into=c("impd","month","year","spp"),
          regex="^(.{4})_([A-Z]+)(\\d{4}).*_(.{6})_.*$", 
          remove=FALSE)

## Clean up a typo and make sure that a fish doesn't have empty and something else
lmb1=lmb1%>%
  mutate(impd=case_when(impd=="OSTL" ~ "OTSL",
                        TRUE ~ impd))%>%
  group_by(id)%>%
  filter(fineitem!="Empty" | all(fineitem=="Empty")) %>%
  ungroup()%>%
  distinct()

## Now we need to isolate what we'll need for modeling
## We're going to use Fish, Shad, Centrarchids, macroinverts, crayfish, and zoops
## Lepomis, Micropterus, and Dorosoma as response vars
## We're going to save tl, size class, month, and shad_cpe for fixed effects
## We're going to save impd as a random effect
lmb2=lmb1%>%
  mutate(fish=case_when(coarseitem=="Fish" ~ 1, TRUE ~ 0),
         dorosoma=case_when(fineitem=="Dorosoma" ~ 1, TRUE ~ 0),
         centrarchid=case_when(fineitem=="Centrarchidae" |
                                 fineitem=="Lepomis" |
                                 fineitem=="Pomoxis" |
                                 fineitem=="Micropterus" ~ 1,
                               TRUE ~ 0),
         aquaticinsect=case_when(coarseitem=="Macroinvertebrate" ~ 1, TRUE ~ 0),
         crayfish=case_when(coarseitem=="Crayfish" ~ 1, TRUE ~ 0),
         zooplankton=case_when(coarseitem=="Zooplankton" |
                                 coarseitem=="Cladocera" ~ 1,
                               TRUE ~ 0))

## Reorganize data so we have one fish per line
lmb3=lmb2%>%
  group_by(id,impd,month,year,tl,size_class)%>%
  summarize(across(fish:zooplankton, ~max(.,na.rm=T)),
            .groups="drop")

## Let's clean up our data a little bit to make modeling easier
lmb4=lmb3%>%
  mutate(impd=factor(impd,
                     levels=c("BUSL","CASL","CLSL","DGSL",
                               "GESL","MESL","NOSL","OTSL",
                               "PTSA","SCSL","SNSL","WSSL"),
                     labels=c("Butler","Clark","Cowley","Douglas",
                              "Geary","Meade","Neosho","Ottawa",
                              "Pottawatomie","Scott","Shawnee","Washington")),
         spp=factor("Largemouth Bass"),
         size_class=factor(size_class,
                           levels=c("small","medium","large")),
         month=factor(month,
                      levels=c("JUNE","AUG","OCT"),
                      labels=c("June","August","October")),
         year=as.numeric(year))

## Finally, we need to add in shad abundance (KDWP gill net CPE)
## First we'll format it like the diet data so we can join
shad=read_csv("shad abundance kdwp/gzs_cpe.csv")%>%
  mutate(impd=factor(impd,
                     levels=c("BUSL","CASL","CLSL","DGSL",
                              "GESL","MESL","NOSL","OTSL",
                              "PTSA","SCSL","SNSL","WSSL"),
                     labels=c("Butler","Clark","Cowley","Douglas",
                              "Geary","Meade","Neosho","Ottawa",
                              "Pottawatomie","Scott","Shawnee","Washington")),
         year=as.numeric(year))

## Combine to create model ready data
lmb5=lmb4%>%
  left_join(shad,
            by=c("impd","year"))%>%
  relocate(gzs_cpe,.after=size_class)

## Check for NAs in predictor variables
lmb5[!complete.cases(lmb5), ]

## Now we need to scale tl and gzs_cpe to promote model convergence
moddat=lmb5%>%
  mutate(tl_z=as.numeric(scale(tl)),
         gzs_cpe_z=as.numeric(scale(gzs_cpe)))

## Finally, we'll apply effects coding to months so comparisons are to the grand mean
contrasts(moddat$month)=contr.sum

################################################################################
## Next step is to create a prediction grid so we can visualize model results
## x-axis is gzs_cpe, y is probability of piscivory, color is size_group, facet is month
## We need to do some transformation magic to get a prediction grid together

## We use the *original* gzs_cpe column from moddat
gzs_mean=mean(moddat$gzs_cpe,na.rm=T)
gzs_sd=sd(moddat$gzs_cpe,na.rm=T)

## Now we need scaling parameters for each size class
size_class_means=moddat%>%
  group_by(size_class)%>%
  summarize(tl_z=mean(tl_z,na.rm=T),
            .groups="drop")

## Create a smooth sequence of REAL, UN-SCALED gzs_cpe values
gzs_grid_real=seq(min(moddat$gzs_cpe,na.rm=T),
                  max(moddat$gzs_cpe,na.rm=T),
                  length.out=100)

## Create the base prediction grid
pred_grid=crossing(month=levels(moddat$month),
                   gzs_cpe=gzs_grid_real,
                   size_class=levels(moddat$size_class))%>%
  left_join(size_class_means,by="size_class")

## Create a new data frame for the model to use
## We add the scaled columns the model
model_input_grid=pred_grid%>%
  mutate(gzs_cpe_z=(gzs_cpe-gzs_mean)/gzs_sd)
################################################################################

################################################################################
## Build fish model
## Occurrence of fish in diet based on month, scaled tl, scaled gzs cpe, and impd
fishmod=glmer(fish~month+tl_z+gzs_cpe_z+(1|impd),
              data=moddat,
              family=binomial(link="logit"))

## Generate predictions using the scaled grid built above
preds_fishmod=predict(fishmod,
                      newdata=model_input_grid,
                      re.form=NA,
                      se.fit=T)

## Add the predictions back to the unscaled grid with confidence intervals
fishmod_pdat=pred_grid%>%
  mutate(size_class=factor(size_class,levels=c("small","medium","large")))%>%
  mutate(log_odds=preds_fishmod$fit,
         log_odds_se=preds_fishmod$se.fit)%>%
  mutate(log_odds_lci=log_odds-1.96*log_odds_se,
         log_odds_uci=log_odds+1.96*log_odds_se)%>%
  mutate(probs=plogis(log_odds),
         probs_lci=plogis(log_odds_lci),
         probs_uci=plogis(log_odds_uci))

## Finally, create the plot
fish_plot=ggplot(fishmod_pdat)+
  geom_line(aes(x=gzs_cpe,y=probs,color=size_class),linewidth=1.2)+
  scale_color_manual(values=c("#FDD0A2","#FD8D3C","#D94701"),
                     labels=c("< 150 mm",
                              "150 to 249 mm",
                              "≥ 250 mm"),
                     name="Size group")+
  scale_x_continuous(breaks=seq(0,30,5),
                     name="")+
  scale_y_continuous(breaks=seq(0,0.50,0.10),
                     name=expression(italic(P)(Fish)))+
  coord_cartesian(xlim=c(-0.3,30.3),
                  ylim=c(0,0.505),
                  expand=F)+
  pubtheme+
  theme(legend.position="none")+
  facet_wrap(~factor(month,
                     levels=c("June","August","October")))

################################################################################
## Build shad model
## Occurrence of gzs in diet based on month, scaled tl, scaled gzs cpe, and impd
gzsmod=glmer(dorosoma~month+tl_z+gzs_cpe_z+(1|impd),
              data=moddat,
              family=binomial(link="logit"))

## Generate predictions using the scaled grid built above
preds_gzsmod=predict(gzsmod,
                      newdata=model_input_grid,
                      re.form=NA,
                      se.fit=T)

## Add the predictions back to the unscaled grid with confidence intervals
gzsmod_pdat=pred_grid%>%
  mutate(size_class=factor(size_class,levels=c("small","medium","large")))%>%
  mutate(log_odds=preds_gzsmod$fit,
         log_odds_se=preds_gzsmod$se.fit)%>%
  mutate(log_odds_lci=log_odds-1.96*log_odds_se,
         log_odds_uci=log_odds+1.96*log_odds_se)%>%
  mutate(probs=plogis(log_odds),
         probs_lci=plogis(log_odds_lci),
         probs_uci=plogis(log_odds_uci))

## Finally, create the plot
gzs_plot=ggplot(gzsmod_pdat)+
  geom_line(aes(x=gzs_cpe,y=probs,color=size_class),linewidth=1.2)+
  scale_color_manual(values=c("#FDD0A2","#FD8D3C","#D94701"),
                     labels=c("< 150 mm",
                              "150 to 249 mm",
                              "≥ 250 mm"),
                     name="Size group")+
  scale_x_continuous(breaks=seq(0,30,5),
                     name="")+
  scale_y_continuous(breaks=seq(0,0.10,0.02),
                     name=expression(italic(P)(`Gizzard Shad`)))+
  coord_cartesian(xlim=c(-0.3,30.3),
                  ylim=c(0,0.101),
                  expand=F)+
  pubtheme+
  facet_wrap(~factor(month,
                     levels=c("June","August","October")))

################################################################################
## Build Centrarchidae model
## Occurrence of centrarchid in diet based on month, scaled tl, scaled gzs cpe, and impd
centrarchidmod=glmer(centrarchid~month+tl_z+gzs_cpe_z+(1|impd),
                  data=moddat,
                  family=binomial(link="logit"))

## Generate predictions using the scaled grid built above
preds_centrarchidmod=predict(centrarchidmod,
                          newdata=model_input_grid,
                          re.form=NA,
                          se.fit=T)

## Add the predictions back to the unscaled grid with confidence intervals
centrarchidmod_pdat=pred_grid%>%
  mutate(size_class=factor(size_class,levels=c("small","medium","large")))%>%
  mutate(log_odds=preds_centrarchidmod$fit,
         log_odds_se=preds_centrarchidmod$se.fit)%>%
  mutate(log_odds_lci=log_odds-1.96*log_odds_se,
         log_odds_uci=log_odds+1.96*log_odds_se)%>%
  mutate(probs=plogis(log_odds),
         probs_lci=plogis(log_odds_lci),
         probs_uci=plogis(log_odds_uci))

## Finally, create the plot
centrarchid_plot=ggplot(centrarchidmod_pdat)+
  geom_line(aes(x=gzs_cpe,y=probs,color=size_class),linewidth=1.2)+
  scale_color_manual(values=c("#FDD0A2","#FD8D3C","#D94701"),
                     labels=c("< 150 mm",
                              "150 to 249 mm",
                              "≥ 250 mm"),
                     name="Size group")+
  scale_x_continuous(breaks=seq(0,30,5),
                     name=expression("Gizzard Shad CPE"))+
  scale_y_continuous(breaks=seq(0,0.25,0.05),
                     name=expression(italic(P)(Centrarchidae)))+
  coord_cartesian(xlim=c(-0.3,30.3),
                  ylim=c(0,0.2525),
                  expand=F)+
  pubtheme+
  theme(legend.position="none")+
  facet_wrap(~factor(month,
                     levels=c("June","August","October")))

################################################################################
## Build Crayfish model
## Occurrence of crayfish in diet based on month, scaled tl, scaled gzs cpe, and impd
crayfishmod=glmer(crayfish~month+tl_z+gzs_cpe_z+(1|impd),
             data=moddat,
             family=binomial(link="logit"))

## Generate predictions using the scaled grid built above
preds_crayfishmod=predict(crayfishmod,
                     newdata=model_input_grid,
                     re.form=NA,
                     se.fit=T)

## Add the predictions back to the unscaled grid with confidence intervals
crayfishmod_pdat=pred_grid%>%
  mutate(size_class=factor(size_class,levels=c("small","medium","large")))%>%
  mutate(log_odds=preds_crayfishmod$fit,
         log_odds_se=preds_crayfishmod$se.fit)%>%
  mutate(log_odds_lci=log_odds-1.96*log_odds_se,
         log_odds_uci=log_odds+1.96*log_odds_se)%>%
  mutate(probs=plogis(log_odds),
         probs_lci=plogis(log_odds_lci),
         probs_uci=plogis(log_odds_uci))

## Finally, create the plot
crayfish_plot=ggplot(crayfishmod_pdat)+
  geom_line(aes(x=gzs_cpe,y=probs,color=size_class),linewidth=1.2)+
  scale_color_manual(values=c("#FDD0A2","#FD8D3C","#D94701"),
                     labels=c("< 150 mm",
                              "150 to 249 mm",
                              "≥ 250 mm"),
                     name="Size group")+
  scale_x_continuous(breaks=seq(0,30,5),
                     name="")+
  scale_y_continuous(breaks=seq(0,0.50,0.10),
                     name=expression(italic(P)(Crayfish)))+
  coord_cartesian(xlim=c(-0.3,30.3),
                  ylim=c(0,0.505),
                  expand=F)+
  pubtheme+
  theme(legend.position="none")+
  facet_wrap(~factor(month,
                     levels=c("June","August","October")))

################################################################################
## Build aquatic insect model
## Occurrence of aquatic insects in diet based on month, scaled tl, scaled gzs cpe, and impd
aquaticinsectmod=glmer(aquaticinsect~month+tl_z+gzs_cpe_z+(1|impd),
                  data=moddat,
                  family=binomial(link="logit"))

## Generate predictions using the scaled grid built above
preds_aquaticinsectmod=predict(aquaticinsectmod,
                               newdata=model_input_grid,
                               re.form=NA,
                               se.fit=T)

## Add the predictions back to the unscaled grid with confidence intervals
aquaticinsectmod_pdat=pred_grid%>%
  mutate(size_class=factor(size_class,levels=c("small","medium","large")))%>%
  mutate(log_odds=preds_aquaticinsectmod$fit,
         log_odds_se=preds_aquaticinsectmod$se.fit)%>%
  mutate(log_odds_lci=log_odds-1.96*log_odds_se,
         log_odds_uci=log_odds+1.96*log_odds_se)%>%
  mutate(probs=plogis(log_odds),
         probs_lci=plogis(log_odds_lci),
         probs_uci=plogis(log_odds_uci))

## Finally, create the plot
aquaticinsect_plot=ggplot(aquaticinsectmod_pdat)+
  geom_line(aes(x=gzs_cpe,y=probs,color=size_class),linewidth=1.2)+
  scale_color_manual(values=c("#FDD0A2","#FD8D3C","#D94701"),
                     labels=c("< 150 mm",
                              "150 to 249 mm",
                              "≥ 250 mm"),
                     name="Size group")+
  scale_x_continuous(breaks=seq(0,30,5),
                     name="")+
  scale_y_continuous(breaks=seq(0,0.60,0.10),
                     name=expression(italic(P)(`Aquatic insect`)))+
  coord_cartesian(xlim=c(-0.3,30.3),
                  ylim=c(0,0.616),
                  expand=F)+
  pubtheme+
  theme(legend.position="none")+
  facet_wrap(~factor(month,
                     levels=c("June","August","October")))

################################################################################
## Build Zooplankton model
## Occurrence of zooplankton in diet based on month, scaled tl, scaled gzs cpe, and impd
zooplanktonmod=glmer(zooplankton~month+tl_z+gzs_cpe_z+(1|impd),
                     data=moddat,
                     family=binomial(link="logit"))

## Generate predictions using the scaled grid built above
preds_zooplanktonmod=predict(zooplanktonmod,
                             newdata=model_input_grid,
                             re.form=NA,
                             se.fit=T)

## Add the predictions back to the unscaled grid with confidence intervals
zooplanktonmod_pdat=pred_grid%>%
  mutate(size_class=factor(size_class,levels=c("small","medium","large")))%>%
  mutate(log_odds=preds_zooplanktonmod$fit,
         log_odds_se=preds_zooplanktonmod$se.fit)%>%
  mutate(log_odds_lci=log_odds-1.96*log_odds_se,
         log_odds_uci=log_odds+1.96*log_odds_se)%>%
  mutate(probs=plogis(log_odds),
         probs_lci=plogis(log_odds_lci),
         probs_uci=plogis(log_odds_uci))

## Finally, create the plot
zoop_plot=ggplot(zooplanktonmod_pdat)+
  geom_line(aes(x=gzs_cpe,y=probs,color=size_class),linewidth=1.2)+
  scale_color_manual(values=c("#FDD0A2","#FD8D3C","#D94701"),
                     labels=c("< 150 mm",
                              "150 to 249 mm",
                              "≥ 250 mm"),
                     name="Size group")+
  scale_x_continuous(breaks=seq(0,30,5),
                     name=expression("Gizzard Shad CPE"))+
  scale_y_continuous(breaks=seq(0,0.25,0.05),
                     name=expression(italic(P)(Zooplankton)))+
  coord_cartesian(xlim=c(-0.3,30.3),
                  ylim=c(0,0.2625),
                  expand=F)+
  pubtheme+
  theme(legend.position="none")+
  facet_wrap(~factor(month,
                     levels=c("June","August","October")))

################################################################################
## Get model outputs together to visualize coefficient estimates

#################################################
## Fish consumption model
## First I'll clean up the slope parameters (tl_z and gzs_cpe_z)
fish1=tidy(fishmod,conf.int=T,conf.method="Wald")%>%
  filter(effect=="fixed",term %in% c("tl_z","gzs_cpe_z"))%>%
  mutate(model="fish",
         conditional_r2=r2_nakagawa(fishmod)[[1]],
         marginal_r2=r2_nakagawa(fishmod)[[2]])%>%
  select(model,term,est=estimate,se=std.error,lci=conf.low,
         uci=conf.high,p=p.value,conditional_r2,marginal_r2)

## Now I'll clean up the categorical parameter (month)
fish2=emmeans(fishmod,specs="month")%>%
  contrast(method="eff")%>%
  summary(infer=T)%>%
  as.data.frame()%>%
  mutate(model="fish",
         conditional_r2=r2_nakagawa(fishmod)[[1]],
         marginal_r2=r2_nakagawa(fishmod)[[2]])%>%
  select(model,term=contrast,est=estimate,se=SE,lci=asymp.LCL,
         uci=asymp.UCL,p=p.value,conditional_r2,marginal_r2)

## Then combine and we have complete model parameter output and r2 values
fishout=bind_rows(fish1,fish2)

#################################################
## Shad consumption model
## First I'll clean up the slope parameters (tl_z and gzs_cpe_z)
gzs1=tidy(gzsmod,conf.int=T,conf.method="Wald")%>%
  filter(effect=="fixed",term %in% c("tl_z","gzs_cpe_z"))%>%
  mutate(model="gzs",
         conditional_r2=r2_nakagawa(gzsmod)[[1]],
         marginal_r2=r2_nakagawa(gzsmod)[[2]])%>%
  select(model,term,est=estimate,se=std.error,lci=conf.low,
         uci=conf.high,p=p.value,conditional_r2,marginal_r2)

## Now I'll clean up the categorical parameter (month)
gzs2=emmeans(gzsmod,specs="month")%>%
  contrast(method="eff")%>%
  summary(infer=T)%>%
  as.data.frame()%>%
  mutate(model="gzs",
         conditional_r2=r2_nakagawa(gzsmod)[[1]],
         marginal_r2=r2_nakagawa(gzsmod)[[2]])%>%
  select(model,term=contrast,est=estimate,se=SE,lci=asymp.LCL,
         uci=asymp.UCL,p=p.value,conditional_r2,marginal_r2)

## Then combine and we have complete model parameter output and r2 values
gzsout=bind_rows(gzs1,gzs2)

#################################################
## Centrarchid consumption model
## First I'll clean up the slope parameters (tl_z and gzs_cpe_z)
centrarchid1=tidy(centrarchidmod,conf.int=T,conf.method="Wald")%>%
  filter(effect=="fixed",term %in% c("tl_z","gzs_cpe_z"))%>%
  mutate(model="centrarchid",
         conditional_r2=r2_nakagawa(centrarchidmod)[[1]],
         marginal_r2=r2_nakagawa(centrarchidmod)[[2]])%>%
  select(model,term,est=estimate,se=std.error,lci=conf.low,
         uci=conf.high,p=p.value,conditional_r2,marginal_r2)

## Now I'll clean up the categorical parameter (month)
centrarchid2=emmeans(centrarchidmod,specs="month")%>%
  contrast(method="eff")%>%
  summary(infer=T)%>%
  as.data.frame()%>%
  mutate(model="centrarchid",
         conditional_r2=r2_nakagawa(centrarchidmod)[[1]],
         marginal_r2=r2_nakagawa(centrarchidmod)[[2]])%>%
  select(model,term=contrast,est=estimate,se=SE,lci=asymp.LCL,
         uci=asymp.UCL,p=p.value,conditional_r2,marginal_r2)

## Then combine and we have complete model parameter output and r2 values
centrarchidout=bind_rows(centrarchid1,centrarchid2)

#################################################
## Crayfish consumption model
## First I'll clean up the slope parameters (tl_z and gzs_cpe_z)
crayfish1=tidy(crayfishmod,conf.int=T,conf.method="Wald")%>%
  filter(effect=="fixed",term %in% c("tl_z","gzs_cpe_z"))%>%
  mutate(model="crayfish",
         conditional_r2=r2_nakagawa(crayfishmod)[[1]],
         marginal_r2=r2_nakagawa(crayfishmod)[[2]])%>%
  select(model,term,est=estimate,se=std.error,lci=conf.low,
         uci=conf.high,p=p.value,conditional_r2,marginal_r2)

## Now I'll clean up the categorical parameter (month)
crayfish2=emmeans(crayfishmod,specs="month")%>%
  contrast(method="eff")%>%
  summary(infer=T)%>%
  as.data.frame()%>%
  mutate(model="crayfish",
         conditional_r2=r2_nakagawa(crayfishmod)[[1]],
         marginal_r2=r2_nakagawa(crayfishmod)[[2]])%>%
  select(model,term=contrast,est=estimate,se=SE,lci=asymp.LCL,
         uci=asymp.UCL,p=p.value,conditional_r2,marginal_r2)

## Then combine and we have complete model parameter output and r2 values
crayfishout=bind_rows(crayfish1,crayfish2)

#################################################
## Aquatic insect consumption model
## First I'll clean up the slope parameters (tl_z and gzs_cpe_z)
aquaticinsect1=tidy(aquaticinsectmod,conf.int=T,conf.method="Wald")%>%
  filter(effect=="fixed",term %in% c("tl_z","gzs_cpe_z"))%>%
  mutate(model="aquaticinsect",
         conditional_r2=r2_nakagawa(aquaticinsectmod)[[1]],
         marginal_r2=r2_nakagawa(aquaticinsectmod)[[2]])%>%
  select(model,term,est=estimate,se=std.error,lci=conf.low,
         uci=conf.high,p=p.value,conditional_r2,marginal_r2)

## Now I'll clean up the categorical parameter (month)
aquaticinsect2=emmeans(aquaticinsectmod,specs="month")%>%
  contrast(method="eff")%>%
  summary(infer=T)%>%
  as.data.frame()%>%
  mutate(model="aquaticinsect",
         conditional_r2=r2_nakagawa(aquaticinsectmod)[[1]],
         marginal_r2=r2_nakagawa(aquaticinsectmod)[[2]])%>%
  select(model,term=contrast,est=estimate,se=SE,lci=asymp.LCL,
         uci=asymp.UCL,p=p.value,conditional_r2,marginal_r2)

## Then combine and we have complete model parameter output and r2 values
aquaticinsectout=bind_rows(aquaticinsect1,aquaticinsect2)

#################################################
## Zooplankton consumption model
## First I'll clean up the slope parameters (tl_z and gzs_cpe_z)
zooplankton1=tidy(zooplanktonmod,conf.int=T,conf.method="Wald")%>%
  filter(effect=="fixed",term %in% c("tl_z","gzs_cpe_z"))%>%
  mutate(model="zooplankton",
         conditional_r2=r2_nakagawa(zooplanktonmod)[[1]],
         marginal_r2=r2_nakagawa(zooplanktonmod)[[2]])%>%
  select(model,term,est=estimate,se=std.error,lci=conf.low,
         uci=conf.high,p=p.value,conditional_r2,marginal_r2)

## Now I'll clean up the categorical parameter (month)
zooplankton2=emmeans(zooplanktonmod,specs="month")%>%
  contrast(method="eff")%>%
  summary(infer=T)%>%
  as.data.frame()%>%
  mutate(model="zooplankton",
         conditional_r2=r2_nakagawa(zooplanktonmod)[[1]],
         marginal_r2=r2_nakagawa(zooplanktonmod)[[2]])%>%
  select(model,term=contrast,est=estimate,se=SE,lci=asymp.LCL,
         uci=asymp.UCL,p=p.value,conditional_r2,marginal_r2)

## Then combine and we have complete model parameter output and r2 values
zooplanktonout=bind_rows(zooplankton1,zooplankton2)

################################################################################
## Combine table outputs and export
out=bind_rows(fishout,gzsout,centrarchidout,
              crayfishout,aquaticinsectout,zooplanktonout)%>%
  mutate(sig=case_when(p<=0.05 ~ "Significant",
                       TRUE ~ "Not significant"),
         term=factor(term,
                     levels=c("tl_z","gzs_cpe_z","June effect",
                              "August effect","October effect"),
                     labels=c("Scaled\nLargemouth\nBass TL","Scaled\nGizzard\nShad CPE","June","August","October")),
         model=factor(model,
                      levels=c("fish","gzs","centrarchid","crayfish","aquaticinsect","zooplankton"),
                      labels=c("Fish","Gizzard Shad","Centrarchidae",
                               "Crayfish","Aquatic insect","Zooplankton")))

## Create data frame with conditional and marginal r2 values so we can put them on plots
r2_labs=out%>%
  distinct(model,conditional_r2,marginal_r2)%>%
  mutate(marg_lab=paste0("'Marginal'~italic(R)^2==",sprintf("%.2f",marginal_r2)), 
         cond_lab=paste0("'Conditional'~italic(R)^2==",sprintf("%.2f",conditional_r2)))

## Create plot
ggplot(out,aes(x=term,y=est,ymin=lci,ymax=uci,color=sig))+
  geom_hline(yintercept=0)+
  geom_pointrange(size=1.2)+
  scale_color_manual(values=c("gray","#D94701"))+
  labs(x="Model parameter",y="Model parameter estimate")+
  geom_text(data=r2_labs,
            aes(label=marg_lab),
            x=5.3,
            y=-2.8,
            hjust=1,
            vjust=0,
            parse=T,
            inherit.aes=F,
            size=6, 
            color="black")+
  geom_text(data=r2_labs,
            aes(label=cond_lab),
            x=5.3,
            y=-3.5,
            hjust=1,
            vjust=0,
            parse=T,
            inherit.aes=F,
            size=6, 
            color="black")+
  facet_wrap(~model,ncol=2,dir="v")+
  pubtheme+
  theme(legend.position="none")
ggsave(plot=last_plot(),"fig3 - parmests.png",height=9,width=16,bg="white")

################################################################################
## Finally, we'll combine prediction plots into one figure
fig4=fish_plot/gzs_plot/centrarchid_plot|crayfish_plot/aquaticinsect_plot/zoop_plot
ggsave(plot=fig4,"fig4 - modpreds.png",height=9,width=16,bg="white")

################################################################################