#--------------------------------------------------------------
#Ben Neely
#11/18/2025
#Generate information for table 1
#--------------------------------------------------------------

## Clear R
cat("\014")  
rm(list=ls())

## Install and load packages
## Checks if package is installed, installs if not, activates for current session
if("tidyverse" %in% rownames(installed.packages()) == FALSE) {install.packages("tidyverse")}
library(tidyverse)

if("gt" %in% rownames(installed.packages()) == FALSE) {install.packages("gt")}
library(gt)

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
## The mollusc line just fixes a typo in the database
lmb=dat%>%
  filter(Species=="MICSAL",
         `Fish TL (mm)`>0)%>%
  select(id=StomachName,
         tl=`Fish TL (mm)`,
         fineitem=`Subcategory`,
         coarseitem=`Overarching Category`)%>%
  mutate(coarseitem=case_when(coarseitem=="Mollusc" ~ "Mollusca",TRUE ~ coarseitem),
         size_class=case_when(tl<150 ~ "small",
                              tl>=150 & tl<250 ~ "medium",
                              tl>=250 ~ "large"))

## Extract impd, month, year, and species from id
##Explanation of the regex:
#^(.{4})/: First 4 characters = impoundment
#_([A-Z]+): All capital letters after first underscore = month
#(\\d{4}): 4-digit year
#.*_: Skip over middle content (like CP2)
#(.{6})$: Last 6 characters = species
lmb1=lmb%>%
  extract(id,into=c("impd","month","year","spp"),
          regex="^(.{4})_([A-Z]+)(\\d{4}).*_(.{6})$",remove=FALSE)

## Clean up a typo
lmb1=lmb1%>%
  mutate(impd=case_when(impd=="OSTL" ~ "OTSL",
                        TRUE ~ impd))

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
         macroinvert=case_when(coarseitem=="Macroinvertebrate" ~ 1, TRUE ~ 0),
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

## Parse out so we just have what we need
lmb6=lmb5%>%
  select(impd,month,tl,gzs_cpe)

## Get min and max LMB TL by site and month
lmb7=lmb6%>%
  group_by(impd,month)%>%
  summarize(min_tl=min(tl),
            max_tl=max(tl),
            .groups="drop")%>%
  pivot_wider(names_from=month,
              values_from=c(min_tl,max_tl),
              names_sep="_")%>%
  select(impd,
         min_tl_June,max_tl_June,
         min_tl_August,max_tl_August,
         min_tl_October,max_tl_October)

## Get min and max Gizzard Shad CPE by site
gzs1=lmb6%>%
  group_by(impd)%>%
  summarize(min_gzs_cpe=round(min(gzs_cpe),2),
            max_gzs_cpe=round(max(gzs_cpe),2),
            .groups="drop")%>%
  mutate(across(.cols=c(min_gzs_cpe,max_gzs_cpe),
                .fns=~sprintf("%.2f",.)))

## Now we'll tie in Secchi data
secchi=read_csv("table 1/secchi.csv")%>%
  mutate(impd=factor(lakecode,
                     levels=c("BUSL","CASL","CLSL","DGSL",
                              "GESL","MESL","NOSL","OTSL",
                              "PTSA","SCSL","SNSL","WSSL"),
                     labels=c("Butler","Clark","Cowley","Douglas",
                              "Geary","Meade","Neosho","Ottawa",
                              "Pottawatomie","Scott","Shawnee","Washington")))%>%
  select(impd,month,secchi_m)

## Get min and max Secchi for each location and month
secchi1=secchi%>%
  group_by(impd)%>%
  summarize(min_secchi=min(secchi_m)*100,
            max_secchi=max(secchi_m)*100)%>%
  select(impd,min_secchi,max_secchi)

################################################################################
out=lmb7%>%
  inner_join(gzs1,by="impd")%>%
  inner_join(secchi1,by="impd")

write_csv(out,"table 1/table data.csv")