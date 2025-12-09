#--------------------------------------------------------------
#Ben Neely
#11/18/2025
#Summary statistics
#--------------------------------------------------------------

## Clear R
cat("\014")  
rm(list=ls())

## Install and load packages
## Checks if package is installed, installs if not, activates for current session
if("janitor" %in% rownames(installed.packages()) == FALSE) {install.packages("janitor")}
library(janitor)

if("tidyverse" %in% rownames(installed.packages()) == FALSE) {install.packages("tidyverse")}
library(tidyverse)

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

################################################################################
################################################################################
## Summarize for first paragraph
## Number of individual fish diets
nrow(lmb5)
length(unique(lmb5$id)) # Double check one fish per line

## Number of fish per impd, year, and month
tmp=xtabs(~impd+year+month,lmb5)%>%
  as_tibble()

## Find the Clark zeroes
print(tmp,n=Inf)%>%
  arrange(n)

## Mean and sd for diets taken during each period
tmp1=tmp%>%
  filter(n>0)

mean(tmp1$n)
sd(tmp1$n)

################################################################################
################################################################################
## Summarize for second paragraph
## Diet item frequency overall and by size class (data for Figure 2)
## I lumped zooplanktons, combined oligochaeta and terrestrial invert, and lumped snail with mollusca
## Note that there are 3500 fish but 4479 rows since fish often had multiple diet items
lmb6=lmb1%>%
  mutate(coarseitem=case_when(coarseitem=="Detritus & Sediment" ~ "Detritus",
                              coarseitem=="Oligochaeta" |
                                coarseitem=="Terrestrial invert matter" ~ "Terrestrial invertebrate",
                              coarseitem=="Ostracod" |
                                coarseitem=="Cladocera" |
                                coarseitem=="Zooplankton" ~ "Zooplankton",
                              coarseitem=="Snail"|
                                coarseitem=="Mollusca" ~ "Mollusk",
                              coarseitem=="Too degraded to analyze" ~ "Unidentifiable",
                              coarseitem=="Terrestrial plant matter" ~ "Terrestrial plant",
                              TRUE ~ coarseitem))
unique(lmb6$coarseitem)
## Now we have 14 coarse items including Empty

## We also have a lot of instances where a fish is listed twice for a coarseitem
## because it had two qualifying fineitems
## For example, a single fish that ate Hemiptera and Diptera would get two lines of macroinvertebrate
lmb6%>%
  group_by(id,coarseitem) %>%
  filter(n()>1)%>%
  arrange(id,coarseitem)

## Keep only unique coarse items per fish
## This converts the data from "Total Count" to "Presence/Absence"
lmb7=lmb6%>%
  distinct(id,coarseitem,.keep_all=T)

## Empty diets
nrow(subset(lmb7,coarseitem=="Empty"))

## Number of fish with one diet item
lmb7%>%
  filter(coarseitem!="Empty")%>%
  count(id)%>%
  filter(n==1)%>%
  nrow()

## Number of fish with two different diet items
lmb7%>%
  count(id)%>% 
  filter(n==2)%>%
  nrow() 

## Number of fish with three or more different diet items
lmb7%>%
  count(id)%>% 
  filter(n>2)%>%
  nrow() 

## Create dataframe with counts of each coarseitem by fish size
lmb8=lmb7%>%
  count(coarseitem,size_class)%>%
  mutate(coarseitem=factor(coarseitem),
         size_class=factor(size_class,levels=c("small","medium","large")))%>%
  complete(coarseitem,size_class,fill=list(n=0))%>%
  group_by(coarseitem)%>%
  mutate(tot_n=sum(n))%>%
  ungroup()

## Look at total number of fish by size
lmb7%>%
  group_by(size_class)%>%
  summarise(n_fish=n_distinct(id))

## Unique items arranged by total abundance (tot_n)
lmb8%>%
  distinct(coarseitem, tot_n)%>%
  arrange(desc(tot_n))

## Small fish
lmb8%>%
  filter(size_class=="small")%>%
  arrange(desc(n))

## Medium fish
lmb8%>%
  filter(size_class=="medium")%>%
  arrange(desc(n))

## Large fish
lmb8%>%
  filter(size_class=="large")%>%
  arrange(desc(n))

################################################################################
################################################################################
## Summarize for third paragraph - Genus/family level consumption
lmb9=lmb6%>%
  filter(coarseitem=="Fish")

## Let's consolidate fineitem a bit and look at the data frame
lmb10=lmb9%>%
  mutate(fineitem=case_when(fineitem=="Fish UNID" |
                              fineitem=="Fish tissue UNID"|
                              fineitem=="Fish"|
                              fineitem=="Fish matter UNID" ~ "Unknown fish",
                            fineitem=="Labidsthes" ~ "Labidesthes",
                            fineitem=="Centrarchidae" ~ "Unknown Centrarchidae",
                            fineitem=="Ictaluridae" ~ "Unknown Ictaluridae",
                            TRUE ~ fineitem))
nrow(lmb10)
xtabs(~fineitem,lmb10)

## We filtered out records so only coarseitem=="Fish" remain
## However, there are 1178 records for 1096 fish
## This is because some fish had multiple fineitems that were "Fish"

## Let's start by finding the number of fish that ONLY have Unknown fish in their diets
lmb10%>%
  group_by(id)%>%
  filter(all(fineitem=="Unknown fish"))%>%
  distinct(id)%>%
  nrow()
## Proportion with only unknown fish
517/1096

## Number of fish with known fish diet items that have 1, 2, or 3+ types
lmb10%>%
  filter(fineitem != "Unknown fish")%>% 
  group_by(id)%>%
  summarize(n_items=n_distinct(fineitem))%>%
  count(diet_complexity=case_when(n_items==1 ~ "1",
                                  n_items==2 ~ "2",
                                  n_items>=3 ~ "3+"))

## Create dataframe with counts of each fineitem by fish size
lmb11=lmb10%>%
  count(fineitem,size_class)%>%
  mutate(fineitem=factor(fineitem),
         size_class=factor(size_class,levels=c("small","medium","large")))%>%
  complete(fineitem,size_class,fill=list(n=0))%>%
  group_by(fineitem)%>%
  mutate(tot_n=sum(n))%>%
  ungroup()

## Look at total number of fish by size
lmb10%>%
  group_by(size_class)%>%
  summarise(n_fish=n_distinct(id))

## Unique items arranged by total abundance (tot_n)
lmb11%>%
  distinct(fineitem, tot_n)%>%
  arrange(desc(tot_n))

## Small fish
lmb11%>%
  filter(size_class=="small")%>%
  arrange(desc(n))

## Medium fish
lmb11%>%
  filter(size_class=="medium")%>%
  arrange(desc(n))

## Large fish
lmb11%>%
  filter(size_class=="large")%>%
  arrange(desc(n))