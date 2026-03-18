#--------------------------------------------------------------
#Ben Neely
#03/18/2026
#Calculate gizzard shad CPE from KDWP gill net data
#--------------------------------------------------------------

## Clear R
cat("\014")  
rm(list=ls())

## Install and load packages
## Checks if package is installed, installs if not, activates for current session
if("FSA" %in% rownames(installed.packages()) == FALSE) {install.packages("FSA")}
library(FSA)

if("tidyverse" %in% rownames(installed.packages()) == FALSE) {install.packages("tidyverse")}
library(tidyverse)

if("lubridate" %in% rownames(installed.packages()) == FALSE) {install.packages("lubridate")}
library(lubridate)

## Read in raw data
fish=read_csv("shad abundance kdwp/fish.csv")
samp=read_csv("shad abundance kdwp/samp.csv")

## Format fish data
fish1=fish%>%
  as.data.frame()%>%
  expandCounts(~group_count)%>%
  select(id=fish_sample_id,
         impd=impoundment_code,
         spp=species,
         tl=length,
         year)

## Format sample data
samp1=samp%>%
  as.data.frame()%>%
  mutate(eff=1)%>%
  select(id=fish_sample_id,
         impd=impoundment_code,
         year,
         eff)

################################################################################
## Get gizzard shad per id from fish1
fish2=fish1%>%
  group_by(id)%>%
  summarize(gzs=n(),
            .groups="drop")

## Add gizzard shad count to each sample
## Replace NA (no gizzard shad observed) with a true 0
## Calculate CPE for each net
out=samp1%>%
  left_join(fish2,by="id")%>%
  replace_na(list(gzs=0))%>%
  mutate(net_cpe=gzs/eff)

################################################################################
## Get size metrics
fish3=fish1%>%
  group_by(impd,year)%>%
  summarize(min_tl_yr=min(tl,na.rm=TRUE),
            max_tl_yr=max(tl,na.rm=TRUE),
            .groups="drop")

################################################################################
## Get gizzard shad CPE for each location
out1=out%>%
  group_by(impd,year)%>%
  summarize(nets=n(),
            gzs_cpe=mean(net_cpe),
            gzs_sd=sd(net_cpe),
            .groups="drop")%>%
  left_join(fish3,by=c("impd","year"))

## Simplify and export so we can tie it into diet data
out2=out1%>%
  select(impd,year,nets,gzs_cpe,min_tl_yr,max_tl_yr)

tmp=out2%>%
  group_by(impd)%>%
  summarize(min_tl=min(min_tl_yr,na.rm=T),
            max_tl=max(max_tl_yr,na.rm=T),
            .groups="drop")

out3=out2%>%
  left_join(tmp,by="impd")

write_csv(out3,"shad abundance kdwp/gzs_cpe.csv")

################################################################################
## Gill net effort
min(out2$nets)
max(out2$nets)
mean(out2$nets)
sd(out2$nets)
