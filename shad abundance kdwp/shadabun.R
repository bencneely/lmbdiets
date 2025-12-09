#--------------------------------------------------------------
#Ben Neely
#11/03/2025
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
         spp=species)

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
## Get gizzard shad CPE for each location/year
out1=out%>%
  group_by(impd,year)%>%
  summarize(nets=n(),
            gzs_cpe=mean(net_cpe),
            gzs_sd=sd(net_cpe),
            .groups="drop")

## Simplify and export so we can tie it into diet data
out2=out1%>%
  select(impd,year,gzs_cpe)

write_csv(out2,"shad abundance kdwp/gzs_cpe.csv")