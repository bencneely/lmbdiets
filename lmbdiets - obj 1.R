#--------------------------------------------------------------
#Ben Neely
#03/18/2026
#Frequency of diet items
#--------------------------------------------------------------

## Clear R
cat("\014")  
rm(list=ls())

## Install and load packages
## Checks if package is installed, installs if not, activates for current session
if("tidyverse" %in% rownames(installed.packages()) == FALSE) {install.packages("tidyverse")}
library(tidyverse)

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
        legend.position.inside=c(0.99,0.01),
        legend.justification=c("right","bottom"),
        legend.title=element_text(size=16),
        legend.text=element_text(size=14),
        legend.background=element_rect(fill=NA))
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

################################################################################
## Top panel will have coarse items in the diet
unique(lmb1$coarseitem)

## We need to pare this down to mostly match the original iteration
## I lumped zooplanktons, combined oligochaeta and terrestrial invert, and lumped snail with mollusca
lmb2=lmb1%>%
  mutate(coarseitem=case_when(coarseitem=="Detritus & Sediment" ~ "Detritus",
                              coarseitem=="Oligochaeta" |
                                coarseitem=="Terrestrial invert matter" ~ "Terrestrial invertebrate",
                              coarseitem=="Ostracod" |
                                coarseitem=="Cladocera" |
                                coarseitem=="Zooplankton" ~ "Zooplankton",
                              coarseitem=="Snail"|
                                coarseitem=="Mollusca" ~ "Mollusk",
                              coarseitem=="Too degraded to analyze" ~ "Unidentified",
                              coarseitem=="Terrestrial plant matter" ~ "Terrestrial plant",
                              coarseitem=="Macroinvertebrate" ~ "Aquatic insect",
                              TRUE ~ coarseitem))
unique(lmb2$coarseitem)
## Now we have 14 coarse items including Empty

## We also have a lot of instances where a fish is listed twice for a coarseitem
## because it had two qualifying fineitems
## For example, a single fish that ate Hemiptera and Diptera would get two lines of macroinvertebrate
lmb2%>%
  group_by(id,coarseitem) %>%
  filter(n()>1)%>%
  arrange(id,coarseitem)

## Keep only unique coarse items per fish
## This converts the data from "Total Count" to "Presence/Absence"
lmb3=lmb2%>%
  distinct(id,coarseitem,.keep_all=T)

## Create dataframe with counts of each coarseitem by fish size
lmb4=lmb3%>%
  count(coarseitem,size_class)%>%
  mutate(coarseitem=factor(coarseitem),
         size_class=factor(size_class,levels=c("small","medium","large")))%>%
  complete(coarseitem,size_class,fill=list(n=0))%>%
  group_by(coarseitem)%>%
  mutate(tot_n=sum(n))%>%
  ungroup()

## Look at how many diet items we saw for each size group
xtabs(~size_class,lmb3)

## Create plot
broad=ggplot(lmb4,aes(x=n,
                      y=fct_reorder(coarseitem,tot_n),
                      fill=fct_relevel(size_class,"small","medium","large")))+
  geom_col(position=position_stack(reverse=T))+
  scale_fill_manual(values=c("#FDD0A2","#FD8D3C","#D94701"),
                    labels=c("< 150 mm: N = 1698",
                             "150 to 249 mm: N = 817",
                             "≥ 250 mm: N = 1583"),
                    name="Size group and count")+
  scale_x_continuous(breaks=seq(0,1200,200),
                     name="Count")+
  scale_y_discrete(name="")+
  coord_cartesian(xlim=c(-5,1250),
                  ylim=c(0.25,14.75),
                  expand=F)+
  pubtheme

################################################################################
## Bottom panel will break out fish species
lmb5=lmb2%>%
  filter(coarseitem=="Fish")

## See what we're working with
xtabs(~fineitem,lmb5)%>%
  as.data.frame()%>%
  arrange(desc(Freq))

## Let's break this down into species if possible, genus if not, and family if necessary
lmb6=lmb5%>%
  mutate(fineitem=case_when(fineitem=="Fish UNID" |
                              fineitem=="Fish tissue UNID"|
                              fineitem=="Fish"|
                              fineitem=="Fish matter UNID" ~ "Unknown fish",
                            fineitem=="Labidsthes" | 
                              fineitem=="Labidesthes" ~ "Brook Silverside",
                            fineitem=="Centrarchidae" ~ "Unknown Centrarchidae",
                            fineitem=="Ictaluridae" ~ "Unknown Ictaluridae",
                            fineitem=="Micropterus" ~ "Largemouth Bass",
                            fineitem=="Dorosoma" ~ "Gizzard Shad",
                            fineitem=="Notemigonus" ~ "Golden Shiner",
                            fineitem=="Gambusia" ~ "Western Mosquitofish",
                            fineitem=="Ictalurus" ~ "Channel Catfish",
                            fineitem=="Pylodictis" ~ "Flathead Catfish",
                            fineitem=="Cyprinus" ~ "Common Carp",
                            TRUE ~ fineitem))

## Look again
xtabs(~fineitem,lmb6)%>%
  as.data.frame()%>%
  arrange(desc(Freq))
## 14 fish genus looks good

## Create dataframe with counts of each fineitem by fish size
lmb7=lmb6%>%
  count(fineitem,size_class)%>%
  mutate(fineitem=factor(fineitem),
         size_class=factor(size_class,levels=c("small","medium","large")))%>%
  complete(fineitem,size_class,fill=list(n=0))%>%
  group_by(fineitem)%>%
  mutate(tot_n=sum(n))%>%
  ungroup()

## Look at how many diets we saw for each size group
xtabs(~size_class,lmb6)

## Create plot
fine=ggplot(lmb7,aes(x=n,
                     y=fct_reorder(fineitem,tot_n),
                     fill=fct_relevel(size_class,"small","medium","large")))+
  geom_col(position=position_stack(reverse=T))+
  scale_fill_manual(values=c("#FDD0A2","#FD8D3C","#D94701"),
                    labels=c("< 150 mm: N = 381",
                             "150 to 249 mm: N = 327",
                             "≥ 250 mm: N = 470"),
                    name="Size group and count")+
  scale_x_continuous(breaks=seq(0,600,100),
                     name="Count")+
  scale_y_discrete(name="")+
  coord_cartesian(xlim=c(-2.5,625),
                  ylim=c(0.25,14.75),
                  expand=F)+
  pubtheme

################################################################################
## Combine plots
out=broad/fine
ggsave(plot=out,"fig2 - dietitems.png",width=10,height=10,bg="white")