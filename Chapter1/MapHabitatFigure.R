rm(list=ls())
library(readr)
library(tidyr)
library(dplyr)
library(broom)
library(ggplot2)
library(ggthemes)
library(mapdata)
library(maps)
library(rgbif)
library(forcats)
library(CoordinateCleaner)
library(ggrepel)
library(png)
library(gridExtra)
library(colourpicker)
library(cowplot)

#########################################
#### Transect map and habitat figure ####
#########################################

#### Site map ####
sites <- read.csv("~/Dropbox/master_data/site/site_details.csv")
sites$'Elevation (m)' <- sites$Mean.Elev

(sites.map <- ggplot(sites, aes(x = Mean.Long, y = Mean.Lat, color=Mean.Elev)) +
  borders("worldHires", ylim = c(55, 58), xlim=c(-6,-2), colour = "gray60", fill = "gray100", size = 0.3) +
  coord_fixed(1.4)+
  geom_point(color = "black", size = 2.5) +
  xlab("Longitude")+
  ylab("Latitude")+
  theme_bw()+
  coord_cartesian(ylim = c(55.8, 58), xlim=c(-6,-2))+
  scale_color_gradient(name  ="Elevation (m)",low="gray70", high="black")+
  geom_point(alpha = 0.9, size = 2.3)+
  theme(text = element_text(size=15)))

#### Habitat per site by latitude ####
## Better annotation in Dataframe.R

habitat <- read.csv("~/Dropbox/2018/Habitats_2018_all.csv")
# combining other deciduous
habitat$X6_oth.decid <- rowSums(habitat[,c("X6_cherry", "X6_elder", "X6_holly", "X6_rose")], na.rm=TRUE)
habitat$X21_oth.decid <- rowSums(habitat[,c("X21_blackthorn", "X21_holly")], na.rm=TRUE)
habitat$s_oth.decid <- rowSums(habitat[,c("s_aspen", "s_cherry", "s_chestnut", "s_elder", "s_hawthorn", "s_holly", "s_lime", "s_whitebeam", "s_other")], na.rm=TRUE)
habitat$m_oth.decid <- rowSums(habitat[,c("m_aspen", "m_cherry", "m_chestnut", "m_holly", "m_lime", "m_whitebeam", "m_horsechestnut")], na.rm=TRUE)
habitat$l_oth.decid <- habitat$l_cherry
habitat$z_oth.decid <- rowSums(habitat[,c("z_blackthorn", "z_cherry", "z_other")], na.rm=TRUE)

# combining conifers
habitat$X6_conifer_all <- habitat$X6_juniper
habitat$s_conifer_all <- rowSums(habitat[,c("s_conifer", "s_pine", "s_yew", "s_juniper")], na.rm=TRUE)
habitat$m_conifer_all <- rowSums(habitat[,c("m_conifer", "m_pine", "m_yew")], na.rm=TRUE)
habitat$l_conifer_all <- rowSums(habitat[,c("l_conifer", "l_pine")], na.rm=TRUE)

# converting 6 stands into small trees
habitat$alder6smt <- habitat$X6_alder*0.5
habitat$ash6smt <- habitat$X6_ash*0.5
habitat$birch6smt <- habitat$X6_birch*0.5
habitat$beech6smt <- habitat$X6_beech*0.5
habitat$hazel6smt <- habitat$X6_hazel*0.5
habitat$rowan6smt <- habitat$X6_rowan*0.5
habitat$syc6smt <- habitat$X6_sycamore*0.5
habitat$willow6smt <- habitat$X6_willow*0.5
habitat$oth.decid6smt <- habitat$X6_oth.decid*0.5
habitat$conifer6smt <- habitat$X6_conifer_all*0.5
# thickets already as small trees

# combining stands, thickets and small trees
habitat$Alder_S <- rowSums(habitat[,c("s_alder", "X21_alder", "alder6smt")], na.rm=TRUE)
habitat$Ash_S <- rowSums(habitat[,c("s_ash", "ash6smt")], na.rm=TRUE)
habitat$Beech_S <- rowSums(habitat[,c("s_beech", "beech6smt")], na.rm=TRUE)
habitat$Birch_S <- rowSums(habitat[,c("s_birch", "X21_birch", "birch6smt")], na.rm=TRUE)
habitat$Hazel_S <- rowSums(habitat[,c("s_hazel", "X21_hazel", "hazel6smt")], na.rm=TRUE)
habitat$Rowan_S <- rowSums(habitat[,c("s_rowan", "rowan6smt")], na.rm=TRUE)
habitat$Sycamore_S <- rowSums(habitat[,c("s_sycamore", "syc6smt")], na.rm=TRUE)
habitat$Willow_S <- rowSums(habitat[,c("s_willow", "willow6smt", "X21_willow", "z_willow")], na.rm=TRUE)
habitat$Conifer_S <- rowSums(habitat[,c("X6_conifer_all", "s_conifer_all")], na.rm=TRUE)
habitat$OthDecid_S <- rowSums(habitat[,c("X6_oth.decid", "X21_oth.decid", "s_oth.decid", "z_oth.decid")], na.rm=TRUE)

##### Use the weighting function to calculate foliage per tree type and size at each NB

## making a function for the weighting equation
weighting <- function(x){return(pi*(x/(2*pi))^2)}  ##weighted by cross secitonal area of the trunk in cm^2 (min possible size for each category)
Sm <- weighting(40)
Med <- weighting(100)
Lar <- weighting(250)
S=1
M=Med/Sm
L=Lar/Sm

Habitat_byNB <- data.frame(Site=habitat$site, 
                           NB=habitat$nestbox, 
                           Alder_S=(S*habitat$Alder_S),
                           Alder_M=(M*habitat$m_alder),
                           Ash_L=(L*habitat$l_ash),
                           Ash_M=(M*habitat$m_ash),
                           Ash_S=(S*habitat$Ash_S),
                           Beech_S=(S*habitat$Beech_S),
                           Beech_M=(M*habitat$m_beech),
                           Beech_L=(L*habitat$l_beech),
                           Birch_S=(S*habitat$Birch_S),
                           Birch_M=(M*habitat$m_birch),
                           Birch_L=(L*habitat$l_birch),
                           Elm_S=(S*habitat$s_elm),
                           Elm_M=(M*habitat$m_elm),
                           Elm_L=(L*habitat$l_elm),
                           Hazel_S=(S*habitat$Hazel_S),
                           Oak_S=(S*habitat$s_oak),
                           Oak_M=(M*habitat$m_oak),
                           Oak_L=(L*habitat$l_oak), 
                           Rowan_S=(S*habitat$Rowan_S),
                           Rowan_M=(M*habitat$m_rowan),
                           Sycamore_S=(S*habitat$Sycamore_S),
                           Sycamore_M=(M*habitat$m_sycamore),
                           Sycamore_L=(L*habitat$l_sycamore),
                           Willow_S=(S*habitat$Willow_S),
                           Willow_M=(M*habitat$m_willow),
                           Willow_L=(L*habitat$l_willow),
                           Conifer_S=(S*habitat$Conifer_S),
                           Conifer_M=(M*habitat$m_conifer_all),
                           Conifer_L=(L*habitat$l_conifer_all),
                           OthDecid_S=(S*habitat$OthDecid_S),  
                           OthDecid_M=(M*habitat$m_oth.decid),
                           OthDecid_L=(L*habitat$l_oth.decid))
Habitat_byNB[is.na(Habitat_byNB)] <- 0

#Foliage score for each tree type at each NB
Habitat_byNB$Alder_FS <- rowSums(Habitat_byNB[,c("Alder_S", "Alder_M")], na.rm=TRUE)
Habitat_byNB$Ash_FS <- rowSums(Habitat_byNB[,c("Ash_S", "Ash_M", "Ash_L")], na.rm=TRUE)
Habitat_byNB$Beech_FS <- rowSums(Habitat_byNB[,c("Beech_S", "Beech_M", "Beech_L")], na.rm=TRUE)
Habitat_byNB$Birch_FS <- rowSums(Habitat_byNB[,c("Birch_S", "Birch_M", "Birch_L")], na.rm=TRUE)
Habitat_byNB$Elm_FS <- rowSums(Habitat_byNB[,c("Elm_S", "Elm_M", "Elm_L")], na.rm=TRUE)
Habitat_byNB$Hazel_FS <- Habitat_byNB$Hazel_S
Habitat_byNB$Oak_FS <- rowSums(Habitat_byNB[,c("Oak_S", "Oak_M", "Oak_L")], na.rm=TRUE)
Habitat_byNB$Rowan_FS <- rowSums(Habitat_byNB[,c("Rowan_S", "Rowan_M")], na.rm=TRUE)
Habitat_byNB$Sycamore_FS <- rowSums(Habitat_byNB[,c("Sycamore_S", "Sycamore_M", "Sycamore_L")], na.rm=TRUE)
Habitat_byNB$Willow_FS <- rowSums(Habitat_byNB[,c("Willow_S", "Willow_M", "Willow_L")], na.rm=TRUE)
Habitat_byNB$Conifer_FS <- rowSums(Habitat_byNB[,c("Conifer_S", "Conifer_M", "Conifer_L")], na.rm=TRUE)
Habitat_byNB$OtherDeciduous_FS <- rowSums(Habitat_byNB[,c("OthDecid_S", "OthDecid_M", "OthDecid_L")], na.rm=TRUE)

## combining NB within each site- mean to account for different number of NBs at sites
site <- read.csv("~/Dropbox/master_data/site/site_details.csv")
Habitat_Site <- Habitat_byNB[,35:46]
Habitat_Site$Site <- Habitat_byNB$Site
Habitat_Site <- aggregate(.~Site, Habitat_Site, mean)


#### Use FS's to make figure ####
#Make mean foliage scores long
Habitat_FS <- Habitat_Site[,1:13]
Habitat_FS_long <- gather(Habitat_FS, key="Taxon", value="FS", 2:13)
Habitat_FS_long$Taxon <- gsub("_FS","", Habitat_FS_long$Taxon)

#FS graph by latitude
site$Site <- site$site
FS_long_siteinfo <- merge(Habitat_FS_long, site, by="Site", duplicates.ok=TRUE)
FS_long_siteinfo <- FS_long_siteinfo[order(FS_long_siteinfo$Mean.Lat),] 
FS_long_siteinfo$Taxon <- as.factor(FS_long_siteinfo$Taxon)
FS_long_siteinfo$Taxon <- factor(FS_long_siteinfo$Taxon , c("Alder","Ash","Beech","Birch","Elm","Hazel","Oak","Rowan","Sycamore","Willow" ,"Conifer" ,"OtherDeciduous"))
AllTaxaCols <- c("darkred", "firebrick3", "chocolate2", "goldenrod", "olivedrab4","darkgreen",  "deepskyblue3", "royalblue4", "slateblue2", "orchid","gray57","gray35")

(habitat <- ggplot(FS_long_siteinfo, aes(fct_inorder(Site), FS))+
  geom_bar(aes(fill=Taxon), stat="identity",alpha = 0.7)+
  theme_bw()+
  ylab("Foliage Score")+
  theme(axis.text.x= element_text(angle=90))+
  scale_fill_manual(values=AllTaxaCols)+
  xlab("Site")+
  scale_y_continuous(expand = c(0,0))+
  theme(text = element_text(size=15)))

#Store legends from figures
legend1 <- get_legend(sites.map)
legend2 <- get_legend(habitat)

#redo figures without legend
sites.map <- ggplot(sites, aes(x = Mean.Long, y = Mean.Lat, color=Mean.Elev)) +
  borders("worldHires", ylim = c(55, 58), xlim=c(-6,-2), colour = "gray60", fill = "gray100", size = 0.3)+
  coord_fixed(1.4)+
  geom_point(color = "black", size = 2.5) +
  xlab("Longitude")+
  ylab("Latitude")+
  theme_bw()+
  coord_cartesian(ylim = c(55.8, 58), xlim=c(-6,-2))+
  scale_color_gradient(name  ="Elevation (m)",low="gray70", high="black")+
  geom_point(alpha = 0.9, size = 2.3)+
  theme(text = element_text(size=15),
        legend.position = "none")

habitat <- ggplot(FS_long_siteinfo, aes(fct_inorder(Site), FS))+
  geom_bar(aes(fill=Taxon), stat="identity",alpha = 0.7)+
  theme_bw()+
  ylab("Foliage Score")+
  theme(axis.text.x= element_text(angle=90))+
  scale_fill_manual(values=AllTaxaCols)+
  xlab("Site")+
  scale_y_continuous(expand = c(0,0))+
  theme(text = element_text(size=15))+
  theme(legend.position="none")

space <- ggplot() + theme_void()

legend <- grid.arrange(legend1, legend2,nrow=2,heights=c(0.9,1.5))
row1 <- grid.arrange(sites.map, space, habitat, space, legend, ncol=5, widths=c(1.7,0.1,2.3,0.1,0.5))
#saved as 6x14.5"