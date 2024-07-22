# Mesh block level, plotting observed SIR and fitted SIR (with and without Covariate)

library(tidyverse)
library(cowplot)

## Load Packages ## -----------------------------------------------------------
library(haven)
library(dplyr)
library(sf)
library(tidyverse)
library(readr)
library(tmap)
library(spdep)

library(CARBayes)
library(coda)

library(rgdal)			# For readOGR()
library(ggplot2)		# For ggplot()
library(rgeos)      # For unionSpatialPolygons(); gIntersection()
library(maptools)   # For unionSpatialPolygons()
library(gridExtra)	# For grid.arrange()
library(scales) 

## Load mesh data ## -----------------------------------------------------------
shape_loc <- "Data/Shapefiles_Mesh_SA1_SA2_SA3_SA4/MB_2016_QLD.shp"
mesh_raw1 <- st_read(shape_loc)
mesh_raw <- st_read(shape_loc) %>% 
  mutate(MB_CODE16 = as.numeric(as.character(MB_CODE16)),
         SA2_MAIN16 = as.numeric(as.character(SA2_MAIN16)),
         SA2_5DIG16 = as.numeric(as.character(SA2_5DIG16)))
mesh_pop <- read.csv("Data/2016 census mesh block counts.csv") %>% 
  mutate(MB_CODE_2016 = as.numeric(MB_CODE_2016)) %>% 
  rename(mesh_population = Person)


## Link mesh block population to mesh_raw (Join mesh_raw with mesh_pop) ##----------------
mesh <- mesh_raw %>% 
  inner_join(., mesh_pop, by = c("MB_CODE16"="MB_CODE_2016"))


## Add mesh population for each of SA2 values (check population count) ## --------------------------------------------------
#calculating sa2 population (sa2_mesh_pop_sum) from mesh totals per SA2 (not using the sa2_pop) file. 

#Calculated mesh proportions where the sum of proportions will now be exactly 1.
sa2_mesh_pop_sum <- mesh %>% 
  st_drop_geometry() %>% 
  group_by(SA2_MAIN16) %>% 
  summarise(mesh_pop_sum = sum(mesh_population))


## Calculate mesh block proportions within each SA2 ## --------------------------
mesh <- mesh %>% 
  inner_join(.,sa2_mesh_pop_sum, by = c("SA2_MAIN16"="SA2_MAIN16"))%>% 
  mutate(mesh_proportion = mesh_population/mesh_pop_sum) 


## mesh_proportion check ##-----------------------------------------------------
mesh_prop_check <- mesh %>%
  st_drop_geometry() %>% 
  group_by(SA2_MAIN16) %>%
  summarise(mesh_prop_sum=sum(mesh_proportion))


## Load simulated cancer data ## ----------------------------------------------
aca <- read_csv("Data/Simulated lung cancer data for QLD.csv")
conc <- read_csv("Data/sa2_2011_to_sa2_2016.csv")


##collapse the sex variables & group by using only sa2_name ##-------------------------
aca_total <- aca %>% 
  dplyr::select(sa2, sa2_name, sex, count, expect) %>% 
  group_by(sa2_name) %>% 
  summarise(sa2_count = sum(count),
            sa2_expect = sum(expect))


## use 2016 SA2 codes instead of 2011 ##-----------------------------------------
aca_SA2_16 <- aca_total %>% 
  inner_join(.,conc, by = c("sa2_name"="SA2_NAME_2011")) %>% 
  dplyr::select(sa2_count, sa2_expect, SA2_MAINCODE_2016, SA2_NAME_2016, sa2_name, RATIO) %>%
  mutate(sa2_count=sa2_count*RATIO,
         sa2_expect=sa2_expect*RATIO) %>%
  group_by(SA2_MAINCODE_2016, SA2_NAME_2016) %>%
  summarise(sa2_count=sum(sa2_count),
            sa2_expect=sum(sa2_expect))


## Calculate cancer count and expect at mesh block level (Link aca to mesh)##---------------------------
DATA <- mesh %>% 
  inner_join(.,aca_SA2_16, by = c("SA2_MAIN16"="SA2_MAINCODE_2016",
                                  "SA2_NAME16"="SA2_NAME_2016")) %>% 
  mutate(mesh_count = sa2_count*mesh_proportion,
         mesh_expect = sa2_expect*mesh_proportion)


## exclude missing values ##--------------------------------------------
DATA<- DATA %>% 
  filter(is.na(mesh_count) == FALSE) %>%
  filter(is.na(mesh_expect) == FALSE) 

## Descriptive statistics ##
summary(DATA$mesh_count)
sd(DATA$mesh_count)


## Creating spatial object ##------------------------------------------
coords <- st_transform(DATA, 3112) #https://epsg.io/?q=australia
DATA_sf <- st_as_sf(DATA, coords)


##SIR Calculation ###
DATA_sf$mesh_expect[DATA_sf$mesh_expect==0]<-0.0001
DATA_sf$SIR<-as.numeric(DATA_sf$mesh_count/DATA_sf$mesh_expect)
count_map<-ggplot(DATA_sf) + geom_sf(aes(fill= mesh_count),colour="black") +
  #scale_fill_continuous(trans="log10")+
  scale_fill_viridis_c(name="",option="inferno", direction = -1, 
                       limits= range(DATA_sf$mesh_count))+
  theme_bw()+#theme(legend.text=element_text(size = 14))
  theme(legend.position = "none")
main_map <- 
  count_map +
  geom_rect(
    xmin = 152.6,
    ymin = -28,
    xmax = 153.6,
    ymax = -27,
    fill = NA, 
    colour = "black",
    size = 0.6
  )

main_map %>% 
  ggdraw() +
  draw_plot(
    {
      main_map + 
        coord_sf(
          xlim = c(152.6, 153.6),
          ylim = c(-28, -27),
          expand = FALSE) +
        labs(title= "Brisbane")+
        theme(legend.position = "none",
              axis.ticks = element_blank(),
              axis.text = element_blank())
    },
    x = 0.65, 
    y = 0,
    width = 0.45, 
    height = 0.45)
SIR_map<-ggplot(DATA_sf) + geom_sf(aes(fill = SIR),color=NA) +
  scale_fill_distiller(limits=c(0,1.6),name="",type = "div", palette = 7, direction = -1)+
  theme_bw()+theme(legend.position = "none")


main_map <- 
  SIR_map +
  geom_rect(
    xmin = 152.6,
    ymin = -28,
    xmax = 153.6,
    ymax = -27,
    fill = NA, 
    colour = "black",
    size = 0.6
  )
main_map %>% 
  ggdraw() +
  draw_plot(
    {
      main_map + 
        coord_sf(
          xlim = c(152.6, 153.6),
          ylim = c(-28, -27),
          expand = FALSE) +
        labs(title= "Brisbane")+
        theme(legend.position = "none",axis.ticks = element_blank(),axis.text = element_blank())
    },
    x = 0.65, 
    y = 0,
    width = 0.45, 
    height = 0.45)

#fitted SIR calculation without covariate
mesh_wocov<-readOGR("mesh_WithoutCov_DATA_sf_model.shp")
mesh_wocov_data<-mesh_wocov@data
DATA_sf$fit_SIR_wocov<-mesh_wocov_data$fitted_SIR
SIR_map2<-ggplot(DATA_sf) + geom_sf(aes(fill = fit_SIR_wocov),color="NA") +
  scale_fill_distiller(limits=c(0,1.6),name="",type = "div", palette = 7, direction = -1)+
  theme_bw()+theme(legend.position = "none")


main_map2 <- 
  SIR_map2 +
  geom_rect(
    xmin = 152.6,
    ymin = -28,
    xmax = 153.6,
    ymax = -27,
    fill = NA, 
    colour = "black",
    size = 0.6
  )
main_map2 %>% 
  ggdraw() +
  draw_plot(
    {
      main_map2 + 
        coord_sf(
          xlim = c(152.6, 153.6),
          ylim = c(-28, -27),
          expand = FALSE) +
        labs(title= "Brisbane")+
        theme(legend.position = "none",axis.ticks = element_blank(),axis.text = element_blank())
    },
    x = 0.65, 
    y = 0,
    width = 0.45, 
    height = 0.45)
#fitted SIR calculation with covariate
mesh_cov<-readOGR("mesh_WithCov_DATA_sf_model.shp")
mesh_cov_data<-mesh_cov@data
mesh_cov_data<-rename(mesh_cov_data,c("MB_CODE16"="MB_CODE"))
MB_codes<-cbind(DATA_sf$MB_CODE16,mesh_cov_data$MB_CODE16)
length(which(MB_codes[,1] != MB_codes[,2]))
MB_codes[192,1:2]
MB_codes2<-matrix(0,nrow=67047,ncol=2)
MB_codes2[,2]<-c(MB_codes[1:191,2],30002180000,MB_codes[192:67046,2])
MB_codes2[,1]<-DATA_sf$MB_CODE16
length(which(MB_codes2[,1] != MB_codes2[,2]))
MB_codes2[1607:1611,]
mesh_cov_data2<-rbind(mesh_cov_data[1:191,c(1,5)],c(30002180000,0))
mesh_cov_data2<-rbind(mesh_cov_data2,mesh_cov_data[192:1608,c(1,5)],c(30023930000,0),
                      mesh_cov_data[1609:67045,c(1,5)])

DATA_sf$ftt_SIR<-mesh_cov_data2$ftt_SIR
SIR_map3<-ggplot(DATA_sf) + geom_sf(aes(fill = ftt_SIR),color="NA") +
  scale_fill_distiller(limits=c(0,2),name="",type = "div", palette = 7, direction = -1)+
  theme_bw()+theme(legend.position = "none")


main_map3 <- 
  SIR_map3 +
  geom_rect(
    xmin = 152.6,
    ymin = -28,
    xmax = 153.6,
    ymax = -27,
    fill = NA, 
    colour = "black",
    size = 0.6
  )
main_map3 %>% 
  ggdraw() +
  draw_plot(
    {
      main_map3 + 
        coord_sf(
          xlim = c(152.6, 153.6),
          ylim = c(-28, -27),
          expand = FALSE) +
        labs(title= "Brisbane")+
        theme(legend.position = "none",axis.ticks = element_blank(),axis.text = element_blank())
    },
    x = 0.65, 
    y = 0,
    width = 0.45, 
    height = 0.45)

# SA1 maps
## Load Mesh & SA1 data ## -----------------------------------------------------------
shape_loc <- "Data/Shapefiles_Mesh_SA1_SA2_SA3_SA4/SA1_2016_AUST.shp"
SA1_raw <- st_read(shape_loc)
SA1_raw_qld <- SA1_raw %>% #Taking records for Queensland only
  filter(STATE_NAME == "Queensland")%>%
  mutate(STATE_CODE = as.numeric(as.character(STATE_CODE)))

shape_loc <- "Data/Shapefiles_Mesh_SA1_SA2_SA3_SA4/MB_2016_QLD.shp"
mesh_raw <- st_read(shape_loc) %>% 
  st_drop_geometry() %>% 
  mutate(MB_CODE16 = as.numeric(as.character(MB_CODE16)),
         SA1_MAIN16 = as.numeric(as.character(SA1_MAIN16)),
         SA2_MAIN16 = as.numeric(as.character(SA2_MAIN16)),
         SA3_CODE16 = as.numeric (as.character(SA3_CODE16)),
         SA4_CODE16 = as.numeric (as.character(SA4_CODE16)))

#Rename 'person' column         
mesh_pop <- read_csv("Data/2016 census mesh block counts.csv") %>% 
  mutate(MB_CODE_2016 = as.numeric(MB_CODE_2016)) %>% 
  rename(mesh_population = Person)

## Link mesh block population to mesh_raw (Join mesh_raw with mesh_pop) ##----------------
mesh <- mesh_raw %>% 
  inner_join(., mesh_pop, by = c("MB_CODE16"="MB_CODE_2016"))%>%
  dplyr::select(MB_CODE16,SA1_MAIN16, SA2_MAIN16, SA3_CODE16, SA4_CODE16, mesh_population)

## Add mesh population for each of SA1 values ## --------------------------------------------------
##Calculate sa1 population
sa1_mesh_pop_sum <- mesh %>% 
  group_by(SA1_MAIN16) %>% 
  summarise(sa1_pop = sum(mesh_population))

## Add mesh population for each of SA2 values ## --------------------------------------------------
##Calculate sa2 population
sa2_mesh_pop_sum <- mesh %>% 
  group_by(SA2_MAIN16) %>% 
  summarise(sa2_pop = sum(mesh_population))

##Calculate sa1 proportions for each SA2
SA1_level <- mesh %>%
  distinct(SA1_MAIN16, SA2_MAIN16, SA3_CODE16, SA4_CODE16) %>%
  inner_join(.,sa1_mesh_pop_sum, by = c("SA1_MAIN16"="SA1_MAIN16"))%>%
  inner_join(.,sa2_mesh_pop_sum, by = c("SA2_MAIN16"="SA2_MAIN16")) %>%
  mutate(sa1_prop = sa1_pop/sa2_pop)

## Load simulated cancer data ## ----------------------------------------------
aca<- read_csv("Data/Simulated lung cancer data for QLD.csv")
conc <- read_csv("Data/sa2_2011_to_sa2_2016.csv")

##collapse the variables & group by using only sa2_name ##-------------------------
aca_total <- aca %>% 
  dplyr::select(sa2, sa2_name,  count, expect) %>% 
  group_by(sa2_name) %>% 
  summarise(sa2_count = sum(count),
            sa2_expect = sum(expect))

## use 2016 SA2 codes instead of 2011 ##-----------------------------------------
aca_SA2_16 <- aca_total %>% 
  inner_join(.,conc, by = c("sa2_name"="SA2_NAME_2011")) %>% 
  dplyr::select(sa2_count, sa2_expect, SA2_MAINCODE_2016, SA2_NAME_2016, sa2_name, RATIO) %>%
  mutate(sa2_count=sa2_count*RATIO,
         sa2_expect=sa2_expect*RATIO) %>%
  group_by(SA2_MAINCODE_2016, SA2_NAME_2016) %>%
  summarise(sa2_count=sum(sa2_count),
            sa2_expect=sum(sa2_expect))


##Join SA1_level with aca_SA2_16 and calculate sa1_count, sa1_expect using sa1_prop
SA1_level_count <- SA1_level %>%
  inner_join(.,aca_SA2_16, by = c("SA2_MAIN16"="SA2_MAINCODE_2016")) %>%
  mutate(sa1_count = sa2_count*sa1_prop,
         sa1_expect = sa2_expect*sa1_prop)
SA1_level_count$SA1_MAIN16<- as.character(SA1_level_count$SA1_MAIN16)
##Join SA1_raw_qld with SA1_level_count
DATA_SA1 <- SA1_raw_qld %>% 
  inner_join(.,SA1_level_count, by = c("SA1_MAIN16"="SA1_MAIN16"))


## Exclude missing values ##--------------------------------------------
DATA<- DATA_SA1 %>% 
  filter(is.na(sa1_count) == FALSE) %>%
  filter(is.na(sa1_expect) == FALSE) 


## Descriptive statistics ##
summary(DATA$sa1_count)
sd(DATA$sa1_count)

## Creating spatial object ##------------------------------------------
coords <- st_transform(DATA, 3112) #https://epsg.io/?q=australia
DATA_sf <- st_as_sf(DATA, coords)

## map SA1 cancer counts ##
DATA_sf$Cancer_counts<- as.numeric(DATA_sf$sa1_count)

## Creating spatial object ##------------------------------------------
coords <- st_transform(DATA, 3112) #https://epsg.io/?q=australia
DATA_sf <- st_as_sf(DATA, coords)

# count map
count_map<-ggplot(DATA_sf) + geom_sf(aes(fill=sa1_count),colour="black") +
  #scale_fill_continuous(trans="log10")+
  scale_fill_viridis_c(name="",option="inferno", direction = -1, 
                       limits= range(DATA_sf$sa1_count))+
  theme_bw()+#theme(legend.text=element_text(size = 14))
  theme(legend.position = "none")
main_map <- 
  count_map +
  geom_rect(
    xmin = 152.6,
    ymin = -28,
    xmax = 153.6,
    ymax = -27,
    fill = NA, 
    colour = "black",
    size = 0.6
  )

main_map %>% 
  ggdraw() +
  draw_plot(
    {
      main_map + 
        coord_sf(
          xlim = c(152.6, 153.6),
          ylim = c(-28, -27),
          expand = FALSE) +
        labs(title= "Brisbane")+
        theme(legend.position = "none",
              axis.ticks = element_blank(),
              axis.text = element_blank())
    },
    x = 0.65, 
    y = 0,
    width = 0.45, 
    height = 0.45)
##SIR Calculation ###
DATA_sf$sa1_expect[DATA_sf$sa1_expect==0]<-0.001
DATA_sf$SIR<-as.numeric(DATA_sf$sa1_count/DATA_sf$sa1_expect)
SIR_map<-ggplot(DATA_sf) + geom_sf(aes(fill = SIR),color=NA) +
  scale_fill_distiller(limits=c(0,1.6),name="",type = "div", palette = 7, direction = -1)+
  theme_bw()#+theme(legend.position = "none")


main_map <- 
  SIR_map +
  geom_rect(
    xmin = 152.6,
    ymin = -28,
    xmax = 153.6,
    ymax = -27,
    fill = NA, 
    colour = "black",
    size = 0.6
  )
main_map %>% 
  ggdraw() +
  draw_plot(
    {
      main_map + 
        coord_sf(
          xlim = c(152.6, 153.6),
          ylim = c(-28, -27),
          expand = FALSE) +
        labs(title= "Brisbane")+
        theme(legend.position = "none",axis.ticks = element_blank(),axis.text = element_blank())
    },
    x = 0.65, 
    y = 0,
    width = 0.45, 
    height = 0.45)

#fitted SIR calculation without covariate
load("Bym_Without_Cov_SA1.Rdata")
fit_SIR<-as.data.frame(cbind("SA1_MAIN16"=DATA_sf_model$SA1_MAIN16,"SIR"=Fit_obs$fitted))
library(dplyr)
DATA_for_join<- as.data.frame(DATA_sf)
DATA_map<-merge(DATA_sf,fit_SIR, by= "SA1_MAIN16",all.x = T)

DATA_map$SIR.y[is.na(DATA_map$SIR.y)]<- DATA_map$SIR.x[is.na(DATA_map$SIR.y)]

library(tidyverse)
SIR_map2<-ggplot(DATA_map) + geom_sf(aes(fill = SIR.y),color="NA") +
  scale_fill_distiller(limits=c(0,1.6),name="",type = "div", palette = 7, direction = -1)+
  theme_bw()+theme(legend.position = "none")


main_map2 <- 
  SIR_map2 +
  geom_rect(
    xmin = 152.6,
    ymin = -28,
    xmax = 153.6,
    ymax = -27,
    fill = NA, 
    colour = "black",
    size = 0.6
  )
main_map2 %>% 
  ggdraw() +
  draw_plot(
    {
      main_map2 + 
        coord_sf(
          xlim = c(152.6, 153.6),
          ylim = c(-28, -27),
          expand = FALSE) +
        labs(title= "Brisbane")+
        theme(legend.position = "none",axis.ticks = element_blank(),axis.text = element_blank())
    },
    x = 0.65, 
    y = 0,
    width = 0.45, 
    height = 0.45)
#fitted SIR calculation with covariate
load("Bym_With_Cov_SA1 (1).Rdata")
fit_SIR<-as.data.frame(cbind("SA1_MAIN16"=DATA_sf_model$SA1_MAIN16,"SIR"=Fit_obs$fitted))
library(dplyr)
DATA_map<-merge(DATA_sf,fit_SIR, by= "SA1_MAIN16",all.x = T)
DATA_map$SIR[is.na(DATA_map$SIR)]<- DATA_map$`SIR (Lung) at SA1`[is.na(DATA_map$SIR)]
SIR_map3<-ggplot(DATA_map) + geom_sf(aes(fill = SIR),color="NA") +
  scale_fill_distiller(limits=c(0,2),name="",type = "div", palette = 7, direction = -1)+
  theme_bw()#+theme(legend.position = "none")
main_map3 <- 
  SIR_map3 +
  geom_rect(
    xmin = 152.6,
    ymin = -28,
    xmax = 153.6,
    ymax = -27,
    fill = NA, 
    colour = "black",
    size = 0.6
  )
main_map3 %>% 
  ggdraw() +
  draw_plot(
    {
      main_map3 + 
        coord_sf(
          xlim = c(152.6, 153.6),
          ylim = c(-28, -27),
          expand = FALSE) +
        labs(title= "Brisbane")+
        theme(legend.position = "none",axis.ticks = element_blank(),axis.text = element_blank())
    },
    x = 0.65, 
    y = 0,
    width = 0.45, 
    height = 0.45)
# SA2 level 

## Load SA2 data ## -----------------------------------------------------------
shape_loc <- "Data/Shapefiles_Mesh_SA1_SA2_SA3_SA4/SA2_2016_AUST.shp"
SA2_raw <- st_read(shape_loc)
SA2_raw_qld <- SA2_raw %>% #Taking records for Queensland only
  filter(STATE_NAME == "Queensland")%>%
  mutate(SA2_MAIN16 = as.numeric(as.character(SA2_MAIN16)))


## Collapse the sex variables & group by using only sa2_name ##-------------------------
aca_total <- aca %>% 
  dplyr::select(sa2, sa2_name, count, expect,irsdscore) %>% 
  group_by(sa2_name,irsdscore) %>% 
  summarise(sa2_count = sum(count),
            sa2_expect = sum(expect))

aca_cov <- aca_total

## Checking for duplicates ##---------------------------------------------------
length(unique(aca_total$sa2_name))

## use 2016 SA2 codes instead of 2011 ##-----------------------------------------
aca_SA2_16 <- aca_total %>% 
  inner_join(.,conc, by = c("sa2_name"="SA2_NAME_2011")) %>% 
  dplyr::select(sa2_count, sa2_expect, SA2_MAINCODE_2016, SA2_NAME_2016,sa2_name, RATIO) %>%
  mutate(sa2_count=sa2_count*RATIO,
         sa2_expect=sa2_expect*RATIO) %>%
  group_by(SA2_MAINCODE_2016, SA2_NAME_2016) %>%
  summarise(sa2_count=sum(sa2_count),
            sa2_expect=sum(sa2_expect))

aca_with_cov <- aca_SA2_16 %>%
  inner_join(.,aca_cov,by = c("SA2_NAME_2016"="sa2_name")) %>%
  dplyr::select(SA2_MAINCODE_2016,SA2_NAME_2016,sa2_count=sa2_count.x, sa2_expect=sa2_expect.x,irsdscore)

##Join SA2_raw_qld with aca_with_cov
DATA_SA2 <- SA2_raw_qld %>% 
  inner_join(.,aca_with_cov, by = c("SA2_MAIN16"="SA2_MAINCODE_2016"))


## exclude missing values ##--------------------------------------------
DATA<- DATA_SA2 %>% 
  filter(is.na(sa2_count) == FALSE) %>%
  filter(is.na(irsdscore)== FALSE)


## Calculate quintile ##----------------------------------------------
SA2_irsd = quantile(DATA$irsdscore, probs = seq(0, 1, .2))
DATA$SA2_irsd <- DATA$irsdscore
DATA$SA2_irsd[DATA$irsdscore<=SA2_irsd[2]] <- 1
DATA$SA2_irsd[DATA$irsdscore>SA2_irsd[2] & DATA$irsdscore<=SA2_irsd[3]] <- 2 
DATA$SA2_irsd[DATA$irsdscore>SA2_irsd[3] & DATA$irsdscore<=SA2_irsd[4]] <- 3
DATA$SA2_irsd[DATA$irsdscore>SA2_irsd[4] & DATA$irsdscore<=SA2_irsd[5]] <- 4
DATA$SA2_irsd[DATA$irsdscore>SA2_irsd[5]] <- 5 


## Descriptive statistics ##
summary(DATA$sa2_count)
sd(DATA$sa2_count)

## Creating spatial object ##------------------------------------------
coords <- st_transform(DATA, 3112) #https://epsg.io/?q=australia
DATA_sf <- st_as_sf(DATA, coords)

# ## map the sa2_count distribution using a quantile classification scheme ##-------------------
# tm_shape(DATA_sf) + tm_polygons(style="quantile", col = "sa2_count") +
#   tm_legend(outside = TRUE, text.size = .8) 


## map SA2 cancer counts ##
DATA_sf$Cancer_counts<- as.numeric(DATA_sf$sa2_count)
# count map
# count map
count_map<-ggplot(DATA_sf) + geom_sf(aes(fill=sa2_count),colour="black") +
  #scale_fill_continuous(trans="log10")+
  scale_fill_viridis_c(name="",option="inferno", direction = -1, 
                       limits= range(DATA_sf$sa2_count))+
  theme_bw()+#theme(legend.text=element_text(size = 14))
  theme(legend.position = "none")
main_map <- 
  count_map +
  geom_rect(
    xmin = 152.6,
    ymin = -28,
    xmax = 153.6,
    ymax = -27,
    fill = NA, 
    colour = "black",
    size = 0.6
  )

main_map %>% 
  ggdraw() +
  draw_plot(
    {
      main_map + 
        coord_sf(
          xlim = c(152.6, 153.6),
          ylim = c(-28, -27),
          expand = FALSE) +
        labs(title= "Brisbane")+
        theme(legend.position = "none",
              axis.ticks = element_blank(),
              axis.text = element_blank())
    },
    x = 0.65, 
    y = 0,
    width = 0.45, 
    height = 0.45)

##SIR Calculation ###
DATA_SA2$SIR<-as.numeric(DATA_SA2$sa2_count/DATA_SA2$sa2_expect)
SIR_map<-ggplot(DATA_SA2) + geom_sf(aes(fill = SIR),color=NA) +
  scale_fill_distiller(limits=c(0,1.6),name="",type = "div", palette = 7, direction = -1)+
  theme_bw()+theme(legend.position = "none")
main_map <- 
  SIR_map +
  geom_rect(
    xmin = 152.6,
    ymin = -28,
    xmax = 153.6,
    ymax = -27,
    fill = NA, 
    colour = "black",
    size = 0.6
  )
main_map %>% 
  ggdraw() +
  draw_plot(
    {
      main_map + 
        coord_sf(
          xlim = c(152.6, 153.6),
          ylim = c(-28, -27),
          expand = FALSE) +
        labs(title= "Brisbane")+
        theme(legend.position = "none",axis.ticks = element_blank(),axis.text = element_blank())
    },
    x = 0.65, 
    y = 0,
    width = 0.45, 
    height = 0.45)

#fitted SIR calculation without covariate
load("BYM_wo_cov_SA2.RData")
library(dplyr)
DATA_SA2$fitted_SIR<-model_bym_wo_cov$fitted.values/DATA_sf$sa2_expect


library(tidyverse)
SIR_map2<-ggplot(DATA_sf) + geom_sf(aes(fill = fitted_SIR),color="NA") +
  scale_fill_distiller(limits=c(0,1.6),name="",type = "div", palette = 7, direction = -1)+
  theme_bw()+theme(legend.position = "none")


main_map2 <- 
  SIR_map2 +
  geom_rect(
    xmin = 152.6,
    ymin = -28,
    xmax = 153.6,
    ymax = -27,
    fill = NA, 
    colour = "black",
    size = 0.6
  )
main_map2 %>% 
  ggdraw() +
  draw_plot(
    {
      main_map2 + 
        coord_sf(
          xlim = c(152.6, 153.6),
          ylim = c(-28, -27),
          expand = FALSE) +
        labs(title= "Brisbane")+
        theme(legend.position = "none",axis.ticks = element_blank(),axis.text = element_blank())
    },
    x = 0.65, 
    y = 0,
    width = 0.45, 
    height = 0.45)
#fitted SIR calculation with covariate
load("BYM_with_cov._SA2.RData")
DATA_sf$fit_SIR_cov<-model_bym$fitted.values/DATA_sf$sa2_expect 
SIR_map3<-ggplot(DATA_sf) + geom_sf(aes(fill = fit_SIR_cov),color="NA") +
  scale_fill_distiller(limits=c(0,2),name="",type = "div", palette = 7, direction = -1)+
  theme_bw()+theme(legend.position = "none")
main_map3 <- 
  SIR_map3 +
  geom_rect(
    xmin = 152.6,
    ymin = -28,
    xmax = 153.6,
    ymax = -27,
    fill = NA, 
    colour = "black",
    size = 0.6
  )
main_map3 %>% 
  ggdraw() +
  draw_plot(
    {
      main_map3 + 
        coord_sf(
          xlim = c(152.6, 153.6),
          ylim = c(-28, -27),
          expand = FALSE) +
        labs(title= "Brisbane")+
        theme(legend.position = "none",axis.ticks = element_blank(),axis.text = element_blank())
    },
    x = 0.65, 
    y = 0,
    width = 0.45, 
    height = 0.45)
# SA3

## Load SA3 data ## -----------------------------------------------------------
shape_loc <- "Data/Shapefiles_Mesh_SA1_SA2_SA3_SA4/SA3_2016_AUST.shp"
SA3_raw <- st_read(shape_loc)
SA3_raw_qld <- SA3_raw %>% #Taking records for Queensland only
  filter(STE_NAME16 == "Queensland")%>%
  mutate(SA3_CODE16 = as.numeric(as.character(SA3_CODE16)))

##Load mesh data
shape_loc <- "Data/Shapefiles_Mesh_SA1_SA2_SA3_SA4/MB_2016_QLD.shp"
mesh_raw <- st_read(shape_loc) %>% 
  st_drop_geometry() %>% 
  mutate(SA1_MAIN16 = as.numeric(as.character(SA1_MAIN16)),
         SA2_MAIN16 = as.numeric(as.character(SA2_MAIN16)),
         SA3_CODE16 = as.numeric (as.character(SA3_CODE16)))%>%
  dplyr::select(SA1_MAIN16,SA2_MAIN16,SA3_CODE16)     

## Collapse the sex variables & group by using only sa2_name ##-------------------------
aca_total <- aca %>% 
  dplyr::select(sa2, sa2_name, sex, count, expect) %>% 
  group_by(sa2_name) %>% 
  summarise(sa2_count = sum(count),
            sa2_expect = sum(expect))

## Use 2016 SA2 codes instead of 2011 ##-----------------------------------------
aca_SA2_16 <- aca_total %>% 
  inner_join(.,conc, by = c("sa2_name"="SA2_NAME_2011")) %>% 
  dplyr::select(sa2_count, sa2_expect, SA2_MAINCODE_2016, SA2_NAME_2016,sa2_name, RATIO) %>%
  mutate(sa2_count=sa2_count*RATIO,
         sa2_expect=sa2_expect*RATIO) %>%
  group_by(SA2_MAINCODE_2016, SA2_NAME_2016) %>%
  summarise(sa2_count=sum(sa2_count),
            sa2_expect=sum(sa2_expect))

##Calculate sa3_count and sa3_expect (link aca_SA2_16 and mesh_raw)
SA3_level_count <- mesh_raw %>%
  distinct(SA2_MAIN16, SA3_CODE16) %>%
  inner_join(.,aca_SA2_16, by = c("SA2_MAIN16"="SA2_MAINCODE_2016"))%>%
  group_by(SA3_CODE16)%>%
  summarise(sa3_count=sum(sa2_count),
            sa3_expect=sum(sa2_expect))

## Load SA1_IRSD data ##
SA1_IRSD <- read_csv("Data/SA1_IRSD.csv")
## Join mesh_raw with SA1_IRSD
SA1_IRSD_full <- SA1_IRSD %>% 
  inner_join(.,mesh_raw, by = c("SA1_MAIN16"="SA1_MAIN16")) 

## Calculate SA3 population ##
sa3_pop <- SA1_IRSD_full %>% 
  group_by(SA3_CODE16) %>% 
  summarise(sa3_pop = sum(SA1_pop))

SA1_IRSD_full <- SA1_IRSD_full %>%
  inner_join(.,sa3_pop,by = c("SA3_CODE16"= "SA3_CODE16"))

SA1_IRSD_calc <- SA1_IRSD_full %>%
  mutate(IRSD_SA3=(SA1_irsdscore*SA1_pop)/sa3_pop)

SA3_IRSD <- SA1_IRSD_calc %>%
  group_by(SA3_CODE16) %>%
  summarise(SA3_irsdscore=round(sum(IRSD_SA3)))

## Calculate quintile ##----------------------------------------------
SA3_irsd = quantile(SA3_IRSD$SA3_irsdscore, probs = seq(0, 1, .2))
SA3_IRSD$SA3_irsd <- SA3_IRSD$SA3_irsdscore
SA3_IRSD$SA3_irsd[SA3_IRSD$SA3_irsdscore<=SA3_irsd[2]] <- 1
SA3_IRSD$SA3_irsd[SA3_IRSD$SA3_irsdscore>SA3_irsd[2] & SA3_IRSD$SA3_irsdscore<=SA3_irsd[3]] <- 2 
SA3_IRSD$SA3_irsd[SA3_IRSD$SA3_irsdscore>SA3_irsd[3] & SA3_IRSD$SA3_irsdscore<=SA3_irsd[4]] <- 3
SA3_IRSD$SA3_irsd[SA3_IRSD$SA3_irsdscore>SA3_irsd[4] & SA3_IRSD$SA3_irsdscore<=SA3_irsd[5]] <- 4
SA3_IRSD$SA3_irsd[SA3_IRSD$SA3_irsdscore>SA3_irsd[5]] <- 5 


write.table(SA3_IRSD)

## Join SA3_IRSD with SA3_level_count and SA3_raw_qld ##
DATA_SA3 <- SA3_level_count %>% 
  inner_join(., SA3_IRSD, by = c("SA3_CODE16"="SA3_CODE16")) %>% 
  inner_join(.,SA3_raw_qld, by = c ("SA3_CODE16"="SA3_CODE16"))


## Descriptive statistics ##
summary(DATA_SA3$sa3_count)
sd(DATA_SA3$sa3_count)


#Create simple feature object
DATA_sf <- st_sf(DATA_SA3)

## Map sa3 count ##
DATA_sf$Cancer_counts<- as.numeric(DATA_sf$sa3_count)

# count map
count_map<-ggplot(DATA_sf) + geom_sf(aes(fill=sa3_count),colour="black") +
  #scale_fill_continuous(trans="log10")+
  scale_fill_viridis_c(name="",option="inferno", direction = -1, 
                       limits= range(DATA_sf$sa3_count))+
  theme_bw()+#theme(legend.text=element_text(size = 14))
  theme(legend.position = "none")
main_map <- 
  count_map +
  geom_rect(
    xmin = 152.6,
    ymin = -28,
    xmax = 153.6,
    ymax = -27,
    fill = NA, 
    colour = "black",
    size = 0.6
  )

main_map %>% 
  ggdraw() +
  draw_plot(
    {
      main_map + 
        coord_sf(
          xlim = c(152.6, 153.6),
          ylim = c(-28, -27),
          expand = FALSE) +
        labs(title= "Brisbane")+
        theme(legend.position = "none",
              axis.ticks = element_blank(),
              axis.text = element_blank())
    },
    x = 0.65, 
    y = 0,
    width = 0.45, 
    height = 0.45)

##SIR Calculation ###
DATA_sf$SIR<-as.numeric(DATA_sf$sa3_count/DATA_sf$sa3_expect)
SIR_map<-ggplot(DATA_sf) + geom_sf(aes(fill = SIR),color=NA) +
  scale_fill_distiller(limits=c(0,1.6),name="",type = "div", palette = 7, direction = -1)+
  theme_bw()+theme(legend.position = "none")
main_map <- 
  SIR_map +
  geom_rect(
    xmin = 152.6,
    ymin = -28,
    xmax = 153.6,
    ymax = -27,
    fill = NA, 
    colour = "black",
    size = 0.6
  )
main_map %>% 
  ggdraw() +
  draw_plot(
    {
      main_map + 
        coord_sf(
          xlim = c(152.6, 153.6),
          ylim = c(-28, -27),
          expand = FALSE) +
        labs(title= "Brisbane")+
        theme(legend.position = "none",axis.ticks = element_blank(),axis.text = element_blank())
    },
    x = 0.65, 
    y = 0,
    width = 0.45, 
    height = 0.45)

#fitted SIR calculation without covariate
load("SA3_BYM_noCov.RData")
library(dplyr)
DATA_sf$fitted_SIR<-model.BYM_noCov$fitted.values/DATA_sf$sa3_expect


library(tidyverse)
SIR_map2<-ggplot(DATA_sf) + geom_sf(aes(fill = fitted_SIR),color="NA") +
  scale_fill_distiller(limits=c(0,1.6),name="",type = "div", palette = 7, direction = -1)+
  theme_bw()+theme(legend.position = "none")


main_map2 <- 
  SIR_map2 +
  geom_rect(
    xmin = 152.6,
    ymin = -28,
    xmax = 153.6,
    ymax = -27,
    fill = NA, 
    colour = "black",
    size = 0.6
  )
main_map2 %>% 
  ggdraw() +
  draw_plot(
    {
      main_map2 + 
        coord_sf(
          xlim = c(152.6, 153.6),
          ylim = c(-28, -27),
          expand = FALSE) +
        labs(title= "Brisbane")+
        theme(legend.position = "none",axis.ticks = element_blank(),axis.text = element_blank())
    },
    x = 0.65, 
    y = 0,
    width = 0.45, 
    height = 0.45)
#fitted SIR calculation with covariate
load("SA3_BYM_withCov.RData")
DATA_sf$fit_SIR_cov<-model_bym$fitted.values/DATA_sf$sa3_expect 
SIR_map3<-ggplot(DATA_sf) + geom_sf(aes(fill = fit_SIR_cov),color="NA") +
  scale_fill_distiller(limits=c(0,2),name="",type = "div", palette = 7, direction = -1)+
  theme_bw()+theme(legend.position = "none")
main_map3 <- 
  SIR_map3 +
  geom_rect(
    xmin = 152.6,
    ymin = -28,
    xmax = 153.6,
    ymax = -27,
    fill = NA, 
    colour = "black",
    size = 0.6
  )
main_map3 %>% 
  ggdraw() +
  draw_plot(
    {
      main_map3 + 
        coord_sf(
          xlim = c(152.6, 153.6),
          ylim = c(-28, -27),
          expand = FALSE) +
        labs(title= "Brisbane")+
        theme(legend.position = "none",axis.ticks = element_blank(),axis.text = element_blank())
    },
    x = 0.65, 
    y = 0,
    width = 0.45, 
    height = 0.45)
# SA4


## Load SA4 data ## -----------------------------------------------------------
shape_loc <-"Data/Shapefiles_Mesh_SA1_SA2_SA3_SA4/SA4_2016_AUST.shp"
SA4_raw <- st_read(shape_loc)
SA4_raw_qld <- SA4_raw %>% #Taking records for Queensland only
  filter(STATE_NAME == "Queensland")%>%
  mutate(SA4_CODE16 = as.numeric(as.character(SA4_CODE16)))



##Load mesh data
shape_loc <- "Data/Shapefiles_Mesh_SA1_SA2_SA3_SA4/MB_2016_QLD.shp"
mesh_raw <- st_read(shape_loc) %>% 
  st_drop_geometry() %>% 
  mutate(SA1_MAIN16 = as.numeric(as.character(SA1_MAIN16)),
         SA2_MAIN16 = as.numeric(as.character(SA2_MAIN16)),
         SA3_CODE16 = as.numeric (as.character(SA3_CODE16)),
         SA4_CODE16 = as.numeric (as.character(SA4_CODE16))) %>%
  dplyr::select(SA1_MAIN16,SA2_MAIN16,SA3_CODE16,SA4_CODE16)     

## Collapse the sex variables & group by using only sa2_name ##-------------------------
aca_total <- aca %>% 
  dplyr::select(sa2, sa2_name, sex, count, expect) %>% 
  group_by(sa2_name) %>% 
  summarise(sa2_count = sum(count),
            sa2_expect = sum(expect))

## Use 2016 SA2 codes instead of 2011 ##-----------------------------------------
aca_SA2_16 <- aca_total %>% 
  inner_join(.,conc, by = c("sa2_name"="SA2_NAME_2011")) %>% 
  dplyr::select(sa2_count, sa2_expect, SA2_MAINCODE_2016, SA2_NAME_2016,sa2_name, RATIO) %>%
  mutate(sa2_count=sa2_count*RATIO,
         sa2_expect=sa2_expect*RATIO) %>%
  group_by(SA2_MAINCODE_2016, SA2_NAME_2016) %>%
  summarise(sa2_count=sum(sa2_count),
            sa2_expect=sum(sa2_expect))


##Calculate sa4_count and sa4_expect (link aca_SA2_16 and mesh_raw)
SA4_level_count <- mesh_raw %>%
  distinct(SA2_MAIN16, SA3_CODE16, SA4_CODE16) %>%
  inner_join(.,aca_SA2_16, by = c("SA2_MAIN16"="SA2_MAINCODE_2016"))%>%
  group_by(SA4_CODE16)%>%
  summarise(sa4_count=sum(sa2_count),
            sa4_expect=sum(sa2_expect))


## Load SA1_IRSD data ##
SA1_IRSD <- read_csv("Data/SA1_IRSD.csv")
## Join mesh_raw with SA1_IRSD
SA1_IRSD_full <- SA1_IRSD %>% 
  inner_join(.,mesh_raw, by = c("SA1_MAIN16"="SA1_MAIN16")) 

## Calculate SA4 population ##
sa4_pop <- SA1_IRSD_full %>% 
  group_by(SA4_CODE16) %>% 
  summarise(sa4_pop = sum(SA1_pop))

SA1_IRSD_full <- SA1_IRSD_full %>%
  inner_join(.,sa4_pop,by = c("SA4_CODE16"= "SA4_CODE16"))

SA1_IRSD_calc <- SA1_IRSD_full %>%
  mutate(IRSD_SA4=(SA1_irsdscore*SA1_pop)/sa4_pop)

SA4_IRSD <- SA1_IRSD_calc %>%
  group_by(SA4_CODE16) %>%
  summarise(SA4_irsdscore=round(sum(IRSD_SA4)))

## Calculate quintile ##----------------------------------------------
SA4_irsd = quantile(SA4_IRSD$SA4_irsdscore, probs = seq(0, 1, .2))
SA4_IRSD$SA4_irsd <- SA4_IRSD$SA4_irsdscore
SA4_IRSD$SA4_irsd[SA4_IRSD$SA4_irsdscore<=SA4_irsd[2]] <- 1
SA4_IRSD$SA4_irsd[SA4_IRSD$SA4_irsdscore>SA4_irsd[2] & SA4_IRSD$SA4_irsdscore<=SA4_irsd[3]] <- 2 
SA4_IRSD$SA4_irsd[SA4_IRSD$SA4_irsdscore>SA4_irsd[3] & SA4_IRSD$SA4_irsdscore<=SA4_irsd[4]] <- 3
SA4_IRSD$SA4_irsd[SA4_IRSD$SA4_irsdscore>SA4_irsd[4] & SA4_IRSD$SA4_irsdscore<=SA4_irsd[5]] <- 4
SA4_IRSD$SA4_irsd[SA4_IRSD$SA4_irsdscore>SA4_irsd[5]] <- 5 

## Join SA4_IRSD with SA4_level_count and SA4_raw_qld ##
DATA_SA4 <- SA4_level_count %>% 
  inner_join(., SA4_IRSD, by = c("SA4_CODE16"="SA4_CODE16")) %>% 
  inner_join(.,SA4_raw_qld, by = c ("SA4_CODE16"="SA4_CODE16"))


## exclude missing values ##--------------------------------------------
DATA<- DATA_SA4 %>% 
  filter(is.na(sa4_count) == FALSE)


## Descriptive statistics ##
summary(DATA$sa4_count)
sd(DATA$sa4_count)

#Create simple feature object
DATA_sf <- st_sf(DATA)

## map SA4 cancer counts ##
DATA_sf$Cancer_counts<- as.numeric(DATA_sf$sa4_count)
# count map
count_map<-ggplot(DATA_sf) + geom_sf(aes(fill=sa4_count),colour="black") +
  #scale_fill_continuous(trans="log10")+
  scale_fill_viridis_c(name="",option="inferno", direction = -1, 
                       limits= range(DATA_sf$sa4_count))+
  theme_bw()+#theme(legend.text=element_text(size = 14))
  theme(legend.position = "none")
main_map <- 
  count_map +
  geom_rect(
    xmin = 152.6,
    ymin = -28,
    xmax = 153.6,
    ymax = -27,
    fill = NA, 
    colour = "black",
    size = 0.6
  )

main_map %>% 
  ggdraw() +
  draw_plot(
    {
      main_map + 
        coord_sf(
          xlim = c(152.6, 153.6),
          ylim = c(-28, -27),
          expand = FALSE) +
        labs(title= "Brisbane")+
        theme(legend.position = "none",
              axis.ticks = element_blank(),
              axis.text = element_blank())
    },
    x = 0.65, 
    y = 0,
    width = 0.45, 
    height = 0.45)
##SIR Calculation ###
DATA_sf$SIR<-as.numeric(DATA_sf$sa4_count/DATA_sf$sa4_expect)
SIR_map<-ggplot(DATA_sf) + geom_sf(aes(fill = SIR),color=NA) +
  scale_fill_distiller(limits=c(0,1.6),name="",type = "div", palette = 7, direction = -1)+
  theme_bw()+theme(legend.position = "none")
main_map <- 
  SIR_map +
  geom_rect(
    xmin = 152.6,
    ymin = -28,
    xmax = 153.6,
    ymax = -27,
    fill = NA, 
    colour = "black",
    size = 0.6
  )
main_map %>% 
  ggdraw() +
  draw_plot(
    {
      main_map + 
        coord_sf(
          xlim = c(152.6, 153.6),
          ylim = c(-28, -27),
          expand = FALSE) +
        labs(title= "Brisbane")+
        theme(legend.position = "none",axis.ticks = element_blank(),axis.text = element_blank())
    },
    x = 0.65, 
    y = 0,
    width = 0.45, 
    height = 0.45)

#fitted SIR calculation without covariate
load("BYM_SA4_nocov.RData")
library(dplyr)
DATA_sf$fitted_SIR<-model.BYM_noCov$fitted.values/DATA_sf$sa4_expect


library(tidyverse)
SIR_map2<-ggplot(DATA_sf) + geom_sf(aes(fill = fitted_SIR),color="NA") +
  scale_fill_distiller(limits=c(0,1.6),name="",type = "div", palette = 7, direction = -1)+
  theme_bw()+theme(legend.position = "none")


main_map2 <- 
  SIR_map2 +
  geom_rect(
    xmin = 152.6,
    ymin = -28,
    xmax = 153.6,
    ymax = -27,
    fill = NA, 
    colour = "black",
    size = 0.6
  )
main_map2 %>% 
  ggdraw() +
  draw_plot(
    {
      main_map2 + 
        coord_sf(
          xlim = c(152.6, 153.6),
          ylim = c(-28, -27),
          expand = FALSE) +
        labs(title= "Brisbane")+
        theme(legend.position = "none",axis.ticks = element_blank(),axis.text = element_blank())
    },
    x = 0.65, 
    y = 0,
    width = 0.45, 
    height = 0.45)
#fitted SIR calculation with covariate
load("BYM_withCov_SA4.RData")
DATA_sf$fit_SIR_cov<-model_bym$fitted.values/DATA_sf$sa4_expect 
SIR_map3<-ggplot(DATA_sf) + geom_sf(aes(fill = fit_SIR_cov),color="NA") +
  scale_fill_distiller(limits=c(0,2),name="",type = "div", palette = 7, direction = -1)+
  theme_bw()+theme(legend.position = "none")
main_map3 <- 
  SIR_map3 +
  geom_rect(
    xmin = 152.6,
    ymin = -28,
    xmax = 153.6,
    ymax = -27,
    fill = NA, 
    colour = "black",
    size = 0.6
  )
main_map3 %>% 
  ggdraw() +
  draw_plot(
    {
      main_map3 + 
        coord_sf(
          xlim = c(152.6, 153.6),
          ylim = c(-28, -27),
          expand = FALSE) +
        labs(title= "Brisbane")+
        theme(legend.position = "none",axis.ticks = element_blank(),axis.text = element_blank())
    },
    x = 0.65, 
    y = 0,
    width = 0.45, 
    height = 0.45)
