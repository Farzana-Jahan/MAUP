### Fitting BYM model without covariate to SA4 levels data in QLD
# codes include data prep and modelling
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

## Load SA4 data ## -----------------------------------------------------------
shape_loc <- "Data/Shapefiles_Mesh_SA1_SA2_SA3_SA4/SA4_2016_AUST.shp"
SA4_raw <- st_read(shape_loc)
SA4_raw_qld <- SA4_raw %>% #Taking records for Queensland only
  filter(STE_CODE16 == "3")%>%
  mutate(SA4_CODE16 = as.numeric(as.character(SA4_CODE16)))


## Load simulated cancer data ## ----------------------------------------------
aca <- read_csv("Data/Simulated lung cancer data for QLD.csv")# not available in the repository, need to be requested to CCQ at: Statistics@cancerqld.org.au
conc <- read_csv("Data/sa2_2011_to_sa2_2016.csv")



## Subset aca data to all cancers ##--------------------------
aca <- aca_full %>% 
  filter(site10group == 23)

##Load mesh data
shape_loc <- "Data/Shapefiles_Mesh_SA1_SA2_SA3_SA4/MB_2016_QLD.shp"
mesh_raw <- st_read(shape_loc) %>% 
  st_drop_geometry() %>% 
  mutate(SA2_MAIN16 = as.numeric(as.character(SA2_MAIN16)),
         SA3_CODE16 = as.numeric (as.character(SA3_CODE16)),
         SA4_CODE16 = as.numeric (as.character(SA4_CODE16)))  %>%
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

# ##SIR Calculation ###

##SIR Calculation ###
DATA_sf$SIR<-as.numeric(DATA_sf$sa4_count/DATA_sf$sa4_expect)

## Fit model when data are aggregated at SA4 level ##--------------------------------

### To handle log0 situation as 0 cannot survive log transformation without being infinite ##----

DATA_sf$observed <- round(DATA_sf$sa4_count)
DATA_sf$expected <- round(DATA_sf$sa4_expect)
DATA_sf$expected[DATA_sf$expected==0] <- 0.001 


### Creating factor variables:covariates ###-----------------------------------------------

DATA_sf$SA4_irsd.f <- factor(DATA_sf$SA4_irsd)
is.factor(DATA_sf$SA4_irsd.f)
levels(DATA_sf$SA4_irsd.f)

## Define neighboring polygons (neighbor as being any contiguous polygon that shares at least one vertex)##-----------
nb <- poly2nb(DATA_sf, queen=FALSE)# nb lists all neighboring polygons
## Create binary, first-order, adjacency spatial weights matrix using neighbours list ##---------------------
nbp_w <- nb2mat(nb, style="B", zero.policy=TRUE)

## Areas with zero counts need to be excluded from the adjacency weights matrix ##---
num_neighbours <- rowSums(nbp_w)
to_remove_index <- which(num_neighbours == 0)
if(sum(to_remove_index)!=0){
  nbp_w <- nbp_w[-to_remove_index, -to_remove_index]
}
DATA_sf_model <- DATA_sf %>%
  dplyr::select(SA4_CODE16,observed,expected,SA4_irsd.f,geometry)
  

###### Fit BYM model with covariate #######
start <- Sys.time()
formula <-observed~offset(log(expected))+SA4_irsd.f

model_bym<- S.CARbym(formula= formula, data=DATA_sf_model, family="poisson", W=nbp_w,
                     burnin=500000, n.sample=1500000, thin=100)#10000 number of iterations

Sys.time() - start
save(model_bym, file="BYM_withCov_SA4.RData")

## Summary of results (summary output of the fitted model) ##---------
print(model_bym)
summary(model_bym)

## Check samples from the model ##----------------------------------------------
summary(model_bym$samples)

###Traces and density plot####
plot(model_bym$samples$tau2)
plot(model_bym$samples$sigma2)
plot(model_bym$samples$beta)


# Combine residuals and fitted values into a data frame
result = data.frame(fitted_values = model_bym$fitted.values, residuals = model_bym$residuals)
DATA_sf$fitted_SIR<-result$fitted_values/DATA_sf$sa4_expect


# Plot observed versus fitted values
Fit_obs <- data.frame(observed=DATA_sf_model$observed, fitted=model_bym$fitted.values)

## Moran's I to the residuals (compute Moran's I statistic and test to assess its significance)##--------------------------------------------
## Define neighboring polygons (neighbor as being any contiguous polygon that shares at least one vertex)##-----------
nb <- poly2nb(DATA_sf, queen=FALSE)# nb lists all neighboring polygons
## Assigning weights to each neighboring polygon ##--------------------------
nbp_w <- nb2listw(nb, style="B", zero.policy=TRUE)
res.model <- residuals(model_bym)
moran.mc(x = res.model,listw = nbp_w, nsim = 10000, zero.policy = TRUE)



