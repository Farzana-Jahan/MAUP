### Fitting BYM model without covariate to SA2 levels data in QLD
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
library(scales)     # For rescale()


## Load SA2 data ## -----------------------------------------------------------
shape_loc <- "Data/Shapefiles_Mesh_SA1_SA2_SA3_SA4/SA2_2016_AUST.shp"
SA2_raw <- st_read(shape_loc)
SA2_raw_qld <- SA2_raw %>% #Taking records for Queensland only
  filter(STE_CODE16 == "3")%>%
  mutate(SA2_MAIN16 = as.numeric(as.character(SA2_MAIN16)))


## Load simulated cancer data ## ----------------------------------------------
aca<- read_csv("Data/Simulated lung cancer data for QLD.csv")# not available in the repository, need to be requested to CCQ at: Statistics@cancerqld.org.au
conc <- read_csv("Data/sa2_2011_to_sa2_2016.csv")

##collapse the sex variables & group by using only sa2_name ##-------------------------
aca_total <- aca %>% 
  dplyr::select(sa2, sa2_name,sex,count, expect) %>% 
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


##Join SA2_raw_qld with aca_SA2_16
DATA_SA2 <- SA2_raw_qld %>% 
  inner_join(.,aca_SA2_16, by = c("SA2_MAIN16"="SA2_MAINCODE_2016"))


## exclude missing values ##--------------------------------------------
DATA<- DATA_SA2 %>% 
  filter(is.na(sa2_count) == FALSE)

## Descriptive statistics ##
summary(DATA$sa2_count)
sd(DATA$sa2_count)

## Creating spatial object ##------------------------------------------
coords <- st_transform(DATA, 3112) #https://epsg.io/?q=australia
DATA_sf <- st_as_sf(DATA, coords)


## map SA2 cancer counts ##
DATA_sf$Cancer_counts<- as.numeric(DATA_sf$sa2_count)

##SIR Calculation ###
DATA_sf$SIR<-as.numeric(DATA_sf$sa2_count/DATA_sf$sa2_expect)


## Define neighboring polygons (neighbor as being any contiguous polygon that shares at least one vertex)##-----------
nb <- poly2nb(DATA_sf, queen=FALSE)# nb lists all neighboring polygons

## Assigning weights to each neighboring polygon ##--------------------------
nbp_w <- nb2listw(nb, style="B", zero.policy=TRUE)


## Moran's I on observed count ##----------------------------------------
## Calculate Moran's I using functions ##-------------------------------- 
# Using moran.test function
moran.test(DATA_sf$sa2_count,nbp_w, zero.policy = TRUE)

## Using MC simulation method to test for significance ##--------------------
MC<- moran.mc(DATA_sf$sa2_count,nbp_w, nsim=999, zero.policy = TRUE)

# View results (including p-value)
MC

#####
## Fit model when data are aggregated at SA2 level ##--------------------------------
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
### To handle log0 situation as 0 cannot survive log transformation without being infinite ##----

DATA_sf$observed <- round(DATA_sf$sa2_count)
DATA_sf$expected<- DATA_sf$sa2_expect
DATA_sf$expected[DATA_sf$expected==0] <- 0.001 


## Fit BYM model ##----------------------------------------------------------
start <- Sys.time()
formula <-observed ~ 1+ offset(log(expected))
model_bym_wo_cov <- S.CARbym(formula=formula, data=DATA_sf, family="poisson", W=nbp_w,
                  burnin=500000, n.sample=1500000, thin=100)#10000 number of iterations

Sys.time() - start
save(model_bym_wo_cov, file="BYM_SA2_noCov.RData")

## Summary of results (summary output of the fitted model) ##---------
print(model_bym_wo_cov)
summary(model.bym_wo_cov)


## The entire posterior distribution can be viewed ##------------
## Traceplots for parameters ##
plot(model.bym_wo_cov$samples$tau2)
plot(model.bym_wo_cov$samples$sigma2)
plot(model.bym_wo_cov$samples$beta)
plot(model.bym_wo_cov$residuals)
plot(model.bym_wo_cov$fitted.values)

# Combine residuals and fitted values into a data frame
result = data.frame(fitted_values =model_bym_wo_cov$fitted.values, residuals = model_bym_wo_cov$residuals)
DATA_sf$fitted_SIR_wo_cov<-result$fitted_values/DATA_sf$expected

# Plot observed versus fitted values
Fit_obs <- data.frame(observed=observed, fitted=model.bym_wo_cov$fitted.values)

## Moran's I to the residuals (compute Moran's I statistic and test to assess its significance)##--------------------------------------------
res.model <- residuals(model.bym_wo_cov)
nb <- poly2nb(DATA_sf, queen=FALSE)
nbp_w <- nb2listw(nb, style="B", zero.policy=TRUE)
moran.mc(x = res.model,listw = nbp_w, nsim = 10000, zero.policy = TRUE)
