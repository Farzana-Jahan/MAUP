### Fitting BYM model without covariate to Mesh block levels data in QLD
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

## Load mesh data ## -----------------------------------------------------------
shape_loc <- "Data/Shapefiles_Mesh_SA1_SA2_SA3_SA4/MB_2016_QLD.shp"
mesh_raw <- st_read(shape_loc) %>% 
  mutate(MB_CODE16 = as.numeric(as.character(MB_CODE16)),
         SA2_MAIN16 = as.numeric(as.character(SA2_MAIN16)),
         SA2_5DIG16 = as.numeric(as.character(SA2_5DIG16)))
mesh_pop <- read_csv("Data/2016 census mesh block counts.csv") %>% 
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

## Load simulated cancer data and concordance file ## ----------------------------------------------
aca<- read_csv("Data/Simulated lung cancer data for QLD.csv") # not available in the repository, need to be requested to CCQ at: Statistics@cancerqld.org.au
conc <- read_csv("Data/sa2_2011_to_sa2_2016.csv")

##collapse the sex variables & group by using only sa2_name ##-------------------------
aca_total <- aca %>% 
  dplyr::select(sa2, sa2_name, sex,count, expect) %>% 
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


## Load SA1_IRSD data ##--------------------------------------------------------
SA1_IRSD <- read_csv("Data/SA1_IRSD.csv")
## Join mesh_raw with SA1_IRSD
SA1_IRSD_full <- SA1_IRSD %>% 
  mutate(SA1_MAIN16 = as.character(SA1_MAIN16)) %>%
  inner_join(.,mesh, by = c("SA1_MAIN16"="SA1_MAIN16")) 

## Calculate mesh level IRSD ##-------------------------------------------------
mesh_IRSD_calc <- SA1_IRSD_full %>%
  mutate(mesh_irsdscore=(SA1_irsdscore*mesh_population)/SA1_pop)


## Calculate quintile ##----------------------------------------------
mesh_irsd_q = quantile(mesh_IRSD_calc$mesh_irsdscore, probs = seq(0, 1, .2))

mesh_IRSD_calc$mesh_irsd <- mesh_IRSD_calc$mesh_irsdscore
mesh_IRSD_calc$mesh_irsd[mesh_IRSD_calc$mesh_irsdscore<=mesh_irsd_q[2]] <- 1
mesh_IRSD_calc$mesh_irsd[mesh_IRSD_calc$mesh_irsdscore>mesh_irsd_q[2] & mesh_IRSD_calc$mesh_irsdscore<=mesh_irsd_q[3]] <- 2 
mesh_IRSD_calc$mesh_irsd[mesh_IRSD_calc$mesh_irsdscore>mesh_irsd_q[3] & mesh_IRSD_calc$mesh_irsdscore<=mesh_irsd_q[4]] <- 3
mesh_IRSD_calc$mesh_irsd[mesh_IRSD_calc$mesh_irsdscore>mesh_irsd_q[4] & mesh_IRSD_calc$mesh_irsdscore<=mesh_irsd_q[5]] <- 4
mesh_IRSD_calc$mesh_irsd[mesh_IRSD_calc$mesh_irsdscore>mesh_irsd_q[5]] <- 5 


## Drop Geometry
mesh_IRSD_calc <- st_sf(mesh_IRSD_calc) %>%
  st_drop_geometry() %>%
  dplyr::select(MB_CODE16, mesh_irsd)


## Join SA1_IRSD with SA1_level_count and SA1_raw_qld ##
DATA_mesh <- mesh_IRSD_calc %>% 
  inner_join(., DATA, by = c("MB_CODE16"="MB_CODE16")) 


## Exclude missing values ##--------------------------------------------
DATA<- DATA_mesh %>% 
  filter(is.na(mesh_count) == FALSE)%>%
  filter(is.na(mesh_expect) == FALSE) 


## Descriptive statistics ##
summary(DATA$mesh_count)
sd(DATA$mesh_count)

DATA_sf<-st_sf(DATA)

## Creating spatial object ##------------------------------------------
coords <- st_transform(DATA_sf, 3112) #https://epsg.io/?q=australia
DATA_sf <- st_as_sf(DATA_sf, coords)


## Map mesh cancer counts ##
DATA_sf$Cancer_counts<- as.numeric(DATA_sf$mesh_count)

##SIR Calculation ###
DATA_sf$SIR<-as.numeric(DATA_sf$mesh_count/DATA_sf$mesh_expect)


#####
## Fit model ##--------------------------------
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
## For 0 count in some meshblocks and integer value, we had to round it and add 1 
# to the measurement in the following calculaton: 
DATA_sf$observed <- round(DATA_sf$mesh_count)
DATA_sf$expected<-DATA_sf$mesh_expect
DATA_sf$expected[DATA_sf$expected==0] <- 0.001

### Creating factor variables ###-----------------------------------------------
DATA_sf$mesh_irsd.f <- factor(DATA_sf$mesh_irsd)
is.factor(DATA_sf$mesh_irsd.f)
levels(DATA_sf$mesh_irsd.f)

DATA_sf_model <- DATA_sf %>%
  dplyr::select(MB_CODE16,observed,expected,mesh_irsd.f, geometry)

## Fit BYM model ##----------------------------------------------------------
start <- Sys.time()
formula <-observed~offset(log(expected))+mesh_irsd.f

model_bym<- S.CARbym(formula= formula, data=DATA_sf_model, family="poisson", W=nbp_w,
                     burnin=500000, n.sample=1500000, thin=100)#10000 number of iterations

Sys.time() - start
save(model_bym, file= "BYM_mesh_cov.RData")
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

# Plot observed versus fitted values
Fit_obs <- data.frame(observed=DATA_sf_model$observed, fitted=model_bym$fitted.values)

## Moran's I to the residuals (compute Moran's I statistic and test to assess its significance)##--------------------------------------------
## Define neighboring polygons (neighbor as being any contiguous polygon that shares at least one vertex)##-----------
nb <- poly2nb(DATA_sf_model, queen=FALSE)# nb lists all neighboring polygons
## Assigning weights to each neighboring polygon ##--------------------------
nbp_w <- nb2listw(nb, style="B", zero.policy=TRUE)
res.model <- residuals(model_bym)
moran.mc(x = res.model,listw = nbp_w, nsim = 10000, zero.policy = TRUE)

