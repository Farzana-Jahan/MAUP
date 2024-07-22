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


## Load Mesh & SA1 data ## -----------------------------------------------------------
shape_loc <- "Data/Shapefiles_Mesh_SA1_SA2_SA3_SA4/SA1_2016_AUST.shp"
SA1_raw <- st_read(shape_loc)
SA1_raw_qld <- SA1_raw %>% #Taking records for Queensland only
  filter(STE_CODE16 == "3")%>%
  mutate(STATE_CODE = as.numeric(as.character(STE_CODE16)))

##Load mesh data
shape_loc <- "Data/Shapefiles_Mesh_SA1_SA2_SA3_SA4/MB_2016_QLD.shp"
mesh_raw <- st_read(shape_loc) %>% 
  st_drop_geometry() %>% 
  mutate(SA2_MAIN16 = as.numeric(as.character(SA2_MAIN16)),
         SA3_CODE16 = as.numeric (as.character(SA3_CODE16)),
         SA4_CODE16 = as.numeric (as.character(SA4_CODE16)),
         SA1_MAIN16= as.numeric(as.character(SA1_MAIN16)),
         MB_CODE16= as.numeric(as.character(MB_CODE16)))

## Load simulated cancer data ## ----------------------------------------------
aca<- read_csv("Data/Simulated lung cancer data for QLD.csv") # not available in the repository, need to be requested to CCQ at: Statistics@cancerqld.org.au
conc <- read_csv("Data/sa2_2011_to_sa2_2016.csv")
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

## Load simulated cancer data and concordance file ## ----------------------------------------------
aca<- read_csv("Data/Simulated lung cancer data for QLD.csv")
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


##Join SA1_level with aca_SA2_16 and calculate sa1_count, sa1_expect using sa1_prop
SA1_level_count <- SA1_level %>%
  inner_join(.,aca_SA2_16, by = c("SA2_MAIN16"="SA2_MAINCODE_2016")) %>%
  mutate(sa1_count = sa2_count*sa1_prop,
        sa1_expect = sa2_expect*sa1_prop)
class(SA1_raw_qld$SA1_MAIN16)
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

##SIR Calculation ###
DATA_sf$SIR<-as.numeric(DATA_sf$sa1_count/DATA_sf$sa1_expect)

# ## Define neighboring polygons (neighbor as being any contiguous polygon that shares at least one vertex)##-----------
nb <- poly2nb(DATA_sf, queen=FALSE)# nb lists all neighboring polygons
 
# ## Assigning weights to each neighboring polygon ##--------------------------
nbp_w <- nb2listw(nb, style="B", zero.policy=TRUE)
# 
# ## Moran's I on observed count ##----------------------------------------
# ## Calculate Moran's I using functions ##-------------------------------- 
# # Using moran.test function
moran.test(DATA_sf$sa1_count,nbp_w, zero.policy = TRUE)
# ## Using MC simulation method to test for significance ##--------------------
MC<- moran.mc(DATA_sf$sa1_count,nbp_w, nsim=999, zero.policy = TRUE)
# # View results (including p-value)
MC

#####
## Fit model when data are aggregated at SA1 level ##--------------------------------
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
DATA_sf$observed <- round(DATA_sf$sa1_count)
DATA_sf$expected <- DATA_sf$sa1_expect
DATA_sf$expected[DATA_sf$expected==0] <- 0.001 

## Fit BYM model ##----------------------------------------------------------
start <- Sys.time()
formula <-observed ~ 1+ offset(log(expected))
model.BYM_wocov<- S.CARbym(formula=formula, data=DATA_sf, family="poisson", W=nbp_w,
                      burnin=500000, n.sample=1500000, thin=100)#10000 number of iterations

Sys.time() - start
save(model.BYM_wocov,file= "BYM_SA1_noCov.RData")

## Summary of results (summary output of the fitted model) ##---------
print(model.BYM_wocov)
summary(model.BYM_wocov)


## Traceplots for parameters ##
plot(model.BYM_wocov$samples$tau2)
plot(model.BYM_wocov$samples$sigma2)
plot(model.BYM_wocov$samples$beta)
plot(model.BYM_wocov$residuals)
plot(model.BYM_wocov$fitted.values)


# Combine residuals and fitted values into a data frame
result = data.frame(fitted_values = model.BYM_wocov$fitted.values, 
                    residuals = model.BYM_wocov$residuals)


# Plot observed versus fitted values
Fit_obs <- data.frame(observed=observed, fitted=model.BYM_wocov$fitted.values)

## Moran's I to the residuals (compute Moran's I statistic and test to assess its significance)##--------------------------------------------
res.model <- residuals(model.BYM_wocov)
nb <- poly2nb(DATA_sf, queen=FALSE)# nb lists all neighboring polygons
nbp_w <- nb2listw(nb, style="B", zero.policy=TRUE)
moran.mc(x = res.model,listw = nbp_w, nsim = 10000)


