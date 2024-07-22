### Fitting BYM model without covariate to SA1 level data in QLD
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


## Load Mesh & SA1 data ## -----------------------------------------------------------
shape_loc <- "Data/Shapefiles_Mesh_SA1_SA2_SA3_SA4/SA1_2016_AUST.shp"
SA1_raw <- st_read(shape_loc)
SA1_raw_qld <- SA1_raw %>% #Taking records for Queensland only
  filter(STE_NAME16 == "Queensland")%>%
  mutate(SA1_MAIN16 = as.numeric(as.character(SA1_MAIN16)))


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

# ## Check the number of SA1_MAIN16 is in this file
# ## It is less than SA1_raw_qld which means there are some SA1 where there are no population data
# length(unique(mesh[["SA1_MAIN16"]]))

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
aca<- read_csv("Data/Simulated lung cancer data for QLD.csv") # not available in the repository, need to be requested to CCQ at: Statistics@cancerqld.org.au
conc <- read_csv("Data/sa2_2011_to_sa2_2016.csv")

##collapse the sex variables & group by using only sa2_name ##-------------------------
aca_total <- aca %>% 
  dplyr::select(sa2, sa2_name, sex, count, expect, SA2irsd) %>% 
  group_by(sa2_name,SA2irsd) %>% 
  summarise(sa2_count = sum(count),
            sa2_expect = sum(expect))

aca_cov <- aca_total

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

## Load SA1_IRSD data ##
SA1_IRSD <- read_csv("Data/SA1_IRSD.csv")
## Join mesh_raw with SA1_IRSD
SA1_IRSD_full <- SA1_IRSD %>% 
  inner_join(.,mesh_raw, by = c("SA1_MAIN16"="SA1_MAIN16")) 

## Calculate SA1 population ##
sa1_pop <- SA1_IRSD_full %>% 
  group_by(SA1_MAIN16) %>% 
  summarise(sa1_pop = sum(SA1_pop))

SA1_IRSD_full <- SA1_IRSD_full %>%
  inner_join(.,sa1_pop,by = c("SA1_MAIN16"= "SA1_MAIN16"))

SA1_IRSD_calc <- SA1_IRSD_full %>%
  mutate(IRSD_SA1=(SA1_irsdscore*SA1_pop)/sa1_pop)

SA1_IRSD <- SA1_IRSD_calc %>%
  group_by(SA1_MAIN16) %>%
  summarise(SA1_irsdscore=round(sum(IRSD_SA1)))

## Calculate quintile ##----------------------------------------------
SA1_irsd = quantile(SA1_IRSD$SA1_irsdscore, probs = seq(0, 1, .2))
SA1_IRSD$SA1_irsd <- SA1_IRSD$SA1_irsdscore
SA1_IRSD$SA1_irsd[SA1_IRSD$SA1_irsdscore<=SA1_irsd[2]] <- 1
SA1_IRSD$SA1_irsd[SA1_IRSD$SA1_irsdscore>SA1_irsd[2] & SA1_IRSD$SA1_irsdscore<=SA1_irsd[3]] <- 2 
SA1_IRSD$SA1_irsd[SA1_IRSD$SA1_irsdscore>SA1_irsd[3] & SA1_IRSD$SA1_irsdscore<=SA1_irsd[4]] <- 3
SA1_IRSD$SA1_irsd[SA1_IRSD$SA1_irsdscore>SA1_irsd[4] & SA1_IRSD$SA1_irsdscore<=SA1_irsd[5]] <- 4
SA1_IRSD$SA1_irsd[SA1_IRSD$SA1_irsdscore>SA1_irsd[5]] <- 5 

# ## Since not all SA1 has IRSD Score, this inner join reduced number of areas
# # ## Join SA1_IRSD with SA1_level_count and SA1_raw_qld ##
 DATA_SA1 <- SA1_raw_qld %>% #
   inner_join(., SA1_IRSD, by = c("SA1_MAIN16"="SA1_MAIN16")) %>%
   inner_join(.,SA1_level_count, by = c ("SA1_MAIN16"="SA1_MAIN16"))

# DATA<- DATA_SA1 %>% 
#    filter(is.na(sa1_count) == FALSE) %>%
#    filter(is.na(sa1_expect) == FALSE) 

## Exclude records with missing values ##--------------------------------------------
DATA<- DATA_SA1%>%
   filter(is.na(sa1_count) == FALSE)%>%
   filter(is.na(sa1_expect) == FALSE)%>%
   dplyr::select(SA1_MAIN16, sa1_pop, sa1_prop, sa1_count,sa1_expect,  SA1_irsdscore, SA1_irsd, geometry )
 
 
## Creating spatial object ##------------------------------------------
coords <- st_transform(DATA, 3112) #https://epsg.io/?q=australia
DATA_sf <- st_as_sf(DATA, coords)

## map SA1 cancer counts ##
DATA_sf$Cancer_counts<- as.numeric(DATA_sf$sa1_count)

##SIR Calculation ###
DATA_sf$SIR<-as.numeric(DATA_sf$sa1_count/DATA_sf$sa1_expect)

## Descriptive statistics ##
summary(DATA_sf$sa1_count)
sd(DATA_sf$sa1_count)

# ### To handle log0 situation as 0 cannot survive log transformation without being infinite ##----
DATA_sf$observed <- round(DATA_sf$sa1_count)
DATA_sf$expected<-DATA_sf$sa1_expect
DATA_sf$expected[DATA_sf$expected==0] <- 0.001 

### Creating factor variables ###-----------------------------------------------
DATA_sf$SA1_irsd.f <- factor(DATA_sf$SA1_irsd)
is.factor(DATA_sf$SA1_irsd.f)
levels(DATA_sf$SA1_irsd.f)

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

DATA_sf_model <- DATA_sf %>%
  dplyr::select(SA1_MAIN16,observed,expected,SA1_irsd.f, geometry)

###### Fit BYM model with covariate #######
start <- Sys.time()
formula <-observed~offset(log(expected))+SA1_irsd.f

model_bym<- S.CARbym(formula= formula, data=DATA_sf_model, family="poisson", W=nbp_w,
                     burnin=500000, n.sample=1500000, thin=100)#10000 number of iterations

Sys.time() - start
save(model_bym, file="BYM_SA1_cov.RDtata")
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
