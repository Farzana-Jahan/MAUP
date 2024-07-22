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

## Collapse the sex variables & group by using only sa2_name ##-------------------------
aca_total <- aca %>% 
  dplyr::select(sa2, sa2_name,sex, count, expect,irsdscore) %>% 
  group_by(sa2_name,irsdscore) %>% 
  summarise(sa2_count = sum(count),
            sa2_expect = sum(expect))

aca_cov <- aca_total

## Checking for duplicates ##---------------------------------------------------
length(unique(aca_total$sa2_name))

## use 2016 SA2 codes instead of 2011 ##-----------------------------------------
aca_SA2_16 <- aca_total %>% 
  inner_join(.,conc, by = c("sa2_name"="SA2_NAME_2011")) %>% 
  dplyr::select(sa2_count, sa2_expect, SA2_MAINCODE_2016, 
                SA2_NAME_2016,sa2_name, RATIO) %>%
  mutate(sa2_count=sa2_count*RATIO,
         sa2_expect=sa2_expect*RATIO) %>%
  group_by(SA2_MAINCODE_2016, SA2_NAME_2016) %>%
  summarise(sa2_count=sum(sa2_count),
            sa2_expect=sum(sa2_expect))

aca_with_cov <- aca_SA2_16 %>%
  inner_join(.,aca_cov,by = c("SA2_NAME_2016"="sa2_name")) %>%
  dplyr::select(SA2_MAINCODE_2016,SA2_NAME_2016,sa2_count=sa2_count.x,
                sa2_expect=sa2_expect.x,irsdscore)

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

## map SA2 cancer counts ##
DATA_sf$Cancer_counts<- as.numeric(DATA_sf$sa2_count)

##SIR Calculation ###
DATA_sf$SIR<-as.numeric(DATA_sf$sa2_count/DATA_sf$sa2_expect)

## Fit model when data are aggregated at SA2 level ##--------------------------------

### To handle log0 situation as 0 cannot survive log transformation without being infinite ##----
DATA_sf$observed <- round(DATA_sf$sa2_count)
DATA_sf$expected <- DATA_sf$sa2_expect
DATA_sf$expected[DATA_sf$expected==0] <- 0.001 

DATA_sf$SA2irsd.f <- factor(DATA_sf$SA2_irsd)
is.factor(DATA_sf$SA2irsd.f)
levels(DATA_sf$SA2irsd.f)

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
  dplyr::select(SA2_MAIN16,observed,expected,SA2irsd.f, geometry)
  

###### Fit BYM model with covariate #######
start <- Sys.time()
formula <-observed~offset(log(expected))+SA2irsd.f

model_bym<- S.CARbym(formula= formula, data=DATA_sf_model, family="poisson", W=nbp_w,
                   burnin=500000, n.sample=1500000, thin=100)#10000 number of iterations

Sys.time() - start
save(model_bym, file="BYM_with_cov.RData")
## Summary of results (summary output of the fitted model) ##---------
print(model_bym)
summary(model_bym)

## Check samples from the model ##----------------------------------------------
summary(model_bym$samples)

###For bym model####
plot(model_bym$samples$tau2)
plot(model_bym$samples$sigma2)
plot(model_bym$samples$beta)


# Combine residuals and fitted values into a data frame
result = data.frame(fitted_values = model_bym$fitted.values, residuals = model_bym$residuals)
fitted_SIR<-result$fitted_values/DATA_sf_model$expected

# Plot observed versus fitted values
Fit_obs <- data.frame(observed=DATA_sf_model$observed, fitted=model_bym$fitted.values)

## Moran's I to the residuals (compute Moran's I statistic and test to assess its significance)##--------------------------------------------
## Define neighboring polygons (neighbor as being any contiguous polygon that shares at least one vertex)##-----------
nb <- poly2nb(DATA_sf, queen=FALSE)# nb lists all neighboring polygons
## Assigning weights to each neighboring polygon ##--------------------------
nbp_w <- nb2listw(nb, style="B", zero.policy=TRUE)
res.model <- residuals(model_bym)
moran.mc(x = res.model,listw = nbp_w, nsim = 10000, zero.policy = TRUE)

