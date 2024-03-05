##################################################################################################################################################
##################################################################################################################################################
## This code belongs to Victor Aguiar de Souza Penha
## Entittled study: Sex, molt, and brood patch: drivers of body condition variation in Atlantic Forest tanagers (Passeriformes: Thraupidae) 
## Scientific journal submitted: Emu Australis
## Year 2023

##################################################################################################################################################
##################################################################################################################################################
## Remove unwanted objects
rm(list=ls())

##################################################################################################################################################
##################################################################################################################################################
## Needed packages
library(phyr)
library(ape)
library(picante)
library(phytools)
library(BBmisc)
library(phyr)
library(devtools)
library(MASS)
library(picante)
library(lme4)
library(phyloch)
library(car)
library(sf)
library(rgdal)
library(raster)  
library(tidyverse)
library(regclass)
library(lmerTest)
library(phangorn)
library(arm)
library(caper)
library(broom.mixed)
library(boot)
library(forestplot)

##################################################################################################################################################
##################################################################################################################################################
## Data: you can find the data here: https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecy.2647
data <- read.csv("./ATLANTIC_BIRD_TRAITS_completed_2018_11_d05.csv")

## Pre-processing the data
## Only tanagers
data <- data[which(data$Order == "Passeriformes"),]
data <- data[which(data$Family == "Thraupidae"),]

## Mass
any(is.na(data$Body_mass.g.))
data <- data[-which(is.na(data$Body_mass.g.)),]
class(data$Body_mass.g.)

## Length
any(is.na(data$Tail_length.mm.))
data <- data[-which(is.na(data$Tail_length.mm.)),]
class(data$Tail_length.mm.)

## Age
table(data$Age)
## Remove unknown
data <- data[-which(data$Age == "Unknown"),]
data <- data[-which(is.na(data$Age)),]
class(data$Age)
data$Age <- as.factor(data$Age)

## Sex
table(data$Sex)
data <- data[-which(data$Sex == "Unknown"),]
any(is.na(data$Sex))
data <- data[-which(is.na(data$Sex)),]
class(data$Sex)
data$Sex <- as.factor(data$Sex)

## Flight Molt
table(data$Flight_molt)
data <- data[-which(is.na(data$Flight_molt)),]
data$Flight_molt <- as.factor(data$Flight_molt)
levels(data$Flight_molt) <- c("Absent", "Present")

## Body molt
table(data$Body_molt)
data <- data[-which(is.na(data$Body_molt)),]
data$Body_molt <- as.factor(data$Body_molt)
levels(data$Body_molt) <- c("Absent", "Present")

## Reproductive stage
table(data$Reproductive_stage)
data <- data[-which(is.na(data$Reproductive_stage)),]
data$Reproductive_stage <- as.factor(data$Reproductive_stage)
levels(data$Reproductive_stage) <- c("Absent", "Present")

##################################################################################################################################################
##################################################################################################################################################
## Remove less than 10 captures
problems <- as.data.frame(sort(table(data$Binomial)))
problems <- problems[-which(problems$Freq < 10),]
data <- data[-which(!data$Binomial %in% problems$Var1),]

##################################################################################################################################################
##################################################################################################################################################
## State
table(data$State)
any(is.na(data$State))
table(data[which(data$State == "BA" | data$State == "CE" | data$State == "SE"),c("Binomial", "State")])
data <- data[-which(data$State == "BA" | data$State == "CE" | data$State == "SE"),]

## Year
table(data$Year)
table(data[which(data$Year == 1975 | data$Year == 1976 | data$Year == 1979 | data$Year == 1980 | data$Year == 1981 | data$Year == 1986 | data$Year == 1995 | data$Year == 2018),c("Binomial", "Year")])
data <- data[-which(data$Year == 1975 | data$Year == 1976 | data$Year == 1979 | 
                      data$Year == 1980 | data$Year == 1981 | data$Year == 1986 | 
                      data$Year == 1995 | data$Year == 2018),]
which(is.na(data$Year))

##################################################################################################################################################
##################################################################################################################################################
## Final check
any(is.na(data$Sex))
table(data$Sex)
any(is.na(data$Flight_molt))
table(data$Flight_molt)
any(is.na(data$Year))
table(data$Year)
table(data$State)
table(data$Age)
table(data$Body_molt)

##################################################################################################################################################
##################################################################################################################################################
## Breeding stage
data$Date2 <- data$Date
data$Date2 <- format(as.Date(data$Date2, format = "%m/%d/%Y"), "%m")
table(data$Date2)
data$BreedingSeason <- NA

## Creating column
data[which(data$Date2 == "10" | data$Date2 == "11" | data$Date2 == "12" | data$Date2 == "01" | 
             data$Date2 == "02" | data$Date2 == "03"),"BreedingSeason"] <- "Breeding"
data[which(!complete.cases(data$BreedingSeason)),"BreedingSeason"] <- "OutsideBreeding"
table(data$BreedingSeason)

##################################################################################################################################################
##################################################################################################################################################
## Estimating the body condition
## Body condition
## Por estado
data$BodyCond <- NA
especies <- unique(data$Binomial)

lisEsp <- c()
lisEsta <- c()
listSMI_mean <- c() 
listSMI_sd <- c() 
listSize <- c()
listBred <- c()
listYear <- c()
listSex <- c()
listAge <- c()

for(i in 1:length(especies)){
  filtro <- data[which(data$Binomial == especies[i]),]
  estados <- unique(filtro$State)
  
  for(j in 1:length(estados)){
    massa <- filtro[which(filtro$State == estados[j]),"Body_mass.g."] 
    comprimento <- filtro[which(filtro$State == estados[j]),"Tail_length.mm."] 
    breed <- filtro[which(filtro$State == estados[j]),"BreedingSeason"] 
    year <- filtro[which(filtro$State == estados[j]),"Year"]
    sex <- filtro[which(filtro$State == estados[j]),"Sex"]
    age <- filtro[which(filtro$State == estados[j]),"Age"]
    
    if(length(massa) < 10){
      data[which(data$State == estados[j] & data$Binomial == especies[i]),"BodyCond"] <- NA
    }
    else{
      if(length(unique(breed)) > 1 & length(unique(year)) > 1 & length(unique(sex)) > 1 & length(unique(age)) > 1){
        modelo <- lmer(massa~comprimento + (1|breed) + (1|year) + (1|sex) + (1|age))
        residuos <- residuals(modelo)
        data[which(data$State == estados[j] & data$Binomial == especies[i]),"BodyCond"] <- round(residuos,3)
      }
      else{
        if(length(unique(breed)) > 1 & length(unique(year)) > 1 & length(unique(sex)) > 1){
          modelo <- lmer(massa~comprimento + (1|breed) + (1|year) + (1|sex))
          residuos <- residuals(modelo)
          data[which(data$State == estados[j] & data$Binomial == especies[i]),"BodyCond"] <- round(residuos,3)
        }
        else{
          if(length(unique(breed)) > 1 & length(unique(year)) > 1 & length(unique(age)) > 1){
            modelo <- lmer(massa~comprimento + (1|breed) + (1|year) + (1|age))
            residuos <- residuals(modelo)
            data[which(data$State == estados[j] & data$Binomial == especies[i]),"BodyCond"] <- round(residuos,3)
          }
          else{
            if(length(unique(year)) > 1 & length(unique(age)) > 1 & length(unique(sex))){
              modelo <- lmer(massa~comprimento + (1|sex) + (1|year) + (1|age))
              residuos <- residuals(modelo)
              data[which(data$State == estados[j] & data$Binomial == especies[i]),"BodyCond"] <- round(residuos,3)
            }
            else{
              if(length(unique(breed)) > 1 & length(unique(sex)) > 1 & length(unique(age)) > 1){
                modelo <- lmer(massa~comprimento + (1|breed) + (1|sex) + (1|age))
                residuos <- residuals(modelo)
                data[which(data$State == estados[j] & data$Binomial == especies[i]),"BodyCond"] <- round(residuos,3)
              }
              else{
                if(length(unique(breed)) > 1 & length(unique(year)) > 1){
                  modelo <- lmer(massa~comprimento + (1|breed) + (1|year))
                  residuos <- residuals(modelo)
                  data[which(data$State == estados[j] & data$Binomial == especies[i]),"BodyCond"] <- round(residuos,3)
                }
                else{
                  if(length(unique(breed)) > 1 & length(unique(age)) > 1){
                    modelo <- lmer(massa~comprimento + (1|breed) + (1|age))
                    residuos <- residuals(modelo)
                    data[which(data$State == estados[j] & data$Binomial == especies[i]),"BodyCond"] <- round(residuos,3)
                  }
                  else{
                    if(length(unique(breed)) > 1 & length(unique(sex)) > 1){
                      modelo <- lmer(massa~comprimento + (1|breed) + (1|sex))
                      residuos <- residuals(modelo)
                      data[which(data$State == estados[j] & data$Binomial == especies[i]),"BodyCond"] <- round(residuos,3)
                    }
                    else{
                      if(length(unique(year)) > 1 & length(unique(age)) > 1){
                        modelo <- lmer(massa~comprimento + (1|year) + (1|age))
                        residuos <- residuals(modelo)
                        data[which(data$State == estados[j] & data$Binomial == especies[i]),"BodyCond"] <- round(residuos,3)
                      }
                      else{
                        if(length(unique(year)) > 1 & length(unique(sex)) > 1){
                          modelo <- lmer(massa~comprimento + (1|year) + (1|sex))
                          residuos <- residuals(modelo)
                          data[which(data$State == estados[j] & data$Binomial == especies[i]),"BodyCond"] <- round(residuos,3)
                        }
                        else{
                          if(length(unique(age)) > 1 & length(unique(sex)) > 1){
                            modelo <- lmer(massa~comprimento + (1|age) + (1|sex))
                            residuos <- residuals(modelo)
                            data[which(data$State == estados[j] & data$Binomial == especies[i]),"BodyCond"] <- round(residuos,3)
                          }
                          else{
                            if(length(unique(breed)) > 1){
                              modelo <- lmer(massa~comprimento + (1|breed))
                              residuos <- residuals(modelo)
                              data[which(data$State == estados[j] & data$Binomial == especies[i]),"BodyCond"] <- round(residuos,3)
                            }
                            else{
                              if(length(unique(year)) > 1){
                                modelo <- lmer(massa~comprimento + (1|year))
                                residuos <- residuals(modelo)
                                data[which(data$State == estados[j] & data$Binomial == especies[i]),"BodyCond"] <- round(residuos,3)
                              }
                              else{
                                if(length(unique(age)) > 1){
                                  modelo <- lmer(massa~comprimento + (1|age))
                                  residuos <- residuals(modelo)
                                  data[which(data$State == estados[j] & data$Binomial == especies[i]),"BodyCond"] <- round(residuos,3)
                                }
                                else{
                                  if(length(unique(sex)) > 1){
                                    modelo <- lmer(massa~comprimento + (1|sex))
                                    residuos <- residuals(modelo)
                                    data[which(data$State == estados[j] & data$Binomial == especies[i]),"BodyCond"] <- round(residuos,3)
                                  }
                                  else{
                                    modelo <- lm(massa~comprimento)
                                    residuos <- residuals(modelo)
                                    data[which(data$State == estados[j] & data$Binomial == especies[i]),"BodyCond"] <- round(residuos,3)
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
      
      lisEsp <- rbind(lisEsp, especies[i])
      lisEsta <- rbind(lisEsta, estados[j])
      listSize <- rbind(listSize, length(massa))
      listSMI_mean <- rbind(listSMI_mean, mean(residuos))
      listSMI_sd <- rbind(listSMI_sd, sd(residuos))
      listBred <- rbind(listBred, paste(as.data.frame(table(breed))$breed, ": ", as.data.frame(table(breed))$Freq, collapse = "; "))
      listYear <- rbind(listYear, paste(as.data.frame(table(year))$year, ": ", as.data.frame(table(year))$Freq, collapse = "; "))
      listAge <- rbind(listAge, paste(as.data.frame(table(age))$age, ": ", as.data.frame(table(age))$Freq, collapse = "; "))
      listSex <- rbind(listSex, paste(as.data.frame(table(sex))$sex, ": ", as.data.frame(table(sex))$Freq, collapse = "; "))
    }
  }
}

##################################################################################################################################################
##################################################################################################################################################
## Removing the ones that the condition was not satisfied
any(is.na(data$BodyCond))
data <- data[-which(!complete.cases(data$BodyCond)),]

##################################################################################################################################################
##################################################################################################################################################
## Phylogeny
## You can find the phylogeny here, either the subset or the complete phylogeny: https://birdtree.org/subsets/
phylogeny <- read.tree("./AllBirdsEricson1.tre")
phylogeny <- maxCladeCred(phylogeny, rooted=T)

## Making the names teh same
data$Binomial <- gsub(" ", "_", data$Binomial)
sort(phylogeny$tip.label)

##################################################################################################################################################
##################################################################################################################################################
## Change names - as in birdtree
unique(data[which(!data$Binomial %in% phylogeny$tip.label),"Binomial"])
data[which(data$Binomial == "Microspingus_cabanisi"),"Binomial"] <- "Poospiza_cabanisi"
data[which(data$Binomial == "Pipraeidea_bonariensis"),"Binomial"] <- "Thraupis_bonariensis"
data[which(data$Binomial == "Tangara_ornata"),"Binomial"] <- "Thraupis_ornata"

## Checking
all(phylogeny$tip.label %in% data$Binomial)
all(data$Binomial %in% phylogeny$tip.label)

## Trimming phylogeny
corte_phylo <- which(!phylogeny$tip.label %in% data$Binomial)
phylo <- drop.tip(phylogeny, corte_phylo)
all(phylo$tip.label %in% data$Binomial)
all(data$Binomial %in% phylo$tip.label)

## Response distribution
hist(normalize(data$BodyCond))

##################################################################################################################################################
##################################################################################################################################################
## Adult only
data <- data[which(data$Age == "Adult"),]
table(data$Age)

##################################################################################################################################################
##################################################################################################################################################
## Remove less than 10 captures
problems <- as.data.frame(sort(table(data$Binomial)))
problems <- problems[-which(problems$Freq < 10),]
data <- data[-which(!data$Binomial %in% problems$Var1),]
sort(table(data$Binomial))

##################################################################################################################################################
##################################################################################################################################################
## Year again, because we removed rows
table(data$Year)
table(data[which(data$Year == 1977 | data$Year == 1999),c("Binomial", "Year")])
data <- data[-which(data$Year == 1977),]


##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################
## Model
## Remove males with brood patch, justifying in the manuscript
data <- data[-which(data$Sex == "Male" & data$Reproductive_stage == "Present"),]

## Checking proportions again
table(data$Sex)
table(data$Flight_molt)
table(data$Year)
table(data$State)
table(data$Age)
table(data$Body_molt)

##################################################################################################################################################
##################################################################################################################################################
## Final numbers
manuscript <- as.data.frame(sort(table(data$Binomial)))
manuscript[order(as.character(manuscript$Var1)),]
sum(manuscript$Freq)
length(manuscript$Var1)


##################################################################################################################################################
##################################################################################################################################################
table(data$Reproductive_stage)
data <- data[-which(data$Year == 1999),]

##################################################################################################################################################
##################################################################################################################################################
row.names(data) <- NULL
## Multicolinearidade
model_multi_haplo <- lmer(normalize(BodyCond) ~   Reproductive_stage + Sex * Flight_molt + Body_molt * Flight_molt + 
                            Reproductive_stage * Body_molt + (1|Longitude_decimal_degrees) + 
                            (1|Latitude_decimal_degrees) + (1|BreedingSeason) + 
                            (1|Year) + (1|Binomial), 
                          data = data, na.action = na.fail)
vif(model_multi_haplo)
## Remove Flight_molt:Body_molt
model_multi_haplo <- lmer(normalize(BodyCond) ~   Reproductive_stage + Sex * Flight_molt + 
                            Reproductive_stage * Body_molt + (1|Longitude_decimal_degrees) + 
                            (1|Latitude_decimal_degrees) + (1|BreedingSeason) + 
                            (1|Year) + (1|Binomial), 
                          data = data, na.action = na.fail)
vif(model_multi_haplo)

## Remove Sex * Flight_molt
model_multi_haplo <- lmer(normalize(BodyCond) ~   Reproductive_stage + Sex + Flight_molt + 
                            Reproductive_stage * Body_molt + (1|Longitude_decimal_degrees) + 
                            (1|Latitude_decimal_degrees) + (1|BreedingSeason) + 
                            (1|Year) + (1|Binomial), 
                          data = data, na.action = na.fail)
vif(model_multi_haplo)

##################################################################################################################################################
##################################################################################################################################################
## PGLMM
## Full model
mod_1 <- pglmm(normalize(BodyCond) ~ Reproductive_stage + Sex + Flight_molt + Reproductive_stage * Body_molt + 
                 (1|Longitude_decimal_degrees) + (1|Latitude_decimal_degrees) + (1|BreedingSeason) + (1|Year) + (1|Species), data, cov_ranef = list(sp = phylo), tree = phylo)

## Null model
mod_1_null <- pglmm(normalize(BodyCond) ~ 1 + 
                      (1|Longitude_decimal_degrees) + (1|Latitude_decimal_degrees) + (1|BreedingSeason) + (1|Year) + (1|Species),  cov_ranef = list(sp = phylo), data, tree = phylo)

## AIC comparision
for(i in 1){
  if(mod_1$AIC < mod_1_null$AIC){
    print("Full model had a lower AIC.")
    print(paste0("Full: ", mod_1$AIC))
    print(paste0("Null: ", mod_1_null$AIC))
  }
  else{
    print("Null model had a lower AIC.")
    print(paste0("Null: ", mod_1_null$AIC))
    print(paste0("Full: ", mod_1$AIC))
  }
}
## Full model was better. Let's evaluate the model.

## Getting predicted values
data$Pred <- as.data.frame(predict(mod_1))[,1]

## Getting residuals
data$Residuals <- as.data.frame(residuals(mod_1))[,1]

## Model summary
summary(mod_1)

## Plotting results
## Plotting confidence intervals
# Getting confidence intervals through bootstrapping for the categoreis
n_bootstrap <- 1000
n_coefs <- length(fixef(mod_1)$Value)
bootstrap_coefs <- matrix(NA, nrow = n_bootstrap, ncol = n_coefs)

## This will take a couple of hours
for (i in 1:n_bootstrap) {
  indices <- sample(nrow(haplos_10_2), replace = TRUE)
  sampled_data <- haplos_10_2[indices, ]
  bootstrap_model <- pglmm(normalize(BodyCond) ~ Reproductive_stage + Sex + Flight_molt + 
                             Reproductive_stage * Body_molt + (1|Longitude_decimal_degrees) + 
                             (1|Latitude_decimal_degrees) + (1|BreedingSeason) + 
                             (1|Year), data = sampled_data, tree = phylo)
  bootstrap_coefs[i, ] <- fixef(bootstrap_model)$Value
}
conf_intervals <- apply(bootstrap_coefs, 2, quantile, c(0.025, 0.975))
print(conf_intervals)

## sRemoving the interaction
conf_intervals_2 <- conf_intervals[,-6]
bootstrap_coefs_2 <- bootstrap_coefs[,-6]

## Plotting results
colnames(conf_intervals_2) <- c("Intercept", "Presence_Brood_Patch_(present)", "Sex_(male)", "Flight_feather_molt_(present)", "Flight_feather_molt_(present)")

## Plotting results
# Create a data frame for forest plot
forest_data <- data.frame(
  Variables = colnames(bootstrap_coefs),
  Estimate = colMeans(bootstrap_coefs),
  CI.Lower = conf_intervals[1, ],
  CI.Upper = conf_intervals[2, ]
)

# Create the forest plot
forestplot(
  mean = forest_data$Estimate,
  lower = forest_data$CI.Lower,
  upper = forest_data$CI.Upper,
  labeltext = forest_data$Variables,
  txt_round = 2,
  is.summary = c(TRUE, rep(FALSE, n_coefs - 1)),
  xticks = c(-2, -1, 0, 1, 2),  # Adjust as needed
  clip = c(-3, 3),  # Adjust as needed
  zero = 0,  # Adjust as needed
  col = fpColors(box = "black", lines = "black", summary = "black")
)

##################################################################################################################################################
##################################################################################################################################################
## Manuscript information
## Molt-breeding overal: results section 1
dim(data[which(data$Reproductive_stage == "Present"),])
dim(data[which(data$Flight_molt == "Present" & data$Reproductive_stage == "Present"),])
dim(data[which(data$Body_molt == "Present" & data$Reproductive_stage == "Present"),])
dim(data[which(data$Body_molt == "Present" & data$Flight_molt == "Present" & data$Reproductive_stage == "Present"),])
table(data[which(data$Body_molt == "Present" & data$Reproductive_stage == "Present"),"Binomial"])
table(data[which(data$Flight_molt == "Present" & data$Reproductive_stage == "Present"),"Binomial"])
table(data[which(data$Body_molt == "Present" & data$Flight_molt == "Present" & data$Reproductive_stage == "Present"),"Binomial"])

##################################################################################################################################################
##################################################################################################################################################
## Mean + sd
round(mean(data$BodyCond),3)
sd(data$BodyCond)
data[which.max(data$BodyCond),]
data[which.min(data$BodyCond),]

##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################
### Makign maps
## https://docs.google.com/file/d/0B__Rg9h09RtfQU9jbEpxcW9xREk/edit?pli=1&resourcekey=0-po0pdYgc-QZ5Pt13twjIsg
my_spdf <- shapefile("./estados_2010.shp",verbose=FALSE)
data.pts <- SpatialPoints(data[,c("Longitude_decimal_degrees", "Latitude_decimal_degrees")])

## Map
lnd <- (SpatialPolygonsDataFrame(Sr = spTransform(my_spdf, CRSobj = CRS("+init=epsg:4326")), data = my_spdf@data))
lnd.f <- fortify(lnd)

## Plotting results
dat_first <- st_as_sf(
  data.pts,
  coords = c("Longitude_decimal_degrees", "Latitude_decimal_degrees"),
  crs = 4326
)
jpeg("./Map.jpeg", width = 8, height = 8, units = 'in', res = 800) 
dat_first %>% 
  ggplot()+
  geom_sf(col = "black", size = 2)+
  xlab("Longitude")+
  ylab("Latitude") + 
  geom_polygon(data = lnd.f, aes(x = long, y = lat, group = group), fill=NA) +
  geom_path(data = lnd.f, aes(x = long, y = lat, group = group), color = "black") +
  theme_bw() + ggtitle(expression(italic("Sampling locations")))
dev.off()

##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################
## SUpplementary material

##################################################################################################################################################
##################################################################################################################################################
## Supplementary table 1
## Paper tables
especies <- unique(data$Binomial)
listespecie <- c()
listestado <- c()
listaanos <- c()
list <- c()
list2 <- c()
list3 <- c()
list4 <- c()
list5 <- c()
list6 <- c()
list7 <- c()
list8 <- c()
list9 <- c()
list10 <- c()
list11 <- c()
list12 <- c()

for(k in 1:length(especies)){
  filtro <- data[which(data$Binomial == especies[k]),]
  estados <- sort(unique(filtro$State))
  
  for(i in 1:length(estados)){
    filtro2 <- filtro[which(filtro$State == estados[i]),]
    anos <- sort(unique(filtro2$Year))
    
    for(j in 1:length(anos)){
      filtro3 <- filtro2[which(filtro2$Year == anos[j]),]
      
      listespecie <- rbind(listespecie, especies[k])
      listestado <- rbind(listestado, estados[i])
      listaanos <- rbind(listaanos, anos[j])
      list <- rbind(list, dim(filtro3[which(filtro3$Sex == "Male" & filtro3$Age == "Adult" & filtro3$Flight_molt == "Present"),])[1])
      list2 <- rbind(list2, dim(filtro3[which(filtro3$Sex == "Male" & filtro3$Age == "Adult" & filtro3$Flight_molt == "Absent"),])[1])
      list3 <- rbind(list3, dim(filtro3[which(filtro3$Sex == "Male" & filtro3$Age == "Adult" & filtro3$Body_molt == "Present"),])[1])
      list4 <- rbind(list4, dim(filtro3[which(filtro3$Sex == "Male" & filtro3$Age == "Adult" & filtro3$Body_molt == "Absent"),])[1])
      list5 <- rbind(list5, dim(filtro3[which(filtro3$Sex == "Male" & filtro3$Age == "Adult" & filtro3$BreedingSeason == "Breeding"),])[1])
      list6 <- rbind(list6, dim(filtro3[which(filtro3$Sex == "Male" & filtro3$Age == "Adult" & filtro3$BreedingSeason == "OutsideBreeding"),])[1])
      list7 <- rbind(list7, dim(filtro3[which(filtro3$Sex == "Female" & filtro3$Age == "Adult" & filtro3$Flight_molt == "Present"),])[1])
      list8 <- rbind(list8, dim(filtro3[which(filtro3$Sex == "Female" & filtro3$Age == "Adult" & filtro3$Flight_molt == "Absent"),])[1])
      list9 <- rbind(list9, dim(filtro3[which(filtro3$Sex == "Female" & filtro3$Age == "Adult" & filtro3$Body_molt == "Present"),])[1])
      list10 <- rbind(list10, dim(filtro3[which(filtro3$Sex == "Female" & filtro3$Age == "Adult" & filtro3$Body_molt == "Absent"),])[1])
      list11 <- rbind(list11, dim(filtro3[which(filtro3$Sex == "Female" & filtro3$Age == "Adult" & filtro3$BreedingSeason == "Breeding"),])[1])
      list12 <- rbind(list12, dim(filtro3[which(filtro3$Sex == "Female" & filtro3$Age == "Adult" & filtro3$BreedingSeason == "OutsideBreeding"),])[1])
    }
  }
}
listespecie <- as.data.frame(listespecie)
listestado <- as.data.frame(listestado)
listaanos <- as.data.frame(listaanos)
list <- as.data.frame(list)
list2 <- as.data.frame(list2)
list3 <- as.data.frame(list3)
list4 <- as.data.frame(list4)
list5 <- as.data.frame(list5)
list6 <- as.data.frame(list6)
list7 <- as.data.frame(list7)
list8 <- as.data.frame(list8)
list9 <- as.data.frame(list9)
list10 <- as.data.frame(list10)
list11 <- as.data.frame(list11)
list12 <- as.data.frame(list12)

final <- cbind(listespecie, listestado, listaanos, list, list2, list3, list4, list5, list6, list7, list8, list9, list10, list11, list12)

##################################################################################################################################################
##################################################################################################################################################
## Supplementary table 3
lisEsp <- as.data.frame(lisEsp)
lisEsta <- as.data.frame(lisEsta)
listSMI_mean <- as.data.frame(listSMI_mean)
listSMI_sd <- as.data.frame(listSMI_sd)
listSize <- as.data.frame(listSize)
listBred <- as.data.frame(listBred)
listYear <- as.data.frame(listYear)
listAge <- as.data.frame(listAge)
listSex <- as.data.frame(listSex)

cbind(lisEsp, lisEsta, listSMI_mean, listSMI_sd, listSize, listBred, listYear, listAge, listSex)

##################################################################################################################################################
##################################################################################################################################################
## Get final data here - contains the unique locations
final_data <- read.csv("./Dados_Finais.csv")
locations <- unique(final_data$Site)
matriz <- as.data.frame(matrix(ncol = 3, nrow = length(locations), NA))
colnames(matriz) <- c("Site", "Species", "Quantity")
matriz$Site <- locations
for(i in 1:length(locations)){
  filtro <- as.data.frame(table(final_data[which(final_data$Site == locations[i]),"Binomial"]))
  matriz$Species[i] <- paste(filtro$Var1, collapse = " ")
  matriz$Quantity[i] <- paste(filtro$Freq, collapse = " ")
}

##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################
##### Figures
## Phylogeny 
plot(phylo)

## Model fit
binnedplot(data$Pred,data$Residuals)

## End of code
##################################################################################################################################################
##################################################################################################################################################
