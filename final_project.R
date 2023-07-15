# Final Project of R discipline

## Chapter I - Introduction to Agent-Based Model (INCOMPLETE)

#In this phase, I am carrying packages and data to use in next steps. 

install.packages("Rcpp")
install.packages("raster")
install.packages("plyr")
install.packages("stringr")
install.packages("wallace")
install.packages("spocc")
install.packages("spThin")
install.packages("dismo")
install.packages("rgeos")
install.packages("ENMeval")
install.packages("dplyr")

library(Rcpp)
library(raster)
library(plyr)
library(stringr)
library(wallace)
library(spocc)
library(spThin)
library(dismo)
library(rgeos)
library(ENMeval)
library(dplyr)

mydata <- read.csv("data/Database.csv", header = TRUE, sep = ",")

## Chapter II - Thermal Plasticity Database Exploration 

#Here, I am going to present a plot and analysis to test if is there relationship between temperature and plasticity. 

#mydata_2 <- na.omit(my.data)
Plasticity <- as.numeric(mydata$mean)
Temperature <- as.numeric(mydata$T)
plot(Plasticity ~ Temperature)
model <- lm(Plasticity ~ Temperature)
summary(model)
Plasticity <- ceiling(Plasticity)

### New plot removing a little values from graph ###
n <- length(Plasticity)
novo <- numeric(n)
for (i in 1:n) {
  ifelse(Plasticity[i] < 1, novo[i] <- NA, novo[i] <- Plasticity[i])
}  

plot(novo ~ Temperature)


#In this phase of database exploration, I am going to do a function to select more frequent data or species from my database.  

### Identifying more frequent species ####

frequent <- function(column) {
  unique(column)
  caracterizando <- table(column)
  caracterizando2 <- as.data.frame(caracterizando)
  ordenado <- sort(caracterizando2$Freq, decreasing = TRUE)
  mais_frequente <- ordenado[1]
  n <- length(caracterizando2$Freq)
  novo_de_novo <- numeric(n)
  for (i in 1:n) {
    if (caracterizando2$Freq[i] == mais_frequente) {
      novo_de_novo[i] <- TRUE
    }
  }
  
  posicao <- which(grepl(1, novo_de_novo))
  names(caracterizando2) <- c("Var1", "Freq")
  resultado <- caracterizando2$Var1[posicao]
  resultado2 <- as.character(resultado)
  return(resultado2)
}
frequent(mydata$species)

especie_nome_final <- frequent(mydata$species)

posicao_sepentina <- which(grepl("serpentina", mydata$species))
posicao_sepentina2 <- posicao_sepentina[1]
genero <- mydata$genus[posicao_sepentina2]

tmp <- cbind(genero, especie_nome_final)

specie_more_frequent <- str_c(tmp, collapse = " ")
specie_more_frequent


#Below, I am going to present a niche model code withdraw of model done utilizing Wallace package. I did it to my more frequent species in previous database.  

source(system.file('shiny/funcs', 'functions.R', package = 'wallace'))

results <- spocc::occ(query = "Chelydra serpentina", from = "gbif", limit = 100, has_coords = TRUE)
results.data <- results[["gbif"]]$data[[formatSpName("Chelydra serpentina")]]
occs.dups <- duplicated(results.data[c('longitude', 'latitude')])
occs <- results.data[!occs.dups,]
occs$latitude <- as.numeric(occs$latitude)
occs$longitude <- as.numeric(occs$longitude)
occs$occID <- row.names(occs)
occs <- occs %>% filter(!(occID %in% 40))
output <- spThin::thin(occs, 'latitude', 'longitude', 'name', thin.par = 1, reps = 100, locs.thinned.list.return = TRUE, write.files = FALSE, verbose = FALSE)
maxThin <- which(sapply(output, nrow) == max(sapply(output, nrow)))
maxThin <- output[[ifelse(length(maxThin) > 1, maxThin[1], maxThin)]]  
occs <- occs[as.numeric(rownames(maxThin)),]  
envs <- raster::getData(name = "worldclim", var = "bio", res = 10, lat = , lon = )
envRes <- 10
if (envRes == 0.5) {
  i <- grep('_', names(envs))
  editNames <- sapply(strsplit(names(envs)[i], '_'), function(x) x[1])
  names(envs)[i] <- editNames
}
i <- grep('bio[0-9]$', names(envs))
editNames <- paste('bio', sapply(strsplit(names(envs)[i], 'bio'), function(x) x[2]), sep='0')
names(envs)[i] <- editNames
envs <- envs[[c('bio01', 'bio12')]]
locs.vals <- raster::extract(envs[[1]], occs[, c('longitude', 'latitude')])
occs <- occs[!is.na(locs.vals), ]  
xmin <- min(occs$longitude)
xmax <- max(occs$longitude)
ymin <- min(occs$latitude)
ymax <- max(occs$latitude)
bb <- matrix(c(xmin, xmin, xmax, xmax, xmin, ymin, ymax, ymax, ymin, ymin), ncol=2)
bgExt <- sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(bb)), 1)))
bgExt <- rgeos::gBuffer(bgExt, width = 0.5)
envsBgCrop <- raster::crop(envs, bgExt)
envsBgMsk <- raster::mask(envsBgCrop, bgExt)
bg.xy <- dismo::randomPoints(envsBgMsk, 10000)
bg.xy <- as.data.frame(bg.xy)  
occs.xy <- occs[c('longitude', 'latitude')]
group.data <- ENMeval::get.jackknife(occ=occs.xy, bg.coords=bg.xy)
occs.grp <- group.data[[1]]
bg.grp <- group.data[[2]]
e <- BioClim_eval(as.data.frame(occs.xy), bg.xy, occs.grp, bg.grp, envsBgMsk)
evalTbl <- e$results
evalMods <- e$models
names(e$predictions) <- "BIOCLIM"
evalPreds <- e$predictions
plot(evalMods[['BIOCLIM']], a = 1, b = 2, p = 0.9)
mod <- evalMods[["BIOCLIM"]]
pred <- evalPreds[["BIOCLIM"]]
plot(pred)
projCoords <- data.frame(x = c(-115.9918, -95.9603, -65.7374, -66.6159, -93.676, -115.9918), y = c(42.1634, 18.1459, 33.8704, 47.9899, 50.5134, 42.1634))
projPoly <- sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(projCoords)), ID=1)))
envsFuture <- raster::getData("CMIP5", var = "bio", res = 10, rcp = 45, model = "CE", year = 50)
predsProj <- raster::crop(envsFuture, projPoly)
predsProj <- raster::mask(predsProj, projPoly)
names(predsProj) <- paste0('bio', sprintf("%02d", 1:19))
predsProj <- raster::subset(predsProj, names(envs))
proj <- dismo::predict(mod, predsProj)
plot(proj)
