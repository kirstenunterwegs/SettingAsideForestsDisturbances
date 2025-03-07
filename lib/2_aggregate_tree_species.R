####################################################
#
# aggreagating the german tree species map from Blickendörfer et al.2023
# to get a higher level forest type map
#
####################################################


#---------- libraries

library(terra)
library(raster)
library(dplyr)
library(ggplot2)

#---------- set working directory

setwd("~/data/") 


#--------- load data

speciesmap<- rast("raw/tree_species_map/species_class_sum_INT1U.tif")

spec.id <- unique(speciesmap$species_class_sum)
spec.id <- spec.id$species_class_sum


# extract species share map for each species respectively

for (i in spec.id) {
  print(i)
  speciesmap.spec<- speciesmap$species_class_sum == i #extract pixels for respective species
  modified_speciesmap <- ifel(speciesmap.spec, 1, 0) #mark spec occurence as 1, other species as 0
  speciesshare.spec<- aggregate(x = modified_speciesmap, 
                                fact = 10, #aggregation factor
                                fun = function(x) sum(x == 1, na.rm = TRUE) / sum(is.finite(x))) #do not need the *100 as the divison by 100 already gives the percentage share per species
  writeRaster(speciesshare.spec,paste("processed/forest_type/species_shares/germany/speciesshare",i,"tif",sep = "."))
  
}



# load generated species share layer

spec.2.birch <- rast("processed/forest_type/species_shares/germany/speciesshare.2.tif")
spec.3.beech <- rast("processed/forest_type/species_shares/germany/speciesshare.3.tif")
spec.4.dougfir <- rast("processed/forest_type/species_shares/germany/speciesshare.4.tif")
spec.5.oak <- rast("processed/forest_type/species_shares/germany/speciesshare.5.tif")
spec.6.alder <- rast("processed/forest_type/species_shares/germany/speciesshare.6.tif")
spec.8.spruce <- rast("processed/forest_type/species_shares/germany/speciesshare.8.tif")
spec.9.pine <- rast("processed/forest_type/species_shares/germany/speciesshare.9.tif")
spec.10.larch <- rast("processed/forest_type/species_shares/germany/speciesshare.10.tif")
spec.14.fir <- rast("processed/forest_type/species_shares/germany/speciesshare.14.tif")
spec.16.odh <- rast("processed/forest_type/species_shares/germany/speciesshare.16.tif")
spec.17.odl <- rast("processed/forest_type/species_shares/germany/speciesshare.17.tif")

spec.share.stack <- c(spec.2.birch, spec.3.beech, spec.4.dougfir, spec.5.oak, spec.6.alder,
                      spec.8.spruce, spec.9.pine, spec.10.larch, spec.14.fir, spec.16.odh, spec.17.odl)
names(spec.share.stack) <- c("birch", "beech", "dougfir", "oak", "alder", "spruce", "pine", "larch", 
                             "fir", "odh", "odl")

#reduce complexity by combining alder(6), birch(2), and other broadleved species (16,17)

spec.otherdecid <- spec.2.birch + spec.6.alder + spec.16.odh + spec.17.odl
spec.otherconif <- spec.4.dougfir + spec.10.larch + spec.14.fir

spec.share.stack.2 <- c(spec.3.beech,  spec.5.oak, spec.8.spruce, spec.9.pine, spec.otherconif, spec.otherdecid)
names(spec.share.stack.2) <- c("beech", "oak", "spruce", "pine", "other.conifers", "other.decidous")


################################################################################
# reclassify species shares into forest types
################################################################################

forest_type_class_all <- function(beech, oak, spruce, pine, o.c, o.d){ #species.share.layer
  forest.type <- rast() #create empty raster to classify
  ext(forest.type) <- ext(beech)
  res(forest.type) <- res(beech)
  crs(forest.type) <- crs(beech)
  # classify change group
  #spruce dominated
  forest.type[spruce >0.66 ] <- 1 #SPRUCE
  forest.type[spruce <0.66  & is.na(oak)& is.na(beech)& is.na(pine) & is.na(o.c) & is.na(o.d)] <- 1
  #beech dominated
  forest.type[beech >0.66 ] <- 2 #BEECH
  forest.type[beech <0.66  & is.na(oak)& is.na(spruce)& is.na(pine) & is.na(o.c) & is.na(o.d) ] <- 2
  #pine dominated
  forest.type[pine >0.66] <- 3 #PINE
  forest.type[pine <0.66  & is.na(oak)& is.na(spruce)& is.na(beech) & is.na(o.c) & is.na(o.d) ] <- 3
  #oak dominated
  forest.type[oak >0.66 ] <- 4 #OAK
  forest.type[oak <0.66  & is.na(pine)& is.na(spruce)& is.na(beech) & is.na(o.c) & is.na(o.d) ] <- 4
 
  #mixed-dominant conifers
  forest.type[spruce <0.66  & pine <0.66 & sum(spruce, pine, o.c) >0.66 ] <- 5
  forest.type[spruce <0.66  & pine <0.66 & sum(spruce, pine, o.c) <0.66 & sum(spruce, pine, o.c) > sum(beech, oak, o.d) ] <- 5
  #mixed-dominant decidous
  forest.type[beech <0.66  & oak <0.66 & sum(beech, oak, o.d) >0.66 ] <- 6
  forest.type[beech <0.66  & oak <0.66 & sum(beech, oak, o.d) <0.66 & sum(beech, oak, o.d) > sum(spruce, pine, o.c) ] <- 6
  
  #mixed- around equal decidous/conifers 
  forest.type[spruce <0.66  & pine <0.66 & beech <0.66 & oak <0.66 & 
                sum(spruce, pine, o.c) >0.4 & sum(beech, oak, o.d) >0.4 ] <- 7

  
  return(forest.type)
} 

forest.type <- forest_type_class_all(spec.3.beech, spec.5.oak, spec.8.spruce, spec.9.pine, spec.otherconif, spec.otherdecid)
writeRaster(forest.type, "processed/forest_type/reclassification/germany/forest.type.class7_germany.tif")

###############################################################################


# check coverage of forest cover map of forest disturbance map

forestcover <- rast("processed/forestcover/forestcover.reclass.tif") # from European forest disturbance map (Senf & Seidl 2021)
forest.type <- rast("processed/forest_type/reclassification/germany/forest.type.class7_germany.tif")

# resample forest type to 30 m resolution for further processing

foresttype.30 <- resample(forest.type, forestcover, method = "ngb") # downsample to same res

# mask out forest cover without forest type information

forest.no.type <- mask(forestcover, foresttype.30, inverse=TRUE, updatevalue=NA)


count <- freq(forest.no.type)
count$ha <- count$count*0.09

#   layer value    count      ha
# 1     1     1 17046622 1534196-> 1,534,196 ha forest cover without forest type information (~ 12%)

count.all.forest <- freq(forestcover)
count.all.forest$ha <- count.all.forest$count*0.09

#   layer value     count       ha
# 1     1     1 141719302 12754737-> 12,754,737 ha forest cover 

###############################################################################

#---- add further forest type for forest cover base map of disturbance mapping,
#---- which is not covered by aggregated forest type map



# assign type 8 (other) forest type to forest type layer, based of where the is NA in forest type and 1 in notype
foresttype.30[notype == 1  & is.na(foresttype.30) ] <- 8


# delete forest type cells which are not covered by forest cover
foresttype.30[is.na(forestcover)] <- NA

writeRaster(foresttype.30, "processed/forest_type/reclassification/forest.type.8.germany.30_covercrop.tif", overwrite=TRUE)


# ---- crop to setaside and control(managed) areas 
setaside<- vect( "data/processed/setaside_forest_sites/all_setaside_forest_sites/setaside_forest.gpkg")
NPs <- vect("data/processed/setaside_forest_sites/NPs_Germany/NPs_Germany_25832.gpkg") # all national parks in Germany


ftype.setaside <- crop(mask(foresttype.30, setaside), setaside)
writeRaster(ftype.setaside, "data/processed/forest_type/ftype.setaside.tif", overwrite=T)

ftype.control <- mask(ftype.fcover, setaside, inverse=TRUE)
ftype.control <- mask(ftype.control, NPs, inverse=TRUE)
writeRaster(ftype.control, "data/processed/forest_type/ftype.control.tif", overwrite=T)




