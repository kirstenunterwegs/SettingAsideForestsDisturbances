##################################################################################
#
# Process raster layer with focal function to get landscape scale information
# for the matching 
#
#(1) process raster layers for matching
#(2) draw samples from setaside areas for matching
#    Draw random points within the setaside areas, which minimum distance apart
#    and where landscape extent covers enough forest area
#
#################################################################################



#---------- libraries

library(terra)
library(dplyr)
library(data.table)
library(spatialEco)
library(sf)

#---------- working directory

setwd("~/data/")


##########################################################################
# -------- (1) process raster data for landscape matching ---------------
##########################################################################

#--------- load data

setaside <- vect( "data/processed/setaside_forest_sites/all_setaside_forest_sites/setaside_forest.gpkg")



fcover.setaside <- rast("data/processed/forestcover/fcover.setaside.tif")
fcover.control <- rast("data/processed/forestcover/fcover.control.tif")

dem.setaside <- rast("data/processed/dem/dem.setaside.tif")
dem.control <- rast("data/processed/dem/dem.control.tif")

slope.setaside <- rast("data/processed/dem/slope.setaside.tif")
slope.control <- rast("data/processed/dem/slope.control.tif")

north.west.setaside <- rast("data/processed/dem/north.west.setaside.tif")
north.west.control <- rast("data/processed/dem/north.west.control.tif")

ftype.setaside <- rast("data/processed/forest_type/ftype.setaside.tif") 
ftype.control <- rast("data/processed/forest_type/ftype.control.tif") 

sitecond.setaside <- rast("data/processed/site_cond/sitecond.setaside.tif")
sitecond.control <- rast("data/processed/site_cond/sitecond.control.tif")



# rasterize set-aside forest areas for sampling

res <-  res(fcover.setaside)
extent <- ext(fcover.setaside)
crs <- crs(fcover.setaside)
raster <- rast(extent=extent, resolution=res, crs=crs)


setaside.r <- rasterize(setaside,raster ,field="newID")
setaside.area <- rasterize(setaside,raster ,field="forest.ha")



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#-------- define focal functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# mean elevation / slope / tri

mean_values <- function(vector) {
  if (all(is.na(vector))) {
    return(NA)
  } else {
    vector <- as.double(vector)
    vector <- vector[!is.na(vector)]
    return(mean(vector))
  }
}

# forest type shares
# 1= spruce, 2= beech, 3= pine, 4= oak, 5= dom. conifer, 6= dom. broadl., 7= mixed, 8= ambigious

share1 <- function(vector) {
  if (all(is.na(vector))) {
    return(NA)
  } else {
    vector <- as.double(vector)
    share <- sum(vector == 1, na.rm = TRUE) / sum(is.finite(vector))
    return(share)
  }
}

share2 <- function(vector) {
  if (all(is.na(vector))) {
    return(NA)
  } else {
    vector <- as.double(vector)
    share <- sum(vector == 2, na.rm = TRUE) / sum(is.finite(vector))
    return(share)
  }
}

share3 <- function(vector) {
  if (all(is.na(vector))) {
    return(NA)
  } else {
    vector <- as.double(vector)
    share <- sum(vector == 3, na.rm = TRUE) / sum(is.finite(vector))
    return(share)
  }
}

share4 <- function(vector) {
  if (all(is.na(vector))) {
    return(NA)
  } else {
    vector <- as.double(vector)
    share <- sum(vector == 4, na.rm = TRUE) / sum(is.finite(vector))
    return(share)
  }
}

share5 <- function(vector) {
  if (all(is.na(vector))) {
    return(NA)
  } else {
    vector <- as.double(vector)
    share <- sum(vector == 5, na.rm = TRUE) / sum(is.finite(vector))
    return(share)
  }
}

share6 <- function(vector) {
  if (all(is.na(vector))) {
    return(NA)
  } else {
    vector <- as.double(vector)
    share <- sum(vector == 6, na.rm = TRUE) / sum(is.finite(vector))
    return(share)
  }
}

share7 <- function(vector) {
  if (all(is.na(vector))) {
    return(NA)
  } else {
    vector <- as.double(vector)
    share <- sum(vector == 7, na.rm = TRUE) / sum(is.finite(vector))
    return(share)
  }
}

share8 <- function(vector) {
  if (all(is.na(vector))) {
    return(NA)
  } else {
    vector <- as.double(vector)
    share <- sum(vector == 8, na.rm = TRUE) / sum(is.finite(vector))
    return(share)
  }
}

# dominant site condition

most <- function(vector) {
  if (all(is.na(vector))) {
    return(NA)
  } else {
    vector <- as.double(vector)
    tbl <- table(vector, useNA = "no")
    return(as.double(names(tbl[which.max(tbl)])))
  }
  
}

# dominant aspect

mostAspect <- function(vector) {
  if (all(is.na(vector))) {
    return(NA)
  } else {
    #vector <- round( as.double(vector)/2,1)*2 # für 10 bins
    vector <- round( as.double(vector),1) # für 20 bins
    tbl <- table(vector, useNA = "no")
    return(as.double(names(tbl[which.max(tbl)])))
  }
  
}


# forestcover share

ncells <- function(vector) {
  length(vector[!is.na(vector)])
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#-------- create matching layers with the respective focal window, so window captures at least 20 ha of forest
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#############################
# --- set aside forests
#############################

# here we find the respective focal window size for each set aside forest area, so we find at least one 
# landscape within a focal window, which covers at least 20 ha of forest, seperately for reserve < 30 ha and > 30 ha

# --- forestcover

# create raster for merging forest focal information
res <-  res(fcover.setaside)
extent <- ext(fcover.setaside)
crs <- crs(fcover.setaside)
setaside_n.fcover<- rast(extent=extent, resolution=res, crs=crs) 

# vector for storing reserves with shape which do not allow sampling
nosample <- array()
window.df <- data.frame(matrix(ncol = 2, nrow = 0))

#i=26
#setaside.aoi.small

for(i in setaside$newID) {
  
  print(i)
  
  setaside.sub <- subset(setaside, setaside$newID == i)
  extent <- ext(setaside.sub)
  fcover.sub <- crop(fcover.setaside, extent)
  
  # check reserve size and adjust focal window
  
  if (setaside.sub$forest.ha < 30) {
    
    neighborhood_matrix <- matrix(1,nrow=17,ncol=17) # 26.01 ha 
    focal.sub <- focal(fcover.sub, w=neighborhood_matrix, fun=ncells)

        #check if min forest cells reached
        if (any(unique(focal.sub) > 220)) {
           setaside_n.fcover <- terra::merge(setaside_n.fcover, focal.sub)
           window.df <- rbind(window.df, c(i, 17))
        }
    
    # if not enough fcover sampled increase focal window 
    
    else {
      print("not enough fcover with focal window, increase focal window 17 -> 18")
      neighborhood_matrix <- matrix(1, nrow = 19, ncol = 19)
      neighborhood_matrix[19, ] <- neighborhood_matrix[, 19] <- NA # 29.16 ha
      focal.sub <- focal(fcover.sub, w=neighborhood_matrix, fun=ncells)
      
      #check if min forest cells reached
      if (any(unique(focal.sub) > 220)) {
        setaside_n.fcover <- terra::merge(setaside_n.fcover, focal.sub)
        window.df <- rbind(window.df, c(i, 18))
      }
    
        # if not enough fcover sampled increase focal window 
        
        else {
          print("not enough fcover with focal window, increase focal window 18 -> 19")
          neighborhood_matrix <- matrix(1,nrow=19,ncol=19)# 32.49 ha 
          focal.sub <- focal(fcover.sub, w=neighborhood_matrix, fun=ncells)
          
          #check if min forest cells reached
          if (any(unique(focal.sub) > 220)) {
            setaside_n.fcover <- terra::merge(setaside_n.fcover, focal.sub)
            window.df <- rbind(window.df, c(i, 19))
          }
          
          # if not enough fcover sampled increase focal window 
          
          else {
            print("not enough fcover with focal window, increase focal window 19 -> 20")
            neighborhood_matrix <- matrix(1, nrow = 21, ncol = 21)
            neighborhood_matrix[21, ] <- neighborhood_matrix[, 21] <- NA # 36 ha
            focal.sub <- focal(fcover.sub, w=neighborhood_matrix, fun=ncells)
            
            #check if min forest cells reached
            if (any(unique(focal.sub) > 220)) {
              setaside_n.fcover <- terra::merge(setaside_n.fcover, focal.sub)
              window.df <- rbind(window.df, c(i, 20))
            }
        

        # if not enough fcover sampled increase focal window 
        
            else {
              print("not enough fcover with focal window, increase focal window 20 -> 21")
              neighborhood_matrix <- matrix(1,nrow=21,ncol=21)# 39.69 ha
              focal.sub <- focal(fcover.sub, w=neighborhood_matrix, fun=ncells)
              
              #check if min forest cells reached
              if (any(unique(focal.sub) > 220)) {
                setaside_n.fcover <- terra::merge(setaside_n.fcover, focal.sub)
                window.df <- rbind(window.df, c(i, 21))
              }
              
              # if not enough fcover sampled increase focal window 
              
              else {
                print("not enough fcover with focal window, increase focal window 21 -> 22")
                neighborhood_matrix <- matrix(1, nrow = 23, ncol = 23)
                neighborhood_matrix[23, ] <- neighborhood_matrix[, 23] <- NA # 43.56 ha
                focal.sub <- focal(fcover.sub, w=neighborhood_matrix, fun=ncells)
                
                #check if min forest cells reached
                if (any(unique(focal.sub) > 220)) {
                  setaside_n.fcover <- terra::merge(setaside_n.fcover, focal.sub)
                  window.df <- rbind(window.df, c(i, 22))
                }
          
                  # if not enough fcover sampled increase focal window
                  
                  else {
                    print("not enough fcover with focal window, increase focal window 22 -> 23")
                    neighborhood_matrix <- matrix(1,nrow=23,ncol=23)# 47.61
                    focal.sub <- focal(fcover.sub, w=neighborhood_matrix, fun=ncells)
                    
                    #check if min forest cells reached
                    if (any(unique(focal.sub) > 220)) {
                      setaside_n.fcover <- terra::merge(setaside_n.fcover, focal.sub)
                      window.df <- rbind(window.df, c(i, 23))
                    }
                    
                    # if not enough fcover sampled increase focal window 
                    
                    else {
                      print("not enough fcover with focal window, increase focal window 23 -> 24")
                      neighborhood_matrix <- matrix(1, nrow = 25, ncol = 25)
                      neighborhood_matrix[25, ] <- neighborhood_matrix[, 25] <- NA # 51.84 ha
                      focal.sub <- focal(fcover.sub, w=neighborhood_matrix, fun=ncells)
                      
                      #check if min forest cells reached
                      if (any(unique(focal.sub) > 220)) {
                        setaside_n.fcover <- terra::merge(setaside_n.fcover, focal.sub)
                        window.df <- rbind(window.df, c(i, 24))
                      }
              
                          # if not enough fcover sampled increase focal window 
                          
                          else {
                            print("not enough fcover with focal window, increase focal window 24 -> 25")
                            neighborhood_matrix <- matrix(1,nrow=25,ncol=25)# 56.25 ha
                            focal.sub <- focal(fcover.sub, w=neighborhood_matrix, fun=ncells)
                            
                            #check if min forest cells reached
                            if (any(unique(focal.sub) > 220)) {
                              setaside_n.fcover <- terra::merge(setaside_n.fcover, focal.sub)
                              window.df <- rbind(window.df, c(i, 25))
                            }
              
            
                              # if not enough fcover sampled increase focal window 
                              
                              else {
                                print("not enough fcover with focal window, increase focal window 25 -> 31")
                                neighborhood_matrix <- matrix(1,nrow=31,ncol=31)# 86.49 ha 
                                focal.sub <- focal(fcover.sub, w=neighborhood_matrix, fun=ncells)
                                setaside_n.fcover <- terra::merge(setaside_n.fcover, focal.sub)
                                window.df <- rbind(window.df, c(i, 31))
                                
                                if(any(unique(focal.sub) > 220)) {next}
                                else{nosample <- c(nosample, i)
                                print("reserve (<30) shape does not allow sampling; fcover < 220") 
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
  
 
  # reserves > 30 ha
     
    if (setaside.sub$forest.ha > 30) {
      neighborhood_matrix <- matrix(1,nrow=15,ncol=15) # 20,25 ha
      focal.sub <- focal(fcover.sub, w=neighborhood_matrix, fun=ncells)

        #check if min forest cells reached
        if (any(unique(focal.sub) > 220)) {
          setaside_n.fcover <- terra::merge(setaside_n.fcover, focal.sub)
          window.df <- rbind(window.df, c(i, 15))
          
        }
      
      # if not enough fcover sampled increase focal window 
      
      else {
        print("not enough fcover with focal window, increase focal window 15 -> 16")
        neighborhood_matrix <- matrix(1, nrow = 17, ncol = 17)
        neighborhood_matrix[17, ] <- neighborhood_matrix[, 17] <- NA # 23.04 ha
        focal.sub <- focal(fcover.sub, w=neighborhood_matrix, fun=ncells)
        
        #check if min forest cells reached
        if (any(unique(focal.sub) > 220)) {
          setaside_n.fcover <- terra::merge(setaside_n.fcover, focal.sub)
          window.df <- rbind(window.df, c(i, 16))
        }
      
          # if not enough fcover sampled increase focal window (first time)
          
          else {
            print("not enough fcover with focal window, increase focal window 16 -> 17")
            neighborhood_matrix <- matrix(1,nrow=17,ncol=17)# 26.01 ha 
            focal.sub <- focal(fcover.sub, w=neighborhood_matrix, fun=ncells)
            
            #check if min forest cells reached
            if (any(unique(focal.sub) > 220)) {
              setaside_n.fcover <- terra::merge(setaside_n.fcover, focal.sub)
              window.df <- rbind(window.df, c(i, 17))
            }
            
            # if not enough fcover sampled increase focal window 
            
            else {
              print("not enough fcover with focal window, increase focal window 17 -> 18")
              neighborhood_matrix <- matrix(1, nrow = 19, ncol = 19)
              neighborhood_matrix[19, ] <- neighborhood_matrix[, 19] <- NA # 29.16 ha
              focal.sub <- focal(fcover.sub, w=neighborhood_matrix, fun=ncells)
              
              #check if min forest cells reached
              if (any(unique(focal.sub) > 220)) {
                setaside_n.fcover <- terra::merge(setaside_n.fcover, focal.sub)
                window.df <- rbind(window.df, c(i, 18))
              }
              
              # if not enough fcover sampled increase focal window 
              
              else {
                print("not enough fcover with focal window, increase focal window 18 -> 19")
                neighborhood_matrix <- matrix(1,nrow=19,ncol=19)# 32.49 ha 
                focal.sub <- focal(fcover.sub, w=neighborhood_matrix, fun=ncells)
                
                #check if min forest cells reached
                if (any(unique(focal.sub) > 220)) {
                  setaside_n.fcover <- terra::merge(setaside_n.fcover, focal.sub)
                  window.df <- rbind(window.df, c(i, 19))
                }
                
                # if not enough fcover sampled increase focal window 
                
                else {
                  print("not enough fcover with focal window, increase focal window 19 -> 20")
                  neighborhood_matrix <- matrix(1, nrow = 21, ncol = 21)
                  neighborhood_matrix[21, ] <- neighborhood_matrix[, 21] <- NA # 36 ha
                  focal.sub <- focal(fcover.sub, w=neighborhood_matrix, fun=ncells)
                  
                  #check if min forest cells reached
                  if (any(unique(focal.sub) > 220)) {
                    setaside_n.fcover <- terra::merge(setaside_n.fcover, focal.sub)
                    window.df <- rbind(window.df, c(i, 20))
                  }
                  
                  
                  # if not enough fcover sampled increase focal window 
                  
                  else {
                    print("not enough fcover with focal window, increase focal window 20 -> 21")
                    neighborhood_matrix <- matrix(1,nrow=21,ncol=21)# 39.69 ha
                    focal.sub <- focal(fcover.sub, w=neighborhood_matrix, fun=ncells)
                    
                    #check if min forest cells reached
                    if (any(unique(focal.sub) > 220)) {
                      setaside_n.fcover <- terra::merge(setaside_n.fcover, focal.sub)
                      window.df <- rbind(window.df, c(i, 21))
                    }
                    
                    # if not enough fcover sampled increase focal window 
                    
                    else {
                      print("not enough fcover with focal window, increase focal window 21 -> 22")
                      neighborhood_matrix <- matrix(1, nrow = 23, ncol = 23)
                      neighborhood_matrix[23, ] <- neighborhood_matrix[, 23] <- NA # 43.56 ha
                      focal.sub <- focal(fcover.sub, w=neighborhood_matrix, fun=ncells)
                      
                      #check if min forest cells reached
                      if (any(unique(focal.sub) > 220)) {
                        setaside_n.fcover <- terra::merge(setaside_n.fcover, focal.sub)
                        window.df <- rbind(window.df, c(i, 22))
                      }
                      
                      # if not enough fcover sampled increase focal window
                      
                      else {
                        print("not enough fcover with focal window, increase focal window 22 -> 23")
                        neighborhood_matrix <- matrix(1,nrow=23,ncol=23)# 47.61
                        focal.sub <- focal(fcover.sub, w=neighborhood_matrix, fun=ncells)
                        
                        #check if min forest cells reached
                        if (any(unique(focal.sub) > 220)) {
                          setaside_n.fcover <- terra::merge(setaside_n.fcover, focal.sub)
                          window.df <- rbind(window.df, c(i, 23))
                        }
                        
                        # if not enough fcover sampled increase focal window 
                        
                        else {
                          print("not enough fcover with focal window, increase focal window 23 -> 24")
                          neighborhood_matrix <- matrix(1, nrow = 25, ncol = 25)
                          neighborhood_matrix[25, ] <- neighborhood_matrix[, 25] <- NA # 51.84 ha
                          focal.sub <- focal(fcover.sub, w=neighborhood_matrix, fun=ncells)
                          
                          #check if min forest cells reached
                          if (any(unique(focal.sub) > 220)) {
                            setaside_n.fcover <- terra::merge(setaside_n.fcover, focal.sub)
                            window.df <- rbind(window.df, c(i, 24))
                          }
                          
                          # if not enough fcover sampled increase focal window 
                          
                          else {
                            print("not enough fcover with focal window, increase focal window 24 -> 25")
                            neighborhood_matrix <- matrix(1,nrow=25,ncol=25)# 56.25 ha
                            focal.sub <- focal(fcover.sub, w=neighborhood_matrix, fun=ncells)
                            
                            #check if min forest cells reached
                            if (any(unique(focal.sub) > 220)) {
                              setaside_n.fcover <- terra::merge(setaside_n.fcover, focal.sub)
                              window.df <- rbind(window.df, c(i, 25))
                            }
                
                          # if not enough fcover sampled increase focal window (fourth time)
                          else {
                            print("not enough fcover with focal window, increase focal window 25 -> 27")
                            neighborhood_matrix <- matrix(1,nrow=27,ncol=27) #  
                            focal.sub <- focal(fcover.sub, w=neighborhood_matrix, fun=ncells)
                            setaside_n.fcover <- terra::merge(setaside_n.fcover, focal.sub)
                            window.df <- rbind(window.df, c(i, 27))
                            
                            if(any(unique(focal.sub) > 220)) { next}
                              else{nosample <- c(nosample, i)
                              print("reserve (>30) shape does not allow sampling; fcover < 220") 
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
 

setaside_n.fcover <- mask(setaside_n.fcover, setaside)
names(window.df) <- c("newID", "window")
saveRDS(window.df, "data/processed/focal/window.df.rds") # stores focal window size per set aside area to secure aggregation of minimum 20 ha forest
writeRaster(setaside_n.fcover, "data/processed/focal/setaside_n.fcover_focal.tif", overwrite=T)

window.df <- readRDS("data/processed/focal/window.df.rds")



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# use focal window per set aside forest to aggregate environmental information 
# with focal window for all possible landscapes
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# --- elevation

# create raster for merging forest focal information
res <-  res(dem.setaside)
extent <- ext(dem.setaside)
crs <- crs(dem.setaside)
setaside_elev.mean <- rast(extent=extent, resolution=res, crs=crs) 

for (i in setaside$newID) {
  
  print(i)
  
  setaside.sub <- subset(setaside, setaside$newID == i)
  extent <- ext(setaside.sub)
  dem.sub <- crop(dem.setaside, extent)
  m <- window.df[window.df$newID == i,]$window # extract necessary window size for further aggregation
  
  # define neighbourhood matrix
  if (m %in% c(15,17,19,21,23,25,27,29,31)){
    neighborhood_matrix <- matrix(1,nrow=m,ncol=m)
  } 
  else{
    neighborhood_matrix <- matrix(1, nrow = m+1, ncol = m+1)
    neighborhood_matrix[m+1, ] <- neighborhood_matrix[, m+1] <- NA 
  }
  
  focal.sub <- focal(dem.sub, w=neighborhood_matrix, fun=mean_values)
  setaside_elev.mean <- terra::merge(setaside_elev.mean, focal.sub)
  
}

writeRaster(setaside_elev.mean, "data/processed/focal/setaside.elev.mean_focal.tif")


# --- slope

# create raster for merging forest focal information
res <-  res(slope.setaside)
extent <- ext(slope.setaside)
crs <- crs(slope.setaside)
setaside.slope.mean <- rast(extent=extent, resolution=res, crs=crs) 

for (i in setaside$newID) {
  
  print(i)
  
  setaside.sub <- subset(setaside, setaside$newID == i)
  extent <- ext(setaside.sub)
  dem.sub <- crop(slope.setaside, extent)
  m <- window.df[window.df$newID == i,]$window # extract necessary window size for further aggregation
  
  # define neighbourhood matrix
  if (m %in% c(15,17,19,21,23,25,27,29,31)){
    neighborhood_matrix <- matrix(1,nrow=m,ncol=m)
  } 
  else{
    neighborhood_matrix <- matrix(1, nrow = m+1, ncol = m+1)
    neighborhood_matrix[m+1, ] <- neighborhood_matrix[, m+1] <- NA 
  }
  
  focal.sub <- focal(dem.sub, w=neighborhood_matrix, fun=mean_values)
  setaside.slope.mean <- terra::merge(setaside.slope.mean, focal.sub)
  
}

writeRaster(setaside.slope.mean, "data/processed/focal/setaside.slope.mean_focal.tif")


# --- North-Westerness

# create raster for merging forest focal information
res <-  res(north.west.setaside)
extent <- ext(north.west.setaside)
crs <- crs(north.west.setaside)
setaside.north.west.most <- rast(extent=extent, resolution=res, crs=crs) 

for (i in setaside$newID) {
  
  print(i)
  
  setaside.sub <- subset(setaside, setaside$newID == i)
  extent <- ext(setaside.sub)
  dem.sub <- crop(north.west.setaside, extent)
  m <- window.df[window.df$newID == i,]$window # extract necessary window size for further aggregation
  
  # define neighbourhood matrix
  if (m %in% c(15,17,19,21,23,25,27,29,31)){
    neighborhood_matrix <- matrix(1,nrow=m,ncol=m)
  } 
  else{
    neighborhood_matrix <- matrix(1, nrow = m+1, ncol = m+1)
    neighborhood_matrix[m+1, ] <- neighborhood_matrix[, m+1] <- NA 
  }
  
  focal.sub <- focal(dem.sub, w=neighborhood_matrix, fun=mostAspect)
  setaside.north.west.most <- terra::merge(setaside.north.west.most, focal.sub)
  
}

writeRaster(setaside.north.west.most, "data/processed/focal/setaside.north.west.most_focal.tif")


# --- sitecondition

# create raster for merging forest focal information
res <-  res(sitecond.setaside)
extent <- ext(sitecond.setaside)
crs <- crs(sitecond.setaside)
setaside.sidecond  <- rast(extent=extent, resolution=res, crs=crs) 

for (i in setaside$newID) {
  
  print(i)
  
  setaside.sub <- subset(setaside, setaside$newID == i)
  extent <- ext(setaside.sub)
  sub <- crop(sitecond.setaside, extent)
  m <- window.df[window.df$newID == i,]$window # extract necessary window size for further aggregation
  
  # define neighbourhood matrix
  if (m %in% c(15,17,19,21,23,25,27,29,31)){
    neighborhood_matrix <- matrix(1,nrow=m,ncol=m)
  } 
  else{
    neighborhood_matrix <- matrix(1, nrow = m+1, ncol = m+1)
    neighborhood_matrix[m+1, ] <- neighborhood_matrix[, m+1] <- NA 
  }
  
  focal.sub <- focal(sub, w=neighborhood_matrix, fun=most)
  setaside.sidecond  <- terra::merge(setaside.sidecond , focal.sub)
  
}

writeRaster(setaside.sidecond, "data/processed/focal/setaside.sitecond.focal.tif")



#----- forest type shares 

#1

# create raster for merging forest focal information
res <-  res(ftype.setaside)
extent <- ext(ftype.setaside)
crs <- crs(ftype.setaside)
setaside.share1  <- rast(extent=extent, resolution=res, crs=crs) 

for (i in setaside$newID) {
  
  print(i)
  
  setaside.sub <- subset(setaside, setaside$newID == i)
  extent <- ext(setaside.sub)
  sub <- crop(ftype.setaside, extent)
  m <- window.df[window.df$newID == i,]$window # extract necessary window size for further aggregation
  
  # define neighbourhood matrix
  if (m %in% c(15,17,19,21,23,25,27,29,31)){
    neighborhood_matrix <- matrix(1,nrow=m,ncol=m)
  } 
  else{
    neighborhood_matrix <- matrix(1, nrow = m+1, ncol = m+1)
    neighborhood_matrix[m+1, ] <- neighborhood_matrix[, m+1] <- NA 
  }
  
  focal.sub <- focal(sub, w=neighborhood_matrix, fun=share1)
  setaside.share1  <- terra::merge(setaside.share1 , focal.sub)
  
}

writeRaster(setaside.share1, "data/processed/focal/setaside.share1.tif")

#2

# create raster for merging forest focal information
res <-  res(ftype.setaside)
extent <- ext(ftype.setaside)
crs <- crs(ftype.setaside)
setaside.share2  <- rast(extent=extent, resolution=res, crs=crs) 

for (i in setaside$newID) {
  
  print(i)
  
  setaside.sub <- subset(setaside, setaside$newID == i)
  extent <- ext(setaside.sub)
  sub <- crop(ftype.setaside, extent)
  m <- window.df[window.df$newID == i,]$window # extract necessary window size for further aggregation
  
  # define neighbourhood matrix
  if (m %in% c(15,17,19,21,23,25,27,29,31)){
    neighborhood_matrix <- matrix(1,nrow=m,ncol=m)
  } 
  else{
    neighborhood_matrix <- matrix(1, nrow = m+1, ncol = m+1)
    neighborhood_matrix[m+1, ] <- neighborhood_matrix[, m+1] <- NA 
  }
  
  focal.sub <- focal(sub, w=neighborhood_matrix, fun=share2)
  setaside.share2  <- terra::merge(setaside.share2 , focal.sub)
  
}

writeRaster(setaside.share2, "data/processed/focal/setaside.share2.tif")

#3

# create raster for merging forest focal information
res <-  res(ftype.setaside)
extent <- ext(ftype.setaside)
crs <- crs(ftype.setaside)
setaside.share3  <- rast(extent=extent, resolution=res, crs=crs) 

for (i in setaside$newID) {
  
  print(i)
  
  setaside.sub <- subset(setaside, setaside$newID == i)
  extent <- ext(setaside.sub)
  sub <- crop(ftype.setaside, extent)
  m <- window.df[window.df$newID == i,]$window # extract necessary window size for further aggregation
  
  # define neighbourhood matrix
  if (m %in% c(15,17,19,21,23,25,27,29,31)){
    neighborhood_matrix <- matrix(1,nrow=m,ncol=m)
  } 
  else{
    neighborhood_matrix <- matrix(1, nrow = m+1, ncol = m+1)
    neighborhood_matrix[m+1, ] <- neighborhood_matrix[, m+1] <- NA 
  }
  
  focal.sub <- focal(sub, w=neighborhood_matrix, fun=share3)
  setaside.share3  <- terra::merge(setaside.share3 , focal.sub)
  
}

writeRaster(setaside.share3, "data/processed/focal/setaside.share3.tif")

#4

# create raster for merging forest focal information
res <-  res(ftype.setaside)
extent <- ext(ftype.setaside)
crs <- crs(ftype.setaside)
setaside.share4  <- rast(extent=extent, resolution=res, crs=crs) 

for (i in setaside$newID) {
  
  print(i)
  
  setaside.sub <- subset(setaside, setaside$newID == i)
  extent <- ext(setaside.sub)
  sub <- crop(ftype.setaside, extent)
  m <- window.df[window.df$newID == i,]$window # extract necessary window size for further aggregation
  
  # define neighbourhood matrix
  if (m %in% c(15,17,19,21,23,25,27,29,31)){
    neighborhood_matrix <- matrix(1,nrow=m,ncol=m)
  } 
  else{
    neighborhood_matrix <- matrix(1, nrow = m+1, ncol = m+1)
    neighborhood_matrix[m+1, ] <- neighborhood_matrix[, m+1] <- NA 
  }
  
  focal.sub <- focal(sub, w=neighborhood_matrix, fun=share4)
  setaside.share4  <- terra::merge(setaside.share4 , focal.sub)
  
}

writeRaster(setaside.share4, "data/processed/focal/setaside.share4.tif")

#5

# create raster for merging forest focal information
res <-  res(ftype.setaside)
extent <- ext(ftype.setaside)
crs <- crs(ftype.setaside)
setaside.share5  <- rast(extent=extent, resolution=res, crs=crs) 

for (i in setaside$newID) {
  
  print(i)
  
  setaside.sub <- subset(setaside, setaside$newID == i)
  extent <- ext(setaside.sub)
  sub <- crop(ftype.setaside, extent)
  m <- window.df[window.df$newID == i,]$window # extract necessary window size for further aggregation
  
  # define neighbourhood matrix
  if (m %in% c(15,17,19,21,23,25,27,29,31)){
    neighborhood_matrix <- matrix(1,nrow=m,ncol=m)
  } 
  else{
    neighborhood_matrix <- matrix(1, nrow = m+1, ncol = m+1)
    neighborhood_matrix[m+1, ] <- neighborhood_matrix[, m+1] <- NA 
  }
  
  focal.sub <- focal(sub, w=neighborhood_matrix, fun=share5)
  setaside.share5  <- terra::merge(setaside.share5 , focal.sub)
  
}

writeRaster(setaside.share5, "data/processed/focal/setaside.share5.tif")

#6

# create raster for merging forest focal information
res <-  res(ftype.setaside)
extent <- ext(ftype.setaside)
crs <- crs(ftype.setaside)
setaside.share6  <- rast(extent=extent, resolution=res, crs=crs) 

for (i in setaside$newID) {
  
  print(i)
  
  setaside.sub <- subset(setaside, setaside$newID == i)
  extent <- ext(setaside.sub)
  sub <- crop(ftype.setaside, extent)
  m <- window.df[window.df$newID == i,]$window # extract necessary window size for further aggregation
  
  # define neighbourhood matrix
  if (m %in% c(15,17,19,21,23,25,27,29,31)){
    neighborhood_matrix <- matrix(1,nrow=m,ncol=m)
  } 
  else{
    neighborhood_matrix <- matrix(1, nrow = m+1, ncol = m+1)
    neighborhood_matrix[m+1, ] <- neighborhood_matrix[, m+1] <- NA 
  }
  
  focal.sub <- focal(sub, w=neighborhood_matrix, fun=share6)
  setaside.share6  <- terra::merge(setaside.share6 , focal.sub)
  
}

writeRaster(setaside.share6, "data/processed/focal/setaside.share6.tif")

#7

# create raster for merging forest focal information
res <-  res(ftype.setaside)
extent <- ext(ftype.setaside)
crs <- crs(ftype.setaside)
setaside.share7  <- rast(extent=extent, resolution=res, crs=crs) 

for (i in setaside$newID) {
  
  print(i)
  
  setaside.sub <- subset(setaside, setaside$newID == i)
  extent <- ext(setaside.sub)
  sub <- crop(ftype.setaside, extent)
  m <- window.df[window.df$newID == i,]$window # extract necessary window size for further aggregation
  
  # define neighbourhood matrix
  if (m %in% c(15,17,19,21,23,25,27,29,31)){
    neighborhood_matrix <- matrix(1,nrow=m,ncol=m)
  } 
  else{
    neighborhood_matrix <- matrix(1, nrow = m+1, ncol = m+1)
    neighborhood_matrix[m+1, ] <- neighborhood_matrix[, m+1] <- NA 
  }
  
  focal.sub <- focal(sub, w=neighborhood_matrix, fun=share7)
  setaside.share7  <- terra::merge(setaside.share7 , focal.sub)
  
}

writeRaster(setaside.share7, "data/processed/focal/setaside.share7.tif")

#8

# create raster for merging forest focal information
res <-  res(ftype.setaside)
extent <- ext(ftype.setaside)
crs <- crs(ftype.setaside)
setaside.share8  <- rast(extent=extent, resolution=res, crs=crs) 

for (i in setaside$newID) {
  
  print(i)
  
  setaside.sub <- subset(setaside, setaside$newID == i)
  extent <- ext(setaside.sub)
  sub <- crop(ftype.setaside, extent)
  m <- window.df[window.df$newID == i,]$window # extract necessary window size for further aggregation
  
  # define neighbourhood matrix
  if (m %in% c(15,17,19,21,23,25,27,29,31)){
    neighborhood_matrix <- matrix(1,nrow=m,ncol=m)
  } 
  else{
    neighborhood_matrix <- matrix(1, nrow = m+1, ncol = m+1)
    neighborhood_matrix[m+1, ] <- neighborhood_matrix[, m+1] <- NA 
  }
  
  focal.sub <- focal(sub, w=neighborhood_matrix, fun=share8)
  setaside.share8  <- terra::merge(setaside.share8 , focal.sub)
  
}

writeRaster(setaside.share8, "data/processed/focal/setaside.share8.tif")




# check if forest shares add up to 1

share.control <- setaside.share1 + setaside.share2 + setaside.share3 + setaside.share4 +
  setaside.share5 + setaside.share6 + setaside.share7 + setaside.share8

# create matching df

setaside.stack <- c(setaside.r, setaside.area, setaside_n.fcover, setaside_elev.mean,
                   setaside.slope.mean, setaside.north.west.most, setaside.sidecond,
                   setaside.share1, setaside.share2 , setaside.share3, setaside.share4 ,
                   setaside.share5, setaside.share6, setaside.share7, setaside.share8)



names(setaside.stack) <- c("reserve", "forest.ha", "fcover", "elev.mean", 
                          "slope.mean", "aspect", "sitecond",
                          "ftype1", "ftype2", "ftype3", "ftype4", 
                          "ftype5", "ftype6","ftype7", "ftype8")

writeRaster(setaside.stack, "data/processed/focal/setaside.stack.tif", overwrite=T)


rm(list = ls()) # empty RAM


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# --- control forest 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# we want to get landscape information of all possible managed forest landscape (here control) to find the best match/ most similar landscape
# as there will be numerous possible control forest landscape for the matching, we only aggregate with on focal window 


# depending on control area, focal function should be computed at a server or chunk up research area!

neighborhood_matrix <- matrix(1,nrow=15,ncol=15) # 15*15*0.09 = 20,25 ha

control_n.fcover <- focal(fcover.control, w=neighborhood_matrix, fun=ncells)
# elevation
control_elev.mean <- focal(dem.control, w=neighborhood_matrix, fun=mean_values)
# slope
control_slope.mean <- focal(slope.control, w=neighborhood_matrix, fun=mean_values)
# expositon (north-westerness)
control_north.west.most <- focal(north.west.control, w=neighborhood_matrix, fun=mostAspect)
# sitecondition
control.sitecond <- focal(sitecond.control, w=neighborhood_matrix, fun=most)
# forest types
control.share1 <- focal(ftype.control, w=neighborhood_matrix, fun= share1)
control.share2 <- focal(ftype.control, w=neighborhood_matrix, fun= share2)
control.share3 <- focal(ftype.control, w=neighborhood_matrix, fun= share3)
control.share4 <- focal(ftype.control, w=neighborhood_matrix, fun= share4)
control.share5 <- focal(ftype.control, w=neighborhood_matrix, fun= share5)
control.share6 <- focal(ftype.control, w=neighborhood_matrix, fun= share6)
control.share7 <- focal(ftype.control, w=neighborhood_matrix, fun= share7)
control.share8 <- focal(ftype.control, w=neighborhood_matrix, fun= share8)


share.control <- control.share1 + control.share2 + control.share3 + control.share4 +
  control.share5 + control.share6 + control.share7 + control.share8


control.stack <- c(control_n.fcover, control_elev.mean, control_slope.mean,
                   control_north.west.most, control.sitecond, 
                   control.share1, control.share2 , control.share3, control.share4 ,
                   control.share5, control.share6, control.share7, control.share8)
names(control.stack) <- c("fcover", "elev.mean", "slope.mean", 
                          "aspect", "sitecond",
                          "ftype1", "ftype2", "ftype3", "ftype4", 
                          "ftype5", "ftype6", "ftype7", "ftype8")


writeRaster(control.stack, "data/processed/focal/control.stack.update.tif", overwrite=T)

rm(list = ls()) # empty RAM


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# --- convert to data frame for matching ---
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# the matching function takes a data frame as input, so we need to convert the 
# focal function aggregated landscape stack to a data frame format
# as the stack for all of Germany is very big, we do that in chunks here:

# NOTE: doing this with data.table is more efficient and recommended!


control.stack <- rast("data/processed/focal/control.stack.tif")

# Set the chunk size (number of cells) for each subset
chunk_size <- 2000000  # Adjust this value based on your memory limitations 

# Get the total number of cells in the raster stack
total_cells <- ncell(control.stack)

# Create an empty data frame to store the final result
control.df <- data.frame(matrix(ncol = 15, nrow = 0))
colnames(control.df) <- c("fcover", "elev.mean","slope.mean", "aspect", "sitecond", "ftype1", "ftype2", "ftype3", "ftype4", "ftype5", "ftype6",    
                          "ftype7", "ftype8", "x", "y"  )

# delete all cells/ rows representing cells with forest cover <= "value assigned"
fcover_min <- 220 #20,25 ha


# -- loop to convert stack to a dataframe format and store chunk df to file

for (i in seq(1, total_cells, by = chunk_size)) {
  
  print(i)
  
  # Get the xy-coordinates of the chunk
  coords <- xyFromCell(control.stack, cell = i:min(i + chunk_size - 1, total_cells))
  
  # Subset the raster stack to a chunk of cells
  chunk <- control.stack[i:min(i + chunk_size - 1, total_cells)]
  
  # Convert the chunk to a data frame
  chunk_df <- as.data.frame(chunk)
  
  # Add the xy-coordinates to the data frame
  chunk_df$x <- coords[, 1]
  chunk_df$y <- coords[, 2]
  
  # Delete rows where fcover is less than or equal to 220
  chunk_df <- chunk_df[chunk_df$fcover > fcover_min, ]
  
  # Delete rows where all columns except 'x' and 'y' are NA
  chunk_df <- chunk_df[!(rowSums(!is.na(chunk_df[, !names(chunk_df) %in% c("x", "y")])) == 0), ]
  
  # Append the chunk data to the result data frame
  #control.df <- rbind(control.df, chunk_df)
  saveRDS(chunk_df, paste0("data/processed/matching/control_chunks/control.df.chunk", i, ".rds"))
}

#rm(list = ls())


# -- load chunks and  create one df.


# Get a list of all .rds files in the folder
rds_files <- list.files("data/data/processed/matching/control_chunks", all.files =TRUE, full.names = TRUE)
# Remove the first two paths
rds_files <- rds_files[-c(1, 2)]

# Read all .rds files and store them in a list of data frames
list_of_dataframes <- lapply(rds_files, readRDS)

# Combine the data frames using rbindlist()
combined_df <- rbindlist(list_of_dataframes)

saveRDS(combined_df, "data/data/processed/matching/control.df.rds")

rm(list = ls()) # empty RAM




##########################################################################
# -------- (2) Draw random landscape in set aside forests ---------------
##########################################################################

# --- load data

setaside.stack <- rast( "data/processed/focal/setaside.stack.tif") # environmental information inkl. forest cover for set aside areas
setaside <- vect( "data/processed/setaside_forest_sites/all_setaside_forest_sites/setaside_forest.gpkg")
sample.points  <-vect("data/processed/focal/init_point.gpkg") # need an initial point with same attribute table structure for upcoming process; I created this one manually
names(sample.points) <- c("reserve",    "forest.ha",  "fcover",     "elev.mean", "slope.mean", "aspect",     "sitecond" ,
                          "ftype1",     "ftype2",     "ftype3", "ftype4",     "ftype5",     "ftype6" ,    "ftype7" ,    "ftype8")


# ---- identify all possible landscape centroids by extracting central landscape pixels from the environmental information stack, 
# ---- which holds the landscape information aggregated with the focal functions, 
# ---- where focal function identified around 20 ha of forest within the reserve (forest cover between 217 - 274 Landsat pixel = 19.53 ha - 24.66 ha)


for (i in setaside$newID) {
  
  print(i)
  
  setaside.sub <- subset(setaside, setaside$newID == i)
  extent <- ext(setaside.sub)
  sub <- crop(setaside.stack, extent)
  points.sub <- as.points(sub, na.rm=TRUE)
  points.sub <- points.sub[points.sub$fcover >= 217 & points.sub$fcover <= 274,] 
  if(length(points.sub) == 0) {next}
  else{
  sample.points  <- rbind(sample.points, points.sub)
  }
}


# filter out initial sample point with fcover = 100
sample.points <- sample.points[sample.points$fcover >= 217,] 

# add id
sample.points$fid <- 1:nrow(sample.points)

# store all possible sample points
writeVector(sample.points, "data/processed/setaside_forest_sites/sublandscapes_setaside/sample.candidats_min_fcover217.gpkg")



# --- randomly draw landscape centroid samples within set aside forest areas for matching (comparable managed landscapes)
# --- perform distance calculation between samples to ensure limited (max 50 %) overlap of landscapes


samples  <-vect("data/processed/focal/init_point.gpkg")  
names(samples) <- c("reserve",    "forest.ha",  "fcover",     "elev.mean", "slope.mean", "aspect",     "sitecond" ,
                      "ftype1",     "ftype2",     "ftype3", "ftype4",     "ftype5",     "ftype6" ,    "ftype7" ,    "ftype8")
sample.candidats <- sample.points
window.df <- readRDS("data/processed/focal/window.df.rds") # window of landscape per reserve to capture enough forest area within landscape
sample.candidats <- merge(sample.candidats, window.df, by.x="reserve",by.y="newID", all.x=TRUE)

# loop through each reserve/set aside forest area and sample depending on the size of this area n amounts of landscape (n progressive increase with reserve size)
# using "subsample.distance" function from the spatial eco package to randomly draw samples with minimum distance within each area



for (i in unique(sample.candidats$reserve)) {
  
  # condition on forest area per reserve
  
  print(i)
  
  # < 30 ha ------ marks the size of set aside forest area/ forest reserve
  
  if (unique(sample.candidats[sample.candidats$reserve == i,]$forest.ha) < 30) 
    
    { n <- 1 # number of landscapes drawn from this reserve size
    
    # subset to reserves
    sample.candid <- sample.candidats[sample.candidats$reserve == i,]
    
    if( length(sample.candid$reserve) == 0) {print(paste0("not enough distance to another sample point in reserve ", i)); next}
    else if( length(sample.candid$reserve) == 1) {samples <- rbind(samples, sample.candid)}
    
    else{
    
    distance <- unique((sample.candid$window*30)/2)
    sample.candid <- sf::st_as_sf(sample.candid) #convert to sf object as input to subsample.distance function
    random_points <- subsample.distance(sample.candid, n, d=distance, d.max=50000, replacement = FALSE)
    random_points <- vect(random_points)

    samples <- rbind(samples, random_points)}

       }
  
  # 30 - 50 ha
  
  else if (unique(sample.candidats[sample.candidats$reserve == i,]$forest.ha) >30 & unique(sample.candidats[sample.candidats$reserve == i,]$forest.ha)  <50) # 30 -50 ha, 2 samples
    
  { n <-  2
  
  # subset to reserves
  sample.candid <- sample.candidats[sample.candidats$reserve == i,]
  
  if( length(sample.candid$reserve) == 0) {print(paste0("not enough distance to another sample point in reserve ", i)); next}
  else if( length(sample.candid$reserve) == 1) {samples <- rbind(samples, sample.candid)}
  
  else if( length(sample.candid$reserve) <= n) {
    n <- 1
    # select randomly further samples
    
    distance <- unique((sample.candid$window*30)/2)
    sample.candid <- sf::st_as_sf(sample.candid) #convert to sf object as input to subsample.distance function
    random_points <- subsample.distance(sample.candid, n, d=distance, d.max=50000, replacement = FALSE)
    random_points <- vect(random_points)
    
    #random_points <- sample.candid[sample(1:nrow(sample.candid), n), ]
    
    samples <- rbind(samples, random_points)
  }
  
  else{
    
    distance <- unique((sample.candid$window*30)/2)
    sample.candid <- sf::st_as_sf(sample.candid) #convert to sf object as input to subsample.distance function
    random_points <- subsample.distance(sample.candid, n, d=distance, d.max=50000, replacement = FALSE)
    random_points <- vect(random_points)
    
    
    samples <- rbind(samples, random_points)} 
  }

  # 50 - 100 ha
  
  else if (unique(sample.candidats[sample.candidats$reserve == i,]$forest.ha) >50 & unique(sample.candidats[sample.candidats$reserve == i,]$forest.ha)  <100) # 30 -50 ha, 2 samples
    
  { n <-  3
  
  # subset to reserves
  sample.candid <- sample.candidats[sample.candidats$reserve == i,]
  
  
  if( length(sample.candid$reserve) == 0) {print(paste0("not enough distance to another sample point in reserve ", i)); next}
  else if( length(sample.candid$reserve) == 1) {samples <- rbind(samples, sample.candid)}
  
  else if( length(sample.candid$reserve) <= n) {
    n <- 1 
    # select randomly further samples
    
    distance <- unique((sample.candid$window*30)/2)
    sample.candid <- sf::st_as_sf(sample.candid) #convert to sf object as input to subsample.distance function
    random_points <- subsample.distance(sample.candid, n, d=distance, d.max=50000, replacement = FALSE)
    random_points <- vect(random_points)
    
    #random_points <- sample.candid[sample(1:nrow(sample.candid), n), ]
    
    samples <- rbind(samples, random_points)
  }
  
  else{
    
    distance <- unique((sample.candid$window*30)/2)
    sample.candid <- sf::st_as_sf(sample.candid) #convert to sf object as input to subsample.distance function
    random_points <- subsample.distance(sample.candid, n, d=distance, d.max=50000, replacement = FALSE)
    random_points <- vect(random_points)
    
    
    samples <- rbind(samples, random_points)} 
  }
  
  # 100 - 200 ha
  
  else if (unique(sample.candidats[sample.candidats$reserve == i,]$forest.ha) >100 & unique(sample.candidats[sample.candidats$reserve == i,]$forest.ha)  <200) 
    
  { n <-  4
  
  # subset to reserves
  sample.candid <- sample.candidats[sample.candidats$reserve == i,]
  
  if( length(sample.candid$reserve) == 0) {print(paste0("not enough distance to another sample point in reserve ", i)); next}
  else if( length(sample.candid$reserve) == 1) {samples <- rbind(samples, sample.candid)}
  
  else if( length(sample.candid$reserve) <= n) {
    n <- 3
    # select randomly  samples
    
    distance <- unique((sample.candid$window*30)/2)
    sample.candid <- sf::st_as_sf(sample.candid) #convert to sf object as input to subsample.distance function
    random_points <- subsample.distance(sample.candid, n, d=distance, d.max=50000, replacement = FALSE)
    random_points <- vect(random_points)
    
    #random_points <- sample.candid[sample(1:nrow(sample.candid), n), ]
    
    samples <- rbind(samples, random_points)
  }
  
  else{
  
  distance <- unique((sample.candid$window*30)/2)
  sample.candid <- sf::st_as_sf(sample.candid) #convert to sf object as input to subsample.distance function
  random_points <- subsample.distance(sample.candid, n, d=distance, d.max=50000, replacement = FALSE)
  random_points <- vect(random_points)
  
  
  samples <- rbind(samples, random_points)} 
  }
  
  # 200 - 400 ha
  
  else if (unique(sample.candidats[sample.candidats$reserve == i,]$forest.ha) >200 & unique(sample.candidats[sample.candidats$reserve == i,]$forest.ha)  <400) # 30 -50 ha, 2 samples
    
  { n <-  5
  
  # subset to reserves
  sample.candid <- sample.candidats[sample.candidats$reserve == i,]
  
  distance <- unique((sample.candid$window*30)/2)
  sample.candid <- sf::st_as_sf(sample.candid) #convert to sf object as input to subsample.distance function
  random_points <- subsample.distance(sample.candid, n, d=distance, d.max=50000, replacement = FALSE)
  random_points <- vect(random_points)
  
  
  samples <- rbind(samples, random_points)} 
  
  # 400 - 800 ha
  
  else if (unique(sample.candidats[sample.candidats$reserve == i,]$forest.ha) >400 & unique(sample.candidats[sample.candidats$reserve == i,]$forest.ha)  <800) # 30 -50 ha, 2 samples
    
  { n <-  6
  
  # subset to reserves
  sample.candid <- sample.candidats[sample.candidats$reserve == i,]
  
  distance <- unique((sample.candid$window*30)/2)
  sample.candid <- sf::st_as_sf(sample.candid) #convert to sf object as input to subsample.distance function
  random_points <- subsample.distance(sample.candid, n, d=distance, d.max=50000, replacement = FALSE)
  random_points <- vect(random_points)
  
  
  samples <- rbind(samples, random_points)} 
  
  # 800 - 1600 ha
  
  else if (unique(sample.candidats[sample.candidats$reserve == i,]$forest.ha) >800 & unique(sample.candidats[sample.candidats$reserve == i,]$forest.ha)  <1600) # 30 -50 ha, 2 samples
    
  { n <-  7
  
  # subset to reserves
  sample.candid <- sample.candidats[sample.candidats$reserve == i,]
  
  distance <- unique((sample.candid$window*30)/2)
  sample.candid <- sf::st_as_sf(sample.candid) #convert to sf object as input to subsample.distance function
  random_points <- subsample.distance(sample.candid, n, d=distance, d.max=50000, replacement = FALSE)
  random_points <- vect(random_points)
  
  
  samples <- rbind(samples, random_points)} 
  
  # 3200 - 6400 ha
  
  else if (unique(sample.candidats[sample.candidats$reserve == i,]$forest.ha) >3200 & unique(sample.candidats[sample.candidats$reserve == i,]$forest.ha)  <6400) # 30 -50 ha, 2 samples
    
  { n <-  9
  
  # subset to reserves
  sample.candid <- sample.candidats[sample.candidats$reserve == i,]
  
  distance <- unique((sample.candid$window*30)/2)
  sample.candid <- sf::st_as_sf(sample.candid) #convert to sf object as input to subsample.distance function
  random_points <- subsample.distance(sample.candid, n, d=distance, d.max=50000, replacement = FALSE)
  random_points <- vect(random_points)
  
  
  samples <- rbind(samples, random_points)} 
  
  # 6400 - 12800 ha
  
  else if (unique(sample.candidats[sample.candidats$reserve == i,]$forest.ha) >6400 & unique(sample.candidats[sample.candidats$reserve == i,]$forest.ha)  <12800) # 30 -50 ha, 2 samples
    
  { n <-  10
  
  # subset to reserves
  sample.candid <- sample.candidats[sample.candidats$reserve == i,]
  
  distance <- unique((sample.candid$window*30)/2)
  sample.candid <- sf::st_as_sf(sample.candid) #convert to sf object as input to subsample.distance function
  random_points <- subsample.distance(sample.candid, n, d=distance, d.max=50000, replacement = FALSE)
  random_points <- vect(random_points)
  
  
  samples <- rbind(samples, random_points)} 
  
  }

# filter out initial sample point with fcover = 100
samples <- samples[samples$fcover >= 217,] #trial with 217 fcover to gain more samples

# storing sublandscape sample candidates for manual editing of sample selection for where code doesn't detetct max amount of samples
writeVector(samples, "data/processed/setaside_forest_sites/sublandscapes_setaside/samples_before_editing_min_fcover217.gpkg")



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# modify samples for final landscape set in set aside forests for matching
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

samples <- vect("data/processed/setaside_forest_sites/sublandscapes_setaside/samples_before_editing_min_fcover217.gpkg")

# In QGis:
# Load the ranomly dawn samples (samples_before_editing_min_fcover217.gpkg) in QGis with the set-aside area boundaries and forest cover layer
# the "subsample.distance" function does not necessarily find the maximum amount of landscapes in one reserve following our minimum distance criteria
# find more samples per reserves where possible and store those in a new point layer: samples_manual.gpkg

# load manually sampled landscape centroids for set aside areas where we can fit more landscapes
new.samples <- vect("data/processed/setaside_forest_sites/sublandscapes_setaside/samples_manual.gpkg")

# identify reserves where I manually collected new samples

reserve.no <- unique(new.samples$reserve)

# delete those samples in reserves in the "original/first" samples object

samples.sub <- samples[!(samples$reserve %in% reserve.no),]

# merge the new.samples object with the original one

samples <- rbind(samples.sub, new.samples)


samples.df <- as.data.frame(samples, geom="XY")


saveRDS(samples.df, "data/processed/matching/setaside.samples.df.rds")
writeVector(samples, "data/processed/matching/setaside_samples_after_editing_min_fcover217.gpkg")



