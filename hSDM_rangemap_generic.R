

###============================================================================================#
###--- Hierarchical Bayesian species distribution modeling using a rangemap as a covariate ----
###============================================================================================#

### Sami Domisch, April 2016

### This script includes: 

### - download and prepare freshwater species and environmental data 
### - move the species point records to the stream network
### - apply a simple linear distance function to the expert range map 
### - create a neighborhood matrix
### - run a GLM
### - run a spatially explicit Bayesian hierarchical model using the hSDM-package


### Please note that this script illustrates a possible workflow. To make the script easy to follow and 
### to speed up the processing, the models use (partly) randomly generated data. The script does therefore 
### not include any model evaluation of further analyses of the outputs. 


### If you find this workflow useful, please cite it as:
###       Domisch, Sami, Adam M. Wilson, and Walter Jetz. "Model-based integration of observed and 
###       expert-based information for assessing the geographic and environmental distribution of 
###       freshwater species." Ecography (2015).
browseURL("http://onlinelibrary.wiley.com/doi/10.1111/ecog.01925/abstract")


### See the hSDM manual, vignette and tutorials for more information regarding the hSDM package
browseURL("http://hsdm.sourceforge.net/")





###===========================#
### Set path and load packages
###===========================#
path="C:/hSDM_rangemap"
dir.create(path)
setwd(path)

if (!require("sp")) { install.packages("sp", dependencies = TRUE) ; library(sp)}
if (!require("raster")) { install.packages("raster", dependencies = TRUE) ; library(raster)}
if (!require("ncdf4")) { install.packages("ncdf4", dependencies = TRUE) ; library(ncdf4)} # see below for Windows
if (!require("rasterVis")) { install.packages("rasterVis", dependencies = TRUE) ; library(rasterVis)}
if (!require("maptools")) { install.packages("maptools", dependencies = TRUE) ; library(maptools)}
if (!require("foreign")) { install.packages("foreign", dependencies = TRUE) ; library(foreign)}
if (!require("rgdal")) { install.packages("rgdal", dependencies = TRUE) ; library(rgdal)}
if (!require("rgeos")) { install.packages("rgeos", dependencies = TRUE) ; library(rgeos)}
if (!require("maps")) { install.packages("maps", dependencies = TRUE) ; library(maps)}
if (!require("dismo")) { install.packages("dismo", dependencies = TRUE) ; library(dismo)}
if (!require("reshape")) { install.packages("reshape", dependencies = TRUE) ; library(reshape)}
if (!require("knitr")) { install.packages("knitr", dependencies = TRUE) ; library(knitr)}
if (!require("hSDM")) { install.packages("hSDM", dependencies = TRUE) ; library(hSDM)}
if (!require("foreach")) { install.packages("foreach", dependencies = TRUE) ; library(foreach)}
if (!require("doParallel")) { install.packages("doParallel", dependencies = TRUE) ; library(doParallel)}
if (!require("snow")) { install.packages("snow", dependencies = TRUE) ; library(snow)}
if (!require("doSNOW")) { install.packages("doSNOW", dependencies = TRUE) ; library(doSNOW)}
if (!require("coda")) { install.packages("coda", dependencies = TRUE) ; library(coda)}
if (!require("scales")) { install.packages("scales", dependencies = TRUE) ; library(scales)}
if (!require("plyr")) { install.packages("plyr", dependencies = TRUE) ; library(plyr)}

### For Windows, download the "ncdf4" library and install locally. 
### Here is an example for Windows 64-bit:
download.file("http://cirrus.ucsd.edu/~pierce/ncdf/win64/ncdf4_1.12.zip", 
               paste0(path, "/ncdf4_1.12.zip"))
install.packages(paste0(path, "/ncdf4_1.12.zip"), repos=NULL) ; library(ncdf4)

### Set color scale for plotting
col <- rev(rainbow(100, start = 0, end = 0.7)) 





###=======================================#
###--- Download examplary species data ---
###=======================================#

### Download point records for the Hardhead, Mylopharodon conocephalus
sp_points <- gbif(genus="Mylopharodon", species="conocephalus", geo=T, sp=T, removeZeros=T, download=T)
sp_points$id <- seq(1:nrow(sp_points)) # add an ID to keep track of the coordinates

### Download example rangemap of the species
if (!file.exists(paste0(path, "/f_CMC01_1.zip"))){
download.file("http://pisces.ucdavis.edu/files/uploads/layers/f_CMC01_1.zip", 
               paste0(path, "/f_CMC01_1.zip"), mode = "wb")
}
### Citation of this example data:
browseURL("http://pisces.ucdavis.edu/content/mylopharodon-conocephalus")
###     Santos, Nicholas R., et al. "A programmable information system for management 
###     and analysis of aquatic species range data in California." Environmental 
###     Modelling & Software 53 (2014): 13-26.


### Unzip file
unzip(paste0(path, "/f_CMC01_1.zip"), exdir=paste0(path, "/f_CMC01_1"), junkpaths=TRUE)

### Read the shape file. Comes in the EPSG:3310 - NAD83 / California Albers projection
### See http://spatialreference.org/ref/epsg/nad83-california-albers/
sp_range <- readShapePoly(paste0(path, "/f_CMC01_1/f_CMC01_1.shp"), verbose=T, 
                          proj4string=CRS("+proj=aea +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 
                                          +x_0=0 +y_0=-4000000 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))
### Project to WGS84
sp_range <- spTransform(sp_range, "+proj=longlat +datum=WGS84")
### Merge polygons
sp_range <- gUnaryUnion(sp_range) 
### Get the extent of the rangemap and add a buffer of 5 degrees (~500km)
ext <- extent(sp_range)
ext@xmin <- ext@xmin - 5
ext@xmax <- ext@xmax + 5
ext@ymin <- ext@ymin - 5
ext@ymax <- ext@ymax + 5



###===================================#
###--- Prepare environmental data ----
###===================================#

### Download exemplary environmental data from earthenv.org

### Upstream elevation
if (!file.exists(paste0(path, "/elevation.nc"))){
download.file("http://data.earthenv.org/streams/elevation.nc", 
               paste0(path, "/elevation.nc"), mode = "wb")
}
### Flow length and flow accumulation
if (!file.exists(paste0(path, "/flow_acc.nc"))){
download.file("http://data.earthenv.org/streams/flow_acc.nc", 
               paste0(path, "/flow_acc.nc"), mode = "wb")
}
### Upstream landcover
if (!file.exists(paste0(path, "/landcover_average.nc"))){
download.file("http://data.earthenv.org/streams/landcover_average.nc", 
               paste0(path, "/landcover_average.nc"), mode = "wb")
}

### Citation: 
browseURL("http://www.nature.com/articles/sdata201573")
###     Domisch, Sami, Giuseppe Amatulli, and Walter Jetz. "Near-global 
###     freshwater-specific environmental variables for biodiversity analyses
###     in 1 km resolution." Scientific data 2 (2015).


### Load layers and rename the single bands
dem <- brick(paste0(path, "/elevation.nc"))
names(dem) <- paste(c("dem"), c("min", "max", "avg", "range"), sep="_")

flow_acc <- brick(paste0(path, "/flow_acc.nc"))
names(flow_acc) <- paste(c("flow"), c("length", "acc"), sep="_")

lc_avg <- brick(paste0(path, "/landcover_average.nc"))
names(lc_avg) <- paste(c("lc_avg"), sprintf("%02d", seq(1:12)), sep="_")



### Crop the environmental layers (in parallel) to the extent of the buffered rangemap 
### Make cluster object
cl <- makePSOCKcluster(3)
# cl <- makePSOCKcluster(1) # if old PC use only 1 core
registerDoParallel(cl) # register parallel backend
getDoParWorkers() # show number of workers

### Crop and scale all layers in the brick
dem_crop <- foreach(i=names(dem), r = unstack(dem),  .final=stack, .packages = c("raster", "ncdf4")) %dopar% {
  options(rasterNCDF4 = TRUE)
  tmp <- crop(r, ext, snap="in")
  scale(tmp)
  }

lc_avg_crop <- foreach(i=names(lc_avg), r = unstack(lc_avg),  .final=stack, .packages = c("raster", "ncdf4")) %dopar% {
  options(rasterNCDF4 = TRUE)
  tmp <- crop(r, ext, snap="in")
  scale(tmp)
  }

flow_acc_crop <- foreach(i=names(flow_acc), r = unstack(flow_acc),  .final=stack, .packages = c("raster", "ncdf4")) %dopar% {
  options(rasterNCDF4 = TRUE)
  tmp <- crop(r, ext, snap="in")
  scale(tmp)
  }


### Create a grid_id -layer, defines the spatial entity in the models
grid_id <- raster(dem_crop[[1]])
grid_id[] <- 1:ncell(grid_id)
grid_id <- mask(grid_id, dem_crop[[1]])
names(grid_id) <- "grid_id"


### Stack the layers
env <- stack(dem_crop, lc_avg_crop, flow_acc_crop, grid_id)


### Plot the data
x11(8,10)
plot(env[[1]], col="grey")
plot(sp_range, border="red", add=T)
points(sp_points, pch=16, col="black")

### Omit one point record far up in the North, potentially a misidentification?
sp_points <- crop(sp_points, ext)


###--------------------------------------------------------------------------#
### Snap the points to the closest pixel based on a distance threshold in km
###--------------------------------------------------------------------------#

### Download the Java-Tool from phycoweb.net
if (!file.exists(paste0(path, "/moveCoordinatesToClosestDataPixel103.jar"))){
download.file("http://www.phycoweb.net/software/rasterGIS/moveCoordinatesToClosestDataPixel103.jar", 
               paste0(path, "/moveCoordinatesToClosestDataPixel103.jar"), mode = "wb")
}
### Citation: 
### Verbruggen, H. (2012) RasterTools: moveCoordinatesToClosestDataPixel.jar version 1.03, available at http://www.phycoweb.net/software
### Write the raster mask to disk (has to be an ASCII file)
writeRaster(lc_avg_crop[[1]], paste0(path, "/raster_mask.asc"), NAflag=-9999, overwrite=TRUE)

### Save the (raw) coordinates for snapping and write to disk
sp_points_df <- as.data.frame(sp_points)[-1] # remove the first column (date)
names(sp_points_df) <- c("id", "longitude", "latitude")
write.csv(sp_points_df, paste0(path, "/points_for_snap.csv"), row.names=FALSE, quote=FALSE)


### Run Java tool: You may need to set the "path" variable in the system settings, 
### see https://www.java.com/en/download/help/path.xml

# system("cmd /c  java -version") # check if Java is installed.  
# system("cmd /c  java -jar moveCoordinatesToClosestDataPixel103.jar") # see options and flags

#    -i   input coordinates file (csv)
#    -r   raster used to determine which pixels have data (esri ascii format)
#    -o   output coordinates file (csv)
# 
# optional parameters
#    -md  maximum distance that new coordinates are allowed to be from original coordinates (in km)

###--- Snapping tolerance of 3 km ----
system("cmd /c  java -jar C:/hSDM_rangemap/moveCoordinatesToClosestDataPixel103.jar  -i  C:/hSDM_rangemap/points_for_snap.csv   -r  C:/hSDM_rangemap/raster_mask.asc    -o  C:/hSDM_rangemap/points_snapped.csv  -md 3")


### Reload the snapped coordinate file. Those points that fell outside the pixels were removed. 
sp_points_snapped <- read.csv(paste0(path, "/points_snapped.csv"), h=T)
head(sp_points_snapped) # contains old and new coordinates
### Remove the old coordinate columns
sp_points_snapped <- subset(sp_points_snapped, select=-c(old_longitude, old_latitude))


### Which points were removed?
"%ni%" <- Negate("%in%") # create a "not in" -function
rows_removed <- which(sp_points$id %ni% sp_points_snapped$id) # get those ID's that were not moved to the stream grids
sp_points_removed <- sp_points[rows_removed,] # subset the raw SpatialPointsDataFrame
### Export these removed points as a shape file
# writeOGR(pts_removed, "points_removed.shp", driver="ESRI Shapefile", layer="points_removed.shp")


### Plot the raw vs. snapped points
plot(flow_acc_crop[[1]], col="grey", main="black = retained \n red = removed ") # stream network
points(sp_points_snapped[c("longitude", "latitude")], pch=16, cex=0.8, col='black') # points that were retained
points(sp_points_removed, pch=16, cex=0.8, col='red') # points that were removed



###=========================================#
###--- Create the expert rangemap decay ----
###=========================================#

### To keep the workflow generic, we only consider the connectivty of the _rangemap_ to other streams. In reality, 
### you would expect to run this also for the point records, and then potentially remove streams that are not 
### connected to any of the points, or where the (coarse) rangemap overlaps with non-connected watersheds.

### Rasterize the rangemap and mask the streams
sp_range$range <- 1
sp_range_r <- rasterize(sp_range, env[[1]], field="range", small=T, na.rm=T, background=0) # slow
sp_range_r <- mask(sp_range_r, env[[1]]) # mask ocean and streams
range_distance <- gridDistance(sp_range_r, origin=1, omit=NA) # get distance from rangemap boundary

# ### Plot with the points
# x11()
# plot(range_distance) # stream network
# points(sp_points_snapped[c("longitude", "latitude")], pch=16, col='blue')


### Calculate the decay. From the rangemap border to 150 within-stream distcance, scaled from 1 to 0
end_buffer=150000 # =150km

fun_buffer <- function(x) { ifelse( x > end_buffer, end_buffer, x)} # beyond 150km --> 0
fun_rescale <- function(x) { scales::rescale(x, to=c(1, 0)) } # within 150km --> scale from 0-1

### Apply the scaling from 0-1
range_distance_linear <- calc(range_distance, fun_buffer)
range_distance_linear <- calc(range_distance_linear, fun_rescale)
names(range_distance_linear) <- "range_distance"
# writeRaster(range_distance_linear, paste0(path, "/range_distance_linear.tif"), overwrite=T)

### Scale the range_distance layer and add to the stack
env <- addLayer(env, scale(range_distance_linear))

### Plot (aggregate cells by factor 2 for better visualization)
plot(aggregate(env[["range_distance"]], fact=2, fun=mean, na.rm=T), col=col, 
     main="Linear within-stream distance \n from rangemap boundary (scaled)") # linear distance decay
plot(sp_range, border="black", add=T)
points(sp_points_snapped[c("longitude", "latitude")], pch=16, cex=0.8, col='black')


###===============================================================================#
###--- Get the range-wide data for predicting the model across the study area ----
###===============================================================================#

data <- as.data.frame(env, na.rm=T, xy=T)
data$stream_id <- seq(1:nrow(data)) # add an id for the stream cells


###=================================================#
###--- Create neighbor-data from spatial models ----
###=================================================#

### Get neighbors for each cell. Only for the study area (=connected streams marked by "range_distance")
env_neigh <- as.data.frame(adjacent(env[["range_distance"]], cells=data$grid_id, 
                                    target=data$grid_id, 
                                    directions=8, sorted=T, pairs=T, id=T))

### Get number of neighbors for each cell
tmp_neigh <- aggregate(env_neigh$from, by=list(env_neigh$from), length)
colnames(tmp_neigh) <- c("grid_id", "n_neighbors")
unique(tmp_neigh$n_neighbors) # any disconnected cells (=0 neighbors)? if yes, need to be removed
n.neighbors <- tmp_neigh$n_neighbors; rm(tmp_neigh) # create a vector
neighbors <- as.numeric(as.factor(env_neigh$to))

### Check data
sum(n.neighbors) == length(neighbors)


###===============================================================================#
###--- Create the presence and non-detection data sets for fitting the models ----
###===============================================================================#

### Note that for illustration purposes, the non-detection as well as the number of 
### trials (=how often a site was visited) are randomly assigned

### Get presences
presence <- cbind.data.frame(sp_points_snapped[c("longitude", "latitude")], presence=1, 
                             raster::extract(env, sp_points_snapped[c("longitude", "latitude")], sp=T, ID=F))
names(presence)[1:2] <- c("x", "y")
### "Non-detections" (random)
ns <- 1000 
set.seed(1234)
sam <- as.data.frame(sampleRandom(env, ns, sp=T, ID=F))
absence <- cbind.data.frame(sam[c("x", "y")], presence=0, subset(sam, select=-c(x, y)))


### Add trials (=number of visits), here only random
presence$trials <- sample(1:3, nrow(presence), replace=T) # up to 3 visists
absence$trials <- sample(3:6, nrow(absence), replace=T) # up to 6 visists

### Merge presences and non-detections
data_fit <- rbind.data.frame(presence, absence)

### Add the stream-id
data_fit <- merge(data_fit, data[c("grid_id", "stream_id")], by="grid_id")

### Sort datafreames by "grid_id" to match neighbor-data
data_fit <- plyr::arrange(data_fit, stream_id)

### Plot all data
plot(aggregate(env[["range_distance"]], fact=2, fun=mean, na.rm=T), 
     main="black = presence \ngrey = non-detection (random)",
     xlab="Longitude", ylab="Latitude", col=col)
points(absence[c("x", "y")], pch=16, cex=0.4, col='grey')
plot(sp_range, border="black", add=T)
points(presence[c("x", "y")], pch=16, cex=0.8, col='black')




###==================#
###--- Run a GLM ----
###==================#

### Use upstream elevation, flow accumulation, landcover (forest cover) and the rangemap as predictors
model <- "~ dem_avg + dem_range + flow_acc + lc_avg_01 + lc_avg_04 + lc_avg_05 + range_distance"
mod_glm <- glm(paste0("presence ", model), data = data_fit, family = "binomial") 
beta_hat_glm <- coef(mod_glm) # get betas as starting values for hSDM later
(as.data.frame(coef(mod_glm)))

### Predict model across the study area
pred_glm <- predict(env, mod_glm, type="response")
### Write to disk
writeRaster(pred_glm, paste0(path, "/pred_glm.tif"), overwrite=T)

### Plot
plot(aggregate(pred_glm, fact=2, fun=mean, na.rm=T), zlim=c(0,1), col=col, 
            main="Prediction GLM, aggregated x 2"); 
plot(sp_range, border="black", add=T)
points(sp_points_snapped[c("longitude", "latitude")], pch=16, cex=0.8, col="black")

###=================================================================================#
###--- Run zero-inflated binomial model with a conditional autoregressive model ---- 
###=================================================================================#

### For demonstration purpose, run one chain with only 1000 iterations. If the species was 
### found at least once at a given site, the environment is considered suitable and any 
### non-detections at that site are due to imperfect detection.

### Note that in this example the non-detections and trials are randomly generated; the 
### resulting models are therefore not useful for any ecological inference but for testing purpose only
mod_hSDM_ZIB_iCAR <- hSDM.ZIB.iCAR(presences = data_fit$presence,
                                  trials = data_fit$trials,
                                  suitability = model,
                                  observability = ~ 1,
                                  spatial.entity = data_fit$stream_id,  
                                  data = data_fit,
                                  n.neighbors = n.neighbors,
                                  neighbors = neighbors,
                                  ## suitability.pred=NULL,
                                  ## spatial.entity.pred=NULL,
                                  suitability.pred = data,
                                  spatial.entity.pred = data$stream_id, 
                                  ### Chains
                                  burnin=1000, mcmc=1000, thin=1, 
                                  ### Starting values
                                  beta.start=beta_hat_glm,  #0
                                  gamma.start=0,
                                  Vrho.start=10, #10
                                  ### CAR process
                                  #  priorVrho="1/Gamma",
                                  priorVrho="Uniform",
                                  #  priorVrho=10,
                                  ### Priors
                                  mubeta=0, Vbeta=100,
                                  mugamma=0, Vgamma=100,
                                  #  shape=2, #0.5, 
                                  # rate=1,  #0.0005, 
                                  Vrho.max=10,
                                  ## Misc
                                  seed=1234, verbose=1,
                                  save.rho=0, save.p=1)


### Check results
str(mod_hSDM_ZIB_iCAR)
summary(mod_hSDM_ZIB_iCAR$mcmc)

### Get traceplots
x11(); plot(mod_hSDM_ZIB_iCAR$mcmc)

### Get model predictions into a dataframe
pred_hbm_df <- data[c("x", "y")]
pred_hbm_df$mean_suitability <- as.numeric(apply(mod_hSDM_ZIB_iCAR$prob.p.pred, 2, mean))
pred_hbm_df$mean_spatial_random <- as.numeric(mod_hSDM_ZIB_iCAR$rho.pred)
### If "save.p=1" in the function, extract also the lower and upper credible intervals
pred_hbm_df$CI_2.5  <- apply(mod_hSDM_ZIB_iCAR$prob.p.pred, 2, quantile, 0.025)
pred_hbm_df$CI_97.5  <- apply(mod_hSDM_ZIB_iCAR$prob.p.pred, 2, quantile, 0.975)

### Create maps
pred_hbm <- stack(
rasterFromXYZ(pred_hbm_df[c(1:2,3)]), # x, y, mean suitability
rasterFromXYZ(pred_hbm_df[c(1:2,4)]),
rasterFromXYZ(pred_hbm_df[c(1:2,5)]),
rasterFromXYZ(pred_hbm_df[c(1:2,6)])
)


### Plot the predictions. Note the very high spatial random effects (scale truncated in the plot), 
### potentially due to the stream network structure 
x11(8,10); par(mfrow=c(2,2))
plot(aggregate(pred_hbm[[1]], fact=2, na.rm=T), col=col, zlim=c(0,1), main="mean suitability")
plot(sp_range, border="black", lwd=2, add=T)
# points(sp_points_snapped[c("longitude", "latitude")], pch=16, cex=0.8, col="black")
plot(aggregate(pred_hbm[[2]], fact=2, na.rm=T), col=col, zlim=c(-20,20), main="spatial random effects")
plot(aggregate(pred_hbm[[3]], fact=2, na.rm=T), col=col, zlim=c(0,1), main="lower credible interval")
plot(aggregate(pred_hbm[[4]], fact=2, na.rm=T), col=col, zlim=c(0,1), main="upper credible interval")

### Estimated detection probability of the species
parameter <- summary(mod_hSDM_ZIB_iCAR$mcmc)$statistics
detection_prob <- parameter[,1]["gamma.(Intercept)"]
inv_logit <- inv.logit(detection_prob)
cat("Detection probability is", round(inv_logit,2))

### Write the 4 raster files to disk
writeRaster(pred_hbm, paste0(path, "/hbm.tif"), overwrite=T, bylayer=T, suffix=names(pred_hbm))

### Stop the cluster object
stopCluster(cl)
# graphics.off()

