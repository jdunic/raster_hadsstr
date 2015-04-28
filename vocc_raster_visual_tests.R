setwd('Meta_analysis_ms')
library(hadsstR)
library(chron)
library(ncdf)
library(sp)
library(beepr)
library(raster)
library(lattice)

dev.new()
dev.new()

hadsst_raster <- loadHadSST1('Data/', hadsstFilename = "HadISST_sst.nc")
sstData <- loadHadSST(directory="Data/", hadsstFilename="HadISST_sst.nc")

cMats <- getClimateChange(sstData, years = 1960:2009)
beep()
AllRasters <- getAllRasters(hadsst_raster, years = 1960:2009)
beep()

#Let's plot some of these matrices
pal <- colorRampPalette(c("blue","white", "red"))

dev.set(2)
with(cMats, image(lon, lat, averageMat, col=pal(80)))
dev.set(3)
plot(subset(AllRasters, 'Average'), col=pal(80), ylim = c(-90, 90), colNA = 'black')

dev.set(2)
with(cMats, image(lon, lat, linearChangeMat, col=pal(80)))
dev.set(3)
plot(subset(AllRasters, 'LinearChange'), col=pal(80), ylim = c(-90, 90), colNA = 'black')

dev.set(2)
with(cMats, image(lon, lat, spatialGradMat, col=pal(80)))
dev.set(3)
plot(subset(AllRasters, 'SpatialGrad'), col=pal(80), ylim = c(-90, 90), colNA = 'black')

dev.set(2)
with(cMats, image(lon, lat, velocityMat, col=pal(80)))
dev.set(3)
plot(subset(AllRasters, 'Velocity'), col=pal(80), ylim = c(-90, 90), colNA = 'black')


#create a velocity matrix where values >200 and < -200 are truncated to those limits
#for easier plotting, as in Burrows et al. 20011
velMatTruncated <- cMats$velocityMat
velMatTruncated[velMatTruncated >200] <- 200
velMatTruncated[velMatTruncated < -200] <- -200

dev.set(2)
latLonGrid <- expand.grid(lon = climateChangeMats$lon, lat = climateChangeMats$lat)
levelplot(velMatTruncated ~ lon * lat, data = latLonGrid, #at = cutpts, 
           pretty = T, 
          col.regions = pal(100),
           at=seq(-200,200,length.out=100))
dev.set(3)

vel_rast <- selectRaster(AllRasters, 'Velocity')
m <- c(-Inf, -200, -200, 200, Inf, 200)
rclmat <- matrix(m, ncol = 3, byrow = TRUE)
velRastTruncated <- raster::reclassify(vel_rast, m)


plot(velRastTruncated, zlim = c(-200, 200), col = pal(100))