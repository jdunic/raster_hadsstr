# Dependencies:
# raster, sp, chron
library(raster)
library(sp)
library(chron)

loadHadSST1 <- function(directory="./", hadsstFilename="HadISST_sst.nc") {
    f1 <- paste0(directory, hadsstFilename)
    b <- raster::brick(f1)
    raster::NAvalue(b) <- -1000
    return(b)
}

# Gets the average temperature for each year, given a single year or a range of 
# years.
# Returns a raster brick object with the average of each year.
# Probably could rewrite this with a stackApply instead of apply wrapping 
# raster::calc()
getSSTAnnualRasters <- function(hadsst_raster, years = 1969:2011) {
    mean_rasts <- 
    apply(matrix(years), 1, function(x) {
    #browser()
        yearIDx <- which(chron::years(hadsst_raster@z$Date) == x)
        subset_x <- raster::subset(hadsst_raster, yearIDx) 
        means <- mean(subset_x, na.rm = TRUE)
        names(means) <- as.character(x)
        return(means)
    })
    #return(mean_rasts)
    mean_brick <- raster::brick(mean_rasts)
    mean_brick <- raster::setZ(mean_brick, as.Date(paste0(years, '-01-01')), 'Date')
    return(mean_brick)
}


# Gets the average temperature over a range of years.
# Returns a single raster object.
getSSTAvgOverYears <- function(hadsst_raster = b, years = 1969:2011) {
    yearIDs <- which(chron::years(hadsst_raster@z$Date) %in% years)
    subset_x <- raster::subset(hadsst_raster, yearIDs)
    SSTAvgOverYears <- mean(subset_x, na.rm = TRUE)
    return(SSTAvgOverYears)
    }


# Gets a raster containing the linear change in temperature over time.
# Returns a single raster object that contains the decadal rate of change of
# temperature change.
getSSTLinChangeRaster <- function(hadsst_raster = b, years = 2000:2011) {
    annual_rasters <- getSSTAnnualRasters(hadsst_raster, years)

    time_ <- I(years - mean(years))
    fun <- function(x) { 
        #browser()
        if (sum(is.na(x)) < length(x)) {
            slope <- 10 * lm(x ~ time_)$coefficients[2]
            return(slope)
        }
        return(NA)
    }
    x2 <- raster::calc(annual_rasters, fun)
    return(x2)
}

# Gets the west to east differences in temperature (corrected for distance by 
# latitude).
# Returns a single raster containing the adjusted differences.
getWEdiffs <- function(avg_raster) {
    f <- matrix(rep(1, 9), nrow=3, ncol=3)
    lat <- sp::coordinates(avg_raster)[, 2]
    WEdiff_raster <- raster::focal(avg_raster, w = f, nrow=3, ncol=3, pad = TRUE,
                          fun = function(x, ...) {
                          ba <- x[6] - x[5]
                          }
    )
    WEdiff_adj <- WEdiff_raster / (111.325 * cos (lat * pi / 180))
}

# Gets the north-south differences in temperature.
# Returns
getNSdiffs <- function(avg_raster) {
    f <- matrix(rep(1, 9), nrow=3, ncol=3)
    NSdiff_raster <- raster::focal(avg_raster, w = f, nrow=3, ncol=3, pad = TRUE,
                          fun = function(x, ...) {
                          #browser()
                          be <- x[2] - x[5]
                          }
    )
    NSdiff_adj <- NSdiff_raster / 111.325
}


# Get spatial gradient raster.
# Returns a raster brick where the first layer contains the magnitude of the 
# spatial gradient and the second layer contains the angle of the vector.
getSpatialGrad_raster <- function(NSraster, WEraster) {
    f <- matrix(rep(1, 9), nrow=3, ncol=3)
    xGrad <- raster::focal(WEraster, w = f, nrow=3, ncol=3, pad = TRUE,
                          fun = function(x, ...) {
                          #browser()
                          wt <- c(1, 2, 1, 2, 0, 2, 1, 2, 1)
                          weighted.mean(x, wt, na.rm = TRUE)
                          })
    yGrad <- raster::focal(NSraster, w = f, nrow = 3, ncol = 3, pad = TRUE, 
                           fun = function(x, ...) {
                           wt <- c(1, 2, 1, 2, 0, 2, 1, 2, 1)
                           weighted.mean(x, wt, na.rm = TRUE)
                           })

    # Get magnitude of resultant vector
    vecSum_rast <- sqrt(xGrad^2 + yGrad^2)

    # Get angle of resultant vector
    vecAngle <- atan2(xGrad, yGrad) * 180 / pi

    # Create correction raster to produce positive angles
    negAngle <- (vecAngle < 0) * 360
#browser()

    vecAngle <- 
        raster::overlay(x = vecAngle, y = negAngle, fun = function(x, y) {x + y})

    spatial_grad_brick <- raster::brick(vecSum_rast, vecAngle)
    names(spatial_grad_brick) <- c('spatialGrad', 'angle') 
    return(spatial_grad_brick)
}


# Gets the magnitude of velocity of climate change raster over a set of years.
# Returns a single raster containing the magnitudes of decadal velocities. To 
# get the direction of velocity get the getSpatialGrad_raster 'angle' raster.
getVelocityMag_raster <- function(hadsst_raster, years = 1969:2009, truncate = TRUE) {
    #browser()
    LinChangeRaster <- getSSTLinChangeRaster(hadsst_raster, years)

    AvgRaster <- getSSTAvgOverYears(hadsst_raster, years)

    WEraster <- getWEdiffs(AvgRaster)
    NSraster <- getNSdiffs(AvgRaster)
    
    SpatGradRaster <- getSpatialGrad_raster(NSraster, WEraster)

    VelocityMagRaster <- 
        LinChangeRaster / raster::subset(SpatGradRaster, 'spatialGrad')

    if (truncate == TRUE) {
        m <- c(-Inf, -200, -200, 200, Inf, 200)
        rclmat <- matrix(m, ncol = 3, byrow = TRUE)
        VelocityMagRaster <- raster::reclassify(VelocityMagRaster, m)
    }

    return(VelocityMagRaster)
}

getAllRasters <- function(hadsst_raster, years = 1969:2009) {
    LinChangeRaster <- getSSTLinChangeRaster(hadsst_raster, years)

    AvgRaster <- getSSTAvgOverYears(hadsst_raster, years)

    WEraster <- getWEdiffs(AvgRaster)
    NSraster <- getNSdiffs(AvgRaster)
    
    SpatGradRaster <- getSpatialGrad_raster(NSraster, WEraster)

    VelocityMagRaster <- 
        LinChangeRaster / raster::subset(SpatGradRaster, 'spatialGrad')

    AllRasters <- raster::brick(AvgRaster, LinChangeRaster, 
                                SpatGradRaster, VelocityMagRaster
                                )
    names(AllRasters) <- c('Average', 'LinearChange', 'SpatialGrad', 
                                   'Angle', 'Velocity') 
    #browser()
    return(AllRasters)
}

selectRaster <- function(AllRasters_obj, raster_name) {
    selected_raster <- raster::subset(AllRasters_obj, raster_name)
    return(selected_raster)
}


# Test
#start1 <- proc.time()
#hadsst_raster <- loadHadSST1('Data/', hadsstFilename = "HadISST_sst.nc")
#velocity <- getVelocityMag_raster(hadsst_raster, years = 1935:2011)
#end1 <- proc.time() - start1
#beep()

# Compare with time to load 
#start2 <- proc.time()
#sstData <- loadHadSST(directory="Data/", hadsstFilename="HadISST_sst.nc")
#cMats <- getClimateChange(sstData, years = 1935:2011)
#end2 <- proc.time() - start2
#beep()


m <- c(-Inf, -200, -200, 200, Inf, 200)
rclmat <- matrix(m, ncol = 3, byrow = TRUE)
velRastTruncated <- raster::reclassify(vel, m)
plot(velRastTruncated, zlim = c(-200, 200), col = pal(100))