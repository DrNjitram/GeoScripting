GetSubNDVI <- function(NDVI_7, NDVI_8, NDVI_delta, ned, munin_name){
  
  "
  This function provides three rasters, of each of the three layers provided, cropped and masked to
  the munincipality provided, extracted from the netherlands shapefile
  
  The three rasterLayers do not need to 100% overlap, but all must overlap the munincipality name
  
  Arguments:
    NDVI_7: Raster*
    NDVI_8: Raster*
    NDVI_delta: Raster*
    ned: SpatialPolygonsDataframe
    munin_name: String

  returns:
    rasters: RasterStack 
  "
  # Retrieve the munincipality file
  west <- ned[ned$NAME_2 == munin_name,  ]
  
  # Mask and crop the file, so that we only have the data for the mask
  N_7_w <- crop(NDVI_L7, west)
  N_7_w <- mask(N_7_w, mask = west)

  N_8_w <- crop(NDVI_L8, west)
  N_8_w <- mask(N_8_w, mask = west)
  
  # Mutually mask them so any gaps are copied
  N_8_w <- mask(N_8_w, mask = N_7_w)
  N_7_w <- mask(N_7_w, mask = N_8_w)
  
  # Do the same for the difference raster
  comp <- crop(NDVI_delta, west)
  comp <- mask(comp, mask = west)
  # We dont have to mutually mask because any gaps will already be in this file
  # since NA -+ x = NA
  
  #Merge them all into the same stack, ensuring that they will have the same extent, otherwise its error town
  rasters <- stack(N_7_w, N_8_w, comp)
  
  return(rasters)
}

GetDifferenceCells <- function(N_s, breaks = c(0.3, 0.7)){
  "
  This function calculates the amount of area that lied above the last break in the layers
  and then examines how much has left and entered this break.

  Default breaks indicate forest (>0.7), grassland (0.3 < x < 0.7) and housing/roads (<0.3)

  Arguments:
    N_s: RasterStack or RasterBrick, only the first two are used
    breaks: Vector consisting of two values with them being 0 < x1 < x2 < 1

  returns:
    df: DataFrame
  "
  
  # Extract the two relevant RasterLayers
  N_7_w <- N_s[[1]]
  N_8_w <- N_s[[2]]
  
  # Create a mask with only values above break 2 and NA them 
  N_7_forest <- N_7_w > breaks[2]
  N_7_forest[N_7_forest == 0] <- NA
  
  # Get the area from the second Layer that matches the mask
  N_8_forest <- mask(N_8_w, mask = N_7_forest)
  
  # Get the cells from the second layer that are still above the break and NA the 0's to only count valid ones
  N_8_maintained_forest <- N_8_forest > breaks[2]
  N_8_maintained_forest[N_8_maintained_forest == 0] <- NA
  
  # Get the cells from the second layer that are below the first break and NA the 0's to only count valid ones
  N_8_removed_forest <- N_8_forest < breaks[1]
  N_8_removed_forest[N_8_removed_forest == 0] <- NA
  
  # We repeat the above, but reverse the layers, to get the amount of cells that went into the top break, or were already there.
  N_8_forest_2 <- N_8_w > breaks[2]
  N_8_forest_2[N_8_forest_2 == 0] <- NA
  
  N_7_forest_2 <- mask(N_7_w, mask = N_8_forest_2)
  
  N_7_created_forest <- N_7_forest_2 < breaks[1]
  N_7_created_forest[N_7_created_forest == 0] <- NA
  
  # Count the cells that have a value
  old_total <- sum(!is.na(N_7_forest[]))
  
  # Count the cells that were already there, passed out of the top break and into the last, and ones that passed out of the top and into the second
  maintained <-sum(!is.na(N_8_maintained_forest[]))
  removed <- sum(!is.na(N_8_removed_forest[]))
  reduced <- old_total - removed - maintained
  
  # Repeat the same
  # We dont recalculate the maintained, since its the same as above by definition
  new_total <-sum(!is.na(N_8_forest_2[]))
  
  created <- sum(!is.na(N_7_created_forest[]))
  improved <- new_total - created - maintained
  
  # Put both sets of values into a vector, and convert it to area in km2
  comp_old <- c(old_total, maintained, reduced, removed)
  comp_old <- (comp_old * (xres(N_8_w) * yres(N_8_w)))/(1000*1000)
  
  comp_new <- c(new_total, maintained, improved, created)
  comp_new <- (comp_new * (xres(N_8_w) * yres(N_8_w)))/(1000*1000)
  
  # Create a Data frame and give it names 
  df <- data.frame(comp_old, comp_new)
  row.names(df) <- c("total", "maintained", "reduced/improved", "removed/created")
  return(df)
}

GetDfDelta <- function(ned, munins, NDVI_delta){
  "
  This function provided information on all munincipalities provided about the NDVI layer.
  Currently only provides the sum divided by area and the standard deviation

  Arguments:
    ned: shapefile
    munins: vector of strings, needs to be in ned$NAME_2 
    NDVI_delta: Rasterlayer

  returns:
    NDVI_df: Dataframe with rows: munincipality_name, diff_per_km2, standard_deviation
  "
  # Set the iterator
  i <- 1
  
  # Create and empty dataframe to hold all values
  NDVI_df <- data.frame(matrix(ncol = 3, nrow = 0))
  
  # Iterate over all munincipalities
  for(munin in munins){
    #Retrieve the shapefile for the munincipality in question
    west <- ned[ned$NAME_2 == munin,  ]
    
    # Gte the delta NDVI specififcally for this municipality
    comp <- crop(NDVI_delta, west)
    comp <- mask(comp, mask = west)
    
    # Get the area by counting all cells with values and multiplying this by the area per pixel
    # We divide by 1000^2 to get the area in square kilometers
    area <- (sum(!is.na(comp[])) * (xres(comp)*yres(comp)))/(1000*1000)
    
    #We extract the sum of all NDVI_deltas to get the net difference in vegetation for this munincipality
    total_ndvi <- sum(comp[!is.na(comp[])])
    
    # Add the difference per km^2 and its standard deviation to a subframe with the munincipality name 
    sub_frame <- data.frame(munin, total_ndvi/area, sd(comp[!is.na(comp[])]))
    # Add the subframe to the end of the dataframe
    NDVI_df <- rbind(NDVI_df, sub_frame)
    i <- i + 1
  }
  return(NDVI_df)
}