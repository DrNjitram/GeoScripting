#Get the filenames
require(raster)
require(ggplot2)

require(rgdal)
require(rasterVis)
require(RColorBrewer)

require("rgeos")
require("maptools")
require("plyr")
require("gridExtra")
require("ggsn")


# Make sure to set the working directory correctly
# check with getwd()
# and set with setwd()
# Also make sure to have at least 20 GB space for temp files
# If needed, set the temp dir with rasterOptions(tmpdir = "a path")

# Load the Scripts
source("./R/scripts.R")

# Create the data dir if it doesnt yet
dir.create("./data", showWarnings = FALSE)

# Retrieve the shapefile for the netherlands and its munincipalities we will use alter
ned <- getData('GADM', country='NLD', level=2, path = './data')
# Set all munincipality names that are NA to be empty strings, this is needed since a function later doesnt accept NA's
ned[is.na(ned$NAME_2), ] <- ""

# This variable controls if you want to remake the processes files again regardless if they already exist
# WARNING: Takes a while to do and will download large files
force_process <- FALSE
# If the prre-processed files are needed, a drive link is in the links.txt file

# Check if our completed files already exists and if so, do not reprocess or redownload them
if((file.exists("./L7_all_masked.grd") == FALSE & file.exists("./L8_complete.grd") == FALSE) | force_process){
  
  # This is a bit of a clunky way to download, but the api to do so is down since the 
  # US government, who runs the USGS, has shut down
  # NOTE: These files are only available for a week after the 28th of januari
  # If not avialable, find link to these files in the link.txt
  if (!file.exists("./data/LC081980242018072601T1-SC20190128104132.tar.gz")){
    download.file(url = "https://edclpdsftp.cr.usgs.gov/orders/espa-martijn.gobes@wur.nl-0101901281910/LC081980242018072601T1-SC20190128104132.tar.gz", destfile = "./data/LC081980242018072601T1-SC20190128104132.tar.gz", method = "auto")
  }
  if (!file.exists("./data/LC081980232018072601T1-SC20190128104019.tar.gz")){
    download.file(url = "https://edclpdsftp.cr.usgs.gov/orders/espa-martijn.gobes@wur.nl-0101901281910/LC081980232018072601T1-SC20190128104019.tar.gz", destfile = "./data/LC081980232018072601T1-SC20190128104019.tar.gz", method = "auto")
  }
  if (!file.exists("./data/LE071990242002072901T1-SC20190128104708.tar.gz")){
    download.file(url = "https://edclpdsftp.cr.usgs.gov/orders/espa-martijn.gobes@wur.nl-0101901281918/LE071990242002072901T1-SC20190128104708.tar.gz", destfile = "./data/LE071990242002072901T1-SC20190128104708.tar.gz", method = "auto")
  }
  if (!file.exists("./data/LE071990232002072901T1-SC20190128105017.tar.gz")){
    download.file(url = "https://edclpdsftp.cr.usgs.gov/orders/espa-martijn.gobes@wur.nl-0101901281918/LE071990232002072901T1-SC20190128105017.tar.gz", destfile = "./data/LE071990232002072901T1-SC20190128105017.tar.gz", method = "auto")
  }
  
  # Retrieves a list of all .tar files that exist in the directory
  files <- list.files('./data', full.names = TRUE, pattern = glob2rx("*tar*"))
  
  # For each of the afformentioned files, untar them in their own subfolder
  for(path in files){
    subfolder <- paste0("./data/", substring(path, 8, 47))
    dir.create(subfolder, showWarnings = FALSE)
    untar(path, exdir = subfolder)
  }
  
  # Create two empty lists to hold the raster stacks
  L7 <- list(rep(0, length(files) - 2))
  L8 <- list(rep(0, 2))
  
  # Set both iterators for L7 and L8
  i <- 1
  j <- 1
  # Initialise the variable that denotes if we have level 2 or level 1 data
  # Level 2 data is marked by having SC in the filename
  level_2 <- FALSE
  
  for(path in files){
    #Retrieve the names of all tifs from their respective folders
    folder <- paste0("./data/", substring(path, 8, 46))
    tifs <- list.files(folder, full.names = TRUE, pattern =  "^.*(\\.tif|\\.TIF)")
    
    # Determine if it is an L7 or L8 file, since they have different band orders
    # Then, depending on if it is level 1 or level 2, order the bands 1 - x and then add quality or other bands
    # Note: From L8 level 1 band 8 is omitted since this is a black and white band which has a different resolution than the others
    if(grepl("LE07", path, fixed = TRUE)){
      if(grepl("SC", path, fixed = TRUE)){ # This is the check for level 2 data
        L7[[i]]  <- stack(tifs[c(4, 5, 6, 7, 8, 9, 10, 3, 2, 1)])
        level_2 <- TRUE
      } else {
        L7[[i]]  <- stack(tifs[c(1, 2, 3, 4, 5, 6, 7, 8, 10)])
      }
      i <- i + 1
    }
    if(grepl("LC08", path, fixed = TRUE)){
      if(grepl("SC", path, fixed = TRUE)){
        L8[[j]] <- stack(tifs[c(4, 5, 6, 7, 8, 9, 10, 3, 2, 1)])
        level_2 <- TRUE
      } else {
        L8[[j]] <- stack(tifs[c(1, 4, 5, 6, 7, 8, 9, 11, 2, 3, 12)])
      }

      j <- j + 1
    }
  }
  
  
  # Check if the L8 file already exists or if force process is marked
  if(file.exists("./L8_level2_complete.grd") == FALSE | force_process){
    # Ensure that the netherlands shapefile is in the same crs
    # The ned object has a lat long projection, so this needs to be changed
    ned <- spTransform(ned, CRSobj = crs(L8[[1]]))
    
    #Crop both rasters to the extent of the netherlands and set all values that are bad (< 1) to be transparent
    L8[[1]] <- crop(L8[[1]], ned)
    L8[[1]][L8[[1]] < 1] <- NA
    
    L8[[2]] <- crop(L8[[2]], ned)
    L8[[2]][L8[[2]] < 1] <- NA  
    
    # If these are level 2 files, the masking is done by different values and on different layers, but the process is the same for both for the rest
    # First get all the pixels to be masked, then mark the rest as NA, so they are transparent
    # Then mask the raster with its respective cloudmask
    if(level_2){
      breakoff <- 400
      layer <- 10
    } else {
      breakoff <- 2750
      layer <- 11
    }

    for(i in 1:length(L8)){
      tmp <- L8[[i]] # Create a tmp layer
      cloud <- tmp[[layer]] > breakoff #Create the cloudmask from the cloud layer
      cloud[cloud == 1] <- NA #Set all values that are 0 to be NA
      L8[[i]] <- mask(x = tmp, mask = cloud) #Do the actual masking and put it back into the list
    }
      
    
    # Merge both raster objects together into one large picture of the netherlands
    # Support for more files can be gotten by iterating over the list
    L8_merged <- merge(L8[[1]], L8[[2]])
    
    #Mask the merged file with the ned shapefile, to limit our pictures to be exactly the netherlands
    L8_masked_merged <- mask(x = L8_merged, mask=ned)
    
    # Set the names of bands of the rasterobject
    if(level_2){
      names(L8_masked_merged) <- c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "BAER", "BRQA", "BPQA")
    } else {
      names(L8_masked_merged) <- c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "B9", "B10", "B11", "BQA")
    }
    
  } else if (!exists("L8_masked_merged")) { # If there is no processing and the object isnt already loaded, load the pre-processed file
    L8_masked_merged <- stack("./L8_level2_complete.grd") 
  }
  

  # Check if the L7 file already exists or if force process is marked
  if(file.exists("./L7_level2_complete.grd") == FALSE | force_process){
    
    # Make sure that the ned object is in the same transformation as the L7 objects
    ned <- spTransform(ned, CRSobj = crs(L7[[1]]))
      
    
    for(i in 1:length(L7)){
      tmp <- L7[[i]] # Set the temp object
      tmp <- crop(tmp, ned) # Crop the bands to the dutch extent
      tmp[tmp < 1] <- NA  # set all values that arent needed to NA
      
      
      # NOTE: Due to errors in the cloudmask for level 2 data this part is not executed
      # This has been discussed with the student assistents
      if(level_2 && FALSE){
        # For level 2 data the cloudmask is more complicated
        # If both bit 6 and 8 are set this means there are clouds and thus we discard the pixel
        # Or if bit 3 is set, it means its a bad pixel and we also discard it
        cloud <- calc(tmp[[7]], fun = function(x) ((bitwAnd(x, 32) & bitwAnd(x, 128)) | (bitwAnd(x, 8))))
        # Set all other values to be transparent
        cloud[cloud == 1] <- NA 
        # Mask the tmp raster with the cloud mask
        L7[[i]] <- mask(x = tmp, mask = cloud)
      } else if (!level_2) {
        #Select all pixels with a cloud value over 700
        cloud <- tmp[[9]] > 700
        # Set all other values to be transparent 
        cloud[cloud == 1] <- NA 
        # Mask the tmp raster with the cloud mask
        L7[[i]] <- mask(x = tmp, mask = cloud)
      }
      
    }
    
    # Merge both L7 objects together
    # For more pictures put this part in a loop
    L7_all <- merge(L7[[1]], L7[[2]])
    # Make sure that the complete picture is masked by the netherlands shape
    L7_all_masked <- mask(x = L7_all, mask=ned)
  } else { # If there is no processing and the object isnt already loaded, load the pre-processed file
    L7_all_masked <- stack("./L7_level2_complete.grd")
  }
} else {
  # If no processing needs to be done, then load the files 
  L8_masked_merged <- stack("./L8_level2_complete.grd")
  L7_all_masked <- stack("./L7_level2_complete.grd")
}

#For both L7 and the L8 data create the NDVI layers using the appropiate layers (RED and NIR)
NDVI_L7 <- overlay(L7_all_masked[[4]], L7_all_masked[[3]], fun = function(x, y){return((x-y)/(x+y))})
NDVI_L8 <- overlay(L8_masked_merged[[5]], L8_masked_merged[[4]], fun = function(x, y){return((x-y)/(x+y))})

NDVI_L7 <- crop(NDVI_L7, NDVI_L8)
NDVI_L8 <- crop(NDVI_L8, NDVI_L7)

#Create the difference NDVI by subracting the L7 NDVI from the L8 NDVI
NDVI_delta <- overlay(NDVI_L8, NDVI_L7, fun = function(x, y) {return(x - y)})

# Ensure that even if data hasnt been processed the ned object is in the right transformation
ned <- spTransform(ned, CRSobj = crs(L7_all_masked))


# Since our L7 data is limited in scope, we will limit ourselves to the munincipalities in Noord-Holland and Zuid-Holland
n_z_holland <- ned[ned$NAME_1 == "Zuid-Holland" |  ned$NAME_1 == "Noord-Holland" , ]


#####################################   This is to get a general overview of the two provinces          #################################################

# Get all munincipalities in N and Z holland
munins <- n_z_holland$NAME_2

# Retrieve a dataframe containing the munincipalites and data regarding their NDVI
NDVI_df <- GetDfDelta(ned, munins, NDVI_delta)

# Assign the correct collumn names
# Note: it is normalised by dividing the sum of NDVI by the area (km^2)
colnames(NDVI_df) <- c("Munincipality_Name", "delta_ndvi_normalised", "sd_nvdi")

# Append the data we are interested in to the shapefile data
n_z_holland$de <- NDVI_df[, 2]
n_z_holland$sd <- NDVI_df[, 3]

# Retrieve all the data their id by rowname
n_z_holland@data$id = rownames(n_z_holland@data)
# Melt all polygon points into a points dataframe
n_z_holland.points = fortify(n_z_holland, region="id")
# Merge the points and data into a dataframe more suitable for plotting
n_z_holland.df = join(n_z_holland.points, n_z_holland@data, by="id")

# Plot the frame, and fill it with the difference by color
plot_de <- ggplot(n_z_holland.df, aes(long, lat, group=group, fill = de)) + 
  geom_polygon() + # add the polygons
  geom_path(col = "black") + # Add paths to more easily differentiate munincipalities
  coord_equal() + # Make the coords equal to prevent stretching
  scale_fill_gradient2(low = "red", mid = "white",  high = "Green", na.value = "grey50", guide = "colourbar", aesthetics = "fill") +
  ggtitle("Difference in NDVI per km^2") +
  labs(fill = "delta NDVI\nper km^2") + 
  north(n_z_holland.df) +
  scalebar(n_z_holland.df, dist = 20)  + 
  annotate("text", x = 578000, y = 5855000, label = "Source: USGS\nProjection: WGS84")

# Plot the frame, and fill it with the difference by color
plot_sd <- ggplot(n_z_holland.df, aes(long, lat, group=group, fill = sd)) + 
  geom_polygon() + # add the polygons
  geom_path(col = "black") + # Add paths to more easily differentiate munincipalities
  coord_equal() + # Make the coords equal to prevent stretching
  scale_fill_gradient2(low = "red", mid = "white",  high = "Green", na.value = "grey50", guide = "colourbar", aesthetics = "fill") +
  ggtitle("Standard Deviation of NDVI deltas") + 
  labs(fill = "SD") + 
  north(n_z_holland.df) +
  scalebar(n_z_holland.df, dist = 20) + 
  annotate("text", x = 578000, y = 5855000, label = "Source: USGS\nProjection: WGS84")

grid.arrange(plot_de, plot_sd, ncol=2)

############################################# Here we zoom in to one specific munincipality ###############
# Set a munincipality to plot
# By default the munincipality with a very high loss of forestation
#munin <- "Wieringen"

# Other interesting ones
#munin <- "Middelharnis"
munin <- "Noordwijkerhout"
#munin <- "Alphen aan den Rijn"

# Get the sub rasters of this munincipality
N_s <- GetSubNDVI(NDVI_7, NDVI_8, NDVI_delta, ned, munin)

# Assign the correct name so they show up in the plot
names(N_s) <- c(paste0(munin, "_2002"), paste0(munin, "_2018"), "Comparison")

# Create a color palette
colr <- colorRampPalette(brewer.pal(11, 'RdYlGn'))

# Plot the rasters using levelplot, using the same color scale
levelplot(N_s, 
          margin=FALSE,                       # suppress marginal graphics
          colorkey=list(
            space='bottom',                   # plot legend at bottom
            labels=list(at=-5:5, font=4),      # legend ticks and labels 
            title = "NDVI"
          ),    
          par.settings=list(
            axis.line=list(col='transparent') # suppress axes and legend outline
          ),
          scales=list(draw=FALSE),            # suppress axis labels
          col.regions=colr,                   # colour ramp
          at=seq(-1, 1, len=101))            # colour ramp breaks



############################################################ Here we gather statistics related to that specific munincipality ###########

# Get the differences in amount of foreststation
# Current (default) breaks are 0.3 and 0.7
df <- GetDifferenceCells(N_s)

par(mfrow = c(1, 2))
# Create a Pie char for these values
# Note: the main argument seems to be failing and 
slices <- df[2:4, 1]
lbls <- c("maintained", "reduced", "removed")
pct <- round(slices, 1)
lbls <- paste(lbls, pct) # add percents to labels 
lbls <- paste(lbls,"km^2") # ad % to labels 
pie(slices,labels = lbls, col=rev(rainbow(length(lbls))), main=paste0(munin, ", 2002: ", round(df[1, 1], 1), " km^2"))


# Pie Chart with Percentages
slices <- df[2:4, 2]
lbls <- c("maintained", "improved", "created")
pct <- round(slices, 1)
lbls <- paste(lbls, pct) # add percents to labels 
lbls <- paste(lbls,"km^2") # ad % to labels 
pie(slices,labels = lbls, col=rev(rainbow(length(lbls))), main=paste0(munin, ", 2018: ", round(df[1, 2], 1), " km^2"))

##########################################################