# Load packages
library(raster)
library(ncdf4)

# Import data ### See part 1 of this guide ###
gebco <- raster("D:\\Analysis\\shared\\GEBCO_2022.nc")    


# Create extent (our map area)
pr.e <- extent(-29.7455, -26.499, 37.337, 39.7485) #xmin, xmax, ymin, ymax
# Create a crop of the bathymetric data
pr.gebco <- crop(gebco, pr.e)

rm(gebco) # to save memory 

plot(pr.gebco)


colbr <- function(x, b1=50, b2=50, r1=-2, r2=-2) {
  # Min/max values of the raster (x)
  mi <- cellStats(x, stat="min")-100
  ma <- cellStats(x, stat="max")+100
  # Create sequences, but only use unique numbers
  s1 <- unique(round(seq(mi, 0, 0-mi/b1),r1))
  s2 <- unique(round(seq(0, ma, ma/b2),r2))
  # Combine sequence for our break points, removing duplicate 0
  s3 <- c(s1, s2[-1])
  # Create a list with the outputs
  # [[1]] = length of the first sequence minus 1 (water)
  # [[2]] = length of the second sequence minus 1 (land)
  # [[3]] = The break points
  x <- list(length(s1)-1, length(s2)-1, s3)
}

pr.br <- colbr(pr.gebco)#, b1=1000, b2=10, r1=-1, r2=-2)
pr.br
pr = getData("GADM", country="PRT", level=0)

blue.col <- colorRampPalette(c("darkblue", "lightblue"))

ocean.pal <- colorRampPalette(
  c("#000000", "#000413", "#000728", "#002650", "#005E8C",
    "#0096C8", "#45BCBB", "#8AE2AE", "#BCF8B9", "#DBFBDC")
)

land.pal <- colorRampPalette(
  c("#467832", "#887438", "#B19D48", "#DBC758", "#FAE769",
    "#FAEB7E", "#FCED93", "#FCF1A7", "#FCF6C1", "#FDFAE0")
)

zbreaks <- seq(-8000, 8000, by=10)
cols <-c(ocean.pal(sum(zbreaks<=0)-1), land.pal(sum(zbreaks>0)))


plot(pr.gebco, col=c(blue.col(pr.br[[1]]), terrain.colors(pr.br[[2]])), breaks=pr.br[[3]])


#plot(pr, add=TRUE)
bpy.colors()


###Plot
#make palette



plot(pr.gebco, col=cols)





# now with ggplot 
library(rasterVis)

#CTD = read.csv("CTD_log.csv")
#TAGS =  read.csv("tags.csv")
#UAV = read.csv("flight log.csv")


nc_bathy = "C:\\Users\\Alec Burslem\\OneDrive - University of St Andrews\\Supervision meetings\\FW debrief\\GEBCO\\bathy4.nc"
bathy2 = readGEBCO.bathy(nc_bathy)


colours = c("#467832", "#887438", "#B19D48", "#DBC758", "#FAE769","#FAEB7E", "#FCED93", "#FCF1A7", "#FCF6C1", "#FDFAE0"),
breaks = pr.br[[3]]


bathy_effort = ggplot(bathy2, aes(x=x, y=y))+
  coord_quickmap() +
  geom_raster(aes(fill=z))+
  scale_fill_gradientn(breaks = pr.br[[3]], colours = c("#467832", "#887438", "#B19D48", "#DBC758", "#FAE769","#FAEB7E", "#FCED93", "#FCF1A7", "#FCF6C1", "#FDFAE0"))+
  labs(fill = "Topo. (m)") 


bathy_effort


###required packages
library(RNetCDF)
library(maps)
library(mapdata)

###Data
#data locations
#bathy_fname <- "galapagos_gebco_08_-92_-2_-88_2.nc" # from https://www.bodc.ac.uk/data/online_delivery/gebco/gebco_08_grid/
#coast_fname <- "galapagos_18563.dat" # from




