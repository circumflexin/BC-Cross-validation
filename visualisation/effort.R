# Load packages
library(marmap)
library(ggplot2)
library(scales)

data_dir = "D:\\Data\\BC_crossval\\Effort\\22"
plot_dir = "D:\\Analysis\\BC Cross validation\\outputs"

az = readGEBCO.bathy(paste(data_dir,"\\study-area_bathy.nc", sep = ""))
CTD = read.csv(paste(data_dir, "\\CTD_log.csv", sep = ""))
TAGS =  read.csv(paste(data_dir,"\\tags.csv", sep = ""))
UAV = read.csv(paste(data_dir,"\\flight log.csv", sep = ""))
  
names(UAV)[1] = "date"

UAV$date2 = as.Date(UAV$date, "%d/%m/%y")


rescale(c(min(az), -1, 0 , max(az)))

myscale = scale_fill_gradientn(colours=c("darkblue","white","darkgreen", "#eab676"), 
                     values = rescale(c(min(az), -0.1, 0 , max(az)))) 

bathy_effort = ggplot(az, aes(x=x, y=y))+
  coord_quickmap() +
  geom_raster(mapping = aes(fill=z)) +
  myscale +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  xlab("Lat.")+
  ylab("Long.")+
  labs(fill = "Elev. (m)")+
  geom_contour(aes(z=z),
               breaks=c(-300, -500, -1000, -1500, -2000), colour="grey", size=0.2, alpha = 1) +
  geom_contour(aes(z=z),
               breaks=c(0), colour="black", size=1, alpha = 1)+
  geom_point(data = CTD, aes(x = Long, y = Lat), shape = 21, colour = "black", size = 2, stroke = 1.5)+
  geom_point(data = TAGS, aes(x = Lon, y = Lat), size = 2, color = "yellow")+
  geom_point(data = UAV, aes(x = Lon, y = Lat), shape = 25, size = 2, colour = "#c72f17", stroke = 1, inherit.aes = FALSE)


ggsave(paste(plot_dir, "effort_22.tif", sep = "\\"), bathy_effort, height = 210, width = 297/2, units = "mm", device = "tiff", dpi = 300)


