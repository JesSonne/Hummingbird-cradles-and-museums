############## historical climate stability
require(climateStability) # for data see Fordham, et al. 2017, Ecography
require(raster)
require(maptools)
require(RColorBrewer)

#all the data is in 2.5 x 2.5 degree scale
data(temperatureDeviation)

#calculating climate stability
t_stab=1/temperatureDeviation

#loading background map
shape=readShapeSpatial("data/Map files/ne_50m_land.shp")
proj4string(shape)=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
shape=crop(shape,extent(bg_map))

andes=readShapeSpatial("data/Map files/Mountains.shp")
proj4string(andes)=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

#color palette for temperature
reds=colorRampPalette(brewer.pal(9,"Reds"))(200)

SpP_ras <- rasterize(andes, t_stab, getCover=TRUE)
SpP_ras[SpP_ras==0] <- NA
Andes_clim=mask(t_stab,SpP_ras)
Andes_clim=trim(Andes_clim)


#aggregating the distribution data for young and old endemic species
source("Scripts/additional functions.R")
#reading richness data
dat=read.csv("data/grid_cell_data.csv",h=T,sep=";")

#constructing identity raster
id=Andes_clim
id[]=1:length(id[])

point=extract(id,dat[,3:2])
point=point[!duplicated(point)]

old=subset(dat[,3:2],dat$r_st_old>0.5)
point_old=extract(id,old)
point_old=point_old[!duplicated(point_old)]

old=subset(dat[,3:2],dat$r_st_old>0.5)
point_old=extract(id,old)
point_old=point_old[!duplicated(point_old)]

raster_dat=id
raster_dat[]=NA

sub=subset(dat[,3:2], dat$r_st_young>0.5)
point=extract(id,sub)
point=point[!duplicated(point)]
raster_dat[point]=2
dub=point

sub=subset(dat[,3:2],dat$r_st_old>0.5 )
point=extract(id,sub)
point=point[!duplicated(point)]
raster_dat[point]=3
dub=append(dub,point)

point=as.numeric(names(which(table(dub)==2)))
raster_dat[point]=1

cols=c("#FAF285", "#6FCCDD" ,"#EE2424")

pos1=which(raster_dat[]==1)
pos2=which(raster_dat[] %in% c(3,2))

clim_dat=data.frame(t_stab=Andes_clim[pos1])
clim_dat$cat="Both"
clim_dat$id=pos1

add=data.frame(t_stab=Andes_clim[pos2])
add$cat="one"
add$id=pos2
clim_dat=rbind(clim_dat,add)

#plotting
windows(12,6);par(mfrow=c(1,3))
plot(Andes_clim,col=reds,legend=F);plot(andes,add=T,border="darkgray",lwd=0.01);plot(shape,border="gray40",add=T)
plot(raster_dat,col=cols,legend=F);plot(andes,add=T,border="darkgray",lwd=0.01);plot(shape,border="gray40",add=T)
boxplot(clim_dat$t_stab~clim_dat$cat,xlab="",ylab="Temperature Stability")

kruskal.test(clim_dat$t_stab~clim_dat$cat)
