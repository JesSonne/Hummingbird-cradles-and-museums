require(raster)
require(maptools)

source("Scripts/additional functions.R")

#reading map files
bg_map=raster("data/Map files/background map.tif")
shape=readShapeSpatial("data/Map files/ne_50m_land.shp")
proj4string(shape)=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
shape=crop(shape,extent(bg_map))

#reading data on species habitat preferences
hab=read.csv("data/habitat_preferences.csv",h=T,sep=";")

#reading richness data
dat=read.csv("data/grid_cell_data.csv",h=T,sep=";")

#Empirical richness of endemic species (Fig. 1) ----
plot(bg_map,col="lightgray",border="lightgray",legend=F, axes=FALSE,box=FALSE)

plot_data=dat
coordinates(plot_data)=~lon+lat
plot_data@proj4string <- crs(shape)
pixels <- SpatialPixelsDataFrame(plot_data, tolerance = 0.0040003, plot_data@data)

richness_end=raster(pixels[,'r_end'])# richness of all endemic species 
richness_young = raster(pixels[,'r_young'])#richness of young endemic species
richness_old = raster(pixels[,'r_old']) #richness of old endemic species

#panel A
plot(bg_map,col="lightgray",border="lightgray",legend=F, axes=FALSE,box=FALSE)
plot(richness_end,col=colfunc(200),add=T,legend.args=list(text='Richness', side=2, font=1, line=0.5, cex=1.5),
     legend.width=1.2, legend.shrink=0.2,zlim=c(1,11))
plot(shape,add=T,border="gray40")


#panel B
plot(bg_map,col="lightgray",border="lightgray",legend=F, axes=FALSE,box=FALSE)
plot(richness_young,col=colfunc(200),add=T,legend.args=list(text='Richness', side=2, font=1, line=0.5, cex=1.5),
     legend.width=1.2, legend.shrink=0.2,zlim=c(1,11))
plot(shape,add=T,border="gray40")


#panel C
plot(bg_map,col="lightgray",border="lightgray",legend=F, axes=FALSE,box=FALSE)
plot(richness_old,col=colfunc(200),add=T,legend.args=list(text='Richness', side=2, font=1, line=0.5, cex=1.5),
     legend.width=1.2, legend.shrink=0.2,zlim=c(1,11))
plot(shape,add=T,border="gray40")


habitat_old=subset(hab[,3:7],hab$class=="old")
habitat_young=subset(hab[,3:7],hab$class=="young")

#constructing habitat by age classe frequency matrix (Supplementary table S3) ----
habitat_mat=rbind(colSums(habitat_young),colSums(habitat_old))

#Applying Fisher's exact test to examine if young and old endemic species have different habitat preferences
fisher.test(habitat_mat)

#plotting the habitat preferences of young and old endemic species
stacked_pie_plot(habitat_old)
stacked_pie_plot(habitat_young)



#visualizing the standardized richness of young and old endemic species using 2D colors (figure 2) ----
df=plot2d(data=dat,vec1="r_st_young", vec2="r_st_old", lon="lon", lat="lat",choose_cols = cols_2D,pix=0.02)    

#plotting results producing Figure 2
plot(bg_map,col="lightgray",border="lightgray",legend=F, axes=FALSE,box=FALSE)
plotRGB(df[[1]],r=1,g=2,b=3,colNA='NA',bgalpha=0,add=T)
plot(shape,add=T,border="gray40")


plot(bg_map,col="lightgray",border="lightgray",legend=F, axes=FALSE,box=FALSE)
plotRGB(df[[1]],r=1,g=2,b=3,colNA='NA',bgalpha=0,add=T)
plot(shape2,add=T,border="gray40")

#plotting legend of Figure 2
fig_2_legend=df[[2]]
plot(fig_2_legend$vec1,fig_2_legend$vec2,
     col=fig_2_legend$col,
     pch=15,
     xlab="standardized richness (young endemic species)",
     ylab="standardized richness (old endemic species)")


#chi square test ont the species' range overlap ----
dat$r_chi=rep(NA,nrow(dat))
dat$r_chi[which(dat$r_young>0|dat$r_old>0)]=1
dat$r_chi[which(dat$r_st_young>0.5&dat$r_st_old>0.5)]=2
dat$r_chi[which(dat$r_st_young<=0.5&dat$r_st_old>0.5)]=3
dat$r_chi[which(dat$r_st_young>0.5&dat$r_st_old<=0.5)]=4

geo_data=dat
coordinates(geo_data)=~lon+lat
geo_data@proj4string <- crs(shape)
pixels <- SpatialPixelsDataFrame(geo_data, tolerance = 0.000750188, geo_data@data)
raster_dat <- raster(pixels[,'r_chi'])
raster_dat=projectRaster(raster_dat,crs=crs(shape))

discrete_colors=c("#0A91A2","#F47C54","#74C9B9","#F9A953")

#hotspot of young and old endemic species. See supplementary fig. S1b
plot(bg_map,col="lightgray",border="lightgray",legend=F, axes=FALSE,box=FALSE)
plot(raster_dat,col=discrete_colors,add=T,legend.args=list(text='Richness', side=2, font=1, line=0.5, cex=1.5),
     legend.width=1.2, legend.shrink=0.2,zlim=c(1,4))
plot(shape,add=T,border="gray40")


#assembling table for chi square test (Supplementary table S1a)
observed_frequencies=table(dat$r_chi)

#YE - standardized richness of young endemic species
#OE - standardized richness of old endemic species 
names(observed_frequencies)=c("YE & OE <0.5","YE & OE >0.5","YE < 0.5 & OE >0.5","YE > 0.5 & OE <0.5")

n=length(which(dat$r_young>0 |dat$r_old>0 ))
z=length(which(dat$r_st_young>0.5))/n #proportion of grid cells with standardized richness of young endemic species exceeding 0.5 
q=length(which(dat$r_st_old>0.5))/n#proportion of grid cells with standardized richness of old endemic species exceeding 0.5

expected_frequencies=rep(NA,4)
expected_frequencies[1]=(1-z)*(1-q)*n #YE & OE <0.5
expected_frequencies[2]=z*q*n #YE & OE >0.5
expected_frequencies[3]=(1-z)*q*n #YE < 0.5 & OE >0.5
expected_frequencies[4]=z*(1-q)*n #YE > 0.5 & OE <0.5


#Chi square formula
Chi=sum((observed_frequencies-expected_frequencies)^2/expected_frequencies)

#find p-value for the Chi-Square test statistic
pchisq(q=Chi, df=3, lower.tail=FALSE)

#Richness-frequency distribution of young and old endemic species (figure 3)----

#richness frequency of young endemic species
sub=subset(dat,dat$r_st_young>0.5)
tab=table(sub$r_young)
dat_young=data.frame(Richness=as.numeric(names(tab)),freq=c(tab),cat=rep("young",length(tab)))
add=data.frame(Richness=5,freq=0,cat="young") # added line to indicate zero frequency
dat_young=rbind(dat_young,add)

#richness frequency of old endemic species
sub=subset(dat,dat$r_st_old>0.5)
tab=table(sub$r_old)
dat_old=data.frame(Richness=as.numeric(names(tab)),freq=c(tab),cat=rep("old",length(tab)))
dat_plot=rbind(dat_young,dat_old)
dat_plot$cat=factor(dat_plot$cat,c("young","old"))

ggplot(dat_plot, aes(x=Richness, y=freq, fill=cat)) + 
  geom_bar(stat="identity", position=position_dodge())+
  theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_fill_manual(values=c("cyan3","red"))+
  theme(legend.position="top")























