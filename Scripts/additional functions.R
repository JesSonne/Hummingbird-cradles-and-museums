#function constructing transpatent colors
t_col <- function(color, percent = 50, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color
  
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  
  ## Save the color
  invisible(t.col)
}



#specifying colors for the richness maps
colfunc<-colorRampPalette(rev(c("#D77B62","#E9A553","#F8D945","#CBDC60","#78C371","#55BE99","#4B96A5","#586698")))

#constructing color vector
n_col=200
vec1=colfunc(n_col)
col_vector=rep(NA,n_col)
for(i in 1:n_col){
  col_vector[i]=t_col(vec1[i])
}

#function to plot 2d colour maps
plot2d=function(data,vec1,vec2,lon,lat,shapefile=NULL,choose_cols=c("green", "yellow", "black", "blue"),pix=0.00080003){
  data=data[complete.cases(data),]
  data$lon=data[,lon]
  data$lat=data[,lat]
  
  data1=data
  coordinates(data)=~lon+lat
  
  
  if(!is.null(shapefile)){
    data@proj4string <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
    pixels <- SpatialPixelsDataFrame(data, tolerance = pix, data@data)
    
    raster <- raster(pixels[,vec1])
    r_A=projectRaster(raster,crs=crs(shapefile))
    raster <- raster(pixels[,vec2])
    r_B=projectRaster(raster,crs=crs(shapefile))
    rr3=brick(r_A,r_B)
    pos=is.na(rr3)
    rr3[is.na(rr3)] <- 0
    v=values(rr3)
    colors <- colormap::colors2d(v, choose_cols)
  }
  
  if(is.null(shapefile)){
    data@proj4string <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
    pixels <- SpatialPixelsDataFrame(data, tolerance = pix, data@data)
    r_A <- raster(pixels[,vec1])
    r_B <- raster(pixels[,vec2])
    rr3=brick(r_A,r_B)
    pos=is.na(rr3)
    rr3[is.na(rr3)] <- 0
    v=values(rr3)
    colors <- colormap::colors2d(v, choose_cols)}
  
  
  col <- rr3[[c(1,1,1)]]
  col[] <- t(col2rgb(colors))
  col[rowSums(v[])==0]=NA
  
  result=list()
  result[[1]]=col
  
  
  
  ss1=seq(min(data1[,vec1]),max(data1[,vec1]),length.out = 200)
  ss2=seq(min(data1[,vec2]),max(data1[,vec2]),length.out = 200)
  mat1=matrix(NA,200,200)
  mat2=matrix(NA,200,200)
  for(i in 1:200){
    for(j in 1:200){
      mat1[i,j]=i
      mat2[i,j]=j
    }}
  
  col_lab=colormap::colors2d(cbind(c(mat1),c(mat2)), choose_cols)
  label=data.frame(vec1=ss1[c(mat1)],vec2=ss2[c(mat2)],col=as.character(col_lab))
  
  result[[2]]=label
  
  
  
  
  return(result)
}

cols_2D=c( rgb(250,242,133, maxColorValue=255,alpha = 0.5),     
       rgb(111,204,221, maxColorValue=255,alpha=0.5),        
       rgb(43,44,119, maxColorValue=255,alpha=0.5),
       rgb(238,36,36, maxColorValue=255,alpha = 0.5))      








