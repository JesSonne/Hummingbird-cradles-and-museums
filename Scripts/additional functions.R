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





#function to plot pie diagrams
stacked_pie_plot=function(habitat_matrix){
  
  tot=nrow(habitat_matrix)
  hl_t=nrow(subset(habitat_matrix,habitat_matrix$humid.lowland==1& rowSums(habitat_matrix )>1))/tot*100
  hl_s=nrow(subset(habitat_matrix,habitat_matrix$humid.lowland==1& rowSums(habitat_matrix )==1))/tot*100
  cl_t=nrow(subset(habitat_matrix,habitat_matrix$cloud.forest==1& rowSums(habitat_matrix )>1))/tot*100
  cl_s=nrow(subset(habitat_matrix,habitat_matrix$cloud.forest==1& rowSums(habitat_matrix )==1))/tot*100
  tl_t=nrow(subset(habitat_matrix,habitat_matrix$humid.higland==1& rowSums(habitat_matrix )>1))/tot*100
  tl_s=nrow(subset(habitat_matrix,habitat_matrix$humid.higland==1& rowSums(habitat_matrix )==1))/tot*100
  Al_t=nrow(subset(habitat_matrix,habitat_matrix$paramo_puna==1& rowSums(habitat_matrix )>1))/tot*100
  Al_s=nrow(subset(habitat_matrix,habitat_matrix$paramo_puna==1& rowSums(habitat_matrix )==1))/tot*100
  Arm_t=nrow(subset(habitat_matrix,habitat_matrix$arid.shrubland==1& rowSums(habitat_matrix )>1))/tot*100
  Arm_s=nrow(subset(habitat_matrix,habitat_matrix$arid.shrubland==1& rowSums(habitat_matrix )==1))/tot*100
  
  
  
  
  habs1=data.frame(habs=rep(unique(colnames(habitat_matrix)),2),
                   share=c(hl_t,cl_t,tl_t,Al_t,Arm_t,hl_s,cl_s,tl_s,Al_s,Arm_s),cat=c(rep("t",5),rep("s",5)))
  
  habs1=habs1[order(habs1$share),]
  
  specific_data <- habs1[order(habs1$habs,habs1$cat), ]
  
  habs_data <-
    aggregate(habs1$share,
              by = list(habs = habs1$habs),
              FUN = sum)
  
  
  habs_data$habs_colors <-  brocolors("crayons")[c("Gray","Green Blue","Midnight Blue","Yellow Orange","Burnt Orange")]
  habs_data$labs=c("Paramo","Arid highlands","Cloud forest","Lowland rainforest","Treeline")
  habs_data=subset(habs_data,habs_data$x>0)
  
  
  
  # adjust these as desired (currently colors all specifics the same as habs)
  specific_data$specific_colors <-brocolors("crayons")[c("Gray","White","Green Blue","White","Midnight Blue","White","Yellow Orange","White","Burnt Orange","White")]
  specific_data$specific_colors1 <-brocolors("crayons")[c("White","White","White","White","White","White","White","White","White","White")]
  specific_data=subset(specific_data,specific_data$share>0)
  
  
  
  specific_data=subset(specific_data,specific_data$share>0)
  
  # format labels to display specific and % market share
  specific_labels <- paste(specific_data$specific, ": ", specific_data$share, "%", sep = "")
  
  # coordinates for the center of the chart
  center_x <- 0.5
  center_y <- 0.5
  
  windows(11,15)
  plot.new()
  
  
  title(main=paste0(names(sort(table(habitat_matrix$category)))[3]," n= ", tot ))
  # draw specific pie chart first
  specific_chart <-
    floating.pie(
      xpos = center_x,
      ypos = center_y,
      x = specific_data$share,
      radius = 0.25,
      border = "white",
      col = specific_data$specific_colors
    )
  habs_chart <-
    floating.pie(
      xpos = center_x,
      ypos = center_y,
      x = habs_data$x,
      radius = 0.20,
      border = "white",
      col = habs_data$habs_colors
    )
  
  
  specific_labels <- paste(specific_data$specific_colors, ": ", specific_data$share, "%", sep = "")
  
  
  pie.labels(
    x = center_x,
    y = center_y,
    angles = habs_chart,
    labels = habs_data$habs,
    radius = 0.125,
    bg = NULL,
    cex = 0.8,
    font = 2,
    col = "black")
}




