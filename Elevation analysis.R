require(ggplot2)
require(quantreg)
require(splines)
require(ggpubr)

source("additional functions.R") # additional functions used for plotting

elev=read.csv("data/species elevations_IOC10,1.csv",h=T,sep=";") # data frame listing species elevations
refs=read.csv("data/taxa reference_IOC10,1.csv",h=T,sep=";") # data frame matching species names in phylogeny and in the elevational dataset   
refs$phyloname=gsub("'","",refs$phyloname)

#Matrix listing the presence and absence of all Andean hummingbird species (n = 229) from sea level to 5000 metres. 
load("data/elv_mat.RData")
colnames(elv_mat)=elv=seq(0,5000,10)


#Matrix listing terminal branch length of endemic species (25% range size quantile) over 1000 permutations; missing species in the phylogeny were added based on a random sample from all R-values  
load("data/Terminal branch lengths/TBL_random.RData")
dim(TBL_random) # 58 species x 1000 permutations

#Matrix listing terminal branch length of endemic species (25% range size quantile) over 1000 permutations; missing species in the phylogeny were added based on a randomly sampled R-value following an 'early burst model'  
load("data/Terminal branch lengths/TBL_early.RData")
dim(TBL_early) # 58 species x 1000 permutations

#Matrix listing terminal branch length of endemic species (25% range size quantile) over 1000 permutations; missing species in the phylogeny were added based on a randomly sampled R-value following an 'late burst model'  
load("data/Terminal branch lengths/TBL_late.RData")
dim(TBL_late) # 58 species x 1000 permutations

#Vector listing terminal branch length of endemic species (25% range size quantile); missing species in the phylogeny were added based on the median R-value  
load("data/Terminal branch lengths/TBL_median.RData")


######################################## Elevational analysis, example with Fig. 4----
Choose_TBL=TBL_random # select one of TBL_random, TBL_early and TBL_late


#Looping to generate the richness of young and old endemic species in each elevational band (j) for each of the 1000 permutations (i) ----
n_rep=1000 
res_mat_young=matrix(0,ncol=ncol(elv_mat),nrow=n_rep)
res_mat_old=matrix(0,ncol=ncol(elv_mat),nrow=n_rep)
res_mat_elv=matrix(0,ncol=ncol(elv_mat),nrow=n_rep)
for(i in 1:n_rep){
  qua1=quantile(Choose_TBL[,i],0.25)
  qua2=quantile(Choose_TBL[,i],0.75)
  
  young=rownames(Choose_TBL)[which(Choose_TBL[,i]<=qua1)] #determining young endemic species
  old=rownames(Choose_TBL)[which(Choose_TBL[,i]>=qua2)] #determining old endemic species
  
  old1=refs$Andes_species[match(old,refs$phyloname)]
  young1=refs$Andes_species[match(young,refs$phyloname)]
  
  for(j in 1:ncol(elv_mat)){
    sps=as.character(rownames(elv_mat)[which(elv_mat[,j] ==1)])
    
    #Extracting the richness of young and old endemic species within each elevational band
    obs_old=length(which(sps %in% old1))
    obs_young=length(which(sps %in% young1))
    
    res_mat_old[i,j]=obs_old
    res_mat_young[i,j]=obs_young
    res_mat_elv[i,j]=elv[j]
  } #end j
  
  print(i)
} # end i


#ordering the results for plotting
x=vector();y1=vector();y2=vector()
for(i in 1:ncol(res_mat_old)){
  x=append(x, res_mat_elv[,i])
  y1=append(y1, res_mat_old[,i])
  y2=append(y2, res_mat_young[,i])
}


#using quantile regression to compute the 95% percent confidence interval ----
n_splines=25
dat_old=data.frame(x=x,y=y1,group="Old")
quan = c(0.025,0.5,0.975) #compute median and the 95% confidence interval 
m1 = rq(y ~ ns(x,n_splines), data=dat_old, tau=quan)
xvals = seq(min(dat_old$x), max(dat_old$x), length.out=100)
rqs_old = data.frame(x=xvals, predict(m1, newdata=data.frame(x=xvals)))
names(rqs_old) = c("x", paste0("p",100*quan))

dat_new=data.frame(x=x,y=y2,group="New")
quan = c(0.025,0.5,0.975) #compute median and the 95% confidence interval 
m1 = rq(y ~ ns(x,n_splines), data=dat_new, tau=quan)
xvals = seq(min(dat_new$x), max(dat_new$x), length.out=100)
rqs_new = data.frame(x=xvals, predict(m1, newdata=data.frame(x=xvals)))
names(rqs_new) = c("x", paste0("p",100*quan))

rqs_old[rqs_old < 0.1] <- 0
rqs_new[rqs_new < 0.1] <- 0

#plotting panel A
panel_A=ggplot() +
  geom_ribbon(data = rqs_old,aes(ymin = p2.5, ymax = p97.5,x=x),fill=t_col("red",50))+
  theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  geom_line(data = rqs_old,aes(x = x, y = p50),color="red",lwd=2) +
  
  geom_ribbon(data = rqs_new,aes(ymin =p2.5, ymax = p97.5,x=x),fill=t_col("cyan3",50))+
  theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  geom_line(data = rqs_new,aes(x = x, y = p50),color="cyan3",lwd=2) +ylim(0,15)+
  
  ylab("Species Richness")+xlab("")+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14))+xlim(0,5000)


#Panel B: elevational variability in the endemic species' branch lengths ----
TBL_median=apply(Choose_TBL,1,median)

data=data.frame(elev=NA,length=NA)
qua1=quantile(TBL_median,0.25)
qua2=quantile(TBL_median,0.75)

young=names(TBL_median)[which(TBL_median<=qua1)] #determining young endemic species
old=names(TBL_median)[which(TBL_median>=qua2)] #determining old endemic species

#Determining the terminal branch length of species in each elevational band
for(i in 1:ncol(elv_mat)){
  sps=as.character(rownames(elv_mat)[which(elv_mat[,i] ==1)])
  phylonam=refs$phyloname[which(refs$Andes_species %in% sps)]
  TBL_elevation=TBL_median[which(names(TBL_median) %in% phylonam)] #branch length of species in a focal elevational band
  
  sub_data=data.frame(elev=elv[i],length=TBL_elevation)
  data=rbind(data,sub_data)
  
  print(i)
}

bands=c(500,1000,1500,2000,2500,3000,3500,4000) #choosing bands to plot
plot_data=subset(data,data$elev %in% bands)
band_median=vector()

panel_B=ggplot(plot_data, aes(x=elev, y=length,group=elev)) + 
  geom_boxplot(fill="gray90")+
  geom_hline(yintercept=quantile(TBL_median,0.25),lty=2,lwd=2,col="cyan3")+
  geom_hline(yintercept=quantile(TBL_median,0.75),lty=2,lwd=2,col="red")+
  theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  xlab("Elevation (meters)") +ylab("Terminal branch length (MY)")+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14))+xlim(0,5000)+xlab("")

windows(18,25)
ggarrange(panel_A,panel_B,labels = c("A", "B"),widths = c(0.6,0.5),
          ncol = 1, nrow = 2)

#Comparing the richness of young and old endemic species to a mid-domain effect model


#function to compute the mid domain effect
#The function shuffles the species along the gradient while maintaining the coherency in the species elevational range limits

#shufling of a single species elevational range limits
shuf=function(string # vector of a species presence (1) and absence (0) along the elevational gradient
){
  tot=sum(string) #total number og elevational bands covered
  max_elv=length(string)-tot+1 #maximum elevation to sample for the species lower boundary
  min_elv=sample(1:max_elv,1) # sampling a random min elevation
  fills=min_elv:(min_elv+tot-1) #add the remaining presences
  
  #add absences
  out=rep(0,length(string))
  out[fills]=1
  return(out)
}

mid_domain=function(com # a species-elevation data matrix. species as rows, elevation bands as columns
){
  x1=t(com)
  #applying the shuf function to each species
  species_list <- split(x1, rep(1:ncol(x1), each = nrow(x1)))
  out=lapply(species_list,shuf)
  null_matrix=do.call("rbind",out)
  
  return(null_matrix)
}


#for the mid-domain effect analysis, we use the median branch length for each species
TBL_median=apply(TBL_random,1,median)

## running the mid domain analysis 1000 times
res_mat_young_mid=matrix(0,ncol=ncol(elv_mat),nrow=1000)
res_mat_old_mid=matrix(0,ncol=ncol(elv_mat),nrow=1000)
res_mat_elv_mid=matrix(0,ncol=ncol(elv_mat),nrow=1000)
for(i in 1:1000){
  
  qua1=quantile(TBL_median,0.25)
  qua2=quantile(TBL_median,0.75)
  
  young=names(TBL_median)[which(TBL_median<=qua1)]
  old=names(TBL_median)[which(TBL_median>=qua2)]
  
  old1=refs$Andes_species[match(old,refs$phyloname)]
  young1=refs$Andes_species[match(young,refs$phyloname)]
  
  elv_mat_null=mid_domain(elv_mat)
  rownames(elv_mat_null)=rownames(elv_mat)
  
  r_old=rep(NA,ncol(elv_mat))
  r_young=rep(NA,ncol(elv_mat))
  for(j in 1:ncol(elv_mat)){
    sps=as.character(rownames(elv_mat_null)[which(elv_mat_null[,j] ==1)])
    
    obs_old=length(which(sps %in% old1))
    obs_young=length(which(sps %in% young1))
    
    r_old[j]=obs_old
    r_young[j]=obs_young
    
  } #end j
  res_mat_old_mid[i,]=r_old
  res_mat_young_mid[i,]=r_young
  res_mat_elv_mid[i,]=elv
  print(i)
} 

x=vector();y=vector()
for(i in 1:ncol(res_mat_young)){
  x=append(x, res_mat_elv_mid[,i])
  y1=append(y, res_mat_old_mid[,i])
  y2=append(y, res_mat_young_mid[,i])
}

n_splines=7
dat_young_mid=data.frame(x=x,y=y2,group="young")
quan = c(0.025,0.5,0.975) #compute median and the 95% confidence interval 
m1 = rq(y ~ ns(x,n_splines), data=dat_young_mid, tau=quan)
xvals = seq(min(dat_young_mid$x), max(dat_young_mid$x), length.out=100)
rqs_young_mid = data.frame(x=xvals, predict(m1, newdata=data.frame(x=xvals)))
names(rqs_young_mid) = c("x", paste0("p",100*quan))


dat_old_mid=data.frame(x=x,y=y1,group="Old")
m1 = rq(y ~ ns(x,n_splines), data=dat_old_mid, tau=quan)
xvals = seq(min(dat_old_mid$x), max(dat_old_mid$x), length.out=100)
rqs_old_mid = data.frame(x=xvals, predict(m1, newdata=data.frame(x=xvals)))
names(rqs_old_mid) = c("x", paste0("p",100*quan))


#rqs_old_mid[rqs_old_mid < 0.1] <- 0
#rqs_new_mid[rqs_new_mid < 0.1] <- 0

p1=ggplot() +
  geom_ribbon(data = rqs_new_mid,aes(ymin = p5, ymax = p95,x=x),fill=t_col("gray",50))+
  theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  geom_line(data = rqs_new_mid,aes(x = x, y = p50),color="gray",lwd=2) +
  
  geom_line(data = rqs_new,aes(x = x, y = p50),color="cyan3",lwd=2) +ylim(0,10)+
  
  #geom_line(data = rqs_new,aes(x = x, y = p50),color="cyan3",lwd=2) +ylim(0,15)+
  
  ylab("Species Richness")+xlab("")+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14))+xlim(0,5000)

p2=ggplot() +
  geom_ribbon(data = rqs_old_mid,aes(ymin = p5, ymax = p95,x=x),fill=t_col("darkgray",50))+
  theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  geom_line(data = rqs_old_mid,aes(x = x, y = p50),color="darkgray",lwd=2) +ylim(0,10)+
  
  geom_line(data = rqs_old,aes(x = x, y = p50),color="red",lwd=2) +ylim(0,10)+
  
  #geom_line(data = rqs_new,aes(x = x, y = p50),color="cyan3",lwd=2) +ylim(0,15)+
  
  ylab("Species Richness")+xlab("")+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14))+xlim(0,5000)

windows(18,25)
ggarrange(p1,p2,labels = c("A", "B"),widths = c(0.6,0.5),
          ncol = 1, nrow = 2)






















