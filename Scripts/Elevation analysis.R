require(ggplot2)
require(quantreg)
require(splines)
require(ggpubr)

setwd("H:/Documents/GITHUB repositories/Hummingbird-cradles-and-museums")
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


######################################## Elevational analysis, example with Fig. 4
Choose_TBL=TBL_random # select one of TBL_random, TBL_early and TBL_late


#Looping to generate the richness of young and old endemic species in each elevational band (j) for each of the 1000 permutations (i) 
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


#using quantile regression to compute the 95% percent confidence interval
n_splines=25
dat_old=data.frame(x=x,y=y1,group="Old")
quan = c(0.025,0.5,0.975) #compute median and the 95% confidence interval 
m1 = rq(y1 ~ ns(x,n_splines), data=dat_old, tau=quan)
xvals = seq(min(dat_old$x), max(dat_old$x), length.out=100)
rqs_old = data.frame(x=xvals, predict(m1, newdata=data.frame(x=xvals)))
names(rqs_old) = c("x", paste0("p",100*quan))

dat_new=data.frame(x=x,y=y2,group="New")
quan = c(0.025,0.5,0.975) #compute median and the 95% confidence interval 
m1 = rq(y2 ~ ns(x,n_splines), data=dat_new, tau=quan)
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


#Panel B: elevational variability in the endemic species' branch lengths
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






plot(res_mat_old[2,]~elv,type="l",col="red")
points(res_mat_young[2,]~elv,type="l",col="blue")


res_mat_old[i,j]=obs_old
res_mat_young[i,j]=obs_young
res_mat_elv[i,j]=elv[j]


































