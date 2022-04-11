age_classes=read.csv("Hummingbird age classes.csv",h=T,sep=";")
young=age_classes$sp[which(age_classes$Quartile=="young")]
olds=age_classes$sp[which(age_classes$Quartile=="old")]

res=read.csv("data.csv",sep=";",h=T)

res$r=rep(NA,nrow(res)) ### richness of all hummingbirds
res$r_young=NA ### richness of young endemic species
res$r_old=NA ### richness of old endemic species

res$r_st_young=NA  ### Standardized richness of young endemic species
res$r_st_old=NA  ### Standardized richness of old endemic species

# loading a named list of 229 vector strings (i.e. one for each hummingbird species found in the Central and Northern Andes). Each vector string contains ID numbers of 0.25o grid cells occupied by each species 
load("dist_0.25degrees_Andes.RData")

#constructing community data matrix, used in the biogeographical null model
mat=matrix(0,ncol=length(dist_0.25degrees_Andes),nrow=nrow(res))
colnames(mat)=names(dist_0.25degrees_Andes)
rownames(mat)=res$ID
for(i in 1:nrow(mat)){
  ID=rownames(mat)[i]
  sp=vector()
  for(j in 1:length(dist_0.25degrees_Andes)){
    if(ID %in% dist_0.25degrees_Andes[[j]]){
      sp=append(sp,names(dist_0.25degrees_Andes)[j])
    }}
  mat[i,which(colnames(mat) %in% sp)]=1
  print(ID)
}

nrep=10 # Number of repititions, set to 1000 in the manuscript
for(i in 1:nrow(res)){
  sp=vector() # list of hummingbird species occuring in a focal grid cell
  ID=res$ID[i]
  for(j in 1:length(dist_0.25degrees_Andes)){
    if(ID %in% dist_0.25degrees_Andes[[j]]){
      sp=append(sp,names(dist_0.25degrees_Andes)[j])
    }
  } # end j
  
  r=length(sp)
  r_young=length(which(sp %in% young))
  r_old=length(which(sp %in% olds))
  
  
  #running the null model
  if(r_young>0|r_old>0){   
    r_old_null=rep(NA,nrep)
    r_young_null=rep(NA,nrep)
    for(k in 1:nrep){
      mat_sam=mat
      mat_sam=mat_sam[-i,]
      if(0 %in% rowSums(mat_sam)){mat_sam=mat_sam[-which(rowSums(mat_sam)==0),]}
      sp_null=vector()
      while(!length(sp_null)==r){
        loc=sample(1:nrow(mat_sam),1) # sample a random grid cell
        ssp=sample(colnames(mat_sam)[which(mat_sam[loc,]==1)],1) # sample a random species from this grid cell
        pps=which(colnames(mat_sam)==ssp) # remove the focal species from the community data matrix
        mat_sam=mat_sam[,-pps]
        # Add species to the null assemblage if not already sampled
        if(0 %in% rowSums(mat_sam)){mat_sam=mat_sam[-which(rowSums(mat_sam)==0),]} 
        sp_null=append(sp_null,ssp)
      } # end while
      
      
      
      r_old_null[k]=length(which(sp_null %in% olds))
      r_young_null[k]=length(which(sp_null %in% young))
      
    } # end k
    
    res$r_st_old[i]=length(which(r_old_null<r_old))/nrep
    res$r_st_young[i]=length(which(r_young_null<r_young))/nrep
    
  } #end if
  
  print(i)
}
