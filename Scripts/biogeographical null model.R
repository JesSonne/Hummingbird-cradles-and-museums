#Matrix listing terminal branch length of endemic species (25% range size quartile) over 1000 permutations; missing species in the phylogeny were added based on a random sample from all R-values  
load("data/Terminal branch lengths/TBL_random.RData")
dim(TBL_random) # 58 species x 1000 permutations

young=c(
  "Aglaeactis castelnaudii",
  "Aglaeactis pamela",
  "Amazilia castaneiventris",
  "Chaetocercus astreans",
  "Chaetocercus berlepschi",
  "Coeligena helianthea",
  "Coeligena orina",
  "Heliangelus micraster",
  "Lepidopyga lilliae",
  "Metallura baroni",
  "Metallura eupogon",
  "Metallura odomae",
  "Metallura theresiae",
  "Oreotrochilus chimborazo",
  "Oreotrochilus cyanolaemus"
)


olds=c(
  "Chalcostigma heteropogon",
  "Coeligena prunellei",
  "Eriocnemis mirabilis",
  "Heliangelus mavors",
  "Heliangelus regalis",
  "Heliangelus strophianus",
  "Heliodoxa imperatrix",
  "Hylonympha macrocerca",
  "Loddigesia mirabilis",
  "Ocreatus addae",
  "Ocreatus peruanus",
  "Oreonympha nobilis",
  "Sternoclyta cyanopectus",
  "Urochroa bougueri",
  "Urochroa leucura",
)


#reading richness data
res=read.csv("grid_cell_data.csv",sep=";",h=T)

res$r_st_young=NA  ### Standardized richness of young endemic species
res$r_st_old=NA  ### Standardized richness of old endemic species


#loading grid cell data 
load("null_model_dat.RData")
total_richness=null_model_dat[[1]]
n_grids=null_model_dat[[2]]

#number of null model repetitions
nrep=1000


#Running the null model 
#mat is a binary species by grid cell data matrix
for(i in 1:nrow(res)){
#Species richness in a focal grid cell
r_young=res$r_young[i] 
r_old=res$r_old[i]
r=res$total_richness[i]
  
#running the null model
if(r_young>0|r_old>0){   
  r_old_null=rep(NA,nrep)
  r_young_null=rep(NA,nrep)
  for(k in 1:nrep){
    mat_sam=mat
    mat_sam=mat_sam[-i,]
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

























