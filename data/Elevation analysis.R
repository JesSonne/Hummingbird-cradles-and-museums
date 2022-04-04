setwd("H:/Documents/GITHUB repositories/Hummingbird-cradles-and-museums/data")
elev=read.csv("species elevations_IOC10,1.csv",h=T,sep=";") # data frame listing species elevations
refs=read.csv("taxa reference_IOC10,1.csv",h=T,sep=";") # data frame matching species names in phylogeny and in the elevational dataset   

#Matrix listing terminal branch length of endemic species (25% range size quantile) over 1000 permutations; missing species in the phylogeny were added based on a random sample from all R-values  
load("TBL_random.RData")
dim(TBL_random) # 58 species x 1000 permutations

#Matrix listing terminal branch length of endemic species (25% range size quantile) over 1000 permutations; missing species in the phylogeny were added based on a randomly sampled R-value following an 'early burst model'  
load("TBL_early.RData")
dim(TBL_early) # 58 species x 1000 permutations

#Matrix listing terminal branch length of endemic species (25% range size quantile) over 1000 permutations; missing species in the phylogeny were added based on a randomly sampled R-value following an 'late burst model'  
load("TBL_late.RData")
dim(TBL_late) # 58 species x 1000 permutations

#Vector listing terminal branch length of endemic species (25% range size quantile); missing species in the phylogeny were added based on the median R-value  
load("TBL_median.RData")

