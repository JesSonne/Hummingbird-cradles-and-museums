#plotting the species habitat preferences as pie diagrams (Fig. 1)
source("Scripts/additional functions.R")

#reading data on species habitat preferences
hab=read.csv("data/habitat_preferences.csv",h=T,sep=";")

habitat_old=subset(hab[,3:7],hab$class=="old")
habitat_young=subset(hab[,3:7],hab$class=="young")

#constructing habitat by age classe frequency matrix (Supplementary table S3)
habitat_mat=rbind(colSums(habitat_young),colSums(habitat_old))

#Applying Fisher's exact test to examine if young and old endemic species have different habitat preferences
fisher.test(habitat_mat)

#plotting the habitat preferences of young and old endemic species
stacked_pie_plot(habitat_old)
stacked_pie_plot(habitat_young)
