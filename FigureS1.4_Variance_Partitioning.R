library("ape")
library("stringr")
library("geosphere")
library("vegan")
library('caper')
library('phytools')
library('ape')

tree_choice<-'CO1'

meta<-read.csv('metadata.csv',header=TRUE)

data_path<-"./Soil_tifs"
files2 <- list.files(path = data_path, pattern = "*.tif", full.names = TRUE)
bioclim_layers2 <- stack(files2)

data_path2<-"./Soil_tifs2"
files3 <- list.files(path = data_path2, pattern = "*.tif", full.names = TRUE)
bioclim_layers3 <- stack(files3)

data_path3<-"./Soil_tifs3"
files4 <- list.files(path = data_path3, pattern = "*.tif", full.names = TRUE)
bioclim_layers4 <- stack(files4)


set.seed(1)


#Read in host (Aphaenogaster) phylogeny
if (tree_choice == 'CO1'){
  host_tree<-read.tree("Full_Trees/COI_TREE_FULL/iqt.contree") ## tree
  host_tree<-root(host_tree, outgroup = "A_umphreyi", resolve.root = TRUE)
}else if (tree_choice == 'CAD'){
  host_tree<-read.tree("Full_Trees/CAD_TREE_FULL/iqt.contree") ## tree
  host_tree<-root(host_tree, outgroup = "A_umphreyi", resolve.root = TRUE)
}else{
  print('WARNING')
}
host_tree<-drop.tip(host_tree,tip='A_umphreyi')

#Relabel the Aphaenogaster tree tips
tip_names <- host_tree$tip.label
simplified_tip_names <- gsub("([0-9]+)$", "", tip_names)
simplified_tip_names <- gsub("(A)(carolinensis|rudis|fulva|tennesseeinsis|picea)", "\\2", simplified_tip_names) # Remove the leading "A" from specific taxa, including picea
simplified_tip_names <- gsub("_$", "", simplified_tip_names)
simplified_tip_names <- gsub("picea", "PICEA", simplified_tip_names, ignore.case = FALSE)
simplified_tip_names <- gsub("rudis", "RUDIS", simplified_tip_names, ignore.case = FALSE)
simplified_tip_names <- gsub("carolinensis", "CARO", simplified_tip_names, ignore.case = FALSE)
simplified_tip_names <- gsub("fulva", "FULVA", simplified_tip_names, ignore.case = FALSE)
simplified_tip_names <- gsub("OCO_4ASP_tennesseeinsis", "OCO_4_ASP_2_TENN", simplified_tip_names, ignore.case = FALSE)
simplified_tip_names <- gsub("DM_53_PK7A_PICEA", "PK_7A_PICEA", simplified_tip_names, ignore.case = FALSE)
simplified_tip_names <- gsub("DM_49_AG7E_RUDIS", "AG_7E_RUDIS", simplified_tip_names, ignore.case = FALSE)
simplified_tip_names <- gsub("FR_6ASP_PICEA", "FR_6_ASP_PICEA", simplified_tip_names, ignore.case = FALSE)
simplified_tip_names <- gsub("DM_51_CCD_CARO", "CC_1D_CARO", simplified_tip_names, ignore.case = FALSE)
simplified_tip_names <- gsub("DM_54_TCI_CARO", "TC_1I_CARO", simplified_tip_names, ignore.case = FALSE)
simplified_tip_names <- gsub("DM_50_BMC_PICEA", "BM_1C_PICEA", simplified_tip_names, ignore.case = FALSE)
simplified_tip_names <- gsub("DM_68_TCI_CARO", "TC_1I_CARO", simplified_tip_names, ignore.case = FALSE)
simplified_tip_names <- gsub("DM_63_AG7E_RUDIS", "AG_7E_RUDIS", simplified_tip_names, ignore.case = FALSE)
simplified_tip_names <- gsub("DM_65_CCD_CARO", "CC_1D_CARO", simplified_tip_names, ignore.case = FALSE)
simplified_tip_names <- gsub("DM_64_BMC_PICEA", "BM_1C_PICEA", simplified_tip_names, ignore.case = FALSE)
simplified_tip_names <- gsub("DM_67_PK7A_PICEA", "PK_7A_PICEA", simplified_tip_names, ignore.case = FALSE)
simplified_tip_names <- gsub("DM_66_CCB_RUDIS", "CC_1B_RUDIS", simplified_tip_names, ignore.case = FALSE)
host_tree$tip.label <- simplified_tip_names



#Normalized Aphaenogaster tree
host_tree2<-host_tree
# Apply sqrt transformation, adding a small constant 
host_tree2$edge.length <- sqrt(host_tree2$edge.length+1e-2)

tipper<-host_tree$tip.label

long<-c()
lat<-c()
elevation<-c()
rownamer<-c()
for (k in 1:length(tipper)){
  long<-c(long,meta$Long[which(grepl(paste0(str_split(tipper[k],'_')[[1]][1],'_'),meta$sample.id))[1]])
  lat<-c(lat,meta$Lat[which(grepl(paste0(str_split(tipper[k],'_')[[1]][1],'_'),meta$sample.id))[1]])
  elevation<-c(elevation,meta$Elevation[which(grepl(paste0(str_split(tipper[k],'_')[[1]][1],'_'),meta$sample.id))[1]])
  rownamer<-c(rownamer,tipper[k])
}

lonlat<-data.frame(long,lat)
colnames(lonlat)=c('longitude','latitude')

values2 <- raster::extract(bioclim_layers2, lonlat)
bioclim_values2 <- as.data.frame(values2)

values3 <- raster::extract(bioclim_layers3, lonlat)
bioclim_values3 <- as.data.frame(values3)

values4 <- raster::extract(bioclim_layers4, lonlat)
bioclim_values4 <- as.data.frame(values4)

bioclim_values<-bind_cols(bioclim_values2, bioclim_values3, bioclim_values4,as.data.frame(elevation),)
rownames(bioclim_values)<-tipper
colnames(bioclim_values) <- c(" min_temp_soil", "soil_moisture_summer", "soil_temp_range","elevation")

# Update specific rows and columns by index
bioclim_values[which(grepl('OCO_',rownames(bioclim_values))==TRUE), 1] <- -5.193484
bioclim_values[which(grepl('OCO_',rownames(bioclim_values))==TRUE), 2] <-0.071459
bioclim_values[which(grepl('OCO_',rownames(bioclim_values))==TRUE), 3] <-32.10958
bioclim_values[which(grepl('PK_',rownames(bioclim_values))==TRUE),2]<- 0.160193 

bioclim_values<-scale(bioclim_values)
df_bioclim_values<-data.frame(bioclim_values)

e1<-c()
e2<-c()
e3<-c()
e4<-c()
rownamer2<-c()
for (k in 1:length(tipper)){
  e1<-c(e1,df_bioclim_values$X.min_temp_soil[which(grepl(paste0(str_split(tipper[k],'_')[[1]][1],'_'),rownames(df_bioclim_values)))[1]])
  e2<-c(e2,df_bioclim_values$soil_temp_range[which(grepl(paste0(str_split(tipper[k],'_')[[1]][1],'_'),rownames(df_bioclim_values)))[1]])
  e3<-c(e3,df_bioclim_values$soil_moisture_summer[which(grepl(paste0(str_split(tipper[k],'_')[[1]][1],'_'),rownames(df_bioclim_values)))[1]])
  e4<-c(e4,df_bioclim_values$elevation[which(grepl(paste0(str_split(tipper[k],'_')[[1]][1],'_'),rownames(df_bioclim_values)))[1]])
  rownamer2<-c(rownamer2,tipper[k])
}
edf<-data.frame(e1,e1+e2,e3,e4)
colnames(edf)<-c('min_soil_temp','max_soil_temp','soil_moisture','elevation')

#Phylogenetic distance matrix
pd<-cophenetic.phylo(host_tree2)

#Spatial distance matrix
sd<-distm(data.frame(long,lat), fun=distHaversine)

ed<-as.matrix(vegdist(edf,method='euclidean',diag=TRUE,upper=TRUE))

dbmem = dbmem(as.dist(sd))

Y<-cbind(dbmem,edf)

if (tree_choice=='CAD'){
  rda1<-dbrda(pd~min_soil_temp+max_soil_temp+soil_moisture+elevation+MEM1+MEM2+MEM3+MEM4+MEM5+MEM6,data=Y)  
  vif.cca(rda1)
  rda1<-dbrda(pd~min_soil_temp+max_soil_temp+elevation+MEM1+MEM2+MEM3+MEM4+MEM5+MEM6,data=Y)  
  vif.cca(rda1)
  rda1<-dbrda(pd~min_soil_temp+max_soil_temp+MEM1+MEM2+MEM3+MEM4+MEM5+MEM6,data=Y)  
  vif.cca(rda1)
  rda1<-dbrda(pd~min_soil_temp+MEM1+MEM2+MEM3+MEM4+MEM5+MEM6,data=Y)  
  vif.cca(rda1)
  vp<-varpart(pd,dbmem,edf$min_soil_temp)
  yellow<-clade.members(x=MRCA(host_tree2,'CL_5C_CARO', 'PK_6B_PICEA'),phy = host_tree2,tip.labels=TRUE)
  green<-clade.members(x=MRCA(host_tree2,'TM_2E_RUDIS', 'AG_5A_PICEA'),phy = host_tree2,tip.labels=TRUE)
  anova.cca(rda1)
  
  dev.off()
  op<-ordiplot(rda1,cex=1.5,display = 'sites',cex.lab=1.5,scaling=1,xlim=c(-0.5,0.5))
  text(rda1, display="bp", col='black',cex=1)
  points(op$sites[which(rownames(op$sites) %in% green),],cex=1.5,pch=16,col='green')
  points(op$sites[which(rownames(op$sites) %in% yellow),],cex=1.5,pch=16,col='yellow')
  points(op$sites,cex=1.5,col='black')
  
  }else{
  rda1<-dbrda(pd~min_soil_temp+max_soil_temp+soil_moisture+elevation+MEM1+MEM2+MEM3+MEM4+MEM5+MEM6+MEM7,data=Y)  
  vif.cca(rda1)
  rda1<-dbrda(pd~min_soil_temp+max_soil_temp+soil_moisture+MEM1+MEM2+MEM3+MEM4+MEM5+MEM6+MEM7,data=Y)  
  vif.cca(rda1)
  rda1<-dbrda(pd~min_soil_temp+max_soil_temp+MEM1+MEM2+MEM3+MEM4+MEM5+MEM6+MEM7,data=Y)  
  vif.cca(rda1)
  rda1<-dbrda(pd~min_soil_temp+MEM1+MEM2+MEM3+MEM4+MEM5+MEM6+MEM7,data=Y)  
  vif.cca(rda1)
  vp<-varpart(pd,dbmem,edf$min_soil_temp)
  purple<-clade.members(x=MRCA(host_tree2,'EM_9D_PICEA', 'RC_5B_PICEA'),phy = host_tree2,tip.labels=TRUE)
  blue<-clade.members(x=MRCA(host_tree2,'BMM_5B_PICEA', 'TC_6A_PICEA'),phy = host_tree2,tip.labels=TRUE)
  red<-setdiff(host_tree2$tip.label,c('A_umphreyi',purple,blue))
  anova.cca(rda1)
  
  dev.off()
  op<-ordiplot(rda1,cex=1.5,display = 'sites',cex.lab=1.5,scaling=2,xlim=c(-0.5,0.5))
  text(rda1, display="bp", col='black',cex=1)
  points(op$sites[which(rownames(op$sites) %in% red),],cex=1.5,pch=16,col='red')
  points(op$sites[which(rownames(op$sites) %in% purple),],cex=1.5,pch=16,col='purple')
  points(op$sites[which(rownames(op$sites) %in% blue),],cex=1.5,pch=16,col='blue')
  points(op$sites,cex=1.5,col='black')
  
}

par(mar=c(3,3,3,3))
plot(vp, bg = c(3, 7), Xnames=c('',''),cex=1.5)





