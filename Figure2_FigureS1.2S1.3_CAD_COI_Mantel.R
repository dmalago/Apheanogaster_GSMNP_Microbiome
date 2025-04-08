library("ape")
library("stringr")
library("geosphere")
library("vegan")
library("gridExtra")
library("grid")
library("cowplot")
library(dplyr)
library(ggplot2)

setwd("~/Desktop/Final Aphaenogaster Code/Data")

meta<-read.csv('metadata.csv',header=TRUE)

data_path<-"./Soil_tifs"
files2 <- list.files(path = data_path, pattern = "*.tif", full.names = TRUE)
bioclim_layers2 <- raster(files2)
bioclim_layers2<-stack(bioclim_layers2)

data_path2<-"./Soil_tifs2"
files3 <- list.files(path = data_path2, pattern = "*.tif", full.names = TRUE)
bioclim_layers3 <- raster(files3)
bioclim_layers3<-stack(bioclim_layers3)

data_path3<-"./Soil_tifs3"
files4 <- list.files(path = data_path3, pattern = "*.tif", full.names = TRUE)
bioclim_layers4 <- stack(files4)
bioclim_layers4<-stack(bioclim_layers4)

tree_choice_list<-c('CO1','CAD')
set.seed(1)

counter<-0
myplots <- vector("list",8)
pvalues<-c()
rstatistic<-c()

for (gene in 1:2){
  tree_choice<-tree_choice_list[gene]
  
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
  
  
  #Phylogenetic distance matrix
  pd<-cophenetic.phylo(host_tree2)
  
  #Spatial distance matrix
  sd<-distm(data.frame(long,lat), fun=distHaversine)
  
  ed<-as.matrix(vegdist(edf,method='euclidean',diag=TRUE,upper=TRUE))
  
  #Perform Mantel tests using spearman and pearson correlation
  print(paste('Tree choice:',tree_choice_list[gene]))
  Mtsp<-mantel(pd,sd,method='spearman')
  pvalues<-c(pvalues,Mtsp$signif)
  rstatistic<-c(rstatistic,Mtsp$statistic)
  Mtpr<-mantel(pd,sd,method='pearson')
  rstatistic<-c(rstatistic,Mtpr$statistic)
  pvalues<-c(pvalues,Mtpr$signif)
  eMtsp<-mantel(pd,ed,method='spearman')
  pvalues<-c(pvalues,eMtsp$signif)
  rstatistic<-c(rstatistic,eMtsp$statistic)
  eMtpr<-mantel(pd,ed,method='pearson')
  pvalues<-c(pvalues,eMtpr$signif)
  rstatistic<-c(rstatistic,eMtpr$statistic)

  temp<-mantel.correlog(pd,sd)
  frepeat<-c(0.000001,temp$break.pts)
  sd_mantel_correlog_sp<-mantel.correlog(D.eco=pd,D.geo=sd,r.type='spearman',cutoff=FALSE)
  sd_mantel_correlog_pr<-mantel.correlog(D.eco=pd,D.geo=sd,r.type='pearson',cutoff=FALSE)
  temp<-mantel.correlog(pd,ed)
  frepeat<-c(temp$break.pts)
  ed_mantel_correlog_sp<-mantel.correlog(D.eco=pd,D.geo=ed,r.type='spearman',cutoff=FALSE)
  ed_mantel_correlog_pr<-mantel.correlog(D.eco=pd,D.geo=ed,r.type='pearson',cutoff=FALSE)
  par(cex=1.5)
  plot(sd_mantel_correlog_sp)
  plot(ed_mantel_correlog_sp)
  
  for (metric in 1:4){
    counter<-counter+1
    if (metric==1){
      mantel_correlog<-sd_mantel_correlog_sp
      xname<-'Spatial Distance Class'
    }else if (metric==2){
      mantel_correlog<-sd_mantel_correlog_pr
      xname<-'Spatial Distance Class'
    }else if (metric==3){
      mantel_correlog<-ed_mantel_correlog_sp
      xname<-'Environmental Distance Class'
    }else{
      mantel_correlog<-ed_mantel_correlog_pr
      xname<-'Environmental Distance Class'
    }
    mantel_df<-as.data.frame(mantel_correlog$mantel.res)
    sig<-c()
    for (k in 1:length(mantel_correlog$mantel.res[,5])){
      if (is.na(mantel_correlog$mantel.res[k,5])){
        sig<-c(sig,100)
      }else{
        if (mantel_correlog$mantel.res[k,5]<0.05){
          sig<-c(sig,1)
        }else{
          sig<-c(sig,0)
        }
      }
    }
    
    mantel_df<-data.frame(mantel_df,sig)
    if (length(which(mantel_df$sig==100)>0)){
    mantel_df<-mantel_df[-which(mantel_df$sig==100),]
    }
    mantel_df$sig<-as.factor(mantel_df$sig)
    color_palette <- c("white", "black")
    
    myplots[[counter]]<-ggplot(mantel_df, aes(x = class.index, y = Mantel.cor, fill = sig, group=1)) +
      geom_line()+guides(fill = "none", color = "none", linetype = "none", shape = "none")+
      geom_point(size = 2, shape = 21) +  # Point plot with fill based on Petal.Width
      scale_fill_manual(values = color_palette) +  # Manual fill colors
      labs(title = "",
           x = xname, y = "Mantel Correlation",
           fill = "Significance") +
      geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 1)+theme_classic()+theme(text = element_text(size=12))+theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1))
    
  }
}

plot_layout <- rbind(c(1,2),c(3,4))

CO1_plots <- list(myplots[[1]], myplots[[2]], myplots[[3]], myplots[[4]])
CAD_plots <- list(myplots[[5]], myplots[[6]], myplots[[7]], myplots[[8]])
CO1CAD_plots<-list(myplots[[1]],myplots[[5]],myplots[[3]],myplots[[7]])

CO1_grid_layout <- plot_grid(myplots[[1]], myplots[[2]], myplots[[3]], myplots[[4]], align = "v", nrow = 2, rel_heights = c(1/4, 1/4, 1/4, 1/4), rel_widths = c(1/4, 1/4, 1/4, 1/4))
CAD_grid_layout <- plot_grid(myplots[[5]], myplots[[6]], myplots[[7]], myplots[[8]], align = "v", nrow = 2, rel_heights = c(1/4, 1/4, 1/4, 1/4), rel_widths = c(1/4, 1/4, 1/4, 1/4))
CO1CAD_grid_layout<-plot_grid(myplots[[1]],myplots[[5]],myplots[[3]],myplots[[7]])


plot(CO1_grid_layout)
plot(CAD_grid_layout)
plot(CO1CAD_grid_layout)



