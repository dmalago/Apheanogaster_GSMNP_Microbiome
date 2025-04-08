library(ape)
library(phytools)
library(stringr)
library(ggtree)
library(rstatix)
library(caper)
library(raster)
library(dplyr)

setwd("~/Desktop/Final Aphaenogaster Code/Data")



tree_choice<-'CO1'
set.seed(1)

#Read in metadata
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

host_tree<-drop.tip(host_tree, 'A_umphreyi')

#Normalized Aphaenogaster tree
host_tree2<-host_tree
host_tree
# Apply sqrt transformation, adding a small constant 
host_tree2$edge.length <- sqrt(host_tree2$edge.length+1e-2)

tipper<-host_tree2$tip.label

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



sedf<-data.frame(e1,e1+e2,e3,scale(e4))

colnames(sedf)<-c('minT','maxT','moisture','elevation')
rownames(sedf)<-rownames(bioclim_values)

pca<-prcomp(sedf)



environment<-pca$x[,1]
names(environment)<-host_tree2$tip.label

fit<-fastAnc(host_tree2,environment,vars=TRUE,CI=TRUE)

obj<-contMap(host_tree2,environment,plot=FALSE)

pdf(file = "HOST_ASR_ENV.pdf",   
    width = 5, # The width of the plot in inches
    height = 5) # The height of the plot in inches
plot(obj,legend=0.7*max(nodeHeights(host_tree2)),
     fsize=c(0.5,0.9))
dev.off()

#Make Screeplot of PCA

# Load necessary libraries
library(vegan)   # For RDA
library(factoextra)  # For visualization


# Visualize the eigenvalues using factoextra's fviz_eig function
fviz_eig(pca)

# Load necessary libraries
library(vegan)   # For RDA
library(ggplot2) # For plotting

# Assuming 'pca' is your prcomp object
# Extract loadings (also known as rotation in the prcomp output) for PC1
loadings_pc1 <- pca$rotation[, "PC1"]

# Create a data frame for plotting
loadings_df <- data.frame(
  terms = rownames(pca$rotation),
  value = loadings_pc1
)



# Plot the loadings using ggplot2 with enhanced aesthetics
ggplot(loadings_df, aes(x = value, y = reorder(terms, value))) +
  geom_col(fill = ifelse(loadings_df$value > 0, "#0A537D", "#b6dfe2")) +
  labs(
    title = "Envrionment: PC1 Loadings",  # Title of the plot
    x = "Loadings value",       # Label for the x-axis
    y = "Variable"              # Label for the y-axis
  ) +
  theme_linedraw() +  # Use linedraw theme
  theme(
    plot.title = element_text(size = 30, face = "bold", hjust = 0.5),  # Make title larger and bold
    axis.title = element_text(size = 25),  # Increase size of axis titles
    axis.text = element_text(size = 18),   # Increase size of axis texts
    legend.position = "none",              # Remove legend
    panel.grid.major = element_blank(),    # Remove major grid lines
    panel.grid.minor = element_blank()     # Remove minor grid lines
  )



# Assuming `data_rarified` is your phyloseq object and `environment` is your named vector
# Extract clade information from sample_data
clades <- sample_data(data_rarified)$Clade

# Ensure that the names of `clades` match the names of `environment`
# This is crucial for accurate matching
clades <- clades[match(names(environment), rownames(sample_data(data_rarified)))]

# Create a new data frame with environment values and their corresponding clade
env_df <- data.frame(
  Environment = environment,
  Clade = clades,
  row.names = names(environment)
)

# Check the structure of the new data frame
str(env_df)



clade_env_kruskal<-kruskal.test(Environment ~ Clade, data = env_df)
kw_pairwise<-pairwise.wilcox.test(env_df$Environment, env_df$Clade,p.adjust.method= "BH")
p_values <- as.data.frame(as.table(kw_pairwise$p.value))
names(p_values) <- c("group1", "group2", "p_value")
# Optionally, sort p_values alphabetically by Group1 then Group2 to maintain consistent ordering
p_values <- p_values[order(p_values$group1, p_values$group2),]

# Print again to verify order
print(p_values)

# Assume the base plot is correctly defined as before
a <- ggplot(env_df, aes(x = Clade, y = Environment, fill = Clade)) +
  geom_violin() +
  geom_point(position = position_jitter(seed = 1, width = 0.25), size = 3) +
  stat_summary(fun.y = mean, geom = "point", size = 7, color = "green") +
  scale_fill_manual(values = c("Blue" = "blue", "Red" = "red", "Purple" = "purple")) +
  labs(x = "Clade", y = "Environment") +
  theme_linedraw() +
  theme(
    axis.title = element_text(size = 25),
    axis.text = element_text(size = 25),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
base_y_pos <- max(env_df$Environment, na.rm = TRUE) * 1.01
p_values<-na.omit(p_values)

# Adding significance annotations with corrected indexing
for (i in 1:nrow(p_values)) {
  y_pos <- base_y_pos + (0.5 * i)  # Adjust y increment as needed
  a <- a + geom_signif(
    annotation = paste("p =", format.pval(p_values$p_value[i], digits = 2)),
    xmin = as.character(p_values$group1[i]),
    xmax = as.character(p_values$group2[i]),
    y_position = y_pos,
    tip_length = 0.02,
    textsize = 6,
    vjust = 0
  )
}

# Print the plot
print(a)



pdf(file = "Clade vs Environment.pdf", width = 8, height = 8)
print(a)
dev.off()


