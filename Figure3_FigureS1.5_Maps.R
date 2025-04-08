library(qiime2R)
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



tree_choice<-'CO1'

meta<-read.csv('metadata.csv',header=TRUE)

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
for (k in 1:length(tipper)){
  if (tipper[k]!='A_umphreyi'){
    long<-c(long,meta$Long[which(grepl(paste0(str_split(tipper[k],'_')[[1]][1],'_'),meta$sample.id))[1]])
    lat<-c(lat,meta$Lat[which(grepl(paste0(str_split(tipper[k],'_')[[1]][1],'_'),meta$sample.id))[1]])
  }else{
    long<-c(long,0.988*min(meta$Long[1:75]))
    lat<-c(lat,mean(meta$Lat[1:75]))
  }
}


lonlat<-data.frame(lat,long)
colnames(lonlat)=c('latitude','longitude')
rownames(lonlat)<-tipper
lonlat




x<-phylo.to.map(host_tree2,lonlat, rotate = TRUE, plot=TRUE)


library(sf)
shapefile_path <- "GRSM_BOUNDARY_POLYGON.shp"  
shape <- st_read(shapefile_path)
shape <- st_transform(shape, crs = st_crs("+proj=longlat +datum=WGS84"))

# Simplify the geometry with a small tolerance
shape <- st_simplify(shape, dTolerance = 0.001)  # Adjust dTolerance as needed
# Filter for only the largest polygon or relevant geometries
shape <- shape[st_area(shape) == max(st_area(shape)), ]
#Extract Coordinates 
shape_coords <- st_coordinates(shape)

#plot(shape_coords,asp=1,pch=19,cex=0.7)

x$map$x <- c(shape_coords[, 1])  # Combine existing map x-coordinates with shapefile's x-coordinates
x$map$y <- c( shape_coords[, 2])  # Combine existing map y-coordinates with shapefile's y-coordinates

if (tree_choice=='CAD'){
  yellow<-clade.members(x=MRCA(host_tree2,'CL_5C_CARO', 'PK_6B_PICEA'),phy = host_tree2,tip.labels=TRUE)
  green<-clade.members(x=MRCA(host_tree2,'TM_2E_RUDIS', 'AG_5A_PICEA'),phy = host_tree2,tip.labels=TRUE)
  clader<-rep('grey',length(tipper))
  clader[which(tipper %in% yellow)]<-'yellow'
  clader[which(tipper %in% green)]<-'green'
  clader[which(tipper=='A_umphreyi')]<-'black'
}else{
  purple<-clade.members(x=MRCA(host_tree2,'EM_9D_PICEA', 'RC_5B_PICEA'),phy = host_tree2,tip.labels=TRUE)
  blue<-clade.members(x=MRCA(host_tree2,'BMM_5B_PICEA', 'TC_6A_PICEA'),phy = host_tree2,tip.labels=TRUE)
  red<-setdiff(host_tree2$tip.label,c('A_umphreyi',purple,blue))
  clader<-rep('blue',length(tipper))
  clader[which(tipper %in% red)]<-'red'
  clader[which(tipper %in% purple)]<-'purple'
  clader[which(tipper=='A_umphreyi')]<-'black'
}



# Determine the total number of tips and nodes in the tree
num_tips <- Ntip(host_tree)
num_nodes <- Nnode(host_tree) # Internal nodes
total_nodes <- num_tips + num_nodes

# Create a colors matrix for all nodes (tips + internal nodes)
#colors <- matrix("red", nrow = num_tips, ncol = 2)



# Assign red to the descendant tip names

cvv<-clader
colors<-matrix(nrow=length(cvv),ncol=2)
for (tip in host_tree2$tip.label) {
  temp<-which(rownames(lonlat)==tip)
  cv<-cvv[temp]
  print(temp)
  if (cv=='blue'){
    colors[which(rownames(lonlat) == tip), 1:2] <- c("blue", "blue")
  }else if (cv=='red'){
    colors[which(rownames(lonlat) == tip), 1:2] <- c("red", "red")
  }else if (cv=='purple'){
    colors[which(rownames(lonlat) == tip), 1:2] <- c("purple", "purple")
  }else if (cv=='green'){
    colors[which(rownames(lonlat) == tip), 1:2] <- c("green", "green")
  }else if (cv=='yellow'){
    colors[which(rownames(lonlat) == tip), 1:2] <- c("yellow", "yellow")
  }else if (cv=='grey'){
    colors[which(rownames(lonlat) == tip), 1:2] <- c("grey", "grey")
  }else{
    colors[which(rownames(lonlat) == tip), 1:2] <- c("white", "white")
    
  }
}

# Assign edge colors based on the `colors` matrix
edge_colors <- rep("red", nrow(host_tree$edge)) # Default to black
for (i in seq_len(nrow(host_tree$edge))) {
  child <- host_tree$edge[i, 2] # Child node
  if (child <= num_tips) { # Check if it's a tip
    edge_colors[i] <- colors[child, 1] # Use the color of the tip
  }
}

# Plot the tree with edge colors matching the tip colors
plot.phylo.to.map(x, 
                  type = "phylogram", 
                  colors = colors, 
                  edge.color = colors, # Match edge colors to tip colors
                  ftype = "off",
                  xlim = c(min(lonlat$longitude) - 0.15, max(lonlat$longitude) + 0.15),
                  ylim = c(min(lonlat$latitude) - 0.15, max(lonlat$latitude) + 0.15),
                  lty = "dashed", cex.points =2)


sites <- unique(lonlat[,1])

library(plotrix)
for (i in 1:length(sites)) {
  # Indices for tips in the current site
  site_indices <- which(lonlat[, 1] == sites[i])
  
  # Count red, blue, and purple tips at the current site
  nred <- length(which(cvv[site_indices] == 'red'))  # Count red tips
  nblue <- length(which(cvv[site_indices] == 'blue'))  # Count blue tips
  npurple <- length(which(cvv[site_indices] == 'purple'))  # Count purple tips
  nyellow <- length(which(cvv[site_indices] == 'yellow'))  # Count red tips
  ngreen <- length(which(cvv[site_indices] == 'green'))  # Count blue tips
  ngrey <- length(which(cvv[site_indices] == 'grey'))  # Count blue tips
  
  if (tree_choice=='CO1'){
  # Ensure valid frequencies
  if (sum(nred, nblue, npurple) > 0) {
    x <- c(nred, nblue, npurple) / sum(nred, nblue, npurple)
  } else {
    x <- c(1/3, 1/3, 1/3)  # Default proportions if no tips exist
  }
  }else{
    if (sum(ngreen, nyellow,ngrey) > 0) {
      x <- c(ngreen, nyellow,ngrey) / sum(ngreen, nyellow,ngrey)
    } else {
      x <- c(1/2, 1/2)  # Default proportions if no tips exist
    }
  }
  # Get site's latitude and longitude
  lat <- lonlat[site_indices[1], 1]  # Latitude, assuming column 2 is lat
  lon <- lonlat[site_indices[1], 2]  # Longitude, assuming column 3 is lon
  
  # Determine if the pie is single-colored
  is_single_color <- sum(x > 0) == 1
  
  if (tree_choice=='CO1'){
  # Plot the pie chart
  if (sum(x > 0) > 0) {  # Ensure the pie chart has valid values
    if (is_single_color) {
      # For single-color pies, draw a filled circle manually
      symbols(
        x = lon, y = lat, circles = 0.03, inches = FALSE,
        bg = c("red", "blue", "purple")[which(x > 0)],  # Use vector indexing for color
        add = TRUE
      )
    } else {
      # For multi-color pies, use floating.pie
      floating.pie(
        lon, lat, x, edges = 200, radius = 0.03, 
        col = c("red", "blue", "purple"),
        border = "black"  # Add outer border for multi-color pies
      )
    }
  }
  }else{
    if (sum(x > 0) > 0) {  # Ensure the pie chart has valid values
      if (is_single_color) {
        # For single-color pies, draw a filled circle manually
        symbols(
          x = lon, y = lat, circles = 0.03, inches = FALSE,
          bg = c("green", "yellow","grey")[which(x > 0)],  # Use vector indexing for color
          add = TRUE
        )
      } else {
        # For multi-color pies, use floating.pie
        floating.pie(
          lon, lat, x, edges = 200, radius = 0.03, 
          col = c("green", "yellow","grey"),
          border = "black"  # Add outer border for multi-color pies
        )
      }
    }
    
  }
}














