library(qiime2R)

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

mapme2<-host_tree$tip.label

mapme3<-mapme

# Simplify names by removing trailing numbers and certain prefixes
mapme3 <- gsub("([0-9]+)$", "", mapme3)
mapme3 <- gsub("^(A)(carolinensis|rudis|fulva|tennesseeinsis|picea)", "\\2", mapme3)
mapme3 <- gsub("_$", "", mapme3)

mapme3 <- gsub("DM_[0-9]+_", "", mapme3)
mapme3 <- gsub("_A", "_", mapme3)

# Standardize certain species names to uppercase
mapme3 <- gsub("picea$", "PICEA", mapme3, ignore.case = FALSE)
mapme3 <- gsub("rudis$", "RUDIS", mapme3, ignore.case = FALSE)
mapme3 <- gsub("carolinensis$", "CARO", mapme3, ignore.case = FALSE)
mapme3 <- gsub("fulva$", "FULVA", mapme3, ignore.case = FALSE)
mapme3 <- gsub("OCO_4ASP_tennesseeinsis", "OCO_4_ASP_2_TENN", mapme3, ignore.case = FALSE)








#Make phyloseq object from qza files
data<-qza_to_phyloseq(features="table.qza", tree="rooted-tree.qza","taxonomy.qza",metadata = "metadata.tsv")

#Remove uncharacterized phyla
data <- subset_taxa(data, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))

#Remove phyla that aren't bacterial
data <- subset_taxa(data, Phylum!="Chytridiomycota"& Phylum !="Ascomycota" & Phylum != "Chlorophyta" & Phylum != "Apicomplexa" & Phylum != "Ciliophora"  & Phylum != "Phragmoplastophyta" & Phylum != "Zoopagomycota")

#remove taxa that have chloroplast and mitochondria at level 5, 6,7 
data<-subset_taxa(data, Family != "Chloroplast")
data<-subset_taxa(data, Genus != "Mitochondria")

#Rarefy your data, here I am rarefying to 22533 sequences per sample
data_rarified = rarefy_even_depth(data, rngseed=1, sample.size=22533, replace=F, verbose = TRUE)

#Removing any taxa which aren't present following rarefying
data_rarified<- prune_taxa(taxa_sums(data_rarified) > 0, data_rarified) 

#Remove samples that are contols
data_rarified<-subset_samples(data_rarified,Type != "Negative")

#Drop tree tips that aren't in the microbiome dataset
#host_tree_preserve <- drop.tip(host_tree, 'A_umphreyi')

tips_to_drop <- setdiff(host_tree$tip.label, colnames(otu_table(data_rarified)))
host_tree <- drop.tip(host_tree, setdiff(host_tree$tip.label, colnames(otu_table(data_rarified))))
tree_tips<-host_tree$tip.label

#Drop microbiome samples that aren't in the Aphaenogaster tree
phyloseq_samples <- sample_names(data_rarified)
matched_samples <- intersect(tree_tips, phyloseq_samples)
data_rarified <- prune_samples(matched_samples, data_rarified)
data_rarified<- prune_taxa(taxa_sums(data_rarified) > 0, data_rarified) 

# Find tip labels for the Aphaenogaster tree
tree_tips <- host_tree$tip.label

# Find sample names from the phyloseq object
phyloseq_samples <- sample_names(data_rarified)

# Find the intersection of tree tips and phyloseq sample names
matched_tips <- intersect(tree_tips, phyloseq_samples)
if (length(matched_tips)!=length(tree_tips)){print('WARNING!')}


#Normalized Aphaenogaster tree
host_tree2<-host_tree
# Apply sqrt transformation, adding a small constant
host_tree2$edge.length <- sqrt(host_tree2$edge.length+1e-2)

tstrains<-c('a035ba5b60c48674e539656fb55b565b','81e7c84a2b9ab8698de9418f4dfdf41d','a5993b81dfe1d1904cdcade915df3116')
#sstrains<-c('dcfef6e645c7c69988f72add2b740629','28950999d0058ddadc3451b7887c960c','5a5c61e50854a84e22b697fd3b94d7f4','d916fbfd7b16ac03ff6630390df17f6a','b5053e2b00f7cb84c24b422f4dece287')
#estrains<-'958df2f535a4d002fc6732dbbfb6fb03'

Wstrains<-which(rownames(otu_table(data_rarified)) %in% tstrains)

Wvi1<-as.vector(otu_table(data_rarified)[Wstrains[1],])
Wvi2<-as.vector(otu_table(data_rarified)[Wstrains[2],])
Wvi3<-as.vector(otu_table(data_rarified)[Wstrains[3],])

Wv1<-Wvi1/max(Wvi1)
Wv2<-Wvi2/max(Wvi2)
Wv3<-Wvi3/max(Wvi3)

LU<-data.frame(Wv1,Wv2,Wv3)
rownames(LU)<-sample_names(data_rarified)
colnames(LU)<-c('W1','W2','W3')

threshold<-0.005

Wvt1<-Wvi1
Wvt2<-Wvi2
Wvt3<-Wvi3

Wvt1[Wvt1<threshold*colSums(otu_table(data_rarified))[1]]<-0
Wvt2[Wvt2<threshold*colSums(otu_table(data_rarified))[1]]<-0
Wvt3[Wvt3<threshold*colSums(otu_table(data_rarified))[1]]<-0

#Strains renumbered to match phylogeny representation

#Strain 1
Wvs1<-sign(Wvt1)
#Strain 3
Wvs2<-sign(Wvt2)
#Strain 2
Wvs3<-sign(Wvt3)

ccf<-Wvs1+10*Wvs2+100*Wvs3
ccf <- ccf[match(host_tree2$tip.label, row.names(LU))]

tipper<-host_tree$tip.label

all3<-tipper[which(ccf==111)]
just2<-tipper[which(ccf==10)]
none<-tipper[which(ccf==0)]
just12<-tipper[which(ccf==11)]
just23<-tipper[which(ccf==110)]
  
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

  clader<-rep('grey',length(tipper))
  clader[which(tipper %in% all3)]<-'red'
  clader[which(tipper %in% none)]<-'white'
  clader[which(tipper %in% just2)]<-'pink'
  clader[which(tipper %in% just12)]<-'violetred'
  clader[which(tipper %in% just23)]<-'magenta'
  clader[which(tipper=='A_umphreyi')]<-'black'



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
  }else if (cv=='pink'){
    colors[which(rownames(lonlat) == tip), 1:2] <- c("pink", "pink")
  }else if (cv=='magenta'){
    colors[which(rownames(lonlat) == tip), 1:2] <- c("magenta", "magenta")
  }else if (cv=='violetred'){
    colors[which(rownames(lonlat) == tip), 1:2] <- c("violetred", "violetred")
  }else if (cv=='white'){
    colors[which(rownames(lonlat) == tip), 1:2] <- c("white", "white")
  }else{
    colors[which(rownames(lonlat) == tip), 1:2] <- c("black", "black")
    
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
  npink <- length(which(cvv[site_indices] == 'pink'))  # Count purple tips
  nmagenta <- length(which(cvv[site_indices] == 'magenta'))  # Count red tips
  nvr <- length(which(cvv[site_indices] == 'violetred'))  # Count blue tips
  nwhite <- length(which(cvv[site_indices] == 'white'))  # Count blue tips
  nblack <- length(which(cvv[site_indices] == 'black'))  # Count blue tips
  
  if (tree_choice=='CO1'){
  # Ensure valid frequencies
  if (sum(nred, nblue, npink,nmagenta,nvr,nwhite,nblack) > 0) {
    x <- c(nred, nblue, npink,nmagenta,nvr,nwhite,nblack) / sum(nred, nblue, npink,nmagenta,nvr,nwhite,nblack)
  } else {
    x <- c(0, 0, 0,0,0,0,1)  # Default proportions if no tips exist
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
        bg = c("red", "blue", "pink","magenta",'violetred','white','black')[which(x > 0)],  # Use vector indexing for color
        add = TRUE
      )
    } else {
      # For multi-color pies, use floating.pie
      floating.pie(
        lon, lat, x, edges = 200, radius = 0.03, 
        col = c("red", "blue", "pink",'magenta','violetred','white','black'),
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














