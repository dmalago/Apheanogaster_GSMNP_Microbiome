
tree_choice<-'CO1'


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

#Read in host (Aphaenogaster) phylogeny
if (tree_choice == 'CO1'){
  host_tree<-read.tree("Full_Trees/COI_TREE_FULL/iqt.contree") ## tree
  #host_tree<-read.tree("COI_Sanger_contigs_aligned.fa.contree") ## tree
  host_tree<-root(host_tree, outgroup = "A_umphreyi", resolve.root = FALSE)
}else if (tree_choice == 'CAD'){
  host_tree<-read.tree("Full_Trees/CAD_TREE_FULL/iqt.contree") ## tree
  host_tree<-root(host_tree, outgroup = "A_umphreyi", resolve.root = FALSE)
}else{
  print('WARNING')
}

mapme<-host_tree$tip.label

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
host_tree
# Apply sqrt transformation, adding a small constant 
host_tree2$edge.length <- sqrt(host_tree2$edge.length+1e-2)

purple<-clade.members(x=MRCA(host_tree2,'EM_9D_PICEA', 'RC_5B_PICEA'),phy = host_tree2,tip.labels=TRUE)
blue<-clade.members(x=MRCA(host_tree2,'BMM_5B_PICEA', 'TC_6A_PICEA'),phy = host_tree2,tip.labels=TRUE)
red<-setdiff(host_tree2$tip.label,c('A_umphreyi',purple,blue))


#Make a new tree so that you can label internal nodes as NULL
host_tree3<-host_tree
host_tree3$node.label<-NULL

tstrains<-c('a035ba5b60c48674e539656fb55b565b','81e7c84a2b9ab8698de9418f4dfdf41d','a5993b81dfe1d1904cdcade915df3116')
#sstrains<-c('dcfef6e645c7c69988f72add2b740629','28950999d0058ddadc3451b7887c960c','5a5c61e50854a84e22b697fd3b94d7f4','d916fbfd7b16ac03ff6630390df17f6a','b5053e2b00f7cb84c24b422f4dece287')
#estrains<-'958df2f535a4d002fc6732dbbfb6fb03'

db<-data_rarified
data_rarified<-prune_samples(sample_names(data_rarified) %in% blue,data_rarified)

Wstrains<-which(rownames(otu_table(data_rarified)) %in% tstrains)

Wvi1<-as.vector(otu_table(data_rarified)[Wstrains[1],])
Wvi2<-as.vector(otu_table(data_rarified)[Wstrains[2],])
Wvi3<-as.vector(otu_table(data_rarified)[Wstrains[3],])

latlong<-data.frame(sample_data(data_rarified)$Long,sample_data(data_rarified)$Lat)
colnames(latlong)=c('longitude','latitude')


nb <- chooseCN(coordinates(latlong), type = 5, d1 = 0, d2 = 10000, plot.nb = FALSE)
distnb <- nbdists(nb, latlong)
fdist <- lapply(distnb, function(x) 1 / x^1)
lw <- nb2listw(nb, style = 'W', glist = fdist, zero.policy = TRUE)  

set.seed(1)
moran.mc(Wvi3, lw,alternative='two.sided',nsim=10000)   #spatial autocorrelation



