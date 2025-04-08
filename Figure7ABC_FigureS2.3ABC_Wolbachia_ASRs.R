
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

host_tree_full<-host_tree

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

#Make a new tree so that you can label internal nodes as NULL
host_tree3<-host_tree
host_tree3$node.label<-NULL

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

ccf<-Wvs1
ccf <- ccf[match(host_tree$tip.label, row.names(LU))]

fit<-ace(ccf, host_tree3, type = "discrete", model = "ER")

plotTree(host_tree2,fsize=0.6,ftype="i")
nodelabels(node=1:host_tree2$Nnode+Ntip(host_tree2), pie=fit$lik.anc,piecol=c('white','red'),cex=0.5)
tiplabels(pie=to.matrix(ccf,sort(unique(ccf))),piecol=c('white','red'),cex=0.3)
add.simmap.legend(colors=c('white','red'),prompt=FALSE,x=0.9*par()$usr[1],y=-max(nodeHeights(host_tree2)),fsize=0.8)



metric<-Wvs1+10*Wvs2+100*Wvs3

namer<-metric
namer[namer==111]<-'all'
namer[namer==10]<-'just2'
namer[namer==0]<-'none'
namer[namer==11]<-'just12'
namer[namer==110]<-'just23'

lat<-namer

lat[lat=='all']<-500
lat[lat=='just23']<-400
lat[lat=='just12']<-300
lat[lat=='just2']<-200
lat[lat=='none']<-100
lat<-as.numeric(lat)

pruned<-setdiff(host_tree_preserve$tip.label, colnames(otu_table(data_rarified)))

df<-data.frame(c(sample_names(data_rarified),pruned),c(lat,rep(600,length(pruned))),c(lat,rep(600,length(pruned))))
colnames(df)<-c('name','lat','long')

if (tree_choice=='CO1'){
loseit<-c('AC_5B_CARO',
          'BMM_4C_PICEA',
          'CC_2C_CARO',
          'CC_5D_PICEA',
          'CCD_CARO',
          'OCO_4ASP_FULVA',
          'OCO_6E_RUDIS',
          'PK_5D_CARO',
          'RC_6B_PICEA',
          'SD_6J_PICEA',
          'TC_6A_PICEA',
          'TC_7A_PICEA',
          'TC_9C_PICEA',
          'TG_4B_PICEA')
}

for (k in 1:length(df[,1])){
  temp<-which(mapme2==df[k,1])
  df[k,1]<-mapme3[temp]
}

df<-df[-which(df[,1] %in% loseit),]
df<-df[which(df[,1] %in% loseit),]
write.csv(df,'Wolbachstatus46.csv')


Wp<-sample_names(data_rarified)[which(Wvs1+Wvs2+Wvs3>0)]
hitter<-c('A_umphreyi',intersect(red,Wp))
red_temp<-drop.tip(host_tree_full, setdiff(host_tree_full$tip.label,hitter))
red_temp<-root(red_temp, outgroup = "A_umphreyi", resolve.root = TRUE)
red_temp<-drop.tip(red_temp, 'A_umphreyi')


