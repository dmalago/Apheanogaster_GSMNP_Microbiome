
tree_choice<-'CO1'
corr<-'spearman'

#Make phyloseq object from qza files
data<-qza_to_phyloseq(features="table.qza", tree="rooted-tree.qza","taxonomy.qza",metadata = "metadata.tsv")

#Remove uncharacterized phyla
data <- subset_taxa(data, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))

#Remove phyla that aren't bacterial
data <- subset_taxa(data, Phylum!="Chytridiomycota"& Phylum !="Ascomycota" & Phylum != "Chlorophyta" & Phylum != "Apicomplexa" & Phylum != "Ciliophora"  & Phylum != "Phragmoplastophyta" & Phylum != "Zoopagomycota")

#remove taxa that have chloroplast and mitochondria at level 5, 6,7 
data<-subset_taxa(data, Family != "Chloroplast")
data<-subset_taxa(data, Genus != "Mitochondria")
data<-subset_taxa(data, Genus != "Wolbachia")
data<-subset_taxa(data, Genus != "Spiroplasma")
data<-subset_taxa(data, Genus != "Entomoplasma")
data<-subset_taxa(data, Genus != "Candidatus_Sulcia")

if (tree_choice=='CO1'){
host_tree<-read.tree("Full_Trees/COI_TREE_FULL/iqt.contree") ## tree
host_tree<-root(host_tree, outgroup = "A_umphreyi", resolve.root = TRUE)
}else if (tree_choice=='CAD'){
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

purple<-clade.members(x=MRCA(host_tree,'EM_9D_PICEA', 'RC_5B_PICEA'),phy = host_tree,tip.labels=TRUE)
blue<-clade.members(x=MRCA(host_tree,'BMM_5B_PICEA', 'TC_6A_PICEA'),phy = host_tree,tip.labels=TRUE)
red<-setdiff(host_tree$tip.label,c('A_umphreyi',purple,blue))


#Drop tree tips that aren't in the microbiome dataset
tips_to_drop <- setdiff(host_tree$tip.label, colnames(otu_table(data)))
host_tree <- drop.tip(host_tree, setdiff(host_tree$tip.label, colnames(otu_table(data))))
tree_tips<-host_tree$tip.label

#Drop microbiome samples that aren't in the Aphaenogaster tree
phyloseq_samples <- sample_names(data)
matched_samples <- intersect(tree_tips, phyloseq_samples)
data <- prune_samples(matched_samples, data)
data<- prune_taxa(taxa_sums(data) > 0, data) 

# Find tip labels for the Aphaenogaster tree
tree_tips <- host_tree$tip.label

# Find sample names from the phyloseq object
phyloseq_samples <- sample_names(data)

# Find the intersection of tree tips and phyloseq sample names
matched_tips <- intersect(tree_tips, phyloseq_samples)
if (length(matched_tips)!=length(tree_tips)){print('WARNING!')}


#Rarefy your data, here I am rarefying to 22533 sequences per sample
data_rarified = rarefy_even_depth(data, rngseed=1, sample.size=10000, replace=F, verbose = TRUE)

#Removing any taxa which aren't present following rarefying
data_rarified<- prune_taxa(taxa_sums(data_rarified) > 0, data_rarified) 

#Remove samples that are contols
data_rarified<-subset_samples(data_rarified,Type != "Negative")


#Drop tree tips that aren't in the microbiome dataset
tips_to_drop <- setdiff(host_tree$tip.label, colnames(otu_table(data_rarified)))
host_tree <- drop.tip(host_tree, setdiff(host_tree$tip.label, colnames(otu_table(data_rarified))))
tree_tips<-host_tree$tip.label


#Find OTU table in biom file
OTU_biom<-otu_table(data_rarified)


#Read in the phylogenetic tree that you made using Qiime2 (from the sv.seqs.fna file given to you by Zymo)
tree_file<-phy_tree(data)

#Rename the tree tips to match the ASV table names
cut_names<-c()
for (k in 1:length(taxa_names(tree_file))){
  #Depending on the Qiime2 output, you may need to split the names with a space or with an underscore
  #cut_names<-c(cut_names,strsplit(taxa_names(tree_file)[k],' ')[[1]][1])
  cut_names<-c(cut_names,strsplit(taxa_names(tree_file)[k],'_')[[1]][1])
}
taxa_names(tree_file)<-cut_names

#Find the taxonomy table in the biom file (this exists for the ASV biom file!)
TAX<-tax_table(data)
colnames(TAX)<-c('Kingdom','Phylum','Class','Order','Family','Genus','Species')
TAX<-gsub('p__', '', TAX)
TAX<-gsub('c__', '', TAX)
TAX<-gsub('o__', '', TAX)
TAX<-gsub('f__', '', TAX)
TAX<-gsub('g__', '', TAX)
TAX<-gsub('s__', '', TAX)

#Normalized Aphaenogaster tree
host_tree2<-host_tree
host_tree
# Apply sqrt transformation, adding a small constant 
host_tree2$edge.length <- sqrt(host_tree2$edge.length+1e-2)

#Find spatial coordinates of each sample
lonlat<-data.frame(sample_data(data_rarified)$Lat,sample_data(data_rarified)$Long)
colnames(lonlat)=c('latitude','longitude')

#Phylogenetic distance matrix
pd<-cophenetic.phylo(host_tree2)
pnames<-rownames(pd)

#Community distance matrix

mantel_p<-c()
myplots <- vector("list",8)

for (metric in 1:4){
if (metric==1){
  cd<-as.matrix(vegdist(sign(t(otu_table(data_rarified))),metric='jaccard',diag=TRUE,upper=TRUE))
}else if (metric==2){
  cd<-as.matrix(vegdist(t(otu_table(data_rarified)),metric='bray',diag=TRUE,upper=TRUE))
}else if (metric==3){
  cd<-as.matrix(UniFrac(data_rarified,weighted=FALSE))
}else{
  cd<-as.matrix(UniFrac(data_rarified,weighted=TRUE))
}
dnames<-sample_names(data_rarified)


dm<-matrix(nrow=length(dnames),ncol=length(dnames))
for (k in 1:length(cd[1,])){
  for (j in 1:length(cd[1,])){
    r1<-which(dnames==pnames[[k]])
    r2<-which(dnames==pnames[[j]])
    dm[k,j]<-cd[r1,r2]
  }
}

if (corr=='pearson'){
  print('gotp')
  mantelv<-mantel(dm,pd,method='pearson')
  print(mantelv)
  mantel_correlog<-mantel.correlog(dm,pd,n.class=20,r.type='pearson')
  mantel_p<-c(mantel_p,mantelv$signif)
}else if (corr=='spearman'){
  print('gots')
  mantelv<-mantel(dm,pd,method='spearman')
  print(mantelv)
  mantel_correlog<-mantel.correlog(dm,pd,n.class=20,r.type='spearman')
  mantel_p<-c(mantel_p,mantelv$signif)
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
mantel_df<-mantel_df[-which(mantel_df$sig==100),]
mantel_df$sig<-as.factor(mantel_df$sig)
color_palette <- c("white", "black")

myplots[[metric]]<-ggplot(mantel_df, aes(x = class.index, y = Mantel.cor, fill = sig, group=1)) +
  geom_line()+guides(fill = "none", color = "none", linetype = "none", shape = "none")+
  geom_point(size = 2, shape = 21) +  # Point plot with fill based on Petal.Width
  scale_fill_manual(values = color_palette) +  # Manual fill colors
  labs(title = "",
       x = 'Host Genetic Distance Class', y = "Mantel Correlation",
       fill = "Significance") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 1)+theme_classic()+theme(text = element_text(size=12))+theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1))

}

plot_layout <- rbind(c(1,2),c(3,4))

comm_plots <- list(myplots[[1]], myplots[[2]], myplots[[3]], myplots[[4]])
comm_grid_layout <- plot_grid(myplots[[1]], myplots[[2]], myplots[[3]], myplots[[4]], align = "v", nrow = 2, rel_heights = c(1/4, 1/4, 1/4, 1/4), rel_widths = c(1/4, 1/4, 1/4, 1/4))

plot(comm_grid_layout)







