
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

#Rarefy your data, here I am rarefying to 22533 sequences per sample
data_rarified = rarefy_even_depth(data, rngseed=1, sample.size=10000, replace=F, verbose = TRUE)

#Removing any taxa which aren't present following rarefying
data_rarified<- prune_taxa(taxa_sums(data_rarified) > 0, data_rarified) 

#Remove samples that are contols
data_rarified<-subset_samples(data_rarified,Type != "Negative")

host_tree<-read.tree("Full_Trees/COI_TREE_FULL/iqt.contree") ## tree
host_tree<-root(host_tree, outgroup = "A_umphreyi", resolve.root = TRUE)

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

db<-data_rarified
data_rarified<-prune_samples(sample_names(data_rarified) %in% blue,data_rarified)


#Spatial distance matrix
sd<-distm(data.frame(sample_data(data_rarified)$Long,sample_data(data_rarified)$Lat), fun=distHaversine)

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

if (corr=='pearson'){
  mantelv<-mantel(cd,sd,method='pearson')
  print(mantelv)
  mantel_correlog<-mantel.correlog(cd,sd,r.type='pearson',cutoff=FALSE)
  mantel_p<-c(mantel_p,mantelv$signif)
}else if (corr=='spearman'){
  mantelv<-mantel(cd,sd,method='spearman')
  print(mantelv)
  mantel_correlog<-mantel.correlog(cd,sd,r.type='spearman',cutoff=FALSE)
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
if (length(which(mantel_df$sig==100))>0){
mantel_df<-mantel_df[-which(mantel_df$sig==100),]
}
mantel_df$sig<-as.factor(mantel_df$sig)
color_palette <- c("white", "black")

myplots[[metric]]<-ggplot(mantel_df, aes(x = class.index, y = Mantel.cor, fill = sig, group=1)) +
  geom_line()+guides(fill = "none", color = "none", linetype = "none", shape = "none")+
  geom_point(size = 2, shape = 21) +  # Point plot with fill based on Petal.Width
  scale_fill_manual(values = color_palette) +  # Manual fill colors
  labs(title = "",
       x = 'Distance Class', y = "Mantel Correlation",
       fill = "Significance") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 1)+theme_classic()+theme(text = element_text(size=12))+theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1))

}

plot_layout <- rbind(c(1,2),c(3,4))

comm_plots <- list(myplots[[1]], myplots[[2]], myplots[[3]], myplots[[4]])
comm_grid_layout <- plot_grid(myplots[[1]], myplots[[2]], myplots[[3]], myplots[[4]], align = "v", nrow = 2, rel_heights = c(1/4, 1/4, 1/4, 1/4), rel_widths = c(1/4, 1/4, 1/4, 1/4))

plot(comm_grid_layout)







