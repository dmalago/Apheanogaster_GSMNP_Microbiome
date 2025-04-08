library('vegan')
library('phyloseq')
library('ape')
library('microshades')
library('speedyseq')
library('ggplot2')
library('cowplot')
library('forcats')
library('tidyverse')


#Make phyloseq object from qza files
data<-qza_to_phyloseq(features="table.qza", tree="rooted-tree.qza","taxonomy.qza",metadata = "metadata.tsv")

#Remove uncharacterized phyla
data <- subset_taxa(data, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))

#Remove phyla that aren't bacterial
data <- subset_taxa(data, Phylum!="Chytridiomycota"& Phylum !="Ascomycota" & Phylum != "Chlorophyta" & Phylum != "Apicomplexa" & Phylum != "Ciliophora"  & Phylum != "Phragmoplastophyta" & Phylum != "Zoopagomycota")

#remove taxa that have chloroplast and mitochondria at level 5, 6,7 
data<-subset_taxa(data, Family != "Chloroplast")
data<-subset_taxa(data, Genus != "Mitochondria")

host_tree<-read.tree("Full_Trees/CAD_TREE_FULL/iqt.contree") ## tree
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

yellow<-clade.members(x=MRCA(host_tree,'CL_5C_CARO', 'PK_6B_PICEA'),phy = host_tree,tip.labels=TRUE)
green<-clade.members(x=MRCA(host_tree,'TM_2E_RUDIS', 'AG_5A_PICEA'),phy = host_tree,tip.labels=TRUE)


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
data_rarified = rarefy_even_depth(data, rngseed=1, sample.size=22533, replace=F, verbose = TRUE)

#Removing any taxa which aren't present following rarefying
data_rarified<- prune_taxa(taxa_sums(data_rarified) > 0, data_rarified) 

#Remove samples that are contols
data_rarified<-subset_samples(data_rarified,Type != "Negative")

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



metadata<-sample_data(data_rarified)

temp_phylo<-phyloseq(OTU_biom,TAX,tree_file)
Mphylo<-temp_phylo

M2phylo<-phyloseq(otu_table(Mphylo),tax_table(Mphylo),sample_data(metadata),phy_tree(Mphylo))


Group<-sample_names(M2phylo)
Group[which(Group %in% yellow)]<-'Yellow'
Group[which(Group %in% green)]<-'Green'
Group[intersect(which(Group!='Green'),which(Group!='Yellow'))]<-'White'

rMphylo<-M2phylo
sample_data(rMphylo)$Group<-Group
rMGphylo<-tax_glom(rMphylo,taxrank = 'Genus')

Wreads<-otu_table(rMGphylo)[which(tax_table(rMGphylo)[,6]=='Wolbachia'),]

diversity_df<-data.frame(sample_data(rMGphylo)$Group,t(Wreads))
colnames(diversity_df)<-c('Group','Wolbachia')
diversity_df$Group<-factor(diversity_df$Group,levels=c('White','Green','Yellow'))

kw_wolbachia<-kruskal.test(Wolbachia ~ Group, data = diversity_df)
if (kw_wolbachia$p.value<0.05){
  pw_wolbachia<-pairwise.wilcox.test(diversity_df$Wolbachia, diversity_df$Group,p.adjust.method = "BH")
  print(pw_wolbachia)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Do you want to show p-values for non-significant pairwise comparisons? (yes or no)
p_value_ns<-'no'
pshow<-'num'

##########################wolbachia Violin Plot#########################################################################################################################################

#Define violin plot
p_wolbachia <- ggplot(diversity_df, aes(x=Group, y=Wolbachia)) + geom_violin(aes(fill=Group),scale='width') +labs(x ="Group", y = "Wolbachia")
#Choose the size of font for the axes titles and labels
p_wolbachia<-p_wolbachia+theme(axis.title = element_text(size = 20))+theme(axis.text = element_text(size = 13))
#Choose the size of font for the legend title and lables
p_wolbachia<-p_wolbachia+theme(legend.title = element_text(size = 20))+theme(legend.text = element_text(size = 15))
#Choose the violin colors for each group
p_wolbachia<-p_wolbachia+scale_fill_manual(values=c("white","green","yellow"))
#Add boxplots inside the violins
p_wolbachia<-p_wolbachia+geom_boxplot(aes(fill=Group),width=0.1)
#Add the p-value for the Kruskal-Wallis test somewhere on your figure (you may have to change the x and y positions of this label)

#If there are significant differences in wolbachia between groups, make a tibble of the pairwise p-values for plotting and add the brackets to the plot
if (kw_wolbachia$p.value<0.02){
  #Names of the groups you're comparing
  group1<-c()
  group2<-c()
  #p-value for the pairwise comparisons
  p.adj<-c()
  #locations of the p-value brackets
  ypos<-c()
  new_y<-23000
  y_step<-3000
  #For each pairwise comparison...
  for (k in 1:length(rownames(pw_wolbachia$p.value))){
    for (j in 1:k){
      if (rownames(pw_wolbachia$p.value)[k]!=colnames(pw_wolbachia$p.value)[j]){
        #If there is a significant difference or you want to also show non-significant differences...
        if (pw_wolbachia$p.value[k,j]<0.1 || p_value_ns=='yes'){
          #Add an entry to your tibble include the names of the two groups being compared and the p-value for the comparison
          group1<-c(group1,rownames(pw_wolbachia$p.value)[k])
          group2<-c(group2,colnames(pw_wolbachia$p.value)[j])
          p.adj<-round(c(p.adj,as.numeric(pw_wolbachia$p.value[k,j])),5)
          #Add the y-position of the bracket and bump the y-location of the bracket up one for next time
          ypos<-c(ypos,new_y)
          new_y<-new_y+y_step
        }
      }
    }
  }
  pstar<-c()
  for (k in 1:length(p.adj)){
    if (p.adj[k]<=0.001){
      pstar<-c(pstar,'***')
    }
    else if (p.adj[k]<=0.01){
      pstar<-c(pstar,'**')
    }
    else if (p.adj[k]<=0.05){
      pstar<-c(pstar,'*')
    }
    else if (p.adj[k]<=0.1){
      pstar<-c(pstar,'.')
    }
    else{
      pstar<-c(pstar,'ns')
    }
  }
  if (pshow=='star'){
    pdisplay<-"{pstar}"
  }
  else{
    pdisplay<-"p = {p.adj}"
  }
  #Create your tibble (what is needed for the stat_pvalue_manual function)
  stat.test_wolbachia<-as_tibble(data.frame(group1,group2,p.adj))
  #Add the pairwise comparisons to your plot
  p_wolbachia<-p_wolbachia+stat_pvalue_manual(stat.test_wolbachia,label=pdisplay,y.position=ypos,size=5)
}
#Make your plot
p_wolbachia 

