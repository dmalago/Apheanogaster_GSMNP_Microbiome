library('vegan')
library('phyloseq')
library('ape')
library('microshades')
library('speedyseq')
library('ggplot2')
library('cowplot')
library('forcats')
library('tidyverse')

#%%%%%%%%%%%%%%%%%%%%%%   FUNCTIONS FOR MAKING NICELY SCALED LEGENDS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xendosymbionts<-'TRUE'

custom_legend2column <- function (mdf, cdf, group_level = "Phylum", subgroup_level = "Genus", x = "Sample",
                                  y = "Abundance", legend_key_size = 0.4, legend_text_size = 10)
{
  if (is.null(mdf[[group_level]])) {
    stop("mdf 'group_level' does not exist")
  }
  
  if (is.null(mdf[[subgroup_level]])) {
    stop("mdf 'subgroup_level' does not exist")
  }
  
  if (is.null(cdf$hex)) {
    stop("cdf 'hex' does not exist")
  }
  
  col_name_group <- paste0("Top_", group_level)
  col_name_subgroup <- paste0("Top_", subgroup_level)
  group_level_names <- unique(cdf[[col_name_group]])
  
  for (i in 1:length(group_level_names))
  {
    if( i == 1)
    {
      complete_legend <-individual_legend2(mdf, cdf, group_level_names[i], col_name_group, col_name_subgroup, legend_key_size = legend_key_size, legend_text_size = legend_text_size)
      tracer<-length(unique(mdf[which(mdf[,28]==group_level_names[i]),29]))+2
    }
    else
    {
      new_legend <-individual_legend2 (mdf, cdf, group_level_names[i], col_name_group, col_name_subgroup, legend_key_size = legend_key_size, legend_text_size =legend_text_size)
      
      complete_height <- tracer
      new_height <-length(unique(mdf[which(mdf[,28]==group_level_names[i]),29]))+2
      tracer<-tracer+new_height
      
      complete_legend <-plot_grid(complete_legend, new_legend, ncol=1, rel_heights = c(complete_height,new_height))
    }
  }
  complete_legend
}

individual_legend2 <- function (mdf,
                                cdf,
                                group_name,
                                col_name_group = "Top_Phylum",
                                col_name_subgroup = "Top_Genus",
                                x = "Sample",
                                y = " Abundance",
                                legend_key_size = 0.4,
                                legend_text_size = 10)
{
  select_mdf <- mdf %>% filter(!!sym(col_name_group) == group_name)
  select_cdf <- cdf %>% filter(!!sym(col_name_group) == group_name)
  select_plot <- ggplot(select_mdf,
                        aes_string(x = x, y = y, fill = col_name_subgroup, text = col_name_subgroup)) +
    geom_col( position="fill") +
    scale_fill_manual(name = group_name,
                      values = select_cdf$hex,
                      breaks = select_cdf[[col_name_subgroup]]) +
    theme(legend.justification = "left") +
    theme(legend.title = element_text(face = "bold")) +
    theme(legend.key.size = unit(legend_key_size, "lines"), text=element_text(size=legend_text_size))
  legend <- get_legend(select_plot)
}



#Make phyloseq object from qza files
data<-qza_to_phyloseq(features="table.qza", tree="rooted-tree.qza","taxonomy.qza",metadata = "metadata.tsv")

#Remove uncharacterized phyla
data <- subset_taxa(data, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))

#Remove phyla that aren't bacterial
data <- subset_taxa(data, Phylum!="Chytridiomycota"& Phylum !="Ascomycota" & Phylum != "Chlorophyta" & Phylum != "Apicomplexa" & Phylum != "Ciliophora"  & Phylum != "Phragmoplastophyta" & Phylum != "Zoopagomycota")

#remove taxa that have chloroplast and mitochondria at level 5, 6,7 
data<-subset_taxa(data, Family != "Chloroplast")
data<-subset_taxa(data, Genus != "Mitochondria")

if (xendosymbionts=='TRUE'){
  data<-subset_taxa(data, Genus != "Wolbachia")
  data<-subset_taxa(data, Genus != "Spiroplasma")
  data<-subset_taxa(data, Genus != "Entomoplasma")
  data<-subset_taxa(data, Genus != "Candidatus_Sulcia")
}

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


for (k in 1:length(TAX[,6])){
  if (is.na(TAX[k,6])){
    if (!(is.na(TAX[k,5]))){
      TAX[k,6]<-paste0('unclassified_',TAX[k,5])
    }else{
      if (!(is.na(TAX[k,4]))){
        TAX[k,6]<-paste0('unclassified_',TAX[k,4])
      }else{
        if (!(is.na(TAX[k,3]))){
          TAX[k,6]<-paste0('unclassified_',TAX[k,3])
        }
        else{
          TAX[k,6]<-paste0('unclassified_',TAX[k,2])
        }
      }
    }
  }
}
for (k in 1:length(TAX[,6])){
  if (TAX[k,6]=='NA'){
    if (TAX[k,5]!='NA'){
      TAX[k,6]<-paste0('unclassified_',TAX[k,5])
    }else{
      if (TAX[k,4]!='NA'){
        TAX[k,6]<-paste0('unclassified_',TAX[k,4])
      }else{
        if (TAX[k,3]!='NA'){
          TAX[k,6]<-paste0('unclassified_',TAX[k,3])
        }
        else{
          TAX[k,6]<-paste0('unclassified_',TAX[k,2])
        }
      }
    }
  }
}
for (k in 1:length(TAX[,6])){
  if (TAX[k,6]=='uncultured'){
    if (TAX[k,5]!='uncultured'){
      TAX[k,6]<-paste0('unclassified_',TAX[k,5])
    }else{
      if (TAX[k,4]!='uncultured'){
        TAX[k,6]<-paste0('unclassified_',TAX[k,4])
      }else{
        if (TAX[k,3]!='uncultured'){
          TAX[k,6]<-paste0('unclassified_',TAX[k,3])
        }
        else{
          TAX[k,6]<-paste0('unclassified_',TAX[k,2])
        }
      }
    }
  }
}
for (k in 1:length(TAX[,6])){
  if (TAX[k,6]==TAX[k,5]){
    TAX[k,6]<-paste0('unclassified_',TAX[k,6])
  }else{
    if (TAX[k,6]==TAX[k,4]){
      TAX[k,6]<-paste0('unclassified_',TAX[k,6])
    }else{
      if (TAX[k,6]==TAX[k,3]){
        TAX[k,6]<-paste0('unclassified_',TAX[k,6])
      }
    }
  }
}

metadata<-sample_data(data_rarified)

temp_phylo<-phyloseq(OTU_biom,TAX,tree_file)
Mphylo<-temp_phylo

M2phylo<-phyloseq(otu_table(Mphylo),tax_table(Mphylo),sample_data(metadata),phy_tree(Mphylo))

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

host_tree<-drop.tip(host_tree, 'A_umphreyi')

yellow<-clade.members(x=MRCA(host_tree,'CL_5C_CARO', 'PK_6B_PICEA'),phy = host_tree,tip.labels=TRUE)
green<-clade.members(x=MRCA(host_tree,'TM_2E_RUDIS', 'AG_5A_PICEA'),phy = host_tree,tip.labels=TRUE)

Group<-sample_names(M2phylo)
Group[which(Group %in% yellow)]<-'Yellow'
Group[which(Group %in% green)]<-'Green'
Group[intersect(which(Group!='Green'),which(Group!='Yellow'))]<-'Wh'

rMphylo<-M2phylo
sample_data(rMphylo)$Group<-Group

#%%%%%%%%%%%%%%%%%%%%%%   MAKE FANCY BAR GRAPHS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Calculate the total abundance of each phylum
summed_phyla <- rowSums(data.frame(otu_table(tax_glom(rMphylo, taxrank="Phylum"))))

# Get the names of the different phyla
phyla_names <- data.frame(tax_table(tax_glom(rMphylo, taxrank="Phylum")))[,2]

# Sort the phylum names based on their abundance
sorted_phyla_list <- phyla_names[order(-summed_phyla)]

# Prepare the OTU table for the plotting program, assuming this must be done at Genus level
mdf_prep <- prep_mdf(rMphylo)

# Create color assignments for the five most abundant phyla
color_objs_GP <- create_color_dfs(mdf_prep, selected_groups = sorted_phyla_list[1:5], cvd = TRUE)

# Extract the modified OTU table and color choices
mdf_GP <- color_objs_GP$mdf
cdf_GP <- color_objs_GP$cdf

# Assign colors; the colors are assigned to the top 5 sorted phyla
new_cdf_GP <- color_reassign(cdf_GP,
                             group_assignment = sorted_phyla_list[1:5],
                             color_assignment = c("micro_cvd_green", "micro_purple", "micro_orange", "micro_brown", "micro_cvd_blue"))



# Generate the legend
GP_legend_new <- custom_legend2column(mdf_GP, new_cdf_GP, legend_key_size=0.6, legend_text_size = 10)

# Define and format the plot
plot2 <- plot_microshades(mdf_GP, new_cdf_GP)
plot_diff2 <- plot2 + 
  scale_y_continuous(labels = scales::percent, expand = expansion(0)) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(size=6)) +
  facet_grid(~fct_relevel(Group, 'Wh','Green', 'Yellow'), scale="free_x", space = "free_x") +
  theme(strip.text.x = element_text(size=15)) +
  theme(axis.text.x = element_text(size=4)) +
  theme(plot.margin = margin(6, 20, 6, 6)) +
  theme(axis.title.x = element_text(size=20)) +
  theme(axis.title.y = element_text(size=20))

# Display the plot with the legend
plot_grid(plot_diff2, GP_legend_new, rel_widths = c(1, .25))

pdf(file = "Microshades_Barplot_W_Wolbachia.pdf", width = 10, height = 6)
plot_grid(plot_diff2, GP_legend_new, rel_widths = c(1, .25))
dev.off()

#%%%%%%%%%%%%%%%%%%%%%%   FIND PHYLUM PERCENTAGES ACROSS ALL SPECIES/POPULATIONS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Sort the read tallies for each phylum names so that the first is the most abundant, the second the second most abundant...
summed_phyla_order<-summed_phyla[order(-summed_phyla)]
#Turn the vector of sums into a dataframe
abundant_phyla<-data.frame(summed_phyla_order[1:5]*100/sum(summed_phyla))
#Name the rows of the dataframe according to their phylum names
rownames(abundant_phyla)<-sorted_phyla_list[1:5]



