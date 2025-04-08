library('vegan')
library('phyloseq')
library('ape')
library('microshades')
library('speedyseq')
library('ggplot2')
library('cowplot')
library('forcats')
library('tidyverse')
library(qiime2R)
library(ggpubr)
library(betareg)
setwd("~/Desktop/Final Aphaenogaster Code/Data")





tree_choice<-'CO1'
set.seed(1)

#Read in metadata
meta<-read.csv('metadata.csv',header=TRUE)
bioclim_values<-read.csv("scaled_bioclim_values.csv")

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

tipper<-host_tree$tip.label
e1<-c()
e2<-c()
e3<-c()
e4<-c()
for (k in 1:length(tipper)){
  e1<-c(e1,bioclim_values$X.min_temp_soil[which(grepl(paste0(str_split(tipper[k],'_')[[1]][1],'_'),bioclim_values$X))[1]])
  e2<-c(e2,bioclim_values$soil_temp_range[which(grepl(paste0(str_split(tipper[k],'_')[[1]][1],'_'),bioclim_values$X))[1]])
  e3<-c(e3,bioclim_values$soil_moisture_summer[which(grepl(paste0(str_split(tipper[k],'_')[[1]][1],'_'),bioclim_values$X))[1]])
  e4<-c(e4,bioclim_values$sample_data.data_rarified..Elevation[which(grepl(paste0(str_split(tipper[k],'_')[[1]][1],'_'),bioclim_values$X))[1]])
}

edf<-data.frame(e1,e1+e2,e3,scale(e4))
host_tree2 <- drop.tip(host_tree2, host_tree2$tip.label[which(is.na(edf[,1])=='TRUE')])
host_tree <- drop.tip(host_tree, host_tree2$tip.label[which(is.na(edf[,1])=='TRUE')])
host_tree <- drop.tip(host_tree, 'A_umphreyi')
if (length(which(is.na(edf[,1])=='TRUE'))>0){
  edf<-edf[-which(is.na(edf[,1])=='TRUE'),]
}
sedf<-edf
tipper2<-tipper[-which(tipper=='FR_6_ASP_PICEA')]
rownames(sedf)<-tipper2


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
tips_to_drop <- setdiff(rownames(sedf), colnames(otu_table(data_rarified)))
sedf2<-sedf[which(!(rownames(sedf) %in% tips_to_drop)),]

#Drop microbiome samples that aren't in the Aphaenogaster tree
phyloseq_samples <- sample_names(data_rarified)
matched_samples <- intersect(rownames(sedf2), phyloseq_samples)
data_rarified <- prune_samples(matched_samples, data_rarified)
data_rarified<- prune_taxa(taxa_sums(data_rarified) > 0, data_rarified) 


colnames(sedf2)<-c('minT','maxT','moisture','elevation')

pca<-prcomp(sedf2)

environment<-pca$x[,1]


redder<-c('TM_7N_PICEA','TM_10C_CARO','TC_2B_RUDIS','CC_9C_RUDIS','CC_2D_CARO','AC_9B_CARO','TC_4E_FULVA','CC_2B_RUDIS','TM_2E_RUDIS','TC_1I_CARO','TM_7M_CARO','TM_9E_RUDIS','CC_6E_RUDIS','TC_6E_RUDIS','TC_9B_PICEA','OCO_6E_RUDIS','SD_4E_PICEA','OCO_7A_RUDIS','OCO_3ASP_RUDIS','CC_1D_CARO','OCO_7G_RUDIS','TM_9E_PICEA','OCO_4C_CARO','AC_7E_CARO','OCO_4_ASP_2_TENN','OCO_4ASP_FULVA','TM_7J_RUDIS','CC_2C_CARO','TC_ASP1_FULVA','AC_6ASP_FULVA','TG_4B_PICEA')
purpler<-c('EM_9D_PICEA','CL_4A_PICEA','TM_2D_PICEA','RC_2B_PICEA','EM_9B_PICEA','RC_5B_PICEA','OCO_5C_PICEA','FR_6_ASP_PICEA','BMM_4B_PICEA','BM_6ASP_CARO','BM_5E_PICEA')

clademe<-rep('blue',length(environment))

dwr<-tax_glom(data_rarified,taxrank='Genus')
wolbachia<-c()
site<-c()
for (k in 1:length(environment)){
  temp<-names(environment)[k]
  if (temp %in% redder){
    clademe[k]<-'red'
  }else if (temp %in% purpler){
    clademe[k]<-'purple'
  }
  hit<-which(sample_names(data_rarified)==temp)
  wolbachia<-c(wolbachia,otu_table(dwr)[which(tax_table(dwr)[,6]=='Wolbachia'),hit])
  site<-c(site,sample_data(dwr)$Site[hit])
}

df<-data.frame(clademe,wolbachia,environment,site)
colnames(df)<-c('Clade','wolbachia','Environment','Site')
df$wolbachia <- ifelse(df$wolbachia == 0, 0.00001, df$wolbachia)
df$wolbachia<-df$wolbachia/max(df$wolbachia)

dfred<-df[which(df$Clade=='red'),]
dfblue<-df[which(df$Clade=='blue'),]
dfpurple<-df[which(df$Clade=='purple'),]

df$Site<-factor(df$Site,levels=c('Cades Cove','Oconaluftee','Abram\'s Creek','Twin Creeks','Tremont','Trillium Gap','Brushy Mountain','Brushy Mountain Myrtle','Ramsey Cascades','Elkmont','Snake Den Ridge','Cataloochee','Purchase Knob','Albright Grove'))
df$Clade<-factor(df$Clade,levels=c('red','purple','blue'))


Q1 <- quantile(df$wolbachia, .25)
Q3 <- quantile(df$wolbachia, .75)
IQR <- IQR(df$wolbachia)

#subset data where points value is outside 1.5*IQR of Q1 and Q3
outliers <- subset(df, df$wolbachia<(Q1 - 1.5*IQR) | df$wolbachia>(Q3 + 1.5*IQR))

Q1 <- quantile(dfred$wolbachia, .25)
Q3 <- quantile(dfred$wolbachia, .75)
IQR <- IQR(dfred$wolbachia)

#subset data where points value is outside 1.5*IQR of Q1 and Q3
outliers <- subset(dfred, dfred$wolbachia<(Q1 - 1.5*IQR) | dfred$wolbachia>(Q3 + 1.5*IQR))

Q1 <- quantile(dfblue$wolbachia, .25)
Q3 <- quantile(dfblue$wolbachia, .75)
IQR <- IQR(dfblue$wolbachia)

#subset data where points value is outside 1.5*IQR of Q1 and Q3
outliers <- subset(dfblue, dfblue$wolbachia<(Q1 - 1.5*IQR) | dfblue$wolbachia>(Q3 + 1.5*IQR))


# Fit the beta regression model
model <- betareg(wolbachia ~ Environment, data = df)
predicted_data <- data.frame(Environment = seq(min(df$Environment), max(df$Environment), length.out = 100))
predicted_data$Wolbachia = predict(model, newdata = predicted_data, type = "response")



p <- ggplot(df, aes(x = Environment, y = wolbachia)) +
  geom_point(aes(fill=Clade), color='black',size = 3, shape=21) +scale_fill_manual(values=c('red','purple','blue'))+# Add points
  geom_line(data = predicted_data, aes(x = Environment, y = Wolbachia), color = "black", size = 1,linetype="dashed") +labs(x = "Environment", y = "Wolbachia Relative Abundance")+
  theme_linedraw() +
  theme(
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 22),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )




# Display the plot
print(p)

model <- betareg(wolbachia ~ Environment, data = dfred)

predicted_data <- data.frame(Environment = seq(min(dfred$Environment), max(dfred$Environment), length.out = 100))
predicted_data$Wolbachia = predict(model, newdata = predicted_data, type = "response")

p <- ggplot(dfred, aes(x = Environment, y = wolbachia)) +
  geom_point(aes(shape=Site,color=Site),size = 3, fill='red') + scale_shape_manual(values=c(21,22,23,24,8,25)) + scale_color_manual(values=c('black','black','black','black','red','black'))+# Add points
  geom_line(data = predicted_data, aes(x = Environment, y = Wolbachia), color = "black", size = 1,linetype="dashed") +labs(x = "Environment", y = "Wolbachia Relative Abundance")+
  theme_linedraw() +
  theme(
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 22),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )




# Display the plot
print(p)

#First without RC blue removed
model <- betareg(wolbachia ~ Environment, data = dfblue)

predicted_data <- data.frame(Environment = seq(min(dfblue$Environment), max(dfblue$Environment), length.out = 100))
predicted_data$Wolbachia = predict(model, newdata = predicted_data, type = "response")

p <- ggplot(dfblue, aes(x = Environment, y = wolbachia)) +
  geom_point(aes(shape=Site,color=Site),size = 3, fill='blue') +  scale_shape_manual(values=c(6,1,21,22,24,23,0,5,2,8,25))+scale_color_manual(values=c('blue','blue','black','black','black','black','blue','blue','blue','blue','black'))+# Add points
  geom_line(data = predicted_data, aes(x = Environment, y = Wolbachia), color = "black", size = 1,linetype="dashed") +labs(x = "Environment", y = "Wolbachia Relative Abundance")+
  theme_linedraw() +
  theme(
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 22),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )




# Display the plot
print(p)


dfblue2<-dfblue[-which(rownames(dfblue)=='RC_6B_PICEA'),]
model <- betareg(wolbachia ~ Environment, data = dfblue2)

predicted_data <- data.frame(Environment = seq(min(dfblue2$Environment), max(dfblue2$Environment), length.out = 100))
predicted_data$Wolbachia = predict(model, newdata = predicted_data, type = "response")

p <- ggplot(dfblue2, aes(x = Environment, y = wolbachia)) +
  geom_point(aes(shape=Site,color=Site),size = 3, fill='blue') +  scale_shape_manual(values=c(6,1,21,22,24,23,0,2,8,25))+scale_color_manual(values=c('blue','blue','black','black','black','blue','blue','blue','blue','black'))+# Add points
  geom_line(data = predicted_data, aes(x = Environment, y = Wolbachia), color = "black", size = 1,linetype="dashed") +labs(x = "Environment", y = "Wolbachia Relative Abundance")+
  theme_linedraw() +
  theme(
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 22),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )




# Display the plot
print(p)

model <- betareg(wolbachia ~ Environment, data = dfpurple)
