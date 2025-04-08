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

hitters<-which(rowSums(otu_table(Wdata))>0.0001*sum(otu_table(data_rarified)))

#Rename to keep in the outgroup
tax_table(data_rarified)[which(rownames(tax_table(data_rarified))=='28950999d0058ddadc3451b7887c960c'),6]<-'Entomoplasma'

Wdata<-subset_taxa(data_rarified,Genus =='Entomoplasma')

Wtree<-phy_tree(Wdata)

Wtree$edge.length <- sqrt(Wtree$edge.length+1e-2)

Wtree<-root(Wtree, outgroup = "28950999d0058ddadc3451b7887c960c", resolve.root = FALSE)


f <- colorRamp(c("white","red"))

tipper<-Wtree$tip.label

summer<-rowSums(otu_table(Wdata))

tstrains<-c('958df2f535a4d002fc6732dbbfb6fb03')



load<-c()
tipname<-c()
for (k in 1:length(tipper)){
  temp<-which(names(summer)==tipper[k])
  load<-c(load,summer[temp])
  if (tipper[k]==tstrains[1]){
    tipname<-c(tipname,'ASV 1')
  }else{
    tipname<-c(tipname,'')
  }
}


shapeclade<-load
shapeclade[which(load>1000)]<-23
shapeclade[which(load<1000)]<-21
shapeclade[which(tipper=='28950999d0058ddadc3451b7887c960c')]<-19

load_old<-load

load<-log(load)/log(200000)

LU<-data.frame(load,load_old)
rownames(LU)<-tipper
colnames(LU)<-c('W','oW')


cols <- rgb(f(LU$W)/255)
cols <- cols[match(Wtree$tip.label, row.names(LU))]

Wtree$tip.label<-tipname

ggtree(Wtree,layout = 'circular')+ geom_tippoint(size=3,color='black',fill=cols,shape=as.numeric(shapeclade))+geom_tiplab(offset = 0.1)





