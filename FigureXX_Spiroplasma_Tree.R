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
tax_table(data_rarified)[which(rownames(tax_table(data_rarified))=='3bf936500f387d1d3714d0f4fa23b939'),6]<-'Spiroplasma'

Wdata<-subset_taxa(data_rarified,Genus =='Spiroplasma')

Wtree<-phy_tree(Wdata)

Wtree$edge.length <- sqrt(Wtree$edge.length+1e-2)

Wtree<-root(Wtree, outgroup = "3bf936500f387d1d3714d0f4fa23b939", resolve.root = FALSE)


f <- colorRamp(c("white","red"))

tipper<-Wtree$tip.label

summer<-rowSums(otu_table(Wdata))

tstrains<-c('dcfef6e645c7c69988f72add2b740629','28950999d0058ddadc3451b7887c960c','5a5c61e50854a84e22b697fd3b94d7f4','d916fbfd7b16ac03ff6630390df17f6a','b5053e2b00f7cb84c24b422f4dece287')



load<-c()
tipname<-c()
for (k in 1:length(tipper)){
  temp<-which(names(summer)==tipper[k])
  load<-c(load,summer[temp])
  if (tipper[k]==tstrains[1]){
    tipname<-c(tipname,'ASV 2')
  }else if (tipper[k]==tstrains[2]){
    tipname<-c(tipname,'ASV 1')
  }else if (tipper[k]==tstrains[3]){
    tipname<-c(tipname,'ASV 5')
  }else if (tipper[k]==tstrains[4]){
    tipname<-c(tipname,'ASV 4')
  }else if (tipper[k]==tstrains[5]){
    tipname<-c(tipname,'ASV 3')
  }else{
    tipname<-c(tipname,'')
  }
}


shapeclade<-load
shapeclade[which(load>1000)]<-23
shapeclade[which(load<1000)]<-21
shapeclade[which(tipper=='3bf936500f387d1d3714d0f4fa23b939')]<-19

load_old<-load

load<-log(load)/log(200000)

LU<-data.frame(load,load_old)
rownames(LU)<-tipper
colnames(LU)<-c('W','oW')


cols <- rgb(f(LU$W)/255)
cols <- cols[match(Wtree$tip.label, row.names(LU))]

Wtree$tip.label<-tipname

ggtree(Wtree,layout = 'circular')+ geom_tippoint(size=3,color='black',fill=cols,shape=as.numeric(shapeclade))+geom_tiplab(offset = 0.1)





