


#Make phyloseq object from qza files
data<-qza_to_phyloseq(features="table.qza", tree="rooted-tree.qza","taxonomy.qza",metadata = "metadata.tsv")

#Remove uncharacterized phyla
data <- subset_taxa(data, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))

#Remove phyla that aren't bacterial
data <- subset_taxa(data, Phylum!="Chytridiomycota"& Phylum !="Ascomycota" & Phylum != "Chlorophyta" & Phylum != "Apicomplexa" & Phylum != "Ciliophora"  & Phylum != "Phragmoplastophyta" & Phylum != "Zoopagomycota")

#remove taxa that have chloroplast and mitochondria at level 5, 6,7 
data<-subset_taxa(data, Family != "Chloroplast")
data<-subset_taxa(data, Genus != "Mitochondria")

#Remove samples that are contols
data_rarified<-subset_samples(data,Type != "Negative")



tdatar = rarefy_even_depth(data, rngseed=1, sample.size=20000, replace=F, verbose = TRUE)
tstrains<-c('a035ba5b60c48674e539656fb55b565b','81e7c84a2b9ab8698de9418f4dfdf41d','a5993b81dfe1d1904cdcade915df3116')
wvs<-which(tax_table(tdatar)[,6]=='Wolbachia')
Wbs<-which(colSums(otu_table(tdatar)[wvs,])>0.005*20000)

#wvs<-which(rownames(tax_table(data_rarified))==tstrains[3])
#Wbs<-which(colSums(otu_table(data_rarified)[wvs,])>0.005*20000)

data<-subset_taxa(data, Genus != "Wolbachia")
data<-subset_taxa(data, Genus != "Spiroplasma")
data<-subset_taxa(data, Genus != "Entomoplasma")
data<-subset_taxa(data, Genus != "Candidatus_Sulcia")


#Rarefy your data, here I am rarefying to 22533 sequences per sample
data_rarified = rarefy_even_depth(data, rngseed=1, sample.size=10000, replace=F, verbose = TRUE)

#Removing any taxa which aren't present following rarefying
data_rarified<- prune_taxa(taxa_sums(data_rarified) > 0, data_rarified) 




metric<-'unifrac'

if (metric=='jaccard'){
  #metric<-'jaccard'
  edist<-vegdist(sign(t(otu_table(data_rarified))),method='jaccard',upper=TRUE,diag=TRUE, binary=FALSE)
  #edist<-as.matrix(edist)
}else if (metric=='bray'){
  #metric<-'bray'
  edist<-vegdist((t(otu_table(data_rarified))),method='bray',upper=TRUE,diag=TRUE, binary=FALSE)
  #edist<-as.matrix(edist)
}else if (metric=='unifrac'){
  metric<-'unifrac'
  #edist <- as.matrix(UniFrac(data_rarified,weighted=FALSE))
  edist <- UniFrac(data_rarified,weighted=FALSE)
}else{
  metric<-'weighted unifrac'
  #edist <- as.matrix(UniFrac(data_rarified,weighted=TRUE))
  edist <- UniFrac(data_rarified,weighted=TRUE)
}

cladal<-rep('W-',length(sample_names(data_rarified)))
cladal[Wbs]<-'W+'

wtest<-adonis2(edist~cladal)
print(wtest)

pcoa_edist<-pcoa(edist)$vectors[,1:2]
types<-c('W-','W+')
colvec<-c('grey','red')
no_type<-cladal
no_type[cladal=='W-']<-1
no_type[cladal=='W+']<-2
no_type<-as.numeric(no_type)
type<-cladal

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#           Plot the PCA

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(pcoa_edist,type='n',cex.lab=1.5,xlab='PCoA1',ylab='PCoA2')
ordihull(pcoa_edist,groups=no_type,draw="polygon",col= colvec[1:length(types)], label=F,alpha=0.25)
points(pcoa_edist, display = "sites", pch=16, col = colvec[no_type])

#Find the centroid (based on mean) of each treatment group in the NMDS plot
xmeans_edist<-c()    #Value of the centroid along the first NMDS axis
ymeans_edist<-c()    #Value of the centroid along the second NMDS axis

#For each treatment group...
for (k in 1:length(types)){
  #Find the points corresponding to that group
  pts<-which(type==types[k])
  #Find the averages values of the points in that treatment group along each axis
  xmeans_edist<-c(xmeans_edist,mean(pcoa_edist[pts,1]))
  ymeans_edist<-c(ymeans_edist,mean(pcoa_edist[pts,2]))
}

#Plot the centroid for each group overtop of your NMDS scatterplot
for (k in 1:length(types)){
  points(xmeans_edist[k],ymeans_edist[k],col=colvec[k],pch=16,cex=2)
  points(xmeans_edist[k],ymeans_edist[k],col='black',pch=1,cex=2)
}

bd<-betadisper(as.dist(edist),as.factor(cladal))
dispersion_anova<-anova(bd)
print(dispersion_anova)


