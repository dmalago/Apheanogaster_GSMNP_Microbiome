library(ape)
library(phytools)
library(stringr)
library(ggtree)
library(rstatix)
library(caper)

map<-'blueto'

meta<-read.csv('metadata.csv',header=TRUE)

host_treeCO1<-read.tree("Full_Trees/COI_TREE_FULL/iqt.contree") ## tree
host_treeCO1<-root(host_treeCO1, outgroup = "A_umphreyi", resolve.root = FALSE)


host_treeCAD<-read.tree("Full_Trees/CAD_TREE_FULL/iqt.contree") ## tree
host_treeCAD<-root(host_treeCAD, outgroup = "A_umphreyi", resolve.root = FALSE)

for (k in 1:2){
  if (k==1){host_tree<-host_treeCO1}else{host_tree<-host_treeCAD}

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
  
  
  
  #Normalized Aphaenogaster tree
  host_tree2<-host_tree
  # Apply sqrt transformation, adding a small constant 
  host_tree2$edge.length <- sqrt(host_tree2$edge.length+1e-2)
  

if (k==1){
  host_treeCO11<-host_tree
  host_treeCO12<-host_tree2
}else{
  host_treeCAD1<-host_tree
  host_treeCAD2<-host_tree2
}
}


cophylo1<-cophylo(host_treeCO11,host_treeCAD1)
cophylo2<-cophylo(host_treeCO12,host_treeCAD2)

col_list<-cophylo2$assoc[,1]
col_list[which(cophylo2$assoc[,1] %in% 'A_umphreyi')]<-'black'

f <- colorRamp(c("purple","lavender","peachpuff","darkorange"))

col_long<-c()
for (k in 1:length(col_list)){
  col_long<-c(col_long,meta$Long[which(grepl(paste0(str_split(col_list[k],'_')[[1]][1],'_'),meta$sample.id))[1]])
}


col_long[which(cophylo2$assoc[,1] %in% 'A_umphreyi')]<-col_long[1]

long<-(col_long-min(col_long))/max(col_long-min(col_long))

col_list[which(cophylo2$assoc[,1] %in% 'A_umphreyi')]<-'black'
for (k in 1:length(col_list)){
  if (k!=which(cophylo2$assoc[,1] %in% 'A_umphreyi')){
    col_list[k]<-rgb(f(long[k])/255)
  }
}

if (map=='blueto'){
  clade_list1<-cophylo2$assoc[,1]
  clade_list2<-cophylo2$assoc[,2]
  col_clade<-c()
  for (k in 1:length(col_list)){
    if (clade_list1[k] %in% blue && clade_list2[k] %in% yellow){
      col_clade<-c(col_clade,'yellow')}
    else if (clade_list1[k] %in% blue && clade_list2[k] %in% green){
      col_clade<-c(col_clade,'green')}
    else{
      col_clade<-c(col_clade,'white')}
  }
}else if (map =='redto'){
  clade_list1<-cophylo2$assoc[,1]
clade_list2<-cophylo2$assoc[,2]
col_clade<-c()
for (k in 1:length(col_list)){
  if (clade_list1[k] %in% red && clade_list2[k] %in% yellow){
    col_clade<-c(col_clade,'yellow')}
  else if (clade_list1[k] %in% red && clade_list2[k] %in% green){
    col_clade<-c(col_clade,'green')}
  else{
    col_clade<-c(col_clade,'white')}
}
}else{
  
  clade_list1<-cophylo2$assoc[,1]
  clade_list2<-cophylo2$assoc[,2]
  col_clade<-c()
  for (k in 1:length(col_list)){
    if (clade_list1[k] %in% purple && clade_list2[k] %in% yellow){
      col_clade<-c(col_clade,'yellow')}
    else if (clade_list1[k] %in% purple && clade_list2[k] %in% green){
      col_clade<-c(col_clade,'green')}
    else{
      col_clade<-c(col_clade,'white')}
  }
}

lwd<-rep(3,length(col_clade))
lwd[which(col_clade=='white')]<-0
plot.cophylo(cophylo2,fsize=c(0.3,0.3),
     link.lty="solid",link.col=col_clade, link.type='curved',link.lwd=lwd)



lightyellow<-clade.members(x=MRCA(host_treeCAD2,'CC_1D_CARO', 'PK_6B_PICEA'),phy = host_treeCAD2,tip.labels=TRUE)
darkyellow<-clade.members(x=MRCA(host_treeCAD2,'CL_5C_CARO', 'BM_7I_PICEA'),phy = host_treeCAD2,tip.labels=TRUE)

lightgreen<-clade.members(x=MRCA(host_treeCAD2,'AG_5A_PICEA', 'OCO_4C_CARO'),phy = host_treeCAD2,tip.labels=TRUE)
darkgreen<-clade.members(x=MRCA(host_treeCAD2,'CL_4A_PICEA', 'TM_2E_RUDIS'),phy = host_treeCAD2,tip.labels=TRUE)



purple<-clade.members(x=MRCA(host_treeCO12,'EM_9D_PICEA', 'RC_5B_PICEA'),phy = host_treeCO12,tip.labels=TRUE)
blue<-clade.members(x=MRCA(host_treeCO12,'BMM_5B_PICEA', 'TC_6A_PICEA'),phy = host_treeCO12,tip.labels=TRUE)
red<-setdiff(host_treeCO12$tip.label,c('A_umphreyi',purple,blue))




#Count rudis
rudis_ly<-length(which(grepl('RUDIS',lightyellow)==TRUE))
rudis_dy<-length(which(grepl('RUDIS',darkyellow)==TRUE))
rudis_lg<-length(which(grepl('RUDIS',lightgreen)==TRUE))
rudis_dg<-length(which(grepl('RUDIS',darkgreen)==TRUE))
rudis_y<-rudis_ly+rudis_dy
rudis_g<-rudis_lg+rudis_dg

rudis_p<-length(which(grepl('RUDIS',purple)==TRUE))
rudis_b<-length(which(grepl('RUDIS',blue)==TRUE))
rudis_r<-length(which(grepl('RUDIS',red)==TRUE))

#Count picea
picea_ly<-length(which(grepl('PICEA',lightyellow)==TRUE))
picea_dy<-length(which(grepl('PICEA',darkyellow)==TRUE))
picea_lg<-length(which(grepl('PICEA',lightgreen)==TRUE))
picea_dg<-length(which(grepl('PICEA',darkgreen)==TRUE))
picea_y<-picea_ly+picea_dy
picea_g<-picea_lg+picea_dg

picea_p<-length(which(grepl('PICEA',purple)==TRUE))
picea_b<-length(which(grepl('PICEA',blue)==TRUE))
picea_r<-length(which(grepl('PICEA',red)==TRUE))


#carolinensis
caro_ly<-length(which(grepl('CARO',lightyellow)==TRUE))
caro_dy<-length(which(grepl('CARO',darkyellow)==TRUE))
caro_lg<-length(which(grepl('CARO',lightgreen)==TRUE))
caro_dg<-length(which(grepl('CARO',darkgreen)==TRUE))
caro_y<-caro_ly+caro_dy
caro_g<-caro_lg+caro_dg

caro_p<-length(which(grepl('CARO',purple)==TRUE))
caro_b<-length(which(grepl('CARO',blue)==TRUE))
caro_r<-length(which(grepl('CARO',red)==TRUE))


#fulva
fulva_ly<-length(which(grepl('FULVA',lightyellow)==TRUE))
fulva_dy<-length(which(grepl('FULVA',darkyellow)==TRUE))
fulva_lg<-length(which(grepl('FULVA',lightgreen)==TRUE))
fulva_dg<-length(which(grepl('FULVA',darkgreen)==TRUE))
fulva_y<-fulva_ly+fulva_dy
fulva_g<-fulva_lg+fulva_dg

fulva_p<-length(which(grepl('FULVA',purple)==TRUE))
fulva_b<-length(which(grepl('FULVA',blue)==TRUE))
fulva_r<-length(which(grepl('FULVA',red)==TRUE))

#Not enough fulva for testing
#data <- as.matrix(rbind(c(rudis_y, picea_y, caro_y, fulva_y), c(rudis_g, picea_g, caro_g, fulva_g)))
#dimnames(data) <- list(c("yellow", "green"), c("rudis", "picea", "caro", "fulva"))

data <- as.matrix(rbind(c(rudis_y, picea_y, caro_y,fulva_y), c(rudis_g, picea_g, caro_g,fulva_g)))
dimnames(data) <- list(c("yellow", "green"), c("rudis", "picea", "caro","fulva"))

fisher.test(data)


data <- as.matrix(rbind(c(rudis_r, picea_r, caro_r,fulva_r), c(rudis_p, picea_p, caro_p,fulva_p), c(rudis_b, picea_b, caro_b,fulva_b)))
dimnames(data) <- list(c("red", "purple","blue"), c("rudis", "picea", "caro","fulva"))

fisher.test(data)

data <- as.matrix(rbind(c(rudis_r, picea_r, caro_r,fulva_r), c(rudis_p, picea_p, caro_p,fulva_p)))
dimnames(data) <- list(c("red", "purple"), c("rudis", "picea", "caro","fulva"))
p1<-pairwise_fisher_test(data, p.adjust.method = "none")
data <- as.matrix(rbind(c(rudis_r, picea_r, caro_r,fulva_r), c(rudis_b, picea_b, caro_p,fulva_b)))
dimnames(data) <- list(c("red", "blue"), c("rudis", "picea", "caro","fulva"))
p2<-pairwise_fisher_test(data, p.adjust.method = "none")
data <- as.matrix(rbind(c(rudis_p, picea_p, caro_p,fulva_p), c(rudis_b, picea_b, caro_p,fulva_b)))
dimnames(data) <- list(c("purple", "blue"), c("rudis", "picea", "caro","fulva"))
p3<-pairwise_fisher_test(data, p.adjust.method = "none")

plist<-c(p1$p,p2$p,p3$p)
clist<-c(rep('red-purple',length(p1$p)),rep('red-blue',length(p2$p)),rep('purple-blue',length(p3$p)))
c1list<-c(p1$group1,p2$group1,p3$group1)
c2list<-c(p1$group2,p2$group2,p3$group2)


plistadj<-p.adjust(plist,method='BH')

pdf<-data.frame(clist,c1list,c2list,plistadj)

