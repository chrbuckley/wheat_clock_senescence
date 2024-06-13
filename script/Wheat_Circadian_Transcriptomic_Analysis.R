#################################################################################
######### Analyses of wheat senescence circadian transcriptomic dataset #########
#################################################################################

#####Introduction#####

#This script reproduces analyses and figures from Buckley et al., 2024

#The script uses the outputs of circadian transcript rhythm algorithms BIO_CYCLE and MetaCycle
#Please see methodology description within the journal article for further details

#####Loading BIO_CYCLE and MetaCycle data and filtering#####

#setwd("~/Wheat_Circadian_Transcriptomic_Analysis/")

#Required packages
library(dplyr)
library(tibble)
library(ggplot2)
library(patchwork)
library(viridis)
library(data.table)
library(ggfortify)
library(VennDiagram)
library(coin)
library(ggdist)
library(topGO)
library(grid)
library(gtable)
library(limorhyde)
library(qs)
library(limma)
library(foreach)
library(GENIE3)

#Mace and Chinese Spring gene orthologies
CS_IDs=read.csv("Mace_to_CS_gene_IDs.csv")

#Load BIO_CYCLE output
mat_bio=read.table("mature_bio_results.tsv", h=T)
sen_bio=read.table("senescence_bio_results.tsv", h=T)

#Order alphabetically
mat_bio=mat_bio[order(mat_bio$ID),]
sen_bio=sen_bio[order(sen_bio$ID),]

#Adjust phase to 24 h cycle
mat_bio$LAG = mat_bio$LAG*(24/mat_bio$PERIOD)
sen_bio$LAG = sen_bio$LAG*(24/sen_bio$PERIOD)

#Load MetaCycle output
mat_meta=read.csv("meta2d_mature_fpkm.csv")
sen_meta=read.csv("meta2d_senescence_fpkm.csv")

#Order alphabetically
mat_meta=mat_meta[order(mat_meta$CycID),]
sen_meta=sen_meta[order(sen_meta$CycID),]

#Adjust phase to 24 h cycle
mat_meta$meta2d_phase = mat_meta$meta2d_phase*(24/mat_meta$meta2d_period)
sen_meta$meta2d_phase = sen_meta$meta2d_phase*(24/sen_meta$meta2d_period)

#Implement filtering for expression. Remove transcripts that are not expressed at least once in each 24 hour period

#MATURE
remove_rows_mat=NULL
for (i in 1:nrow(mat_bio)){ #for all rows
  test=NULL
  for (x in 10:ncol(mat_bio)){ #and all columns with expression data
    test_val = mat_bio[i,x] > 0 #test whether values are above 0
    test=c(test,test_val)
  }
  if (length(which(grepl("TRUE", test[1:24]))) == "0" | length(which(grepl("TRUE", test[25:48]))) == "0"){
    remove_rows_mat = c(remove_rows_mat, i)  
  }
}
mat_bio=mat_bio[-remove_rows_mat,]
mat_meta=mat_meta[-remove_rows_mat,]

#SENESCENT
remove_rows_sen=NULL
for (i in 1:nrow(sen_bio)){ #for all rows
  test=NULL
  for (x in 10:ncol(sen_bio)){ #and all columns with expression data
    test_val = sen_bio[i,x] > 0 #test whether values are above 0
    test=c(test,test_val)
  }
  if (length(which(grepl("TRUE", test[1:24]))) == "0" | length(which(grepl("TRUE", test[25:48]))) == "0"){
    remove_rows_sen = c(remove_rows_sen, i)  
  }
}
sen_bio=sen_bio[-remove_rows_sen,]
sen_meta=sen_meta[-remove_rows_sen,]

#####Quality control of expression data via PCA#####

#Assessment of BIO_CYCLE's estimate of rhythmic transcripts

#Retain rhythmic transcripts only
mat_bio_rhy=filter(mat_bio, Q_VALUE < 0.05)
sen_bio_rhy=filter(sen_bio, Q_VALUE < 0.05)

#By timepoint
trans_mat_rhy = data.frame(t(mat_bio_rhy[,9:56]))
trans_sen_rhy = data.frame(t(sen_bio_rhy[,9:56]))

trans_mat_rhy= mutate(trans_mat_rhy, Timepoint=paste("ZT",substr(rownames(trans_mat_rhy),4,5),sep=""),
                      Stage="Mature")
trans_sen_rhy= mutate(trans_sen_rhy, Timepoint=paste("ZT",substr(rownames(trans_sen_rhy),4,5),sep=""),
                      Stage="Senescent")

pca_mat_tp = prcomp(trans_mat_rhy[,1:(ncol(trans_mat_rhy)-2)], center=TRUE,scale. = T)
pca_sen_tp = prcomp(trans_sen_rhy[,1:(ncol(trans_sen_rhy)-2)], center=TRUE,scale. = T)


#Plot by timepoint for mature tissue
autoplot(pca_mat_tp, x=1,y=2,data=trans_mat_rhy, colour="Timepoint", size=2)+
  scale_colour_viridis("Time point",discrete=T)+
  labs(x="\nPC1 (34.23%)", y="PC2 (14.14%)\n")+
  scale_y_continuous(expand=c(0.02,0.02))+
  scale_x_continuous(expand=c(0.02,0.02))+
  theme_bw()+
  theme(axis.line = element_line(size = 0.5, colour = "black"), 
        axis.ticks = element_blank(), panel.border = element_blank(), panel.grid = element_blank(),
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"), axis.title = element_text(size = 16),
        legend.text = element_text(size=12), legend.title = element_text(size=14))

#Plot by timepoint for senescent tissue
autoplot(pca_sen_tp, x=1,y=2,data=trans_sen_rhy, colour="Timepoint",size=2)+
  scale_colour_viridis(discrete=T)+
  labs(x="\nPC1 (34.8%)", y="PC2 (12.89%)\n")+
  scale_y_continuous(expand=c(0.02,0.02))+
  scale_x_continuous(expand=c(0.02,0.02))+
  theme_bw()+
  theme(axis.line = element_line(size = 0.5, colour = "black"), 
        axis.ticks = element_blank(), panel.border = element_blank(), panel.grid = element_blank(),
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"), axis.title = element_text(size = 16),
        legend.text = element_text(size=12), legend.title = element_text(size=14))

#Now do the same but combine mature and senescent to visualise separation by tissue type

#Firstly need to make sure all the genes are shared
overlap = intersect(mat_bio_rhy$ID, sen_bio_rhy$ID)

mat_shared = which(mat_bio_rhy$ID %in% overlap)
sen_shared = which(sen_bio_rhy$ID %in% overlap)
trans_mat_rhy=trans_mat_rhy[,mat_shared]
trans_sen_rhy=trans_sen_rhy[,sen_shared]
trans_mat_rhy= mutate(trans_mat_rhy, Timepoint=paste("ZT",substr(rownames(trans_mat_rhy),4,5),sep=""),
                      Stage="Mature")
trans_sen_rhy= mutate(trans_sen_rhy, Timepoint=paste("ZT",substr(rownames(trans_sen_rhy),4,5),sep=""),
                      Stage="Senescent")

rownames(trans_mat_rhy)=paste("M", rownames(trans_mat_rhy), sep=" ")
rownames(trans_sen_rhy)=paste("S", rownames(trans_sen_rhy), sep=" ")
colnames(trans_mat_rhy)[1:(ncol(trans_mat_rhy)-2)]=1:(ncol(trans_mat_rhy)-2)
colnames(trans_sen_rhy)[1:(ncol(trans_sen_rhy)-2)]=1:(ncol(trans_sen_rhy)-2)

trans_all_rhy = rbind(trans_mat_rhy,trans_sen_rhy)
trans_all_rhy$Stage = as.factor(trans_all_rhy$Stage)
pca_all_tp = prcomp(trans_all_rhy[,1:(ncol(trans_all_rhy)-2)], center=TRUE, scale. = TRUE)


autoplot(pca_all_tp, x=1,y=2,data=trans_all_rhy, colour="Stage")+
  scale_colour_manual(values=c("#197000","#BD9B00"))+
  labs(x="\nPC1 (33.93%)", y="PC2 (17.31%)\n")+
  scale_y_continuous(expand=c(0.02,0.02))+
  scale_x_continuous(expand=c(0.02,0.02))+
  theme_bw()+
  theme(axis.line = element_line(size = 0.5, colour = "black"), 
        axis.ticks = element_blank(), panel.border = element_blank(), panel.grid = element_blank(),
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"), axis.title = element_text(size = 16),
        legend.text = element_text(size=12), legend.title = element_text(size=14))


#Visualise the overlap of rhythmicity estimates of BIO_CYCLE and MetaCycle

#Filter MetaCycle output to retain rhythmic transcripts
mat_meta_rhy=filter(mat_meta, meta2d_BH.Q < 0.05)
sen_meta_rhy=filter(sen_meta, meta2d_BH.Q < 0.05)

#Intersect of BIO_CYCLE and MetaCycle
mat_rhy=subset(mat_bio_rhy, ID %in% mat_meta_rhy$CycID)
sen_rhy=subset(sen_bio_rhy, ID %in% sen_meta_rhy$CycID)

#Make a figure of the intersect!

venn.diagram(
  x=list(mat_bio_rhy$ID, sen_bio_rhy$ID,
         mat_meta_rhy$CycID, sen_meta_rhy$CycID),
  category.names = c("Mature\nBIO_CYCLE" , "Senescent\nBIO_CYCLE" , "Mature\nMetaCycle","Senescent\nMetaCycle"),
  output=T,
  filename = 'venn_diagram_bio_meta.png',
  imagetype = "png",
  fill=c("#197000","#BD9B00","#197000","#BD9B00"),
  colour=c("#197000","#BD9B00","#197000","#BD9B00"),
  lwd=2,
  lty="blank",
  fontfamily="helvetica",
  cat.fontfamily="helvetica",
  cat.cex=1.5,
  cex=1.5,
  width=5050,
  height=4000,
  resolution=600)

#####Analysis of transcripts that are uniquely rhythmic in senescent tissue#####

overlap = intersect(mat_rhy[,1], sen_rhy[,1])

mat_shared = subset(mat_rhy, ID %in% overlap)
sen_shared = subset(sen_rhy, ID %in% overlap)

#Uniquely rhythmic in mature tissue
mat_unique = mat_rhy[-which(mat_rhy$ID %in% overlap),]
mat_unique = mat_unique[-which(mat_unique$ID %in% sen_bio_rhy$ID),]
mat_unique = mat_unique[-which(mat_unique$ID %in% sen_meta_rhy$CycID),]

#Uniquely rhythmic in senescent tissue
sen_unique = sen_rhy[-which(sen_rhy$ID %in% overlap),]
sen_unique = sen_unique[-which(sen_unique$ID %in% mat_bio_rhy$ID),]
sen_unique = sen_unique[-which(sen_unique$ID %in% mat_meta_rhy$CycID),]

#Make a simple bar plot

ggplot()+
  geom_bar(mapping=aes(x=c(1,2,3),y=c(2072,8593,5709)),stat="identity",fill=c("#197000","grey","#BD9B00"))+
  geom_text(mapping=aes(x=c(1,2,3),y=c(1472,7993,5109),label=c("Mat\n(2072)","Both\n(8593)","Sen\n(5709)")),size=6)+
  theme_bw()+
  labs(y="Gene count\n")+
  theme(axis.line.x = element_line(size = 0.5, colour = "black"),axis.line.y = element_line(size = 0.5, colour = "black"),
        axis.ticks = element_line(colour="NA"), panel.border = element_blank(), panel.grid = element_blank(),
        axis.text.x = element_blank(),axis.title.x = element_blank(),
        axis.text.y = element_text(size = 16, colour = "black"), axis.title = element_text(size = 20))+
  scale_y_continuous(expand=c(0,0))


#GO enrichment analysis of uniquely rhythmic transcripts in senescent tissue

go_terms = read.csv("go_terms.txt", stringsAsFactors = F)
go_terms=unique(go_terms)

#Remove genes without GO terms
go_terms = go_terms[which(go_terms$GO.term.accession != ""),]

#Create list with element for each gene, containing vectors with all terms for each gene
gene2GO = tapply(go_terms$GO.term.accession, go_terms$Gene.stable.ID, function(x)x)

#Define gene list
#Assign 1 if gene is in test list, 0 if not
#sen unique
tmp <- ifelse(go_terms$Gene.stable.ID %in% sen_unique$ID, 1, 0)
genelist_sen_unique=tmp
names(genelist_sen_unique)=go_terms$Gene.stable.ID

#Create the GO data classes
GOdata_sen_unique <- new("topGOdata",
                         ontology = "BP",
                         allGenes = genelist_sen_unique,
                         geneSelectionFun = function(x)(x == 1),
                         annot = annFUN.gene2GO, gene2GO = gene2GO)
#Fisher's exact test for GO enrichment!
Fisher_sen_unique <- runTest(GOdata_sen_unique, algorithm = "classic", statistic = "fisher")
go_results_sen_unique <- GenTable(GOdata_sen_unique, raw.p.value = Fisher_sen_unique, topNodes = length(Fisher_sen_unique@score))
#Clean up the results for processing
go_results_sen_unique$raw.p.value = gsub("e", "E", go_results_sen_unique$raw.p.value)
go_results_sen_unique$raw.p.value = gsub("< ", "", go_results_sen_unique$raw.p.value)
go_results_sen_unique$raw.p.value = as.numeric(go_results_sen_unique$raw.p.value)
go_results_sen_unique$Observed.Expected = go_results_sen_unique$Significant/go_results_sen_unique$Expected

#Curate a top results list to remove duplicate terms
#Order of this might vary slightly from run to run
su=go_results_sen_unique[c(1,2,3,4,5,6,7,8,9,10,11,13,16,18,19,20,21,22,23,24),]
su$Term=with(su, reorder(Term , Observed.Expected, mean , na.rm=T))

ggplot(su,aes(x=log(Observed.Expected),y=Term,fill=-log(raw.p.value)))+
  geom_bar(stat="identity",width=0.15,colour="transparent")+
  geom_point(su,mapping=aes(log(Observed.Expected),Term,fill=-log(raw.p.value),size=Significant),pch=21,stroke=0,show.legend = T)+
  scale_fill_gradientn(expression(paste("-log(",italic(p),")")),colours=c("#f7d7a6","#fcd69d","darkorange","red","red1"),values=c(0,0.15,0.4,0.75,1), breaks=seq(30,70,8))+
  labs(x="\nlog(Observed/Expected)",y="")+
  theme_bw()+
  theme(panel.background = element_blank(),axis.line = element_line(size = 0.5, colour = "black"), 
        axis.ticks = element_line(colour="NA"), panel.border = element_blank(), panel.grid = element_blank(),
        axis.text.x = element_text(size = 10, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"), axis.title = element_text(size = 12),
        legend.text = element_text(size=10), legend.title = element_text(size=12))+
  guides(size=guide_legend(override.aes = list(fill="black")))

#A different method of assessing gain or loss of rhythmicity is to look at average 'periodicity'. This is BIO_CYCLE's measure of rhythmicity
#Do this for different families of TFs

#Load TFs from Ramirez-Gonzalez et al 2018
tfs = read.csv("transcription_factors_to_use_high_confidence.csv")
#These are from CS. Filter against Mace.
tfs_Mace=subset(CS_IDs, Triticum.aestivum.gene.stable.ID %in% tfs$locus)
tfs_Mace=merge(tfs_Mace, tfs, by.x="Triticum.aestivum.gene.stable.ID", by.y="locus")
tfs_Mace$identifier=1:nrow(tfs_Mace)
tfs_Mace=tfs_Mace%>%group_by(Gene.stable.ID)%>%filter(identifier== max(identifier)) %>% 
  distinct

#Get average periodicity of superfamilies of TFs

#For this, need all data before filtering for rhythmicity or expression.
#This is because transcripts that are called rhythmic by BIO_CYCLE (q < 0.05) have periodicity typically close to 1
#so the comparison would be less meaningful
mat=read.table("mature_bio_results.tsv", h=T)
sen=read.table("senescence_bio_results.tsv", h=T)

sen_tfs = subset(sen, ID %in% tfs_Mace$Gene.stable.ID)
mat_tfs = subset(mat, ID %in% tfs_Mace$Gene.stable.ID)
sen_tfs=mutate(merge(sen_tfs,tfs_Mace, by.x="ID",by.y="Gene.stable.ID"),stage="Senescent")
mat_tfs=mutate(merge(mat_tfs,tfs_Mace, by.x="ID",by.y="Gene.stable.ID"),stage="Mature")

#Summarise periodicity by family
tfs_fams = rbind(sen_tfs,mat_tfs) %>% group_by(superfamily,stage)%>%
  summarise(mean_periodicity = mean(AVG_REPS_PERIODICITY),
            mean_q = mean(Q_VALUE),n=n(),
            sd=sd(AVG_REPS_PERIODICITY))

#Test for significant differences in mean periodicity between mature and senescent by t-test
tfs_fams$superfamily=as.character(tfs_fams$superfamily)
tfs_fams$sig=NA
tfs_fams$level=NA
#Retain TF families with >10 members
tfs_fams=tfs_fams[-(which(tfs_fams$n < 10)),]

for (i in seq(2,nrow(tfs_fams),2)){
  x=sen_tfs[sen_tfs$superfamily == unlist(tfs_fams[i,1]),]$AVG_REPS_PERIODICITY
  y=mat_tfs[mat_tfs$superfamily == unlist(tfs_fams[i,1]),]$AVG_REPS_PERIODICITY
  t=t.test(x,y)
  print(unlist(tfs_fams[i,1]))
  print(t)
  tfs_fams$sig[c(i,i-1)]=t$p.value
  if(tfs_fams$sig[i] < 0.05){tfs_fams$level[c(i,i+1)]="*"}
  if(tfs_fams$sig[i] < 0.01){tfs_fams$level[c(i,i+1)]="**"}
  if(tfs_fams$sig[i] < 0.001){tfs_fams$level[c(i,i+1)]="***"}
}

all_tfs_fams=rbind(sen_tfs,mat_tfs)
tfs_fams$superfamily =with(tfs_fams, reorder(superfamily , mean_periodicity, mean , na.rm=T))

#Calculate differences in periodicity between mature and senescent
diffs_periodicity=data.frame(fam=tfs_fams[seq(1,94,2),1],diffs=tfs_fams[seq(1,94,2),3]-tfs_fams[seq(2,94,2),3])
arrange(diffs_periodicity, mean_periodicity)

ggplot(tfs_fams, aes(y=superfamily, x=mean_periodicity))+
  geom_line()+
  geom_point(mapping=aes(colour=stage,size=n))+
  scale_size_continuous("Gene count",range=c(1,5),breaks=c(10,50,250,500))+
  scale_color_manual("Stage",values=c("#197000","#BD9B00"))+
  scale_fill_manual(values=c("#197000","#BD9B00"))+
  geom_text(tfs_fams[tfs_fams$stage=="Senescent",c(-2)], mapping=aes(y=superfamily,x=mean_periodicity+0.02,label=level),nudge_y = -0.25,hjust=0,size=5)+
  labs(x="\nMean periodicity",y="")+
  theme_bw()+
  theme(panel.background = element_blank(),axis.line = element_line(size = 0.5, colour = "black"),
        axis.ticks = element_line(colour="NA"), panel.border = element_blank(), panel.grid = element_blank(),
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"), axis.title = element_text(size = 14),
        legend.text = element_text(size=12), legend.title = element_text(size=14))+
  guides(size=guide_legend(override.aes = list(fill="black")))+
  guides(colour=guide_legend(override.aes = list(size=4)))


#This looks particularly interesting for WRKY TFs. Now we can investigate this for individual WRKYs

WRKY=filter(all_tfs_fams, superfamily == "WRKY")
WRKY_sen=filter(WRKY, stage == "Senescent")
WRKY_mat=filter(WRKY, stage == "Mature")

#Retain WRKYs that are the most rhythmic in senescent flag leaves. This is a sanity check.
keep=filter(WRKY, stage == "Senescent" & AVG_REPS_PERIODICITY > 0.99)
WRKY=WRKY[WRKY$ID %in% keep$ID,]
WRKY_ID=WRKY$ID
WRKY=data.frame(ID=WRKY_ID[1:(nrow(WRKY)/2)],
                periodicity_diff=(WRKY[WRKY$stage=="Mature",]$AVG_REPS_PERIODICITY-WRKY[WRKY$stage=="Senescent",]$AVG_REPS_PERIODICITY),
                periodicity_sen=WRKY[WRKY$stage=="Senescent",]$AVG_REPS_PERIODICITY,
                periodicity_mat=WRKY[WRKY$stage=="Mature",]$AVG_REPS_PERIODICITY)

#Retain largest periodicity differences. Arbritrary cut off as I want to plot 20.
wrky=as.character(filter(WRKY,periodicity_diff < -0.905)$ID)

#Plot expression comparisons
#Firstly make a new df for plotting expression data easily. This is handy for lots of future analyses.
mat_names=paste("mat_",colnames(mat[,9:56]),sep="")
sen_names=paste("sen_",colnames(sen[,9:56]),sep="")
all=cbind(mat[,9:56],sen[,9:56])
colnames(all)=c(mat_names,sen_names)
rownames(all)=mat$ID
metadata=data.table(sample_id=c(mat_names,sen_names),
                    time=as.numeric(substr(c(mat_names,sen_names),8,9)),
                    cond=substr(c(mat_names,sen_names),1,3))

#Filter for target WRKYs
df = data.frame(t(all[which(rownames(all) %in% wrky),]))
df$sample_id= metadata$sample_id
df = merge(df, metadata[, .(sample_id, cond, time)], by = 'sample_id')

df_=NULL
for (x in 2:(ncol(df)-2)){
  y=data.frame(df[,c(x,22,23)])
  colnames(y)[1]="V1"
  df_=rbind(df_,y)
}
df_$cond=rep(c(rep("Mature",48),rep("Senescent",48)),20)
df_$locus=rep(wrky,each=96)

df_means = df_ %>%
  group_by(time,cond,locus)%>%
  summarise(mean=mean(V1),
            sd=sd(V1)/sqrt(n()))


ggplot(df_means, aes(y=mean,x=time,colour=cond))+
  facet_wrap(~locus, scales="free_y",nrow=5)+
  geom_point(size=1)+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), size = 0.5, width = 0)+
  geom_line()+
  scale_colour_manual("",labels =c("Mature","Senescent"),values=c("#197000","#BD9B00"))+
  #scale_shape_discrete("",labels=)+
  labs(x = "\nTime (ZT)", y="Mean expression (FPKM)\n")+
  theme_bw()+
  theme(axis.ticks = element_line(colour="NA"), panel.border = element_rect(colour="light grey",size=1), 
        panel.grid.major = element_line(colour="light grey", size=0.1),
        panel.grid.minor = element_blank(),axis.text.x = element_text(size = 10, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"), axis.title = element_text(size = 12),
        strip.background = element_rect(colour=NA,fill="NA"), strip.text = element_text(size=9),
        legend.text = element_text(size=10), legend.position = "top")+
  guides(shape = guide_legend(override.aes = list(size = 3)))


#Now plot a metagene (average expression across all genes) for each TF family for further visualisation
fams_df=NULL

for (i in unique(all_tfs_fams$superfamily)){
  x=filter(all_tfs_fams, superfamily == i)
  x=as.character(unique(x$ID)) 
  df = data.frame(t(all[which(rownames(all) %in% x),]))
  df$sample_id= metadata$sample_id
  df = merge(df, metadata[, .(sample_id, cond, time)], by = 'sample_id')
  
  df_=NULL
  for (z in 2:(ncol(df)-2)){
    y=data.frame(df[,c(z,ncol(df)-1,ncol(df))])
    colnames(y)[1]="V1"
    df_=rbind(df_,y)
  }
  
  df_$cond=rep(c(rep("Mature",48),rep("Senescent",48)),length(x))
  df_$superfamily=rep(i,each=(length(x)*96))
  
  fams_df=rbind(fams_df,df_)
  
}

df_means = fams_df %>%
  group_by(time,cond,superfamily)%>%
  summarise(mean=mean(V1),
            sd=sd(V1)/sqrt(n()))

#Subset to a few relevant families
#df_means=subset(df_means, superfamily %in% c("NAC", "WRKY", "MYB-related","AP2/EREBP","Pseudo ARR-B", "Sigma70-like"))

ggplot(df_means, aes(y=mean,x=time,colour=cond))+
  facet_wrap(~superfamily, scales="free_y",nrow=9)+
  geom_point(size=1)+
  #geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), size = 0.5, width = 0)+
  geom_line()+
  scale_colour_manual("",labels =c("Mature","Senescent"),values=c("#197000","#BD9B00"))+
  labs(x = "\nTime (ZT)", y="Mean expression (FPKM)\n")+
  theme_bw()+
  theme(axis.ticks = element_line(colour="NA"), panel.border = element_rect(colour="light grey",size=1), 
        panel.grid.major = element_line(colour="light grey", size=0.1),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"), axis.title = element_text(size = 14),
        strip.background = element_rect(colour=NA,fill=NA), strip.text = element_text(size=10),
        legend.text = element_text(size=12),legend.position = "top")+
  guides(shape = guide_legend(override.aes = list(size = 2)))



#####Analysis of rhythmic parameters of gloabl rhythmic genes#####

mat_rhy_ = mutate(mat_rhy, Stage="Mature")
sen_rhy_ = mutate(sen_rhy, Stage="Senescence")
all_rhy = rbind(mat_rhy_, sen_rhy_)

#Circadian period
#Use permutation test of the distribution
oneway_test(PERIOD ~ as.factor(Stage), data=all_rhy,distribution = approximate(nresample = 10000))
#Or t-test
t.test(mat_rhy_$PERIOD,sen_rhy_$PERIOD)

#Period vs developmental stage plot

ggplot(all_rhy, aes(y=Stage, x=PERIOD,colour=Stage,fill=Stage))+
  stat_halfeye(.width=0,point_colour=NA,justification=-0.105,alpha=0.6,adjust=0.5,height=0.9)+
  geom_jitter(height=0.075,size=0.2,alpha=0.2,show.legend = F)+
  geom_boxplot(width=0.075, fill="transparent",colour="black",show.legend = F,outlier.colour = NA)+
  scale_colour_manual(values=c("#197000","#BD9B00"),guide="none")+
  scale_fill_manual(values=c("#197000","#BD9B00"))+
  scale_y_discrete(expand=c(0.075,0.075),labels=c("",""))+
  labs(y = "", x="\nCircadian period (h)")+
  theme_bw()+
  theme(axis.line.x = element_line(size = 0.5, colour = "black"),axis.line.y = element_blank(),
        axis.ticks = element_line(colour="NA"), panel.border = element_blank(), panel.grid = element_blank(),
        axis.text.x = element_text(size = 16, colour = "black"),
        axis.text.y = element_text(size = 16, colour = "black"), axis.title = element_text(size = 18),
        legend.text = element_text(size=16), legend.title = element_text(size=16))+
  guides(fill = guide_legend(override.aes = list(colour=c("#197000","#BD9B00"),fill=c("#197000","#BD9B00"),stroke=0,linetype="blank")))

#Phase vs developmental stage plot

#Use permutation test of the distribution. This is complicated by the non-linear nature of phase
oneway_test(LAG ~ as.factor(Stage), data=all_rhy,distribution = approximate(nresample = 10000))
#Or t-test
t.test(mat_rhy_$LAG,sen_rhy_$LAG)

ggplot(all_rhy, aes(y=Stage, x=LAG,colour=Stage,fill=Stage))+
  stat_halfeye(.width=0,point_colour=NA,justification=-0.095,alpha=0.6,adjust=0.5)+
  geom_jitter(height=0.075,size=0.2,alpha=0.2,show.legend = F)+
  geom_boxplot(width=0.075, fill="transparent",colour="black",show.legend = F)+
  scale_colour_manual(values=c("#197000","#BD9B00"),guide="none")+
  scale_fill_manual(values=c("#197000","#BD9B00"))+
  scale_y_discrete(expand=c(0.075,0.075),labels=c("",""))+
  labs(y = "", x="\nCircadian phase (h)")+
  theme_bw()+
  theme(axis.line.x = element_line(size = 0.5, colour = "black"),axis.line.y = element_blank(),
        axis.ticks = element_line(colour="NA"), panel.border = element_blank(), panel.grid = element_blank(),
        axis.text.x = element_text(size = 16, colour = "black"),
        axis.text.y = element_text(size = 16, colour = "black"), axis.title = element_text(size = 18),
        legend.text = element_text(size=16), legend.title = element_text(size=16))+
  scale_x_continuous(limits = c(0,24),breaks=c(0,4,8,12,16,20,24))+
  guides(fill = guide_legend(override.aes = list(colour=c("#197000","#BD9B00"),fill=c("#197000","#BD9B00"),stroke=0,linetype="blank",size=2)))

#Now compare circadian parameters at the gene level. I.e. what is the difference from mature to senescent tissue for each gene.

phase_diff = mat_shared$LAG-sen_shared$LAG
period_diff = mat_shared$PERIOD-sen_shared$PERIOD
amp_diff = log(mat_shared$AMPLITUDE/sen_shared$AMPLITUDE)

diffs_df = tibble(phase_diff,period_diff,amp_diff)
diffs_df = cbind(mat_shared[,1:2], diffs_df)

#Adjust phase differences below and above 12 by adding 24 or subtracting 24
diffs_df[diffs_df$phase_diff > 12,]$phase_diff = diffs_df[diffs_df$phase_diff > 12,]$phase_diff-24
diffs_df[diffs_df$phase_diff < -12,]$phase_diff = diffs_df[diffs_df$phase_diff < -12,]$phase_diff+24

#Check mean gene-level differences
mean(diffs_df$phase_diff)
mean(diffs_df$period_diff)
mean(diffs_df$amp_diff)

#Check distribution of differences. Are they shifted away from 0?
#Period
ggplot(diffs_df, aes(x=period_diff)) + 
  geom_density(alpha = 0.5, size = 0.5,fill=viridis(4)[1],colour=viridis(4)[1])+
  geom_vline(linetype=1,size=0.5,xintercept=0,colour="red")+
  geom_vline(linetype=2,size=0.5,xintercept=mean(diffs_df$period_diff))+
  scale_colour_manual("Mean",values="black",labels="")+
  geom_text(data.frame(x=(0.52)+(12*0.03),y=(0.34)*0.98,label="0.52 h"),mapping=aes(x=x,y=y,label=label),hjust=0)+
  guides(colour=guide_legend(keywidth = 4))+
  scale_y_continuous(expand=c(0,0),limits=c(0,0.34))+
  labs(x="\nPeriod difference (h)",y="Density\n")+
  theme_bw()+
  theme(panel.background = element_blank(),axis.line = element_line(size = 0.5, colour = "black"), 
        axis.ticks = element_line(colour="NA"), panel.border = element_blank(), panel.grid = element_blank(),
        axis.text.x = element_text(size = 10, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"), axis.title = element_text(size = 11),
        legend.text = element_text(size=10), legend.title = element_text(size=12),legend.position = "none")

#Phase
ggplot(diffs_df, aes(x=phase_diff)) + 
  geom_density(alpha = 0.5, size = 0.5,fill=viridis(4)[2],colour=viridis(4)[2])+
  geom_vline(linetype=1,size=0.5,xintercept=0,colour="red")+
  geom_vline(linetype=2,size=0.5,xintercept=mean(diffs_df$phase_diff))+
  scale_colour_manual("Mean",values="black",labels="")+
  geom_text(data.frame(x=(0.44)+(24*0.03),y=(0.25)*0.98,label="0.44 h"),mapping=aes(x=x,y=y,label=label),hjust=0)+
  guides(colour=guide_legend(keywidth = 4))+
  scale_y_continuous(expand=c(0,0),limits=c(0,0.25))+
  labs(x="\nPhase difference (h)",y="Density\n")+
  theme_bw()+
  theme(panel.background = element_blank(),axis.line = element_line(size = 0.5, colour = "black"), 
        axis.ticks = element_line(colour="NA"), panel.border = element_blank(), panel.grid = element_blank(),
        axis.text.x = element_text(size = 10, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"), axis.title = element_text(size = 11),
        legend.text = element_text(size=10), legend.title = element_text(size=12),legend.position = "none",axis.title.y= element_blank())

#Amplitude
ggplot(diffs_df, aes(x=amp_diff)) + 
  geom_density(alpha = 0.5, size = 0.5,fill=viridis(4)[3],colour=viridis(4)[3])+
  geom_vline(linetype=1,size=0.5,xintercept=0,colour="red")+
  geom_vline(linetype=2,size=0.5,xintercept=mean(diffs_df$amp_diff))+
  scale_colour_manual("Mean",values="black",labels="")+
  geom_text(data.frame(x=(0.10)+(4*0.03),y=(1.6)*0.98,label="0.10"),mapping=aes(x=x,y=y,label=label),hjust=0)+
  guides(colour=guide_legend(keywidth = 4))+
  scale_y_continuous(expand=c(0,0),limits=c(0,1.6))+
  scale_x_continuous(limits=c(-2,2))+
  labs(x="\nAmplitude change (log(Mature:senescence))",y="Density\n")+
  theme_bw()+
  theme(panel.background = element_blank(),axis.line = element_line(size = 0.5, colour = "black"), 
        axis.ticks = element_line(colour="NA"), panel.border = element_blank(), panel.grid = element_blank(),
        axis.text.x = element_text(size = 10, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"), axis.title = element_text(size = 11),
        legend.text = element_text(size=10), legend.title = element_text(size=12),legend.position = "none",axis.title.y= element_blank())

#####Categorise rhythmic genes by differences in their parameters#####
#Then analyse within these categories

#Categorise by changes in period into short period transcripts, long and unchanged
short_period = filter(diffs_df, period_diff > 0.5)
long_period = filter(diffs_df, period_diff < -0.5)
unchanged_period = filter(diffs_df, period_diff < 0.5 & period_diff > -0.5)

#Perform GO enrichment on these categories as before
tmp <- ifelse(go_terms$Gene.stable.ID %in% short_period$ID, 1, 0)
genelist_short_period=tmp
names(genelist_short_period)=go_terms$Gene.stable.ID

tmp <- ifelse(go_terms$Gene.stable.ID %in% long_period$ID, 1, 0)
genelist_long_period=tmp
names(genelist_long_period)=go_terms$Gene.stable.ID

tmp <- ifelse(go_terms$Gene.stable.ID %in% unchanged_period$ID, 1, 0)
genelist_unchanged_period=tmp
names(genelist_unchanged_period)=go_terms$Gene.stable.ID


GOdata_short_period <- new("topGOdata",
                           ontology = "BP",
                           allGenes = genelist_short_period,
                           geneSelectionFun = function(x)(x == 1),
                           annot = annFUN.gene2GO, gene2GO = gene2GO)

GOdata_long_period <- new("topGOdata",
                          ontology = "BP",
                          allGenes = genelist_long_period,
                          geneSelectionFun = function(x)(x == 1),
                          annot = annFUN.gene2GO, gene2GO = gene2GO)

GOdata_unchanged_period <- new("topGOdata",
                               ontology = "BP",
                               allGenes = genelist_unchanged_period,
                               geneSelectionFun = function(x)(x == 1),
                               annot = annFUN.gene2GO, gene2GO = gene2GO)

Fisher_short_period <- runTest(GOdata_short_period, algorithm = "classic", statistic = "fisher")
Fisher_long_period <- runTest(GOdata_long_period, algorithm = "classic", statistic = "fisher")
Fisher_unchanged_period <- runTest(GOdata_unchanged_period, algorithm = "classic", statistic = "fisher")

go_results_short_period <- GenTable(GOdata_short_period, raw.p.value = Fisher_short_period, topNodes = length(Fisher_short_period@score),
                                    numChar = 120)

go_results_long_period <- GenTable(GOdata_long_period, raw.p.value = Fisher_long_period, topNodes = length(Fisher_long_period@score),
                                   numChar = 120)

go_results_unchanged_period <- GenTable(GOdata_unchanged_period, raw.p.value = Fisher_unchanged_period, topNodes = length(Fisher_unchanged_period@score),
                                        numChar = 120)

go_results_short_period$raw.p.value = gsub("e", "E", go_results_short_period$raw.p.value)
go_results_short_period$raw.p.value = gsub("< ", "", go_results_short_period$raw.p.value)
go_results_short_period$raw.p.value = as.numeric(go_results_short_period$raw.p.value)

go_results_long_period$raw.p.value = gsub("e", "E", go_results_long_period$raw.p.value)
go_results_long_period$raw.p.value = gsub("< ", "", go_results_long_period$raw.p.value)
go_results_long_period$raw.p.value = as.numeric(go_results_long_period$raw.p.value)

go_results_unchanged_period$raw.p.value = gsub("e", "E", go_results_unchanged_period$raw.p.value)
go_results_unchanged_period$raw.p.value = gsub("< ", "", go_results_unchanged_period$raw.p.value)
go_results_unchanged_period$raw.p.value = as.numeric(go_results_unchanged_period$raw.p.value)

#PLOTTING PERIOD GO ENRICHMENTS
#First calculate observed/expected and add groupings
go_results_short_period$Observed.Expected = (go_results_short_period$Significant-go_results_short_period$Expected)/go_results_short_period$Expected
go_results_long_period$Observed.Expected = (go_results_long_period$Significant-go_results_long_period$Expected)/go_results_long_period$Expected
go_results_unchanged_period$Observed.Expected = (go_results_unchanged_period$Significant-go_results_unchanged_period$Expected)/go_results_unchanged_period$Expected
go_results_short_period$Change = "Short"
go_results_long_period$Change = "Long"
go_results_unchanged_period$Change = "Unchanged"

go_period = rbind(go_results_short_period,go_results_long_period,go_results_unchanged_period)
go_period$Enrichment=NA

#Curate top results lists to remove duplicate terms
#Order of this might vary slightly from run to run
top10_short = head(go_results_short_period$Term,20)[c(1:4,6:8,10:12)]
positions_s=NULL
for (i in 1:length(top10_short)){
  x=grep(noquote(paste('^',top10_short[i],'$',sep="")), go_period$Term)
  positions_s=c(positions_s,x)
}
go_period$Enrichment[positions_s]="Short period-enriched"

top10_long = head(go_results_long_period$Term,20)[c(1,3,4,5,6,7,9,11,12,14)]
positions_l=NULL
for (i in 1:length(top10_long)){
  x=grep(noquote(paste('^',top10_long[i],'$',sep="")), go_period$Term)
  positions_l=c(positions_l,x)
}
go_period$Enrichment[positions_l]="Long period-enriched"

top10_unchanged_per = head(go_results_unchanged_period$Term,20)[c(1,2,3,4,6,7,8,10,11,12)]
positions_u=NULL
for (i in 1:length(top10_unchanged_per)){
  x=grep(noquote(paste('^',top10_unchanged_per[i],'$',sep="")), go_period$Term)
  positions_u=c(positions_u,x)
}
go_period$Enrichment[positions_u]="Unchanged period-enriched"

#filter by these top terms
go_period=na.omit(go_period)#remove rows without annotated enrichment
go_period$raw.p.value=as.numeric(go_period$raw.p.value)

#Plot this to summarise
go_period$Change = factor(go_period$Change, levels=c("Short","Unchanged","Long"))
go_period$Enrichment = factor(go_period$Enrichment, levels=c("Short period-enriched","Unchanged period-enriched","Long period-enriched"))

ggplot(go_period, aes(x=Change, y=Term,size=ifelse(Observed.Expected==-1, NA, Observed.Expected), colour=-log(raw.p.value)))+
  facet_wrap(~Enrichment, nrow=3, scales="free_y")+
  geom_point(alpha=0.75)+
  scale_x_discrete(labels=c("Short\n(n=4139)","Unchanged\n(n=2301)","Long\n(n=2153)"))+
  scale_color_gradientn(expression(paste("-log(",italic(p),")value")),colours=c("blue","blue","purple","red","red"),values=c(0,0.1,0.5,0.9,1))+
  scale_size_continuous("(O-E)/E", range=c(0.5,8), breaks=c(1,5,10,25,50))+
  labs(y="",x="\nPeriod change")+
  theme_bw()+
  theme(panel.background = element_blank(),axis.line = element_line(size = 0.5, colour = "black"), 
        axis.ticks = element_line(colour="NA"), panel.border = element_blank(), panel.grid = element_blank(),
        axis.text.x = element_text(size = 11, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"), axis.title = element_text(size = 14), 
        legend.text = element_text(size=12), legend.title = element_text(size=14),
        strip.background = element_rect(fill=NA,colour=NA), strip.text = element_text(size=11))


#Now provide example circadian rhythm changes of genes in the significantly enriched GO categories

#short - carbohydrate metabolic process
carbohydrate=unique(readLines("carbohydrate_metabolic_process.txt")[-1])
carbohydrate_=subset(short_period, ID %in% carbohydrate)
#select transcripts that best demonstrate the change in period in question
carbohydrate_=subset(diffs_df, ID %in% carbohydrate_$ID)[c(14,31,84),]
carbohydrate_$term="carbohydrate\n metabolic process"
#unchanged - circadian rhythm
rhythmic=unique(readLines("circadian_rhythm.txt"))[-1]
rhythmic_=subset(unchanged_period, ID %in% rhythmic)
rhythmic_=subset(diffs_df, ID %in% rhythmic_$ID)[c(7,8,9),]
rhythmic_$term="rhythmic\n process"
#long - photosynthesis, light harvesting
photo=unique(readLines("photosynthesis.txt")[-1])
photo_=subset(long_period, ID %in% photo)
photo_=subset(diffs_df, ID %in% photo_$ID)[c(14,17,18),]
photo_$term="photosynthesis,\n light harvesting"

period_combine=rbind(carbohydrate_,rhythmic_,photo_)

#Get expression data for these genes
df=data.frame(0)
for (i in 1:nrow(period_combine)){
  x = data.frame(t(all[which(rownames(all) %in% as.character(period_combine$ID[i])),]))
  colnames(x)=period_combine$ID[i]
  df=cbind(df,x)
}
df=df[,-1]
df$sample_id= metadata$sample_id
df = merge(df, metadata[, .(sample_id, cond, time)], by = 'sample_id')

df_=NULL
for (x in 2:(ncol(df)-2)){
  y=data.frame(df[,c(x,11,12)])
  colnames(y)[1]="V1"
  df_=rbind(df_,y)
}

df_$cond=rep(c(rep("Mature",48),rep("Senescent",48)),9)
df_$locus=rep(as.character(period_combine$ID),each=96)
df_$term=rep(period_combine$term, each=96)
df_$term=factor(df_$term, levels=unique(df_$term))


df_means = df_ %>%
  group_by(time,cond,locus,term)%>%
  summarise(mean=mean(V1),
            sd=sd(V1)/sqrt(2))

plot1=ggplot(df_means, aes(y=mean,x=time,colour=cond))+
  facet_wrap(~factor(locus,levels=unique(df_$locus)), scales="free_y",ncol =3,dir="h",strip.position = "top")+
  geom_point(size=1.5)+
  geom_line()+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), size = 0.5, width = 0)+
  scale_colour_manual("",labels =c("Mature","Senescent"),values=c("#197000","#BD9B00"))+
  scale_x_continuous(breaks=c(seq(48,93,12)))+
  labs(x = "\nTime (ZT)", y="Mean expression (FPKM)\n")+
  theme_bw()+
  theme(axis.ticks = element_line(colour="NA"), panel.border = element_rect(colour="light grey",size=1), 
        panel.grid.major = element_line(colour="light grey", size=0.1),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"), axis.title = element_text(size = 12),
        strip.background = element_rect(colour=NA,fill=NA), strip.text = element_text(size=10,colour="black"),
        legend.text = element_text(size=12), legend.position = "top",
        plot.margin = unit(c(0,0,0,0),"cm"),plot.background = element_rect(fill="white"))+
  guides(colour = guide_legend(override.aes = list(size = 2.5)))

gt1 = ggplot_gtable(ggplot_build(plot1))

#A bit of manipulation to get horizontal strip labels for the GO terms. Need to plot a different figure for this purpose.
df_$locus=rep(seq(1,3,1),each=96)
df_means = df_ %>%
  group_by(time,cond,locus,term)%>%
  summarise(mean=mean(V1),
            sd=sd(V1))

plot2=ggplot(df_means, aes(y=mean,x=time,colour=cond))+
  facet_grid(term~locus, scales="free_y")+
  geom_point(size=1)+
  geom_line()+
  scale_colour_manual("",labels =c("Mature","Senescent"),values=c("#197000","#BD9B00"), guide="none")+
  #scale_shape_discrete("",labels=)+
  labs(x = "\nTime (ZT)", y="Mean expression (FPKM)\n")+
  theme_bw()+
  theme(axis.ticks = element_line(colour="NA"), panel.border = element_rect(colour="light grey",size=1), 
        panel.grid.major = element_line(colour="light grey", size=0.1),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"), axis.title = element_text(size = 12,colour="black"),
        strip.background = element_rect(colour=NA,fill="light grey"), strip.text.y = element_text(size=11, angle=270,colour="black"),
        legend.key.size = unit(0.75, "cm"), legend.text = element_text(size=9),plot.background = element_rect(fill="white"))+
  guides(shape = guide_legend(override.aes = list(size = 3)))

gt2 = ggplot_gtable(ggplot_build(plot2))

gt.side1 = gtable_filter(gt2, 'strip-r-1')
gt.side2 = gtable_filter(gt2, 'strip-r-2')
gt.side3 = gtable_filter(gt2, 'strip-r-3')

gt1 = gtable_add_cols(gt1, widths=gt.side1$widths[1], pos = -1)
gt1 = gtable_add_grob(gt1, zeroGrob(), t = 1, l = ncol(gt1), b=nrow(gt1))

panel_id <- gt1$layout[grep('panel-.+1$', gt1$layout$name),]
gt1 = gtable_add_grob(gt1, gt.side1, t = panel_id$t[1], l = ncol(gt1))
gt1 = gtable_add_grob(gt1, gt.side2, t = panel_id$t[2], l = ncol(gt1))
gt1 = gtable_add_grob(gt1, gt.side3, t = panel_id$t[3], l = ncol(gt1))

grid.newpage()
plot(gt1)

#Next, plot the proportions of genes in these GO terms that fall in short, long, unchanged
#rhythmic process
rhythmic=subset(diffs_df, ID %in% rhythmic)
rhythmic$category = NA
rhythmic[rhythmic$period_diff > 0.5,]$category = "> 0.5 h"
#rhythmic[rhythmic$period_diff < -0.5,]$category = "< -0.5 h" #there are no long period genes in this category
rhythmic[is.na(rhythmic$category),]$category= "< 0.5 h & > -0.5 h"
#carbohydrate metabolic process
carbohydrate=subset(diffs_df, ID %in% carbohydrate)
carbohydrate$category = NA
carbohydrate[carbohydrate$period_diff > 0.5,]$category = "> 0.5 h"
carbohydrate[carbohydrate$period_diff < -0.5,]$category = "< -0.5 h"
carbohydrate[is.na(carbohydrate$category),]$category= "< 0.5 h & > -0.5 h"
#photosynthesis, light harvesting
photo=subset(diffs_df, ID %in% photo)
photo$category = NA
photo[photo$period_diff > 0.5,]$category = "> 0.5 h"
photo[photo$period_diff < -0.5,]$category = "< -0.5 h"
photo[is.na(photo$category),]$category= "< 0.5 h & > -0.5 h"
#all transcripts
diffs_df$category = NA
diffs_df[diffs_df$period_diff > 0.5,]$category = "> 0.5 h"
diffs_df[diffs_df$period_diff < -0.5,]$category = "< -0.5 h"
diffs_df[is.na(diffs_df$category),]$category= "< 0.5 h & > -0.5 h"

#Order the categories within the plot
diffs_df$which="a"
carbohydrate$which="b"
rhythmic$which="c"
photo$which="d"

#Proportion plot
ggplot(rbind(rhythmic,diffs_df,photo,carbohydrate),aes(x=which,fill=category))+
  geom_bar(stat="count",position="fill",width=0.75)+
  geom_text(data.frame(label=c("47.4%","52.6%","48.0%","26.8%","25.2%","55.8%","25.6%","18.6%","19.1%","12.8%","68.1%"),
                       position=c(0.237,0.737,0.24,0.614,0.874,0.279,0.686,0.907,0.096,0.255,0.659),
                       which=c(rep("c",2),rep("a",3),rep("b",3),rep("d",3)), category=rep("> 0.5 h",11)),
            mapping=aes(x=which,y=position,label=label),colour=c("black","black",rep(c("black","black","white"),3)))+
  scale_fill_viridis("Period difference",discrete=T)+
  labs(x="",y="Proportion (%)\n")+
  scale_y_continuous(expand=c(0,0))+
  scale_x_discrete(labels=c("All","carbohydrate\nmetabolic process\n(GO:0005975)","circadian rhythm\n(GO:0007623)",
                            "photosynthesis,\nlight harvesting\n(GO:0009765)"))+
  theme_bw()+
  theme(panel.background = element_blank(),axis.line = element_line(size = 0.5, colour = "black"), 
        axis.ticks = element_line(colour="NA"), panel.border = element_blank(), panel.grid = element_blank(),
        axis.text.x = element_text(size = 10, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"), axis.title = element_text(size = 12), 
        legend.text = element_text(size=10), legend.title = element_text(size=12))


#####Differential rhythmicity and differential expression using package LimoRhyde#####
#Code adapted from vignette available at https://github.com/hugheylab/limorhyde

#Define some parameters for LimoRhyde
period = 24
qvalRhyCutoff = 0.05
qvalDrCutoff = 0.05

#Decompose timepoint into sin and cos
metadata = cbind(metadata, limorhyde(metadata$time, 'time_'))


#LimoRhyde independent analysis of rhythmicity
#Please note, this analysis of rhythmicity was not used for this study, but was necessary to calculate for the subsequent calculation
#of 'differential rhythmicity'. Only BIO_CYCLE and MetaCycle were used to analyse rhythmicity of transcripts in this study.

rhyLimma = foreach(condNow = unique(metadata$cond), .combine = rbind) %do% {
  design = model.matrix(~ time_cos + time_sin, data = metadata[cond == condNow])
  vm= voom(all[, metadata$cond == condNow],design,plot=T)
  assign(paste("vm_rhy",condNow,sep=""),vm)
  fit = lmFit(vm, design)
  fit = eBayes(fit, trend = TRUE)
  rhyNow = data.table(topTable(fit, coef = 2:3, number = Inf), keep.rownames = TRUE)
  setnames(rhyNow, 'rn', 'gene_id')
  rhyNow[, cond := condNow]}

norm_rhy = cbind(vm_rhymat$E,vm_rhysen$E) # voom gives normalised counts!

rhyLimmaSummary = rhyLimma[, .(P.Value = min(P.Value)), by = gene_id]
rhyLimmaSummary[, adj.P.Val := p.adjust(P.Value, method = 'BH')]
setorderv(rhyLimmaSummary, 'adj.P.Val')


#LimoRhyde assessment of differential rhythmicity

design = model.matrix(~ cond * (time_cos + time_sin), data = metadata)
vm = voom(all,design,plot=T)
fit = lmFit(vm, design)
fit = eBayes(fit, trend = TRUE)
drLimma = data.table(topTable(fit, coef = 5:6, number = Inf), keep.rownames = TRUE)
setnames(drLimma, 'rn', 'gene_id')
drLimma = drLimma[gene_id %in% rhyLimmaSummary[adj.P.Val <= qvalRhyCutoff]$gene_id]
drLimma[, adj.P.Val := p.adjust(P.Value, method = 'BH')]
setorderv(drLimma, 'adj.P.Val')

norm_dr = vm$E

#LimoRhyde assessment of differential expression

design = model.matrix(~ cond + time_cos + time_sin, data = metadata)
#vm=voom(all,design,plot=T)
fit = lmFit(all, design)
fit = eBayes(fit, trend = TRUE)
deLimma = data.table(topTable(fit, coef = 2, number = Inf), keep.rownames = TRUE)
setnames(deLimma, 'rn', 'gene_id')
de_dr = deLimma[(gene_id %in% drLimma[adj.P.Val <= qvalDrCutoff]$gene_id)]
deLimma = deLimma[!(gene_id %in% drLimma[adj.P.Val <= qvalDrCutoff]$gene_id)]
deLimma[, adj.P.Val := p.adjust(P.Value, method = 'BH')]
setorderv(deLimma, 'adj.P.Val')

norm_de = vm$E


#Plot a comparison of expression of the top result of rhythmicity, differential rhythmicity and differential expression to illustrate
#differences between the statistical techniques

geneIdsNow = c(rhyLimmaSummary$gene_id[1L], drLimma$gene_id[1L], deLimma$gene_id[1L])
df = data.table(cbind(norm_rhy[geneIdsNow[1],],norm_dr[geneIdsNow[2],],norm_de[geneIdsNow[3],]))
df[, sample_id := colnames(all[geneIdsNow, ])]
df = merge(df, metadata[, .(sample_id, cond, time)], by = 'sample_id')

df_=rbind(df[,-c(3,4)], df[,-c(2,4)],df[,-c(2,3)], use.names=FALSE)
df_$name=c(rep("Rhythmicity",96),rep("Differential\nrhythmicity",96),rep("Differential\nexpression",96))
df_$name=factor(df_$name,levels=c("Rhythmicity","Differential\nrhythmicity","Differential\nexpression"))
df_$cond=rep(c(rep("Mature",48),rep("Senescent",48)),3)

df_means = df_ %>%
  group_by(time,cond,name)%>%
  summarise(mean=mean(V1),
            sd=sd(V1))

ggplot(df_means, aes(y=mean,x=time,colour=cond, shape=name))+
  facet_grid(name~cond, scales="free_y")+
  geom_point(size=2)+
  geom_line()+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), size = 0.5, width = 0)+
  scale_colour_manual("",labels =c("Mature","Senescent"),values=c("#197000","#BD9B00"), guide="none")+
  scale_shape_discrete("",labels=geneIdsNow)+
  labs(x = "\nTime (ZT)", y="Mean expression (FPKM)\n")+
  theme_bw()+
  theme(axis.ticks = element_line(colour="NA"), panel.border = element_rect(colour="light grey",size=1), 
        panel.grid.major = element_line(colour="light grey", size=0.1),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 16, colour = "black"),
        axis.text.y = element_text(size = 16, colour = "black"), axis.title = element_text(size = 18),
        strip.background = element_rect(colour=NA,fill="light grey"), strip.text = element_text(size=16),
        legend.key.size = unit(0.75, "cm"), legend.text = element_text(size=16), legend.position="top")+
  guides(shape = guide_legend(override.aes = list(size = 3)))+
  scale_x_continuous(breaks=c(seq(48,84,12)))

#####Further investigation of differentially rhythmic transcripts#####

#Isolate significantly differentially rhythmic genes (q < 0.05)
dr_genes=filter(drLimma,adj.P.Val < 0.05)
#Subset these to retain only transcripts also identified as rhythmic in both conditions by BIO_CYCLE and MetaCycle
dr_diffs = subset(diffs_df, ID %in% dr_genes$gene_id)

#Do a GO term enrichment of these DR transcripts, as before
tmp <- ifelse(go_terms$Gene.stable.ID %in% dr_diffs$ID, 1, 0)
genelist_dr=tmp
names(genelist_dr)=go_terms$Gene.stable.ID

GOdata_dr <- new("topGOdata",
                 ontology = "BP",
                 allGenes = genelist_dr,
                 geneSelectionFun = function(x)(x == 1),
                 annot = annFUN.gene2GO, gene2GO = gene2GO)


Fisher_dr <- runTest(GOdata_dr, algorithm = "classic", statistic = "fisher")

go_results_dr <- GenTable(GOdata_dr, raw.p.value = Fisher_dr, topNodes = length(Fisher_dr@score))

go_results_dr$raw.p.value = gsub("e", "E", go_results_dr$raw.p.value)
go_results_dr$raw.p.value = gsub("< ", "", go_results_dr$raw.p.value)
go_results_dr$raw.p.value = as.numeric(go_results_dr$raw.p.value)

go_results_dr$Observed.Expected = go_results_dr$Significant/go_results_dr$Expected
go_results_dr$Stage = "DR"

#Curate the list to remove duplicate terms
#List may vary slightly run to run
g=go_results_dr[c(1,3:12,14,17,19,23:25,26,27,29),]
g$Term =with(g, reorder(Term , Observed.Expected, mean , na.rm=T))

ggplot(g,aes(x=log(Observed.Expected),y=Term,fill=-log(raw.p.value)))+
  geom_bar(stat="identity",width=0.15,colour="transparent")+
  geom_point(g,mapping=aes(log(Observed.Expected),Term,fill=-log(raw.p.value),size=Significant),pch=21,stroke=0,show.legend = T)+
  scale_fill_gradientn(expression(paste("-log(",italic(p),")")),colours=c("#f7d7a6","#fcd69d","darkorange","red","red1"),values=c(0,0.15,0.4,0.75,1),limits=c(7,41),
                       breaks=seq(8,43,8))+
  scale_size_continuous("Gene count",range=c(3,8),breaks=c(5,25,50,100))+
  labs(x="\nlog(Observed/Expected)",y="")+
  theme_bw()+
  theme(panel.background = element_blank(),axis.line = element_line(size = 0.5, colour = "black"), 
        axis.ticks = element_line(colour="NA"), panel.border = element_blank(), panel.grid = element_blank(),
        axis.text.x = element_text(size = 10, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"), axis.title = element_text(size = 12),
        legend.text = element_text(size=10), legend.title = element_text(size=12))
  guides(size=guide_legend(override.aes = list(fill="black")))

#Check the circadian rhythm parameter changes of DR genes
  
ggplot()+  
  geom_density(diffs_df, mapping=aes(x=period_diff), colour="black",fill="black",alpha=0.5)+
  geom_density(dr_diffs, mapping=aes(x=period_diff), colour="blue",fill="blue",alpha=0.5)+
  labs(y="Density\n", x="\nPeriod difference (h)")+
  scale_y_continuous(expand=c(0,0))+
  theme_bw()+
  theme(panel.background = element_blank(),axis.line = element_line(size = 0.5, colour = "black"), 
        axis.ticks = element_line(colour="NA"), panel.border = element_blank(), panel.grid = element_blank(),
        axis.text.x = element_text(size = 10, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"), axis.title = element_text(size = 12),
        legend.text = element_text(size=10), legend.title = element_text(size=12),
        strip.background = element_blank(), strip.text = element_text(size=10))

ggplot()+
  geom_density(diffs_df, mapping=aes(x=phase_diff), colour="black",fill="black",alpha=0.5)+
  geom_density(dr_diffs, mapping=aes(x=phase_diff), colour="blue",fill="blue",alpha=0.5)+
  scale_x_continuous(limits=c(-8,8))+
  labs(y="", x="\nPhase difference (h)")+
  scale_y_continuous(expand=c(0,0))+
  theme_bw()+
  theme(panel.background = element_blank(),axis.line = element_line(size = 0.5, colour = "black"), 
        axis.ticks = element_line(colour="NA"), panel.border = element_blank(), panel.grid = element_blank(),
        axis.text.x = element_text(size = 10, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"), axis.title = element_text(size = 12),
        legend.text = element_text(size=10), legend.title = element_text(size=12),
        strip.background = element_blank(), strip.text = element_text(size=10))

ggplot(rbind(mutate(dr_diffs,Type="Differentially rhythmic"),mutate(diffs_df,Type="All")),aes(x=amp_diff,fill=Type,colour=Type))+
  geom_density(alpha=0.5)+
  scale_x_continuous(limits=c(-2,2))+
  labs(y="", x="\nAmplitude change (log(Mature:senescence))")+
  scale_y_continuous(expand=c(0,0))+
  scale_fill_manual("",values=c("black","blue"))+
  scale_colour_manual("",values=c("black","blue"))+
  theme_bw()+
  theme(panel.background = element_blank(),axis.line = element_line(size = 0.5, colour = "black"), 
        axis.ticks = element_line(colour="NA"), panel.border = element_blank(), panel.grid = element_blank(),
        axis.text.x = element_text(size = 10, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"), axis.title = element_text(size = 12),
        legend.text = element_text(size=10), legend.title = element_text(size=12),
        strip.background = element_blank(), strip.text = element_text(size=10))

#Now statistically test how different the distributions are with Fisher-Pitman Permutation test
oneway_test(period_diff ~ as.factor(Term), data=rbind(mutate(dr_diffs, Term="DR"),mutate(diffs_df,Term="All")),distribution = approximate(nresample = 10000))#*
oneway_test(phase_diff ~ as.factor(Term), data=rbind(mutate(dr_diffs, Term="DR"),mutate(diffs_df,Term="All")),distribution = approximate(nresample = 10000))#***
oneway_test(amp_diff ~ as.factor(Term), data=rbind(mutate(dr_diffs, Term="DR"),mutate(diffs_df,Term="All")),distribution = approximate(nresample = 10000))#*

#Now look at DR TFs - are there any interesting candidates that might be driving changes?

#Get DR TFs
tfs_dr = subset(dr_genes, gene_id %in% tfs_Mace$Gene.stable.ID)

#Filter against transcripts that are rhythmic in both tissue types. Get circadian parameter changes
dr_diffs = subset(diffs_df, ID %in% dr_genes$gene_id)
tfs_diffs = subset(diffs_df, ID %in% tfs_Mace$Gene.stable.ID)
tfs_dr_diffs = subset(diffs_df, ID %in% tfs_dr$gene_id)

#Differentially rhythmic TFs are more advanced in phase than all differentially rhythmic
#and much more advanced than all rhythmic
mean(filter(tfs_dr_diffs, phase_diff > -12 & phase_diff < 12)$phase_diff)
mean(filter(dr_diffs, phase_diff > -12 & phase_diff < 12)$phase_diff)
mean(filter(tfs_diffs, phase_diff > -12 & phase_diff < 12)$phase_diff)
mean(filter(diffs_df, phase_diff > -12 & phase_diff < 12)$phase_diff)

#Smaller differences between period and amplitude. Period change seemingly reduced for DR and DR TFs
mean(tfs_dr_diffs$period_diff)
mean(dr_diffs$period_diff)
mean(tfs_diffs$period_diff)
mean(diffs_df$period_diff)

mean(tfs_dr_diffs$amp_diff)
mean(dr_diffs$amp_diff)
mean(tfs_diffs$amp_diff)
mean(diffs_df$amp_diff)

#Statistical tests of these parameter differences. Test within pairs, i.e. DR TFs vs TFs, DR vs All
t.test(tfs_dr_diffs$phase_diff,tfs_diffs$phase_diff)
t.test(diffs_df$phase_diff,dr_diffs$phase_diff)

t.test(tfs_dr_diffs$period_diff,tfs_diffs$period_diff)
t.test(diffs_df$period_diff,dr_diffs$period_diff)

#Visualise these differences by plotting some boxplots

#Combine the data from the different genesets
tfs_phase_comp=rbind(mutate(tfs_diffs,Type="Rhythmic\nTFs"),mutate(tfs_dr_diffs,Type="Differentially\nrhythmic\nTFs"),mutate(dr_diffs,Type="Differentially\nrhythmic"),mutate(diffs_df,Type="All\nrhythmic"))
tfs_phase_comp$Type=factor(tfs_phase_comp$Type, levels=c("All\nrhythmic","Rhythmic\nTFs","Differentially\nrhythmic",
                                                         "Differentially\nrhythmic\nTFs"))


#Phase
ggplot(tfs_phase_comp,aes(y=phase_diff,x=Type,colour=Type,fill=Type))+
  geom_jitter(size=1.5,alpha=0.5,width=0.25)+
  geom_boxplot(outlier.shape = NA,colour="black",alpha=0.5)+
  labs(x="", y="Phase difference (h)\n")+
  scale_fill_viridis(discrete=T,begin=0.15,end=0.85)+
  scale_colour_viridis(discrete=T,begin=0.15,end=0.85)+
  theme_bw()+
  theme(panel.background = element_blank(),axis.line = element_line(size = 0.5, colour = "black"), 
        axis.ticks = element_line(colour="NA"), panel.border = element_blank(), panel.grid = element_blank(),
        axis.text.x = element_text(size = 10, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"), axis.title = element_text(size = 12),
        legend.text = element_text(size=10), legend.title = element_text(size=12),legend.position = "none")

#Period
ggplot(tfs_phase_comp,aes(y=period_diff,x=Type,colour=Type,fill=Type))+
  geom_jitter(size=1.5,alpha=0.5,width=0.25)+
  geom_boxplot(outlier.shape = NA,colour="black",alpha=0.5)+
  labs(x="", y="Period difference (h)\n")+
  scale_fill_viridis(discrete=T,begin=0.15,end=0.85)+
  scale_colour_viridis(discrete=T,begin=0.15,end=0.85)+
  theme_bw()+
  theme(panel.background = element_blank(),axis.line = element_line(size = 0.5, colour = "black"), 
        axis.ticks = element_line(colour="NA"), panel.border = element_blank(), panel.grid = element_blank(),
        axis.text.x = element_text(size = 10, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"), axis.title = element_text(size = 12),
        legend.text = element_text(size=10), legend.title = element_text(size=12),legend.position = "none")


#Check whether any family of TFs has significantly more members in the set of DR TFs than would be expected by chance
tfs_dr_diffs = merge(tfs_dr_diffs, tfs_Mace, by.x="ID",by.y="Gene.stable.ID")

#Get proportions of each TF family in DR TFs set
proportion_dr = tfs_dr_diffs %>% group_by(superfamily)%>%
  summarise(count=n(),
            proportion=(n()/nrow(tfs_dr_diffs))*100)

#Need to compare these proportions against proportions of each TF family among all transcripts rhythmic in both tissues
#Get superfamily info for tfs_diffs
tfs_diffs=merge(tfs_diffs, tfs_Mace[,c(2,7)], by.x="ID",by.y="Gene.stable.ID")
proportion_all = tfs_diffs %>% group_by(superfamily)%>%
  summarise(count=n(),
            proportion=(n()/nrow(tfs_diffs))*100)

#Make a contingency table and run a Fisher's exact test
fishers_tf = proportion_dr
fishers_tf$P.value=NA

for (i in 1:nrow(proportion_dr)){
  #fishers
  x=proportion_dr$count[i] #superfamily in dr
  y=proportion_all$count[which(proportion_all$superfamily==proportion_dr$superfamily[i])]#superfamily in all
  
  contingency.table <- matrix(c(x, sum(proportion_dr$count)-x, y, sum(proportion_all$count)-y), nrow = 2)
  fishers_tf$P.value[i]=fisher.test(contingency.table)$p.value
  #fold enrichment
  a=proportion_dr$proportion[i]
  b=proportion_all$proportion[which(proportion_all$superfamily==proportion_dr$superfamily[i])]
  fishers_tf$proportion_diff[i]=log(a/b)
}

#Reorder by proportion difference
fishers_tf$superfamily =with(fishers_tf, reorder(superfamily , proportion_diff, mean , na.rm=T))  

#Exclude if family count is just 1 as this is a bit meaningless
fishers_tf=filter(fishers_tf, count > 1)

ggplot(fishers_tf,aes(x=proportion_diff,y=superfamily,fill=-log(P.value)))+
  geom_bar(stat="identity",width=0.15,colour="transparent")+
  geom_point(fishers_tf,mapping=aes(proportion_diff,superfamily,fill=-log(P.value),size=count),pch=21,stroke=0,show.legend = T)+
  scale_fill_gradientn(expression(paste("-log(",italic(p),")")),colours=c("blue","blue","purple","red","red"),values=c(0,0.1,0.5,0.9,1))+
  scale_size_continuous("Gene count",range=c(3,8))+
  labs(x="\nlog(Fold enrichment)",y="")+
  theme_bw()+
  theme(panel.background = element_blank(),axis.line = element_line(size = 0.5, colour = "black"), 
        axis.ticks = element_line(colour="NA"), panel.border = element_blank(), panel.grid = element_blank(),
        axis.text.x = element_text(size = 10, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"), axis.title = element_text(size = 12),
        legend.text = element_text(size=10), legend.title = element_text(size=12))+
  guides(size=guide_legend(override.aes = list(fill="black")))


#####String DB protein-protein network#####
#Rationale was to assess whether there is evidence that this set of phase-advanced DR TFs work in a network

#I extracted a list of DR TFs, converted their Mace IDs to CS IDs (where one-to-one orthology was available)
#Then downloaded sequences of these genes and BLASTed them in STRING to find the old CSS (2014) IDs that STRING uses
#See methods for further details

#Import resulting network map files from STRING
mapping=read.csv("string_mapping.csv")
mapping$stringId=substr(mapping$stringId, 1,nchar(mapping$stringId)-2)
#Add TF superfamily info
mapping=merge(mapping,tfs,by.x="queryItem",by.y="locus",all.x=T)

coords=read.csv("string_network_coordinates.csv")
graphs=read.csv("string_interactions_short.csv")

#Rejig the dataframes to make the network plottable
graphs$order=1:614

node1=merge(graphs, coords, by.x="node1_string_id", by.y="identifier")
node1=node1[order(node1$order),]
node2=merge(graphs, coords, by.x="node2_string_id", by.y="identifier")
node2=node2[order(node2$order),]

node1$number=1:614
node2$number=1:614

merged_nodes=rbind(node1,node2)
merged_nodes$number=as.factor(merged_nodes$number)
merged_nodes=merged_nodes[order(merged_nodes$number),]
merged_nodes$pch="A"

colours=rep("black",615)

#Read in some additional info about the broad processes that each network protein is involved in
convert_ids=read.csv("differentially_rhythmic_tf_candidates_IWGSC_to_CSS.csv")
#Subset
convert_ids=subset(convert_ids, CSS.ID %in% substr(merged_nodes$node1_string_id,6,(nchar(as.character(merged_nodes$node1_string_id))-2))
                   | CSS.ID %in% substr(merged_nodes$node2_string_id,6,(nchar(as.character(merged_nodes$node2_string_id))-2)))
coords$new=substr(coords$identifier,6,(nchar(as.character(coords$identifier))-2))
convert_ids=merge(convert_ids,coords,by.x="CSS.ID", by.y="new")
convert_ids$combined_score=0
convert_ids$number=as.factor(0)
convert_ids$Refined_process=as.character(convert_ids$Refined_process)
convert_ids$Refined_process[which(convert_ids$Refined_process == "")]="AA"
convert_ids$Refined_process=as.factor(convert_ids$Refined_process)
convert_ids$Genome=substr(convert_ids$Gene.stable.ID,9,9)

#Plot
ggplot(merged_nodes, aes(x=x_position,y=y_position, linewidth=combined_score,linetype=number))+
  geom_line(show.legend=T, alpha=0.6, colour="dark grey")+
  scale_linetype_manual(values=rep(1, 615), guide="none")+
  geom_point(merged_nodes, mapping=aes(x=x_position,y=y_position,pch=pch),size=28,colour="transparent",fill="transparent",stroke=0)+
  geom_point(convert_ids, mapping=aes(x=x_position,y=y_position, fill=Refined_process,colour=Refined_process),size=15.5, pch=21,stroke=0.5)+
  scale_fill_manual(values=c("grey92","#968EBA","#8BB2E1","#71CE94","#feec50"),labels=c("Undetermined",levels(convert_ids$Refined_process)[-1]))+
  scale_colour_manual(values=c("grey","#645e80","#5c7899","#457a59","#9c902f"),labels=c("Undetermined",levels(convert_ids$Refined_process)[-1]),guide="none")+
  scale_shape_manual(values=c(21),labels="")+
  scale_linewidth_continuous(range=c(0.15,1.5), breaks=c(0.4,0.6,0.8,1), limits=c(0.4,1))+
  scale_x_continuous(expand=c(0,0.1),limits=c(min(coords$x_position),max(coords$x_position)))+
  scale_y_continuous(limits=c(min(coords$y_position)-0.02,max(coords$y_position)+0.01))+
  geom_text(convert_ids[convert_ids$Wheat_name == "",],
            mapping=aes(x=x_position,y=y_position+0.0075, label=Superfamily),show.legend=F,size=2.25,colour="black")+
  geom_text(convert_ids[convert_ids$Wheat_name == "",],
            mapping=aes(x=x_position,y=y_position-0.015, label=Genome),show.legend=F,size=2.25,colour="black")+
  geom_text(convert_ids[convert_ids$Wheat_name != "",],
            mapping=aes(x=x_position,y=y_position+0.015, label=Superfamily),show.legend=F,size=2.25,colour="black")+
  geom_text(convert_ids[convert_ids$Wheat_name != "",],
            mapping=aes(x=x_position,y=y_position, label=Wheat_name),show.legend=F,size=2.25, fontface="bold",colour="black")+
  geom_text(convert_ids[convert_ids$Wheat_name != "",],
            mapping=aes(x=x_position,y=y_position-0.015, label=Genome),show.legend=F,size=2.25,colour="black")+
  guides(linewidth = guide_legend(override.aes = list(size=8,label=""), title="Combined score"))+
  geom_text(convert_ids,#dummy!
            mapping=aes(x=2,y=2, label=Genome),show.legend=T,size=2.9,colour="black")+
  guides(linewidth = guide_legend(override.aes = list(size=8,label=""), title="Combined score",order=1,keywidth = 3))+
  guides(fill = guide_legend(override.aes = list(size=12,stroke=0.5,label="",colour=c("grey","#645e80","#5c7899","#457a59","#9c902f"),
                                                 linetype="blank"),title="Related process",order=2,title.hjust = 0))+
  guides(shape = guide_legend(override.aes = list(fill="grey92",
                                                stroke=1, colour="grey",pch=21,label="",
                                                linetype="blank"),title="Annotation",order=3, title.hjust=0.3))+
  theme_bw()+
  theme(axis.text = element_blank(),axis.line = element_blank(), axis.title = element_blank(),
        axis.ticks = element_blank(), panel.grid = element_blank(),legend.text=element_text(size=10),
        legend.title=element_text(size=12),panel.border = element_blank())

#Analysis of putative targets of the subnetwork of DR TFs

#GENIE3 GRN Analysis
#This is the process for using GENIE3 for predicting all putative interactions of genes that are rhythmic in both mature and senescent tissue
#I decided to only focus on the interactions of these genes in senescent tissue, however.

#The analysis itself is commented out as it is computationally expensive

#set.seed(123)
#sen_matrix=as.matrix(sen_shared[,c(1,9:56)])
#rownames(sen_matrix)=sen_matrix[,1]
#sen_matrix=sen_matrix[,-1]
#genie_senescent=GENIE3(sen_matrix)
#write.csv(genie_senescent, "GENIE3_Output_Matrix_Senescent.csv", q=F)

#Get all links
#links_senescent=getLinkList(genie_senescent)

#Keep top 1% of links for rigour
#toplinks=head(links_senescent,n=(nrow(links_senescent)*0.01))
#write.csv(toplinks, "GENIE3_sen_rhy_toplinks.csv", q=F,row.names=F)

#Read in toplinks file
toplinks=read.csv("GENIE3_sen_rhy_toplinks.csv")

#Keep links that are associated with DR TFs (i.e. the protein subnetwork)
DR_links=subset(toplinks,regulatoryGene %in% tfs_dr_diffs$ID)

#Check circadian parameter changes of the target genes
mean(diffs_df[diffs_df$ID %in% unique(DR_links$targetGene),]$phase_diff)
mean(diffs_df$phase_diff)
mean(diffs_df[diffs_df$ID %in% unique(DR_links$targetGene),]$period_diff)
mean(diffs_df$period_diff)
#Parameter changes are nearly identical to the global means 

#Plot these comparisons
DR_links_df=diffs_df[diffs_df$ID %in% unique(DR_links$targetGene),]

ggplot()+  
  geom_density(diffs_df, mapping=aes(x=period_diff), colour="black",fill="black",alpha=0.5)+
  geom_density(DR_links_df, mapping=aes(x=period_diff), colour="blue",fill="blue",alpha=0.5)+
  geom_vline(xintercept=mean(diffs_df$period_diff),linetype="dashed")+
  labs(y="Density\n", x="\nCircadian period difference (h)")+
  scale_y_continuous(expand=c(0,0))+
  theme_bw()+
  theme(panel.background = element_blank(),axis.line = element_line(size = 0.5, colour = "black"), 
        axis.ticks = element_line(colour="NA"), panel.border = element_blank(), panel.grid = element_blank(),
        axis.text.x = element_text(size = 10, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"), axis.title = element_text(size = 12),
        legend.text = element_text(size=10), legend.title = element_text(size=12),
        strip.background = element_blank(), strip.text = element_text(size=10))

ggplot()+  
  geom_density(diffs_df, mapping=aes(x=phase_diff), colour="black",fill="black",alpha=0.5)+
  geom_density(DR_links_df, mapping=aes(x=phase_diff), colour="blue",fill="blue",alpha=0.5)+
  geom_vline(xintercept=mean(diffs_df$phase_diff),linetype="dashed")+
  labs(y="Density\n", x="\nCircadian phase difference (h)")+
  scale_y_continuous(expand=c(0,0))+
  theme_bw()+
  theme(panel.background = element_blank(),axis.line = element_line(size = 0.5, colour = "black"), 
        axis.ticks = element_line(colour="NA"), panel.border = element_blank(), panel.grid = element_blank(),
        axis.text.x = element_text(size = 10, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"), axis.title = element_text(size = 12),
        legend.text = element_text(size=10), legend.title = element_text(size=12),
        strip.background = element_blank(), strip.text = element_text(size=10))

#Do a GO enrichment with the DR links

tmp <- ifelse(go_terms$Gene.stable.ID %in% DR_links$targetGene, 1, 0)
genelist_dr=tmp
names(genelist_dr)=go_terms$Gene.stable.ID

GOdata_dr <- new("topGOdata",
                 ontology = "BP",
                 allGenes = genelist_dr,
                 geneSelectionFun = function(x)(x == 1),
                 annot = annFUN.gene2GO, gene2GO = gene2GO)

Fisher_dr <- runTest(GOdata_dr, algorithm = "classic", statistic = "fisher")

go_results_dr <- GenTable(GOdata_dr, raw.p.value = Fisher_dr, topNodes = length(Fisher_dr@score))

go_results_dr$raw.p.value = gsub("e", "E", go_results_dr$raw.p.value)
go_results_dr$raw.p.value = gsub("< ", "", go_results_dr$raw.p.value)
go_results_dr$raw.p.value = as.numeric(go_results_dr$raw.p.value)

go_results_dr$Observed.Expected = go_results_dr$Significant/go_results_dr$Expected
go_results_dr$Stage = "DR"

#Curate the list to remove duplicate terms

g=go_results_dr[c(1,2,3,4,5,6,9,10,11,12,13,15,16,18,19,20,23,24,25,27),]

g$Term =with(g, reorder(Term , Observed.Expected, mean , na.rm=T))

ggplot(g,aes(x=log(Observed.Expected),y=Term,fill=-log(raw.p.value)))+
  geom_bar(stat="identity",width=0.15,colour="transparent")+
  geom_point(g,mapping=aes(log(Observed.Expected),Term,fill=-log(raw.p.value),size=Significant),pch=21,stroke=0,show.legend = T)+
  scale_fill_gradientn(expression(paste("-log(",italic(p),")")),colours=c("#f7d7a6","#fcd69d","darkorange","red","red1"),values=c(0,0.15,0.4,0.75,1),limits=c(20,70),
                       breaks=seq(20,70,10))+
  scale_size_continuous("Gene count",range=c(3,8),breaks=c(5,25,50,100))+
  labs(x="\nlog(Observed/Expected)",y="")+
  theme_bw()+
  theme(panel.background = element_blank(),axis.line = element_line(size = 0.5, colour = "black"), 
        axis.ticks = element_line(colour="NA"), panel.border = element_blank(), panel.grid = element_blank(),
        axis.text.x = element_text(size = 10, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"), axis.title = element_text(size = 12),
        legend.text = element_text(size=10), legend.title = element_text(size=12))+
  guides(size=guide_legend(override.aes = list(fill="black")))

#####Further analysis of differentially expressed transcripts#####

#First, upregulated genes
deup_genes=filter(deLimma,adj.P.Val < 0.05 & logFC > 1.2)

#GO enrichment analysis
tmp <- ifelse(go_terms$Gene.stable.ID %in% deup_genes$gene_id, 1, 0)
genelist_deup=tmp
names(genelist_deup)=go_terms$Gene.stable.ID

GOdata_deup <- new("topGOdata",
                   ontology = "BP",
                   allGenes = genelist_deup,
                   geneSelectionFun = function(x)(x == 1),
                   annot = annFUN.gene2GO, gene2GO = gene2GO)

Fisher_deup <- runTest(GOdata_deup, algorithm = "classic", statistic = "fisher")

go_results_deup <- GenTable(GOdata_deup, raw.p.value = Fisher_deup, topNodes = length(Fisher_deup@score))

go_results_deup$raw.p.value = gsub("e", "E", go_results_deup$raw.p.value)
go_results_deup$raw.p.value = gsub("< ", "", go_results_deup$raw.p.value)
go_results_deup$raw.p.value = as.numeric(go_results_deup$raw.p.value)

#Second, downregulated genes
dedown_genes=filter(deLimma,adj.P.Val < 0.05 & logFC < -1.2)

#GO enrichment analysis
tmp <- ifelse(go_terms$Gene.stable.ID %in% dedown_genes$gene_id, 1, 0)
genelist_dedown=tmp
names(genelist_dedown)=go_terms$Gene.stable.ID

GOdata_dedown <- new("topGOdata",
                     ontology = "BP",
                     allGenes = genelist_dedown,
                     geneSelectionFun = function(x)(x == 1),
                     annot = annFUN.gene2GO, gene2GO = gene2GO)


Fisher_dedown <- runTest(GOdata_dedown, algorithm = "classic", statistic = "fisher")

go_results_dedown <- GenTable(GOdata_dedown, raw.p.value = Fisher_dedown, topNodes = length(Fisher_dedown@score))

go_results_dedown$raw.p.value = gsub("e", "E", go_results_dedown$raw.p.value)
go_results_dedown$raw.p.value = gsub("< ", "", go_results_dedown$raw.p.value)
go_results_dedown$raw.p.value = as.numeric(go_results_dedown$raw.p.value)

go_results_deup$Stage = "DE_UP"
go_results_dedown$Stage = "DE_DOWN"

#Curate the lists to remove duplicate terms
#Precise list order might vary slightly from run-to-run
deup_head=go_results_deup[c(1,4,5,6,7,8,9,10,11,15,18,21,22,26,28),]
dedown_head=go_results_dedown[c(1,2,4,5,7,8,9,13,16,17,18,19,20,22,26),]

#Combine into one
combine=rbind(deup_head,dedown_head)
combine$proportion=combine$Significant/combine$Expected

#Trick to make a gap in the plot
combine[31,]=0
combine$Stage[31]="DE_DOWN"
combine$Term[31]=""
combine$Significant[31]=NA
combine$proportion[31]=NA

#Reorder the levels for the plot
combine$Term=factor(combine$Term, levels=(arrange(combine,Stage,proportion)$Term))

ggplot(combine,aes(x=proportion,y=Term,fill=-log(raw.p.value),colour=Stage))+
  geom_bar(stat="identity",width=0.35,linewidth=0.75)+
  geom_point(combine,mapping=aes(proportion,Term,fill=-log(raw.p.value),size=Significant),pch=21,stroke=0.75,show.legend = T)+
  scale_fill_gradientn(expression(paste("-log(",italic(p),")")),colours=c("grey25","grey30","grey60","grey90"),values=c(0,0.1,0.5,1))+
  scale_size_continuous("Gene count",range=c(2,6.5),breaks=c(50,250,500,1000))+
  scale_color_manual("",values=c("dark blue","dark red"),labels=c("Upregulated","Downregulated"))+
  scale_x_continuous(expand=c(0,0),limits=c(0,12))+
  labs(x="\nObserved/Expected",y="")+
  theme_bw()+
  theme(panel.background = element_blank(),axis.line = element_line(size = 0.5, colour = "black"), 
        axis.ticks = element_line(colour=NA), panel.border = element_blank(), panel.grid = element_blank(),
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"), axis.title = element_text(size = 14),
        legend.text = element_text(size=12), legend.title = element_text(size=14))+
  guides(size=guide_legend(override.aes = list(fill="black")))+
  guides(colour=guide_legend(override.aes = list(fill="transparent",colour=c("dark red","dark blue"))))


#Now make a simple log fold and p-value plot

#Categorise into significance bands
deLimma$cat="Not significant\n"
deLimma[deLimma$adj.P.Val < 0.05 & deLimma$logFC > 1.2,]$cat="Up-regulated\n(n=6310)\n"
deLimma[deLimma$adj.P.Val < 0.05 & deLimma$logFC < -1.2,]$cat="Down-regulated\n(n=8785)\n"

#Label the RBCS genes that are separated on the plot
deLimma$RBCS=""
deLimma[deLimma$logFC < -4000 & -log10(deLimma$adj.P.Val) > 8 ,]$RBCS="RBCS"

#Outer plot
p1=ggplot(deLimma, aes(x=logFC,y=-log10(adj.P.Val),colour=cat))+
  geom_point(alpha=0.5,show.legend = T,size=1)+
  stat_ellipse(deLimma[deLimma$logFC < -4000 & -log10(deLimma$adj.P.Val) > 8 ,],mapping=aes(x=logFC,y=-log10(adj.P.Val)),colour="dark blue",show.legend = F)+
  scale_colour_manual(values=c("dark blue","grey15","dark red"))+
  geom_text(mapping=aes(x=-4500,y=30,label="RBCS"),size=3,colour="dark blue",show.legend=F,hjust=0)+
  geom_rect(aes(xmin = -500, xmax = 500, ymin = 0, ymax = 65), color = "black", alpha = 0,size=0.5)+
  geom_path(aes(x,y,group=grp), 
            data=data.frame(x = c(-500,-7800,-1250,-500), y=c(0,35,70,65),grp=c(1,1,2,2),cat=c(rep("",4))),
            linetype='dashed',colour="black",size=0.5)+
  labs(y=expression(paste("-log10(",italic(p),")","\n")), x="\nlog2(Fold change)")+
  scale_y_continuous(expand=c(0,0))+
  theme_bw()+
  theme(panel.background = element_blank(),axis.line = element_line(size = 0.5, colour = "black"), 
        axis.ticks = element_line(colour="NA"), panel.border = element_blank(), panel.grid = element_blank(),
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"), axis.title = element_text(size = 14),
        legend.text = element_text(size=12), legend.title = element_text(size=14))+
  guides(colour=guide_legend(override.aes = list(linetype="blank",size=2,label=""),title=""))

#Window plot
p2=ggplot(deLimma, aes(x=logFC,y=-log10(adj.P.Val),colour=cat))+
  geom_point(alpha=0.5,show.legend = F,size=0.5)+
  scale_x_continuous(limits=c(-500,500))+
  scale_colour_manual(values=c("dark blue","grey15","dark red"))+
  labs(y=expression(paste("-log10(",italic(p),")")), x="log2(Fold change)")+
  scale_y_continuous(expand=c(0,0),limits=c(0,62))+
  theme_bw()+
  theme(panel.background = element_blank(),axis.line = element_line(size = 0.5, colour = "black"), 
        axis.ticks = element_line(colour="NA"), panel.border = element_blank(), panel.grid = element_blank(),
        axis.text.x = element_text(size = 10, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"), axis.title = element_text(size = 12),
        legend.text = element_text(size=10), legend.title = element_text(size=12))

p1 + 
  annotation_custom(ggplotGrob(p2), xmin = -7800, xmax = -1250, ymin = 35, ymax = 70) +
  geom_rect(aes(xmin = -7800, xmax = -1250, ymin = 35, ymax = 70), color='black', alpha=0,size=0.4)

#####Analysis of changes within the oscillator#####

#Import wheat clock gene IDs
clock_Mace=read.csv("Clock_gene_IDs.csv")

#Merge to get circadian parameter differences
clock_diffs=subset(diffs_df, ID %in% clock_Mace$Gene.stable.ID)
clock=merge(clock_diffs,clock_Mace,by.x="ID",by.y="Gene.stable.ID")
clock$Name=gsub("LUX BOA", "LUX/BOA", clock$Name)
clock$Name=gsub("TaPRR37","Ppd",clock$Name)

#Phase in the oscillator itself is more advanced
t.test(clock$phase_diff, diffs_df$phase_diff)

#Attempt to summarise elements of the oscillator's response at senescence
#Trying some clustering of circadian parameter changes
clock=clock[order(clock$Name),]
set.seed(2)#set seed for reproducible results
pca_clock_diffs = prcomp(scale(clock[,c(3:5)]), center=TRUE, scale. = TRUE)
#K-means clustering
#Find optimal K using the elbow plot method
km1=kmeans(pca_clock_diffs$x,1)
km2=kmeans(pca_clock_diffs$x,2)
km3=kmeans(pca_clock_diffs$x,3)
km4=kmeans(pca_clock_diffs$x,4)
km5=kmeans(pca_clock_diffs$x,5)
km6=kmeans(pca_clock_diffs$x,6)
km7=kmeans(pca_clock_diffs$x,7)
km8=kmeans(pca_clock_diffs$x,8)

ratio_ss = c(km1$tot.withinss/km1$totss,km2$tot.withinss/km2$totss, km3$tot.withinss/km3$totss, 
             km4$tot.withinss/km4$totss,km5$tot.withinss/km5$totss, km6$tot.withinss/km6$totss, 
             km7$tot.withinss/km7$totss,km8$tot.withinss/km8$totss)
nclus = c(1:8)

plot(nclus,ratio_ss)
#elbow is at k=6

#Visualise parameter changes (period, phase, amplitude) as PCA and colour by cluster
clock$cluster=km6$cluster
means=clock%>%group_by(cluster)%>%summarise(period=mean(period_diff),phase=mean(phase_diff),amp=mean(amp_diff),
                                            sd_period=sd(period_diff),sd_phase=sd(phase_diff),sd_amp=sd(amp_diff))
#Rejig for plotting
clock$row=1:59
clock=arrange(clock,cluster)
clock$cluster=c(rep(2,6),rep(4,11),rep(3,7),1,rep(5,14),rep(6,20))
clock$cluster=as.factor(clock$cluster)
clock=arrange(clock,row)

ggplot()+
  geom_point(pca_clock_diffs$x, mapping=aes(PC1,PC2,fill=as.factor(clock$cluster)),colour="transparent",size=2.75,pch=21)+
  scale_fill_manual("Cluster",values=c("#A881C5","#9FC783","#FFD45B","#FF8F57","#70BAFF","grey"))+
  labs(x="PC1 (54.20%)",y="PC2 (36.01%)")+
  theme_bw()+
  theme(panel.background = element_blank(),axis.line = element_line(size = 0.5, colour = "black"), 
        axis.ticks = element_line(colour="NA"), panel.border = element_blank(), panel.grid = element_blank(),
        axis.text.x = element_text(size = 8, colour = "black"),
        axis.text.y = element_text(size = 8, colour = "black"), axis.title = element_text(size = 10),
        legend.position = "none")

#Make summaries of what the clusters represent. I.e. plot phase, period and ampllitude changes per cluster

means$cluster=as.factor(c(2,4,3,1,5,6))
means$cluster=as.factor(means$cluster)

ggplot()+
  geom_bar(means, mapping=aes(y=period,x=cluster,fill=cluster),stat="identity",width=0.75)+
  geom_errorbar(means, mapping=aes(ymax=period+sd_period,ymin=period-sd_period,x=cluster),width=0.15)+
  scale_fill_manual("Cluster",values=c("#A881C5","#9FC783","#FFD45B","#FF8F57","#70BAFF","grey"))+
  labs(y="\nMean difference (h)",title="Period")+
  geom_hline(yintercept=0,size=0.25)+
  theme_bw()+
  theme(panel.background = element_blank(),axis.line = element_blank(), 
        axis.ticks = element_line(colour="NA"), panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 8, colour = "black"), axis.title = element_text(size = 10),
        axis.title.x = element_blank(), axis.line.x=element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5,size=12),panel.border=element_rect())

ggplot()+
  geom_bar(means, mapping=aes(y=phase,x=cluster,fill=cluster),stat="identity",width=0.75)+
  geom_errorbar(means, mapping=aes(ymax=phase+sd_phase,ymin=phase-sd_phase,x=cluster),width=0.15)+
  scale_fill_manual("Cluster",values=c("#A881C5","#9FC783","#FFD45B","#FF8F57","#70BAFF","grey"))+
  labs(y="Mean difference (h)",title="Phase")+
  geom_hline(yintercept=0,size=0.25)+
  theme_bw()+
  theme(panel.background = element_blank(),axis.line = element_blank(), 
        axis.ticks = element_line(colour="NA"), panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 8, colour = "black"), axis.title = element_text(size = 10),
        axis.title.x = element_blank(), axis.line.x=element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5,size=12),panel.border=element_rect())

ggplot()+
  geom_bar(means, mapping=aes(y=amp,x=cluster,fill=cluster),stat="identity",width=0.75)+
  geom_errorbar(means, mapping=aes(ymax=amp+sd_amp,ymin=amp-sd_amp,x=cluster),width=0.15)+
  scale_fill_manual("Cluster",values=c("#A881C5","#9FC783","#FFD45B","#FF8F57","#70BAFF","grey"))+
  labs(y="\nAmplitude change (log(Mature:senescence))", title="Amplitude")+
  geom_hline(yintercept=0,size=0.25)+
  theme_bw()+
  theme(panel.background = element_blank(),axis.line = element_blank(), 
        axis.ticks = element_line(colour="NA"), panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 8, colour = "black"), axis.title = element_text(size = 10),
        axis.title.x = element_blank(), axis.line.x=element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5,size=12),panel.border=element_rect())

#Plot representative gene expression from each transcript cluster

#Compile representatives
clock_shapes=c("TraesMAC2D03G01101810","TraesMAC3B03G01584430","TraesMAC3D03G02002020","TraesMAC6B03G03549430","TraesMAC7B03G04136500","TraesMAC4B03G02333170")

#Get gene expression data for these transcripts
df = data.frame(t(all[as.character(clock_shapes),]))
df$sample_id= metadata$sample_id
df = merge(df, metadata[, .(sample_id, cond, time)], by = 'sample_id')

df_=NULL
for (x in 2:(ncol(df)-2)){
  y=data.frame(df[,c(x,8,9)])
  colnames(y)[1]="V1"
  df_=rbind(df_,y)
}

df_$cond=rep(c(rep("Mature",48),rep("Senescent",48)),(ncol(df)-3))
df_$locus=rep(c("Ppd_2D","TaGI_3B","TaLUX/BOA_3D","TaTOC1_6B","TaLHY_7B","TaPRR73_4B"),each=96)#Naming system has been updated. These are old


df_means = df_ %>%
  group_by(time,cond,locus)%>%
  summarise(mean=mean(V1),
            sd=sd(V1))

df_means$locus=factor(df_means$locus,levels=c("Ppd_2D","TaGI_3B","TaLUX/BOA_3D","TaTOC1_6B","TaLHY_7B","TaPRR73_4B"))

#Scaling function for y axis
scaleFUN <- function(x) sprintf("%.0f", x)

ggplot(df_means, aes(y=mean,x=time,colour=cond))+
  facet_wrap(~locus, scales="free_y",nrow=1)+
  geom_point(size=1)+
  geom_line()+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), size = 0.5, width = 0)+
  scale_colour_manual("",labels =c("Mature","Senescent"),values=c("#197000","#BD9B00"))+
  labs(x = "\nTime (ZT)", y="Mean expression (FPKM)\n")+
  theme_bw()+
  theme(axis.ticks = element_line(colour="NA"), panel.border = element_rect(colour="light grey",size=1), 
        panel.grid.major = element_line(colour="light grey", size=0.1),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"), axis.title = element_text(size = 11),
        strip.background = element_rect(colour=NA,fill=NA), strip.text = element_text(size=11),
        legend.text = element_text(size=10),legend.position = "right")+
  guides(shape = guide_legend(override.aes = list(size = 3)))+
  scale_x_continuous(breaks=seq(48,84,12))+
  scale_y_continuous(labels=scaleFUN)


#Now make a list of each cluster's members
ggplot()+
  geom_point(data.frame(x=seq(1,6,1),y=rep(6.3,6)),mapping=aes(x=x,y=y),size=4,
             colour=c("#A881C5","#9FC783","#FFD45B","#FF8F57","#70BAFF","grey"))+
  geom_text(data.frame(x=seq(1,6,1),y=rep(7.5,6),text=paste("Cluster ",seq(1,6,1),sep="")),mapping=aes(x=x,y=y,label=text),fontface="bold")+
  geom_text(arrange(clock,cluster),mapping=aes(x=c(rep(1,1),rep(2,6),rep(3,7),rep(4,11),rep(5,14),rep(6,20)),y=c(5,rev(seq(0,5,1)),rev(seq(-1,5,1)),rev(seq(-5,5,1)),rev(seq(-8,5,1)),rev(seq(-14,5,1))),label=Name),size=4)+
  geom_vline(xintercept=seq(1.5,5.5,1),linewidth=0.25)+
  scale_x_continuous(limits=c(0.85,6.15))+
  theme(plot.background = element_blank(),axis.line = element_blank(),axis.title = element_blank(),
        axis.ticks = element_blank(),axis.text = element_blank(),panel.background =element_blank())

#Expression plots for all putative clock genes

clock_Mace = filter(clock_Mace, nchar(Name) > 3)
clock_Mace=clock_Mace[order(clock_Mace$Name),]
clock_Mace=clock_Mace[-c(1:12,46:48),]#filter to most important oscillator components
clock_Mace=clock_Mace[order(clock_Mace$Gene.stable.ID),]
df = data.frame(t(all[which(rownames(all) %in% unlist(clock_Mace[,2])),]))
df$sample_id= metadata$sample_id
df = merge(df, metadata[, .(sample_id, cond, time)], by = 'sample_id')


df_=NULL
for (x in 2:(ncol(df)-2)){
  y=data.frame(df[,c(x,66,67)])
  colnames(y)[1]="V1"
  df_=rbind(df_,y)
}

df_$cond=rep(c(rep("Mature",48),rep("Senescent",48)),(ncol(df)-3))
df_$locus=rep(clock_Mace$Name,each=96)


df_means = df_ %>%
  group_by(time,cond,locus)%>%
  summarise(mean=mean(V1),
            sd=sd(V1))

#Order by approximate phase
#Mat - getting mean phase
x=mat_bio[mat_bio$ID %in% clock_Mace$Gene.stable.ID,]
x=merge(x[,1:10],clock_Mace[,c(2,3)], by.x="ID", by.y="Gene.stable.ID")
m=x%>%group_by(Name)%>%summarise(mean_phase=mean(LAG))
#Sen - getting mean phase
z=mat_bio[mat_bio$ID %in% clock_Mace$Gene.stable.ID,]
z=merge(z[,1:10],clock_Mace[,c(2,3)], by.x="ID", by.y="Gene.stable.ID")
s=z%>%group_by(Name)%>%summarise(mean_phase=mean(LAG))

#Average for mat and sen across homoeologue groups
m_s=rbind(m,s)
m_s=m_s%>%group_by(Name=substr(Name, 1,nchar(m_s$Name)-1))%>%summarise(phase=mean(mean_phase))
order=m_s[order(m_s$phase),]$Name
order=paste(rep(order,each=3),c("A","B","D"),sep="")
order=order[order %in% clock_Mace$Name]

df_means$locus=factor(df_means$locus, levels = order)

ggplot(df_means, aes(y=mean,x=time,colour=cond))+
  facet_wrap(~locus, scales="free_y",nrow=11)+
  geom_point(size=0.5)+
  geom_line(size=0.25)+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), size = 0.25, width = 0)+
  scale_colour_manual("",labels =c("Mature","Senescent"),values=c("#197000","#BD9B00"))+
  labs(x = "\nTime (ZT)", y="Mean expression (FPKM)\n")+
  theme_bw()+
  theme(axis.ticks = element_line(colour="NA"), panel.border = element_rect(colour="light grey",size=0.5), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"), axis.title = element_text(size = 16),
        strip.background = element_rect(colour=NA,fill="white"), strip.text = element_text(size=11),
        legend.text = element_text(size=14), legend.title = element_text(size=16),legend.position = "top")+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_x_continuous(breaks=seq(48,84,12))

