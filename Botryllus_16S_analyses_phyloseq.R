library("qiime2R", lib.loc="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")
library("ggplot2", lib.loc="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")
library("RColorBrewer", lib.loc="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")
library("vegan", lib.loc="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")
library("dplyr")
library("scales", lib.loc="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")

physeq_Bs_only<-qza_to_phyloseq(
  features="table-sans-Blanks_R1.qza",
  tree="rooted-tree.qza",
  "taxonomy_Bs.qza",
  metadata = "Bs_metadata.txt"
)

physeq_Bs<-qza_to_phyloseq(
  features="merged_table_Bs.qza",
  tree="rooted_tree.qza",
  "merged_taxonomy_Bs.qza",
  metadata = "Bs_metadata.txt"
)

physeq_BsR2<-qza_to_phyloseq(
  features="table-sans-Blanks_R1.qza",
  tree="rooted-tree.qza",
  "taxonomy_Bs.qza",
  metadata = "Bs_metadata_R2.txt"
)


physeq_BsR2

physeq_Bs2=prune_taxa(taxa_sums(physeq_BsR2)>0,physeq_BsR2)
physeq_Bs2
physeq_norm = transform_sample_counts(physeq_Bs2, function(x) 6933 * x/sum(x))

wunifrac_norm_all = phyloseq::distance(physeq_norm, method="unifrac", weighted=T)
ordination_all_norm = ordinate(physeq_norm, method="PCoA", distance=wunifrac_norm_all)

#TO get color palettes
# https://medialab.github.io/iwanthue/

coul <- brewer.pal(8, "Dark2") 
coul <- colorRampPalette(coul)(20)

custom<-c("#b92300",
           "#0073f5",
           "#feba49",
           "#0178c7",
           "#718a00",
           "#e1007e",
           "#00d27c",
           "#dc0055",
           "#00c798",
           "#94226d",
           "#bfce79",
           "#ff73b1",
           "#385d00",
           "#9a2547",
           "#009882",
           "#784704",
           "#00663e",
           "#fbb699",
           "#719161")

custom2<-c("#016dbb",
           "#58b946",
           "#792295",
           "#007f2e",
           "#7763da",
           "#d2db83",
           "#c40066",
           "#e84ef7",
           "#6b0080",
           "#1d4c00",
           "#47a5ff",
           "#f07139",
           "#f5ff00",
           "#825900",
           "#8d000c",
           "#793c19","#dc0055", "pink","grey")

p_nmds_norm1=plot_ordination(physeq_norm, ordination_all_norm, color="Sample")+
     geom_point(size=5)+
     geom_point(aes(fill=Sample), colour="black", size=5,shape=21)+ 
     stat_ellipse(geom = "polygon",aes(fill = Sample), alpha = 0.15)+
     theme(legend.title = element_text(size=20),legend.text = element_text(size=18),axis.text = element_text(size=12),legend.key.size = unit(2,"line"))+theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13))

p_nmds_norm1

ps2.rarefied = rarefy_even_depth(physeq_Bs, rngseed=1, sample.size=min(sample_sums(physeq_Bs)), replace=F)
wunifrac_r_all = phyloseq::distance(ps2.rarefied, method="unifrac", weighted=T)
ordination_all_r = ordinate(ps2.rarefied, method="NMDS", distance=wunifrac_r_all)

p_nmds_rare1=plot_ordination(ps2.rarefied, ordination_all_r, color="Sample")+
  geom_point(size=5)+
  geom_point(aes(fill=Sample), colour="black", size=5,shape=21)+ 
  stat_ellipse(geom = "polygon",aes(fill = Sample), alpha = 0.15)+
  theme(legend.title = element_text(size=20),legend.text = element_text(size=18),axis.text = element_text(size=12),legend.key.size = unit(2,"line"))+theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13))+
  scale_fill_manual(values=c("#f9af1a","#3ec6ef"))

p_nmds_rare1
library("ggpubr")
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "ns", "ns"))
a_my_comparisons <- list(c("Bschlosseri","seawater"))

p1=plot_richness(ps2.rarefied, x="Sample", measures=c("Observed","Shannon","InvSimpson"), color = "Sample")+
    geom_boxplot(aes(fill=Sample))+
     theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12),strip.text.x = element_text(size=13),axis.text.y = element_text(size=12),axis.title = element_text(size=13))+
     scale_fill_manual(values=c("#f9af1a","#3ec6ef"))+scale_color_manual(values=color3)+
  stat_compare_means(comparisons = a_my_comparisons, symnum.args = symnum.args, method = "t.test")

p1


speciesList_ID<-tapply(sample_names(physeq_norm), get_variable(physeq_norm, "Description"), c)

speciesPhyseq <- lapply(speciesList_ID, prune_samples, physeq_norm)

speciesOTUtable <- lapply(speciesPhyseq,otu_table)
speciesAvg <- lapply(speciesOTUtable,rowMeans)
pooledOTUtable = t(do.call(rbind,speciesAvg))
pooledOTUtable = data.frame(OTU=row.names(pooledOTUtable),pooledOTUtable)
View(pooledOTUtable)

TT = tax_table(physeq_norm)
TT = TT[, which(apply(!apply(TT, 2, is.na), 2, any))]

tdf = data.frame(TT, OTU = taxa_names(physeq_norm))
pOTUtax = merge(pooledOTUtable, tdf, by.x = "OTU")
write.table(pOTUtax,"pOTUtax_R2.txt",sep="\t")
write.table(pOTUtax,"pOTUtax_merged.txt",sep="\t")

pOTU = data.frame(pOTUtax,SeqTotal = rowSums(pOTUtax[,2:12]))
pOTU_class=read.csv("pOTU_Class.csv",header=TRUE)
pOTU.avg.c =pOTU_class[,c(2:12,15)]
melt.class = melt(pOTU.avg.c,id.vars="Class")
colnames(melt.class)[2]="species"

agg.class=aggregate(.~Class+species,melt.class,sum)
write.csv(agg.class,"agg_class.csv")
agg.class=read.csv("agg.class2.csv",header=TRUE)

bar2=ggplot(agg.class,aes(x=species,y=value,fill=Class)) +
  geom_bar(stat="identity",position="fill",color="black") +
  scale_y_continuous(labels = percent_format())+
  xlab("Species") +
  ylab("Relative Abudance")+
  theme(axis.text.y  = element_text(size=12))+
  theme(axis.text.x = element_text(angle=30,vjust=1, hjust = 1, size=12),legend.key.size =  unit(1.5,"line"))+
  facet_grid(.~sample,drop = TRUE,space = "free",scales="free")+
  scale_fill_manual(values=custom2)
bar2
#Core
core=read.table("core-features-1.tsv",sep="/t")
View(core)

pOTU_core=read.csv("pOTU_core.csv",header=TRUE)
pOTU.avg.c =pOTU_core[,c(2:12,15)]
pOTU.avg.c =core_potu[,c(2:12,15)]

#To make box plot for Taxas ==
box<-read.table("per_rel_otu_box.txt",header=TRUE)
melt_box<-melt(box)
b1<-ggplot(melt_box,aes(x=OTU,y=value,color=OTU))+geom_boxplot()+
  geom_jitter(position=position_jitter(0.2))+theme_bw()+
theme(axis.text.x = element_text(angle=30,hjust=0.95))+
xlab("Taxa")+ylab("% Relative Abundance")

#To add title facet ==
melt_box$title<-"Botryllus"
melt_box_sw$title<-"Seawater"

# To run Deseq ==
diagdds = phyloseq_to_deseq2(ps2.rarefied, ~ Sample)
res2 <- results(diagdds, contrast=c("Sample","Bschlosseri","seawater"))
alpha = 0.05
sigtab = res2[which(res2$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps2.rarefied)[rownames(sigtab), ], "matrix"))
dim(sigtab)
top<-sigtab[(sigtab$baseMean>50) & (abs(sigtab$log2FoldChange)>4),]
top2<-
dsq_p1<-ggplot(top, aes(x=Class, y=log2FoldChange, shape=Phylum,color=Order)) + geom_point(size=6) + 
   theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+scale_colour_manual(values = coul)
dsq_p1

dsq_p2<-ggplot(top, aes(x=Class, y=log2FoldChange, shape=Phylum,color=Order)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+scale_colour_manual(values = coul)
dsq_p2

sigtab_l2f = cbind(as(sigtab_l2f, "data.frame"), as(tax_table(ps2.rarefied)[rownames(sigtab_l2f), ], "matrix"))

sigtab2 = res2[which(res2$pvalue < alpha), ]
View(sigtab2)
sigtab2 = cbind(as(sigtab2, "data.frame"), as(tax_table(ps2.rarefied)[rownames(sigtab2), ], "matrix"))
top2<-sigtab2[(sigtab2$baseMean>5) & (abs(sigtab$log2FoldChange)>2),]
write.table(top2,"DE_top2.txt",sep="\t") #I manually blasted unannotate Gammas
top2<-read.table("DE_top2.txt",header=TRUE,row.names = "X")
top2_ed<-read.table("DE_top2_ed.txt",header=TRUE,row.names = "X")

#dsq_p3<-ggplot(top2_ed, aes(x=Class, y=log2FoldChange, shape=Phylum,color=Order)) + geom_point(size=6) + 
#  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+scale_colour_manual(values = custom2)+scale_shape_manual(values=c(15,16,17,18,9,7))

dsq_p3<-ggplot(top2_ed, aes(x=Class, y=log2FoldChange,color=Order))  + 
      theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+scale_fill_manual(values = custom2)+geom_point(aes(fill=Order), colour="black", size=5,shape=21)+
  theme(axis.text = element_text(size=12),legend.text = element_text(size=12,face="italic"),legend.title = element_text(size=14),axis.title.y = element_text(size=12))

dsq_p3
fig2=ggarrange(p1, p_nmds_norm1,legend="bottom")
fig2

tiff("Fig1.tiff",width=10,height=5,units ="in",type=c("cairo"),res=300)
plot(fig2)
dev.off()

tiff("Fig1A.tiff",width=8,height=5,units ="in",type=c("cairo"),res=300)
plot(p1)
dev.off()

tiff("Fig1B.tiff",width=8,height=5,units ="in",type=c("cairo"),res=300)
plot(p_nmds_norm1)
dev.off()
#Merging deseq and otu tabel to make heatmap
top2_ed<-read.table("DE_top2_ed.txt",header=TRUE,row.names = "X")

library(tibble)
top2_ed<-rownames_to_column(top2_ed,"X")

otu_top2<-merge(top2_ed,pooledOTUtable,by.x="X",by.y="OTU")

otu_top2<-column_to_rownames(otu_top2,"X")
filtered_top <- subset(otu_top2, padj < 0.05)
View(filtered_top)
# Top 20 values (sort column 'A' in descending order)
filtered_top <- filtered_top[order(-filtered_top$log2FoldChange), ][1:25, ]

# Bottom 20 values (sort column 'A' in ascending order)
filtered_bottom <- filtered_top[order(filtered_top$log2FoldChange), ][1:25, ]

# Combine both top and bottom 20 into one data frame
combined_df <- rbind(filtered_top, filtered_bottom)
combined_df$Class <- make.unique(as.character(combined_df$Class))
rownames(combined_df) <- NULL
combined_df<-column_to_rownames(combined_df,"Class")
combined_df<-combined_df[-c(1:12)]

metadata = read.delim("Bs_metadata.txt")
meta_df<-column_to_rownames(metadata,"Description")
meta_df<-meta_df[-c(1:2)]
View(metadata)
anno_color<- list(
  Sample = c(
    "Bschlosseri_CT" = "#f9af1a",
    "seawater"="#3ec6ef"))
pheatmap(combined_df,show_rownames = TRUE, cluster_rows = TRUE,cluster_cols = TRUE, scale = 'row',annotation_colors = anno_color,annotation_col =meta_df,colorRampPalette(c(
  "#6a1185", "white", "#d4061a"
))(25
))


#Run t-test of class between seawater and Bs ===
# I added the relative abundance of all classes show in Fig. 2C and made a data frame 
per_rel<-read.delim("per_rel_otu_all.txt",row.names=1)
View(per_rel)

# Subset the data frame for each Class

alpha_values <- per_rel["b_Alpha", ]
View(alpha_values)
AVP_samples <- as.numeric(alpha_values[grep("AVP", colnames(alpha_values))])
seawater_samples <- as.numeric(alpha_values[grep("seawater", colnames(alpha_values))])
View(AVP_samples)
t.test(AVP_samples, seawater_samples, alternative="two.sided", var.equal=FALSE)
# Welch Two Sample t-test
# 
# data:  AVP_samples and seawater_samples
# t = -1.5935, df = 8.9938, p-value = 0.1664
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -22.989521  4.946197
# sample estimates:
#   mean of x mean of y 
# 19.33687  38.30473 

gamma_values <- per_rel["c_Gamma", ]
View(gamma_values)
AVP_samples <- as.numeric(gamma_values[grep("AVP", colnames(gamma_values))])
seawater_samples <- as.numeric(gamma_values[grep("seawater", colnames(gamma_values))])
View(AVP_samples)
t.test(AVP_samples, seawater_samples, alternative="two.sided", var.equal=FALSE)
# Welch Two Sample t-test
# 
# data:  AVP_samples and seawater_samples
# t = 4.9877, df = 8.4978, p-value = 0.0007227
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.3008674 13.3860887
# sample estimates:
#   mean of x mean of y 
# 16.531444  9.687966 

ep_values <- per_rel["d_Epsilonproteobacteria", ]
View(ep_values)
AVP_samples <- as.numeric(ep_values[grep("AVP", colnames(ep_values))])
seawater_samples <- as.numeric(ep_values[grep("seawater", colnames(ep_values))])
View(AVP_samples)
t.test(AVP_samples, seawater_samples, alternative="two.sided", var.equal=FALSE)

# Welch Two Sample t-test
# 
# data:  AVP_samples and seawater_samples
# t = 2.0806, df = 7.0176, p-value = 0.05192
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.6087351  9.5658965
# sample estimates:
#   mean of x mean of y 
# 4.8968700 0.4182893

fl_values <- per_rel["e_Flavobacteriia", ]
View(fl_values)
AVP_samples <- as.numeric(fl_values[grep("AVP", colnames(fl_values))])
seawater_samples <- as.numeric(fl_values[grep("seawater", colnames(fl_values))])
View(AVP_samples)
t.test(AVP_samples, seawater_samples, alternative="two.sided", var.equal=FALSE)
# Welch Two Sample t-test
# 
# data:  AVP_samples and seawater_samples
# t = -17.951, df = 8.8655, p-value = 2.827e-08
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -32.89475 -25.51671
# sample estimates:
#   mean of x mean of y 
# 6.728689 35.934420 

ver_values <- per_rel["f_Verruco", ]
View(ver_values)
AVP_samples <- as.numeric(ver_values[grep("AVP", colnames(ver_values))])
seawater_samples <- as.numeric(ver_values[grep("seawater", colnames(ver_values))])
View(AVP_samples)
t.test(AVP_samples, seawater_samples, alternative="two.sided", var.equal=FALSE)


# Welch Two Sample t-test
# 
# data:  AVP_samples and seawater_samples
# t = 1.754, df = 8.7708, p-value = 0.1142
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.2386857  1.8577539
# sample estimates:
#   mean of x mean of y 
# 2.497115  1.687581

#PERMDISP==

grouping_variable <- sample_data(physeq_norm)$Sample
# Run betadisper
bd <- betadisper(wunifrac_norm_all, grouping_variable)
summary(bd)
# Perform a significance test for differences in dispersion
anova(bd)
# Analysis of Variance Table
# 
# Response: Distances
# Df   Sum Sq   Mean Sq F value Pr(>F)  
# Groups     1 0.014873 0.0148734  5.5902 0.0423 *
#   Residuals  9 0.023946 0.0026606 
# Permutation test for significance
permutest(bd)

# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
# Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)  
# Groups     1 0.014873 0.0148734 5.5902    999  0.048 *
#   Residuals  9 0.023946 0.0026606                       
# Visualize the results
plot(bd)

