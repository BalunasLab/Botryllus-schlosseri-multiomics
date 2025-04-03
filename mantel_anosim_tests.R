install.packages("vegan","geosphere","dplyr")
library(vegan)
library(geosphere)
library(dplyr)

data<- read.csv("Ametabolomics.csv",header = TRUE, row.names = 1)
head(data)
data_m<- as.matrix(data)

groups<- c("botryllus","botryllus","botryllus","botryllus","botryllus","botryllus","botryllus","botryllus", "seawater", "seawater","seawater")

#####Now for ANOSIM, first between seawater and botryllus...
#Adding the groups
data$groups<-groups
anosim_result<- anosim(data_m, data$groups, distance = "bray", permutations = 999)
anosim_result
#When interpreting these results you want to look at the ANOSIM statistic R and the Significance values. A Significance value less than 0.05 is generally considered to be statistically significant, and means the null hypothesis can be rejected.

#“The ANOSIM statistic “R” compares the mean of ranked dissimilarities between groups to the mean of ranked dissimilarities within groups. An R value close to “1.0” suggests dissimilarity between groups while an R value close to “0” suggests an even distribution of high and low ranks within and between groups”. In other words, the higher the R value, the more dissimilar your groups are in terms of microbial community composition.

#Data is significant

##PERMANOVA
##calculate the distance first
set.seed(123)
data_bray<- vegdist(data_m, method = "bray")
perm.met<-adonis2(data_bray~as.factor(data$groups), data=data, permutations = 999, method ="bray")
perm.met


####Now for individuals!!!#### 
#Can't do for individuals because there is no replication
#for botryllus we will use the panmetabolome data!
panmet<- read.csv("panmet_binary.csv", header = TRUE, row.names = 1)
panmet_abundance<- read.csv("panmet_abundance.csv",header = TRUE, row.names = 1)
indiv<- c("B1", "B2","B3","B4","B5","B6","B7","B8")
panmet_m <- as.matrix(panmet)
panmet_abundance_m<- as.matrix(panmet_abundance)

panmet$indiv<-indiv
panmet_abundance$indiv<-indiv

###First with braycurtis distance 
Pan_anosimbray<- anosim(panmet_m, panmet$indiv, distance = "bray", permutations = 999)
Pan_anosimbray

###Euclidean
Pan_anosimE<- anosim(panmet_m, panmet$indiv, distance = "euclidean", permutations = 999)
Pan_anosimbray



###Now for the for the mantel test we will use the community, metals and metabolomics###
library(vegan)
library(geosphere)
library(dplyr)

metabo<- read.csv("Bmetabolomics.csv", header = TRUE,row.names = 1)
metals<- read.csv("Bmetals.csv", header = TRUE, row.names = 1)
community<- read.csv("Ccommunity.csv",  row.names = 1)
sets<- c("botryllus","botryllus","botryllus","seawater","seawater")
##For botryllus

Bmetabo<- metabo[1:3,1:ncol(metabo)]
Bcom<- community[1:3,1:ncol(community)]
Bmet<- metals[1:3,1:ncol(metals)]

#Now let's calculate the distance(only for Botryllus samples)
dist.metabo<- vegdist(Bmetabo, method = "bray")
dist.com<- vegdist(Bcom, method = "bray")
dist.met <- vegdist(Bmet, method = "euclidean")##MAybe euclidean might be better? I'll try that one 
##I'm doing three pairwise mantel tests

#metabolomics vs community
metabo_comMantel<- mantel(dist.metabo,dist.com, method = "spearman", permutations = 999, na.rm = TRUE)
metabo_comMantel
#metabolomics vs metals
metabo_metalsMantel<- mantel(dist.metabo, dist.met, method = "spearman", permutations = 999, na.rm = TRUE
)
metabo_metalsMantel

#community vs metals
com_metMantel<- mantel(dist.com, dist.met, method = "spearman", permutations = 999, na.rm = TRUE)
com_metMantel

#Permanova 
metabo$sets<- sets
community$sets<-sets
metals$sets<-sets
###Now let's do a permanova for botryllus and seawater
dist.metabo_BS<- vegdist(metabo, method = "bray")
dist.com_BS<- vegdist(community, method = "bray")
dist.met_BS <- vegdist(metals, method = "euclidean")
set.seed(123)
perm.metabo<-adonis2(dist.metabo_BS~as.factor(metabo$sets), data=metabo, permutations = 999, method ="bray")
perm.metabo


perm.comm<- adonis2(dist.com_BS~as.factor(community$sets), data= community, permutations =999, method ="bray")
perm.comm
