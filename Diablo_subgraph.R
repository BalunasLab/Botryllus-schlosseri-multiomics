
###Diablo model with mixomics package version 6.23.2###
#I'll run two different models one that will be global, including all the metabolome, microbiome and metallome and one that will be decided by the model in the tunning option
# 04.14.23 First trial
install.packages(c("mixOmics", "IMIFA", "compositions","phyloseq","microbiome","igraph","dplyr","ggplot2"))

library(mixOmics)
library(IMIFA)
library(compositions)
library(phyloseq)
library(microbiome)
library(igraph)
library(dplyr)
library(ggplot2)

#Also let's try to be more organized. Data will be Dmetabolomics, Dcommunity and Emetals

#Reading the data

M<- read.csv("Dmetabolomics.csv", header = TRUE,row.names = 1)
Me<- read.csv("Emetals.csv", header = TRUE, row.names = 1)
Mi<- read.csv("Dcommunity.csv",  row.names = 1)

#Now transforming the data....

#For metals I'll use scale
SMe<- scale(Me, center= TRUE, scale = TRUE)
#For metabolomics I'll use square rootand pareto scale but first let's see distribution
sqM<- sqrt(M)
#pmetabo<- pareto_scale(sqmetabo, centering= TRUE)
Mdf <- as.data.frame(sqM)
#For community now with center log ratio (following what alyssa did) but from her code!
#community<- clr(community)
Mim <- as.matrix(Mi)
Miotu = otu_table(Mim, taxa_are_rows = FALSE)

# CLR transform
MiotuT <- microbiome::transform(Miotu, 'clr')
MiotuTdf<-as.data.frame(MiotuT) 


#remove where stdev = 0 ;for now only for community and I'll keep the transformed metabolomics data
#pmetabo2<-Filter(function(x) sd(x) != 0, pmetabo)
MituTdf <- Filter(function(x) sd(x) != 0, MiotuTdf)
Mdf <- Filter(function(x) sd(x) !=0, Mdf)
pMdf<- pareto_scale(Mdf, centering = 
                      TRUE)


###CONVERT Files for DIABLO
MM<-data.matrix(pMdf, rownames.force = TRUE)
MiM<-data.matrix(MituTdf, rownames.force = TRUE)
MeM <-SMe

##Now let's make the list
MMiMe <- list(metabolome=MM, microbiome=MiM, metallome=MeM)

types<- factor(c("tunnicate", "tunnicate", "tunnicate", "seawater","seawater", "seawater"))

#Let's do the design
designM2 = matrix(0.1, ncol = length(MMiMe), nrow = length(MMiMe), 
                  dimnames = list(names(MMiMe), names(MMiMe)))
diag(designM2) = 0 # set diagonal to 0s
designM2

#Now the tunning
tune.MMiMe <- tune.block.splsda(X= MMiMe, Y= types , ncomp = 3, validation = "Mfold" , folds = 3, nrepeat = 10, dist = "centroids.dist", near.zero.var = TRUE, design = designM2)


#this is to extract which features I will be using...
list.keepX<- tune.MMiMe$choice.keepX

#Now let's run the tuned model

result.diablo.MMiMe.tuned <- block.splsda(MMiMe,types,ncom=3, keepX = list.keepX, design = designM2, near.zero.var = TRUE)

plotIndiv(result.diablo.MMiMe.tuned,ind.names = FALSE, style = 'graphics')

plotDiablo(result.diablo.MMiMe.tuned, ncomp = 1)

plotDiablo(result.diablo.MMiMe.tuned, ncomp = 2)

plotDiablo(result.diablo.MMiMe.tuned, ncomp = 3)

xaxis<-list(c(-2.5,2.5),c(-2.5,2.5),c(-2,2))
yaxis<-list(c(-4,3),c(-15,10),c(-0.7,1.0))

##now with the axis and everything else for the block graphs
plotIndiv(result.diablo.MMiMe.tuned, style = 'graphics', study = "global", col.per.group = c('#3ec6ef','#f35b04'), subtitle = c("Metabolome","Microbiome","Metallome"), xlim = xaxis, ylim = yaxis)

plotDiablo(result.diablo.MMiMe.tuned, ncomp = 1, col.per.group = c('#3ec6ef','#f35b04'))

plotDiablo(result.diablo.MMiMe.tuned, ncomp = 2, col.per.group = c('#3ec6ef','#f35b04'))

plotDiablo(result.diablo.MMiMe.tuned, ncomp = 3,col.per.group = c('#3ec6ef','#f35b04'))



##Let's do the circos with cutoff of 0.7

uno<-1
dos<-2
tres<-3
##Cutoff 0.95
circosPlot(result.diablo.MMiMe.tuned, comp = 1:2, cutoff = 0.9, line = TRUE, color.Y = c('#3ec6ef','#f35b04'),linkWidth = c(0.2:1.0),
           color.blocks= c('#8d6a9f','#d3fad6', '#d7907b'), 
           color.cor = c("navyblue", "brown1"), size.labels = 1.5,size.variables = 0.6, size.legend = 1.1)

circosPlot(result.diablo.MMiMe.tuned, comp = 1:2, cutoff = 0.9, line = TRUE, color.Y = c('#3ec6ef','#f35b04'),linkWidth = c(0.2:1.0),
           color.blocks= c('#8d6a9f','#d3fad6', '#d7907b'), 
           color.cor = c("#394ba0", "#cd0000"), size.labels = 1.5,size.variables = 0.6, size.legend = 1.1)


circosPlot(result.diablo.MMiMe.tuned, cutoff = 0.9, line = TRUE, color.Y = c('#00cffc','#fbaf00'),
           color.blocks= c('#459268', '#8fb339', '#ddfff7'), 
           color.cor = c("navyblue", "brown1"), size.labels = 1.5,size.variables = 0.6, size.legend = 1.1)

###Now let's do the networks
##all together
par(mar=c(1, 1, 1, 1))
my.networkmime_tunned<- network(result.diablo.MMiMe.tuned,blocks = c(1,2,3), color.node = c('#459268', '#8fb339', '#ddfff7'), cutoff = 0.7)



write.graph(my.networkmime_tunned$gR, file = "myNetwork07MiMeT.gml", format = "gml")

# circular plot with correlations for tuned model
#component 1 and 2
plotVar(result.diablo.MMiMe.tuned, style= 'lattice', legend = TRUE, legend.title = "")
#no overlapping
plotVar(result.diablo.MMiMe.tuned, overlap = FALSE, style= 'lattice', legend = TRUE, legend.title = "")
#component 1 and 3
plotVar(result.diablo.MMiMe.tuned,comp = c(1,3), style= 'lattice', legend = TRUE, legend.title = "")
#no overlapping
plotVar(result.diablo.MMiMe.tuned,comp = c(1,3), overlap = FALSE, style= 'lattice', legend = TRUE, legend.title = "")

##Now let's run the full model

result.diablo.MMiMe<- block.splsda(MMiMe, types, ncomp = 3, design = designM2, near.zero.var = TRUE)

plotIndiv(result.diablo.MMiMe,ind.names = FALSE, style = 'graphics')

xaxis<-list(c(-30,30),c(-40,60),c(-4,5))
yaxis<-list(c(-15,10),c(-25,30),c(-1.0,0.8))

plotIndiv(result.diablo.MMiMe, style = 'graphics', study = "global", col.per.group = c('#00cffc','#fbaf00'), subtitle = c("Metabolome","Microbiome","Metallome"), xlim = xaxis, ylim = yaxis)


plotDiablo(result.diablo.MMiMe, ncomp = 1, col.per.group = c('#00cffc','#fbaf00'))

plotDiablo(result.diablo.MMiMe, ncomp = 2, col.per.group = c('#00cffc','#fbaf00'))

plotDiablo(result.diablo.MMiMe, ncomp = 3,col.per.group = c('#00cffc','#fbaf00'))

uno<-1
dos<-2
tres<-3
#circos plot cutoff 0.9
circosPlot(result.diablo.MMiMe, cutoff = 0.9,comp = uno, line = TRUE, color.Y = c('#00cffc','#fbaf00'),
           color.blocks= c('#459268', '#8fb339', '#ddfff7'), 
           color.cor = c("navyblue", "brown1"), size.labels = 1.5,size.variables = 0.6, size.legend = 1.0)

###Now let's do the networks
##all together
par(mar=c(1, 1, 1, 1))
my.networkmime<- network(result.diablo.MMiMe,blocks = c(1,2,3), color.node = c('#459268', '#8fb339', '#ddfff7'), cutoff = 0.8)


write.graph(my.networkmime$gR, file = "test2_myNetwork08MiMe.gml", format = "gml")

write.table()
#component 1 and 2
plotVar(result.diablo.MMiMe, style= 'lattice', legend = TRUE, legend.title = "")
#no overlapping
plotVar(result.diablo.MMiMe, overlap = FALSE, style= 'lattice', legend = TRUE, legend.title = "")
#component 1 and 3
plotVar(result.diablo.MMiMe,comp = c(1,3), style= 'lattice', legend = TRUE, legend.title = "")
#no overlapping
plotVar(result.diablo.MMiMe,comp = c(1,3), overlap = FALSE, style= 'lattice', legend = TRUE, legend.title = "")


#Now to delimit the interactions and find a good threshold for the full network, I'll use associationsubgraphs
install.packages(c("associationsubgraphs","tidyr","stringr","dplyr"))
library(associationsubgraphs)
library(tidyr)
library(stringr)
library(dbplyr)

network_bot<- read.csv("myNetwork08MiMe._impOTU_Bot.csv", header = TRUE)

network_bot[c("a","b")] <- str_split_fixed(network_bot$shared.name, ".(interacts with).", n=2)
head(network_bot)

network_bot2<- select(network_bot, c("weight", "a", "b"))
head(network_bot2)
network_bot2<- network_bot2 [c("a","b","weight")]
head(network_bot2)

#Let's do positive first
network_Bpos<- filter(network_bot2, weight > 0)
head(network_Bpos)
network_Bpos<-as.data.frame(network_Bpos)
network_Bpos<- rename(network_Bpos,strength=weight)

visualize_subgraph_structure(network_Bpos)
#I think I'll settle for 0.980? or 0.981
networkBposf<- filter(network_Bpos, strength>=0.94)
write.csv(networkBposf, "network_Bpf_094.csv")

#I'll use 0.96 for both!


##Now for negative

network_Bneg<- filter(network_bot2, weight<0)
head(network_Bneg)

network_Bneg<- rename(network_Bneg,strength=weight)
network_Bneg<- mutate(network_Bneg, strength= abs(strength))
head(network_Bneg)

visualize_subgraph_structure(network_Bneg)
networkBnegf<- filter(network_Bneg, strength>=0.94)
networkBnegf<-mutate(networkBnegf, strength = -1*(strength))
write.csv(networkBnegf, "network_Bnf_094.csv")
head(networkBnegf)

#For 0.98
networkBnegf_096<- filter(network_Bneg, strength>=0.96)
networkBnegf_096<-mutate(networkBnegf_096, strength = -1*(strength))
write.csv(networkBnegf_096, "network_Bnf_096.csv")
