##Now I'll try for upset plots but with stacked bars##
install.packages(c("UpSetR","mltools","data.table","VennDiagram","grid","futile.logger","dplyr"))
library(UpSetR)
library(mltools)
library(data.table)
library(VennDiagram)
library(grid)
library(futile.logger)
library(dplyr)
#first to read the table that will have the info I need

origins<- read.csv("features_annot_origin.csv")
#Now let's get the info into the botryllus table first

bot_data <- read.csv("botryllus_t5_updated.csv")
bot_data <- rename(bot_data,compound = Compound)
bot_data_origin<- left_join(bot_data, origins, by="compound")
write.csv(bot_data_origin, file = "botryllus_data_withOrigins.csv")
pdf(file="~/Documents/UMICH/Research_projects/botryllus/R_data/upset&venn/upset_wstackedbar_updated.pdf",
    width = 10,
    height = 8)
upset(bot_data_origin, nsets = 8, mb.ratio = c(0.65, 0.35), order.by = "freq", point.size = 3.5, line.size = 1.5, mainbar.y.label = "Shared features", sets.x.label = "Number of features", text.scale = c(1.8, 1.4, 1.3, 1.3, 1.5, 1.4), queries = list (
  list (query= elements, params = list("origin2", c("mixed", "microbial", "not annotated", "invertebrate")), color = "#ef476f" , active = T),# this is for microbial attributes
  list(query= elements, params = list("origin2", c("mixed", "not annotated", "invertebrate")), color = "#ffd166", active =T), #this is for invertebrates attributes
  list(query=elements, params = list ("origin2", "not annotated", "mixed"), color = "#06d6a0", active = T),#this is for mixed attributes
  list(query=elements, params = list ("origin2", "not annotated"), color = "#264653", active = T)), query.legend = "top") #this is for not annotated attributes
dev.off()

write.csv(bot_data_origin, file = "botpanmet.csv")

#Now for an upset plot from seawater vs botryllus and annotation
botsea<- read.csv("binary_botsea.csv")

botsea_origin <- left_join(botsea, origins, by= "compound")
## Will write the output of the joints as a csv 
write.csv(botsea_origin, file = "botsea&originupdated.csv")
test<- read.csv("botsea&origin.csv")

#Now the upset plot
pdf(file="~/Documents/UMICH/Research_projects/botryllus/R_data/upset&venn/upset_botsea3.pdf",
    width = 10,
    height = 8)
upset(botsea_origin, nsets=2, order.by = "freq", mb.ratio = c(0.65, 0.35), point.size = 5.5, line.size = 2, mainbar.y.label = "Shared features", sets.x.label = "Number of features", text.scale = c(1.8, 1.6, 1.5, 1.5, 1.7, 1.6), queries = list (
  list(query = elements, params =list ("origin.3", c("annotated", "not annotated")), color = "#c39ed5", active = T), #this is for annotated
  list(query =elements, params= list ("origin.3", "not annotated"), color = "#264653", active =T)), query.legend = "top") #this is for unnanotated
dev.off()