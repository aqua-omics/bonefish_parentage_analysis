#Scripts below have been adapted from Miller-Crews et al., 2021 and the R package dendextend's vignettes
#https://github.com/imillercrews/ParentageAnalysis
#https://cran.r-project.org/web/packages/dendextend/vignettes/

#Author: David J. Bradshaw II, Ph.D.
#Email: dbradshaw3366@gmail.com

#Install packages
#install.packages(c("tidyverse", "ggfortify", "pvclust", "gplots", "igraph"))


#Load libraries
library(tidyverse)
library(ggfortify)
library(pvclust)
library(gplots)
library(igraph)
library(dendextend)

##Load in metadata
data.bam.names.plinked <- read.csv('C:/Users/dbrad/Documents/HBOI/Bonefish_Revised_Manuscript/Bonefish_PA/bonefish_metadata_cd_plinked.csv')

#Create a sample id column
data.bam.names.plinked$sample_id <-
  gsub(".fq.bt2.bam", "", data.bam.names.plinked$bam_id)



##Load in Full IBS Matrix
adults.spawns.ibs <- read.delim("C:/Users/dbrad/Documents/HBOI/Bonefish_Revised_Manuscript/Bonefish_PA/WL_HBOI_R_1_hwe_0.1/plink.mibs", header = FALSE)


#Use metadata file to name rows and columns of the IBS dataframe
colnames(adults.spawns.ibs)=data.bam.names.plinked$pit_tag
rownames(adults.spawns.ibs)=data.bam.names.plinked$pit_tag

#Store order numbers for each type of sample
Spawn_1_Positions <- filter(data.bam.names.plinked, sex=="S1")$order
Spawn_2_Positions <- filter(data.bam.names.plinked, sex=="S2")$order
Spawn_Positions <- sort(c(Spawn_1_Positions, Spawn_2_Positions))
Female_Positions <- filter(data.bam.names.plinked, sex=="F")$order
Male_Positions <- filter(data.bam.names.plinked, sex=="M")$order

filter(data.bam.names.plinked, sex=="F")$order

####Females Prep####

#Filter Male samples from metadata
females.spawns.data.bam.names.plinked <- filter(data.bam.names.plinked, sex!="M")

#Subset main mibs for just females and spawns
females.spawns.ibs <- adults.spawns.ibs[c(1:8, Female_Positions), c(1:8, Female_Positions)]

####Spawn 1 + Females####

#Create Spawn 1 and Female focused dendrograms (all replicates, and separated replicates)
m.ibs.spawn1a.females.plinked.V2 <- females.spawns.ibs[c(1,9:31),c(1,9:31)]
m.ibs.spawn1b.females.plinked.V2 <- females.spawns.ibs[c(3,9:31),c(3,9:31)]
m.ibs.spawn1c.females.plinked.V2 <- females.spawns.ibs[c(4,9:31),c(4,9:31)]
m.ibs.spawn1d.females.plinked.V2 <- females.spawns.ibs[c(5,9:31),c(5,9:31)]
m.ibs.spawn1.females.plinked.V2 <- females.spawns.ibs[c(Spawn_1_Positions,9:31), 
                                                      c(Spawn_1_Positions,9:31)]

#Run each replicate separately and all replicates together to get idea of dendrogram
plot(pvclust(m.ibs.spawn1a.females.plinked.V2, method.hclust = "average",method.dist = 'euclidian'), main='Spawn 1a and Females')
plot(pvclust(m.ibs.spawn1b.females.plinked.V2, method.hclust = "average",method.dist = 'euclidian'), main='Spawn 1b and Females')
plot(pvclust(m.ibs.spawn1c.females.plinked.V2, method.hclust = "average",method.dist = 'euclidian'), main='Spawn 1c and Females')
plot(pvclust(m.ibs.spawn1d.females.plinked.V2, method.hclust = "average",method.dist = 'euclidian'), main='Spawn 1d and Females')
plot(pvclust(m.ibs.spawn1.females.plinked.V2, method.hclust = "average",method.dist = 'euclidian'), main='Spawn 1 and Females')

#Generally, a-d separately range from 0.42-0.24 in Height
#a-d together ranges from 0.5-0.05 in Height

#Prepare for plotting by making curated versions of metadata that only has one "offspring" sample each
sample.data.spawn1a.females.plinked <- females.spawns.data.bam.names.plinked[c(1,9:31),]
sample.data.spawn1b.females.plinked <- females.spawns.data.bam.names.plinked[c(3,9:31),]
sample.data.spawn1c.females.plinked <- females.spawns.data.bam.names.plinked[c(4,9:31),]
sample.data.spawn1d.females.plinked <- females.spawns.data.bam.names.plinked[c(5,9:31),]
sample.data.spawn1.females.plinked <- females.spawns.data.bam.names.plinked[c(Spawn_1_Positions,9:31),]

##Plot for Spawn 1a Females
#Save pvclust result, run with more permutations (nboot = 10000 instead of = 1000) than test above, as a variable and convert to a dendextend object
result <- pvclust(m.ibs.spawn1a.females.plinked.V2, method.hclust = "average",method.dist = 'euclidian', nboot = 10000)
dend <- as.dendrogram(result)

#By default, the dend has no colors to the labels so let's add some color:
#Save desired metadata category as a list of factors
status <- as.factor(sample.data.spawn1a.females.plinked$status)
status
#Save the desired metadata as a list of numbers (Alphabetical by default) associated with each metadata subtype
colors_to_use <- as.numeric(as.factor(sample.data.spawn1a.females.plinked[,11]))
colors_to_use
#Sort them based on their order in dend:
colors_to_use <- colors_to_use[order.dendrogram(dend)]
colors_to_use
#Convert the numbers to colors of your choice
colors_to_use[colors_to_use==1] <- "black"
colors_to_use[colors_to_use==2] <- "purple"
colors_to_use[colors_to_use==3] <- "blue"
colors_to_use[colors_to_use==4] <- "orange"
colors_to_use
#Save this list as the labels_colors for your dendextend object
labels_colors(dend) <- colors_to_use
# Now each state has a color
labels_colors(dend) 

#Play around with the ylim to get a view you like and add in the legend based upon your colors above and the results from pvclust
dend %>% hang.dendrogram()%>%  plot(ylim=c(0.20,0.44))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "purple", "blue", "orange"))

#Add a red rectangle around branchings with au pvalues > 80
pvrect(result, alpha=.80, pv="au", type="geq")

#Make pvclust statistics text smaller
result %>% text(cex=.75)

#Save the perfected image above
tiff(file="Spawn 1a with Females default.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram()%>%  plot(ylim=c(0.20,0.44))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "purple", "blue", "orange"))
result %>% text
dev.off()

#Image without statistics, but boxed around branchings with au p-values >80
tiff(file="Spawn 1a with Females boxed.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram()%>%  plot(ylim=c(0.20,0.44))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "purple", "blue", "orange"))
pvrect(result, alpha=.80, pv="au", type="geq")
dev.off()

#Image without statistics, but with bolded significant branchings with au p-values >80
tiff(file="Spawn 1a with Females bolded.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram() %>%
  pvclust_show_signif(result, signif_type =c("au"), alpha = 0.2)%>% 
  plot(ylim=c(0.20,0.44))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "purple", "blue", "orange"))
dev.off()

#Image without statistics, but with bolded significant branchings with au p-values >80, without a legend
tiff(file="Spawn 1a with Females bolded no legend.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram() %>%
  pvclust_show_signif(result, signif_type =c("au"), alpha = 0.2)%>% 
  plot(ylim=c(0.20,0.44))
dev.off()

##Plot for Spawn 1b Females
##Coloring labels based upon metadata column
#Save pvclust result as a variable and convert to a dendextend object
result <- pvclust(m.ibs.spawn1b.females.plinked.V2, method.hclust = "average",method.dist = 'euclidian', nboot = 10000)
dend <- as.dendrogram(result)

#By default, the dend has no colors to the labels so let's add some color:
#Save desired metadata category as a list of factors
status <- as.factor(sample.data.spawn1b.females.plinked$status)
#Save the desired metadata as a list of numbers (Alphabetical by default) associated with each metadata subtype
colors_to_use <- as.numeric(as.factor(sample.data.spawn1b.females.plinked[,11]))
colors_to_use
#Sort them based on their order in dend:
colors_to_use <- colors_to_use[order.dendrogram(dend)]
colors_to_use
#Convert the numbers to colors of your choice
colors_to_use[colors_to_use==1] <- "black"
colors_to_use[colors_to_use==2] <- "purple"
colors_to_use[colors_to_use==3] <- "blue"
colors_to_use[colors_to_use==4] <- "orange"
#Save thsi list as the labels_colors for your dendextend object
labels_colors(dend) <- colors_to_use
# Now each state has a color
labels_colors(dend) 

#Play around with the ylim to get a view you like and add in the legend based upon your colors above and the results from pvclust
dend %>% hang.dendrogram()%>%  plot(ylim=c(0.20,0.44))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "purple", "blue", "orange"))
result %>% text

#Save the perfected image above
tiff(file="Spawn 1b with Females.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram()%>%  plot(ylim=c(0.20,0.44))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "purple", "blue", "orange"))
result %>% text
dev.off()

#Image without statistics, but boxed around branchings with au p-values >80
tiff(file="Spawn 1b with Females boxed.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram()%>%  plot(ylim=c(0.20,0.44))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "purple", "blue", "orange"))
pvrect(result, alpha=.80, pv="au", type="geq")
dev.off()

#Image without statistics, but with bolded significant branchings with au p-values >80
tiff(file="Spawn 1b with Females bolded.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram() %>%
  pvclust_show_signif(result, signif_type =c("au"), alpha = 0.2)%>% 
  plot(ylim=c(0.20,0.44))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "purple", "blue", "orange"))
dev.off()

#Image without statistics, but with bolded significant branchings with au p-values >80, without a legend
tiff(file="Spawn 1b with Females bolded no legend.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram() %>%
  pvclust_show_signif(result, signif_type =c("au"), alpha = 0.2)%>% 
  plot(ylim=c(0.20,0.44))
dev.off()

##Plot for Spawn 1c Females
##Coloring labels based upon metadata column
#Save pvclust result as a variable and convert to a dendextend object
result <- pvclust(m.ibs.spawn1c.females.plinked.V2, method.hclust = "average",method.dist = 'euclidian', nboot = 10000)
dend <- as.dendrogram(result)

#By default, the dend has no colors to the labels so let's add some color:
#Save desired metadata category as a list of factors
status <- as.factor(sample.data.spawn1c.females.plinked$status)
#Save the desired metadata as a list of numbers (Alphabetical by default) associated with each metadata subtype
colors_to_use <- as.numeric(as.factor(sample.data.spawn1c.females.plinked[,11]))
colors_to_use
#Sort them based on their order in dend:
colors_to_use <- colors_to_use[order.dendrogram(dend)]
colors_to_use
#Convert the numbers to colors of your choice
colors_to_use[colors_to_use==1] <- "black"
colors_to_use[colors_to_use==2] <- "purple"
colors_to_use[colors_to_use==3] <- "blue"
colors_to_use[colors_to_use==4] <- "orange"
#Save thsi list as the labels_colors for your dendextend object
labels_colors(dend) <- colors_to_use
# Now each state has a color
labels_colors(dend) 

#Play around with the ylim to get a view you like and add in the legend based upon your colors above and the results from pvclust
dend %>% hang.dendrogram()%>%  plot(ylim=c(0.20,0.44))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "purple", "blue", "orange"))
result %>% text

#Save the perfected image above
tiff(file="Spawn 1c with Females.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram() %>% pvclust_show_signif(result, signif_type =c("au"))%>%  plot(ylim=c(0.20,0.44))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "purple", "blue", "orange"))
result %>% text
dev.off()

#Image without statistics, but boxed around branchings with au p-values >80
tiff(file="Spawn 1c with Females boxed.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram()%>%  plot(ylim=c(0.20,0.44))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "purple", "blue", "orange"))
pvrect(result, alpha=.80, pv="au", type="geq")
dev.off()

#Image without statistics, but with bolded significant branchings with au p-values >80
tiff(file="Spawn 1c with Females bolded.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram() %>%
  pvclust_show_signif(result, signif_type =c("au"), alpha = 0.2)%>% 
  plot(ylim=c(0.20,0.44))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "purple", "blue", "orange"))
dev.off()

#Image without statistics, but with bolded significant branchings with au p-values >80 without a legend
tiff(file="Spawn 1c with Females bolded no legend.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram() %>%
  pvclust_show_signif(result, signif_type =c("au"), alpha = 0.2)%>% 
  plot(ylim=c(0.20,0.44))
dev.off()

##Plot for Spawn 1d Females
##Coloring labels based upon metadata column
#Save pvclust result as a variable and convert to a dendextend object
result <- pvclust(m.ibs.spawn1d.females.plinked.V2, method.hclust = "average",method.dist = 'euclidian', nboot = 10000)
dend <- as.dendrogram(result)

#By default, the dend has no colors to the labels so let's add some color:
#Save desired metadata category as a list of factors
status <- as.factor(sample.data.spawn1d.females.plinked$status)
#Save the desired metadata as a list of numbers (Alphabetical by default) associated with each metadata subtype
colors_to_use <- as.numeric(as.factor(sample.data.spawn1d.females.plinked[,11]))
colors_to_use
#Sort them based on their order in dend:
colors_to_use <- colors_to_use[order.dendrogram(dend)]
colors_to_use
#Convert the numbers to colors of your choice
colors_to_use[colors_to_use==1] <- "black"
colors_to_use[colors_to_use==2] <- "purple"
colors_to_use[colors_to_use==3] <- "blue"
colors_to_use[colors_to_use==4] <- "orange"
#Save thsi list as the labels_colors for your dendextend object
labels_colors(dend) <- colors_to_use
# Now each state has a color
labels_colors(dend) 

#Play around with the ylim to get a view you like and add in the legend based upon your colors above and the results from pvclust
dend %>% hang.dendrogram()%>%  plot(ylim=c(0.20,0.44))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "purple", "blue", "orange"))
result %>% text

#Save the perfected image above
tiff(file="Spawn 1d with Females.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram()%>%  plot(ylim=c(0.20,0.44))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "purple", "blue", "orange"))
result %>% text
dev.off()

#Image without statistics, but boxed around branchings with au p-values >80
tiff(file="Spawn 1d with Females boxed.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram()%>%  plot(ylim=c(0.20,0.44))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "purple", "blue", "orange"))
pvrect(result, alpha=.80, pv="au", type="geq")
dev.off()

#Image without statistics, but with bolded significant branchings with au p-values >80
tiff(file="Spawn 1d with Females bolded.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram() %>%
  pvclust_show_signif(result, signif_type =c("au"), alpha = 0.2)%>% 
  plot(ylim=c(0.20,0.44))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "purple", "blue", "orange"))
dev.off()

#Image without statistics, but with bolded significant branchings with au p-values >80 without a legend
tiff(file="Spawn 1d with Females bolded no legend.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram() %>%
  pvclust_show_signif(result, signif_type =c("au"), alpha = 0.2)%>% 
  plot(ylim=c(0.20,0.44))
dev.off()

##Plot for Spawn 1 Females
##Coloring labels based upon metadata column
#Save pvclust result as a variable and convert to a dendextend object
result <- pvclust(m.ibs.spawn1.females.plinked.V2, method.hclust = "average",method.dist = 'euclidian', nboot = 10000)
dend <- as.dendrogram(result)

#By default, the dend has no colors to the labels so let's add some color:
#Save desired metadata category as a list of factors
status <- as.factor(sample.data.spawn1.females.plinked$status)
status
#Save the desired metadata as a list of numbers (Alphabetical by default) associated with each metadata subtype
colors_to_use <- as.numeric(as.factor(sample.data.spawn1.females.plinked[,11]))
colors_to_use
#Sort them based on their order in dend:
colors_to_use <- colors_to_use[order.dendrogram(dend)]
colors_to_use
#Convert the numbers to colors of your choice
colors_to_use[colors_to_use==1] <- "black"
colors_to_use[colors_to_use==2] <- "purple"
colors_to_use[colors_to_use==3] <- "blue"
colors_to_use[colors_to_use==4] <- "orange"
colors_to_use
#Save thsi list as the labels_colors for your dendextend object
labels_colors(dend) <- colors_to_use
# Now each state has a color
labels_colors(dend) 

#Play around with the ylim to get a view you like and add in the legend based upon your colors above and the results from pvclust
dend %>% hang.dendrogram()%>% pvclust_show_signif(result, signif_type =c("au"), alpha = 0.2)%>%  plot(ylim=c(0.1,0.53))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "purple", "blue", "orange"))
result %>% text

#Save the perfected image above
tiff(file="Spawn 1 with Females.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram()%>%  plot(ylim=c(0.1,0.53))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "purple", "blue", "orange"))
result %>% text
dev.off()

#Image without statistics, but boxed around branchings with au p-values >80
tiff(file="Spawn 1 with Females boxed.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram()%>%  plot(ylim=c(0.1,0.53))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "purple", "blue", "orange"))
pvrect(result, alpha=.80, pv="au", type="geq")
dev.off()

#Image without statistics, but with bolded significant branchings with au p-values >80
tiff(file="Spawn 1 with Females bolded.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram() %>%
  pvclust_show_signif(result, signif_type =c("au"), alpha = 0.2)%>% 
  plot(ylim=c(0.1,0.53))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "purple", "blue", "orange"))
dev.off()


####Spawn 2 + Females####

#Create Spawn 2 and Female focused dendrograms (all replicates, and separated replicates)
m.ibs.spawn2a.females.plinked.V2 <- females.spawns.ibs[c(2,9:31),c(2,9:31)]
m.ibs.spawn2b.females.plinked.V2 <- females.spawns.ibs[c(6,9:31),c(6,9:31)]
m.ibs.spawn2c.females.plinked.V2 <- females.spawns.ibs[c(7,9:31),c(7,9:31)]
m.ibs.spawn2d.females.plinked.V2 <- females.spawns.ibs[c(8,9:31),c(8,9:31)]
m.ibs.spawn2.females.plinked.V2 <- females.spawns.ibs[c(Spawn_2_Positions,9:31),c(Spawn_2_Positions,9:31)]

#Run each replicate separately and all replicates together to get idea of dendrogram
plot(pvclust(m.ibs.spawn2a.females.plinked.V2, method.hclust = "average",method.dist = 'euclidian'), main='Spawn 2a and Females')
plot(pvclust(m.ibs.spawn2b.females.plinked.V2, method.hclust = "average",method.dist = 'euclidian'), main='Spawn 2b and Females')
plot(pvclust(m.ibs.spawn2c.females.plinked.V2, method.hclust = "average",method.dist = 'euclidian'), main='Spawn 2c and Females')
plot(pvclust(m.ibs.spawn2d.females.plinked.V2, method.hclust = "average",method.dist = 'euclidian'), main='Spawn 2d and Females')
plot(pvclust(m.ibs.spawn2.females.plinked.V2, method.hclust = "average",method.dist = 'euclidian'), main='Spawn 2 and Females')

#Generally, a-d separately range from 0.42-0.33 in Height
#a-d together ranges from 0.5-0.05 in Height

#Prepare for plotting by making curated versions of metadata that only has one "offspring" sample each
sample.data.spawn2a.females.plinked <- females.spawns.data.bam.names.plinked[c(2,9:31),]
sample.data.spawn2b.females.plinked <- females.spawns.data.bam.names.plinked[c(6,9:31),]
sample.data.spawn2c.females.plinked <- females.spawns.data.bam.names.plinked[c(7,9:31),]
sample.data.spawn2d.females.plinked <- females.spawns.data.bam.names.plinked[c(8,9:31),]
sample.data.spawn2.females.plinked <- females.spawns.data.bam.names.plinked[c(Spawn_2_Positions,9:31),]

##Plot for Spawn 2a Females
#Save pvclust result, run with more permutations (nboot = 10000 instead of = 1000) than test above, as a variable and convert to a dendextend object
result <- pvclust(m.ibs.spawn2a.females.plinked.V2, method.hclust = "average",method.dist = 'euclidian', nboot = 10000)
dend <- as.dendrogram(result)

#By default, the dend has no colors to the labels so let's add some color:
#Save desired metadata category as a list of factors
status <- as.factor(sample.data.spawn2a.females.plinked$status)
status
#Save the desired metadata as a list of numbers (Alphabetical by default) associated with each metadata subtype
colors_to_use <- as.numeric(as.factor(sample.data.spawn2a.females.plinked[,11]))
colors_to_use
#Sort them based on their order in dend:
colors_to_use <- colors_to_use[order.dendrogram(dend)]
colors_to_use
#Convert the numbers to colors of your choice
colors_to_use[colors_to_use==1] <- "black"
colors_to_use[colors_to_use==2] <- "purple"
colors_to_use[colors_to_use==3] <- "blue"
colors_to_use[colors_to_use==4] <- "orange"
colors_to_use
#Save this list as the labels_colors for your dendextend object
labels_colors(dend) <- colors_to_use
# Now each state has a color
labels_colors(dend) 

#Play around with the ylim to get a view you like and add in the legend based upon your colors above and the results from pvclust
dend %>% hang.dendrogram()%>%  plot(ylim=c(0.30,0.41))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "purple", "blue", "orange"))
result %>% text

#Save the perfected image above
tiff(file="Spawn 2a with Females.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram()%>%  plot(ylim=c(0.30,0.41))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "purple", "blue", "orange"))
result %>% text
dev.off()

#Image without statistics, but boxed around branchings with au p-values >80
tiff(file="Spawn 2a with Females boxed.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram()%>%  plot(ylim=c(0.30,0.41))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "purple", "blue", "orange"))
pvrect(result, alpha=.80, pv="au", type="geq")
dev.off()

#Image without statistics, but with bolded significant branchings with au p-values >80
tiff(file="Spawn 2a with Females bolded.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram() %>%
  pvclust_show_signif(result, signif_type =c("au"), alpha = 0.2)%>% 
  plot(ylim=c(0.30,0.41))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "purple", "blue", "orange"))
dev.off()

#Image without statistics, but with bolded significant branchings with au p-values >80 without a legend
tiff(file="Spawn 2a with Females bolded no legend.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram() %>%
  pvclust_show_signif(result, signif_type =c("au"), alpha = 0.2)%>% 
  plot(ylim=c(0.30,0.41))
dev.off()

##Plot for Spawn 2b Females
##Coloring labels based upon metadata column
#Save pvclust result as a variable and convert to a dendextend object
result <- pvclust(m.ibs.spawn2b.females.plinked.V2, method.hclust = "average",method.dist = 'euclidian', nboot = 10000)
dend <- as.dendrogram(result)

#By default, the dend has no colors to the labels so let's add some color:
#Save desired metadata category as a list of factors
status <- as.factor(sample.data.spawn2b.females.plinked$status)
#Save the desired metadata as a list of numbers (Alphabetical by default) associated with each metadata subtype
colors_to_use <- as.numeric(as.factor(sample.data.spawn2b.females.plinked[,11]))
colors_to_use
#Sort them based on their order in dend:
colors_to_use <- colors_to_use[order.dendrogram(dend)]
colors_to_use
#Convert the numbers to colors of your choice
colors_to_use[colors_to_use==1] <- "black"
colors_to_use[colors_to_use==2] <- "purple"
colors_to_use[colors_to_use==3] <- "blue"
colors_to_use[colors_to_use==4] <- "orange"
#Save thsi list as the labels_colors for your dendextend object
labels_colors(dend) <- colors_to_use
# Now each state has a color
labels_colors(dend) 

#Play around with the ylim to get a view you like and add in the legend based upon your colors above and the results from pvclust
dend %>% hang.dendrogram()%>%  plot(ylim=c(0.27,0.42))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "purple", "blue", "orange"))
result %>% text

#Save the perfected image above
tiff(file="Spawn 2b with Females.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram()%>%  plot(ylim=c(0.27,0.42))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "purple", "blue", "orange"))
result %>% text
dev.off()

#Image without statistics, but boxed around branchings with au p-values >80
tiff(file="Spawn 2b with Females boxed.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram()%>%  plot(ylim=c(0.27,0.42))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "purple", "blue", "orange"))
pvrect(result, alpha=.80, pv="au", type="geq")
dev.off()

#Image without statistics, but with bolded significant branchings with au p-values >80
tiff(file="Spawn 2b with Females bolded.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram() %>%
  pvclust_show_signif(result, signif_type =c("au"), alpha = 0.2)%>% 
  plot(ylim=c(0.27,0.42))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "purple", "blue", "orange"))
dev.off()

#Image without statistics, but with bolded significant branchings with au p-values >80 without a legend
tiff(file="Spawn 2b with Females bolded no legend.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram() %>%
  pvclust_show_signif(result, signif_type =c("au"), alpha = 0.2)%>% 
  plot(ylim=c(0.27,0.42))
dev.off()

##Plot for Spawn 2c Females
##Coloring labels based upon metadata column
#Save pvclust result as a variable and convert to a dendextend object
result <- pvclust(m.ibs.spawn2c.females.plinked.V2, method.hclust = "average",method.dist = 'euclidian', nboot = 10000)
dend <- as.dendrogram(result)

#By default, the dend has no colors to the labels so let's add some color:
#Save desired metadata category as a list of factors
status <- as.factor(sample.data.spawn2c.females.plinked$status)
#Save the desired metadata as a list of numbers (Alphabetical by default) associated with each metadata subtype
colors_to_use <- as.numeric(as.factor(sample.data.spawn2c.females.plinked[,11]))
colors_to_use
#Sort them based on their order in dend:
colors_to_use <- colors_to_use[order.dendrogram(dend)]
colors_to_use
#Convert the numbers to colors of your choice
colors_to_use[colors_to_use==1] <- "black"
colors_to_use[colors_to_use==2] <- "purple"
colors_to_use[colors_to_use==3] <- "blue"
colors_to_use[colors_to_use==4] <- "orange"
#Save thsi list as the labels_colors for your dendextend object
labels_colors(dend) <- colors_to_use
# Now each state has a color
labels_colors(dend) 

#Play around with the ylim to get a view you like and add in the legend based upon your colors above and the results from pvclust
dend %>% hang.dendrogram()%>%  plot(ylim=c(0.27,0.42))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "purple", "blue", "orange"))
result %>% text

#Save the perfected image above
tiff(file="Spawn 2c with Females.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram()%>%  plot(ylim=c(0.27,0.42))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "purple", "blue", "orange"))
result %>% text
dev.off()

#Image without statistics, but boxed around branchings with au p-values >80
tiff(file="Spawn 2c with Females boxed.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram()%>%  plot(ylim=c(0.27,0.42))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "purple", "blue", "orange"))
pvrect(result, alpha=.80, pv="au", type="geq")
dev.off()

#Image without statistics, but with bolded significant branchings with au p-values >80
tiff(file="Spawn 2c with Females bolded.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram() %>%
  pvclust_show_signif(result, signif_type =c("au"), alpha = 0.2)%>% 
  plot(ylim=c(0.27,0.42))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "purple", "blue", "orange"))
dev.off()

#Image without statistics, but with bolded significant branchings with au p-values >80 without a legend
tiff(file="Spawn 2c with Females bolded no legend.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram() %>%
  pvclust_show_signif(result, signif_type =c("au"), alpha = 0.2)%>% 
  plot(ylim=c(0.27,0.42))
dev.off()

##Plot for Spawn 2d Females
##Coloring labels based upon metadata column
#Save pvclust result as a variable and convert to a dendextend object
result <- pvclust(m.ibs.spawn2d.females.plinked.V2, method.hclust = "average",method.dist = 'euclidian', nboot = 10000)
dend <- as.dendrogram(result)

#By default, the dend has no colors to the labels so let's add some color:
#Save desired metadata category as a list of factors
status <- as.factor(sample.data.spawn2d.females.plinked$status)
#Save the desired metadata as a list of numbers (Alphabetical by default) associated with each metadata subtype
colors_to_use <- as.numeric(as.factor(sample.data.spawn2d.females.plinked[,11]))
colors_to_use
#Sort them based on their order in dend:
colors_to_use <- colors_to_use[order.dendrogram(dend)]
colors_to_use
#Convert the numbers to colors of your choice
colors_to_use[colors_to_use==1] <- "black"
colors_to_use[colors_to_use==2] <- "purple"
colors_to_use[colors_to_use==3] <- "blue"
colors_to_use[colors_to_use==4] <- "orange"
#Save thsi list as the labels_colors for your dendextend object
labels_colors(dend) <- colors_to_use
# Now each state has a color
labels_colors(dend) 

#Play around with the ylim to get a view you like and add in the legend based upon your colors above and the results from pvclust
dend %>% hang.dendrogram()%>%  plot(ylim=c(0.27,0.42))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "purple", "blue", "orange"))
result %>% text

#Save the perfected image above
tiff(file="Spawn 2d with Females.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram()%>%  plot(ylim=c(0.27,0.42))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "purple", "blue", "orange"))
result %>% text
dev.off()

#Image without statistics, but boxed around branchings with au p-values >80
tiff(file="Spawn 2d with Females boxed.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram()%>%  plot(ylim=c(0.27,0.42))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "purple", "blue", "orange"))
pvrect(result, alpha=.80, pv="au", type="geq")
dev.off()

#Image without statistics, but with bolded significant branchings with au p-values >80 without a legend
tiff(file="Spawn 2d with Females bolded no legend.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram() %>%
  pvclust_show_signif(result, signif_type =c("au"), alpha = 0.2)%>% 
  plot(ylim=c(0.27,0.42))
dev.off()

##Plot for Spawn 2 Females
##Coloring labels based upon metadata column
#Save pvclust result as a variable and convert to a dendextend object
result <- pvclust(m.ibs.spawn2.females.plinked.V2, method.hclust = "average",method.dist = 'euclidian', nboot = 10000)
dend <- as.dendrogram(result)

#By default, the dend has no colors to the labels so let's add some color:
#Save desired metadata category as a list of factors
status <- as.factor(sample.data.spawn2.females.plinked$status)
status
#Save the desired metadata as a list of numbers (Alphabetical by default) associated with each metadata subtype
colors_to_use <- as.numeric(as.factor(sample.data.spawn2.females.plinked[,11]))
colors_to_use
#Sort them based on their order in dend:
colors_to_use <- colors_to_use[order.dendrogram(dend)]
colors_to_use
#Convert the numbers to colors of your choice
colors_to_use[colors_to_use==1] <- "black"
colors_to_use[colors_to_use==2] <- "purple"
colors_to_use[colors_to_use==3] <- "blue"
colors_to_use[colors_to_use==4] <- "orange"
colors_to_use
#Save thsi list as the labels_colors for your dendextend object
labels_colors(dend) <- colors_to_use
# Now each state has a color
labels_colors(dend) 

#Play around with the ylim to get a view you like and add in the legend based upon your colors above and the results from pvclust
dend %>% hang.dendrogram()%>%  plot(ylim=c(0.1,0.53))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "purple", "blue", "orange"))
result %>% text

#Save the perfected image above
tiff(file="Spawn 2 with Females.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram()%>%  plot(ylim=c(0.1,0.53))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "purple", "blue", "orange"))
result %>% text
dev.off()

#Image without statistics, but boxed around branchings with au p-values >80
tiff(file="Spawn 2 with Females boxed.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram()%>%  plot(ylim=c(0.1,0.53))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "purple", "blue", "orange"))
pvrect(result, alpha=.80, pv="au", type="geq")
pvclust_show_signif(dend, result, signif_type = "au")
dev.off()

#Image without statistics, but with bolded significant branchings with au p-values >80
tiff(file="Spawn 2 with Females bolded.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram() %>%
  pvclust_show_signif(result, signif_type =c("au"), alpha = 0.2)%>% 
  plot(ylim=c(0.1,0.53))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "purple", "blue", "orange"))
dev.off()

####Males Prep####

#Filter Female samples from metadata
males.spawns.data.bam.names.plinked <- filter(data.bam.names.plinked, sex!="F")

#Subset main mibs for just females and spawns
male.spawns.ibs <- adults.spawns.ibs[c(1:8, Male_Positions), c(1:8, Male_Positions)]

####Spawn 1 + Males####

#Create Spawn 1 and Female focused dendrograms (all replicates, and separated replicates)
m.ibs.spawn1a.males.plinked.V2 <- male.spawns.ibs[c(1,9:17),c(1,9:17)]
m.ibs.spawn1b.males.plinked.V2 <- male.spawns.ibs[c(3,9:17),c(3,9:17)]
m.ibs.spawn1c.males.plinked.V2 <- male.spawns.ibs[c(4,9:17),c(4,9:17)]
m.ibs.spawn1d.males.plinked.V2 <- male.spawns.ibs[c(5,9:17),c(5,9:17)]
m.ibs.spawn1.males.plinked.V2 <- male.spawns.ibs[c(Spawn_1_Positions,9:17),c(Spawn_1_Positions,9:17)]


#Run each replicate separately and all replicates together to get idea of dendrogram
plot(pvclust(m.ibs.spawn1a.males.plinked.V2, method.hclust = "average",method.dist = 'euclidian'), main='Spawn 1a and Males')
plot(pvclust(m.ibs.spawn1b.males.plinked.V2, method.hclust = "average",method.dist = 'euclidian'), main='Spawn 1b and Males')
plot(pvclust(m.ibs.spawn1c.males.plinked.V2, method.hclust = "average",method.dist = 'euclidian'), main='Spawn 1c and Males')
plot(pvclust(m.ibs.spawn1d.males.plinked.V2, method.hclust = "average",method.dist = 'euclidian'), main='Spawn 1d and Males')
plot(pvclust(m.ibs.spawn1.males.plinked.V2, method.hclust = "average",method.dist = 'euclidian'), main='Spawn 1 and Males')

#Generally, a-d separately range from 0.42-0.30 in Height
#a-d together ranges from 0.5-0.05 in Height

#Prepare for plotting by making curated versions of metadata that only has one "offspring" sample each
sample.data.spawn1a.males.plinked <- males.spawns.data.bam.names.plinked[c(1,9:17),]
sample.data.spawn1b.males.plinked <- males.spawns.data.bam.names.plinked[c(3,9:17),]
sample.data.spawn1c.males.plinked <- males.spawns.data.bam.names.plinked[c(4,9:17),]
sample.data.spawn1d.males.plinked <- males.spawns.data.bam.names.plinked[c(5,9:17),]
sample.data.spawn1.males.plinked <- males.spawns.data.bam.names.plinked[c(Spawn_1_Positions,9:17),]

##Plot for Spawn 1a Males
##Coloring labels based upon metadata column
#Save pvclust result as a variable and convert to a dendextend object
result <- pvclust(m.ibs.spawn1a.males.plinked.V2, method.hclust = "average",method.dist = 'euclidian', nboot = 10000)
dend <- as.dendrogram(result)

#By default, the dend has no colors to the labels so let's add some color:
#Save desired metadata category as a list of factors
status <- as.factor(sample.data.spawn1a.males.plinked$status)
status
#Save the desired metadata as a list of numbers (Alphabetical by default) associated with each metadata subtype
colors_to_use <- as.numeric(as.factor(sample.data.spawn1a.males.plinked[,11]))
colors_to_use
#Sort them based on their order in dend:
colors_to_use <- colors_to_use[order.dendrogram(dend)]
colors_to_use
#Convert the numbers to colors of your choice
colors_to_use[colors_to_use==1] <- "black"
colors_to_use[colors_to_use==2] <- "blue"
colors_to_use[colors_to_use==3] <- "orange"
colors_to_use
#Save this list as the labels_colors for your dendextend object
labels_colors(dend) <- colors_to_use
# Now each state has a color
labels_colors(dend) 

#Play around with the ylim to get a view you like and add in the legend based upon your colors above and the results from pvclust
dend %>% hang.dendrogram()%>%  plot(ylim=c(0.28,0.42))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "blue", "orange"))
result %>% text

#Save the perfected image above
tiff(file="Spawn 1a with Males.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram()%>%  plot(ylim=c(0.28,0.42))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "blue", "orange"))
result %>% text
dev.off()

#Image without statistics, but boxed around branchings with au p-values >80
tiff(file="Spawn 1a with Males boxed.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram()%>%  plot(ylim=c(0.28,0.42))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "blue", "orange"))
pvrect(result, alpha=.80, pv="au", type="geq")
dev.off()

#Image without statistics, but with bolded significant branchings with au p-values >80 without a legend
tiff(file="Spawn 1a with Males bolded no legend.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram() %>%
  pvclust_show_signif(result, signif_type =c("au"), alpha = 0.2)%>% 
  plot(ylim=c(0.28,0.42))
dev.off()

##Plot for Spawn 1b Males
##Coloring labels based upon metadata column
#Save pvclust result as a variable and convert to a dendextend object
result <- pvclust(m.ibs.spawn1b.males.plinked.V2, method.hclust = "average",method.dist = 'euclidian', nboot = 10000)
dend <- as.dendrogram(result)

#By default, the dend has no colors to the labels so let's add some color:
#Save desired metadata category as a list of factors
status <- as.factor(sample.data.spawn1b.males.plinked$status)
#Save the desired metadata as a list of numbers (Alphabetical by default) associated with each metadata subtype
colors_to_use <- as.numeric(as.factor(sample.data.spawn1b.males.plinked[,11]))
colors_to_use
#Sort them based on their order in dend:
colors_to_use <- colors_to_use[order.dendrogram(dend)]
colors_to_use
#Convert the numbers to colors of your choice
colors_to_use[colors_to_use==1] <- "black"
colors_to_use[colors_to_use==2] <- "blue"
colors_to_use[colors_to_use==3] <- "orange"
#Save thsi list as the labels_colors for your dendextend object
labels_colors(dend) <- colors_to_use
# Now each state has a color
labels_colors(dend) 

#Play around with the ylim to get a view you like and add in the legend based upon your colors above and the results from pvclust
dend %>% hang.dendrogram()%>%  plot(ylim=c(0.29,0.42))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "blue", "orange"))
result %>% text

#Save the perfected image above
tiff(file="Spawn 1b with Males.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram()%>%  plot(ylim=c(0.29,0.42))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "blue", "orange"))
result %>% text
dev.off()

#Image without statistics, but boxed around branchings with au p-values >80
tiff(file="Spawn 1b with Males boxed.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram()%>%  plot(ylim=c(0.29,0.42))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "blue", "orange"))
pvrect(result, alpha=.80, pv="au", type="geq")
dev.off()

#Image without statistics, but with bolded significant branchings with au p-values >80 without a legend
tiff(file="Spawn 1b with Males bolded no legend.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram() %>%
  pvclust_show_signif(result, signif_type =c("au"), alpha = 0.2)%>% 
  plot(ylim=c(0.29,0.42))
dev.off()

##Plot for Spawn 1c Males
##Coloring labels based upon metadata column
#Save pvclust result as a variable and convert to a dendextend object
result <- pvclust(m.ibs.spawn1c.males.plinked.V2, method.hclust = "average",method.dist = 'euclidian', nboot = 10000)
dend <- as.dendrogram(result)

#By default, the dend has no colors to the labels so let's add some color:
#Save desired metadata category as a list of factors
status <- as.factor(sample.data.spawn1c.males.plinked$status)
#Save the desired metadata as a list of numbers (Alphabetical by default) associated with each metadata subtype
colors_to_use <- as.numeric(as.factor(sample.data.spawn1c.males.plinked[,11]))
colors_to_use
#Sort them based on their order in dend:
colors_to_use <- colors_to_use[order.dendrogram(dend)]
colors_to_use
#Convert the numbers to colors of your choice
colors_to_use[colors_to_use==1] <- "black"
colors_to_use[colors_to_use==2] <- "blue"
colors_to_use[colors_to_use==3] <- "orange"
#Save thsi list as the labels_colors for your dendextend object
labels_colors(dend) <- colors_to_use
# Now each state has a color
labels_colors(dend) 

#Play around with the ylim to get a view you like and add in the legend based upon your colors above and the results from pvclust
dend %>% hang.dendrogram()%>%  plot(ylim=c(0.27,0.42))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "blue", "orange"))
result %>% text

#Save the perfected image above
tiff(file="Spawn 1c with Males.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram()%>%  plot(ylim=c(0.27,0.42))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "blue", "orange"))
result %>% text
dev.off()

#Image without statistics, but boxed around branchings with au p-values >80
tiff(file="Spawn 1c with Males boxed.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram()%>%  plot(ylim=c(0.27,0.42))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "blue", "orange"))
pvrect(result, alpha=.80, pv="au", type="geq")
dev.off()

#Image without statistics, but with bolded significant branchings with au p-values >80 without a legend
tiff(file="Spawn 1c with Males bolded no legend.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram() %>%
  pvclust_show_signif(result, signif_type =c("au"), alpha = 0.2)%>% 
  plot(ylim=c(0.27,0.42))
dev.off()

##Plot for Spawn 1d Males
##Coloring labels based upon metadata column
#Save pvclust result as a variable and convert to a dendextend object
result <- pvclust(m.ibs.spawn1d.males.plinked.V2, method.hclust = "average",method.dist = 'euclidian', nboot = 10000)
dend <- as.dendrogram(result)

#By default, the dend has no colors to the labels so let's add some color:
#Save desired metadata category as a list of factors
status <- as.factor(sample.data.spawn1d.males.plinked$status)
#Save the desired metadata as a list of numbers (Alphabetical by default) associated with each metadata subtype
colors_to_use <- as.numeric(as.factor(sample.data.spawn1d.males.plinked[,11]))
colors_to_use
#Sort them based on their order in dend:
colors_to_use <- colors_to_use[order.dendrogram(dend)]
colors_to_use
#Convert the numbers to colors of your choice
colors_to_use[colors_to_use==1] <- "black"
colors_to_use[colors_to_use==2] <- "blue"
colors_to_use[colors_to_use==3] <- "orange"
#Save thsi list as the labels_colors for your dendextend object
labels_colors(dend) <- colors_to_use
# Now each state has a color
labels_colors(dend) 

#Play around with the ylim to get a view you like and add in the legend based upon your colors above and the results from pvclust
dend %>% hang.dendrogram()%>%  plot(ylim=c(0.27,0.42))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "blue", "orange"))
result %>% text

#Save the perfected image above
tiff(file="Spawn 1d with Males.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram()%>%  plot(ylim=c(0.27,0.42))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "blue", "orange"))
result %>% text
dev.off()

#Image without statistics, but boxed around branchings with au p-values >80
tiff(file="Spawn 1d with Males boxed.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram()%>%  plot(ylim=c(0.27,0.42))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "blue", "orange"))
pvrect(result, alpha=.80, pv="au", type="geq")
dev.off()

#Image without statistics, but with bolded significant branchings with au p-values >80 without a legend
tiff(file="Spawn 1d with Males bolded no legend.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram() %>%
  pvclust_show_signif(result, signif_type =c("au"), alpha = 0.2)%>% 
  plot(ylim=c(0.27,0.42))
dev.off()

##Plot for Spawn 1 Males
##Coloring labels based upon metadata column
#Save pvclust result as a variable and convert to a dendextend object
result <- pvclust(m.ibs.spawn1.males.plinked.V2, method.hclust = "average",method.dist = 'euclidian', nboot = 10000)
dend <- as.dendrogram(result)

#By default, the dend has no colors to the labels so let's add some color:
#Save desired metadata category as a list of factors
status <- as.factor(sample.data.spawn1.males.plinked$status)
status
#Save the desired metadata as a list of numbers (Alphabetical by default) associated with each metadata subtype
colors_to_use <- as.numeric(as.factor(sample.data.spawn1.males.plinked[,11]))
colors_to_use
#Sort them based on their order in dend:
colors_to_use <- colors_to_use[order.dendrogram(dend)]
colors_to_use
#Convert the numbers to colors of your choice
colors_to_use[colors_to_use==1] <- "black"
colors_to_use[colors_to_use==2] <- "blue"
colors_to_use[colors_to_use==3] <- "orange"
colors_to_use
#Save thsi list as the labels_colors for your dendextend object
labels_colors(dend) <- colors_to_use
# Now each state has a color
labels_colors(dend) 

#Play around with the ylim to get a view you like and add in the legend based upon your colors above and the results from pvclust
dend %>% hang.dendrogram()%>%  plot(ylim=c(0.1,0.53))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "blue", "orange"))
result %>% text

#Save the perfected image above
tiff(file="Spawn 1 with Males.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram()%>%  plot(ylim=c(0.1,0.53))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "blue", "orange"))
result %>% text
dev.off()

#Image without statistics, but boxed around branchings with au p-values >80
tiff(file="Spawn 1 with Males boxed.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram()%>%  plot(ylim=c(0.1,0.53))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "blue", "orange"))
pvrect(result, alpha=.80, pv="au", type="geq")
dev.off()

#Image without statistics, but with bolded significant branchings with au p-values >80 without a legend
tiff(file="Spawn 1 with Males bolded no legend.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram() %>%
  pvclust_show_signif(result, signif_type =c("au"), alpha = 0.2)%>% 
  plot(ylim=c(0.1,0.53))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "blue", "orange"))
dev.off()

####Spawn 2 + Males####

m.ibs.spawn2a.males.plinked.V2 <- male.spawns.ibs[c(2,9:17),c(2,9:17)]
m.ibs.spawn2b.males.plinked.V2 <- male.spawns.ibs[c(6,9:17),c(6,9:17)]
m.ibs.spawn2c.males.plinked.V2 <- male.spawns.ibs[c(7,9:17),c(7,9:17)]
m.ibs.spawn2d.males.plinked.V2 <- male.spawns.ibs[c(8,9:17),c(8,9:17)]
m.ibs.spawn2.males.plinked.V2 <- male.spawns.ibs[c(Spawn_2_Positions,9:17),c(Spawn_2_Positions,9:17)]

plot(pvclust(m.ibs.spawn2a.males.plinked.V2, method.hclust = "average",method.dist = 'euclidian'), main='Spawn 1a and Males')
plot(pvclust(m.ibs.spawn2b.males.plinked.V2, method.hclust = "average",method.dist = 'euclidian'), main='Spawn 1b and Males')
plot(pvclust(m.ibs.spawn2c.males.plinked.V2, method.hclust = "average",method.dist = 'euclidian'), main='Spawn 1c and Males')
plot(pvclust(m.ibs.spawn2d.males.plinked.V2, method.hclust = "average",method.dist = 'euclidian'), main='Spawn 1d and Males')
plot(pvclust(m.ibs.spawn2.males.plinked.V2, method.hclust = "average",method.dist = 'euclidian'), main='Spawn 2 and Males')

#Generally, a-d separately range from 0.42-0.26 in Height
#a-d together ranges from 0.5-0.05 in Height

#Prepare for plotting by making curated versions of metadata that only has one "offspring" sample each
sample.data.spawn2a.males.plinked <- males.spawns.data.bam.names.plinked[c(2,9:17),]
sample.data.spawn2b.males.plinked <- males.spawns.data.bam.names.plinked[c(6,9:17),]
sample.data.spawn2c.males.plinked <- males.spawns.data.bam.names.plinked[c(7,9:17),]
sample.data.spawn2d.males.plinked <- males.spawns.data.bam.names.plinked[c(8,9:17),]
sample.data.spawn2.males.plinked <- males.spawns.data.bam.names.plinked[c(Spawn_2_Positions,9:17),]

##Plot for Spawn 2a Males
##Coloring labels based upon metadata column
#Save pvclust result as a variable and convert to a dendextend object
result <- pvclust(m.ibs.spawn2a.males.plinked.V2, method.hclust = "average",method.dist = 'euclidian', nboot = 10000)
dend <- as.dendrogram(result)

#By default, the dend has no colors to the labels so let's add some color:
#Save desired metadata category as a list of factors
status <- as.factor(sample.data.spawn2a.males.plinked$status)
status
#Save the desired metadata as a list of numbers (Alphabetical by default) associated with each metadata subtype
colors_to_use <- as.numeric(as.factor(sample.data.spawn2a.males.plinked[,11]))
colors_to_use
#Sort them based on their order in dend:
colors_to_use <- colors_to_use[order.dendrogram(dend)]
colors_to_use
#Convert the numbers to colors of your choice
colors_to_use[colors_to_use==1] <- "black"
colors_to_use[colors_to_use==2] <- "blue"
colors_to_use[colors_to_use==3] <- "orange"
colors_to_use
#Save this list as the labels_colors for your dendextend object
labels_colors(dend) <- colors_to_use
# Now each state has a color
labels_colors(dend) 

#Play around with the ylim to get a view you like and add in the legend based upon your colors above and the results from pvclust
dend %>% hang.dendrogram()%>%  plot(ylim=c(0.24,0.44))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "blue", "orange"))
result %>% text

#Save the perfected image above
tiff(file="Spawn 2a with Males.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram()%>%  plot(ylim=c(0.24,0.44))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "blue", "orange"))
result %>% text
dev.off()

#Image without statistics, but boxed around branchings with au p-values >80
tiff(file="Spawn 2a with Males boxed.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram()%>%  plot(ylim=c(0.24,0.4))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "blue", "orange"))
pvrect(result, alpha=.80, pv="au", type="geq")
dev.off()

#Image without statistics, but with bolded significant branchings with au p-values >80
tiff(file="Spawn 2a with Males bolded.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram() %>%
  pvclust_show_signif(result, signif_type =c("au"), alpha = 0.2)%>% 
  plot(ylim=c(0.24,0.4))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "blue", "orange"))
dev.off()

#Image without statistics, but with bolded significant branchings with au p-values >80 without a legend
tiff(file="Spawn 2a with Males bolded no legend.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram() %>%
  pvclust_show_signif(result, signif_type =c("au"), alpha = 0.2)%>% 
  plot(ylim=c(0.24,0.4))
dev.off()

##Plot for Spawn 2b Males
##Coloring labels based upon metadata column
#Save pvclust result as a variable and convert to a dendextend object
result <- pvclust(m.ibs.spawn2b.males.plinked.V2, method.hclust = "average",method.dist = 'euclidian', nboot = 10000)
dend <- as.dendrogram(result)

#By default, the dend has no colors to the labels so let's add some color:
#Save desired metadata category as a list of factors
status <- as.factor(sample.data.spawn2b.males.plinked$status)
#Save the desired metadata as a list of numbers (Alphabetical by default) associated with each metadata subtype
colors_to_use <- as.numeric(as.factor(sample.data.spawn2b.males.plinked[,11]))
colors_to_use
#Sort them based on their order in dend:
colors_to_use <- colors_to_use[order.dendrogram(dend)]
colors_to_use
#Convert the numbers to colors of your choice
colors_to_use[colors_to_use==1] <- "black"
colors_to_use[colors_to_use==2] <- "blue"
colors_to_use[colors_to_use==3] <- "orange"
#Save thsi list as the labels_colors for your dendextend object
labels_colors(dend) <- colors_to_use
# Now each state has a color
labels_colors(dend) 

#Play around with the ylim to get a view you like and add in the legend based upon your colors above and the results from pvclust
dend %>% hang.dendrogram()%>%  plot(ylim=c(0.27,0.42))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "blue", "orange"))
result %>% text

#Save the perfected image above
tiff(file="Spawn 2b with Males.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram()%>%  plot(ylim=c(0.27,0.42))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "blue", "orange"))
result %>% text
dev.off()

#Image without statistics, but boxed around branchings with au p-values >80
tiff(file="Spawn 2b with Males boxed.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram()%>%  plot(ylim=c(0.27,0.42))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "blue", "orange"))
pvrect(result, alpha=.80, pv="au", type="geq")
dev.off()

#Image without statistics, but with bolded significant branchings with au p-values >80
tiff(file="Spawn 2b with Males bolded.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram() %>%
  pvclust_show_signif(result, signif_type =c("au"), alpha = 0.2)%>% 
  plot(ylim=c(0.27,0.42))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "blue", "orange"))
dev.off()

#Image without statistics, but with bolded significant branchings with au p-values >80 without a legend
tiff(file="Spawn 2b with Males bolded no legend.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram() %>%
  pvclust_show_signif(result, signif_type =c("au"), alpha = 0.2)%>% 
  plot(ylim=c(0.27,0.42))
dev.off()

##Plot for Spawn 2c Males
##Coloring labels based upon metadata column
#Save pvclust result as a variable and convert to a dendextend object
result <- pvclust(m.ibs.spawn2c.males.plinked.V2, method.hclust = "average",method.dist = 'euclidian', nboot = 10000)
dend <- as.dendrogram(result)

#By default, the dend has no colors to the labels so let's add some color:
#Save desired metadata category as a list of factors
status <- as.factor(sample.data.spawn2c.males.plinked$status)
#Save the desired metadata as a list of numbers (Alphabetical by default) associated with each metadata subtype
colors_to_use <- as.numeric(as.factor(sample.data.spawn2c.males.plinked[,11]))
colors_to_use
#Sort them based on their order in dend:
colors_to_use <- colors_to_use[order.dendrogram(dend)]
colors_to_use
#Convert the numbers to colors of your choice
colors_to_use[colors_to_use==1] <- "black"
colors_to_use[colors_to_use==2] <- "blue"
colors_to_use[colors_to_use==3] <- "orange"
#Save thsi list as the labels_colors for your dendextend object
labels_colors(dend) <- colors_to_use
# Now each state has a color
labels_colors(dend) 

#Play around with the ylim to get a view you like and add in the legend based upon your colors above and the results from pvclust
dend %>% hang.dendrogram()%>%  plot(ylim=c(0.25,0.44))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "blue", "orange"))
result %>% text

#Save the perfected image above
tiff(file="Spawn 2c with Males.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram()%>%  plot(ylim=c(0.25,0.43))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "blue", "orange"))
result %>% text
dev.off()

#Image without statistics, but boxed around branchings with au p-values >80
tiff(file="Spawn 2c with Males boxed.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram()%>%  plot(ylim=c(0.25,0.43))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "blue", "orange"))
pvrect(result, alpha=.80, pv="au", type="geq")
dev.off()

#Image without statistics, but with bolded significant branchings with au p-values >80
tiff(file="Spawn 2c with Males bolded.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram() %>%
  pvclust_show_signif(result, signif_type =c("au"), alpha = 0.2)%>% 
  plot(ylim=c(0.25,0.43))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "blue", "orange"))
dev.off()

#Image without statistics, but with bolded significant branchings with au p-values >80 without a legend
tiff(file="Spawn 2c with Males bolded no legend.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram() %>%
  pvclust_show_signif(result, signif_type =c("au"), alpha = 0.2)%>% 
  plot(ylim=c(0.25,0.43))
dev.off()

##Plot for Spawn 2d Males
##Coloring labels based upon metadata column
#Save pvclust result as a variable and convert to a dendextend object
result <- pvclust(m.ibs.spawn2d.males.plinked.V2, method.hclust = "average",method.dist = 'euclidian', nboot = 10000)
dend <- as.dendrogram(result)

#By default, the dend has no colors to the labels so let's add some color:
#Save desired metadata category as a list of factors
status <- as.factor(sample.data.spawn2d.males.plinked$status)
#Save the desired metadata as a list of numbers (Alphabetical by default) associated with each metadata subtype
colors_to_use <- as.numeric(as.factor(sample.data.spawn2d.males.plinked[,11]))
colors_to_use
#Sort them based on their order in dend:
colors_to_use <- colors_to_use[order.dendrogram(dend)]
colors_to_use
#Convert the numbers to colors of your choice
colors_to_use[colors_to_use==1] <- "black"
colors_to_use[colors_to_use==2] <- "blue"
colors_to_use[colors_to_use==3] <- "orange"
#Save thsi list as the labels_colors for your dendextend object
labels_colors(dend) <- colors_to_use
# Now each state has a color
labels_colors(dend) 

#Play around with the ylim to get a view you like and add in the legend based upon your colors above and the results from pvclust
dend %>% hang.dendrogram()%>%  plot(ylim=c(0.25,0.43))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "blue", "orange"))
result %>% text

#Save the perfected image above
tiff(file="Spawn 2d with Males.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram()%>%  plot(ylim=c(0.27,0.42))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "blue", "orange"))
result %>% text
dev.off()

#Image without statistics, but boxed around branchings with au p-values >80
tiff(file="Spawn 2d with Males boxed.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram()%>%  plot(ylim=c(0.27,0.42))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "blue", "orange"))
pvrect(result, alpha=.80, pv="au", type="geq")
dev.off()

#Image without statistics, but with bolded significant branchings with au p-values >80
tiff(file="Spawn 2d with Males bolded.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram() %>%
  pvclust_show_signif(result, signif_type =c("au"), alpha = 0.2)%>% 
  plot(ylim=c(0.27,0.42))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "blue", "orange"))
dev.off()

#Image without statistics, but with bolded significant branchings with au p-values >80 without a legend
tiff(file="Spawn 2d with Males bolded no legend.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram() %>%
  pvclust_show_signif(result, signif_type =c("au"), alpha = 0.2)%>% 
  plot(ylim=c(0.27,0.42))
dev.off()

##Plot for Spawn 2 Males
##Coloring labels based upon metadata column
#Save pvclust result as a variable and convert to a dendextend object
result <- pvclust(m.ibs.spawn2.males.plinked.V2, method.hclust = "average",method.dist = 'euclidian', nboot = 10000)
dend <- as.dendrogram(result)

#By default, the dend has no colors to the labels so let's add some color:
#Save desired metadata category as a list of factors
status <- as.factor(sample.data.spawn2.males.plinked$status)
status
#Save the desired metadata as a list of numbers (Alphabetical by default) associated with each metadata subtype
colors_to_use <- as.numeric(as.factor(sample.data.spawn2.males.plinked[,11]))
colors_to_use
#Sort them based on their order in dend:
colors_to_use <- colors_to_use[order.dendrogram(dend)]
colors_to_use
#Convert the numbers to colors of your choice
colors_to_use[colors_to_use==1] <- "black"
colors_to_use[colors_to_use==2] <- "blue"
colors_to_use[colors_to_use==3] <- "orange"
colors_to_use
#Save thsi list as the labels_colors for your dendextend object
labels_colors(dend) <- colors_to_use
# Now each state has a color
labels_colors(dend) 

#Play around with the ylim to get a view you like and add in the legend based upon your colors above and the results from pvclust
dend %>% hang.dendrogram()%>%  plot(ylim=c(0.1,0.53))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "blue", "orange"))
result %>% text

#Save the perfected image above
tiff(file="Spawn 2 with Males.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram()%>%  plot(ylim=c(0.1,0.53))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "blue", "orange"))
result %>% text
dev.off()

#Image without statistics, but boxed around branchings with au p-values >80
tiff(file="Spawn 2 with Males boxed.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram()%>%  plot(ylim=c(0.1,0.53))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "blue", "orange"))
pvrect(result, alpha=.80, pv="au", type="geq")
dev.off()

#Image without statistics, but with bolded significant branchings with au p-values >80
tiff(file="Spawn 2 with Males bolded.tiff", units = "in", height = 10, width = 10, res=300)
dend %>% hang.dendrogram() %>%
  pvclust_show_signif(result, signif_type =c("au"), alpha = 0.2)%>% 
  plot(ylim=c(0.1,0.53))
legend("top", horiz="TRUE", legend = levels(status), fill = c("black", "blue", "orange"))
dev.off()












####Plink - Mandel Test Prep####

##Females

#Extract ids
female_ids <- filter(data.bam.names.plinked, sex=="F")$bam_id 

#Remove extra bits
female_ids <- gsub(".fq.bt2.bam", "", female_ids)

#Check on length
length(female_ids)
#23/23

##Females

#Extract ids
male_ids <- filter(data.bam.names.plinked, sex=="M")$bam_id 

#Remove extra bits
male_ids <- gsub(".fq.bt2.bam", "", male_ids)

#Check on length
length(male_ids)
#9/9

#Create a dataframe of all male and female combos
male_female_combos <- 
  expand.grid(male_ids, female_ids)
#9*23 = 207

#Change columnn names
colnames(male_female_combos) <- c("male_ids", "female_ids")


#Load the draft mandel dataframe
update_parents_S11_S9 <-
  read_delim("../mendel_parentage/update_parents_S11_S9.txt")

#Create a loop to make all verisions of the parentage trios
for(row in 4:nrow(male_female_combos)){
  male_id <- as.character(male_female_combos[row, "male_ids"])
  female_id <- as.character(male_female_combos[row, "female_ids"])
  temp_df <- update_parents_S11_S9
  temp_df$`Within-family ID of father` <-
    gsub("S11" , male_id, temp_df$`Within-family ID of father`)
  temp_df$`Within-family ID of mother` <-
    gsub("S9" , female_id, temp_df$`Within-family ID of mother`)
  write_delim(temp_df, 
              paste0("../mendel_parentage/update_parents_", 
                     male_id, "_", female_id, ".txt"))
}

ls(pattern, "update_parents_", "../mendel_parentage")




####Mendelian Errors Summarizing###

#Upload curated mendel results
all_families_imendel <- read_tsv("../WL_HBOI_R_1_hwe_0.1/all_families_imendel.tsv")

#Add in paternal pit tag id
all_families_imendel_pit_tag <-
  merge(select(data.bam.names.plinked, sample_id, pit_tag),
        all_families_imendel,
        by.x = "sample_id",
        by.y = "paternal_id")

#Change column names
colnames(all_families_imendel_pit_tag) <- c("paternal_id", "paternal_pit_tag", "maternal_id", "offspring_id", "num_mendelian_errors")

#Add in maternal pit tag id
all_families_imendel_pit_tag <-
  merge(select(data.bam.names.plinked, sample_id, pit_tag),
        all_families_imendel_pit_tag,
        by.x = "sample_id",
        by.y = "maternal_id")

#Change column names
colnames(all_families_imendel_pit_tag) <- c("maternal_id", "maternal_pit_tag", "paternal_id", "paternal_pit_tag", "offspring_id", "num_mendelian_errors")

#Create a column that combines paternal and maternal ids into one column
all_families_imendel_pit_tag <-
  all_families_imendel_pit_tag %>%
  unite("parentage_pair_id", c(paternal_id, maternal_id), sep = "_", remove = FALSE)

#Create a column that combines paternal and maternal pit tag ids into one column
all_families_imendel_pit_tag <-
  all_families_imendel_pit_tag %>%
  unite("parentage_pair_pit_tag", c(paternal_pit_tag, maternal_pit_tag), sep = "_", remove = FALSE)

#Add in offspring "pit tag id"
all_families_imendel_pit_tag <-
  merge(select(data.bam.names.plinked, sample_id, pit_tag),
        all_families_imendel_pit_tag,
        by.x = "sample_id",
        by.y = "offspring_id")

#Change column names
colnames(all_families_imendel_pit_tag) <- 
  c("offspring_id", "offspring_pit_tag", "parentage_pair_id", "maternal_id", "parentage_pair_pit_tag", "maternal_pit_tag", "paternal_id", "paternal_pit_tag", "num_mendelian_errors")

#Change order of columns
all_families_imendel_pit_tag <- 
  all_families_imendel_pit_tag %>%
  subset(select = c("parentage_pair_id", "parentage_pair_pit_tag", 
                    "paternal_id", "paternal_pit_tag",
                    "maternal_id", "maternal_pit_tag",
                    "offspring_id", "offspring_pit_tag",
                    "num_mendelian_errors"))


#Create blank dataframes to hold loop data
all_families_imendel_summ_mat <- data.frame()
all_families_imendel_summ_pat <- data.frame()

#For each offspring get the mean and standard deviation for all males and females
for(offspring in unique(all_families_imendel$offspring_id)){
  all_families_imendel_summ_pat <-
    all_families_imendel_pit_tag %>%
    filter(offspring_id==offspring) %>%
    group_by(paternal_pit_tag, offspring_id, offspring_pit_tag) %>%
    summarise(mean = mean(num_mendelian_errors),
              sd = sd(num_mendelian_errors))  %>%
    rbind(all_families_imendel_summ_pat)
  
  all_families_imendel_summ_mat <-
    all_families_imendel_pit_tag %>%
    filter(offspring_id==offspring) %>%
    group_by(maternal_pit_tag, offspring_id, offspring_pit_tag) %>%
    summarise(mean = mean(num_mendelian_errors),
              sd = sd(num_mendelian_errors))  %>%
    rbind(all_families_imendel_summ_mat)
}

#Get the value and pit tag of the minimum average mendelian errors for males
all_families_imendel_summ_pat_mins <-
  all_families_imendel_summ_pat %>%
  group_by(offspring_pit_tag) %>%
  summarize(min=min(mean)) %>%
  merge(select(all_families_imendel_summ_pat, paternal_pit_tag, mean, sd),
        by.x = "min",
        by.y = "mean") %>%
  select(offspring_id, offspring_pit_tag, paternal_pit_tag, min, sd)

colnames(all_families_imendel_summ_pat_mins) <- c("offspring_id", "offspring_pit_tag", "paternal_pit_tag", "min_paternal_me_mean", "min_paternal_me_sd")

#Get the value and pit tag of the minimum average mendelian errors for females
all_families_imendel_summ_mat_mins <-
  all_families_imendel_summ_mat %>%
  group_by(offspring_pit_tag) %>%
  summarize(min=min(mean)) %>%
  merge(select(all_families_imendel_summ_mat, maternal_pit_tag, mean, sd),
        by.x = "min",
        by.y = "mean") %>%
  select(offspring_id, offspring_pit_tag, maternal_pit_tag, min, sd)

colnames(all_families_imendel_summ_mat_mins) <- c("offspring_id", "offspring_pit_tag", "maternal_pit_tag", "min_maternal_me_mean", "min_maternal_me_sd")

#Get the value and pit tag of the minimum average mendelian errors for parentage pairs
all_families_imendel_summ_parentage_pair_mins <-
  all_families_imendel_pit_tag %>%
  group_by(offspring_pit_tag) %>%
  summarize(min=min(num_mendelian_errors)) %>%
  merge(select(all_families_imendel_pit_tag, parentage_pair_pit_tag, num_mendelian_errors),
        by.x = "min",
        by.y = "num_mendelian_errors") %>%
  select(offspring_pit_tag, parentage_pair_pit_tag, min)

colnames(all_families_imendel_summ_parentage_pair_mins) <- c("offspring_pit_tag", "parentage_pair_pit_tag", "min_parentage_pair_me")

#Merge all three summary dataframes
all_families_imendel_summ_mins <-
merge(all_families_imendel_summ_pat_mins,
      all_families_imendel_summ_mat_mins)

clipr::write_clip(all_families_imendel_summ_mins)




