#Load the required libraries
library(factoextra)
library(clv)
library(ggplot2)
library(pheatmap)

#Import the dataset
codon_usage <- read.csv("C:/Users/DARKO/Desktop/Projects/Codon usage taxa project/codon_usage.csv")
View(codon_usage)

#Make the data object adenodat which will be used for data wrangling
adenodat=codon_usage

#Create advir dataset containing only virus columns
advir <- adenodat[grep(c("vrl"), adenodat$Kingdom), ]

#use a matchpattern object containing "Human" and "adenovirus" strings to extract 
#Human adenovirus data
matchpattern <- c("Human", "adenovirus")

#Use another grepl to select rows containing "adenovirus" 
advir<- advir[grepl("Human adenovirus", advir$SpeciesName),]

#Make Rownames => Species
rownames(advir)<-advir$SpeciesName

#Finalized adenovirus dataframe -> check
str(advir)

#Using the silhouette score to evaluate k, keep in mind this is statistically optimal
#but biologically theoretical part must also be evaluated
fviz_nbclust(advir[,6:69], hcut, method = "silhouette", k.max = 25) +
  labs(subtitle = "Find optimal number of clusters")

#This doesnt mean that other k cound not be better due to Separation Height Variability

  
#Creating the distance matric of all codon usage values
dmatrix=dist(advir[,6:69], method = "euclidean")
  
#Will be using the complete linkage to perform aglomerative clustering
ahc=hclust(dmatrix, method = "complete")



#This will  create the algomerative hierarchical clustering dendrogram for k=4
fviz_dend(ahc, cex = 0.7, k = 4,
          k_colors = c("blue", "green3", "red", "aquamarine"),
          rect = TRUE, lower_rect = -0.1, rect_lty = 4 )


#To perform the inter/intra-cluster distance analysis the hierarchy should be cut into k parts
#According to the second elbow the dendrogram will be cut into 4 parts
cutdend<-cutree(ahc, k=4)


#To visualizae the heatmap and both codon based and species based dendrograms 
#make sure data is numeric and its log transformed for better clustering and visualization
data<- as.data.frame(sapply(advir[,6:69],as.numeric))
data<-log(data+1)

#rownames should correspond to the species names
rownames(data)<-advir$SpeciesName

#Using the pheatmap package a publication ready visualization is created
pheatmap(data, cutree_rows = 4)


#cll.scatt.data() is used to evaluate and output inter/intra-cluster distances
#and cls.attrib() to evaluate cluster size
#ggplot will be used to plot the intercluster distance medians and IQRs

hid=cls.scatt.data(data,cutdend, dist = "euclidean")
dfinter=as.data.frame(hid$intercls.complete)
intradim<-hid$intracls.complete
att<- cls.attrib(data,cutdend)
clsize<-att$cluster.size
dfinter=as.data.frame(hid$intercls.complete)
ggplot(stack(dfinter), aes(x = ind,y = values)) +
  geom_boxplot()




#This will  create the alliterative hierarchical clustering dendrogram
#but this time using k as 7, much better biologically speaking
fviz_dend(ahc, cex = 0.7, k = 7,
          k_colors = c("blue", "green3", "red", "aquamarine","orange","purple","gray"),
          rect = TRUE, lower_rect = -0.1, rect_lty = 4 )


#To perform the inter/intra-cluster distance analysis the hierarchy should be cut into k parts
cutdend2<-cutree(ahc, k=7)




#Repeat the same log transformation process for second dataset meant for k=7
data2<- as.data.frame(sapply(advir[,6:69],as.numeric))
data2<-log(data+1)

#rownames again assigned as in advir
rownames(data2)<-advir$SpeciesName

#Using the pheatmap() a publication ready visualization is created for k=7 structure
pheatmap(data, cutree_rows = 7)


#Using this code, hierarchical clustering intracluster dimensions,
#intercluster distances will be calculated and plotted for the new design

hid2=cls.scatt.data(data2,cutdend2, dist = "euclidean")
dfinter2=as.data.frame(hid2$intercls.complete)
intradim2<-hid2$intracls.complete
att2<- cls.attrib(data2,cutdend2)
clsize2<-att2$cluster.size
dfinter2=as.data.frame(hid2$intercls.complete)
ggplot(stack(dfinter2), aes(x = ind,y = values)) +
  geom_boxplot()


#Visualize the dendrogram with k=6
fviz_dend(ahc, cex = 0.7, k = 6,
          k_colors = c("blue", "green3", "red", "aquamarine","orange","purple","gray"),
          rect = TRUE, lower_rect = -0.1, rect_lty = 4 )



#Cut the dendrogram for the hierarchy  k =6 parts
cutdend3<-cutree(ahc, k=6)


#log transformation process for second dataset meant for k=6
data3<- as.data.frame(sapply(advir[,6:69],as.numeric))
data3<-log(data+1)

#rownames again assigned as in advir k=6
rownames(data)<-advir$SpeciesName

#Pheatmap() for a publication ready visualization is created for k=6 structure
pheatmap(data, cutree_rows = 6)


#Last setting for k=6 to check the silhouette differences between k=7 and k-6
#Same steps as in previous two setting except for k=6
hid3=cls.scatt.data(data3,cutdend3, dist = "euclidean")
dfinter3=as.data.frame(hid3$intercls.complete)
intradim3<-hid3$intracls.complete
att3<- cls.attrib(data3,cutdend3)
clsize3<-att3$cluster.size
dfinter3=as.data.frame(hid3$intercls.complete)
ggplot(stack(dfinter3), aes(x = ind,y = values)) +
  geom_boxplot()




