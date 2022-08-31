#R code for analysis of Gene expression profile of CD22 CAR T-cell products from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE200296
#See PDF for citation and figure output

#***Initial data loading and cleaning***

#read in data (removed metadata from txt file)
dat <- read.csv("data.csv", row.name=1)

#how many negative values and how many rows with negative values
length(dat[dat < 0])
#441
sum(rowSums(dat<0) > 0)
#170

#I don't love deleting 170 genes but have to do something with them for log transform
#Will start with deleting and can try something else later if needed

#delete rows with negatives
dat[dat < 0] <- NA
dat <- na.omit(dat)

#read in annotations (manually extracted from metadata)
ann <- read.csv("ann.csv", row.name=1, header=F)

#get class from annotations
cl <- as.numeric(ann[2,])


#***Outlier inspection and removal***

#correlation heatmap
library(gplots)
dat.cor <- cor(dat)

layout(matrix(c(1,1,1,1,1,1,1,1,2,2), 5, 2, byrow = TRUE))
par(oma=c(5,7,1,1))
cx <- rev(colorpanel(25,"yellow","black","blue"))
leg <- seq(min(dat.cor,na.rm=T),max(dat.cor,na.rm=T),length=10)
image(dat.cor,main="Correlation heat map of samples in CD22 CAR T-cell product dataset",axes=F,col=cx)
axis(1,at=seq(0,1,length=ncol(dat.cor)),label=dimnames(dat.cor)[[2]],cex.axis=0.9,las=2)
axis(2,at=seq(0,1,length=ncol(dat.cor)),label=dimnames(dat.cor)[[2]],cex.axis=0.9,las=2)

image(as.matrix(leg),col=cx,axes=F)
tmp <- round(leg,2)
axis(1,at=seq(0,1,length=length(leg)),labels=tmp,cex.axis=1)

#cluster dendogram
dat.trans <- t(dat)
dat.dist <- dist(dat.trans,method="euclidean")
dat.clust <- hclust(dat.dist,method="single")
plot(dat.clust,labels=names(dat),cex=0.75)


#average correlation plot
dat.avg <- apply(dat.cor,1,mean)
par(oma=c(3,0.1,0.1,0.1))
plot(c(1,length(dat.avg)),range(dat.avg),type="n",xlab="",ylab="Avg r",main="Avg correlation of CD22 CAR T-cell product samples",axes=F)
points(dat.avg,bg="red",col=1,pch=21,cex=1.25)
axis(1,at=c(1:length(dat.avg)),labels=dimnames(dat)[[2]],las=2,cex.lab=0.4,cex.axis=0.6)
axis(2)
abline(v=seq(0.5,62.5,1),col="grey")


#Decided to drop GSM6030549
#Figure out sample index
which(colnames(dat)=='GSM6030549')


#drop is from dataframe and class labels
dat2 <- subset(dat, select=-GSM6030549)
cl2 <- cl[c(1:13, 15:43)]


#***Differential gene expression analysis***

#log2 the data
dat.log <- log2(dat2)


#To keep this analysis simple, I’m going to make this into a two-class problem and look at differential gene expression 
#I’m going to Group no CRS (class 0) with grade 1 CRS (class 1), and group Grade 2, 3, 4 CRS (class 2,3,4) together. 
#The reason for drawing the line here is that class 1 CRS results in mild symptoms, no effect on cancer treatment, and no need for specific treatment to modulate CRS. 
#At grade 2 CRS is when it start to become a problem clinically causing interruption in treatment and requiring separate intervention. 


#create t-test function to get test statistic for each gene
t.test.all.genes <- function(x,s1,s2) {
       x1 <- x[s1]
x2 <- x[s2]
x1 <- as.numeric(x1)
x2 <- as.numeric(x2)
t.out <- t.test(x1,x2, alternative="two.sided",var.equal=T)
out <- as.numeric(t.out$p.value)
return(out)
} 

#run the function to calc p-values for each gene
pv <- apply(dat.log, 1, t.test.all.genes, s1=(cl2<2), s2=(cl2>1))

#look at genes with a p-value less than 0.05
pv[pv < 0.05]


#calc fold change of those genes

#mean of each group
dat.log.1 <- apply(dat.log[cl2<2],1,mean)
dat.log.2 <- apply(dat.log[cl2>1],1,mean)

fold <- dat.log.1 - dat.log.2

# get fold change where p-value is less than 0.05
fold[names(pv[pv<0.05])]
