#setwd("C:\\Users\\Anjana Thimmaiah\\Documents\\Research 2017\\paper1&2_data\\paper1")
getwd()
setwd(getwd())
#-----------------------------------------------------------------------------------------

#install required packages

packages <- c("bnlearn","prospectr","pracma","bclust", "clusterSim","WriteXLS","dplyr","ChemoSpec")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}
#load the libararies for the session
library(MASS)
library(Rgraphviz)
library(graph)
library(BiocGenerics)
library(parallel)
library(prospectr)
library(pracma)
library(dplyr)
library(base)
library(grid)
library(WriteXLS);
library(clusterSim);
library(bclust);
library(prospectr)
library(bnlearn)
library(seewave)
library(stats)
library(nlme)
library(gclus)
library(HSAUR)
library(tools)
library(Rgraphviz)

#load raw data to the work space

#XXX<-read.csv("raw_low_highanddata_Oc21no380.csv",header=T,row.names=1) #bclustresult_Oct26raw_alldat.rda to be used with this data
XXX<-read.csv("H1optb_Bio.csv")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

#BiocManager::install(("Rgraphviz"),header=F,row.names=NULL) 

XXX1<-XXX[1:1764,1:5]#choosing your rows and columns
#XXX1<-XXX[1:400,1:42]


dim(XXX1)


i=ncol(XXX1)# calculate the number of columns of XXX1

for(j in 1:i)
  #using sprintf or paste, create string:for automatic renaming of clusters
  
  names(XXX1)<-sprintf("sample %d",1:i)

#------------------------------------------------------------------------------------------------------

# Bayesian Learning strats from here, if you have your components from SMCR then
# that data will be you 'Xfinal'

# for this case we had to prepare the data from the forginig section to 'Xfinal' 

#prepare the new data for BN learning 

#nn=ncol(m1)
#n=nn-1
#xcust=t(aggregate(m1[, 1:n], list(m1$groups), rms))#mean
#Xfinal=as.data.frame(xcust[2:nn,])
Xfinal=as.data.frame(XXX[1:1764,1:5])

#change names of the variables

j=ncol(Xfinal)# calculate the number of columns/clusters
for(i in 1:j)
  #using sprintf or paste, create string:for automatic renaming of clusters
  names(Xfinal)<-sprintf("group %d",1:j)

#View(Xfinal)




#==========================================================================================================
#assesing the normality of the new dataset


#hist(Xfinal$group5, prob=TRUE)            # prob=TRUE for probabilities not counts

#lines(density(Xfinal$group2))             # add a density estimate with defaults
#lines(density(Xfinal$group5, adjust=2), lty="dotted")
# a

#============================================================================================================
#correllation plot of the new variables 


my.abs=abs(cor(Xfinal))
my.colors=dmat.color(my.abs)
my.ordered=order.single(cor(Xfinal))
cpairs(Xfinal, my.ordered,panel.colors=my.colors,gap=.5)

#============================================================================================================
#

#
#                                 3. PERFORM BAYESIAN STRUCTURE LEARNING                                       #
#-------------------------------------------------------------------------------------------------------------##



#score beased, hillclibing algorithm(greed sreach)
#

bh = hc(Xfinal,optimized = F,score="bic-g")

#Plot the network and format its arcs according to the strength of the dependencies they represent.
#
hlight <- list(nodes = nodes(bh),arcs = arcs(bh),col = "blue", textCol = "blue")#,lwd =2

strength1= arc.strength(bh, Xfinal,criterion = "bic-g")
strength1


strength.plot(bh, strength1, shape = "circle", sub = " Hill-Climbing",highlight = hlight) 
#possible layouts; dots, neato, twopi, circo and fdp
score(bh, data = Xfinal, type = "bic-g")
#----------------------------------------------------------------------------------------------------------
#tabu

mtabu=tabu(Xfinal,optimized = F,score="bic-g")

#graphviz.plot(mtabu)

strengthq2= arc.strength(mtabu, Xfinal,criterion = "bic-g")

strengthq2

score(mtabu, data = Xfinal, type = "bic-g")

hlight <- list(nodes = nodes(mtabu), arcs = arcs(mtabu),col = "green", textCol = "green")#,lwd =2

strength.plot(mtabu, strengthq2, shape = "circle", sub = "Tabu Search",highlight = hlight )



#---------------------------------------------------------------------------

#hybrid methods

hybr_mmpc=mmhc(Xfinal)

strength3 = arc.strength(hybr_mmpc, Xfinal)
strength3

score(hybr_mmpc, data = Xfinal, type = "bic-g")

hlight <- list(nodes = nodes(hybr_mmpc), arcs = arcs(hybr_mmpc),col = "red", textCol = "red")#,lwd =2



strength.plot(hybr_mmpc, strength3, shape = "circle", sub = "Max-Min Hill-Climbing",highlight = hlight )
#--------------------------------------------------------------------------------------------------
#model selection using cross validation and score
bn.cv(Xfinal, 'hc', k=3, loss = "logl-g")#, loss.args = list(target = "F"))
bn.cv(Xfinal, 'tabu', k=3, loss = "logl-g")#, loss.args = list(target = "F"))
bn.cv(Xfinal, 'mmhc', k=3, loss = "logl-g")#, loss.args = list(target = "F"))
#=================================================================================================


