# CURRENT PROBLEM - TRUNCATION OF TIME POINTS NOT WORKING IN SIDE THE ASSIGN DATA FUNCTION, SO TVEC STARTS AT T = 14 BAAAAD


workDir<-"/home/graeme/Dropbox/Thymic Development/BusulfanProjectFiles/"
#workDir<-"C:/Users/Graeme/Dropbox/Thymic Development/BusulfanProjectFiles"
setwd(workDir)
E<-exp(1)
source("BusulfanPackages.R")


##########################################################################################################################################################################################################

#Code for bootstrapping naive busulfan data
#User must ensure the following is correct
cellA<-"gdT44lo"
cellB<-"gdT44hi"
#precursor<-paste("sp",substr(cell,1,1),sep="")
precursor<-"dn4"
qZero<-1

#requires folder ROutputs/ exists
##########################################################################################################################################################################################################

#This code will take in the outputs (except ki67 stuff) of DataPrepare.R and output the best fit and some number of bootstrap replicates.

#The bootstrapping is setup as follows:
#First the precursor (SP4/SP8) is fitted to extract nu. This is then used as an input into the peripheral fits where the donor fraction and total counts are fitted simultaneously using 
#the GenSA to minimize the sum of squares.
#Currently the counts are log transformed and the donor fraction is transformed using ArcSinSqrt.
#This is all editable below.

##########################################################################################################################################################################################################

#Find number of cores, leave one spare
nLogicalCores<-detectCores(all.tests = FALSE, logical = TRUE) -1
#setwd("C:/Users/Graeme/Dropbox/Thymic Development/Current/Bootstrap and analysis code/SelfContained")

########################################################################################################################
# In this case we actually need TWO peripheral cell populations at the same time, so we just define a funciton to get 
# the data in a more modular way. This should probably be included into the other codes to save re-using loads of this 
# same code shit a million times, but meh.
get.cell.data<-function(cell,precursor){
  #import data
  dataFolder<-"PreparedDataR/"
  countData <- read.table(paste(dataFolder,cell,"TotalCounts.txt",sep=""), sep="\t")
  ratioData <-read.table(paste(dataFolder,cell,"SPRatio.txt",sep=""), sep="\t")
  precursorData <-read.table(paste(dataFolder,precursor,"TotalCounts.txt",sep=""), sep="\t")
  
  #get rid of all data prior to reconstitution and shift everything so that we start at t = 0
  reconstitution.time<-42 #days
  apply.reconstitution.time<-function(set){
    temp<-set[set$V1>=reconstitution.time,]
    temp[,1]<-temp[,1]-reconstitution.time

    temp # This code differs from other apply.reconstitution() definitions because it is INSIDE A FUNCTION and therefore i had trouble with assign() .enviroments, so structure here different
  }
  
  countData<-apply.reconstitution.time(countData)
  ratioData<-apply.reconstitution.time(ratioData)
  precursorData<-apply.reconstitution.time(precursorData)
  
  #define timepoint vectors
  tVec<-countData[,1]
  precursor.tVec<-precursorData[,1]
  
  #scale N0 by 10^7 to make fitting easier
  countData[,2] <-countData[,2]/10^7
  peripheral.data<-data.frame(t=tVec,counts=countData[,2],ratio=ratioData[,2])
  return(list(precursorData,peripheral.data))
}
# Assuming for now they have the same precursor - basically just what data we use to define the level of chimerism chi,
# Might affect peripheral plots - just be aware of what you are plotting! (i.e. [B/chi]/[A/chi] = B/A only if chi's are same)
cellA.data<-get.cell.data(cellA,precursor)
cellB.data<-get.cell.data(cellB,precursor)
tVec.A<-cellA.data[[2]]["t"][,1]
tVec.B<-cellB.data[[2]]["t"][,1]
tVec.merged<-unique(c(tVec.A,tVec.B))

index.times.A<-as.data.frame(table(tVec.A))["Freq"][,1]
index.times.B<-as.data.frame(table(tVec.B))["Freq"][,1]


indexes.A<-foreach(i = 1:length(index.times.A),.combine='c') %do% {
  rep(i,index.times.A[i])
}
indexes.B<-foreach(i = 1:length(index.times.B),.combine='c') %do% {
  rep(i,index.times.B[i])
}

########################################################################################################################
#define transforms
Tform1 <- function(x) {
  log(x)
}

Tform2 <- function(x) {
  asin(sqrt(x))
}

transform.frames<-function(frame){
  # assumes counts, ratios (no time) structure
  data.frame(Tcounts=Tform1(frame["counts"]),Tratio=Tform2(frame["ratio"]))
}
ln<-function(x){
  log(x, base = exp(1))
}

ssr.from.frame<-function(frame){
  ssrC<-sum(frame["counts"]^2)
  ssrR<-sum(frame["ratio"]^2)

  return(nrow(frame)*ln(ssrC*ssrR))
}

source("gdTfunctions/CppODE.R", local = TRUE)
################################################
SSRfn<-function(inVec,dVec,nuVal){
  
  inVec<-as.data.frame(t(as.numeric(inVec)))
  names(inVec)<-header
  
  inVec<-cbind(inVec,chi=0.5,nu=nuVal)
  residuals<-residFn(inVec,dVec)

  ssrA<-ssr.from.frame(residuals[[1]])
  ssrB<-ssr.from.frame(residuals[[2]])
  SSR<-ssrA+ssrB # here is a sum because we have logs, see definition of 'ssr.from.frame()'.   

  return(ssrA+ssrB)
}
################################################
plot.functions<-function(inVec,nuVal){
  plot.times<-seq(1:max(tVec.merged))
  inVec<-as.data.frame(t(as.numeric(inVec)))
  names(inVec)<-header
  
  inVec<-cbind(inVec,chi=0.5,nu=nuVal)
  IC<-NULL
  for(i in 1:length(odeICs)){
    IC[[i]]<-with(as.list(inVec),eval(odeICs[[i]]))
  }
  
  sol<-NDSolve(plot.times,inVec,IC)
  sol.A<-sol[[1]]
  sol.B<-sol[[2]]
  return(list(sol.A,sol.B))
}
################################################
residFn<-function(inVec,dVec){
  # all of this assumes good ordering of time points in merged/cell specific time vectors and data
  IC<-NULL
  for(i in 1:length(odeICs)){
    IC[[i]]<-with(as.list(inVec),eval(odeICs[[i]]))
  }
  sol<-NDSolve(tVec.merged,inVec,IC)
  sol.A<-sol[[1]][sol[[1]]$t %in% tVec.A,]
  sol.B<-sol[[2]][sol[[2]]$t %in% tVec.B,]
  
  # calculate counts and ratio on transformed scale
  finalA.vals<-transform.frames(sol.A[indexes.A,][,-1])
  finalB.vals<-transform.frames(sol.B[indexes.B,][,-1])
  
  # function - data
  residA<-finalA.vals-dVec[[1]]
  residB<-finalB.vals-dVec[[2]]

  return(list(residA,residB))
}

################################################

findFit<-function(dVec,nuVal){
  out <- GenSA(par=start,lower = lower, upper = upper, fn = SSRfn, control=list(verbose=TRUE,temperature=temp,maxit=maxit),dVec,nuVal)
  fitPars<<-out["par"][[1]]
  return(fitPars)
  
}
temp<-10^5
maxit<-2*10^3
################################################################################################
################################################################################################
# cell A is defined as being thymically sourced, so we use it's precursor to estimate nu
precursor.data<-cellA.data[[1]][,2]
precursor.times<-cellA.data[[1]][,1]
precursorFit<-nls(Tform1(precursor.data)~logA-precursor.times*nu, start=list(logA=10,nu=0.1))
thymusBF<-predict(precursorFit) #defined on transformed scale!
precursorResids<-thymusBF-Tform1(precursor.data)
nuVal<-coef(precursorFit)[[2]]

# transform data
data.vector<-list(
  transform.frames(cellA.data[[2]]),
  transform.frames(cellB.data[[2]]))

header<-c("lambdaA","lambdaB","TonA","beta","A0","B0","muA","muB")#fit lambda separately, hence why dimension is only 9 (set chi to fixed value)
dimension <- length(header)
#c("lambdaA","lambdaB","TonA","beta","A0","B0","muA","muB")
start<-rep(0.1, dimension)
lower <-c(-1,-1,10^-5,10^-5,10^-5,10^-5,10^-5,10^-5)
upper <- c(1,1,1,5,5,5,1,1)

bestFit<-findFit(data.vector,nuVal)
names(bestFit)<-header
#expHeader<-c("lnL/2","nu",header,"temp","maxit")
#data.vector<-c(SSRfn(bestFit,data.vector,nuVal),nuVal,bestFit,temp,maxit)
#names(data.vector)<-expHeader

#write.xlsx(bestFit,paste(cell,"bestFit.xlsx",sep=""))

################################################################################################
predicted<-plot.functions(bestFit,nuVal)
cellACountsData<-cellA.data[[2]][c("t","counts")]
cellACountsFn<-predicted[[1]][c("t","counts")]

#cellBCountsData<-cellB.data[[2]][c("t","counts")]
#cellBCountsFn<-predicted[[2]][c("t","counts")]

cellARatioData<-cellA.data[[2]][c("t","ratio")]
cellARatioFn<-predicted[[1]][c("t","ratio")]

#cellBRatioData<-cellB.data[[2]][c("t","ratio")]
#cellBRatioFn<-predicted[[2]][c("t","ratio")]

ggplot(cellACounts,aes(x=t,y=counts))+geom_point()+geom_line(data=cellACountsFn)
ggplot(cellBCountsData,aes(x=t,y=counts))+geom_point()+geom_line(data=cellBCountsFn)


ggplot(cellARatioData,aes(x=t,y=ratio))+geom_point()+geom_line(data=cellARatioFn)
ggplot(cellBRatioData,aes(x=t,y=ratio))+geom_point()+geom_line(data=cellBRatioFn)
