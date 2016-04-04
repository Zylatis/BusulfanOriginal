#directory of file
#workDir<-"/home/graeme/Dropbox/Thymic Development/BusulfanProjectFiles/"
workDir<-"C:/Users/Graeme/Dropbox/Thymic Development/BusulfanProjectFiles"
setwd(workDir)
dataFolder<-"PreparedDataR/"
source("BusulfanPackages.R")

#working with q = 0 here. See BusulfanNaiveFit.R for further discussion about q != 0 stuff.
q<-0

##########################################################################################################################################################################################################
##########################################################################################################################################################################################################

#Code for analysing output of bootstrap and best fit data (from BusulfanNaiveFit.R)
#Currently only tabular data is printed to the screen, but can easily be plotted/exported if required.
#Only works for one cell at a time, but could be looped and paired with export functions.

##########################################################################################################################################################################################################

#choose cell type
cell<-"4nai"

#import experiment data
countData <- read.table(paste(dataFolder,cell,"TotalCounts.txt",sep=""), sep="\t")
ratioData <-read.table(paste(dataFolder,cell,"SPRatio.txt",sep=""), sep="\t")
precursorData <-read.table(paste(dataFolder,"sp",substr(cell,1,1),"TotalCounts.txt",sep=""), sep="\t")

t.Shift<-42

#shift back to days post BMT rather than post-reconstitution
countData[,1]<-countData[,1]#+t.Shift
ratioData[,1]<-ratioData[,1]#+t.Shift
precursorData[,1]<-precursorData[,1]#+t.Shift

data.names<-c("t","fn")
names(countData)<-data.names
names(ratioData)<-data.names
names(precursorData)<-data.names

#define timepoint vectors
tVec<-countData[,1]

#get best and bootfit data
bestFitData<- read.xlsx(paste("ROutputs/RCD",cell,"best.xlsx",sep=""),sheetName="Results")
bootFitData<- read.xlsx(paste("ROutputs/RCD",cell,"Boot.xlsx",sep=""),sheetName="Results")

if(ncol(bestFitData)!=ncol(bootFitData)||(sum(names(bootFitData)==names(bestFitData))!=ncol(bestFitData))){
  print("Problem! Headers in boot and best fit files don't match, will cause headaches!")
}
data.header<-names(bootFitData)
#restore 10^7 factor (can do in previous boot outputs if easier, or remove for ease of plotting etc)
bestFitData[,3]<-bestFitData[,3]*10^7
bootFitData[,3]<-bootFitData[,3]*10^7

#make header of parameters of interest
header<-c("log(2)/nu","lambda","N0","theta","I0","NT=I0+N0","I0/NT","theta/NT")
#define useful parameter functions
theta<-function(data){
  data[,"TonN"]*data[,"N0"] #TonN = Theta on N
}
I0<-function(data){
  data[,"alpha"]*data[,"N0"] #I_0 = alpha*N_0
}
#define peripheral ratio and total count functions
ratioFn<-function(t,pars){
  x1<-with(pars,1/(1 - ((-1 - alpha/exp(lambda*(-1 + q)*t) + mu)*(lambda - nu))/(mu*(lambda - nu) + (-1 + exp((lambda - nu)*t))*TonN)))
  return(x1)
}

totalCountFn<-function(t,pars){
  with(pars,(N0*(alpha*exp(lambda*t)*(lambda - nu) + exp(lambda*q*t)*(lambda - nu - TonN) + exp((lambda - nu + lambda*q)*t)*TonN))/(exp(lambda*(1 + q)*t)*(lambda - nu)))
}

#define function to take fitted parameters and turn into 'physically interesting' parameters (see 'header')
makePhysicalBits<-function(data){
  out<-as.data.frame(cbind(
    log(2)/data[,"nu"],                                 #thymic involution halflife
    data[,"lambda"],                                    # net loss rate
    data["N0"],                                         # initial displacable numbers (donor + disp. host)
    theta(data),                                        # thymic output numbers
    I0(data),                                           # incumbent numbers
    data[,"N0"]+I0(data),                               # NT_0=I_0+N_0
    I0(data)/( data[,"N0"]+I0(data)),                   # initial incumbent fraction
    theta(data)/( data[,"N0"]+I0(data))))               # initial thymic output fraction of total pool
  names(out)<-header
  return(out)
}

#function to calculate quantiles of columns in a dataframe
calcQuantiles<-function(data){  
  lower<-list()
  upper<-list()
 for(i in 1:length(data[1,])){
   lower[i]<-quantile(data[,i],0.025)
   upper[i]<-quantile(data[,i],0.975)
 }
 out<-as.data.frame(cbind(lower,upper))
 return(out)
}

#PHYSICAL PARAMETER DATA TABLE
##########################################################################################################################################################################################################

#make joint dataframe showing best fit and CI values for these physical parameters
bestFitDataPhysical<-makePhysicalBits(bestFitData)
bootFitDataPhysical<-makePhysicalBits(bootFitData)
names(bestFitDataPhysical)<-header

bootFitDataPhysical<-as.data.frame(t(calcQuantiles(bootFitDataPhysical)))
names(bootFitDataPhysical)<-header

physicalPars<-rbind(bestFitDataPhysical,bootFitDataPhysical)
rownames(physicalPars)<-c("bestFit","lower","upper")
physicalPars<-t(physicalPars)
#RATIO AND TOTAL COUNT PLOTTING DATA WITH CIs
#now also with incumbent frac 22 oct 15
##########################################################################################################################################################################################################

makePlotCI<-function(fun){
  lower<-list()
  upper<-list()
  val<-list()
  tVals<-list()
  for( t in t.Shift:400){
    tVals[t]<-t
    val[t]<-fun(t-t.Shift,bestFitData)
    bootFnVal<-fun(t-t.Shift,bootFitData)
    lower[t]<-quantile(bootFnVal,0.025)
    upper[t]<-quantile(bootFnVal,0.975)
  }
  out<-as.data.frame(cbind(unlist(tVals),unlist(val),unlist(lower),unlist(upper)))
  names(out)<-c("t","fn","lowerCI","upperCI")
  return(out)
}

incumbentFrac<-function(t,pars){
  with(pars,alpha*N0/totalCountFn(t,pars))
}

#need plotting code
countBits<-makePlotCI(totalCountFn)
ratioBits<-makePlotCI(ratioFn)
incumbentBits<-makePlotCI(incumbentFrac)


pointSize<-3
colN="red"
countPlot<-ggplot(countBits[c("t","fn")],aes(t,fn))+geom_line(aes(colour=colN))+geom_ribbon(aes(ymin=countBits$lowerCI,ymax=countBits$upperCI,colour=colN,fill=colN),alpha=0.2)+ylab("Total Counts (10^-7)")  + xlab("t (days)")+ ylim(0,5*10^7)+geom_point(data=countData,size=pointSize)+theme(legend.position="none")+ggtitle(paste(cell, " Busulfan fits to total counts and donor fraction in LN+SP",sep="")) + theme(plot.title = element_text(lineheight=.8, face="bold"))
ratioPlot<-ggplot(ratioBits[c("t","fn")],aes(t,fn))+geom_line(aes(colour=colN))+geom_ribbon(aes(ymin=ratioBits$lowerCI,ymax=ratioBits$upperCI,colour=colN,fill=colN),alpha=0.2)+ylab("Rescaled donor frac")  + xlab("t (days)")+ ylim(0,1)+geom_point(data=ratioData,size=pointSize)+theme(legend.position="none")#+ggtitle(paste(cell, " Busulfan fits to total counts and donor fraction in LN+SP",sep="")) + theme(plot.title = element_text(lineheight=.8, face="bold"))
incumbentFracPlot<-ggplot(incumbentBits[c("t","fn")],aes(t,fn))+geom_line(aes(colour=colN))+geom_ribbon(aes(ymin=incumbentBits$lowerCI,ymax=incumbentBits$upperCI,colour=colN,fill=colN),alpha=0.2)+ylab("Incumbent frac")  + xlab("t (days)")+ ylim(0,0.25)+theme(legend.position="none")#+ggtitle(paste(cell, " Busulfan fits to total counts and donor fraction in LN+SP",sep="")) + theme(plot.title = element_text(lineheight=.8, face="bold"))


grid.arrange(countPlot, ratioPlot,incumbentFracPlot,ncol = 1)


#LIFETIME AND INTERDIVISION TIME ESTIMATION
##########################################################################################################################################################################################################

#rho function as written in paper appendix
#rhoFn<-function(f,sig,p,Tval){
#  (epsilon*p*Tval - f*(1 + f*(-1 + sig) + p*(epsilon + sig - epsilon*sig)*Tval))/((-1 + f)*(2 + f*(-1 + sig))*Tval)
#}
rhoFn<-function(kappa,sig,lambda,Tval){
  ( kappa*(1 + kappa*(-1 + sig) + lambda*(epsilon + sig - epsilon*sig)*Tval)-epsilon*lambda*Tval)/((1 - kappa)*(2 + kappa*(-1 + sig))*Tval)
  }

#20% of thymically derived cells are ki67+
epsilon<-0.2
#vector of sigma values to loop over
sigT<-c(0.1,1,10)
#import ki67 data (taken to be >30wks in these files, see DataPrepare.nb)
ki67Donor <- read.table(paste(dataFolder,cell,"ki67donor.txt",sep=""), sep="\t")[,2]/100
ki67Host <- read.table(paste(dataFolder,cell,"ki67host.txt",sep=""), sep="\t")[,2]/100
ki67Full<-c(unlist(rbind(ki67Donor,ki67Host)))
ki67Host<-unlist(ki67Host)

#function to calculate pooled (disp+inc) interdivision and life-times
#relies on global definition of invIncBoot, invRhoBoot, and invDeltaBoot which is not a great idea!
computePoolFrac<-function(time){
  bootSet<-list()
  
  for(i in 1:length(bootFitData[,2])){
    bootSet[i]<-with(as.list(bootFitData[i,]),incumbentFrac(time,bootFitData[i,]))
  }
  bootSet<-unlist(bootSet)
  
  pooledDelta<-list()
  pooledRho<-list()
  for(i in 1:length(invDeltaBoot)){
    incVal<-invIncBoot[i]
    deltaVal<-invDeltaBoot[i]
    rhoVal<-invRhoBoot[i]
    incFrac<-sample(bootSet,1)
    pooledDelta[i]<-incFrac*incVal+(1-incFrac)*deltaVal
    pooledRho[i]<-incFrac*incVal+(1-incFrac)*rhoVal
  }
  pooledDelta<-unlist(pooledDelta)
  pooledRho<-unlist(pooledRho)
  return<-cbind(pooledDelta,pooledRho)
}
#round values
format<-function(bootSet){
round(c(mean(bootSet),quantile(bootSet,0.025),quantile(bootSet,0.975)))
}
#make lists to dump stuff into. no doubt a way better way to do this!
lifetime<-list()
interdiv<-list()
incumtime<-list()

deltaSet<-list()

pooledT1Delta<-list()
pooledT2Delta<-list()
pooledT1Rho<-list()
pooledT2Rho<-list()
set.seed(0)

monteCarloCount<-10^4
Tlist<-log(rlnorm(monteCarloCount, 4, 0.5))
#loop over sigma values
for(j in 1:length(sigT)){
sig<-sigT[j]

  invRhoBoot<-list()
  invDeltaBoot<-list()
  invIncBoot<-list()
  deltaBoot<-list()
  set.seed(0)
   for(i in 1:monteCarloCount) {
     ki67Choice<-sample(ki67Full,1)
     ki67ChoiceI<-sample(ki67Host,1)
     
     lambdaChoice<-sample(bootFitData[,2],1)
     Tchoice<-sample(Tlist,1)
    
    
     rho<-rhoFn(ki67Choice,sig,lambdaChoice,Tchoice)
     delta<-lambdaChoice+rho
     deltaBoot[i]<-delta
     rhoI<-rhoFn(ki67ChoiceI,sig,0,Tchoice)
       invRhoBoot[i]<-1/rho
       invDeltaBoot[i]<-1/delta
       invIncBoot[i]<-1/rhoI
   }
  
  
  invRhoBoot<-unlist(invRhoBoot)
  invDeltaBoot<-unlist(invDeltaBoot)
  invIncBoot<-unlist(invIncBoot)
  
  joined<-as.data.frame(cbind(invRhoBoot,invDeltaBoot,invIncBoot))
  joined<-joined[(joined$invRhoBoot>0)&(joined$invIncBoot>0),]
  
  interdiv<-rbind(interdiv,format(joined[,1]))
  lifetime<-rbind(lifetime,format(joined[,2]))
  incumtime<-rbind(incumtime,format(joined[,3]))

  #compute pooled population stuff at 14 and 44 wks (post-BMT!)
  reconsTime<-6*7
  t1<-14*7-reconsTime
  t2<-44*7-reconsTime
  
  t1Pooled<-computePoolFrac(t1)
  t2Pooled<-computePoolFrac(t2)

  pooledT1Delta<-rbind(pooledT1Delta,format(t1Pooled[,1]))
  pooledT2Delta<-rbind(pooledT2Delta,format(t2Pooled[,1]))

  pooledT1Rho<-rbind(pooledT1Rho,format(t1Pooled[,2]))
  pooledT2Rho<-rbind(pooledT2Rho,format(t2Pooled[,2]))
}

r.names<-c("sig=0.1","1","10")
c.names<-c("bestFit","lowerCI","upperCI")

 kinetic.Parameters<-list(lifetime,interdiv,incumtime,pooledT1Delta,pooledT2Delta,pooledT1Rho,pooledT2Rho)
 
 for(i in 1:length(kinetic.Parameters)){
   rownames(kinetic.Parameters[[i]])<-r.names
   colnames(kinetic.Parameters[[i]])<-c.names
 }


#THYMIC AND PERIPHERAL DIVISION PLOTTING DATA WITH CIs
##########################################################################################################################################################################################################

E<-exp(1)
thymus<-function(t,pars){
  internal<-pars
  names(pars)<-data.header
  output<-with(pars,TonN*N0*E^(-t*nu))
  return(output)
}

thymicProduction<-makePlotCI(thymus)
thymusPlot<-ggplot(thymicProduction[c("t","fn")],aes(t,fn))+geom_line(aes(colour=colN))+geom_ribbon(aes(ymin=thymicProduction$lowerCI,ymax=thymicProduction$upperCI,colour=colN,fill=colN),alpha=0.2)+
  ylab("Daily export from thymus")  +
  xlab("t (days)")+ 
  ylim(0,600000)+theme(legend.position="none")#+ggtitle(paste(cell, " Busulfan fits to total counts and donor fraction in LN+SP",sep="")) + theme(plot.title = element_text(lineheight=.8, face="bold"))
