switch(Sys.info()[['sysname']],
       Windows= {print("I'm a Windows PC.")
                 setwd("C:/Users/Graeme/Dropbox/Thymic Development/BusulfanProjectFiles/")},
       Linux  = {print("I'm a penguin.")
                 setwd("/home/graeme/Dropbox/Thymic Development/BusulfanProjectFiles/")},
       Darwin = {print("I'm a Mac.")})

########################################################################################################
#CODE FOR SOLVING PDE CASE OF AGE-DEPENDENT LOSS RATES WITHOUT INCUMBENT CELLS
########################################################################################################
args <- commandArgs(trailingOnly = TRUE)
# Change model parameters here - see busulfan PNAS paper for definitions
tau<-42            		# reconstitution time
model<- as.numeric(args[1])     # choice of lambda model
p<- as.numeric(args[2])		# get choice of p (atm looking at 0.5 (not 1/2! seems to mess with 'args'), 1 and 2
treatment.age<-8*7  		# Age at BMT (even though we start t = 0 at treatment time, need treatment age for 'width' of initial distribution of host cells)
N0<-2.78            		# Letting this be free gives crazyness so we fix it based on estimates from WT mice
E<-exp(1)           		# 'E' comes from using FortranForm[] in MM so we define it here
cell<-"4nai"        		# Choose cell type 
TonNMax<- (1 + p)/(treatment.age*p)
########################################################################################################
source("BusulfanPackages.R")
nLogicalCores<-6 #set number of cores to use
########################################################################################################

#define relevant thymic precursor
precursor<-paste("sp",substr(cell,1,1),sep="")
#import data
data.folder<-"PreparedDataR/"
export.folder<-"ROutputs/AgeDependentLossFreeMu/"
countData <- read.table(paste(data.folder,cell,"TotalCounts.txt",sep=""), sep="\t")
ratioData <-read.table(paste(data.folder,cell,"SPRatio.txt",sep=""), sep="\t")
precursor.data <-read.table(paste(data.folder,precursor,"TotalCounts.txt",sep=""), sep="\t")
mu.data <-read.table(paste(data.folder,precursor,"muData.txt",sep=""), sep="\t")
mu.time.vector<-mu.data[,1]
########################################################################################################
# We want to compare this model to that involving incumbents, meaning we need to consider the exact same PERIPHERAL dataset.
# Thus, while we may use pre-reconstitution data in the thymus to estimate the reconstitution rate, we must only consider t>tau data in the periphery.
# (there are other reasons too, basically pre-reconstitution chimerism 'chi' is actually an effective chi = chi(t) so we can't use it in the normal way)

apply.reconstitution.time<-function(set){
  temp<-set[set$V1>=tau,]
  temp[,1]<-temp[,1]
  
  assign(deparse(substitute(set)),temp,envir=.GlobalEnv)
}

apply.reconstitution.time(countData)
apply.reconstitution.time(ratioData)
apply.reconstitution.time(precursor.data)

names(ratioData)<-c("t","ratio")
names(countData)<-c("t","counts")

#define timepoint vectors
tVec<-countData[,1]

#re-scaled counts - helps with numerics a bit
countData[,2] <-countData[,2]/10^7


# Read in the PDE pde.functions as created by the MM file MakePDEFunctions.nb
options(warn=-1) # stops the importer complaining about missing line ends etc
pde.functions<-read.csv(paste("PDEfunctions/functionsModel",toString(model),".csv",sep=""), header = FALSE,quote = "", stringsAsFactors=F)[[1]]
pde.functions<-parse(text=pde.functions)

lambda.function<-read.csv(paste("PDEfunctions/lambdaFunctionModel",toString(model),".csv", sep=""),header = FALSE,quote = "", stringsAsFactors=F)[[1]]
lambda.function<-parse(text=lambda.function)
options(warn=0)

# function to call when we want to evaluate a particular function with a given set of parameters, and at a given 'a' and 't'
eval.fn<-function(a,t,all.pars,fn,extra.fn){
  with(all.pars,eval(fn)*eval(extra.fn))
}

########################################################################################################

# Main function that calculates the rescaled donor fraction and the total counts at a given time (integrates over 'a')
# Currently the pde.functions[[i]] calls are hardcoded to correspond to the order in which the pde.functions are exported from MakePDEFunctions.nb
calculate.outputs<-function(t,all.pars,extra.fn){
  #get integrals of thymically derived distribution PRE-RECONSTITUTION i.e. between treatment time and reconstitution time
  thymic.donor.pre<-integrate(eval.fn,lower = max(0,t-tau),upper = t+10^-10,t,all.pars,pde.functions[[1]],extra.fn)[[1]]
  thymic.total.pre<-integrate(eval.fn,lower = max(0,t-tau),upper = t+10^-10,t,all.pars,pde.functions[[2]],extra.fn)[[1]]
  
  #get integrals of thymically derived distribution POST-RECONSTITUTION
  thymic.donor.post<-integrate(eval.fn,lower = 0, upper = max(0,t-tau)+10^-10,t,all.pars,pde.functions[[3]],extra.fn)[[1]]
  thymic.total.post<-integrate(eval.fn,lower = 0, upper = max(0,t-tau)+10^-10,t,all.pars,pde.functions[[4]],extra.fn)[[1]]
  
  #get integral of initial age distribution at this time
  initial.total<-integrate(eval.fn,lower = t, upper = t+treatment.age,t,all.pars,pde.functions[[5]],extra.fn)[[1]]
  
  combined.donor<-thymic.donor.pre+thymic.donor.post
  combined.total<-initial.total+thymic.total.pre+thymic.total.post
  
  donor.frac<-combined.donor/combined.total

  return(list(donor.frac,combined.total))
}

calc.mean.loss<-function(t,pars){
  calculate.outputs(t,pars,lambda.function)[[2]]/calculate.outputs(t,pars,1)[[2]]
}

########################################################################################################

findFit<-function(transformed.data.vector,nuVal){
  
  #internal numeric toggles
  maxit<-1.5*10^3
  temp<-10^8
  # range and start values for PERIPHERALLY fitted parameters - see 'header' (defined later, not a great idea but meh...)
  lower<-c(10^-5,10^-5,10^-5,0.01)
  upper<-c(1,1,TonNMax,1)
  start<-c(0.1,0.01,0.01,0.01)
  
  out <- GenSA(par=start,lower = lower, upper = upper, fn = cost.function, control=list(verbose=TRUE,temperature=temp,maxit=maxit),transformed.data.vector,nuVal)
  fitPars<-out["par"][[1]]
  return(fitPars)
  
}

########################################################################################################
# apply approprirate transforms to the donor fraction and count inputs
transform.values<-function(current.pars){
  vals<-sapply(tVec,calculate.outputs,current.pars,1)
  vals<-cbind(Tform2(unlist(t(vals)[,1])),Tform(unlist(t(vals)[,2])))
  return(vals)
}

########################################################################################################

calculate.residuals<-function(current.pars,transformed.data){
  transformed.fits<-transform.values(current.pars)
  residuals<-transformed.fits-transformed.data
  return(residuals)
}
SSRbest<<-10^10
########################################################################################################
# cost function defined as the product of the SSRs of the two peripheral sets we are fitting to
# see busulfan PNAS paper SI for derivation. 
cost.function<-function(inVec,transformed.data,nuVal){
  pars<-c(inVec,nuVal) #add back nu - need to be in header in correct order!
  names(pars)<-header
  pars<-as.data.frame(t(pars))
  residuals<-calculate.residuals(pars,transformed.data)
  
  SSR<-(sum(residuals[,1]^2))*(sum(residuals[,2]^2))
  SSRbest<<-min(SSR,SSRbest)
# print(cbind(SSRbest,SSR,pars))
# if(SSR=="NaN"){
#   browser()
# }
  return(SSR)
}

plot.times<-1:(53*7)
########################################################################################################
# make plots of data and pde.functions for a given set of peripheral and thymic parameters
make.plot<-function(inVec,nuVal){
  pars<-c(inVec,nuVal) #add back nu and mu - need to be in header in correct order!
  names(pars)<-header
  pars<-as.data.frame(t(pars))
  
  vals<-sapply(plot.times,calculate.outputs,pars,1)
  ratio.vals<-data.frame(t=plot.times,ratio=unlist(t(vals)[,1]))
  count.vals<-data.frame(t=plot.times,counts=unlist(t(vals)[,2]))
  
  ratio.plot<-ggplot(ratio.vals,aes(t,ratio))+geom_line(aes(colour="red"))+ylim(0,1)+geom_point(data=ratioData)+theme(legend.position="none")
  count.plot<-ggplot(count.vals,aes(t,counts))+geom_line(aes(colour="red"))+ylim(0,4)+geom_point(data=countData)+theme(legend.position="none")
}

########################################################################################################

Tform<-function(x){
  log(x)
}

Tform2<-function(x){
  asin(sqrt(x))
}

fit.thymic.parameters<-function(precursor.counts,donor.precursor.dp1){
  #hardcode nu function 
  nu.fit<-nls(precursor.counts ~ logA-tVec*nu, start=list(logA=10,nu=0.1))
  nuVal<-coef(nu.fit)[[2]]
  
  #hardcode mu function (time function defined globally)
  #mu.fit<-nls(donor.precursor.dp1 ~ Tform((10^6*exp(logA)*(mu.time.vector/tau)**(1/mu))/E**(nuVal*mu.time.vector)), start=list(logA=10,mu=0.1))
  #muVal<-coef(mu.fit)[[2]]
  
  return(list(nuVal,predict(nu.fit)))
}
########################################################################################################
#stuff we need regardless of if best fit exists already
#we make it so that the fit.thymic.pars function takes TRANSFORMED data so we can use it more easily in the bootstrap part
thymic.pars<-fit.thymic.parameters(Tform(precursor.data[,2]),Tform(mu.data[,2]))
nu.val<-thymic.pars[[1]]
#mu.val<-thymic.pars[[2]]
aBMT<-treatment.age

header<-c("l","r","TonN","mu","nu")

#ratio data first in here!!!!
transformed.data<-cbind(Tform2(ratioData[,2]),Tform(countData[,2]))
########################################################################################################
#see if best fit exists
best.file<-paste(export.folder,"RCD",cell,"Model",toString(model),"P",toString(p),"best.xlsx",sep="")
if(!file.exists(best.file)){
  print("Best fit doesn't exist, calculating:")
  
  best.fit.pars<-c(findFit(transformed.data,nu.val),nu.val)
  names(best.fit.pars)<-header
  best.fit.pars<-as.data.frame(t(best.fit.pars))
  lambda.mean<-mean(sapply(plot.times,calc.mean.loss,best.fit.pars))
  
  best.SSR<-data.frame(SSR=cost.function(as.numeric(best.fit.pars),transformed.data,nu.val))
  lambda.mean<-data.frame(lambda.mean=lambda.mean)
  best.output<-cbind(best.SSR,best.fit.pars,lambda.mean)
  
  
  write.xlsx(best.output,best.file,row.names=FALSE)
  print(best.output)
} else if(file.exists(best.file)){
  print("Best fit exists, importing:")
  best.fit.pars<-read.xlsx(best.file,sheetIndex=1)
  
}
########################################################################################################
#Bootstrap stuff
########################################################################################################

#get TRANSFORMED best-fit values of all three cases (precursor, peripheral count and ratio)
TnuBF<-thymic.pars[[2]]
#TmuBF<-thymic.pars[[4]]
TperipheralBF<-transform.values(best.fit.pars)

#use this to calculate precursor residuals to be bootstrapped along with the count data
nuResidData<-TnuBF-Tform(precursor.data[,2])
#muResidData<-TmuBF-Tform(mu.data[,2])

#have to sample mu separately as different number of points
residuals<-calculate.residuals(best.fit.pars,transformed.data)
pairedMice<-data.frame(nuResidData,residuals[,1],residuals[,2]) 

#setting seed not required, more for debugging purposes (only set once so each replicate is still randomly chosen, just the start point is the same for the whole lot)
set.seed(0)
#define bootstrap functon
doBoot<-function(pairedMice,i){
  #choose residuals with replacement
  boot.choice<-pairedMice[sample(1:length(pairedMice[,1]),replace=TRUE),]
  mu.residual.choice<-muResidData[sample(1:length(muResidData),replace=TRUE)]
  #make new precursor, count, and ratio data (best fit + residuals)
  new.nu.data<-TnuBF+boot.choice[[1]]
  new.mu.data<-TmuBF+mu.residual.choice
  new.ratio.data<-TperipheralBF[[1]]+boot.choice[[2]]
  new.count.data<-TperipheralBF[[2]]+boot.choice[[3]]
  #fit new precursor to get new nu
  thymic.pars<-fit.thymic.parameters(new.nu.data,new.mu.data)
  new.nu.val<-thymic.pars[[1]]
  #new.mu.val<-thymic.pars[[2]]
  
  new.data<-c(new.ratio.data,new.count.data)
  #return new nu, and also fit on count and ratio with this new nu
  output<-c(findFit(new.data,new.nu.val),new.nu.val)
  
  temp.frame<-as.data.frame(t(output))
  names(temp.frame)<-header
  
  function.vals<-t(sapply(tVec,calculate.outputs,temp.frame,1))
  colnames(function.vals)<-c("ratio.vals","count.vals")
  write.xlsx(function.vals,paste(export.folder,"replicates/",toString(i),".xlsx",sep=""),row.names=FALSE)
  return(output)
  
}

# #number of bootstrap replicates
# n.Replicates<-1
# print("Start bootstrap")
# #make parallel cluster thing
# cl <- makeCluster(nLogicalCores)
# registerDoParallel(cl)
# boot.vals<-foreach(i=1:n.Replicates, .combine='rbind',.export=c("N0","aBMT","mu.time.vector","p")) %dopar% { #already exporting tau and treatment time apparently
#   E<-exp(1)
#   library(GenSA)
#   library(xlsx)
#   doBoot(pairedMice,i)
# }
# colnames(boot.vals) <- header
# stopCluster(cl)
# 
# #output bootstrap data
# write.xlsx(boot.vals, file = paste(export.folder,"RCD",cell,"Boot.xlsx",sep=""),
#            sheetName = "Results", row.names = FALSE)