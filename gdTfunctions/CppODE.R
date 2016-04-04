library(odeintr)
library(xlsx)
############################################################################################################

#get c++ files
cFnFileTemplate<-"gdTfunctions/CppODEs.txt"
options(warn=-1)
odeSystem <- paste(readLines(cFnFileTemplate), collapse=" ")
odeICs<-parse(text=read.csv("gdTfunctions/CppICs.csv", header = FALSE,quote = "", stringsAsFactors=F)[[1]])
options(warn=0)

############################################################################################################


cppHeaderCode<-"
#include <math.h>
"

cppDef<-"
double E = 2.718282;
"
fns<-c("Ad","Ah","Bd","Bh")
nFns<-length(fns)

#compile solver - takes time!
internalHeader<-c("chi","TonA","beta","lambdaA","lambdaB","nu","muA","muB","A0") # fixed here, should be made more polymorphic for other types of het considered! i.e. add kappa for tempHet
compile_sys("compiledSystem", odeSystem, internalHeader,sys_dim=nFns,headers=cppHeaderCode,globals=cppDef)

############################################################################################################
NDSolve <- function(tau,parameters,IC){
  #set pars
  parL<-as.list(parameters)
  compiledSystem_set_params(chi = with(parL,chi),
                            TonA = with(parL,TonA),
                            beta = with(parL,beta),
                            lambdaA = with(parL,lambdaA),
                            lambdaB = with(parL,lambdaB),
                            nu=with(parL,nu),
                            muA=with(parL,muA),
                            muB=with(parL,muB),
                            A0=with(parL,A0))
  
  resultRaw<-compiledSystem_at(IC, tau)
  result<-resultRaw[,-1] #drop time column from dataframe, not needed really.
  colnames(result)<-fns
  chiVal<-with(parL,chi)

  Ad<-result$Ad
  Ah<-result$Ah
  Bd<-result$Bd
  Bh<-result$Bh
  
  total.A<-Ad+Ah
  total.B<-Bd+Bh

  A.data<-data.frame(t=tau,counts=total.A,ratio=Ad/(chiVal*total.A))
  B.data<-data.frame(t=tau,counts=total.B,ratio=Bd/(chiVal*total.B))

  return(list(A.data,B.data))
}

############################################################################################################
# #testbed
# impPars<-read.xlsx("4naibestFit (copy).xlsx",sheetIndex=1)
# impParsN<-as.numeric(t(impPars)[-1,])
# impParsH<-t(impPars)[-2,]
# names(impParsN)<-impParsH
# impParsN<-as.data.frame(t(impParsN))
# impParsN<-cbind(impParsN,chi=0.5)
# 
# testPars<-data.frame(nu=0.004289,chi=0.5,TonN=0.014,k=8.71,mu=0.484,lambdaR=0.034,lambdaN=0.005689,R0=1.8313,betaN=0.7712,betaR=0.001,NT=3.005)
# #testPars<-impParsN
# IC<-NULL
# for(i in 1:length(odeICs)){
#   IC[[i]]<-with(as.list(testPars),eval(odeICs[[i]]))
# }
# tauList<-seq(from =0, to = 50, by=0.01)
# 
# out<-NDSolve(tauList,testPars,IC)
# totalDat<-cbind(tauList+6*7,out[[1]])
# ratioDat<-cbind(tauList+6*7,out[[2]])
# 
# #write.xlsx(totalDat,"cd8TestTotal.xlsx")
