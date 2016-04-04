#workDir<-"/home/graeme/Dropbox/Thymic Development/BusulfanProjectFiles/"
workDir<-"C:/Users/Graeme/Dropbox/Thymic Development/BusulfanProjectFiles"
setwd(workDir)
########################################################################################################################################################################
#workDir<-"C:/Users/Graeme/Dropbox/Thymic Development/BusulfanProjectFiles"

source("BusulfanPackages.R")
raw.Data<-read.xlsx("RawData/ChimeraCountsCompiled NO WT.xlsx",sheetIndex=1,header=TRUE)
n.Mice<-nrow(raw.Data)

#go through and remove cols which only contain NA (i.e. missing data - read.xlsx also seems to pick up phantom columns for some reason)
del.Cols<-colSums(is.na(raw.Data))
del.Cols.Names<-as.vector(names(del.Cols[del.Cols!=n.Mice]))

#re-define the raw data dataframe.
raw.Data<-raw.Data[,del.Cols.Names]

#order by age post-BMT
raw.Data<-raw.Data[order(raw.Data$days.post.bmt),]
header<-names(raw.Data)
mice<-raw.Data[,1]
########################################################
########################################################
#things you might want to edit here!
population.List<-c("donor","host")

cell.Types.nai<-c("4nai","8nai","gdT44lo","gdT44hi")
cell.Types.others<-c("4mem","8Tem","8Tcm")

#thymic precursor order must match the order of the joined list c(cell.Types.nai,cell.Types.others)
thymic.precursors<-c("sp4","sp8","dn4","dn4","sp4","sp8","sp8")
peripheral.precursors<-c("4nai","8nai","8nai")


cell.Types<-c(cell.Types.nai,cell.Types.others)

if(length(thymic.precursors)!=length(cell.Types)){
  stop("thymic precursor number not same as master cell number, stopping.")
}
organ.Set<-c("LN","SP")
expFolder<-"PreparedDataR/"

#defined reconstitution time NOT USED ANYMORE LEAVE THIS ALONE  truncation done in later codes
t.Shift<-0.

#mice to delete
naughty.mice<-c("TH108-30","TH108-29","TH108-40")

########################################################
########################################################

raw.Data<-raw.Data[!(raw.Data$mouse.id  %in% naughty.mice),]
raw.Data<-raw.Data[raw.Data$days.post.bmt >=t.Shift ,]
n.Mice<-nrow(raw.Data)
time.vector<-raw.Data[,"days.post.bmt"]
names(time.vector)="t"
########################################################
collect.data<-function(cell,organ.list,set){
  donor.counts<-rep(0,n.Mice)
  host.counts<-rep(0,n.Mice)
  
  for(i in 1:length(organ.list)){
    organ<-organ.list[i]
    donor<-paste(organ,".donor.",cell,sep="")
    host<-paste(organ,".host.",cell,sep="")
    
    donor.counts<-donor.counts+set[,donor]
    host.counts<-host.counts+set[,host]
  }  
  return(list(donor.counts,host.counts))
}
########################################################
#this makes the UNSCALED peripheral donor fraction
make.donor.frac<-function(cell,organ.list,set){
  pooled.data<-collect.data(cell,organ.list,set)
  
  donor.counts<-pooled.data[[1]]
  host.counts<-pooled.data[[2]]
  
  donor.frac<-donor.counts/(donor.counts+host.counts)
  return(donor.frac)
}
########################################################
make.counts<-function(cell,organ.list,set){
  pooled.data<-collect.data(cell,organ.list,set)
  
  donor.counts<-pooled.data[[1]]
  host.counts<-pooled.data[[2]]
  
  totals<-donor.counts+host.counts
  return(totals)
}
########################################################

dp1.frac<-make.donor.frac("dp1",c("TH"),raw.Data)
dp1.frac<-data.frame(t=time.vector,dp1.frac=dp1.frac)

thymic.precursor.fracs<-make.donor.frac(thymic.precursors,c("TH"),raw.Data)  
thymic.precursor.counts<-make.counts(thymic.precursors,c("TH"),raw.Data)  

names(thymic.precursor.fracs)<-thymic.precursors
names(thymic.precursor.counts)<-thymic.precursors

#unscaled donor frac in periphery for all cells
peripheral.fracs<-make.donor.frac(cell.Types,organ.Set,raw.Data)  
peripheral.counts<-make.counts(cell.Types,organ.Set,raw.Data)  
names(peripheral.fracs)<-cell.Types
names(peripheral.counts)<-cell.Types

#deal with missing data - replace all NA's with -1 so that we can easily remove from both counts and donor frac later
peripheral.fracs[is.na(peripheral.fracs)]<- -1
peripheral.counts[is.na(peripheral.counts)]<- -1

#make peripheral donor frac rescaled by thymic precursors
rescaled.peripheral.fracs<-foreach(i = 1:ncol(peripheral.fracs), .combine=cbind) %do% {
  peripheral.fracs[i]/thymic.precursor.fracs[i]
}

#now also take any memory cells and make the peripheral donor frac rescaled by appropraite naive frac
#this is done independent of any scaling to thymic populations so changing those precursors *shouldn't* affect this...
rescaled.peripheral.fracs.mem<-foreach(i = 1:length(cell.Types.others), .combine=cbind) %do% {
  peripheral.fracs[cell.Types.others[i]]/peripheral.fracs[peripheral.precursors[[i]]]
}

# while considering pre-reconstitution data we pickup some small mismatches where we divide by zero
# can also just happen because we divide by (say) a zero value of DP1 but the naive value was small but non zero so get infinity - assume these should all be zero and KEEP rather than delete
# (can't have non-zero donor in periph without non-zero in thymus)
for(i in 1:ncol(rescaled.peripheral.fracs)){
  set<- rescaled.peripheral.fracs[is.nan(rescaled.peripheral.fracs[[i]]),]
  if(nrow(set)>0){
    rescaled.peripheral.fracs[is.nan(rescaled.peripheral.fracs[[i]]),][,i]<-0
  }
  set<- rescaled.peripheral.fracs[is.infinite(rescaled.peripheral.fracs[[i]]),]
  if(nrow(set)>0){
    rescaled.peripheral.fracs[is.infinite(rescaled.peripheral.fracs[[i]]),][,i]<-0
  }
}

for(i in 1:ncol(rescaled.peripheral.fracs.mem)){
  set<- rescaled.peripheral.fracs.mem[is.nan(rescaled.peripheral.fracs.mem[[i]]),]
  if(nrow(set)>0){
    rescaled.peripheral.fracs.mem[is.nan(rescaled.peripheral.fracs.mem[[i]]),][,i]<-0
  }
  set<- rescaled.peripheral.fracs.mem[is.infinite(rescaled.peripheral.fracs.mem[[i]]),]
  if(nrow(set)>0){
    rescaled.peripheral.fracs.mem[is.infinite(rescaled.peripheral.fracs.mem[[i]]),][,i]<-0
  }
}

########################################################
#get pre-reconstitution donor data in the thymus for use in PDE model, see definition of 'mu' in age dependent loss stuff.
mu.cells<-unique(thymic.precursors)
mu.data<-raw.Data["days.post.bmt"]
for(i in 1:length(mu.cells)){
  mu.data<-cbind(mu.data,raw.Data[paste("TH.donor.",mu.cells[i],sep="")]/dp1.frac[,2])
  
}
mu.data<-mu.data[mu.data$days.post.bmt<=42,]
names(mu.data)<-c("t",mu.cells)
########################################################
# It's possible (i.e. gdT44lo) that the data sometimes gives rescaled fractions >1 so we replace those slightly larger than 1 with 1. This is just due to small errors in the denominator
# Anything larger than 1.X is kept so it throws errors as larger than 1.1 means a big problem probably!
rescaled.peripheral.fracs[rescaled.peripheral.fracs>1]<-1.
#rescaled.peripheral.fracs[rescaled.peripheral.fracs>1&rescaled.peripheral.fracs<1.1]<-1.
#rescaled.peripheral.fracs[rescaled.peripheral.fracs>1.1]<--1.

#7th Dec need to revisit this, should really drop points instead of making them all 1 but need to learn how to remove rows based on logical check of any column

rescaled.peripheral.fracs<-cbind(time.vector,rescaled.peripheral.fracs)
peripheral.counts<-cbind(time.vector,peripheral.counts)
thymic.precursor.counts<-cbind(time.vector,thymic.precursor.counts)
rescaled.peripheral.fracs.mem<-cbind(time.vector,rescaled.peripheral.fracs.mem)
########################################################
export.Table<-function(frame,file.name){
  for(i in 2:length(frame)){
    temp<-cbind(frame[,1]-t.Shift,frame[,i])
  
    #here we remove the missing data in the 8Tem/cm case (manifests in the form of negative values, see replacement done on peripheral.fracs and peripheral.counts)
    temp<-temp[temp[,2]>=0,]                             
    file<-paste(expFolder,names(frame)[[i]],file.name,sep="")
    write.table(temp,file,sep="\t",row.names=F,col.names=F)
  }
}
########################################################
export.Table(rescaled.peripheral.fracs,"SPRatio.txt")
export.Table(peripheral.counts,"TotalCounts.txt")
export.Table(thymic.precursor.counts,"TotalCounts.txt")
export.Table(rescaled.peripheral.fracs.mem,"NaiRatio.txt")

export.Table(mu.data,"muData.txt")
########################################################
#ki67
########################################################
ki67.data<-read.xlsx("RawData/ki67Data NO WT.xlsx",sheetIndex=1,header=TRUE)
ki67.data<-ki67.data[ki67.data$days.post.bmt>=30*7,] #only keeping ki67 data post 30wks
ki67.data<-ki67.data[order(ki67.data$days.post.bmt),]
ki67.mice<-ki67.data[,1]
ki67.mice.counts<-raw.Data[raw.Data$mouse.id %in% ki67.mice,]
ki67.time.vector<-ki67.data[,"days.post.bmt"]
#get number of Ki67+ cells as well as total number of cells
get.ki67.numbers<-function(cell,organ,population){
  label<-paste(organ,".",population,".",cell,sep="")
  ki67perc<-ki67.data[,label]/100
  counts<-ki67.mice.counts[,label]

  return(data.frame(ki67pos=ki67perc*counts,ki67total=counts))
}
#now pool these numbers across organs and calculate final (weighted average) ki67+%
for(k in 1:length(population.List)){
  population<-toString(population.List[k])
  result<-ki67.time.vector
  
  for(j in 1:length(cell.Types)){
    cell<- cell.Types[j]
    pooledTemp<-rep(0,length(ki67.mice))
    for(i in 1:length(organ.Set)){
      pooledTemp<- pooledTemp + get.ki67.numbers(cell,organ.Set[i],population)
    } #END I LOOP
    
    #exported as percentages NOT fractions
    result<-cbind(result,100*pooledTemp[,1]/pooledTemp[,2])
  } #END J LOOP
  result<-as.data.frame(result)
  names(result)<-c("t",cell.Types)
  fileE<-paste("ki67",population,".txt",sep="")

  export.Table(result,fileE)
} #END K LOOP

########################################################################################################################################################################