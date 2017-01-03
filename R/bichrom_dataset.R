#chromwoodherb=read.csv("csomewh.csv",header=TRUE, row.names=1) #Removed chromosome numbers larger than 50 for first analysis
#' Cleans and organizes dataset for BiChroM analyses
#' @author Rosana Zenil-Ferguson
#'
#' @details bichrom_dataset organizes chromosome numbers in three ways. First, it transforms chromosome numbers into haploid chromosome numbers, and for chromosome numbers that are not divisible by 2, floor value of chromosome number divided by 2 is taken. Second, changes haploid values larger than specified size to value size+1. And third, for taxa that are in binary state 1, redefines haploid chromosome numbers from size+2 to 2*(size+1). That is if a taxon has 6 chromosomes, state.1=1, and size=50 bichrom_dataset function will recode this individual with value 54 (3 haploid chromosome number, state 0 ranges from 1 to 51, and state 1 ranges from 52 to 102) for analyses with Q_bichrom
#' @param chromwoodherb dataframe with two columns. First column contains chromosome numbers, second column has values of binary state. Row.names of data frame are taxon names
#' @param size maximum haploid chromosome number allowed in analyses. See details.
#' @param state.0 value of binary state 0 in second column of chromwoodherb by default "H" indicating herbaceous.
#' @param state.1 value of binary state 1 in second column of chromwoodherb by default "W" indicating woody
#' @return Returns a data frame with two columns. First column taxon name and second column haploid chromosome numbers codified by state. See details.
#' @seealso Q_bichrom
#' @export
bichrom_dataset<-function(chromwoodherb, size=50, state.0="H",state.1="W"){
  print("Make sure the taxa in your tree and dataset are in the same order")
  # make the space.state thing.
  last.state<-size
  csomenumber<-as.numeric(chromwoodherb[,1])
  aux1<-csomenumber/2
  aux2<-which(aux1%%1!=0)
  haploidnumber<-aux1
  haploidnumber[aux2]<-floor(aux1[aux2])
  largerthan<-which(haploidnumber >last.state)
  taxa.morethanlaststate<-chromwoodherb[largerthan,] #Taxa with a lot more of
  haploidnumber[largerthan]<-last.state+1
  type<-chromwoodherb[,2]
  herbaceous<-which(type==state.0)
  woody<-which(type==state.1)
  bichromvals<-rep(0,dim(chromwoodherb)[1])
  bichromvals[herbaceous]<-haploidnumber[herbaceous]
  bichromvals[woody]<-last.state+haploidnumber[woody]+1
  ########################################################################
  ##Formal input
  bichrom.dataset<-data.frame(Taxon=row.names(chromwoodherb),bichrom_values=bichromvals)
  write.table(bichrom.dataset,file="bichromvals.txt",sep=",",row.names=FALSE,col.names=FALSE)
  return(bichrom.dataset)
}
