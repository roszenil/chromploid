#' Cleans and organizes dataset for PloidEvol analyses
#'
#' @details This function cleans and organizes a vector of ploidies assuming ploidy values range from 2x to 12x. It deletes values larger than 12x and redefines ploidy values as value-1 since Q_ploidevol indices range from 1 to 11. If larger ploidy values are needed user needs to specify own Q.FUN to use in negloglikelihood.
#' @seealso \code{\link{negloglikelihood}}
#' @param ploidies a vector with ploidy values, and taxon row.names
#' @param size maximum ploidy value allow for analyses, by default is 12. See details for more information.
#' @return ploidy.dataset a data frame whose first column is taxon name and second column the ploidy values-1. Saves the data frame as .txt file for future references. See details.
#' @export
ploidevol_dataset<-function(ploidies, size=12){
  larger.than12<-which(ploidies[ ,1]>12)
  ploidyvalues<-ploidies[-larger.than12, ]
  ploidy.dataset<-data.frame(Taxon=row.names(ploidies)[-larger.than12] ,ploidvals=as.numeric(ploidyvalues)-1)
  write.table(ploidy.dataset,file="ploidvals.txt",sep=",",row.names=FALSE,col.names=FALSE)
  return(ploidy.dataset)
}
