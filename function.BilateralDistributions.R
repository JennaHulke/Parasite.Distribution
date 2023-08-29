### function to analyze observed data of parasite intensities that are distributed between 2 discrete habitats (paired organs)
BilateralDisturbutions <- function(file){
  mydata <- file  
  intensity <- mydata[,3]
  #return(index = mean(mydata.sum[,8])) ##index of dispersion
  mydata.sum <-matrix(data = NA, nrow = nrow(mydata), ncol= 9)
  colnames(mydata.sum) <- c("Group A", "Group B", "total.intensity", "diff", "abs.diff", "Variance", "mean", "ID", "single")
  mydata.sum[,1] <- mydata[,1] 
  mydata.sum[,2] <- mydata[,2]
  mydata.sum[,3] <- mydata.sum[,1]+ mydata.sum[,2] ##gives total intensity per host
  mydata.sum[,4] <- mydata.sum[,1] - mydata.sum[,2] ## difference between glands (this is directional)
  mydata.sum[,5] <- abs(mydata.sum[,4]) ##gives absolute value of difference between glands to make it non-directional
  
  for (i in 1:nrow(mydata.sum)){
    mydata.sum[i,6] <- var(x=c(mydata.sum[i,1], mydata.sum[i,2])) ##variance between the glands
    mydata.sum[i,7] <- mean(c(mydata.sum[i,1], mydata.sum[i,2])) ##mean of glands 
    mydata.sum[i,8] <- mydata.sum[i,6]/mydata.sum[i,7] ##Index of dispersion
    mydata.sum[i,9] <- sum(mydata.sum[i,1] ==1, mydata.sum[i,2]==1) ##number of single infections
  }
  index_of_dispersion <- mean(mydata.sum[,8])
  print(paste("Index of Dispersion:", index_of_dispersion))
}
