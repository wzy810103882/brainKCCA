partition<-function(imageDat){
  position<-which(colnames(imageDat)=="V1")
  partImgdat<-list()
  for(i in 1:(length(position)-1)){
    partImgdat[[i]]<-imageDat[,position[i]:position[i+1]-1]
  }
  partImgdat[[length(position)]]<-imageDat[,utils::tail(position,n=1):dim(imageDat)[2]]
  return(partImgdat)
}
