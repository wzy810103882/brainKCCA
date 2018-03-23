#read nii region file provided by user and transform it into RData. 
#oh<-regionList("AAL_MNI_2mm.nii","RegionList.txt")
#example<-regionList("AAL_MNI_2mm.nii", "RegionList.txt")

regionList<-function(regionData, regionCode, resolution="2mm"){
  
  cat("reading and manipulating the regionData...", "\n")
  largeMatrix<-oro.nifti::readNIfTI(regionData)
  longVector<-expand.grid(largeMatrix)
  #Same as below
  #longVector<-0
  #for(i in 1:prod(dim(largeMatrix))) longVector[i] = largeMatrix[[i]]
  
  cat("reading and manipulating the regionCode...", "\n")
  regionCode<-read.table(regionCode)
  
  if(dim(regionCode)[2]!=3) stop("Region list can only have 3 columns.")
  
  center2mm = c(46,64,37)
  if(resolution=="2mm") coords2mm = expand.grid(-2*(1:91-center2mm[1]),2*(1:109-center2mm[2]),2*(1:91-center2mm[3]))
  if(resolution=="3mm") coords2mm = expand.grid(-3*(1:91-center2mm[1]),3*(1:109-center2mm[2]),3*(1:91-center2mm[3]))
  temp<-NULL
  for(i in 1:dim(regionCode)[1])
    temp<- rbind(temp,t(as.matrix(colMeans(coords2mm[which(largeMatrix==regionCode[i,2],arr.ind = F),]))))
  regionCode<-cbind(temp, regionCode)
  return(list(longVector, regionCode))
  
}






