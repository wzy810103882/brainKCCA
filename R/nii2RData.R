#' Read nii image(s) file(s) and save it (them) as Rdata
#' 
#' This function can take nii data file(s) (the file extension is "nii.gz" or "nii") as input 
#' (you only need to specify the name of file before extension) and transform them 
#' into RData file(s) (in the form of .RData).
#' @author Xubo Yue, Chia-Wei Hsu (tester), Jian Kang (maintainer)
#' @param niiFile1 the nii data file(s) you would like to use. 
#' If you would like to read one or more files, type in the name in the string vector. 
#' @param resolution the resolution of your region data. It can take "2mm" as default. 
#' If user would like to use 3mm resolution, type in "3mm".
#' @param saveName whether to save processed imaging data. 
#' If you do not have enough space or do not want to use space to store processed data, just type in "None" (default);
#' Else you need to specify name in this argument. For example, saveName="myName.RData".
#' @param regionCode the region code provided by user or default. 
#' It should have 3 columns with index, region code and region name.
#' @param niiFile2 the nii region file you would like to use. It has default nii file and can be left as blank.
#' @param imgPath the directory where your nii file(s) is (are) located. 
#' It chooses your current working directory as default.
#' @param datPath the directory where you would like to store .RData file(s). 
#' It chooses your current working directory as default.
#' @return the processed imaging data.
#' @details This function must accept the nii imaging data as argument to manipulate the raw data to processed data. 
#' The file extension is "nii.gz" or "nii".
#' You only need to specify the name of file without extension in the function argument. 
#' For example, your interested imaging file is called "brain.nii.gz" or "brain.nii".
#' You only need to write niiFile1="brain" in the first argument rather than "brain.nii.gz" or something else.
#' The package can either uses its default region code and region data or uses user-defined region information. \cr
#' Here are some notes which may be useful: \cr
#' (1) When reading multiple imaging files, make sure they have the same resolution, region code and region data.
#' (2) If you would like to read multiple imaging files, and if you did not use saveName="None" argument, it can 
#' sometimes take a large storage. For example, if you read 100 files, then 100
#' "RData" files will be generated and thus take a large storage space. Thus, saveName="None" is highly recommended.
#' @export
#' @examples  
#' #testcase1 <- nii2RData(niiFile1 = "preproc_con21_rest_MNI_2mm", saveName = "None") 
#' #testcase2 <- nii2RData(niiFile1 = "preproc_con21_rest_MNI_2mm", saveName = "Example1.RData") 
#' #testcase3 <- nii2RData(niiFile1 = c("preproc_con21_rest_MNI_2mm", "preproc_con21_rest_MNI_2mm"), saveName = "Example2.RData") 
#' #nii2RData(niiFile1 = c("data1_3mm", "data2_3mm"), resolution="3mm", saveName = "Example2.RData") 
#' #nii2RData(niiFile1 = "data1_3mm", resolution="3mm", saveName = "Example2.RData")
nii2RData<-function(niiFile1, resolution="2mm", saveName="None", regionCode="", niiFile2="", imgPath=getwd(), datPath=getwd()){
  
  cat("Reading nii files...this progress may take a long time...", "\n")
  
  if(niiFile2=="") {
    temp<-NULL
    for(i in 1:dim(regionIdx)[1]){
      if(i==1) pb <- txtProgressBar(min = 0, max = dim(regionIdx)[1], style = 3)
      Sys.sleep(0.1)
      setTxtProgressBar(pb, i)
      temp<- rbind(temp,t(as.matrix(colMeans(coords2mm[which(aal_region_2mm==regionIdx[i,1],arr.ind = F),]))))
      
    }
    regionIdx<-cbind(temp,regionIdx)
    user<-list(aal_region_2mm, regionIdx)
  }
  else user<-regionList(niiFile2, regionCode, resolution)
  
  region_select = (user[[1]]>0)
  regioncode = user[[1]][region_select]
  
  imgdat=list()
  
  for(i in 1:length(niiFile1)){
    cat("\n", "reading data ", i, "\n" )
    input_img = oro.nifti::readNIfTI(fname=file.path(imgPath,paste(niiFile1[i],".nii.gz", sep="")))
    time_points = dim(input_img)[4]
    imgdat[[i]] = matrix(0,nrow=time_points,ncol=sum(region_select))
    for(j in 1:time_points){
      if(j==1) pb <- txtProgressBar(min = 0, max = time_points, style = 3)
      Sys.sleep(0.1)
      setTxtProgressBar(pb, j)
      temp_img = input_img[,,,j]
      imgdat[[i]][j,] = temp_img[region_select]
    }
    division<- cumsum(table(regioncode))
    colnames(imgdat[[i]])<-imgdat[[i]][1,]
    colnames(imgdat[[i]])[c(1,1+division[-length(division)])]<-"V1"
  }

  imgdat[[length(niiFile1)+1]]<-user[[2]]
  
  if(saveName!="None") { 
    cat("\n", "saving data...","\n")
    save(imgdat,file=file.path(datPath,saveName))  
  }
  
  cat("\n", "Completed...", "\n")
  close(pb)

  return(imgdat)
}

