#' Calculation of Strength of the Connectivity among multiple Brain Regions
#'
#' This function can calculate strength of the connectivity among multiple brain regions via kernel canonical correlation.
#' Permutaiton test is employed to assess the statistical significance.
#' @author Xubo Yue, Chia-Wei Hsu (tester), Jian Kang (maintainer)
#' @param imageDat There are two options for this argument: (1)(vectors of) imaging data (with extension .nii or .nii.gz).
#' You do not need to add extension in your argument, simply type in the name of file is enough. (2) the processed image data
#' produced by nii2RData function. No matter which option you choose, the result will be the same (as long as you use the same
#' dataset). When you would like to read and save nii data first and run kcca test later, you can first save output of nii2RData
#' and, in the future, choose option (2) to run kcca in order to save time (avoid read same dataset multiple times).
#' @param region user-specified multiple brain regions (as vector, for example, c(1,2,30)).
#' @param resolution the resolution of your region data. It can take "2mm" as default.
#' If user would like to use 3mm resolution, type in "3mm".
#' @param saveName whether to save processed imaging data.
#' If you do not have enough space or do not want to use space to store processed data, just type in "None" (default);
#' Else you need to specify name in this argument. For example, saveName="myName.RData".
#' @param kernel the kernel function used in training and predicting.
#' The default kernel is the radial basis kernel function "Gaussian" (rbfdot). Use "?kernlab::kcca" to find more available kernels.
#' @param regionCode the region code provided by user or default. It should have 3 columns with index, region code and region name.
#' @param niiFile2 the nii region file you would like to use. It has default nii file and can be left as blank.
#' @param imgPath the directory where your nii file(s) is (are) located. It chooses your current working directory as default.
#' @param datPath the directory where you would like to store .RData file(s). It chooses your current working directory as default.
#' @param parallel whether to use parallel computing. Type FALSE as not using parallel and TRUE as using parallel. Parallel is not
#' recommended in local computer as it may slow down your system. Use parallel in cluster is preferred
#' @param loc this argument can accept argument "local" or "cluster". if you choose to use parallel computing, please specify whether
#' you run your code in your local computer or cluster.
#' If you did not choose parallel computing, then "local" or "clutser" makes no difference.
#' @param perm number of permutation. Default time is 50.
#' @param saveData whether to save output as R data. Type "None" as not save output.
#' Type name of file if you would like to save. For example, "output_two.RData".
#' @return (lists of) list of brain regions, p-value, region type ("two" or "multiple"), and sregion name.
#' @details (1) Kernel canonical correlation analysis (KCCA) can explore the nonlinear relationship between two variables.
#' It transformed sample vectors into the Hilbert space and maximize correlation coefficient by solving quadratically regularized
#' Lagrangean function.
#' Refer to Kang's paper for more details: Kang J, Bowman FD, Mayberg H, Liu H (2016). "A depression network of functionallyconnected
#' regions discovered via multi-attribute canonical correlation graphs."NeuroImage,141, 431-441.
#' (2) Use rgl.snapshot() function if you would like to save plot but forgot to use TRUE in screenShot argument.
#' @references \url{https://www.ncbi.nlm.nih.gov/pubmed/27474522}
#' @export
#' @examples
#' #result0<-permkCCA_multipleRegion(imageDat=c("preproc_con21_rest_MNI_2mm", "preproc_con23_rest_MNI_2mm"), region=c(1,60,70), regionCode="RegionList.txt", niiFile2="AAL_MNI_2mm.nii") #user can provide region data and region list
#' #result1<-permkCCA_multipleRegion(imageDat="UM_1_0050272_func_preproc", region=c(1,75), resolution="3mm", regionCode="RegionList90.txt", niiFile2="AAL_90_3mm.nii") #for 3mm resolution, user must provide region data and region list
#' #result2<-permkCCA_multipleRegion(imageDat="preproc_con21_rest_MNI_2mm", region=c(1,2,3)) #user can use default region data and region list only for 2mm resolution
#' #testcase1 <- nii2RData(niiFile1 = "preproc_con21_rest_MNI_2mm", saveName = "None") #one can process nii file first,
#' #result1<-permkCCA_multipleRegion(imageDat = testcase1, region = c(1,5,10)) #and then use processed nii file as input.
#' #result4<-permkCCA_multipleRegion(imageDat = testcase1, region = c(1,5)) #and then use processed nii file as input.
#' #result5<-permkCCA_multipleRegion(imageDat=c("preproc_con21_rest_MNI_2mm", "preproc_con23_rest_MNI_2mm"), region=c(1,60,90), regionCode="RegionList.txt", niiFile2="AAL_MNI_2mm.nii")
#' #summary_kcca(result1, 0.05, 1, "markdown")
permkCCA_multipleRegion = function(imageDat, region, resolution="2mm", saveName="None", kernel= "rbfdot",
                                   regionCode="", niiFile2="", imgPath=getwd(), datPath=getwd(), parallel=FALSE,
                                   loc="local", perm=50, saveData="None"){

  cat("\n", "Checking data format...","\n")

  if(length(region)<2) {warning("length of your vector should be at least 2, please check your input vector!", call. = FALSE); return(0)}
  if(sum(table(region)>1)) {warning("your input numbers should be different, please check your input vector!", call. = FALSE); return(0)}
  if(class(imageDat)!="list") imageDat<-nii2RData(imageDat, resolution, saveName, regionCode, niiFile2, imgPath, datPath)
  else cat("\n", "processed image data detected, do not need to read nii data again...","\n")

  if(length(region)==2){

    temporary_variable <- imageDat[[length(imageDat)]]
    if(max(region) > length(temporary_variable[,1])) {warning("region index is out of upper bound.", call. = FALSE); return(0)}
    if(min(region) <= 0) {warning("region index is out of lower bound.", call. = FALSE); return(0)}

    cat("Calculating KCC...","\n")

    if(regionCode!=""){
      regionCode<-read.table(regionCode)
      if(dim(regionCode)[2]==3){
        regionNum<-regionCode[,2]
        names(regionNum)<-regionCode[,3]
      }
      else stop("Region list can only have 3 columns.")
    }
    else{
      regionNum<-c(1:116)
      names(regionNum)<-node[,6]
    }

    final_output<-NULL

    for(i in 1:(length(imageDat)-1)){
      partImgdat<-partition(imageDat[[i]])
      results = perm_kCCA(partImgdat[[region[1]]],partImgdat[[region[2]]],permNum=perm,kernel="rbfdot")
      output<-NULL
      output[[1]]<-c(region[1],region[2])
      output[[2]]<-results[[3]]
      output[[3]]<-"two"
      output[[4]]<-c(names(regionNum[region[1]]),names(regionNum[region[2]]))
      final_output[[i]]<-output
    }

    final_output[[(length(imageDat))]] <- imageDat[[(length(imageDat))]]
    if(saveData!="None") save(final_output, file=saveData)
    return(final_output)
  }


  if(length(region)>2){
    temporary_variable <- imageDat[[length(imageDat)]]
    if(max(region) > length(temporary_variable[,1])) {warning("region index is out of upper bound.", call. = FALSE); return(0)}
    if(min(region) <= 0) {warning("region index is out of lower bound.", call. = FALSE); return(0)}
    if(regionCode!=""){
      regionCode<-read.table(regionCode)
      if(dim(regionCode)[2]==3){
        regionNum<-regionCode[,2]
        names(regionNum)<-regionCode[,3]
      }
      else stop("Region list can only have 3 columns.")
    }
    else{
      regionNum<-c(1:116)
      names(regionNum)<-node[,6]
    }

    cat("Calculating KCC...","\n")

    final_output<-NULL

    for(j in 1:(length(imageDat)-1)){

      partImgdat <- partition(imageDat[[j]])
      pairs <- expand.grid(region,region) #combination of two pairs
      pairs <- unique(t(apply(pairs,1,sort))) #remove same value with reverse order
      pairs <- pairs[pairs[,1]!=pairs[,2],] #remove same value in two column

      if(parallel==TRUE){
        if(loc=="local")  no_cores <- parallel::detectCores() - 1
        else no_cores <- parallel::detectCores()
        if(loc!="local") cl <- parallel::makeCluster(no_cores, type="FORK")
        else cl <- parallel::makeCluster(no_cores)
        if(loc=="local") parallel::clusterExport(cl, list("kcor", "kcca"))
        else {
          parallel::clusterEvalQ(cl, library(kernlab))
          parallel::clusterExport(cl, list("perm_kCCA", "kcor", "kcca", "sapply_pb", "regionList", "partition", "rbfdot"), envir=environment())

        }

        results<-parallel::parLapply(cl, 1:nrow(pairs), function(x)
          perm_kCCA(partImgdat[[pairs[x,1]]],partImgdat[[pairs[x,2]]],permNum=perm,kernel= "rbfdot"))

        parallel::stopCluster(cl)
      }
      else{
        results<-lapply(1:nrow(pairs), function(x) brainKCCA::perm_kCCA(partImgdat[[pairs[x,1]]],partImgdat[[pairs[x,2]]],permNum=perm,kernel="rbfdot"))
      }

      output<-NULL
      output[[1]]<-pairs

      output[[2]]<-results[[1]][[3]]
      for(i in 1:choose(length(region),2)){
        output[[2]][[i]]<-results[[i]][[3]]
      }
      output[[3]]<-"multiple"
      output[[4]]<-c(names(regionNum[pairs[1,1]]),names(regionNum[pairs[1,2]]))
      for(i in 2:nrow(pairs)){
        output[[4]]<-rbind( output[[4]], c(names(regionNum[pairs[i,1]]),names(regionNum[pairs[i,2]]))  )
      }

      tempPValue<-NULL
      for(i in 1:length(output[[2]]))
        tempPValue<-c(tempPValue,output[[2]][[i]])

      adjustedPValue<-p.adjust(tempPValue,"BY")
      for(i in 1:length(output[[2]]))
        output[[2]][[i]]<-adjustedPValue[i]

      final_output[[j]]<-output
    }

    final_output[[(length(imageDat))]] <- imageDat[[(length(imageDat))]]
    if(saveData!="None") save(final_output, file=saveData)
    return(final_output)
  }

}
