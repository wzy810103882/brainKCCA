#' Calculation of Strength of the Connectivity among multiple Brain Regions
#' This function is the core for kernel canonical correlation. Generally you do not need to use this function unless you are famaliar with kcca
#' algorithm.
#' @author Xubo Yue, Chiawei Hsu (tester), Jian Kang (maintainer)
#' @rdname perm_kCCA
#' @param x region 1, a matrix containing data index by row.
#' @param y region 2, a matrix containing data index by row.
#' @param sig inverse kernel width for the Radial Basis kernel function "rbfdot" and the Laplacian kernel "laplacedot".
#' @param gama regularization parameter (default: 0.1).
#' @param ncomps number of canonical components (default: 1).
#' @param permNum number of permutation.
#' @param kernel type of kernel.
#' @return (lists of) list of brain regions, permutation coefficient, p-value, region name and region type ("two" or "multiple").
#' @details Kernel canonical correlation analysis (KCCA) can explore the nonlinear relationship between two variables.
#' It transformed sample vectors into the Hilbert space and maximize correlation coefficient by solving quadratically regularized Lagrangean function.
#' Refer to Kang's paper for more details: Kang J, Bowman FD, Mayberg H, Liu H (2016). "A depression network of functionallyconnected regions discovered via multi-attribute canonical correlation graphs."NeuroImage,141, 431-441.
#' @references \url{https://www.ncbi.nlm.nih.gov/pubmed/27474522}
#' @importFrom kernlab kcor kcca
#' @export
perm_kCCA = function(x,y,sig=0.1,gama=0.1,ncomps=1,permNum = 500,kernel="rbfdot") {
  n = nrow(x)
  res =  kcca(x,y,kpar=list(sigma=sig),gamma = gama, ncomps = ncomps, kernel=kernel)
  rcorcoef = abs(kcor(res)[1])
  permcoef = sapply_pb(1:permNum, function(i) return(abs(kcor(kcca(x[sample(1:n,n),],y,kpar=list(sigma=sig),gamma = gama, ncomps = ncomps)))[1]))
  return(list(permcoef = permcoef,rcorcoef=rcorcoef,pvalue=mean(permcoef>rcorcoef)))
}

#' @rdname perm_kCCA
#' @export
perm_kCCA_par = function(x,y,sig=0.1,gama=0.1,ncomps=1,permNum = 500,kernel="rbfdot") {
  n = nrow(x)
  res =  kcca(x,y,kpar=list(sigma=sig),gamma = gama, ncomps = ncomps, kernel=kernel)
  rcorcoef = abs(kcor(res)[1])
  permcoef = sapply_pb(1:permNum, function(i) return(abs(kcor(kcca(x[sample(1:n,n),],y,kpar=list(sigma=sig),gamma = gama, ncomps = ncomps)))[1]))
  return(list(permcoef = permcoef,rcorcoef=rcorcoef,pvalue=mean(permcoef>rcorcoef)))
}




sapply_pb <- function(X, FUN, ...) {
  env <- environment()
  pb_Total <- length(X)
  counter <- 0
  pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)
  
  wrapper <- function(...){
    curVal <- get("counter", envir = env)
    assign("counter", curVal +1 ,envir=env)
    utils::setTxtProgressBar(get("pb", envir=env), curVal +1)
    FUN(...)
  }
  res <- sapply(X, wrapper, ...)
  close(pb)
  res
}


partition<-function(imageDat){
  position<-which(colnames(imageDat)=="V1")
  partImgdat<-list()
  for(i in 1:(length(position)-1)){
    partImgdat[[i]]<-imageDat[,position[i]:position[i+1]-1]
  }
  partImgdat[[length(position)]]<-imageDat[,utils::tail(position,n=1):dim(imageDat)[2]]
  return(partImgdat)
}


regionList<-function(regionData, regionCode, resolution="2mm"){
  
  cat("reading and manipulating the regionData...", "\n")
  largeMatrix<-oro.nifti::readNIfTI(regionData)
  longVector<-expand.grid(largeMatrix)
  
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


