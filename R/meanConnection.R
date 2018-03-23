#' Calculate percentage of connection in all pairwise brain regions.
#'
#' This function can create a list of significant (threshold is defined by user) region pairs.
#' @author Xubo Yue, Chia-Wei Hsu (tester), Jian Kang (maintainer)
#' @param path the path where csv files located
#' @param threshold the threshold for significance of percentage of 
#' connection (if percentage exceeds threhold, then
#' the region pair is significant). Typically, it can be 15-30\%.
#' @return the object containing significant regions.
#' @details you need to specify the path where csv files 
#' (containing KCCA information)locoated. This function
#' will read all csv files listed in the path. 
#' @export
#' @examples
#' #summary_data <- meanConnection(path = "UM_1_result", threshold=0.25)
#' #multipleRegion_plot(summary_data, significance=NA)

meanConnection<-function(path = getwd(), threshold=0.2){
  
  currentPath <- getwd()
  setwd(path)
  csvFiles <- list.files(pattern='*\\.csv')
  data_file1 <- read.csv(csvFiles[1]) 
  
  data_file1 <- read.csv(csvFiles[1]) 
  
  if(length(csvFiles)>1){
    for(i in 2:length(csvFiles)){
      data_file2 <- read.csv(csvFiles[i])
      data_file1 <- rbind(data_file1,data_file2)
      data_file1[] <- lapply(data_file1, function(x) type.convert(as.character(x)))
      data_file1 <- aggregate(. ~ index1+index2+region1+region2, data_file1, sum)
    }
  }
    
  data_file1$percent <- data_file1$significant/data_file1$count
  data_file2N <- data_file1[data_file1$percent>=threshold,]
  result_noP <- NULL
  result_noP[[1]] <- cbind(data_file2N$index1, data_file2N$index2)
  result_noP[[2]] <- result_noP[[3]] <- NULL
  result_noP[[4]] <- list(cbind(data_file2N$region1, data_file2N$region2))
  final_noP <- NULL
  final_noP[[1]] <- result_noP
  
  setwd(currentPath)
  return(final_noP)
}





