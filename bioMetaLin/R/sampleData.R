#' @title Sample data for R meta analysis package
#' @description A simple list of dataframes that is for simple testing and running.
#' @format A list of 3 data frames, each of which has 101 columns
#' It's a list of 3 data frames, each of which has 101 columns. For each data frame, first column indicating group membership,
#' columns 2 to p+1 representing the biomarker measurements for p biomarkers.
#'  \describe{
#'  \item{Dataframe1}{Data1 consists of 200 continuous data samples}
#'  \item{Dataframe2}{Data2 consists of 150 continous data samples}
#'  \item{Dataframe3}{Data3 consists of 300 categorical data samples}
#'  }
#' @param p  number of columns(by default p=100)
#' @return A list of 3 data frames.
#' @examples
#' sampleData()
#' sampleData(p=200)
#' @export
sampleData <- function(p=100){
  set.seed(123)
  data1 <- data.frame(group = sample(1:3, 200, replace = TRUE),
                      matrix(rnorm(p * 200), ncol = p))
  data2 <- data.frame(group = sample(1:2, 150, replace = TRUE),
                      matrix(rnorm(p * 150), ncol = p))
  data3 <- data.frame(group = sample(1:4, 300, replace = TRUE),
                      matrix(sample(0:2, size = 300 * p, replace = TRUE),
                             ncol = p))
  data3[, -1] <- lapply(data3[, -1], factor)
  colnames(data1) <- c("group", paste("gene", 1:p, sep = ""))
  colnames(data2) <- colnames(data3) <- colnames(data1)
  return(list(Data1 = data1, Data2 = data2, Data3 = data3))
}
