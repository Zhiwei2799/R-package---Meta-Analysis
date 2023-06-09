#' @title Weighted Stouffer's Method: a technique for meta-analysis.
#' @description weight is the square root of the sample size, wi
#' \deqn{Z_{weighted~Stouffer} = \frac{\sum_{i=1}^{k} w_i Z_i}{\sqrt{\sum_{i=1}^{k} w_i^2}}}
#' @param datalist input datalist
#' @param stat to show stat, by default stat=FALSE
#' @param p.adj  p values adjustment, by default p.adj = FALSE
#' @return p_value
#' @return statistic if stat=TRUE
#' @examples
#' dl <- sampleData(p=100)
#' weighted_stouffer(dl, stat = FALSE,p.adj = FALSE)
#' @examples
#' dl <- sampleData(p=100)
#' weighted_stouffer(dl, stat = FALSE,p.adj = 'fdr')
#' @export
weighted_stouffer <- function(datalist, stat = FALSE, p.adj=FALSE) {
  pmat<-p_matrix(datalist)
  w <- matrix(0, nrow = nrow(pmat), ncol = ncol(pmat))
  statistic <- c()
  p_value <- c()
  sample_size <- matrix(NA,nrow=nrow(pmat),ncol=ncol(pmat))
  for (i in 1:ncol(pmat)){
    sample_size[,i] <- unname(sapply(datalist[[i]], function(df) length(complete.cases(df))))[-1]
  }
  w <- sqrt(sample_size)
  for (i in 1:nrow(pmat)){
    statistic[i] <- sum(w[i,]*qnorm(1-pmat[i,])) /sqrt(sum(w[i,]^2))
    p_value[i] <- pnorm(statistic[i], lower.tail=F)
  }
  if (p.adj !=FALSE){
    p_value <- p.adjust(p_value,p.adj)
  }
  if (stat){
    return(list(statistic = statistic,p_value=p_value))
  }else{return(list(p_value=p_value))}
}
