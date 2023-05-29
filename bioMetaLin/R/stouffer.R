#' @title Stoufferer's Method: a technique for meta-analysis.
#' @description Stouffer's method is used to combine the results from several independence tests
#' bearing upon the same overall hypothesis (H0). The method sums the inverse normal transformed
#' p-values.
#' \deqn{Z_{Stouffer} = \sum\limits_{i=1}^{k} Z_{i}/\sqrt{k}}
#'
#' @param datalist input datalist
#' @param stat to show stat, by default stat=FALSE
#' @param p.adj  p values adjustment, by default p.adj = FALSE
#' @return statistic
#' @return p_value
#' @examples
#' dl <- sampleData(p=100)
#' stouffer(dl, stat = FALSE,p.adj = FALSE)
#' @examples
#' dl <- sampleData(p=100)
#' stouffer(dl, stat = FALSE,p.adj = 'fdr')
#' @export
stouffer<- function(datalist,stat=FALSE,p.adj=FALSE){
  pmat<-p_matrix(datalist)
  statistic<-apply(pmat, 1, function(p_values){
    sum(qnorm(p_values, lower.tail = FALSE)) / sqrt(length(p_values))})
  p_value <- pnorm(statistic,lower.tail = FALSE)
  if(p.adj != FALSE){
      p_value <- p.adjust(p_value,p.adj)
  }
  if(stat){return(list(statistic = statistic, p_value = p_value))
    }else{return(list(p_value = p_value))}
}
