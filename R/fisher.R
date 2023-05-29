#' @title Fisher's method: a technique for meta-analysis.
#' @description Fisher's method is used to combine the results from several independence tests
#' bearing upon the same overall hypothesis (H0). The method sums up the log-transformed p-values
#' obtained from individual studies.
#' \deqn{X_{fisher}^2 = \sum\limits_{i=1}^{k}\ln(p_i)}
#' @param datalist input datalist
#' @param stat to show stat, by default stat=FALSE
#' @param p.adj  p values adjustment, by default p.adj = FALSE
#' @return p_value
#' @return statistic if stat=TRUE
#' @examples
#' dl <- sampleData(p=100)
#' fisher(dl, stat = FALSE,p.adj = FALSE)
#' @examples
#' dl <- sampleData(p=100)
#' fisher(dl, stat = FALSE,p.adj = 'fdr')
#' @export
fisher <- function(datalist, stat = FALSE, p.adj=FALSE) {
  pmat<-p_matrix(datalist)
  statistic <- apply(pmat,1,FUN = function(p_values) -2 * sum(log(p_values)))
  df <- 2 * ncol(pmat)
  p_value <- pchisq(statistic, df,lower.tail = FALSE)
  if (p.adj !=FALSE){
  p_value <- p.adjust(p_value,p.adj)
  }
  if (stat){
    return(list(statistic = statistic,p_value=p_value))
  }else{return(list(p_value=p_value))}
}
