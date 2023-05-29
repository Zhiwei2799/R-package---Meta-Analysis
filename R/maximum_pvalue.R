#' @title Maximum p-value (maxO): a technique for meta analysis
#' @description The maxP method takes maximum p-value as the test statistic.
#' It follows a beta distribution with degrees of freedom \deqn{\alpha = K} and \deqn{\beta = 1}
#' under the null hypothesis
#'
#' @param datalist input datalist
#' @param stat to show stat, by default stat=FALSE
#' @param p.adj  p values adjustment, by default p.adj = FALSE
#' @return p_value
#' @return statistic if stat=TRUE
#' @examples
#' dl <- sampleData(p=100)
#' maxP(dl, stat = FALSE,p.adj = FALSE)
#' @examples
#' dl <- sampleData(p=100)
#' maxP(dl, stat = FALSE,p.adj = 'fdr')
#' @export
maxP<- function(datalist, stat = FALSE, p.adj=FALSE){
  pmat <- p_matrix(datalist)
  max_p <- apply(pmat, 1 ,max)
  alpha <- apply(pmat, 1, length)
  beta <- rep(1, length(alpha))
  p_value <- pbeta(max_p, alpha, beta,lower.tail = FALSE)
  if (p.adj !=FALSE){
    p_value <- p.adjust(p_value,p.adj)
  }
  if (stat){
    return(list(statistic = max_p,p_value=p_value))
  }else{return(list(p_value=p_value))}
}



