#' @title Minimum p-value (minP): a technique for meta analysis
#' @description The minP method takes the minimum p-vlaue among the K studies as the test statistic.
#' It follows a beta distribution with degrees of freedom  \deqn{\alpha = 1} and \deqn{\beta = k} under the null hypothesis.
#' @param datalist input datalist
#' @param stat to show stat, by default stat=FALSE
#' @param p.adj  p values adjustment, by default p.adj = FALSE
#' @return p_value
#' @return statistic if stat=TRUE
#' @examples
#' dl <- sampleData(p=100)
#' minP(dl, stat = FALSE,p.adj = FALSE)
#' @examples
#' dl <- sampleData(p=100)
#' minP(dl, stat = FALSE,p.adj = 'fdr')
#' @export
minP <- function(datalist, stat = FALSE, p.adj=FALSE) {
  pmat <- p_matrix(datalist)
  min_p <- apply(pmat, 1, min)
  beta <- apply(pmat, 1, length)
  alpha <- rep(1, length(beta))
  p_value <- pbeta(min_p, alpha, beta,lower.tail = FALSE)
  if (p.adj !=FALSE){
    p_value <- p.adjust(p_value,p.adj)
  }
  if (stat){
    return(list(statistic = min_p,p_value=p_value))
  }else{return(list(p_value=p_value))}
}

