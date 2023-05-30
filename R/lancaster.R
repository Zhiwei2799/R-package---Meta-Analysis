#' @title Lancaster's Method: a technique for meta analysis
#' @description Lancaster's method generalized Fisher's method by assignning different
#' weights using the additivity of chi-squared distribution.
#' \deqn{T = \sum\limits_{i=1}^{k} [X_{n_i}^{2}]^{-1}(1-p_i)},
#' where \deqn{[X_{n_i}^{2}]^{-1}} is the inverse cumulative chi-square distributon function with ni degrees of freedom
#' @param datalist input datalist
#' @param stat to show stat, by default stat=FALSE
#' @param p.adj  p values adjustment, by default p.adj = FALSE
#' @return p_value
#' @return statistic if stat=TRUE
#' @examples
#' dl <- sampleData(p=100)
#' lancaster(dl, stat = FALSE,p.adj = FALSE)
#' @examples
#' dl <- sampleData(p=100)
#' lancaster(dl, stat = FALSE,p.adj = 'fdr')
#' @export
lancaster <- function(datalist, stat=FALSE, p.adj=FALSE) {
  pmat <- p_matrix(datalist)
  sample_size <- matrix(NA,nrow=nrow(pmat),ncol=ncol(pmat))
  p_value <- c()
  statistic <- c()
  w <- matrix(NA,nrow=nrow(pmat),ncol=ncol(pmat))
  for (i in 1:ncol(pmat)){
    sample_size[,i] <- unname(sapply(datalist[[i]], function(df) length(complete.cases(df))))[-1]
  }
  for (i in 1:nrow(pmat)){
    statistic[i] <- sum(qchisq((1-pmat[i,]), df = sample_size[i,]))
    p_value[i] <- pchisq(statistic[i],df=sum(sample_size[i,]),lower.tail = F)
  }
  if (p.adj !=FALSE){
    p_value <- p.adjust(p_value,p.adj)
  }
  if (stat){
    return(list(statistic = statistic, p_value=p_value))
  }else{return(list(p_value = p_value))}
}





