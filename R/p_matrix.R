#' @title Computing a P values Matrix for a List of Dataframes
#' @description Within a dataset, biomarker measurements has to be all continuous or all categorical.
#' If each biomarker in a dataset has >2 groups, either ANOVA or Kruskal Wallis test will be chosen.
#'
#' If each biomarker in a dataset has =2 groups, either t-tests or Wilcoxon rank sum test will be chosen.
#'
#' If each biomarker in a dataset is categorical, either Fisher's exact test or chi-squared tests will be chosen.
#' @param datalist input datalist
#' @param stat_test statistical tests applied to each biomarker in data frames.Default is false
#' @return A matrix of P-values
#' @return A matrix of statistical tests if stat_test = TRUE
#' @examples
#' dl <-sampleData()
#' p_matrix(dl, stat_test=FALSE)
#' @export
p_matrix <- function(datalist, stat_test = FALSE) {
  if (!is.list(datalist)) {
    stop("Input is not a list")
  }

  if (!all(sapply(datalist, is.data.frame))) {
    stop("Input list must be a list of data frames")
  }

  if (!all(sapply(datalist, function(df) ncol(df) == ncol(datalist[[1]])))) {
    stop("Data frames must have the same number of columns")
  }
  if (any(sapply(datalist, function(df) any(is.na(df))))){
    warning("Data frame has missing values, and they will be ignored when computing p-values")
  }
  tests <- lapply(datalist, p_vector)
  p_value <- sapply(tests, function(x) x$p_value)
  if (stat_test){
    stat_test <- sapply(tests, function(x) x$stat_test)
    return(list(p_value=p_value, stat_test=stat_test))
  } else {
    return(p_value)
  }
}



#' @title Computing a P values Vector for a Dataframe
#' @param df Input dataframe
#' @param stat_test statistical tests applied to each biomarker in a data frame. Default is false
#' @return A vector of P-values
#' @return A vector of statistical tests if stat_test = TRUE
#' @examples
#' dl <-sampleData()
#' df <- dl[[1]]
#' p_vector(df, stat_test=FALSE)
#' @export
p_vector <-function(df,stat_test = TRUE) {
  n_groups <- length(unique(df$group))
  p_value <-c()
  tests <- c()
  if (is.numeric(df[,2])) {
    if (n_groups > 2) {
      for (i in 2:ncol(df)) {
        sp <- c()
        fp <- c()
        for (j in 1:n_groups) {
          sp[j] <- shapiro.test(df[df$group==j,i])$p.value
          if (j > 1) {
            fp[j-1] <- var.test(df[df$group==(j-1),i], df[df$group==j,i])$p.value
          }
        }
        if(all(sp > 0.05) & all(fp > 0.05)) {
          test <- aov(df[,i]~group, data=df,na.action = na.omit)
          p_value[i-1] <- summary(test)[[1]][["Pr(>F)"]][[1]]
          tests[i-1] <- "Anova"
        } else {
          p_value[i-1] <- kruskal.test(df[,i]~group,data=df,na.action = na.omit)$p.value
          tests[i-1] <- "Kruskal"
        }
      }
    } else {
      for (i in 2:ncol(df)) {
        sp1 <- shapiro.test(df[df$group==1,i])$p.value
        sp2 <- shapiro.test(df[df$group==2,i])$p.value
        fp <- var.test(df[df$group==1,i],df[df$group==2,i])$p.value
        if(sp1>0.05&sp2>0.05&fp>0.05){
          p_value[i-1] <- t.test(df[df$group==1,i],df[df$group==2,i],var.equal=TRUE,na.action = na.omit)$p.value
          tests[i-1] <- 'TtestEV'
        } else if(sp1>0.05&sp2>0.05&fp<=0.05){
          p_value[i-1] <- t.test(df[df$group==1,i],df[df$group==2,i],var.equal=FALSE,na.action = na.omit)$p.value
          tests[i-1] <- 'TtestUEV'
        } else {
          p_value[i-1] <- wilcox.test(df[df$group==1,i],df[df$group==2,i],na.action = na.omit)$p.value
          tests[i-1] <- 'Wilcox'
        }
      }
    }
  } else {
      for (i in 2:ncol(df)) {
        if (sum(table(df[, i], df$group) >= 5) / length(table(df[, i], df$group)) >= 0.8){
          chisq <- chisq.test(table(df[, i], df$group))
          tests[i-1] <- "Chisq"
          p_value[i-1] <- chisq$p.value
        }else{
          fisher <- fisher.test(table(df[, i], df$group))
          tests[i-1] <- "Fisher"
          p_value[i-1] <- fisher$p.value
        }
      }
  }
  if (stat_test) return(list(p_value = p_value, stat_test= tests))
  else return(list(p_value = p_value))
}
