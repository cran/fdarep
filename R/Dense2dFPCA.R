#' Two-Dimensional Functional Principal Component Analysis for dense repeated functional data.
#'
#'Note: The code works for dense functional data observed on a regular grid, missing values are allowed.
#'
#' @param X.age.year An n by (num.years x num.ages) input data matrix, such that the ith row of the matrix gives the observed values for the ith individual. The values in each row are sorted first by years (dimension 1) and then by ages (dimension 2) within each year.
#' @param n The sample size.
#' @param num.years Dimension 1
#' @param num.ages Dimension 2
#' @param fpca.op A list of options control parameters specified by \code{list(name=value)} for the two-dimesnional FPCA; check fdapace::FPCA for available control options and default settings.
#' @param pc.num A scalar denoting the maximum number of components to consider for the two-dimensional FPCA; default: chosen by FVE if NULL.
#'
#' @details The code works for dense functional data (with missing values), with density in both the direction of (age) dimension 2 and (year) dimension 1.
#'
#' @return A list containing the following fields:
#' \item{mu}{An num.ages by num.years matrix containing the bivariate mean function estimate.}
#' \item{pc.num}{A scalar denoting the selected number of components for the two-dimensional FPCA.}
#' \item{res.2D.FPCA}{A list containing the FPCA output for the fitted two-dimensional FPCA.}
#' \item{scores}{An n by pc.num matrix of the estimated scores, such that the ith row of the matrix comprises estimated scores for the ith individual.}
#' \item{eig}{An (num.years x num.ages) by pc.num matrix of the estimated product eigen functions. The estimated eigenfunctions in the otput eig are in the form of a vector rather than a matrix. For example, the first column in eig gives the first estimated eigenfunction such that \eqn{gamma(s,t) -> eig[ ( (s-1)*num.ages + t ), 1]} where LHS is the bivariate function in the usual form and RHS gives the corresponding element in the output vector. The rows are sorted first by years (dimension 1) and then by ages (dimension 2) within each year.}
#' \item{FVE}{A vector of length pc.num, indicating the fraction of total variance explained by each product function, with corresponding 'FVEthreshold'.}
#'
#' @examples
#' \donttest{
#' n <- 100 ### sample size
#' N <- 100
#' num.ages <- 20 ### dimension 2
#' num.years <- 15 ### dimension 1
#' dense_grid <- seq(0,1,length=N)
#' Lt <- list()
#' Ly <- list()
#' for (i in 1:n) {
#'   Lt[[i]] <- dense_grid ### dense time grid
#'   y_temp <- matrix(0,num.ages,num.years)
#'   for (s in 1:num.ages) {
#'     for (t in 1:num.years) {
#'       y_temp[s,t] <- y_temp[s,t]+cos(Lt[[i]][t])+rnorm(1,0,0.5)
#'     }
#'   }
#'   Ly[[i]] <- y_temp ### dense functional data
#' }
#' X.age.year <- matrix(0,n,num.years*num.ages)
#' for (i in 1:n) {
#'   X.age.year[i,] <- as.vector(Ly[[i]]) ### data matrix
#' }
#' res <- Dense2dFPCA(X.age.year, n , 15, 20, fpca.op=NULL,pc.num=2)
#' # Basic output
#' res$mu
#' res$pc.num
#' res$res.2D.FPCA
#' res$eig
#' res$FVE
#' res$pc.num
#' cumsum(res$FVE)
#' }
#'
#' @references
#' \itemize{
#' \item \cite{Chen, K., Delicado, P., & Müller, H. G. (2017). Modelling function-valued stochastic processes, with applications to fertility dynamics. Journal of the Royal Statistical Society Series B: Statistical Methodology, 79(1), 177-196.}
#' \item \cite{Chen, K., & Müller, H. G. (2012). Modeling repeated functional observations. Journal of the American Statistical Association, 107(500), 1599-1609.}
#' \item \cite{Hall, P.,  Müller, H.G. and Wang, J.L. (2006). Properties of principal component methods for functional and longitudinal data analysis. Annals of Statistics, 34(3), 1493-1517.}
#' \item \cite{Yao, F., Müller, H. G., & Wang, J. L. (2005). Functional data analysis for sparse longitudinal data. Journal of the American statistical association, 100(470), 577-590.}
#' }
#' @export

Dense2dFPCA <- function(X.age.year, n, num.years, num.ages, fpca.op=list(), pc.num=NULL){

  # Argument checking
  if(!("matrix" %in% class(X.age.year))){
    stop("The input X.age.year is not a matrix.")
  }
  if(!(class(n) %in% c("numeric","integer"))){
    stop("The input n is not an integer.")
  }
  if(as.integer(n)-n != 0 & n <= 0){
    stop("The input n is not a positive integer.")
  }

  if(!(class(num.years) %in% c("numeric","integer"))){
    stop("The input num.years is not an integer.")
  }
  if(as.integer(num.years)-num.years != 0 & num.years<=0){
    stop("The input num.years is not a positive integer.")
  }

  if(!(class(num.ages) %in% c("numeric","integer"))){
    stop("The input num.ages is not an integer.")
  }
  if(as.integer(num.ages)-num.ages != 0 & num.ages<=0){
    stop("The input num.ages is not a positive integer.")
  }

  if(!is.null(pc.num)){

  if(!(class(pc.num) %in% c("numeric","integer"))){
    stop("The input pc.num is not an integer.")
  }
  if(as.integer(pc.num)-pc.num != 0 | pc.num<=0){
    stop("The input pc.num is not a positive integer.")
  }

    }

  # n by (num.years*num.ages)
  if(!all(dim(X.age.year) == c(n, num.years*num.ages))){
    stop("The dimension of input X.age.year does not match with n by (num.years*num.ages)")
  }
  if (is.null(fpca.op)) {
    fpca.op <- list()
  }


  if( is.null(fpca.op$dataType)==FALSE){if(fpca.op$dataType == 'Sparse'){warning(sprintf("Warning: This function works for dense data only"))}}
  fpca.op$dataType = "Dense"

  # s-year; t-age

  # obs: Z(s,t) <-> Z( (s-1)*num.ages + t)
  # cov-fct: Cov(Z(s1,t1), Z(s2,t2)) <-> Cov( Z( (s1-1)*num.ages + t1), Z( (s2-1)*num.ages + t2) )
  # eigen-fct: gamma(s,t) <-> gamma( (s-1)*num.ages + t )

  if (sum(is.na(X.age.year)) == 0){
    tList <- lapply(1:n, function(x) 1:(num.ages*num.years))
    yList <- lapply(1:n, function(x) X.age.year[x,])
  } else {
    fpca.op$dataType <- 'DenseWithMV'
    tList <- vector('list', length  = n)
    yList <- vector('list', length = n)
    tvec <- 1:(num.ages*num.years)
    for (i in 1:n){
      ind <- which(!is.na(X.age.year[i,]))
      yList[[i]] <- X.age.year[i,ind]
      tList[[i]] <- tvec[ind]
    }

  }

  out.res <- fdapace::FPCA(yList, tList, fpca.op)
  if (is.null(pc.num)){
    pc.jk <- dim(out.res$xiEst)[2] ### pc.num value: default value based on FVE
    } else {
    pc.jk <- min(pc.num,dim(out.res$xiEst)[2]) ### pc.num value: based on user input
    if(pc.num > dim(out.res$xiEst)[2]){
      warning(sprintf("Desired number of components %d exceeds the maximum number of components %d returned by FPCA",
                      pc.num, dim(out.res$xiEst)[2]))
    }
  }

  pc.num <- pc.jk
  eig.Age.Year <- out.res$phi[,1:pc.num]

  mean.X <- apply(X.age.year, 2, mean, na.rm = TRUE)
  mean.X.mu <- t(matrix(mean.X,ncol=num.ages,nrow=num.years,byrow=T))

  FVE <- out.res$lambda[1:pc.num]/sum(out.res$lambda)
  scores <- out.res$xiEst[,1:pc.num]

  return(list(mu = mean.X.mu, pc.num = pc.num, res.2D.FPCA = out.res, scores=scores, eig = eig.Age.Year, FVE = FVE))
}


