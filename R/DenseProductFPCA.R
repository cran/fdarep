#' Product Functional Principal Component Analysis for dense repeated functional data.
#'
#'Note: The code works for dense functional data observed on a regular grid, missing values are allowed, written by Kehui Chen 10/09/2015, based on the original code by Pedro Delicado.
#'
#' @param X.age.year An n by (num.years x num.ages) input data matrix, such that the ith row of the matrix gives the observed values for the ith individual. The values in each row are sorted first by years (dimension 1) and then by ages (dimension 2) within each year.
#' @param n The sample size.
#' @param num.years Dimension 1
#' @param num.ages Dimension 2
#' @param fpca.op1 A list of options control parameters specified by \code{list(name=value)} for FPCA in direction of (age) dimension 2; check fdapace::FPCA for available control options and default settings.
#' @param fpca.op2 A list of options control parameters specified by \code{list(name=value)} for FPCA in direction of (year) dimension 1; check fdapace::FPCA for available control options and default settings.
#' @param pc.j A scalar denoting the maximum number of components to consider for FPCA in direction of (age) dimension 2; default: chosen by FVE if NULL.
#' @param pc.k A scalar denoting the maximum number of components to consider for FPCA in direction of (year) dimension 1; default: chosen by FVE if NULL.
#'
#' @details The code works for dense functional data (with missing values), with density in both the direction of (age) dimension 2 and (year) dimension 1.
#'
#' @return A list containing the following fields:
#' \item{Xest}{An n by (num.years x num.ages) estimated matrix, based on the fitted Product FPCA model. The ith row of the matrix gives the observed values for the ith individual. The values in each row are sorted first by years (dimension 1) and then by ages (dimension 2) within each year.}
#' \item{mu}{An num.ages by num.years matrix containing the bivariate mean function estimate.}
#' \item{pc.j}{A scalar denoting the selected number of components for FPCA in direction of (age) dimension 2.}
#' \item{pc.k}{A scalar denoting the selected number of components for FPCA in direction of (year) dimension 1.}
#' \item{scores}{An n by (pc.k x pc.j) matrix of the estimated scores, such that the ith row of the matrix comprises estimated scores \eqn{chi_{1,1}}, \eqn{chi_{1,2}},...\eqn{chi_{1,pc.k}}, \eqn{chi_{2,1}}, \eqn{chi_{2,2}}, ..., \eqn{chi_{2,pc.k}}, ..., \eqn{chi_{pc.j,1}}, \eqn{chi_{pc.j,2}}, ..., \eqn{chi_{pc.j,pc.k}} for the ith individual.}
#' \item{res.psi}{A list containing the FPCA output for FPCA in direction of (age) dimension 2.}
#' \item{res.phi}{A list containing the FPCA output for FPCA in direction of (year) dimension 1.}
#' \item{eig}{An (num.years x num.ages) by (pc.k x pc.j) matrix of the estimated product eigen functions. The rows are sorted first by years (dimension 1) and then by ages (dimension 2) within each year. The columns are sorted as in the scores.}
#' \item{psi}{An num.ages by pc.j matrix containing the estimated eigenfunctions from FPCA in direction of (age) dimension 2.}
#' \item{phi}{An num.years by pc.k matrix containing the estimated eigenfunctions from FPCA in direction of (year) dimension 1.}
#' \item{FVE}{A vector of length (pc.k x pc.j), indicating the fraction of total variance explained by each product function, with corresponding 'FVEthreshold'.}
#' \item{VarOrdered}{Variance explained by each term. The terms are ordered by \eqn{var(chi_{jk})}. One can select the best model by truncating at a desired level of FVE, and use names(VarOrdered) to see the corresponding model terms.}
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
#'  Lt <- lapply(1:n, function(i) dense_grid)
#'  Ly <- lapply(Lt, function(grid) {
#'  y_temp <- matrix(0, num.ages, num.years)
#'   lapply(1:num.ages, function(s) {
#'     lapply(1:num.years, function(t) {
#'       y_temp[s, t] <<- y_temp[s, t] + cos(grid[t]) + rnorm(1, 0, 0.5)
#'     })
#'   })
#'   y_temp
#'  })
#' X.age.year <- matrix(0,n,num.years*num.ages)
#' lapply(1:n, function(i) {
#'   X.age.year[i,] <<- as.vector(Ly[[i]])
#' })
#' res<-DenseProductFPCA(X.age.year, n, 15, 20, fpca.op1=NULL, fpca.op2=NULL, pc.j = NULL, pc.k = NULL)
#' # Basic output
#' res$Xest
#' res$mu
#' res$pc.j
#' res$pc.k
#' res$scores
#' res$res.psi
#' res$psi
#' res$FVE
#' res$VarOrdered
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

DenseProductFPCA <- function(X.age.year, n, num.years, num.ages, fpca.op1=list(), fpca.op2=list(), pc.j = NULL, pc.k = NULL){

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

  if(!is.null(pc.j)){
  if(!(class(pc.j) %in% c("numeric","integer"))){
    stop("The input pc.j is not an integer.")
  }
  if(as.integer(pc.j)-pc.j != 0 | pc.j<=0){
    stop("The input pc.j is not a positive integer.")
  }
    }

  if(!is.null(pc.k)){
  if(!(class(pc.k) %in% c("numeric","integer"))){
    stop("The input pc.k is not an integer.")
  }
  if(as.integer(pc.k)-pc.k != 0 | pc.k<=0){
    stop("The input pc.k is not a positive integer.")
  }
    }

  # n by (num.years*num.ages)
  if(!all(dim(X.age.year) == c(n, num.years*num.ages))){
    stop("The dimension of input X.age.year does not match with n by (num.years*num.ages)")
  }
  if (is.null(fpca.op1)) {
    fpca.op1 <- list()
  }
  if (is.null(fpca.op2)) {
    fpca.op2 <- list()
  }

  # step 1: centering
  mean.X <- apply(X.age.year, 2, mean, na.rm = TRUE)
  mean.X.mu <- t(matrix(mean.X,ncol=num.ages,nrow=num.years,byrow=T))
  mean.Xmat <- matrix(rep(mean.X,n), nrow = n, ncol = num.years*num.ages, byrow = TRUE)
  X.age.year.c <- X.age.year-mean.Xmat
  X.age.c <- matrix(t(X.age.year.c), nrow = n*num.years, ncol = num.ages, byrow = TRUE)

  # step 2: computing marginal eigenfunctions psi

  # use PACE functions so that it can handle missing values
  if (sum(is.na(X.age.c)) == 0){
    tList <- lapply(1:(n*num.years), function(x) 1:num.ages)
    yList <- lapply(1:(n*num.years), function(x) X.age.c[x,])
  } else {
    fpca.op1$dataType <- 'DenseWithMV'
    tList <- vector('list', length  = n*num.years)
    yList <- vector('list', length = n*num.years)
    tvec <- 1:num.ages
    for (i in 1:(n*num.years)){
      ind <- which(!is.na(X.age.c[i,]))
      yList[[i]] <- X.age.c[i,ind]
      tList[[i]] <- tvec[ind]
    }
  }

  res.psi <- fdapace::FPCA(yList, tList, optns = fpca.op1) ## FPCA for (age) dimension 2
  if (is.null(pc.j)){
    pc.s <- dim(res.psi$xiEst)[2] ###pc.j value: default value based on FVE
  }  else {
    pc.s <- min(pc.j,dim(res.psi$xiEst)[2]) ###pc.j value: based on user input
    if(pc.j > dim(res.psi$xiEst)[2]){
      warning(sprintf("Desired number of components %d exceeds the maximum number of components %d returned by FPCA",
                      pc.j, dim(res.psi$xiEst)[2]))
    }
  }

  psi <- res.psi$phi[, 1:pc.s]
  FVE.s <- res.psi$cumFVE[1:pc.s]
  xiEst.s <- res.psi$xiEst[, 1:pc.s]
  if(pc.s==1){xiEst.s <- matrix(res.psi$xiEst[, 1:pc.s],ncol=pc.s)}
  res.psi$selectK <- pc.s
  pc.j <- pc.s
  res.psi$xiEst <- xiEst.s

  # step 3: computing marginal eigenfunctions phi
  X.year.c <- matrix(X.age.year.c,  nrow = n*num.ages, ncol = num.years)
  if (sum(is.na(X.year.c)) == 0){
    tList <- lapply(1:(n*num.ages), function(x) 1:num.years)
    yList <- lapply(1:(n*num.ages), function(x) X.year.c[x,])
  } else {
    fpca.op2$dataType <- 'DenseWithMV'
    tList <- vector('list', length  = n*num.ages)
    yList <- vector('list', length = n*num.ages)
    tvec <- 1:num.years
    for (i in 1:(n*num.ages)){
      ind <- which(!is.na(X.year.c[i,]))
      yList[[i]] <- X.year.c[i,ind]
      tList[[i]] <- tvec[ind]
    }
  }

  res.phi <- fdapace::FPCA(yList, tList, optns = fpca.op2) ##FPCA for (year) dimension 1

  if (is.null(pc.k)){
    pc.t <- dim(res.phi$xiEst)[2] ###pc.k value: default value based on FVE
  }  else {
    pc.t <- min(pc.k,dim(res.phi$xiEst)[2]) ###pc.k value: based on user input
    if(pc.k > dim(res.phi$xiEst)[2]){
      warning(sprintf("Desired number of components %d exceeds the maximum number of components %d returned by FPCA",
                      pc.k, dim(res.phi$xiEst)[2]))
    }
  }

  phi <- res.phi$phi[, 1:pc.t]
  FVE.t <- res.phi$cumFVE[1:pc.t]
  xiEst.t <- res.phi$xiEst[, 1:pc.t]
  if(pc.t==1){xiEst.t <- matrix(res.phi$xiEst[, 1:pc.t],ncol=pc.t)}
  res.phi$selectK <- pc.t
  pc.k <- pc.t
  res.phi$xiEst <- xiEst.t

  eig.Age.Year <- matrix(0, nrow= num.years*num.ages, ncol= pc.j*pc.k)

  i <- 0
  for (j in 1:pc.j){
    for (k in 1:pc.k){
      i <- i +1
      eig.Age.Year[,i] <- as.numeric(outer(psi[,j],phi[, k],"*"))  ## Psi_j(s)*Phi_k(t) where every col represents one (j,k) combination

    }
  }

  lmfull.simu <- mFPCA.lm(t(X.age.year.c),eig.Age.Year)
  scores <- lmfull.simu$coefs
  Xest <- mean.Xmat + scores%*%t(eig.Age.Year) ##estimated X.age.year

  df.eig.Age.Year <- as.data.frame(eig.Age.Year)
  r <- dim(df.eig.Age.Year)[2]
  mFPCA.Age.Year.pctvar <- numeric(r)
  for (i in 1:r){
    mFPCA.Age.Year.pctvar[i] <-
      round(mFPCA.lm(t(X.age.year.c),df.eig.Age.Year,model=i)$perct.model,2)
  }

  best.model <- order(-mFPCA.Age.Year.pctvar)
  VarOrdered <- mFPCA.Age.Year.pctvar[best.model]
  namesV <- as.character(t(outer(paste("Psi",1:pc.j,sep=""),paste("Phi",1:pc.k,sep=""),paste,sep=".")))
  names(VarOrdered) <- namesV[best.model]
  return(list(Xest = Xest,mu = mean.X.mu, scores = scores, res.psi = res.psi, res.phi = res.phi,
              eig = eig.Age.Year, pc.j = pc.j, psi = psi, pc.k = pc.k, phi = phi, FVE = mFPCA.Age.Year.pctvar, VarOrdered = VarOrdered))
}

