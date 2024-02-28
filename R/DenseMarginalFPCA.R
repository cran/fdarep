#' Marginal Functional Principal Component Analysis for dense repeated functional data.
#'
#'Note: The code works for dense functional data observed on a regular grid, missing values are allowed, written by Kehui Chen 10/09/2015, based on the original code by Pedro Delicado.
#'
#' @param X.age.year An n by (num.years*num.ages) input data matrix, such that the ith row of the matrix gives the observed values for the ith individual. The values in each row are sorted first by years (dimension 1) and then by ages (dimension 2) within each year.
#' @param n The sample size.
#' @param num.years Dimension 1
#' @param num.ages Dimension 2
#' @param fpca.op1 A list of options control parameters specified by \code{list(name=value)} for FPCA in direction of (age) dimension 2; check fdapace::FPCA for available control options and default settings.
#' @param fpca.op2 A list of options control parameters specified by \code{list(name=value)} for FPCA in direction of (year) dimension 1; check fdapace::FPCA for available control options and default settings.
#' @param pc.j A scalar denoting the maximum number of components to consider for marginal FPCA in direction of (age) dimension 2; default: chosen by FVE if NULL
#' @param pc.k A vector of length pc.j denoting the maximum number of components to consider for nested FPCA in direction of (year) dimension 1; default: chosen by FVE if NULL
#'
#' @details The code works for dense functional data (with missing values), with density in both the direction of (age) dimension 2 and (year) dimension 1.
#'
#' @return A list containing the following fields:
#' \item{Xest}{An n by (num.years x num.ages) estimated matrix, based on the fitted Marginal FPCA model. The ith row of the matrix gives the observed values for the ith individual. The values in each row are sorted first by years (dimension 1) and then by ages (dimension 2) within each year.}
#' \item{mu}{An num.ages by num.years matrix containing the bivariate mean function estimate.}
#' \item{pc.j}{A scalar denoting the selected number of components for FPCA in direction of (age) dimension 2.}
#' \item{pc.k}{A vector of length pc.j, denoting the selected number of components for FPCA in direction of (year) dimension 1.}
#' \item{scores}{An n by sum(pc.k) matrix of the estimated scores, such that the ith row of the matrix comprises estimated scores \eqn{chi_{1,1}}, \eqn{chi_{1,2}}, ..., \eqn{chi_{1,pc.k[1]}}, \eqn{chi_{2,1}}, \eqn{chi_{2,2}}, ..., \eqn{chi_{2,pc.k[2]}}, ..., \eqn{chi_{pc.j,1}}, \eqn{chi_{pc.j,2}}, ..., \eqn{chi_{pc.j,pc.k[j]}} for the ith individual.}
#' \item{res.psi}{A list containing the FPCA output for FPCA in direction of (age) dimension 2.}
#' \item{res.phi}{A list containing the FPCA output for FPCA in direction of (year) dimension 1.}
#' \item{eig}{An (num.years x num.ages) by sum(pc.k) matrix of the estimated product eigen functions. The rows are sorted first by years (dimension 1) and then by ages (dimension 2) within each year. The columns are sorted as in the scores.}
#' \item{psi}{An num.ages by pc.j matrix containing the estimated eigenfunctions for FPCA in direction of (age) dimension 2.}
#' \item{phi}{An num.years by sum(pc.k) matrix, containing the estimated eigenfunctions for FPCA in direction of (year) dimension 1.}
#' \item{FVE}{A vector of length sum(pc.k), indicating the fraction of total variance explained by each product function, with corresponding 'FVEthreshold'.}
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
#' res <- DenseMarginalFPCA(X.age.year, n, 15, 20,
#' fpca.op1=NULL, fpca.op2=NULL, pc.j = NULL, pc.k = NULL)
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
#' #Additional scores
#' fpca_psi <- res$res.psi
#' xi_i_j <- fpca_psi$xiEst
#' str(fpca_psi$xiEst)
#' #xi scores for 1st individual
#' xi_i_j[1:num.years,]
#' #xi scores for 2nd individual
#' xi_i_j[(num.years+1):(2*num.years),]
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


DenseMarginalFPCA <- function(X.age.year, n, num.years, num.ages, fpca.op1=list(), fpca.op2=list(), pc.j = NULL, pc.k = NULL){

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
  if(as.integer(pc.j)-pc.j != 0 & pc.j<=0){
    stop("The input pc.j is not a positive integer.")
  }
    }

  if(!is.null(pc.k)){
  # To check if pc.k is a vector of length pc.j
  if(length(pc.k)!=pc.j){
    stop("The input pc.k is not a vector of length pc.j! \n")
  }

  # to check if every element of pc.k is a positive integer
  if (!all(sapply(pc.k, function(x) {
    (as.integer(x)-x==0) && x > 0
  }))) {
    stop("One or more elements of pc.k are not positive integers. \n")
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


  Lt <- lapply(1:(n), function(x) 1:num.years) ## ith element of list is a vector from 1 to num.years
  Ly <- lapply(1:(n), function(x) matrix(as.vector(X.age.year.c[x,]), nrow = num.ages)) ## converts to matrix, in each row : each year and all ages and then next year

  # step 2: computing marginal eigenfunctions psi

  # use PACE functions so that it can handle missing values

  if (sum(is.na(X.age.c)) == 0){
    tList <- lapply(1:(n*num.years), function(x) 1:num.ages)
    yList <- lapply(1:(n*num.years), function(x) X.age.c[x,])
  } else {
    tList <- vector('list', length  = n*num.years) ## not clear!
    yList <- vector('list', length = n*num.years)
    tvec <- 1:num.ages
    for (i in 1:(n*num.years)){
      ind <- which(!is.na(X.age.c[i,]))
      yList[[i]] <- X.age.c[i,ind]
      tList[[i]] <- tvec[ind]
    }
  }
  res.psi <- fdapace::FPCA(yList, tList, optns = fpca.op1) ## FPCA along age dimension 2
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

  xiEst <- xiEst.s # xi_i[obs]_j[#eigen;col](t)[year]



  # step 3
  res.phi <- vector('list', length = pc.j)
  pc.k.2 <- numeric(pc.j)
  for (j in 1:pc.j){
    fpc <- matrix(xiEst[,j],nrow=n,ncol=num.years,byrow=TRUE) ## Xi_j(t)
    tList <- lapply(1:n, function(x) 1:num.years)
    yList <- lapply(1:n, function(x) fpc[x,])
    res <- fdapace::FPCA(yList, tList, optns = fpca.op2) ## Karhunen Loeve decomposition of Xi_j(t)
    res.phi[[j]] <- res
    eigenfunctions <- res$phi
    if (is.null(pc.k)){
      pc.k.2[j] <- dim(eigenfunctions)[2]
    } else {
      pc.k.2[j] <- min(pc.k[j],dim(eigenfunctions)[2])
      if(pc.k[j] > dim(eigenfunctions)[2]){
        warning(sprintf("Desired number of components %d exceeds the maximum number of components %d returned by FPCA",
                        pc.k[j], dim(eigenfunctions)[2]))
      }
      eigenfunctions <- eigenfunctions[,1:pc.k.2[j]]
    }

    if (j == 1){
      phi <- eigenfunctions
    } else {
      phi <- cbind(phi,eigenfunctions)
    }
  }

  pc.k <- pc.k.2
  eig.Age.Year <- matrix(0, nrow= num.years*num.ages, ncol= sum(pc.k)) ## Psi_j(s)*Phi_j_k(t) where every col represents one (j,k) combination

  i <- 0
  id <- 1:pc.j
  for (j in 1:pc.j){
    for (k in 1:pc.k[j]){
      i <- i +1
      #eig.Age.Year[,i] <- as.numeric(outer(psi[,j],phi[,(j-1)*pc.k[j]+k],"*"))
      eig.Age.Year[,i] <- as.numeric(outer(psi[,j],phi[,sum(pc.k[id<j])+k],"*")) ##Not very clear!

    }
  }

  lmfull.simu <- mFPCA.lm(t(X.age.year.c),eig.Age.Year)
  scores <- lmfull.simu$coefs
  Xest <- mean.Xmat+scores%*%t(eig.Age.Year) ##estimated X.age.year

  df.eig.Age.Year <- as.data.frame(eig.Age.Year)
  r <- dim(df.eig.Age.Year)[2]
  mFPCA.Age.Year.pctvar <- numeric(r)
  for (i in 1:r){
    mFPCA.Age.Year.pctvar[i] <-
      round(mFPCA.lm(t(X.age.year.c),df.eig.Age.Year,model=i)$perct.model,2) ## Not clear!
  }

  ###
  best.model <- order(-mFPCA.Age.Year.pctvar)
  VarOrdered <- mFPCA.Age.Year.pctvar[best.model]
  namesV <- vector('list', length = dim(eig.Age.Year)[2])
  i <- 0
  for (j in 1:pc.j){
    for (k in 1:pc.k[j]){
      i <- i+1
      namesV[i] <- as.character(paste('Psi',j, 'Phi',j,k, sep = "")) ## to indicate what are the columns of eig.Age.Year
    }
  }
  names(VarOrdered) <- namesV[best.model] ## to select components in a particular order
  return(list(Xest = Xest, mu = mean.X.mu, scores = scores, res.psi = res.psi, res.phi = res.phi, eig = eig.Age.Year,
              pc.j = pc.j, psi = psi, pc.k = pc.k, phi = phi, FVE = mFPCA.Age.Year.pctvar, VarOrdered = VarOrdered))
}


