#' Product Functional Principal Component Analysis for sparse repeated functional data.
#'
#' Note: The code works for sparse functional data.
#'
#' @param sSup  A vector of length num.ages representing the common grid for every individual in the direction of (age) dimension 2.
#' @param Lt A list of n vectors containing the observation time points for every individual, such that the ith element of the list gives the num.years.i points in the direction of (year) dimension 1 at which the functional-valued stochastic process is observed for the ith individual.
#' @param Ly A list of n matrices containing the observed values for every individual, such that the ith element is an num.ages by num.years.i matrix of observed values for the i-th individual.
#' @param fpca.op1 A list of options control parameters specified by \code{list(name=value)} for FPCA in direction of (age) dimension 2; check fdapace::FPCA for available control options and default settings.
#' @param fpca.op2 A list of options control parameters specified by \code{list(name=value)}for FPCA in direction of (year) dimension 1; check fdapace::FPCA for available control options and default settings.
#' @param pc.j A scalar denoting the maximum number of components to consider for FPCA in direction of (age) dimension 2; default: chosen by FVE if NULL
#' @param pc.k A scalar denoting the maximum number of components to consider for FPCA in direction of (year) dimension 1; default: chosen by FVE if NULL.
#' @param bw_mu_min The minimum bandwidth value considered for bandwidth selection in mean function estimation, such that the final bandwidth chosen by 5-fold cross validation is above this minimum value; default:NULL, bandwidth chosen by 5-fold cross validation over the default range.
#' @param bw_mu_max The maximum bandwidth value considered for bandwidth selection in mean function estimation, such that the final bandwidth chosen by 5-fold cross validation is below this value; default:NULL, bandwidth chosen by 5-fold cross validation over the default range.
#'
#' @details This code works for sparse functional data, with the notion of sparsity defined as follows. Sparsity in the year direction (dimension 1) means that the years at which the data are observed for a country (or individual unit) are sparsely distributed. However for the ith county (or individual unit), if the data are available for a particular year (dimension 1), then it is available for all the ages (dimension 2) in sSup corresponding to that specific year. Thus along (age) dimension 2, data type is dense. The 'usergrid' control option in FPCA indicates whether to use observation grid for fitting, if false FPCA will use equidistant grid. logical - default:FALSE. Along (age) dimension 2, FPCA is done for only for sSup as observation grid. Depending on the choice of usergrid for 'fpca.op2', FPCA in (year) dimension 1 is either fitted on the observed (pooled) grid or on the internal regular grid of default length 51.
#'
#' @return A list containing the following fields:
#' \item{age.grid}{A vector of length num.ages, representing the grid used for fitting FPCA in the direction of (age) dimension 2, same as sSup.}
#' \item{year.grid}{A vector of length nWorkGrid, representing the grid used for fitting FPCA in the direction of (year) dimension 1.}
#' \item{mu}{An num.ages by nWorkGrid matrix containing the bivariate mean function estimate.}
#' \item{bwMu}{The selected bandwidth for mean function estimation.}
#' \item{pc.j}{A scalar denoting the selected number of components for FPCA in direction of (age) dimension 2.}
#' \item{pc.k}{A scalar denoting the selected number of components for FPCA in direction of (year) dimension 1.}
#' \item{res.psi}{A list containing the FPCA output for FPCA in direction of (age) dimension 2.}
#' \item{res.phi}{A list containing the FPCA output for FPCA in direction of (year) dimension 1.}
#' \item{scores}{A list of pc.j matrices containing the estimated scores, such that the jth element of the list is an n by pc.k matrix with its ith row comprising the estimated scores \eqn{chi_{j,1}},..., \eqn{chi_{j,pc.k[j]}} for the ith individual.}
#' \item{psi}{An num.ages by pc.j matrix containing the estimated eigenfunctions from FPCA in direction of (age) dimension 2.}
#' \item{phi}{An nWorkGrid by pc.k matrix, containing the estimated eigenfunctions from FPCA in direction of (year) dimension 1.}
#' \item{VarOrdered}{A list of pc.j vectors each of length pc.k, containing the variance explained by each term. The terms are ordered by \eqn{var(chi_{jk})}. One can select the best model by truncating at a desired level of FVE, and use names(VarOrdered) to see the corresponding model terms.}
#'
#' @examples
#' \donttest{
#' Ly <- lapply(1:20, function(i){matrix(rnorm(13*(i)), 13, i)})
#' Lt <- lapply(1:20, function(i){1:(i)})
#' sSup <- c(1:13)
#' pc.j <- 2
#' pc.k <- 3
#' fpca.op1 <- NULL
#' fpca.op2 <- NULL
#' bw_mu_max <- 5.625000/2
#' bw_mu_min <- NULL
#' res <- SparseProductFPCA(sSup, Lt, Ly, fpca.op1, fpca.op2, pc.j, pc.k, bw_mu_min, bw_mu_max)
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


SparseProductFPCA <- function(sSup, Lt, Ly, fpca.op1 = list(), fpca.op2 = list(), pc.j = NULL, pc.k = NULL, bw_mu_min=NULL, bw_mu_max=NULL){

  # Argument checking
  if(!is.list(Ly)){
    stop('Ly should be list. \n')
  }
  if(!is.list(Lt)){
    stop('Lt should be list. \n')
  }
  if(!is.numeric(sSup)){
    stop('sSup should be numeric vector. \n')
  }

  if( length(Lt) != length(Ly)){
    stop('Lt and Ly should have the same length. \n')
  }

  n_sSup = length(sSup)
  ni_y = sapply(Ly,function(x) dim(!is.na(x)))
  ni_tt = unlist(lapply(Lt,function(x) sum(!is.na(x))))
  if(!all(ni_y[2,] == ni_tt)){
    stop("SparseProductFPCA is aborted because the length of Lt and column numbers of Ly do not match!\n");
  }
  if(!all(ni_y[1,] == n_sSup)){
    stop("SparseProductFPCA is aborted because the length of sSup and row numbers of Ly do not match!\n");
  }
  if( !all(unlist(lapply(Ly,function(x) typeof(x) %in% c('integer', 'double') ) ) ) ){
    stop("SparseProductFPCA is aborted because 'Ly' members are not all of type double or integer! Try  \"lapply(y,function(x) typeof(x))\" to see the current types. \n");
  }
  if( !all(unlist(lapply(Lt,function(x) typeof(x) %in% c('integer', 'double'))) ) ){
    stop("SparseProductFPCA is aborted because 'Lt' members are not all of type double or integer! Try  \"lapply(t,function(x) typeof(x))\" to see the current types. \n");
  }

  if(any( unlist( lapply(Lt, function(x) length(x) != length(unique(x))))) ){
    stop("SparseProductFPCA is aborted because within-subject 'Lt' members have duplicated values.  Try  \"which( unlist( lapply(t, function(x) length(x) != length(unique(x)))))\" to see potentially problematic entries. \n");
  }
  if( any(sapply(Lt[seq_len(min(1001, length(Lt)))], is.unsorted, na.rm=TRUE)) ) {
    stop('Each vector in Lt should be in ascending order.')
  }
  if(min(unlist(Ly),na.rm=TRUE)==-Inf){
    stop('There are entries in Ly which are -Inf.')
  }
  if(max(unlist(Ly),na.rm=TRUE)==Inf){
    stop('There are entries in Ly which are Inf.')
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
  if(!(class(pc.k) %in% c("numeric","integer"))){
    stop("The input pc.k is not an integer.")
  }
  if(as.integer(pc.k)-pc.k != 0 & pc.k<=0){
    stop("The input pc.k is not a positive integer.")
  }
   }

  if (is.null(fpca.op1)) {
    fpca.op1 <- list()
  }
  if (is.null(fpca.op2)) {
    fpca.op2 <- list()
  }

  if(!is.null(bw_mu_min)){
  if(!typeof(bw_mu_min) %in% c('integer', 'double')){
    stop("SparseProductFPCA is aborted because bw_mu_minis not the type double or integer! \n");
  }
  }

  if(!is.null(bw_mu_max)){
  if(!typeof(bw_mu_max) %in% c('integer', 'double')){
    stop("SparseProductFPCA is aborted because bw_mu_minis not the type double or integer! \n");
  }
   }

  n <- length(Lt)
  s_t.n <- do.call(cbind, Ly)
  t.n <- do.call(c, Lt) ### Lt[[i]]: years for ith individual
  t.n_cum <- c(0,cumsum(sapply(Lt, function(x) length(x)))) ## n1, n1+n2, n1+n2+n3, ... ###Index to extract individual data
  yList1 <- lapply(1:length(t.n), function(i) s_t.n[,i]) ### all ASFRs for an individual at a particular year for all ages
  tList <- lapply(1:length(t.n), function(i) sSup) ### same ages for all individuals
  ### Sparsity assumed only in year dimension, for every available year of an individual: ASFR for all ages in sSup available


  # step 1: centering and computing the mean function mu(s,t)
  obsGrid <- sort(unique(t.n)) ## To get the pooled range of years observed
  mu_s_t <-  matrix(0,ncol=length(obsGrid),nrow=length(sSup))
  for(i in 1:length(sSup)){ ### ith row corresponding to jth age in sSup
    Ltt <- list()
    Lyy <- list()
    Ltt <- lapply(1:n, function(j) {Lt[[j]]})
    Lyy <- lapply(1:n, function(j) {Ly[[j]][i,]})

    # get the bandwidth using CV to estimnate the mean function (below)
    suppressWarnings({
      bw_mu <- CVLwls1D(Lyy, Ltt, kernel= "epan", npoly=1, nder=0, dataType= "sparse", kFolds = 5,
                                 useBW1SE = FALSE)
    })

    bw_mu_fin <- bw_mu
    bw_mu_final <- bw_mu

    if(is.null(bw_mu_min)==FALSE)
    {
      if(bw_mu_min > 0){bw_mu_fin <- max(bw_mu,bw_mu_min)}
      else{warning(sprintf("The minimum bandwidth should be greater than 0"))}
    }

    if(is.null(bw_mu_max)==FALSE)
    {
      if(bw_mu_max>0){bw_mu_final <- min(bw_mu_fin,bw_mu_max)}
      else{warning(sprintf("The maximum bandwidth should be greater than 0"))}
    }

    # Get the mean function using the bandwith estimated above:
    xin <- numeric()
    yin <- numeric()
    win <- numeric()
    mu <- numeric()
    xin <- unlist(Ltt)
    yin <- unlist(Lyy)[order(xin)]
    xin <- sort(xin)
    win <- rep(1, length(xin))
    mu <- fdapace::Lwls1D(bw_mu_final, kernel_type = "epan", npoly = 1, nder = 0,
                         xin = xin, yin= yin, xout = obsGrid, win = win)
    mu_s_t[i,] <- mu
  }
  #message(sprintf("Bandwidth for mean estimation: The bandwidth chosen for mean function estimation is %f",bw_mu_final))

  ##### centered functional values #####
  ind <- sapply(t.n, function(t){which(obsGrid == t)})
  yList <- lapply(1:length(t.n), function(j) yList1[[j]]-mu_s_t[,ind[j]]) ### all ASFRs for an individual at a particular year for all ages

  # step 2: computing marginal eigenfunctions psi
  # fpca.op1 = NULL
  if( is.null(fpca.op1$usergrid)==FALSE){if(fpca.op1$usergrid == FALSE){warning(sprintf("Warning: Along dimension 2, FPCA is done for only for sSup as observation grid"))}}
  if( is.null(fpca.op1$dataType)==FALSE){if(fpca.op1$dataType == 'Sparse'){warning(sprintf("Warning: Along dimension 2, data type is dense"))}}
  if( is.null(fpca.op1$userMu)==FALSE){warning(sprintf("Warning: Along dimension 2, FPCA does not accept user defined mean since the function internally decentralizes using bivariate mean function"))}
  fpca.op1$usergrid <- TRUE       # FALSE
  fpca.op1$userMu <- list("t"= sSup, "mu"=rep(0, length(sSup))) # if NULL
  # fpca.op1$dataType <- 'Dense'    # Dense dimension age
  #fpca.op1$FVEthreshold <- 0.9999 # 0.99
  #fpca.op1$error <- FALSE         # TRUE
  fpca.op1 <- fdapace::SetOptions(yList, tList, fpca.op1);
  res.s <- fdapace::FPCA(yList, tList, optns = fpca.op1) ##FPCA in age dimension
  if (is.null(pc.j)){
    pc.s <- dim(res.s$xiEst)[2] ###pc.j value: default value based on FVE
  } else {
    pc.s <- min(pc.j,dim(res.s$xiEst)[2]) ###pc.j value: based on user input
    if(pc.j > dim(res.s$xiEst)[2]){
      warning(sprintf("Desired number of components %d exceeds the maximum number of components %d returned by FPCA",
                      pc.j, dim(res.s$xiEst)[2]))
    }
  }
  psi <- res.s$phi[, 1:pc.s]
  FVE.s <- res.s$cumFVE[1:pc.s]
  xiEst.s <- res.s$xiEst[, 1:pc.s]
  if(pc.s==1){xiEst.s <- matrix(res.s$xiEst[, 1:pc.s],ncol=pc.s)}
  res.s$selectK <- pc.s
  res.s$xiEst <- xiEst.s
  res.s$phi <- psi

  # step 3: computing marginal eigenfunctions phi
  xiEst <- vector('list', length = pc.s)
  Lambda <- vector('list', length = pc.s)
  #if( is.null(fpca.op2$dataType)==FALSE){if(fpca.op2$dataType == 'Dense'){warning(sprintf("Along dimension 1, data type is sparse"))}}
  if( is.null(fpca.op2$error)==FALSE){if(fpca.op2$error == 'TRUE'){warning(sprintf("Along dimension 1, we assume measurement error to be FALSE in order to avoid unstable estimates for sparse data. Hence changed to FALSE"))}}
  if(is.null(fpca.op2$error)== TRUE){
    fpca.op2$error <- FALSE
    message(sprintf("Along dimension 1, we assume measurement error to be FALSE in order to avoid unstable estimates for sparse data. Hence error = FALSE is used"))}
  if(is.null(fpca.op2$userMu)==FALSE){warning(sprintf("Along dimension 1, FPCA does not accept user defined mean since the function internally decentralizes using bivariate mean function"))}
  #fpca.op2 <- NULL
  #fpca.op2$dataType <- 'Sparse'       # Sparse dimension year
  #fpca.op2$error <- FALSE             # TRUE
  #fpca.op2$FVEthreshold <- 0.9999     # 0.99
  #fpca.op2$usergrid <- FALSE          # FALSE
  #fpca.op2$methodBwCov <- "GCV"
  # fpca.op2$maxK <- min(length(obsGrid), n)
  #fpca.op2$verbose <- FALSE
  fpca.op2$error <- FALSE
  ran <- range(unlist(Lt))
  # if(fpca.op2$userBwCov == 'NULL'){
  #   fpca.op2$userBwCov <- (ran[2]-ran[1])/10
  #   warning(sprintf("Warning: If left NULL, "))
  #   }
  yList <- list()
  tList <- list()
  for (i in 1:length(Ly)) {
    for (s in 1:length(sSup)) {
      yList <- c(yList, list(Ly[[i]][s,]))
      tList <- c(tList, list(Lt[[i]]))
    }
  }
  #obsGrid1=seq(floor(range(unlist(tList))[1]),ceiling(range(unlist(tList))[2]),length=1000)
  obsGrid <- sort(unique(unlist(Lt)))
  fpca.op2$userMu <- list("t"= obsGrid, "mu"=rep(0, length(obsGrid))) # NULL
  #fpca.op2$userMu <- list("t"= obsGrid1, "mu"=rep(0, length(obsGrid1))) # NULL
  fpca.op2 <- fdapace::SetOptions(yList, tList, fpca.op2);


  ## get raw covariances for every pair of time points
  n <- length(Lt)
  yin <- unlist(lapply(1:n, function(i){
    L <- length(Lt[[i]])
    if(L==1){
      Ly[[i]] <- matrix(Ly[[i]], nrow = length(sSup)) ## For points with 1 observation at a single time point, it contributes only to variance
    }
    ind_of_t <- which(obsGrid %in% Lt[[i]])
    unlist(lapply(1:L, function(j1){
      sapply(1:L, function(j2){
        y <- (Ly[[i]][,j1] - mu_s_t[,ind_of_t[j1]])*(Ly[[i]][,j2] - mu_s_t[,ind_of_t[j2]])
        return(pracma::trapz(sSup, y)) ##intregration of X_c(j1,s)*X_c(j2,s) ds
      })
    }))
  }))
  ## Construction of pairs of time points in the range of the time domain
  xin <- t(do.call(cbind, lapply(1:n, function(i){
    L <- length(Lt[[i]])
    do.call(cbind,lapply(1:L, function(j1){
      sapply(1:L, function(j2){
        return(c(Lt[[i]][j1], Lt[[i]][j2]))
      })
    }))
  })))
  ran <- range(unlist(Lt))
  if(fpca.op2$usergrid == TRUE){
    workGrid <- sort(unique(unlist(Lt)))
  } else {
    if(fpca.op2$dataType == 'Dense' || fpca.op2$dataType == 'DenseWithMV'){
      tt <- unlist(Lt)
      nRegGrid <- length(unique(signif(tt[!is.na(tt)],6)));
    } else { # for Sparse and p>>n
      nRegGrid <- 51;
    }
    workGrid <- seq(ran[1],ran[2],length.out=nRegGrid)
  }
  cov.t <- fdapace::Lwls2D(bw = (ran[2]-ran[1])/10, kern="epan", xin, yin, xout1 = workGrid, xout2 = workGrid)
  cov.t <- (cov.t + t(cov.t))/2  # C(t1, t2) = int cov(X(s,t1), X(s,t2)) ds
  eigObj.t <- GetEigenAnalysisResults(smoothCov = cov.t, regGrid = workGrid, fpca.op2, muWork = rep(0, 51))
  phi <- eigObj.t$phi ## Getting eigen functions from the eigen decomposition of marginal covariance function G_T

  data <- cbind(workGrid, phi)
  # write.table(data,file="adniPhiProduct.txt", append=FALSE, col.names=FALSE, row.names=FALSE)

  tList <- lapply(1:n, function(i) t.n[c((t.n_cum[i]+1):t.n_cum[i+1])]) ##Index to extract individual years but same as Lt
  #obsGrid1 <- seq(floor(range(unlist(tList))[1]),ceiling(range(unlist(tList))[2]),length=1000)
  #obsGrid1 <- seq(0,1,length=1000)
  obsGrid <- sort(unique(unlist(Lt)))
  fpca.op2$userMu <- list("t"= obsGrid, "mu"=rep(0, length(obsGrid))) # NULL
  #fpca.op2$userMu <- list("t"= obsGrid, "mu"=rep(0, length(obsGrid1))) # NULL
  for (j in 1:pc.s){
    yList <- lapply(1:n, function(i) xiEst.s[c((t.n_cum[i]+1): t.n_cum[i+1]),j]) ### Xi_j(t) = integration X_c(s,t) x Psi_j(s) ds
    fpca.op2 <- fdapace::SetOptions(yList, tList, fpca.op2);
    res <- fdapace::FPCA(yList, tList, optns = fpca.op2) ##FPCA in year dimension 2
    if (is.null(pc.k)){
      pc.t <- min(dim(res$xiEst)[2],dim(phi)[[2]]) ###pc.k value: default value based on FVE
    } else {
      pc.t <- min(pc.k,dim(res$xiEst)[2],dim(phi)[[2]]) ###pc.k value: based on user input
      if(pc.k > max(dim(res$xiEst)[2],dim(phi)[[2]])){
        warning(sprintf("Warning: Desired number of components %d exceeds the maximum number of components %d returned by FPCA",
                        pc.t, dim(res$xiEst)[2]))
      }
    }
    fittedCov <- res$fittedCov
    Phi <- matrix(phi[,1:pc.t], ncol = pc.t)

    Lambda[[j]] <- sapply(1:pc.t, function(k){
      v <- sapply(1:length(workGrid), function(l){
        pracma::trapz(workGrid, fittedCov[,l]*phi[,k]) # int cov(Xi_j(t1), Xi_j(t)) * phi_k(t) dt = lambda_k * phi_k(t1)
      })
      pracma::trapz(workGrid, v*phi[,k]) # int lambda_k * phi_k(t1) * phi_k(t1) dt1 = lambda_k

    })

    scoresObj <- GetCEScores(y=yList, t=tList, fpca.op2, mu=rep(0, length(workGrid)), obsGrid=workGrid, fittedCov=fittedCov,
                                      lambda=Lambda[[j]], phi=Phi, sigma2 = NULL)
    xiEst[[j]] <- t(do.call(cbind, scoresObj[1,]))
  }


  mu_s_t <- matrix(0,ncol=length(workGrid),nrow=length(sSup))
  for(i in 1:length(sSup)){ ### ith row corresponding to jth age in sSup
    Ltt <- list()
    Lyy <- list()
    Ltt <- lapply(1:n, function(j) {Lt[[j]]})
    Lyy <- lapply(1:n, function(j) {Ly[[j]][i,]})

    # get the bandwidth using CV to estimnate the mean function (below)
    suppressWarnings({
      bw_mu <- CVLwls1D(Lyy, Ltt, kernel= "epan", npoly=1, nder=0, dataType= "sparse", kFolds = 5,
                                 useBW1SE = FALSE)
    })

    bw_mu_fin <- bw_mu
    bw_mu_final <- bw_mu

    if(is.null(bw_mu_min)==FALSE)
    {
      if(bw_mu_min > 0){bw_mu_fin=max(bw_mu,bw_mu_min)}
      else{warning(sprintf("The minimum bandwidth should be greater than 0"))}
    }

    if(is.null(bw_mu_max)==FALSE)
    {
      if(bw_mu_max>0){bw_mu_final=min(bw_mu_fin,bw_mu_max)}
      else{warning(sprintf("The maximum bandwidth should be greater than 0"))}
    }

    # Get the mean function using the bandwith estimated above on the workGrid
    xin <- numeric()
    yin <- numeric()
    win <- numeric()
    mu <- numeric()
    xin <- unlist(Ltt)
    yin <- unlist(Lyy)[order(xin)]
    xin <- sort(xin)
    win <- rep(1, length(xin))
    mu <- fdapace::Lwls1D(bw_mu_final, kernel_type = "epan", npoly = 1, nder = 0,
                         xin = xin, yin= yin, xout = workGrid, win = win)
    mu_s_t[i,] <- mu
  }

  ## calculate FVE
  TotalVariance <- sum(unlist(Lambda))
  return(list("age.grid" = sSup, "year.grid" = workGrid,
              "mu" = mu_s_t,"bwMu" = bw_mu_final,
              "pc.j"=pc.s,"pc.k"=pc.t,
              "res.psi" = res.s, "res.phi" = res,
              "psi" = psi, "phi"=Phi,
              "scores" = xiEst, "VarOrdered" = Lambda
  ))
}
