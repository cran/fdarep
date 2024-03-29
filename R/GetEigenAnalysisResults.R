# phi: a nRegGrid * no_FVE
# The input smoothCov may be truncated.
# #' \cite{This function is copied from the fdapace package.}

GetEigenAnalysisResults <- function(smoothCov, regGrid, optns, muWork = NULL) {

  maxK <- optns$maxK
  FVEthreshold <- optns$FVEthreshold
  FVEfittedCov <- optns$FVEfittedCov
  verbose <- optns$verbose

  gridSize <- regGrid[2] - regGrid[1]
  numGrids <- nrow(smoothCov)

  eig <- eigen(smoothCov)

  positiveInd <- eig[['values']] >= 0
  if (sum(positiveInd) == 0) {
    stop('All eigenvalues are negative. The covariance estimate is incorrect.')
  }
  d <- eig[['values']][positiveInd]
  eigenV <- eig[['vectors']][, positiveInd, drop=FALSE]

  if (maxK < length(d)) {
    if (optns[['verbose']]) {
      message(sprintf("At most %d number of PC can be selected, thresholded by `maxK` = %d. \n", length(d), maxK))
    }

    d <- d[1:maxK]
    eigenV <- eigenV[, 1:maxK, drop=FALSE]
  }

  # thresholding for corresponding FVE option
  #(not before to avoid not being able to reach the FVEthreshold when pos eigenvalues > maxk)
  # i.e. default FVE 0.9999 outputs all components remained here.
  FVE <- cumsum(d) / sum(d) #* 100  # cumulative FVE for all available eigenvalues from fitted cov
  no_opt <- min(which(FVE >= FVEthreshold #* 100
                      )) # final number of component chosen based on FVE

  # normalization
  if (is.null(muWork)) {
    muWork = 1:dim(eigenV)[1]
  }

  phi <- apply(eigenV, 2, function(x) {
                    x <- x / sqrt(fdapace::trapzRcpp(regGrid, x^2))
                    if ( 0 <= sum(x*muWork) )
                      return(x)
                    else
                      return(-x)
  })
  lambda <- gridSize * d;

  if(!is.null(FVEfittedCov)){
    no_fittedCov <- min(which(FVE >= FVEfittedCov)) # number of components used to construct fittedCov
  }else{
    no_fittedCov <- dim(phi)[2]
  }
  fittedCovUser <- phi[,1:no_fittedCov, drop=FALSE] %*% diag(x=lambda[1:no_fittedCov], nrow = length(lambda[1:no_fittedCov])) %*% t(phi[,1:no_fittedCov, drop=FALSE])
  fittedCov <- phi %*% diag(x=lambda, nrow = length(lambda)) %*% t(phi)
  if(any(diag(fittedCovUser) ==0)){
    fittedCorrUser = NULL
  } else{
    fittedCorrUser <- (diag(1/sqrt(diag(fittedCovUser)))) %*% fittedCovUser %*% (diag(1/sqrt(diag(fittedCovUser))))
    diag(fittedCorrUser) = 1
  }
  return(list(lambda = lambda[1:no_opt], phi = phi[,1:no_opt, drop=FALSE],
              cumFVE = FVE, kChoosen=no_opt,
              fittedCov=fittedCov, fittedCovUser=fittedCovUser,
              fittedCorrUser=fittedCorrUser))
}
