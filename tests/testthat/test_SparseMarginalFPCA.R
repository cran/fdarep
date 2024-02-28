require(testthat)
require(fdapace)

testthat::test_that('mismatched input dimension for Ly and Lt',{
  set.seed(1234)
  Ly <- lapply(1:20, function(i){matrix(rnorm(13*(i)), 13, i)})
  Lt <- lapply(1:20, function(i){1:(i)})
  sSup <- c(1:13)
  pc.j <- 2
  pc.k <- c(2,3)
  fpca.op1 <- NULL
  fpca.op2 <- NULL
  bw_mu_max <- 5.625000/2
  bw_mu_min <- NULL
  testthat::expect_error(SparseMarginalFPCA(sSup, Lt[[1:10]], Ly, fpca.op1, fpca.op2, pc.j, pc.k, bw_mu_min, bw_mu_max))
})




testthat::test_that('check output dimension',{
  set.seed(1234)
  Ly <- lapply(1:20, function(i){matrix(rnorm(13*(i)), 13, i)})
  Lt <- lapply(1:20, function(i){1:(i)})
  sSup <- c(1:13)
  pc.j <- 2
  pc.k <- c(2,3)
  fpca.op1 <- NULL
  fpca.op2 <- NULL
  bw_mu_max <- 5.625000/2
  bw_mu_min <- NULL
  res <- SparseMarginalFPCA(sSup, Lt, Ly, fpca.op1, fpca.op2, pc.j, pc.k, bw_mu_min, bw_mu_max)

  testthat::expect_equal(res$pc.j, 2)
  testthat::expect_equal(res$pc.k, c(2,3))
})


testthat::test_that('check output eigenvalues',{
  set.seed(1234)
  n <- 50
  J <- 2
  K <- 2
  N <- 50
  sGrid <- seq(0, 1, length=N)
  tGrid <- seq(0, 1, length=N)
  Psi <- fdapace::CreateBasis(J+1, pts=sGrid, type="fourier")[,2:(J+1)]
  Phi1 <- fdapace::CreateBasis(K, pts=tGrid, type='sin')
  Phi2 <- fdapace::CreateBasis(K, pts=tGrid, type='poly')
  Lambda1 <- c(5, 2)
  Lambda2 <- c(2, 1)

  Xlist <- lapply(1:n, function(o){
    xi <- rnorm(J * K, 0, 1)
    X11 <- xi[1] * sqrt(Lambda1[1]) * kronecker(Psi[,1], t(Phi1[,1]))
    X12 <- xi[2] * sqrt(Lambda1[2]) * kronecker(Psi[,1], t(Phi1[,2]))
    X21 <- xi[3] * sqrt(Lambda2[1]) * kronecker(Psi[,2], t(Phi1[,1]))
    X22 <- xi[4] * sqrt(Lambda2[2]) * kronecker(Psi[,2], t(Phi1[,2]))
    return(X11 + X12 + X21 + X22)
  })

  Lt <- lapply(1:n, function(i){
    1:(i)
  })
  Ly <- lapply(1:n, function(i){
    return(t(Xlist[[i]][1:(i),,drop = FALSE]))
  })

  res <- SparseMarginalFPCA(sGrid, Lt, Ly, fpca.op1=list(error=FALSE), fpca.op2=list(error=FALSE), pc.j=2, pc.k=c(2,2), bw_mu_min=5, bw_mu_max=50)

  testthat::expect_equal(c(5,2)/7, as.numeric(res$VarOrdered[[1]])/sum(res$VarOrdered[[1]]), tolerance = 0.3)
  testthat::expect_equal(c(2,1)/3, as.numeric(res$VarOrdered[[2]])/sum(res$VarOrdered[[2]]), tolerance = 0.3)
})




