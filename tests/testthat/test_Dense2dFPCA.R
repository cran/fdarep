require(testthat)
require(pracma)
require(fdapace)

testthat::test_that('mismatched input dimension for sample size',{
  set.seed(1234)
  n <- 100 ### sample size
  N <- 100
  num.ages <- 20 ### dimension 2
  num.years <- 15 ### dimension 1
  dense_grid <- seq(0,1,length=N)
  Lt <- list()
  Ly <- list()
  for (i in 1:n) {
    Lt[[i]] <- dense_grid ### dense time grid
    y_temp <- matrix(0,num.ages,num.years)
    for (s in 1:num.ages) {
      for (t in 1:num.years) {
        y_temp[s,t] <- y_temp[s,t]+cos(Lt[[i]][t])+rnorm(1,0,0.5)
      }
    }
    Ly[[i]] <- y_temp ### dense functional data
  }
  X.age.year <- matrix(0,n,num.years*num.ages)
  for (i in 1:n) {
    X.age.year[i,] <- as.vector(Ly[[i]]) ### data matrix
  }
  testthat::expect_error(Dense2dFPCA(X.age.year, 200 , 15, 20, fpca.op=NULL,pc.num=2))
})


testthat::test_that('mismatched input dimension for dimension1 and dimension2',{
  set.seed(1234)
  n <- 100 ### sample size
  N <- 100
  num.ages <- 20 ### dimension 2
  num.years <- 15 ### dimension 1
  dense_grid <- seq(0,1,length=N)
  Lt <- list()
  Ly <- list()
  for (i in 1:n) {
    Lt[[i]] <- dense_grid ### dense time grid
    y_temp <- matrix(0,num.ages,num.years)
    for (s in 1:num.ages) {
      for (t in 1:num.years) {
        y_temp[s,t] <- y_temp[s,t]+cos(Lt[[i]][t])+rnorm(1,0,0.5)
      }
    }
    Ly[[i]] <- y_temp ### dense functional data
  }
  X.age.year <- matrix(0,n,num.years*num.ages)
  for (i in 1:n) {
    X.age.year[i,] <- as.vector(Ly[[i]]) ### data matrix
  }
  testthat::expect_error(Dense2dFPCA(X.age.year, n, 20, 20, fpca.op=NULL, pc.num=2))
  testthat::expect_error(Dense2dFPCA(X.age.year, n, 15, 15, fpca.op=NULL, pc.num=2))
})



testthat::test_that('check output dimension',{
  set.seed(1234)
  n <- 100 ### sample size
  N <- 100
  num.ages <- 20 ### dimension 2
  num.years <- 15 ### dimension 1
  dense_grid <- seq(0,1,length=N)
  Lt <- list()
  Ly <- list()
  for (i in 1:n) {
    Lt[[i]] <- dense_grid ### dense time grid
    y_temp <- matrix(0,num.ages,num.years)
    for (s in 1:num.ages) {
      for (t in 1:num.years) {
        y_temp[s,t] <- y_temp[s,t]+cos(Lt[[i]][t])+rnorm(1,0,0.5)
      }
    }
    Ly[[i]] <- y_temp ### dense functional data
  }
  X.age.year <- matrix(0,n,num.years*num.ages)
  for (i in 1:n) {
    X.age.year[i,] <- as.vector(Ly[[i]]) ### data matrix
  }
  res <- Dense2dFPCA(X.age.year, 100, 15, 20, fpca.op=NULL, pc.num=2)
  testthat::expect_equal(res$pc.num, 2)
})


testthat::test_that('check output eigenvalues',{
  set.seed(1234)
  n <- 100
  J <- 2
  K <- 2
  N <- 20
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

  Xmat <- t(sapply(Xlist, function(mat) as.vector(t(mat))))
  res <- Dense2dFPCA(Xmat, n, N, N, pc.num = 4)
  testthat::expect_equal(c(5,2,2,1)/10, as.numeric(res$FVE), tolerance = 0.2)
})




