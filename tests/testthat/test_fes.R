context("Testing fes function")
test_that("Testing that fes calculates correctly sum of points", {
  # fes function 2D
  myfes<-fes(acealanme, imax=2000)
  summyfes<-sum(myfes$fes)
  expect_equal(object=summyfes, expected=-1127743, tolerance=1, scale=1)
  
  # fes function 2D
  myfes<-fes2(acealanme, imax=2000)
  summyfes<-sum(myfes$fes)
  expect_equal(object=summyfes, expected=-1127743, tolerance=1, scale=1)
  
  # fes function 1D
  myfes<-fes(acealanme1d, imax=2000)
  summyfes<-sum(myfes$fes)
  expect_equal(object=summyfes, expected=-14790, tolerance=1, scale=1)
  
  # fes function 1D
  myfes<-fes2(acealanme1d, imax=2000)
  summyfes<-sum(myfes$fes)
  expect_equal(object=summyfes, expected=-14790, tolerance=1, scale=1)

  #+-*/
  expect_equal(object=summyfes+summyfes, expected=-29581, tolerance=1, scale=1)
  expect_equal(object=summyfes-summyfes, expected=0, tolerance=1, scale=1)
  expect_equal(object=2*summyfes, expected=-29581, tolerance=1, scale=1)
  expect_equal(object=summyfes/0.5, expected=-29581, tolerance=1, scale=1)

  # fes2d21d
  myfes<-fes2d21d(acealanme, imax=2000)
  summyfes<-sum(myfes$fes)
  expect_equal(object=summyfes, expected=-9091, tolerance=1, scale=1)
  
  # fes function 2D with per=c(T,F)
  acealanme2 <- acealanme
  acealanme2$per <- c(T,F)
  myfes<-fes(acealanme2, imax=2000)
  summyfes<-sum(myfes$fes)
  expect_equal(object=summyfes, expected=-1016915, tolerance=1, scale=1)
  
  # fes function 2D with per=c(F,T)
  acealanme2 <- acealanme
  acealanme2$per <- c(F,T)
  myfes<-fes(acealanme2, imax=2000)
  summyfes<-sum(myfes$fes)
  expect_equal(object=summyfes, expected=-1019403, tolerance=1, scale=1)

  # fes function 2D with per=c(F,F)
  acealanme2 <- acealanme
  acealanme2$per <- c(F,F)
  myfes<-fes(acealanme2, imax=2000)
  summyfes<-sum(myfes$fes)
  expect_equal(object=summyfes, expected=-919238, tolerance=1, scale=1)

  # fes2 function 2D with per=c(T,F)
  acealanme2 <- acealanme
  acealanme2$per <- c(T,F)
  myfes<-fes2(acealanme2, imax=2000)
  summyfes<-sum(myfes$fes)
  expect_equal(object=summyfes, expected=-1016897, tolerance=1, scale=1)
  
  # fes2 function 2D with per=c(F,T)
  acealanme2 <- acealanme
  acealanme2$per <- c(F,T)
  myfes<-fes2(acealanme2, imax=2000)
  summyfes<-sum(myfes$fes)
  expect_equal(object=summyfes, expected=-1019409, tolerance=1, scale=1)

  # fes2 function 2D with per=c(F,F)
  acealanme2 <- acealanme
  acealanme2$per <- c(F,F)
  myfes<-fes2(acealanme2, imax=2000)
  summyfes<-sum(myfes$fes)
  expect_equal(object=summyfes, expected=-919226, tolerance=1, scale=1)

})

