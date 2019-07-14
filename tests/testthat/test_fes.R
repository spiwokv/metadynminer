context("Testing fes function")
test_that("Testing that fes calculates correctly sum of points", {
  # fes function 2D
  myfes<-fes(acealanme, imax=2000)
  summyfes<-sum(myfes$fes)
  expect_equal(object=summyfes, expected=-1118950, tolerance=1, scale=1)
  
  # fes function 2D
  myfes<-fes2(acealanme, imax=2000)
  summyfes<-sum(myfes$fes)
  expect_equal(object=summyfes, expected=-1118950, tolerance=1, scale=1)
  
  # fes function 1D
  myfes<-fes(acealanme1d, imax=2000)
  summyfes<-sum(myfes$fes)
  expect_equal(object=summyfes, expected=-14733, tolerance=1, scale=1)
  
  # fes function 1D
  myfes<-fes2(acealanme1d, imax=2000)
  summyfes<-sum(myfes$fes)
  expect_equal(object=summyfes, expected=-14733, tolerance=1, scale=1)

  #+-*/
  summyfes2<-summyfes+summyfes
  expect_equal(object=summyfes2, expected=-29466, tolerance=1, scale=1)
  summyfes2<-summyfes-summyfes
  expect_equal(object=summyfes2, expected=0, tolerance=1, scale=1)
  summyfes2<-2*summyfes
  expect_equal(object=summyfes2, expected=-29466, tolerance=1, scale=1)
  summyfes2<-summyfes/0.5
  expect_equal(object=summyfes2, expected=-29466, tolerance=1, scale=1)

  # fes2d21d
  myfes<-fes2d21d(acealanme, imax=2000)
  summyfes<-sum(myfes$fes)
  expect_equal(object=summyfes, expected=-9042, tolerance=1, scale=1)
  myfes<-fes2d21d(acealanme, imax=2000, remdim=1)
  summyfes<-sum(myfes$fes)
  expect_equal(object=summyfes, expected=-10791, tolerance=1, scale=1)
  
  # fes function 2D with per=c(T,F)
  acealanme2 <- acealanme
  acealanme2$per <- c(T,F)
  myfes<-fes(acealanme2, imax=2000)
  summyfes<-sum(myfes$fes)
  expect_equal(object=summyfes, expected=-1009738, tolerance=1, scale=1)
  
  # fes function 2D with per=c(F,T)
  acealanme2 <- acealanme
  acealanme2$per <- c(F,T)
  myfes<-fes(acealanme2, imax=2000)
  summyfes<-sum(myfes$fes)
  expect_equal(object=summyfes, expected=-1011742, tolerance=1, scale=1)

  # fes function 2D with per=c(F,F)
  acealanme2 <- acealanme
  acealanme2$per <- c(F,F)
  myfes<-fes(acealanme2, imax=2000)
  summyfes<-sum(myfes$fes)
  expect_equal(object=summyfes, expected=-913007, tolerance=1, scale=1)

  # fes2 function 2D with per=c(T,F)
  acealanme2 <- acealanme
  acealanme2$per <- c(T,F)
  myfes<-fes2(acealanme2, imax=2000)
  summyfes<-sum(myfes$fes)
  expect_equal(object=summyfes, expected=-1009724, tolerance=1, scale=1)
  
  # fes2 function 2D with per=c(F,T)
  acealanme2 <- acealanme
  acealanme2$per <- c(F,T)
  myfes<-fes2(acealanme2, imax=2000)
  summyfes<-sum(myfes$fes)
  expect_equal(object=summyfes, expected=-1011740, tolerance=1, scale=1)

  # fes2 function 2D with per=c(F,F)
  acealanme2 <- acealanme
  acealanme2$per <- c(F,F)
  myfes<-fes2(acealanme2, imax=2000)
  summyfes<-sum(myfes$fes)
  expect_equal(object=summyfes, expected=-912991, tolerance=1, scale=1)
  
  # fes function 1D with per=c(F)
  acealanme1d2 <- acealanme1d
  acealanme1d2$per <- c(F)
  myfes<-fes(acealanme1d2, imax=2000)
  summyfes<-sum(myfes$fes)
  expect_equal(object=summyfes, expected=-13303, tolerance=1, scale=1)

  # fes2 function 1D with per=c(F)
  acealanme1d2 <- acealanme1d
  acealanme1d2$per <- c(F)
  myfes<-fes2(acealanme1d2, imax=2000)
  summyfes<-sum(myfes$fes)
  expect_equal(object=summyfes, expected=-13303, tolerance=1, scale=1)
  
  # fes2d21d with per=(T,F)
  acealanme2 <- acealanme
  acealanme2$per <- c(T,F)
  myfes<-fes2d21d(acealanme2, imax=2000)
  summyfes<-sum(myfes$fes)
  expect_equal(object=summyfes, expected=-8890, tolerance=1, scale=1)
  
  # fes2d21d with per=(F,T)
  acealanme2 <- acealanme
  acealanme2$per <- c(F,T)
  myfes<-fes2d21d(acealanme2, imax=2000)
  summyfes<-sum(myfes$fes)
  expect_equal(object=summyfes, expected=-8490, tolerance=1, scale=1)

  # fes2d21d with per=(F,F)
  acealanme2 <- acealanme
  acealanme2$per <- c(F,F)
  myfes<-fes2d21d(acealanme2, imax=2000)
  summyfes<-sum(myfes$fes)
  expect_equal(object=summyfes, expected=-8363, tolerance=1, scale=1)

})

