context("Testing fesminima and related functions")
test_that("Testing that fesminima correctly identifies energy minima", {
  # fesminima 2D
  myfes<-fes(acealanme, imax=2000)
  mins<-fesminima(myfes)
  minA <- mins$minima[1,6]
  expect_equal(object=minA, expected=-41.50, tolerance=0.02, scale=1)

  # summary
  minA <- summary(mins)[1,8]
  expect_equal(object=minA, expected=38.3, tolerance=0.1, scale=1)
  
  # fesminima 1D
  myfes<-fes(acealanme1d, imax=2000)
  mins<-fesminima(myfes)
  minA <- mins$minima[1,4]
  expect_equal(object=minA, expected=-72.12, tolerance=0.02, scale=1)

  # summary
  minA <- summary(mins)[1,6]
  expect_equal(object=minA, expected=40.8, tolerance=0.1, scale=1)
 
  # feprof and summary 2D
  myfes<-fes(acealanme, imax=2000)
  mins<-fesminima(myfes)
  profs<-feprof(mins)
  minA<-summary(profs)[2,7]
  expect_equal(object=minA, expected=0.10, tolerance=0.01, scale=1)

  # feprof and summary 1D
  myfes<-fes(acealanme1d, imax=2000)
  mins<-fesminima(myfes)
  profs<-feprof(mins)
  minA<-summary(profs)[2,6]
  expect_equal(object=minA, expected=5.89, tolerance=0.01, scale=1)

})

