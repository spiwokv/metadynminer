context("Testing fesminima and related functions")
test_that("Testing that fesminima correctly identifies energy minima", {
  # fesminima 2D
  myfes<-fes(acealanme, imax=2000)
  mins<-fesminima(myfes)
  minA <- summary(mins)[1,8]
  expect_equal(object=minA, expected=38.3, tolerance=0.1, scale=1)
  
  # fesminima 1D
  myfes<-fes(acealanme1d, imax=2000)
  mins<-fesminima(myfes)
  minA <- summary(mins)[1,6]
  expect_equal(object=minA, expected=40.8, tolerance=0.1, scale=1)
  
  # fesminima 2D with per=c(T,F)
  acealanme2<-acealanme
  acealanme2$per<-c(T,F)
  myfes<-fes(acealanme, imax=2000)
  mins<-fesminima(myfes)
  minA <- summary(mins)[1,8]
  expect_equal(object=minA, expected=39.0, tolerance=0.1, scale=1)

  # fesminima 2D with per=c(F,T)
  acealanme2<-acealanme
  acealanme2$per<-c(F,T)
  myfes<-fes(acealanme, imax=2000)
  mins<-fesminima(myfes)
  minA <- summary(mins)[1,8]
  expect_equal(object=minA, expected=42.0, tolerance=0.1, scale=1)

  # fesminima 2D with per=c(F,F)
  acealanme2<-acealanme
  acealanme2$per<-c(F,F)
  myfes<-fes(acealanme, imax=2000)
  mins<-fesminima(myfes)
  minA <- summary(mins)[1,8]
  expect_equal(object=minA, expected=48.6, tolerance=0.1, scale=1)
  
  # fesminima 1D with per=c(F)
  acealanme1d2<-acealanme1d
  acealanme2$per<-c(F)
  myfes<-fes(acealanme1d2, imax=2000)
  mins<-fesminima(myfes)
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
  
  # feprof and summary 2D with pertc(T,F)
  acealanme2<-acealanme
  acealanme2$per<-c(T,F)
  myfes<-fes(acealanme2, imax=2000)
  mins<-fesminima(myfes)
  profs<-feprof(mins)
  minA<-summary(profs)[2,7]
  expect_equal(object=minA, expected=-0.60, tolerance=0.01, scale=1)

  # feprof and summary 2D with pertc(F,T)
  acealanme2<-acealanme
  acealanme2$per<-c(F,T)
  myfes<-fes(acealanme2, imax=2000)
  mins<-fesminima(myfes)
  profs<-feprof(mins)
  minA<-summary(profs)[2,7]
  expect_equal(object=minA, expected=0.13, tolerance=0.01, scale=1)

  # feprof and summary 2D with pertc(F,F)
  acealanme2<-acealanme
  acealanme2$per<-c(F,F)
  myfes<-fes(acealanme2, imax=2000)
  mins<-fesminima(myfes)
  profs<-feprof(mins)
  minA<-summary(profs)[2,7]
  expect_equal(object=minA, expected=-0.67, tolerance=0.01, scale=1)
  
  # feprof and summary 1D with pertc(F)
  acealanme1d2<-acealanme1d
  acealanme1d2$per<-c(F)
  myfes<-fes(acealanme1d2, imax=2000)
  mins<-fesminima(myfes)
  profs<-feprof(mins)
  minA<-summary(profs)[2,6]
  expect_equal(object=minA, expected=5.91, tolerance=0.01, scale=1)

})

