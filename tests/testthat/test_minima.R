context("Testing fesminima and related functions")
test_that("Testing that fesminima correctly identifies energy minima", {
  # fesminima 2D
  myfes<-fes(acealanme, imax=2000)
  mins<-fesminima(myfes)
  minA <- summary(mins)[1,8]
  expect_equal(object=minA, expected=38.9, tolerance=0.1, scale=1)
  
  # fesminima 1D
  myfes<-fes(acealanme1d, imax=2000)
  mins<-fesminima(myfes)
  minA <- summary(mins)[1,6]
  expect_equal(object=minA, expected=42.2, tolerance=0.1, scale=1)
  
  # fesminima 2D with per=c(T,F)
  acealanme2<-acealanme
  acealanme2$per<-c(T,F)
  myfes<-fes(acealanme2, imax=2000)
  mins<-fesminima(myfes)
  minA <- summary(mins)[1,8]
  expect_equal(object=minA, expected=37.7, tolerance=0.1, scale=1)

  # fesminima 2D with per=c(F,T)
  acealanme2<-acealanme
  acealanme2$per<-c(F,T)
  myfes<-fes(acealanme2, imax=2000)
  mins<-fesminima(myfes)
  minA <- summary(mins)[1,8]
  expect_equal(object=minA, expected=44.1, tolerance=0.1, scale=1)

  # fesminima 2D with per=c(F,F)
  acealanme2<-acealanme
  acealanme2$per<-c(F,F)
  myfes<-fes(acealanme2, imax=2000)
  mins<-fesminima(myfes)
  minA <- summary(mins)[1,8]
  expect_equal(object=minA, expected=49.2, tolerance=0.1, scale=1)
  
  # fesminima 1D with per=c(F)
  acealanme1d2<-acealanme1d
  acealanme2$per<-c(F)
  myfes<-fes(acealanme1d2, imax=2000)
  mins<-fesminima(myfes)
  minA <- summary(mins)[1,6]
  expect_equal(object=minA, expected=42.2, tolerance=0.1, scale=1)

  # feprof and summary 2D
  myfes<-fes(acealanme, imax=2000)
  mins<-fesminima(myfes)
  profs<-feprof(mins)
  minA<-summary(profs)[2,7]
  expect_equal(object=minA, expected=-0.030, tolerance=0.01, scale=1)

  # feprof and summary 1D
  myfes<-fes(acealanme1d, imax=2000)
  mins<-fesminima(myfes)
  profs<-feprof(mins)
  minA<-summary(profs)[2,6]
  expect_equal(object=minA, expected=5.97, tolerance=0.01, scale=1)
  
  # feprof and summary 2D with pertc(T,F)
  acealanme2<-acealanme
  acealanme2$per<-c(T,F)
  myfes<-fes(acealanme2, imax=2000)
  mins<-fesminima(myfes)
  profs<-feprof(mins)
  minA<-summary(profs)[2,7]
  expect_equal(object=minA, expected=-0.787, tolerance=0.01, scale=1)

  # feprof and summary 2D with pertc(F,T)
  acealanme2<-acealanme
  acealanme2$per<-c(F,T)
  myfes<-fes(acealanme2, imax=2000)
  mins<-fesminima(myfes)
  profs<-feprof(mins)
  minA<-summary(profs)[2,7]
  expect_equal(object=minA, expected=0.068, tolerance=0.01, scale=1)

  # feprof and summary 2D with pertc(F,F)
  acealanme2<-acealanme
  acealanme2$per<-c(F,F)
  myfes<-fes(acealanme2, imax=2000)
  mins<-fesminima(myfes)
  profs<-feprof(mins)
  minA<-summary(profs)[2,7]
  expect_equal(object=minA, expected=-0.72, tolerance=0.01, scale=1)
  
  # feprof and summary 1D with pertc(F)
  acealanme1d2<-acealanme1d
  acealanme1d2$per<-c(F)
  myfes<-fes(acealanme1d2, imax=2000)
  mins<-fesminima(myfes)
  profs<-feprof(mins)
  minA<-summary(profs)[2,6]
  expect_equal(object=minA, expected=6.01, tolerance=0.01, scale=1)

  # fespoint 2D
  point<-fespoint(acealanme, c(3,3), imax=2000)
  expect_equal(object=minA, expected=-17.41, tolerance=0.01, scale=1)

  # fespoint 1D
  point<-fespoint(acealanme1d, c(3), imax=2000)
  expect_equal(object=minA, expected=-51.43, tolerance=0.01, scale=1)

  # fespoint 2D with pertc(T,F)
  acealanme2<-acealanme
  acealanme2$per<-c(T,F)
  point<-fespoint(acealanme2, c(3,3), imax=2000)
  expect_equal(object=minA, expected=-13.26, tolerance=0.01, scale=1)

  # fespoint 2D with pertc(F,T)
  acealanme2<-acealanme
  acealanme2$per<-c(F,T)
  point<-fespoint(acealanme2, c(3,3), imax=2000)
  expect_equal(object=minA, expected=-7.35, tolerance=0.01, scale=1)

  # fespoint 2D with pertc(F,F)
  acealanme2<-acealanme
  acealanme2$per<-c(F,F)
  point<-fespoint(acealanme2, c(3,3), imax=2000)
  expect_equal(object=minA, expected=-6.61, tolerance=0.01, scale=1)

  # fespoint 1D with pertc(F)
  acealanme1d2<-acealanme1d
  acealanme1d2$per<-c(F)
  point<-fespoint(acealanme1d2, c(3), imax=2000)
  expect_equal(object=minA, expected=-29.46, tolerance=0.01, scale=1)


  
})

