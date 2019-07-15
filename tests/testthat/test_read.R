context("Testing read.hills and related functions")
test_that("Testing that read.hills correctly loads hills and ralated functions analyze them", {
  # read.hills 2D
  l1<-"1 -1.409  2.808 0.3 0.3 1.111 10"
  l2<-"2 -2.505  2.791 0.3 0.3 1.111 10"
  l3<-"3 -2.346  2.754 0.3 0.3 1.069 10"
  l4<-"4 -1.198  2.872 0.3 0.3 1.074 10"
  fourhills<-c(l1,l2,l3,l4)
  tf <- tempfile()
  writeLines(fourhills, tf)
  h<-read.hills(tf, per=c(TRUE,TRUE))
  myfes<-fes(h)
  summyfes<-sum(myfes$fes)
  expect_equal(object=summyfes, expected=-4098, tolerance=1, scale=1)
  
  # +
  h<-read.hills(tf, per=c(TRUE,TRUE))
  h2<-h+h
  myfes<-fes(h2)
  summyfes<-sum(myfes$fes)
  expect_equal(object=summyfes, expected=-8195, tolerance=1, scale=1)
  
  # read.hills 1D
  l1<-"1 -2.600 0.3 1.111 10"
  l2<-"2 -1.970 0.3 1.106 10"
  l3<-"3 -1.287 0.3 1.107 10"
  l4<-"4 -0.877 0.3 1.092 10"
  fourhills<-c(l1,l2,l3,l4)
  tf <- tempfile()
  writeLines(fourhills, tf)
  h<-read.hills(tf, per=c(TRUE))
  myfes<-fes(h)
  summyfes<-sum(myfes$fes)
  expect_equal(object=summyfes, expected=-135.3, tolerance=0.1, scale=1)

  # +
  h<-read.hills(tf, per=c(TRUE))
  h2<-h+h
  myfes<-fes(h2)
  summyfes<-sum(myfes$fes)
  expect_equal(object=summyfes, expected=-270, tolerance=1, scale=1)

})

