context("Testing read.hills and related functions")
test_that("Testing that read.hills correctly loads hills and ralated functions analyze them", {
  # read.hills
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
  expect_equal(object=summyfes, expected=-4097, tolerance=1, scale=1)

})

