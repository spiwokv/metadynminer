context("Testing neb and related functions")
test_that("Testing that fesminima correctly identifies energy minima", {
  # neb
  myfes<-fes(acealanme)
  mins<-fesminima(myfes)
  myneb<-neb(mins, min1="A", min2="B")
  thalf<-summary(myneb)[1,2]
  expect_equal(object=thalf, expected=2.91, tolerance=0.01, scale=1)
  
})

