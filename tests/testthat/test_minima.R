context("Testing fesminima and related functions")
test_that("Testing that fesminima correctly identifies energy minima", {
  # fesminima
  myfes<-fes(acealanme, imax=2000)
  mins<-fesminima(myfes)
  minA <- mins$minima[1,6]
  expect_equal(object=minA, expected=-41.50, tolerance=0.02, scale=1)
  
})

