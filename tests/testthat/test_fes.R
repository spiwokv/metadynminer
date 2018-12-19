context("Testing fes function")
test_that("Testing that fes calculates correctly sum of points", {
  myfes<-fes(acealanme, imax=2000)
  summyfes<-sum(myfes$fes)
  expect_identical(summyfes, -1127743, tolerance=1, scale = 1)
})

