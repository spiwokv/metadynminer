context("Testing fes function")
test_that("Testing that fes calculates correctly sum of points", {
  myfes<-fes(acealanme, imax=2000)
  expect_identical(sum(myfes$fes), -1127743)
})

