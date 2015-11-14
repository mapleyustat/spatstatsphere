library(spatstatsphere)
context("Distances on sphere.")

test_that("Distance on unit sphere is correct", {
  expect_equal(unitspheredist(-pi/2, 0, pi/2, 0), pi)
  expect_equal(unitspheredist(pi/4, 0, pi/4, 0), 0)
})
