sas_data_message <- paste(
  "The original multicentre data from the thesis were not supplied;",
  "the SAS reference value cannot be verified against simulated data."
)

test_that("SAS Drug 1 LS-mean is reproduced", {
  skip(sas_data_message)
  expect_equal(NA_real_, 13.1569, tolerance = 1e-3)
})

test_that("SAS Drug 2 LS-mean is reproduced", {
  skip(sas_data_message)
  expect_equal(NA_real_, 18.14846, tolerance = 1e-3)
})

test_that("SAS Drug 3 LS-mean is reproduced", {
  skip(sas_data_message)
  expect_equal(NA_real_, 16.25425, tolerance = 1e-3)
})

test_that("SAS time LS-means are reproduced", {
  skip(sas_data_message)
  expect_equal(NA_real_, c(14.95442, 15.67009, 17.73251), tolerance = 1e-3)
})

test_that("SAS treatment difference is reproduced", {
  skip(sas_data_message)
  expect_equal(NA_real_, -4.991562, tolerance = 1e-3)
})

test_that("SAS treatment-difference standard error is reproduced", {
  skip(sas_data_message)
  expect_equal(NA_real_, 1.099005, tolerance = 1e-3)
})

test_that("SAS treatment-difference degrees of freedom are reproduced", {
  skip(sas_data_message)
  expect_equal(NA_real_, 46.10, tolerance = 1e-3)
})
