sas_fit <- fit_multicentre_ar1()
sas_covariance <- VarCorr(sas_fit)
sas_fixed <- fixef(sas_fit)
sas_type3 <- anova(sas_fit)
sas_drug_means <- lsmeans(sas_fit, ~Drug)
sas_time_means <- lsmeans(sas_fit, ~Time)
sas_drug_pairs <- lsmeans(sas_fit, pairwise ~ Drug)$contrasts
sas_time_pairs <- lsmeans(sas_fit, pairwise ~ Time)$contrasts

test_that("the multilocation example data are complete", {
  expect_equal(nrow(multicentre), 153L)
  expect_equal(sum(is.na(multicentre$Y)), 28L)
  expect_equal(sum(!is.na(multicentre$Y)), 125L)
})

test_that("SAS covariance parameters are reproduced to 1e-3", {
  expect_absolute_error_below(
    sas_covariance$estimate,
    c(5.1737, 10.6702, 0.9351),
    tolerance = 1e-3
  )
  expect_absolute_error_below(
    deviance(sas_fit),
    490.62570532,
    tolerance = 1e-3
  )
})

test_that("fixed-effect estimates are reproduced to 1e-3", {
  expected <- c(
    12.0524590,
    4.7647059,
    3.9411765,
    0.9210892,
    2.3922237,
    -0.1427475,
    -0.4735052,
    0.8233152,
    0.3342713
  )
  expect_absolute_error_below(sas_fixed, expected, tolerance = 1e-3)
})

test_that("stored type III F statistics are reproduced", {
  expect_equal(round(sas_type3$statistic, 2), c(11.43, 59.27, 1.35))
  expect_absolute_error_below(
    sas_type3$p.value[[3L]],
    0.2597,
    tolerance = 1e-3
  )
})

test_that("SAS Drug 1 LS-mean is reproduced to 1e-3", {
  expect_absolute_error_below(
    sas_drug_means$estimate[[1L]],
    13.1568,
    tolerance = 1e-3
  )
  expect_absolute_error_below(
    sas_drug_means$std.error[[1L]],
    1.5290,
    tolerance = 1e-3
  )
  expect_equal(round(sas_drug_means$df[[1L]], 2), 2.98)
})

test_that("SAS Drug 2 LS-mean is reproduced to 1e-3", {
  expect_absolute_error_below(
    sas_drug_means$estimate[[2L]],
    18.1484,
    tolerance = 1e-3
  )
  expect_absolute_error_below(
    sas_drug_means$std.error[[2L]],
    1.5324,
    tolerance = 1e-3
  )
  expect_equal(round(sas_drug_means$df[[2L]], 2), 3.01)
})

test_that("SAS Drug 3 LS-mean is reproduced to 1e-3", {
  expect_absolute_error_below(
    sas_drug_means$estimate[[3L]],
    17.0516,
    tolerance = 1e-3
  )
  expect_absolute_error_below(
    sas_drug_means$std.error[[3L]],
    1.5327,
    tolerance = 1e-3
  )
  expect_equal(round(sas_drug_means$df[[3L]], 2), 3.01)
})

test_that("SAS Time 1 LS-mean is reproduced to 1e-3", {
  expect_absolute_error_below(
    sas_time_means$estimate[[1L]],
    14.9543,
    tolerance = 1e-3
  )
  expect_absolute_error_below(
    sas_time_means$std.error[[1L]],
    1.3961,
    tolerance = 1e-3
  )
  expect_equal(round(sas_time_means$df[[1L]], 2), 2.08)
})

test_that("SAS Time 2 LS-mean is reproduced to 1e-3", {
  expect_absolute_error_below(
    sas_time_means$estimate[[2L]],
    15.6700,
    tolerance = 1e-3
  )
  expect_absolute_error_below(
    sas_time_means$std.error[[2L]],
    1.3985,
    tolerance = 1e-3
  )
  expect_equal(round(sas_time_means$df[[2L]], 2), 2.09)
})

test_that("SAS Time 3 LS-mean is reproduced to 1e-3", {
  expect_absolute_error_below(
    sas_time_means$estimate[[3L]],
    17.7324,
    tolerance = 1e-3
  )
  expect_absolute_error_below(
    sas_time_means$std.error[[3L]],
    1.4032,
    tolerance = 1e-3
  )
  expect_equal(round(sas_time_means$df[[3L]], 2), 2.12)
})

test_that("SAS Drug 1 minus Drug 2 contrast is reproduced to 1e-3", {
  comparison <- sas_drug_pairs[1L, ]
  expect_absolute_error_below(
    comparison$estimate,
    -4.9916,
    tolerance = 1e-3
  )
  expect_absolute_error_below(
    comparison$std.error,
    1.0990,
    tolerance = 1e-3
  )
  expect_equal(round(comparison$df, 1), 46.1)
  expect_lt(
    abs(as.numeric(comparison$p.value) - 4.0e-5),
    1e-6
  )
})

test_that("SAS Drug 1 minus Drug 3 contrast is reproduced to 1e-3", {
  comparison <- sas_drug_pairs[2L, ]
  expect_absolute_error_below(
    comparison$estimate,
    -3.8948,
    tolerance = 1e-3
  )
  expect_absolute_error_below(
    comparison$std.error,
    1.0987,
    tolerance = 1e-3
  )
  expect_equal(round(comparison$df, 1), 46.1)
  expect_equal(round(comparison$p.value, 4), 0.0009)
})

test_that("SAS Drug 2 minus Drug 3 contrast is reproduced to 1e-3", {
  comparison <- sas_drug_pairs[3L, ]
  expect_absolute_error_below(
    comparison$estimate,
    1.0968,
    tolerance = 1e-3
  )
  expect_absolute_error_below(
    comparison$std.error,
    1.1043,
    tolerance = 1e-3
  )
  expect_equal(round(comparison$df, 1), 46.9)
  expect_equal(round(comparison$p.value, 4), 0.3257)
})

test_that("SAS Time 1 minus Time 2 contrast is reproduced to 1e-3", {
  comparison <- sas_time_pairs[1L, ]
  expect_absolute_error_below(
    comparison$estimate,
    -0.7157,
    tolerance = 1e-3
  )
  expect_absolute_error_below(
    comparison$std.error,
    0.1845,
    tolerance = 1e-3
  )
  expect_equal(round(comparison$df, 1), 68.1)
  expect_equal(round(comparison$p.value, 4), 0.0002)
})

test_that("SAS Time 1 minus Time 3 contrast is reproduced to 1e-3", {
  comparison <- sas_time_pairs[2L, ]
  expect_absolute_error_below(
    comparison$estimate,
    -2.7781,
    tolerance = 1e-3
  )
  expect_absolute_error_below(
    comparison$std.error,
    0.2714,
    tolerance = 1e-3
  )
  expect_equal(round(comparison$df, 1), 73.2)
  expect_lt(comparison$p.value, 0.0001)
})

test_that("SAS Time 2 minus Time 3 contrast is reproduced to 1e-3", {
  comparison <- sas_time_pairs[3L, ]
  expect_absolute_error_below(
    comparison$estimate,
    -2.0624,
    tolerance = 1e-3
  )
  expect_absolute_error_below(
    comparison$std.error,
    0.2056,
    tolerance = 1e-3
  )
  expect_equal(round(comparison$df, 1), 68.6)
  expect_lt(comparison$p.value, 0.0001)
})
