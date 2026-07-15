sas_targets_path <- system.file(
  "validation",
  "sas",
  "sas-stat-14.3-targets.csv",
  package = "lmmix"
)
official_sas_targets <- utils::read.csv(
  sas_targets_path,
  na.strings = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

sas_targets <- function(model, quantity) {
  official_sas_targets[
    official_sas_targets$model == model &
      official_sas_targets$quantity == quantity,
    ,
    drop = FALSE
  ]
}

unstructured_target <- function(targets, dimension) {
  out <- matrix(NA_real_, nrow = dimension, ncol = dimension)
  indices <- regmatches(
    targets$term,
    regexec("UN\\(([0-9]+),([0-9]+)\\)", targets$term)
  )
  for (index in seq_len(nrow(targets))) {
    row <- as.integer(indices[[index]][[2L]])
    column <- as.integer(indices[[index]][[3L]])
    out[row, column] <- targets$estimate[[index]]
    out[column, row] <- targets$estimate[[index]]
  }
  out
}

growth_sas_coefficients <- function(fit) {
  transformation <- rbind(
    Intercept = c(1, 1, 0, 0),
    "Gender F" = c(0, -1, 0, 0),
    Age = c(0, 0, 1, 1),
    "Age*Gender F" = c(0, 0, 0, -1)
  )
  variance <- transformation %*% vcov(fit) %*% t(transformation)
  data.frame(
    term = rownames(transformation),
    estimate = as.vector(transformation %*% fixef(fit)),
    std_error = sqrt(diag(variance))
  )
}

official_split_fit <- lmm(
  sas_split_plot,
  Y ~ A * B,
  random = list(
    Block = ~ 1 | Block,
    A_Block = ~ 1 | Block:A
  ),
  method = "REML",
  ddf = "satterthwaite"
)

official_growth_un_fit <- lmm(
  sas_growth,
  y ~ Gender * Age,
  repeated = ~ Age | Person,
  structure = "un",
  method = "ML",
  ddf = "satterthwaite"
)

official_growth_cs_fit <- lmm(
  sas_growth,
  y ~ Gender * Age,
  repeated = ~ Age | Person,
  structure = "cs",
  method = "ML",
  ddf = "satterthwaite"
)

official_random_fit <- lmm(
  sas_random_coefficients,
  Y ~ Month,
  random = ~ 1 + Month | Batch,
  method = "REML",
  ddf = "satterthwaite"
)

official_line_fit <- lmm(
  sas_line_source,
  Y ~ (Cult + Dir + Irrig)^2,
  random = list(
    Block = ~ 1 | Block,
    Block.Dir = ~ 1 | Block:Dir,
    Block.Irrig = ~ 1 | Block:Irrig
  ),
  repeated = ~ Sbplt | Block:Cult,
  structure = "toep(4)",
  method = "REML",
  ddf = "residual"
)

test_that("official SAS example data are reproduced exactly", {
  expect_equal(nrow(sas_split_plot), 24L)
  expect_equal(nrow(sas_growth), 108L)
  expect_equal(nrow(sas_random_coefficients), 108L)
  expect_equal(sum(is.na(sas_random_coefficients$Y)), 24L)
  expect_equal(nrow(sas_line_source), 108L)
})

test_that("SAS split-plot covariance and likelihood targets are reproduced", {
  covariance_target <- sas_targets("split_plot", "covariance")
  covariance_fit <- VarCorr(official_split_fit)
  expect_absolute_error_below(
    covariance_fit$estimate,
    covariance_target$estimate,
    tolerance = 1e-3
  )

  deviance_target <- sas_targets("split_plot", "deviance")$estimate
  expect_absolute_error_below(
    deviance(official_split_fit),
    deviance_target,
    tolerance = 1e-5
  )
})

test_that("SAS split-plot Type III tests are reproduced", {
  target <- sas_targets("split_plot", "type3")
  result <- anova(official_split_fit)

  expect_equal(result$num.df, target$num_df)
  expect_absolute_error_below(result$den.df, target$den_df, 1e-2)
  expect_absolute_error_below(result$statistic, target$statistic, 1e-2)
  expect_absolute_error_below(
    result$p.value,
    as.numeric(target$p_value),
    1e-4
  )
})

test_that("SAS unstructured repeated covariance is reproduced", {
  target <- sas_targets("growth_un", "covariance")
  target_matrix <- unstructured_target(target, dimension = 4L)
  expect_absolute_error_below(
    official_growth_un_fit$covariance$residual_base,
    target_matrix,
    tolerance = 1e-3
  )

  deviance_target <- sas_targets("growth_un", "deviance")$estimate
  expect_absolute_error_below(
    deviance(official_growth_un_fit),
    deviance_target,
    tolerance = 1e-5
  )
})

test_that("SAS unstructured fixed-effect solution is reproduced", {
  target <- sas_targets("growth_un", "fixed")
  result <- growth_sas_coefficients(official_growth_un_fit)

  expect_identical(result$term, target$term)
  expect_absolute_error_below(result$estimate, target$estimate, 1e-3)
  expect_absolute_error_below(result$std_error, target$std_error, 1e-3)
})

test_that("SAS unstructured Type III statistics are reproduced", {
  target <- sas_targets("growth_un", "type3")
  result <- anova(official_growth_un_fit)

  expect_absolute_error_below(result$statistic, target$statistic, 1e-2)
})

test_that("SAS compound-symmetry covariance is reproduced", {
  target <- sas_targets("growth_cs", "covariance")
  residual_matrix <- official_growth_cs_fit$covariance$residual_base
  observed <- c(
    residual_matrix[1L, 2L],
    residual_matrix[1L, 1L] - residual_matrix[1L, 2L]
  )
  expect_absolute_error_below(observed, target$estimate, tolerance = 1e-3)

  deviance_target <- sas_targets("growth_cs", "deviance")$estimate
  expect_absolute_error_below(
    deviance(official_growth_cs_fit),
    deviance_target,
    tolerance = 1e-5
  )
})

test_that("SAS compound-symmetry fixed effects are reproduced", {
  target <- sas_targets("growth_cs", "fixed")
  result <- growth_sas_coefficients(official_growth_cs_fit)

  expect_identical(result$term, target$term)
  expect_absolute_error_below(result$estimate, target$estimate, 1e-3)
  expect_absolute_error_below(result$std_error, target$std_error, 1e-3)

  type3_target <- sas_targets("growth_cs", "type3")
  type3_result <- anova(official_growth_cs_fit)
  expect_absolute_error_below(
    type3_result$statistic,
    type3_target$statistic,
    1e-2
  )
})

test_that("SAS random-coefficients covariance and likelihood are reproduced", {
  target <- sas_targets("random_coefficients", "covariance")
  observed <- c(
    official_random_fit$covariance$g$Batch[1L, 1L],
    official_random_fit$covariance$g$Batch[2L, 1L],
    official_random_fit$covariance$g$Batch[2L, 2L],
    official_random_fit$covariance$residual_base[1L, 1L]
  )
  expect_absolute_error_below(observed, target$estimate, tolerance = 1e-3)

  deviance_target <- sas_targets("random_coefficients", "deviance")$estimate
  expect_absolute_error_below(
    deviance(official_random_fit),
    deviance_target,
    tolerance = 1e-5
  )
})

test_that("SAS random-coefficients fixed and random solutions are reproduced", {
  fixed_target <- sas_targets("random_coefficients", "fixed")
  expect_absolute_error_below(
    fixef(official_random_fit),
    fixed_target$estimate,
    tolerance = 5e-3
  )
  expect_absolute_error_below(
    sqrt(diag(vcov(official_random_fit))),
    fixed_target$std_error,
    tolerance = 1e-3
  )

  random_target <- sas_targets("random_coefficients", "random")
  random_result <- ranef(official_random_fit)[[1L]]
  observed <- as.vector(t(as.matrix(random_result[c("(Intercept)", "Month")])))
  expect_absolute_error_below(
    observed,
    random_target$estimate,
    tolerance = 1e-3
  )

  type3_target <- sas_targets("random_coefficients", "type3")
  type3_result <- anova(official_random_fit)
  expect_absolute_error_below(type3_result$den.df, type3_target$den_df, 1e-2)
  expect_absolute_error_below(
    type3_result$statistic,
    type3_target$statistic,
    1e-2
  )
  expect_absolute_error_below(type3_result$p.value, 0.0478, 1e-4)
})

test_that("SAS line-source Toeplitz covariance is reproduced", {
  target <- sas_targets("line_source", "covariance")
  residual <- official_line_fit$covariance$residual_base
  observed <- c(
    official_line_fit$covariance$g$Block[1L, 1L],
    official_line_fit$covariance$g$Block.Dir[1L, 1L],
    official_line_fit$covariance$g$Block.Irrig[1L, 1L],
    residual[1L, 2L],
    residual[1L, 3L],
    residual[1L, 4L],
    residual[1L, 1L]
  )
  expect_absolute_error_below(observed, target$estimate, tolerance = 1e-4)

  deviance_target <- sas_targets("line_source", "deviance")$estimate
  expect_absolute_error_below(
    deviance(official_line_fit),
    deviance_target,
    tolerance = 1e-5
  )
})
