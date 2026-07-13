orthodont_data <- function() {
  data <- as.data.frame(nlme::Orthodont)
  data$Occasion <- ave(
    seq_len(nrow(data)),
    data$Subject,
    FUN = seq_along
  )
  data
}

fit_orthodont_intercept <- function(ddf = "satterthwaite") {
  lmm(
    data = orthodont_data(),
    formula = distance ~ age + Sex,
    random = ~ 1 | Subject,
    ddf = ddf
  )
}

fit_multicentre_ar1 <- function() {
  lmm(
    data = multicentre,
    formula = Y ~ Drug * Time,
    random = ~ 1 | Center,
    repeated = ~ Time | Center:Subject,
    structure = "ar1"
  )
}
