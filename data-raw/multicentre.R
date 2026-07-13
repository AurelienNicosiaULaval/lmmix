# Generate the documented example data.
#
# These data mimic the design described in Annex B. They are not the original
# data used to obtain the SAS reference values in that annex.

set.seed(20260713)

multicentre <- expand.grid(
  Time = factor(1:3),
  SubjectNumber = 1:6,
  Drug = factor(1:3),
  Center = factor(c("R", "S", "T")),
  KEEP.OUT.ATTRS = FALSE
)
multicentre$Subject <- factor(sprintf(
  "%s-D%s-%02d",
  multicentre$Center,
  multicentre$Drug,
  multicentre$SubjectNumber
))

center_effect <- c(R = -1.4, S = 0.2, T = 1.2)
drug_effect <- c(`1` = 0, `2` = 4.7, `3` = 3.1)
time_effect <- c(`1` = 0, `2` = 0.8, `3` = 2.5)
interaction_effect <- matrix(
  c(0, 0.1, 0.3, 0, 0.4, 0.8, 0, -0.2, 0.2),
  nrow = 3,
  byrow = TRUE
)
subject_effect <- stats::rnorm(nlevels(multicentre$Subject), sd = 1.4)
names(subject_effect) <- levels(multicentre$Subject)

innovation <- stats::rnorm(nrow(multicentre), sd = 1.2)
residual <- numeric(nrow(multicentre))
for (subject in levels(multicentre$Subject)) {
  rows <- which(multicentre$Subject == subject)
  rows <- rows[order(multicentre$Time[rows])]
  residual[rows[[1L]]] <- innovation[rows[[1L]]] / sqrt(1 - 0.55^2)
  for (index in 2:length(rows)) {
    residual[rows[[index]]] <-
      0.55 * residual[rows[[index - 1L]]] + innovation[rows[[index]]]
  }
}

drug_index <- as.integer(multicentre$Drug)
time_index <- as.integer(multicentre$Time)
multicentre$Y <-
  13 +
  center_effect[as.character(multicentre$Center)] +
  drug_effect[as.character(multicentre$Drug)] +
  time_effect[as.character(multicentre$Time)] +
  interaction_effect[cbind(drug_index, time_index)] +
  subject_effect[as.character(multicentre$Subject)] +
  residual

missing_rows <- c(5, 17, 38, 74, 92, 119, 145, 160)
multicentre$Y[missing_rows] <- NA_real_
multicentre$SubjectNumber <- NULL
multicentre <- multicentre[c("Center", "Drug", "Subject", "Time", "Y")]

usethis::use_data(multicentre, overwrite = TRUE)
