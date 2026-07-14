# Reproduce data sets used in official SAS PROC MIXED examples.
#
# Primary source: SAS Institute Inc. (2017), SAS/STAT 14.3 User's Guide,
# Chapter 79, The MIXED Procedure, Examples 79.1, 79.2, 79.5, and 79.6.
# https://go.documentation.sas.com/doc/en/statug/14.3/

sas_split_plot <- data.frame(
  Block = rep(seq_len(4L), each = 6L),
  A = rep(rep(seq_len(3L), each = 2L), times = 4L),
  B = rep(seq_len(2L), times = 12L),
  Y = c(
    56, 41, 50, 36, 39, 35,
    30, 25, 36, 28, 33, 30,
    32, 24, 31, 27, 15, 19,
    30, 25, 35, 30, 17, 18
  )
)
sas_split_plot$Block <- factor(sas_split_plot$Block)
sas_split_plot$A <- factor(sas_split_plot$A)
sas_split_plot$B <- factor(sas_split_plot$B)

growth_wide <- data.frame(
  Person = seq_len(27L),
  Gender = c(rep("F", 11L), rep("M", 16L)),
  y8 = c(
    21, 21, 20.5, 23.5, 21.5, 20, 21.5, 23, 20, 16.5, 24.5,
    26, 21.5, 23, 25.5, 20, 24.5, 22, 24, 23, 27.5, 23, 21.5,
    17, 22.5, 23, 22
  ),
  y10 = c(
    20, 21.5, 24, 24.5, 23, 21, 22.5, 23, 21, 19, 25,
    25, 22.5, 22.5, 27.5, 23.5, 25.5, 22, 21.5, 20.5, 28, 23,
    23.5, 24.5, 25.5, 24.5, 21.5
  ),
  y12 = c(
    21.5, 24, 24.5, 25, 22.5, 21, 23, 23.5, 22, 19, 28,
    29, 23, 24, 26.5, 22.5, 27, 24.5, 24.5, 31, 31, 23.5, 24,
    26, 25.5, 26, 23.5
  ),
  y14 = c(
    23, 25.5, 26, 26.5, 23.5, 22.5, 25, 24, 21.5, 19.5, 28,
    31, 26.5, 27.5, 27, 26, 28.5, 26.5, 25.5, 26, 31.5, 25, 28,
    29.5, 26, 30, 25
  )
)
sas_growth <- data.frame(
  Person = rep(growth_wide$Person, each = 4L),
  Gender = rep(growth_wide$Gender, each = 4L),
  Age = rep(c(8, 10, 12, 14), times = nrow(growth_wide)),
  y = as.vector(t(as.matrix(growth_wide[c("y8", "y10", "y12", "y14")])))
)
sas_growth$Person <- factor(sas_growth$Person)
sas_growth$Gender <- factor(sas_growth$Gender, levels = c("F", "M"))

random_coefficients_wide <- rbind(
  c(1, 0, 101.2, 103.3, 103.3, 102.1, 104.4, 102.4),
  c(1, 1, 98.8, 99.4, 99.7, 99.5, NA, NA),
  c(1, 3, 98.4, 99.0, 97.3, 99.8, NA, NA),
  c(1, 6, 101.5, 100.2, 101.7, 102.7, NA, NA),
  c(1, 9, 96.3, 97.2, 97.2, 96.3, NA, NA),
  c(1, 12, 97.3, 97.9, 96.8, 97.7, 97.7, 96.7),
  c(2, 0, 102.6, 102.7, 102.4, 102.1, 102.9, 102.6),
  c(2, 1, 99.1, 99.0, 99.9, 100.6, NA, NA),
  c(2, 3, 105.7, 103.3, 103.4, 104.0, NA, NA),
  c(2, 6, 101.3, 101.5, 100.9, 101.4, NA, NA),
  c(2, 9, 94.1, 96.5, 97.2, 95.6, NA, NA),
  c(2, 12, 93.1, 92.8, 95.4, 92.2, 92.2, 93.0),
  c(3, 0, 105.1, 103.9, 106.1, 104.1, 103.7, 104.6),
  c(3, 1, 102.2, 102.0, 100.8, 99.8, NA, NA),
  c(3, 3, 101.2, 101.8, 100.8, 102.6, NA, NA),
  c(3, 6, 101.1, 102.0, 100.1, 100.2, NA, NA),
  c(3, 9, 100.9, 99.5, 102.2, 100.8, NA, NA),
  c(3, 12, 97.8, 98.3, 96.9, 98.4, 96.9, 96.5)
)
colnames(random_coefficients_wide) <- c(
  "Batch", "Month", paste0("y", seq_len(6L))
)
sas_random_coefficients <- data.frame(
  Batch = rep(random_coefficients_wide[, "Batch"], each = 6L),
  Month = rep(random_coefficients_wide[, "Month"], each = 6L),
  Replicate = rep(seq_len(6L), times = nrow(random_coefficients_wide)),
  Y = as.vector(t(random_coefficients_wide[, paste0("y", seq_len(6L))]))
)
sas_random_coefficients$Batch <- factor(sas_random_coefficients$Batch)

line_source_wide <- data.frame(
  Block = c(1, 1, 1, 2, 2, 2, 3, 3, 3),
  Cult = c(
    "Luke", "Nugaines", "Bridger",
    "Nugaines", "Bridger", "Luke",
    "Nugaines", "Luke", "Bridger"
  ),
  matrix(c(
    2.4, 2.7, 5.6, 7.5, 7.9, 7.1, 6.1, 7.3, 7.4, 6.7, 3.8, 1.8,
    2.2, 2.2, 4.3, 6.3, 7.9, 7.1, 6.2, 5.3, 5.3, 5.2, 5.4, 2.9,
    2.9, 3.2, 5.1, 6.9, 6.1, 7.5, 5.6, 6.5, 6.6, 5.3, 4.1, 3.1,
    2.4, 2.2, 4.0, 5.8, 6.1, 6.2, 7.0, 6.4, 6.7, 6.4, 3.7, 2.2,
    2.6, 3.1, 5.7, 6.4, 7.7, 6.8, 6.3, 6.2, 6.6, 6.5, 4.2, 2.7,
    2.2, 2.7, 4.3, 6.9, 6.8, 8.0, 6.5, 7.3, 5.9, 6.6, 3.0, 2.0,
    1.8, 1.9, 3.7, 4.9, 5.4, 5.1, 5.7, 5.0, 5.6, 5.1, 4.2, 2.2,
    2.1, 2.3, 3.7, 5.8, 6.3, 6.3, 6.5, 5.7, 5.8, 4.5, 2.7, 2.3,
    2.7, 2.8, 4.0, 5.0, 5.2, 5.2, 5.9, 6.1, 6.0, 4.3, 3.1, 3.1
  ), nrow = 9L, byrow = TRUE)
)
names(line_source_wide)[3:14] <- paste0("Y", seq_len(12L))
sas_line_source <- data.frame(
  Block = rep(line_source_wide$Block, each = 12L),
  Cult = rep(line_source_wide$Cult, each = 12L),
  Sbplt = rep(seq_len(12L), times = nrow(line_source_wide)),
  Irrig = rep(c(seq_len(6L), 6:1), times = nrow(line_source_wide)),
  Dir = rep(rep(c("North", "South"), each = 6L), nrow(line_source_wide)),
  Y = as.vector(t(as.matrix(line_source_wide[paste0("Y", seq_len(12L))])))
)
sas_line_source$Block <- factor(sas_line_source$Block)
sas_line_source$Cult <- factor(sas_line_source$Cult)
sas_line_source$Dir <- factor(sas_line_source$Dir)
sas_line_source$Irrig <- factor(sas_line_source$Irrig)

stopifnot(
  nrow(sas_split_plot) == 24L,
  nrow(sas_growth) == 108L,
  nrow(sas_random_coefficients) == 108L,
  sum(is.na(sas_random_coefficients$Y)) == 24L,
  nrow(sas_line_source) == 108L
)

save(
  sas_split_plot,
  sas_growth,
  sas_random_coefficients,
  sas_line_source,
  file = "data/sas_examples.rda",
  compress = "xz"
)
