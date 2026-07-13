# Reproduce Table 5.20 from the source thesis.
#
# Source: Mahsa Mohseni Bonab, "Programmation R et SAS pour modèles
# linéaires mixtes", Table 5.20, p. 144. The table was transcribed from the
# locally supplied final manuscript and checked against its rendered PDF page.

response_matrices <- list(
  "R.1" = rbind(
    c(17, NA, NA),
    c(12, 14, 15),
    c(12, 11, 14),
    c(13, 13, 17),
    c(12, 13, NA)
  ),
  "R.2" = rbind(
    c(18, 19, 21),
    c(19, NA, NA),
    c(18, 19, NA),
    c(16, 15, 19),
    c(19, 18, NA)
  ),
  "R.3" = rbind(
    c(16, 16, 18),
    c(16, 16, 18),
    c(16, 16, 19),
    c(19, 19, 23),
    c(17, 16, 20)
  ),
  "S.1" = rbind(
    c(18, 18, 21),
    c(15, 14, 16),
    c(6, 6, NA),
    c(16, 17, 18)
  ),
  "S.2" = rbind(
    c(23, 23, 26),
    c(15, 15, 19),
    c(15, 16, 19),
    c(17, 17, 21)
  ),
  "S.3" = rbind(
    c(23, NA, NA),
    c(17, 17, 20),
    c(18, 19, 21),
    c(18, NA, NA)
  ),
  "T.1" = rbind(
    c(8, 10, 11),
    c(5, 7, 7),
    c(2, 5, 6),
    c(9, 11, 11),
    c(8, 9, 9),
    c(13, 14, 16),
    c(11, 12, NA),
    c(19, 20, NA)
  ),
  "T.2" = rbind(
    c(16, 18, NA),
    c(12, NA, NA),
    c(12, NA, NA),
    c(16, 17, 18),
    c(19, 20, 20),
    c(16, NA, NA),
    c(13, 16, 17),
    c(13, 15, NA)
  ),
  "T.3" = rbind(
    c(11, NA, NA),
    c(13, 14, 14),
    c(14, 16, 18),
    c(14, 15, 16),
    c(15, 16, 18),
    c(18, 18, 20),
    c(7, NA, NA),
    c(11, NA, NA)
  )
)

wide_groups <- lapply(names(response_matrices), function(group_name) {
  identifiers <- strsplit(group_name, ".", fixed = TRUE)[[1L]]
  responses <- response_matrices[[group_name]]
  data.frame(
    Center = identifiers[[1L]],
    Drug = identifiers[[2L]],
    Subject = seq_len(nrow(responses)),
    Y1 = responses[, 1L],
    Y2 = responses[, 2L],
    Y3 = responses[, 3L]
  )
})
multicentre_wide <- do.call(rbind, wide_groups)

multicentre <- data.frame(
  Center = rep(multicentre_wide$Center, each = 3L),
  Drug = rep(multicentre_wide$Drug, each = 3L),
  Subject = rep(multicentre_wide$Subject, each = 3L),
  Time = rep(seq_len(3L), nrow(multicentre_wide)),
  Y = as.vector(t(as.matrix(multicentre_wide[c("Y1", "Y2", "Y3")]))),
  stringsAsFactors = FALSE
)
multicentre$Center <- factor(multicentre$Center, levels = c("R", "S", "T"))
multicentre$Drug <- factor(multicentre$Drug, levels = as.character(1:3))
multicentre$Subject <- factor(multicentre$Subject, levels = as.character(1:8))
multicentre$Time <- factor(multicentre$Time, levels = as.character(1:3))

stopifnot(
  nrow(multicentre) == 153L,
  sum(is.na(multicentre$Y)) == 28L,
  sum(!is.na(multicentre$Y)) == 125L
)

save(multicentre, file = "data/multicentre.rda", compress = "xz")
