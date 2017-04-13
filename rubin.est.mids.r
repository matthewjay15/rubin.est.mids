rubin.est.mids <- function(analyses = NULL) {
  # (C) Matthew Jay, March 2017.
  m <- length(analyses)
  predictors <- analyses[[1]]$coefnames
  n.outcomes <- nrow(summary(analyses[[1]])$coefficients)
  logits <- lapply(analyses,
                   function(x) summary(x)$coefficients[, colnames(summary(x)$coefficient) %in% predictors])
  r.total <- matrix(0, nrow = n.outcomes, ncol = length(predictors))
  for (i in 1:m) {
    r.total <- as.matrix(logits[[i]]) + r.total
  }
  qbar <- r.total / m
  se <- lapply(analyses,
               function(x) summary(x)$standard.errors[, colnames(summary(x)$standard.errors) %in% predictors])
  r.total <- matrix(0, nrow = n.outcomes, ncol = length(predictors))
  for (i in 1:m) {
    r.total <- as.matrix(se[[i]]) + r.total
  }
  ubar <- r.total / m
  r.total <- matrix(0, nrow = n.outcomes, ncol = length(predictors))
  for (i in 1:m) {
    r.total <- ((as.matrix(logits[[i]]) - qbar)^2) + r.total
  }
  b <- (1 / (m-1)) * r.total
  t.var <- ubar + (1 + (1/m))*b
  or <- round(exp(qbar), 2)
  lcl <- round(exp(qbar - 1.96*t.var), 2)
  ucl <- round(exp(qbar + 1.96*t.var), 2)
  z <- qbar / t.var
  p <- round((1 - pnorm(abs(z), 0, 1)) * 2, 3)
  output <- matrix(NA, nrow = length(predictors), ncol = n.outcomes)
  if (length(predictors) > 1) {
    colnames(output) <- paste0(rownames(qbar),
                               ": logit (SE); OR (CI); p |")
    for (i in 1:n.outcomes) { # cols
      for (j in 1:length(predictors)) { # rows
        output[j, i] <- paste0(round(qbar[i, j], 3),
                               " (", round(t.var[i, j], 3), "); ",
                               or[i, j],
                               " (", lcl[i, j], ", ", ucl[i, j], "); ",
                               p[i, j], " |"
      } 
    }
  } else {
    colnames(output) <- paste0(names(qbar), ": logit (SE); OR (CI); p |")
    for (i in 1:n.outcomes) { # cols
      output[i] <- paste0(round(qbar[i], 3),
                          " (", round(t.var[i], 3), "); ",
                          or[i],
                          " (", lcl[i], ", ", ucl[i], "); ",
                          p[i], " |"
    }  
  }
  rownames(output) <- predictors
  print(noquote(output), right = T)
}