##########
## Functions to fit unleash
##########


#' Wrapper for \code{\link{calc_loglik}}
#'
#' @inheritParams ashr::calc_loglik
#' @param fitted_g The fitted g
#' @param xi The variance inflation parameter.
#' @param pen The variance inflation penalty.
#'
#' @author David Gerard
#'
var_inflate_obj <- function(xi, fitted_g, data, pen = 1) {
  data$s <- data$s * xi
  llike <- ashr::calc_loglik(g = fitted_g, data = data)
  vpen <- -pen / xi
  return(llike + vpen)
}

#' UNimodally Leveraging Empirical-null with Adaptive SHrinkage
#'
#' This is a wrapper for \code{\link{ash}} with a variance inflation term that
#' expands the variances to as large a value as possible while allowing
#' the prior to be unimodal.
#'
#' @inheritParams ashr::ash
#' @param pen The variance inflation penalty.
#' @param maxiter The maximum number of iterations to run.
#' @param tol The stopping criterion.
#' @param plot_update A logical. Should we plot the updates?
#'
#' @author David Gerard
#'
#' @export
#'
unleash <- function(betahat, sebetahat,
                    mixcompdist = c("uniform", "halfuniform", "normal", "+uniform", "-uniform"),
                    df = NULL, pen = 1, maxiter = 500, tol = 10^-6, plot_update = FALSE, ...) {

  mixcompdist <- match.arg(mixcompdist)
  aout <- ashr::ash.workhorse(betahat = betahat, sebetahat = sebetahat, mixcompdist = mixcompdist,
                        df = df)

  ## Set initial value ---
  xi <- 1
  llike_new <- var_inflate_obj(xi = xi, fitted_g = aout$fitted_g, data = aout$data, pen = pen)
  iterindex <- 1
  err <- tol + 1

  if (plot_update) {
    dat <- data.frame(Index = 0, llike = llike_new)
  }

  while (iterindex < maxiter & err > tol) {
    llike_old <- llike_new
    oout <- stats::optim(par = xi, fn = var_inflate_obj, fitted_g = aout$fitted_g, data = aout$data,
                         method = "Brent", upper = 10, lower = 0, control = list(fnscale = -1),
                         pen = pen)
    xi <- oout$par
    new_sebetahat <- xi * sebetahat
    aout <- ashr::ash.workhorse(betahat = betahat, sebetahat = new_sebetahat, mixcompdist = mixcompdist,
                          df = df)
    llike_new <- var_inflate_obj(xi = xi, fitted_g = aout$fitted_g, data = aout$data, pen = pen)
    err <- abs(llike_new - llike_old)

    if (plot_update) {
      dat <- rbind(dat, c(iterindex, llike_new))
      pl <- ggplot2::ggplot(data = dat, mapping = ggplot2::aes_string(x = "Index", y = "llike")) +
        ggplot2::geom_line() +
        ggplot2::theme_bw() +
        ggplot2::ylab("Penalized Log-Likelihood") +
        ggplot2::ggtitle(paste0("Error = ", format(err, digits = 2)))
      print(pl)
    }

    iterindex <- iterindex + 1
  }

  aout$sd_inflate <- xi
  return(aout)
}

