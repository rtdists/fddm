#' Methods for ddm objects
#' 
#' Implemented S3 methods for objects of class \code{ddm} as returned by
#' function \code{\link{ddm}()}.
#' 
#' @param object,x object of class \code{ddm}
#' @param formula see \code{\link{model.frame}}
#' @param dpar which distributional parameter or DDM parameter to focus on. In
#'   addition to the five DDM parameters \code{c("drift", "boundary", "ndt",
#'   "bias", "sv")}, some methods accept \code{"full"} which returns information
#'   for all estimated parameters.
#' @param digits minimal number of significant digits, see
#'   \code{\link{print.default}}.
#' @param ... further arguments passed to or from other methods.
#' 
#' @details The methods should fail with an informative error if a
#'   distributional parameter is selected in \code{dpar} that is fixed and not
#'   estimated.
#' 
#' @name ddm-methods
NULL

## some of the methods are implemented based on the corresponding methods
## in package betareg: https://cran.r-project.org/package=betareg

## MISSING METHODS: residuals

#' @rdname ddm-methods
#' @export
print.ddm <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("\nCall:",
      deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)),
      "", sep = "\n")

  cat("DDM fit with",  length(x$coefficients), "estimated and",
      length(x$fixed_dpar), "fixed distributional parameters.\n")
  if (length(x$fixed_dpar) > 0) {
    cat("Fixed:",
        paste(names(x$fixed_dpar), x$fixed_dpar, sep = " = ", collapse = ", "),
        "\n")
  }
  cat("\n")

  cat(paste("drift coefficients (", x$link$drift$name,
            " link):\n", sep = ""))
  print.default(format(x$coefficients$drift, digits = digits), print.gap = 2,
                quote = FALSE)
  cat("\n")

  if (length(x$coefficients) > 1) {
    for (par in seq_len(length(x$coefficients) - 1)) {
      cur_par <- names(x$coefficients)[par + 1]
      cat(paste(cur_par, " coefficients (", x$link[[cur_par]]$name,
                " link):\n", sep = ""))
      print.default(format(x$coefficients[[cur_par]], digits = digits),
                    print.gap = 2, quote = FALSE)
      cat("\n")
    }
  }

  if (!x$converged) {
    warning("optimisation algorithm indicates that model did not converge.\n",
            "Proceed with caution!\n", call. = FALSE)
  }

  invisible(x)
}

#' @rdname ddm-methods
#' @export
summary.ddm <- function(object, ...) {
  ## extend coefficient table for each dpar
  dpars <- object$dpar
  tbl <- list()
  for (i in 1:length(dpars)) {
    cf <- object$coefficients[[dpars[i]]]
    se <- sqrt(diag(object$vcov[[dpars[i]]]))
    if (dpars[i] == "drift") {
      tbl[[i]] <- cbind(cf, se, cf/se, 2 * pnorm(-abs(cf/se)))
      colnames(tbl[[i]]) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    } else {
      tbl[[i]] <- cbind(cf, se)
      colnames(tbl[[i]]) <- c("Estimate", "Std. Error")
    }
  }
  names(tbl) <- dpars
  object$coefficients <- tbl

  ## delete some slots
  # object$fitted.values <- object$terms <- object$model <- object$y <-
  #   object$x <- object$levels <- object$contrasts <- object$start <- NULL

  ## return
  class(object) <- "summary.ddm"
  return(object)
}


#' @rdname ddm-methods
#' @export
print.summary.ddm <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("\nCall:",
      deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)),
      "", sep = "\n")

  cat("DDM fit with",  length(x$coefficients), "estimated and",
      length(x$fixed_dpar), "fixed distributional parameters.\n")
  if (length(x$fixed_dpar) > 0) {
    cat("Fixed:",
        paste(names(x$fixed_dpar), x$fixed_dpar, sep = " = ", collapse = ", "),
        "\n")
  }
  cat("\n")

  cat(paste("drift coefficients (", x$link$drift$name,
            " link):\n", sep = ""))
  printCoefmat(x$coefficients$drift, digits = digits)
  cat("\n")

  if (NROW(x$coefficients) > 1) {
    for (par in seq_len(length(x$coefficients) - 1)) {
      cur_par <- names(x$coefficients)[par + 1]
      cat(paste(cur_par, " coefficients (", x$link[[cur_par]]$name,
                " link):\n", sep = ""))
      printCoefmat(x$coefficients[[cur_par]], digits = digits)
      cat("\n")
    }
  }

  if (!x$converged) {
    warning("optimisation algorithm indicates that model did not converge.\n",
            "Proceed with caution!\n", call. = FALSE)
  }

  invisible(x)
}

#' @rdname ddm-methods
#' @export
coef.ddm <- function(
    object,
    dpar = c("drift", "boundary", "ndt", "bias", "sv", "full"),
    ...) {
  cf <- object$coefficients
  dpar <- match.arg(dpar)

  if (dpar == "full") {
    all_names <- lapply(cf, names)
    dpar_names <- names(cf)
    cf <- unlist(cf)
    names(cf) <- paste0(
      unlist(all_names), " (",
      rep(dpar_names, vapply(all_names, length, 0L)), ")"
    )
    return(cf)
  } else {
    if (!(dpar %in% names(cf)))
      stop(dpar, " was fixed and not estimated.", call. = FALSE)
    cf[[dpar]]
  }
}

#' @rdname ddm-methods
#' @export
vcov.ddm <- function(
    object,
    dpar = c("drift", "boundary", "ndt", "bias", "sv"),
    ...) {
  cf <- object$vcov
  dpar <- match.arg(dpar)

  if (!(dpar %in% names(cf))) { # just picks out the right vcov
    stop(dpar, " was fixed and not estimated.", call. = FALSE)
  }
  return(cf[[dpar]])
}

#' @rdname ddm-methods
#' @export
model.frame.ddm <- function(formula, ...) {
  if (!is.null(formula$model)) return(formula$model)
  formula$terms <- formula$terms$full
  formula$call$formula <- formula$formula <- formula(formula$terms)
  NextMethod()
}

#' @rdname ddm-methods
#' @export
model.matrix.ddm <- function(
    object,
    dpar = c("drift", "boundary", "ndt", "bias", "sv"),
    ...) {
  dpar <- match.arg(dpar)
  if (!(dpar %in% object$dpar))
      stop(dpar, " was fixed and not estimated.", call. = FALSE)
  rval <- if (!is.null(object$x[[dpar]])) object$x[[dpar]]
    else model.matrix(object$terms[[dpar]], model.frame(object), contrasts = object$contrasts[[dpar]])
  return(rval)
}

#' @rdname ddm-methods
#' @export
terms.ddm <- function(
    x,
    dpar = c("drift", "boundary", "ndt", "bias", "sv"),
    ...) {
  dpar <- match.arg(dpar)
  if (!(dpar %in% x$dpar))
      stop(dpar, " was fixed and not estimated.", call. = FALSE)
  x$terms[[dpar]]
}

#' @rdname ddm-methods
#' @export
logLik.ddm <- function(object, ...) {
  structure(object$loglik,
            df = sum(sapply(object$coefficients, length)),
            class = "logLik")
}


#' @rdname ddm-methods
#' @export
update.ddm <- function(object, ...) {
  stop("no update method for ddm objects. Refit model using fddm::ddm()",
       call. = FALSE)
}
## note: update method could be structured similar to ddm call with one formula
## per estimated ddm parameter.
