
fit_nlminb <- function(init, 
                       objective, gradient, 
                       lower, upper, 
                       control, 
                       ...) ## ... arguments are ignored
{
  out <- nlminb(start = init, 
                objective = objective,
                gradient = gradient,
                lower = lower, upper = upper,
                control = control)
  list(
    coefficients = out$par,
    loglik = -out$objective,
    converged = out$convergence == 0,
    optim = out
  )
}
