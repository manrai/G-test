# assume beta(1,1)
# computes posterior for beta-binomial given num observations & sample size
# draws trials number of samples from this posterior

sample.beta <- function(freq, n, trials=10000) {
  obs_success <- freq * n
  obs_failure <- n - obs_success
  out <- sapply(1:length(freq), function(i)
    rbeta(n = trials, shape = obs_success[i] + 1, shape2 = obs_failure[i] + 1) %>%
      setNames(1:trials)
  )
  colnames(out) <- names(freq)
  rownames(out) <- NULL
  return(out)
}