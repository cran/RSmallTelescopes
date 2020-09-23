#' Simulate Power 
#'
#' Estimate statistical power for point estimate of effect size plus the lower 
#' and upper bounds of a confidence interval.
#'
#' @param data Dataset (matrix).
#' @param n.original The sample size of the original analysis (scalar).
#' @param analysis Function to produce a p value and an effect size estimate.
#' @param sims.CI The number of simulated samples used to construct CI (scalar); default = 10,000.
#' @param CI.level The confidence level of the interval (scalar); default = .90.
#' @param sims.power The number of samples to be simulated (scalar); default = 10,000.
#' @param sims.samples Number of samples analyzed at upper/lower bounds of CI (scalar); default = 11.
#' @param seed Allows randomly generated numbers to be reproducible (scalar); default = 1.
#' @return Displays statistical power for point estimate of an effect size plus the lower and upper 
#' bounds of a confidence interval. List contains the following components:
#' \item{n.replication}{The sample size of the replication analysis.}
#' \item{n.original}{The sample size of the original analysis.}
#' \item{sims.CI}{The number of simulated samples used to construct CI.}
#' \item{CI.level}{The confidence level of the interval.}
#' \item{sims.power}{The number of samples simulated.}
#' \item{sims.samples}{Number of samples analyzed at upper/lower bounds of CI.}
#' \item{es.estimate}{Point estimate of effect size.}
#' \item{es.power}{Estimated power for the point estimate of effect size.}
#' \item{CI.lower.estimate}{Effect size estimate at the lower bound of the CI.}
#' \item{CI.lower.power}{Estimated power for the lower bound of the CI.}
#' \item{CI.upper.estimate}{Effect size estimate at the upper bound of the CI.}
#' \item{CI.upper.power}{Estimated power for the upper bound of the CI.}
#' @examples 
#' # create or import dataset
#'  example.data <- matrix(rnorm(50), 25, 2)
#'
#' # conduct empirical small telescopes analysis
#'  SimulatePower(
#'    data = example.data, 
#'    n.original = 10, 
#'    analysis = function(data) {
#'      corr <- cor.test(data[,1], data[,2])
#'      return(list(effect.size = corr$estimate, p.value = corr$p.value))
#'    }, 
#'    sims.CI = 100, 
#'    sims.power = 100)
#' @import stats
#' @export
SimulatePower <- function(data, n.original, analysis, sims.CI = 10000, CI.level = .90,
                          sims.power = 10000, sims.samples = 11, seed = 1) {
  n.rep <- dim(data)[1]
  set.seed(seed)
  point.est <- analysis(data)$effect.size
  power.point <- EstimatePower(data, n.original, sims.power, analysis)
  effect.size <- rep(0, sims.CI)
  for (i in 1:sims.CI) {
    set.seed(i * seed)
    data.sample <- data[sample(1:n.rep, n.rep, replace = T), ]
    effect.size[i] <- analysis(data.sample)$effect.size
  }
  ordered.samples <- (1:sims.CI)[sort.list(effect.size)]
  tail.size <- (1 - CI.level) / 2
  effectsize.lci <- rep(0, sims.samples)
  power.lci <- rep(0, sims.samples)
  power.uci <- rep(0, sims.samples)
  for (i in 1:sims.samples) {
    sample.ci <- ordered.samples[round(sims.CI * tail.size + 1 + 
                 (i - sims.samples + (sims.samples - 1) / 2))]
    if (i == ((sims.samples - 1) / 2)) {   
      effectsize.lci <- effect.size[sample.ci]
    }
    set.seed(sample.ci * seed)
    data.ci <- data[sample(1:n.rep, n.rep, replace = T), ]
    power.lci[i] <- EstimatePower(data.ci, n.original, sims.power, analysis)
  }
  for (i in 1:sims.samples) {
    sample.ci <- ordered.samples[round(sims.CI * (1 - tail.size) +
                 (i - sims.samples + (sims.samples - 1) / 2))]
    if (i == ((sims.samples - 1) / 2)) {
      effectsize.uci <- effect.size[sample.ci]
    }
    set.seed(sample.ci * seed)
    data.ci <- data[sample(1:n.rep, n.rep, replace = T), ]
    power.uci[i] <- EstimatePower(data.ci, n.original, sims.power, analysis)
  }
  return(list(
    n.replication = n.rep,
    n.original = n.original,
    sims.CI = sims.CI, 
    CI.level = CI.level, 
    sims.power = sims.power, 
    sims.samples = sims.samples, 
    es.estimate = point.est,
    es.power = power.point, 
    CI.lower.estimate = effectsize.lci, 
    CI.lower.power = median(power.lci), 
    CI.upper.estimate = effectsize.uci, 
    CI.upper.power = median(power.uci)))
}

#' Estimate Power 
#'
#' Estimate statistical power of an effect size parameter by simulation using 
#' original sample size.
#' 
#' @param data Dataset (matrix).
#' @param n.original The sample size of the original analysis (scalar).
#' @param sims.power The number of samples to be simulated (scalar).
#' @param analysis Function to produce a p value and an effect size estimate. 
#' @return Power estimate generated through simulation (scalar).
#' @examples 
#' # create or import dataset
#'  example.data <- matrix(rnorm(50), 25, 2)  
#' 
#' # estimate statistical power
#'  EstimatePower(
#'    data = example.data, 
#'    n.original = 10, 
#'    analysis = function(data) {
#'      corr <- cor.test(data[,1], data[,2])
#'      return(list(effect.size = corr$estimate, p.value = corr$p.value))
#'    }, 
#'    sims.power = 100)
#' @export 
EstimatePower <- function(data, n.original, sims.power, analysis) {
  reject <- 0
  n <- dim(data)[1]
  for (i in 1:sims.power) {
    data.sample <- data[sample(1:n, n.original, replace = T), ]
    p.value <- analysis(data.sample)$p.value
    reject <- reject + (p.value < .05)
  }
  power <- reject / sims.power
  return(power)
}