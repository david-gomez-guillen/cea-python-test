# THALASSA interface of the advanced cancer model.
#
# The model itself is written in Python, in model.py. This file only describes it
# to the app (parameters, strategies, strata, calibration...) and calls it through
# reticulate, which runs Python from R and converts the values between the two:
# an R list becomes a Python dict, a numeric vector a list, and a pandas data
# frame comes back as an R data.frame.
#
# reticulate uses the Python of the RETICULATE_PYTHON environment variable, which
# the app sets to the Python of the model's conda environment (environment.yml).

library(reticulate)

# Import model.py as a Python module. Its functions are then called as
# model$simulate(...), model$base_case(), and so on.
model <- import_from_path('model', path='.')

get.overview <- function() {
  # Markdown shown in the Overview tab. The text lives in overview.md. Everything
  # described there is derived from model.py and from the rest of this interface,
  # so it must be kept in sync with them.
  return(paste(readLines('overview.md'), collapse='\n'))
}

get.model.settings <- function() {
  return(list(cost.unit='€', effectiveness.unit='QALY', run.initial.guess=TRUE))
}

get.strategies <- function() {
  # The names must match model.STRATEGIES in model.py. The descriptions and
  # attributes describe strategy_arms() and simulate_strategy() in model.py, and
  # must be kept in sync with them.
  return(list(
    list(
      name='chemotherapy',
      display.name='Chemotherapy for all',
      description='The comparator. Every patient receives chemotherapy for chemo.months months while progression-free.',
      attributes=list(immunotherapy='None', chemotherapy='All patients', biomarker.test='No')
    ),
    list(
      name='immunotherapy',
      display.name='Immunotherapy for all',
      description='Every patient receives the immunotherapy until progression, for at most immuno.max.months months. It multiplies the rate of progression by hr.progression.positive in biomarker-positive patients and by hr.progression.negative in the rest.',
      attributes=list(immunotherapy='All patients', chemotherapy='None', biomarker.test='No')
    ),
    list(
      name='biomarker_immunotherapy',
      display.name='Biomarker-guided immunotherapy',
      description='Every patient is tested for the biomarker, at cost.biomarker.test each. Biomarker-positive patients (p.biomarker.positive) receive the immunotherapy and the rest receive chemotherapy.',
      attributes=list(immunotherapy='Biomarker-positive', chemotherapy='Biomarker-negative', biomarker.test='Yes')
    )
  ))
}

get.strategy.attributes <- function() {
  # What describes the strategies, shown as columns of the Strategies tab. Only
  # who receives the immunotherapy is drawn on the base case plot, as the shape
  # of the point.
  return(list(
    immunotherapy=list(
      label='Immunotherapy for',
      plot='shape',
      values=c(None='x', `All patients`='circle', `Biomarker-positive`='square')
    ),
    chemotherapy='Chemotherapy for',
    biomarker.test='Biomarker test'
  ))
}

get.strata <- function() {
  # Years since the start of first-line treatment. They match the ten years of the
  # horizon in model.py (HORIZON_YEARS).
  return(paste('Year', 1:10))
}

# ---- Parameters ---------------------------------------------------------------

# The base value of every parameter is read from model.py, so that the Python
# model stays the one place where they are set. Python names use underscores and
# R names dots: p_death_pfs in Python is p.death.pfs here.
python.name <- function(name) gsub('.', '_', name, fixed=TRUE)
base.value <- function(name) model$base_case()[[python.name(name)]]

# Constraints: each returns TRUE when the value is acceptable, or the message
# shown for it when it is not. A constraint that ties two parameters together is
# declared on both, so that the app highlights whichever of them the user
# changed. A parameter split into strata arrives as a list, hence the unlist().
is.probability <- function(par.name, params) {
  x <- unlist(params[[par.name]])
  if (any(x < 0 | x > 1)) sprintf('%s must be between 0 and 1', par.name) else TRUE
}
death.pfs.below.progressed <- function(par.name, params) {
  if (any(unlist(params[['p.death.pfs']]) >= unlist(params[['p.death.progressed']])))
    'Death while progression-free must be less likely than after progression'
  else TRUE
}
utility.pfs.above.progressed <- function(par.name, params) {
  if (params[['utility.pfs']] <= params[['utility.progressed']])
    sprintf('Utility while progression-free (%g) must be higher than after progression (%g)',
            params[['utility.pfs']], params[['utility.progressed']])
  else TRUE
}
positive.benefit.more <- function(par.name, params) {
  if (params[['hr.progression.positive']] >= params[['hr.progression.negative']])
    sprintf(paste('Biomarker-positive patients must benefit more than biomarker-negative ones',
                  '(hazard ratios %g and %g)'),
            params[['hr.progression.positive']], params[['hr.progression.negative']])
  else TRUE
}

# distribution is the one the PSA draws the parameter from. The treatment
# schedules and the discount rate declare none, which keeps them out of the PSA.
# min.value and max.value bound the values the user can enter: 0 to 1 for
# probabilities, utilities and the discount rate, and above 0 for costs.
parameter <- function(name, display.name, class, constraints=NULL, distribution=NULL,
                      min.value=NULL, max.value=NULL) {
  p <- list(name=name, display.name=display.name, base.value=base.value(name), class=class)
  if (!is.null(constraints)) p$constraints <- constraints
  if (!is.null(distribution)) p$distribution <- distribution
  if (!is.null(min.value)) p$min.value <- min.value
  if (!is.null(max.value)) p$max.value <- max.value
  return(p)
}

get.parameters <- function() {
  return(list(
    parameter('p.progression', 'Monthly probability of progression on chemotherapy', 'Natural history',
              list(is.probability), distribution='beta',
              min.value=0, max.value=1),
    parameter('p.death.pfs', 'Monthly probability of death while progression-free', 'Natural history',
              list(is.probability, death.pfs.below.progressed), distribution='beta',
              min.value=0, max.value=1),
    parameter('p.death.progressed', 'Monthly probability of death after progression', 'Natural history',
              list(is.probability, death.pfs.below.progressed), distribution='beta',
              min.value=0, max.value=1),

    parameter('p.biomarker.positive', 'Share of biomarker-positive patients', 'Biomarker and immunotherapy',
              list(is.probability), distribution='beta',
              min.value=0, max.value=1),
    parameter('hr.progression.positive', 'Hazard ratio of progression, immunotherapy vs chemotherapy, biomarker-positive', 'Biomarker and immunotherapy',
              list(positive.benefit.more), distribution='lognormal'),
    parameter('hr.progression.negative', 'Hazard ratio of progression, immunotherapy vs chemotherapy, biomarker-negative', 'Biomarker and immunotherapy',
              list(positive.benefit.more), distribution='lognormal'),

    parameter('chemo.months', 'Months of chemotherapy', 'Treatment schedules'),
    parameter('immuno.max.months', 'Maximum months of immunotherapy (stopped earlier at progression)', 'Treatment schedules'),

    parameter('cost.chemo', 'Monthly cost of chemotherapy', 'Costs', distribution='gamma', min.value=0),
    parameter('cost.immuno', 'Monthly cost of immunotherapy', 'Costs', distribution='gamma', min.value=0),
    parameter('cost.biomarker.test', 'Cost of a biomarker test', 'Costs', distribution='gamma', min.value=0),
    parameter('cost.pfs.care', 'Monthly cost of care while progression-free', 'Costs', distribution='gamma', min.value=0),
    parameter('cost.progressed.care', 'Monthly cost of care after progression, including later lines', 'Costs', distribution='gamma', min.value=0),
    parameter('cost.end.of.life', 'One-off cost of end-of-life care', 'Costs', distribution='gamma', min.value=0),

    parameter('utility.pfs', 'Utility while progression-free', 'Utilities',
              list(utility.pfs.above.progressed), distribution='beta',
              min.value=0, max.value=1),
    parameter('utility.progressed', 'Utility after progression', 'Utilities',
              list(utility.pfs.above.progressed), distribution='beta',
              min.value=0, max.value=1),

    parameter('discount', 'Annual discount rate for costs and QALYs', 'General',
              min.value=0, max.value=1)
  ))
}

# ---- Overview tab ---------------------------------------------------------------

get.model.states <- function() {
  # State diagram shown in the Overview tab, next to the text of overview.md. It
  # describes transition_matrix() in model.py and must be kept in sync with it.
  return(list(
    title='Model states',
    description='The three states of the cohort and the monthly transitions between them. Hover a state or a transition for what drives it.',
    nodes=data.frame(
      id=c('pfs', 'progressed', 'dead'),
      label=c('Progression-free', 'Progressed', 'Dead'),
      description=c(
        'Alive and progression-free, on first-line treatment. The whole cohort starts here.',
        'Alive after the disease has progressed, receiving later lines of treatment.',
        'Absorbing state. The end-of-life cost is paid on entry.'
      )
    ),
    edges=data.frame(
      from=c('pfs', 'pfs', 'progressed'),
      to=c('progressed', 'dead', 'dead'),
      label=c('Progression', 'Death', 'Death'),
      description=c(
        'p.progression, its rate multiplied by hr.progression.positive or hr.progression.negative under immunotherapy.',
        'p.death.pfs',
        'p.death.progressed'
      )
    )
  ))
}

get.code.sample <- function() {
  # Code shown in the Code panel of the Overview tab. model.py is written to be
  # read, so it is shown as it is.
  return(list(
    `model.py`=list(code=readLines('model.py'), language='python')
  ))
}

# ---- Simulation -----------------------------------------------------------------

to.python <- function(pars) {
  # The app gives a parameter split into strata as a list keyed by stratum. The
  # model expects one value per year, in order, so the list is turned into a
  # plain vector in the order of the strata. reticulate passes a vector of
  # several values as a Python list and a single value as a number.
  args <- lapply(pars, function(value) {
    if (is.list(value)) unname(unlist(value[get.strata()])) else value
  })
  names(args) <- python.name(names(pars))
  return(args)
}

plain.data.frame <- function(df) {
  attr(df, 'pandas.index') <- NULL
  rownames(df) <- NULL
  return(df)
}

run.simulation <- function(strategies, pars) {
  results <- do.call(model$simulate, c(list(as.list(strategies)), to.python(pars)))
  # pandas data frames come back as R data frames that still point to their
  # pandas index. The app copies results between processes, and a pointer to a
  # Python object cannot be copied, so only the plain data frame is kept.
  results$summary <- plain.data.frame(results$summary)
  results$outcomes <- plain.data.frame(results$outcomes)
  return(results)
}

# ---- Calibration ----------------------------------------------------------------

get.calibration.schemes <- function() {
  return(list(
    standard=list(
      description='Progression-free survival of the chemotherapy arm of a pivotal trial',
      parameters='p.progression',
      target=list(
        # Progression-free survival at the end of each year of follow-up in the
        # chemotherapy arm of a (hypothetical) trial. The trial followed patients
        # for five years only, so years 6 to 10 are not calibrated and keep the
        # base value of p.progression.
        PFS=list(
          `Year 1`=.30,
          `Year 2`=.12,
          `Year 3`=.07,
          `Year 4`=.05,
          `Year 5`=.04
        )
      ),
      strata=paste('Year', 1:5),
      initial_guess=rep(.10, 5),
      error_function=calibration.error,
      latent_space_training_set=generate.training.dataset,
      evaluation.plots=NULL
    )))
}

calibration.error <- function(pars, target) {
  calibration.strategy <- 'chemotherapy'
  # The target is a named list with one entry per stratum, so it is flattened into
  # a named numeric vector to match the simulated values by stratum name.
  target.pfs <- unlist(target$PFS)
  result <- tryCatch({
    results <- run.simulation(calibration.strategy, pars)
    pfs <- unlist(results$pfs[[calibration.strategy]])
    names(pfs) <- get.strata()
    # Only the strata present in the target contribute to the error.
    list(error=sum((pfs[names(target.pfs)] - target.pfs)^2),
         output=list(PFS=pfs))
  }, error=function(e) {
    # Parameter sets the model rejects (e.g. probabilities above 1) get an
    # infinite error, which the app records as a failed attempt.
    pfs <- rep(NA, length(get.strata()))
    names(pfs) <- get.strata()
    list(error=Inf, output=list(PFS=pfs))
  })
  return(result)
}

generate.training.dataset <- function(initial_guess, n, ...) {
  # Samples for the latent space methods: the initial guess scaled by a uniform
  # random factor per parameter, kept between 0 and 0.5 so that every sample is a
  # valid monthly probability of progression.
  f.pars <- list(...)
  variation <- f.pars$variation

  n_params <- length(initial_guess)
  dataset <- matrix(NA, nrow=n, ncol=n_params)
  for (i in 1:n) {
    factors <- runif(n_params, min=1-variation, max=1+variation)
    dataset[i, ] <- pmin(.5, pmax(0, initial_guess * factors))
  }
  dataset <- dataset[sample(nrow(dataset)), ]
  return(dataset)
}
