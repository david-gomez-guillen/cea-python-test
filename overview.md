# Advanced cancer model

A minimal cost-effectiveness model of first-line treatment for a hypothetical advanced
(metastatic) solid tumour, written as a teaching example for an oncology audience new
to cost-effectiveness analysis. The model is written in Python (`model.py`, shown in
the Code panel below), and this app runs it from R through `reticulate`.

It answers the following question: in patients newly diagnosed with advanced cancer,
is a new immunotherapy good value for money compared with standard chemotherapy, and
does it become better value if only biomarker-positive patients receive it?

## Strategies

- **Chemotherapy for all** (`chemotherapy`): the comparator. Chemotherapy is given
  for `chemo.months` months to patients who remain progression-free.
- **Immunotherapy for all** (`immunotherapy`): everyone receives the immunotherapy
  until progression, for at most `immuno.max.months` months.
- **Biomarker-guided immunotherapy** (`biomarker_immunotherapy`): everyone is tested
  for the biomarker. Biomarker-positive patients receive the immunotherapy and the
  rest receive chemotherapy.

## Model structure

The model is a three-state Markov cohort model, the structure used by most oncology
economic evaluations:

- **Progression-free (PFS)**: alive, with the disease controlled by first-line
  treatment. The whole cohort starts here.
- **Progressed**: alive after the disease has progressed. Patients do not return to
  PFS.
- **Dead**: absorbing state, reached from both of the other states.

Cycles last one month and the cohort is followed for ten years (120 cycles). Costs
and QALYs of a cycle are counted on the cohort at the start of the cycle, before it
moves on, and both are discounted.

Some assumptions are worth keeping in mind when reading the results:

- **The hazard ratio applies to the rate, not the probability.** Trials report the
  effect of a treatment as a hazard ratio (HR). The model converts the monthly
  probability of progression into a rate, multiplies it by the HR and converts it
  back.
- **Durable response.** The immunotherapy lowers the rate of progression for the
  whole horizon, even after it is stopped at `immuno.max.months`. In real appraisals,
  whether and how fast the effect wanes is often the assumption that matters most.
- **The biomarker is predictive, not prognostic.** It changes how well the
  immunotherapy works (`hr.progression.positive` compared with
  `hr.progression.negative`), but under chemotherapy positive and negative patients do
  equally well.
- **Survival after progression is the same in every strategy.** The immunotherapy
  gains life years only by delaying progression.

## Parameters

| Parameter | Base value | Class | Meaning |
|---|---|---|---|
| `p.progression` | 0.10 | Natural history | Monthly probability of progressing while progression-free on chemotherapy. Can take one value per year of follow-up. |
| `p.death.pfs` | 0.01 | Natural history | Monthly probability of dying while progression-free. |
| `p.death.progressed` | 0.07 | Natural history | Monthly probability of dying after progression. |
| `p.biomarker.positive` | 0.35 | Biomarker and immunotherapy | Share of patients who are biomarker-positive. |
| `hr.progression.positive` | 0.30 | Biomarker and immunotherapy | Hazard ratio of progression, immunotherapy versus chemotherapy, in biomarker-positive patients. |
| `hr.progression.negative` | 0.90 | Biomarker and immunotherapy | The same hazard ratio in biomarker-negative patients. |
| `chemo.months` | 4 | Treatment schedules | Months of chemotherapy for patients who remain progression-free. |
| `immuno.max.months` | 24 | Treatment schedules | Maximum months of immunotherapy, which stops earlier at progression. |
| `cost.chemo` | 2,500 | Costs | Monthly cost of chemotherapy. |
| `cost.immuno` | 2,000 | Costs | Monthly cost of the immunotherapy, at its net price. |
| `cost.biomarker.test` | 400 | Costs | Cost of testing one patient, paid for everyone under `biomarker_immunotherapy`. |
| `cost.pfs.care` | 600 | Costs | Monthly cost of follow-up while progression-free. |
| `cost.progressed.care` | 2,000 | Costs | Monthly cost of care after progression, including later lines of therapy. |
| `cost.end.of.life` | 8,000 | Costs | One-off cost of end-of-life care, paid on death. |
| `utility.pfs` | 0.75 | Utilities | Utility of a year lived progression-free. |
| `utility.progressed` | 0.55 | Utilities | Utility of a year lived after progression. |
| `discount` | 0.03 | General | Annual discount rate for costs and QALYs. |

The values are hypothetical but in a plausible range. A constant `p.progression` of
0.10 gives a median progression-free survival of about six months on chemotherapy,
and `p.death.progressed` of 0.07 a median survival after progression of about nine
and a half months. The class is only used to group the parameters in the interface.

## Strata

Results are reported by year since the start of treatment, from `Year 1` to `Year 10`.
The probabilities of progression and death can be given one value per year, which is
what the calibration uses.

## Outputs

- `summary`: one row per strategy with the cost (`C`, the discounted total cost per
  patient over the horizon) and the effect (`E`, the discounted QALYs per patient).
  The ratio of their differences between two strategies is the ICER, in € per QALY.
- `outcomes`: the details behind the summary for each strategy: undiscounted life
  years and progression-free life years, discounted QALYs, and the discounted cost
  split into drugs, care, end-of-life care and testing.
- `pfs` and `os`: progression-free and overall survival at the end of each year, that
  is, the share of the cohort still progression-free and still alive. These are the
  curves a trial reports.

In the base case, biomarker-guided immunotherapy costs about €9,100 more than
chemotherapy and gains about 0.27 QALYs, an ICER of about €34,000 per QALY. That is
above a willingness-to-pay threshold of €25,000 per QALY, so at that threshold
chemotherapy remains the cost-effective choice. Immunotherapy for all gains only
another 0.03 QALYs over the biomarker-guided strategy, at about €204,000 per extra
QALY: in biomarker-negative patients the immunotherapy costs a lot and does little.

## Calibration

The `standard` scheme calibrates `p.progression` against the progression-free survival
of the chemotherapy arm of a (hypothetical) pivotal trial, using the `chemotherapy`
strategy. It estimates one value per year for the five years the trial followed
patients, starting from 0.10 in each year.

| | Year 1 | Year 2 | Year 3 | Year 4 | Year 5 |
|---|---|---|---|---|---|
| PFS target | 0.30 | 0.12 | 0.07 | 0.05 | 0.04 |

The error is the sum of squared differences between the simulated and the target PFS
over the five calibrated years. Parameter sets the model rejects, such as
probabilities that add up to more than 1, get an infinite error.

PFS at the end of a year depends on the probabilities of progression in all the
earlier years, so the years cannot be fitted one at a time. The targets are met
exactly, with monthly probabilities of progression that fall from about 0.085 in the
first year to about 0.008 in the fifth. This is the long tail of patients whose
disease stays controlled for years, which a constant probability cannot reproduce.

Years 6 to 10 are beyond the follow-up of the trial, so they are not calibrated and
keep the base value of 0.10. That is an extrapolation choice, and it shows on the
calibration plot as a sudden drop in PFS after year 5. Such choices matter: with the
calibrated values, the ICER of biomarker-guided immunotherapy falls to about €27,000
per QALY, still above a €25,000 threshold. It falls to about €23,000, below the
threshold, if the fifth-year probability is carried forward instead of reverting to
0.10. The extrapolation choice alone changes the decision.

For the latent space methods, the training set is drawn by scaling the initial guess
by a uniform random factor per parameter, kept between 0 and 0.5.

## Things to try

- Find the price of the immunotherapy (`cost.immuno`) at which biomarker-guided
  immunotherapy becomes cost-effective at €25,000 per QALY. How does that price
  change at a threshold of your own country?
- Shorten `immuno.max.months`. How do costs and QALYs change, given that the effect
  is assumed to persist?
- Make the biomarker less discriminating by moving `hr.progression.positive` and
  `hr.progression.negative` closer together. When is testing no longer worth it?
- Run the deterministic sensitivity analysis and see which parameters drive the
  ICER.
