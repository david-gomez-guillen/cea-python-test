"""A minimal cost-effectiveness model of first-line treatment for advanced cancer.

This file is meant to be read from top to bottom by someone learning
cost-effectiveness analysis (CEA). Every step of the model is a small function
with a docstring that explains both what it does and why.

The clinical question
---------------------
Patients are newly diagnosed with an advanced (metastatic) solid tumour. Standard
first-line care is chemotherapy. A new immunotherapy delays progression, but it
is far more expensive and is given for much longer. Part of the patients carry a
biomarker that predicts a better response to the immunotherapy. Which of these
strategies gives the best value for money?

    chemotherapy               everyone gets chemotherapy (the comparator)
    immunotherapy              everyone gets the immunotherapy
    biomarker_immunotherapy    everyone is tested; biomarker-positive patients get
                               the immunotherapy, the rest get chemotherapy

The model
---------
It is a Markov cohort model with the three states used in almost every oncology
economic evaluation:

    PFS          alive and progression-free, on first-line treatment
    Progressed   alive after the disease has progressed
    Dead         absorbing state

    Cycle:    one month
    Horizon:  ten years (120 cycles)

Rather than following patients one by one, a cohort model keeps track of the
share of a group of patients (the cohort) that is in each state, and moves those
shares from state to state once per cycle. Everyone starts in PFS.

For each strategy the model adds up what the cohort costs and how many
quality-adjusted life years (QALYs) it lives. Comparing those totals between
strategies is what a cost-effectiveness analysis does: the extra cost of a
strategy divided by the extra QALYs it brings is its incremental
cost-effectiveness ratio (ICER), in euros per QALY.

How to use it
-------------
From Python::

    import model
    results = model.simulate(model.STRATEGIES, **model.base_case())
    print(results["summary"])

From R (and from the THALASSA app), through thalassa_interface.R, which calls
``simulate()`` with reticulate. See example.py for a full worked example that
also computes the ICERs.
"""

from dataclasses import asdict, dataclass, fields
from typing import Sequence, Union

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Structure of the model
# ---------------------------------------------------------------------------

STATES = ("pfs", "progressed", "dead")
STRATEGIES = ("chemotherapy", "immunotherapy", "biomarker_immunotherapy")

MONTHS_PER_YEAR = 12
HORIZON_YEARS = 10
N_CYCLES = HORIZON_YEARS * MONTHS_PER_YEAR  # one cycle per month

# A probability that may change with time since the start of treatment: either
# one value for the whole horizon, or one value per year of follow-up.
TimeVarying = Union[float, Sequence[float]]


# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------

@dataclass
class Parameters:
    """Every input of the model, with its base-case value.

    The numbers are hypothetical but in a plausible range for an advanced solid
    tumour. All probabilities are per month (one cycle), and all costs are in
    euros: per month for ongoing costs, per event for one-off costs.
    """

    # ---- Natural history under chemotherapy --------------------------------
    p_progression: TimeVarying = 0.10
    """Monthly probability of progressing while in PFS on chemotherapy. Either one
    value or one value per year of follow-up (10 values). 0.10 gives a median
    progression-free survival of about six months."""

    p_death_pfs: TimeVarying = 0.01
    """Monthly probability of dying while progression-free (any cause)."""

    p_death_progressed: TimeVarying = 0.07
    """Monthly probability of dying after progression. 0.07 gives a median
    survival after progression of about nine and a half months."""

    # ---- Biomarker and immunotherapy effect --------------------------------
    p_biomarker_positive: float = 0.35
    """Share of patients who are biomarker-positive."""

    hr_progression_positive: float = 0.30
    """Hazard ratio of progression, immunotherapy versus chemotherapy, in
    biomarker-positive patients. Below 1 means the immunotherapy is better."""

    hr_progression_negative: float = 0.90
    """Hazard ratio of progression, immunotherapy versus chemotherapy, in
    biomarker-negative patients."""

    # ---- Treatment schedules -----------------------------------------------
    chemo_months: float = 4
    """Months of chemotherapy given to patients who stay progression-free."""

    immuno_max_months: float = 24
    """Maximum months of immunotherapy. It is stopped earlier at progression."""

    # ---- Costs ---------------------------------------------------------------
    cost_chemo: float = 2500
    """Monthly cost of chemotherapy (drug and administration)."""

    cost_immuno: float = 2000
    """Monthly cost of the immunotherapy (drug and administration, at the net
    price paid after confidential discounts)."""

    cost_biomarker_test: float = 400
    """One-off cost of testing one patient for the biomarker."""

    cost_pfs_care: float = 600
    """Monthly cost of follow-up care while progression-free (visits, scans)."""

    cost_progressed_care: float = 2000
    """Monthly cost of care after progression, including later lines of therapy."""

    cost_end_of_life: float = 8000
    """One-off cost of end-of-life care, paid when a patient dies."""

    # ---- Quality of life ----------------------------------------------------
    utility_pfs: float = 0.75
    """Utility of a year lived progression-free (1 = full health, 0 = dead)."""

    utility_progressed: float = 0.55
    """Utility of a year lived after progression."""

    # ---- General -------------------------------------------------------------
    discount: float = 0.03
    """Annual discount rate, applied to both costs and QALYs."""


def base_case() -> dict:
    """The base-case value of every parameter, as a dictionary.

    It can be unpacked straight into ``simulate()``::

        simulate(STRATEGIES, **base_case())
    """
    return asdict(Parameters())


# ---------------------------------------------------------------------------
# Building blocks
# ---------------------------------------------------------------------------
# Small helpers for the conversions every oncology model needs.

def probability_to_rate(p: float) -> float:
    """Convert a probability over one cycle into a rate (hazard).

    A probability says what share of the patients have the event within the
    cycle; a rate says how fast events happen at any moment. Assuming the rate is
    constant within the cycle, the share that does NOT have the event is
    exp(-rate), hence ``rate = -ln(1 - p)``.
    """
    return -np.log(1 - p)


def rate_to_probability(rate: float) -> float:
    """Convert a rate (hazard) back into a probability over one cycle."""
    return 1 - np.exp(-rate)


def apply_hazard_ratio(p: float, hazard_ratio: float) -> float:
    """Apply a hazard ratio to a probability.

    Trials report treatment effects as hazard ratios (HR), and a hazard ratio
    multiplies a rate, never a probability. Multiplying the probability directly
    is a common mistake: it gives the wrong answer, and with an HR above 1 it can
    even give a probability above 1. The correct way is to convert to a rate,
    multiply, and convert back.
    """
    return rate_to_probability(probability_to_rate(p) * hazard_ratio)


def discount_factor(cycle: int, annual_rate: float) -> float:
    """How much a cost or a QALY of a given cycle is worth today.

    Money and health count for less the further in the future they are. At an
    annual rate of 3%, what happens one year from now is worth 1 / 1.03 of the
    same today. Cycles are months, so month ``cycle`` is ``cycle / 12`` years
    away.
    """
    return 1 / (1 + annual_rate) ** (cycle / MONTHS_PER_YEAR)


def value_in_cycle(value: TimeVarying, cycle: int) -> float:
    """The value of a time-varying probability in a given cycle.

    A probability is either one number for the whole horizon, or a list with one
    number per year of follow-up, in which case the one of the current year is
    used. Year-specific values are what the calibration estimates (see
    thalassa_interface.R).
    """
    if np.ndim(value) == 0:
        return float(value)
    if len(value) != HORIZON_YEARS:
        raise ValueError(f"Expected one value per year ({HORIZON_YEARS}), got {len(value)}")
    year = cycle // MONTHS_PER_YEAR
    return float(value[year])


def transition_matrix(p_progression: float, p_death_pfs: float,
                      p_death_progressed: float) -> np.ndarray:
    """The transition matrix of one cycle.

    Row: the state a patient is in now. Column: the state one cycle later, in
    the order of STATES (PFS, Progressed, Dead).

                     PFS                 Progressed        Dead
        PFS        [ stay               progression       death in PFS       ]
        Progressed [ 0                  stay              death progressed   ]
        Dead       [ 0                  0                 1                  ]

    Each row adds up to 1, because everyone has to end up somewhere, so the
    probability of staying is whatever the other transitions leave. Nobody goes
    back to PFS once progressed, and nobody leaves Dead.
    """
    p_stay_pfs = 1 - p_progression - p_death_pfs
    p_stay_progressed = 1 - p_death_progressed
    matrix = np.array([
        [p_stay_pfs, p_progression,     p_death_pfs],
        [0.0,        p_stay_progressed, p_death_progressed],
        [0.0,        0.0,               1.0],
    ])
    # A matrix with negative entries means the inputs are not valid probabilities
    # (for instance, progression and death in PFS adding up to more than 1).
    if np.any(matrix < 0) or np.any(matrix > 1):
        raise ValueError("Invalid transition probabilities: "
                         f"progression={p_progression}, death in PFS={p_death_pfs}, "
                         f"death after progression={p_death_progressed}")
    return matrix


# ---------------------------------------------------------------------------
# One treatment arm
# ---------------------------------------------------------------------------

def simulate_arm(treatment: str, hazard_ratio: float, params: Parameters) -> dict:
    """Follow a cohort that all receives the same first-line treatment.

    Args:
        treatment: ``"chemotherapy"`` or ``"immunotherapy"``.
        hazard_ratio: hazard ratio of progression versus chemotherapy, applied
            for the whole horizon. Ignored for chemotherapy.
        params: the model inputs.

    Returns:
        A dictionary with the cohort trace (the share in each state at the start
        of each cycle, plus the end of the last one), and the discounted costs
        by component, the discounted QALYs and the undiscounted life years
        accumulated over the horizon.

    Costs and QALYs of a cycle are counted on the cohort at the start of that
    cycle, before it moves on (no half-cycle correction; with monthly cycles the
    difference is small).
    """
    if treatment not in ("chemotherapy", "immunotherapy"):
        raise ValueError(f"Unknown treatment: {treatment}")

    # The share of the cohort in each state: everyone starts progression-free.
    cohort = np.array([1.0, 0.0, 0.0])
    trace = [cohort]

    totals = {"cost_drug": 0.0, "cost_care": 0.0, "cost_end_of_life": 0.0,
              "qalys": 0.0, "life_years": 0.0, "pfs_life_years": 0.0}

    chemo_months = int(round(params.chemo_months))
    immuno_max_months = int(round(params.immuno_max_months))

    for cycle in range(N_CYCLES):
        pfs, progressed, _ = cohort

        # ---- 1. Transition probabilities of this cycle ----------------------
        # Under chemotherapy, the natural-history probabilities apply as they
        # are. The immunotherapy lowers the probability of progressing through
        # its hazard ratio. The model assumes a durable response: the effect
        # lasts for the whole horizon, even after the treatment is stopped.
        # Whether and how fast it wanes is one of the key assumptions of real
        # immunotherapy models.
        p_progression = value_in_cycle(params.p_progression, cycle)
        if treatment == "immunotherapy":
            p_progression = apply_hazard_ratio(p_progression, hazard_ratio)
        matrix = transition_matrix(p_progression,
                                   value_in_cycle(params.p_death_pfs, cycle),
                                   value_in_cycle(params.p_death_progressed, cycle))

        # ---- 2. Costs of this cycle -----------------------------------------
        # Drug costs are paid only for patients still on treatment: those in PFS
        # during the treatment period, which is fixed for chemotherapy and
        # capped for the immunotherapy. Care costs depend on the state. The
        # end-of-life cost is a one-off cost paid on the transition to Dead, so
        # it is counted on the patients who die during the cycle.
        on_chemo = treatment == "chemotherapy" and cycle < chemo_months
        on_immuno = treatment == "immunotherapy" and cycle < immuno_max_months
        drug = pfs * (params.cost_chemo * on_chemo + params.cost_immuno * on_immuno)
        care = pfs * params.cost_pfs_care + progressed * params.cost_progressed_care
        new_deaths = pfs * matrix[0, 2] + progressed * matrix[1, 2]
        end_of_life = new_deaths * params.cost_end_of_life

        # ---- 3. Health of this cycle ----------------------------------------
        # A month is 1/12 of a year, so a month in a state with utility u is
        # worth u / 12 QALYs. Life years are the same without the utility weight.
        qalys = (pfs * params.utility_pfs
                 + progressed * params.utility_progressed) / MONTHS_PER_YEAR

        # ---- 4. Add them up, discounted -------------------------------------
        # Life years are reported undiscounted, as survival is in a trial.
        factor = discount_factor(cycle, params.discount)
        totals["cost_drug"] += drug * factor
        totals["cost_care"] += care * factor
        totals["cost_end_of_life"] += end_of_life * factor
        totals["qalys"] += qalys * factor
        totals["life_years"] += (pfs + progressed) / MONTHS_PER_YEAR
        totals["pfs_life_years"] += pfs / MONTHS_PER_YEAR

        # ---- 5. Move the cohort one cycle forward ---------------------------
        # Multiplying the shares by the matrix sends each share where its row
        # says, which gives the shares at the start of the next cycle.
        cohort = cohort @ matrix
        trace.append(cohort)

    return {"trace": np.array(trace), **totals}


# ---------------------------------------------------------------------------
# One strategy
# ---------------------------------------------------------------------------

def strategy_arms(strategy: str, params: Parameters) -> list:
    """What each subgroup of patients receives under a strategy.

    A strategy is a mix of treatment arms. The cohort is split by biomarker
    status, and each subgroup is given a treatment. The biomarker is assumed to
    be predictive but not prognostic: it changes how well the immunotherapy
    works, but under chemotherapy both subgroups do equally well.

    Returns:
        A list of ``(share of the cohort, treatment, hazard ratio)``.
    """
    positive = params.p_biomarker_positive
    negative = 1 - positive
    if strategy == "chemotherapy":
        return [(1.0, "chemotherapy", 1.0)]
    if strategy == "immunotherapy":
        return [(positive, "immunotherapy", params.hr_progression_positive),
                (negative, "immunotherapy", params.hr_progression_negative)]
    if strategy == "biomarker_immunotherapy":
        return [(positive, "immunotherapy", params.hr_progression_positive),
                (negative, "chemotherapy", 1.0)]
    raise ValueError(f"Unknown strategy: {strategy}")


def simulate_strategy(strategy: str, params: Parameters) -> dict:
    """Simulate a strategy as the weighted mix of its treatment arms.

    Because the subgroups never mix, the result of the whole cohort is the
    average of the results of each subgroup, weighted by its share. The
    biomarker test, when the strategy uses one, is paid once for every patient
    at the start.
    """
    result = None
    for share, treatment, hazard_ratio in strategy_arms(strategy, params):
        arm = simulate_arm(treatment, hazard_ratio, params)
        weighted = {key: share * value for key, value in arm.items()}
        result = weighted if result is None else {
            key: result[key] + weighted[key] for key in result}

    tested = strategy == "biomarker_immunotherapy"
    result["cost_test"] = params.cost_biomarker_test if tested else 0.0
    result["cost_total"] = (result["cost_drug"] + result["cost_care"]
                            + result["cost_end_of_life"] + result["cost_test"])
    return result


# ---------------------------------------------------------------------------
# The whole model
# ---------------------------------------------------------------------------

def simulate(strategies: Sequence[str], **parameters) -> dict:
    """Run the model for the given strategies.

    Args:
        strategies: names of the strategies to simulate (see STRATEGIES).
        **parameters: the model inputs, named as the fields of ``Parameters``.
            Any input left out keeps its base-case value.

    Returns:
        A dictionary with:

        ``summary``
            One row per strategy, with the discounted total cost per patient
            (``C``) and the discounted QALYs per patient (``E``). These are the
            two numbers a cost-effectiveness analysis compares.
        ``outcomes``
            One row per strategy with the details behind them: undiscounted life
            years and progression-free life years, discounted QALYs, and the
            discounted cost split by component.
        ``pfs`` and ``os``
            For each strategy, the progression-free survival and the overall
            survival at the end of each year of follow-up, that is, the share of
            the cohort still in PFS and still alive. These are the curves a
            trial reports, and what the model is calibrated against.
        ``cohort_info``
            For each strategy, the full cohort trace: one row per cycle (plus the
            end of the last one), one column per state.
    """
    if isinstance(strategies, str):
        strategies = [strategies]
    params = Parameters(**parameters)

    # Only the probabilities may vary by year of follow-up.
    time_varying = {"p_progression", "p_death_pfs", "p_death_progressed"}
    for field in fields(Parameters):
        if field.name not in time_varying and np.ndim(getattr(params, field.name)) != 0:
            raise ValueError(f"{field.name} cannot vary by year")

    summary, outcomes = [], []
    pfs, os, cohort_info = {}, {}, {}
    year_ends = [MONTHS_PER_YEAR * year for year in range(1, HORIZON_YEARS + 1)]

    for strategy in strategies:
        result = simulate_strategy(strategy, params)
        trace = result["trace"]

        summary.append({"strategy": strategy,
                        "C": result["cost_total"],
                        "E": result["qalys"]})
        outcomes.append({"strategy": strategy,
                         **{key: value for key, value in result.items() if key != "trace"}})

        # Survival at the end of each year: the share in PFS, and the share in
        # PFS or Progressed (that is, not dead).
        pfs[strategy] = trace[year_ends, 0].tolist()
        os[strategy] = (trace[year_ends, 0] + trace[year_ends, 1]).tolist()
        cohort_info[strategy] = trace

    return {"summary": pd.DataFrame(summary),
            "outcomes": pd.DataFrame(outcomes),
            "pfs": pfs,
            "os": os,
            "cohort_info": cohort_info}
