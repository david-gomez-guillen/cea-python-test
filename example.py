"""Worked example: run the base case and work out which strategy is cost-effective.

Run it with::

    python example.py

It prints the health outcomes and costs of each strategy, and then the
incremental analysis: the strategies sorted by cost, the dominated ones ruled out,
and the incremental cost-effectiveness ratio (ICER) of each remaining strategy
against the previous one.
"""

import pandas as pd

import model

# Willingness-to-pay thresholds: the most the health system is prepared to pay
# for one QALY. They differ between countries; €25,000 is the reference used
# here, and €50,000 shows how the decision changes under a more generous one.
THRESHOLDS = (25_000, 50_000)


def incremental_analysis(summary: pd.DataFrame) -> pd.DataFrame:
    """Compute ICERs along the efficiency frontier.

    Comparing every strategy with the cheapest one is not enough when there are
    more than two. The usual steps are:

    1. Sort the strategies from cheapest to most expensive.
    2. Rule out strongly dominated strategies: those that cost more and give
       fewer (or equal) QALYs than another one. Nobody should choose them.
    3. Compute the ICER of each strategy against the previous one that is still
       in: (cost difference) / (QALY difference).
    4. Rule out extendedly dominated strategies: those whose ICER is higher than
       the ICER of the next, more effective, strategy. A mix of their neighbours
       would give more QALYs for the same money. Repeat until ICERs only grow.

    The strategies left form the efficiency frontier. The cost-effective one is
    the most effective whose ICER is below the threshold.
    """
    table = summary.sort_values("C").reset_index(drop=True)
    table["status"] = ""

    def frontier():
        return table[table["status"] == ""].index

    # Step 2: strong dominance.
    for i in table.index:
        cheaper = table.loc[: i - 1]
        if ((cheaper["E"] >= table.loc[i, "E"]) & (cheaper["status"] == "")).any():
            table.loc[i, "status"] = "dominated"

    # Steps 3 and 4: ICERs and extended dominance.
    while True:
        rows = frontier()
        table["ICER"] = float("nan")
        for previous, current in zip(rows[:-1], rows[1:]):
            table.loc[current, "ICER"] = (
                (table.loc[current, "C"] - table.loc[previous, "C"])
                / (table.loc[current, "E"] - table.loc[previous, "E"]))
        icers = table.loc[rows, "ICER"]
        # An ICER higher than the next one along the frontier.
        extended = icers[icers > icers.shift(-1)]
        if extended.empty:
            break
        table.loc[extended.index[0], "status"] = "extendedly dominated"

    return table


def cost_effective(table: pd.DataFrame, threshold: float) -> str:
    """The most effective strategy of the frontier whose ICER is below the threshold."""
    affordable = table[(table["status"] == "") & ~(table["ICER"] > threshold)]
    return affordable.iloc[-1]["strategy"]


def main():
    pd.set_option("display.width", 120)
    results = model.simulate(model.STRATEGIES, **model.base_case())

    print("Health outcomes and costs per patient (costs and QALYs discounted)\n")
    print(results["outcomes"].set_index("strategy").T.round(3).to_string())

    print("\nProgression-free (PFS) and overall survival (OS) at the end of each year\n")
    years = [f"Year {y}" for y in range(1, model.HORIZON_YEARS + 1)]
    for curve in ("pfs", "os"):
        print(pd.DataFrame(results[curve], index=years).head(5)
              .round(3).to_string(header=[f"{curve.upper()} {s}" for s in model.STRATEGIES]))
        print()

    print("Incremental analysis\n")
    table = incremental_analysis(results["summary"])
    print(table.round({"C": 0, "E": 3, "ICER": 0}).to_string(index=False))

    print()
    for threshold in THRESHOLDS:
        print(f"At €{threshold:,} per QALY, the cost-effective strategy is: "
              f"{cost_effective(table, threshold)}")


if __name__ == "__main__":
    main()
