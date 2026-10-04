# Advanced cancer CEA model (Python)

A minimal cost-effectiveness model of first-line treatment for advanced cancer,
written in Python as a teaching example for oncologists who are new to
cost-effectiveness analysis. It compares chemotherapy, immunotherapy for all, and
biomarker-guided immunotherapy with a three-state Markov cohort model
(progression-free, progressed, dead).

The model runs on its own in Python, and in the THALASSA app
through `thalassa_interface.R`, which calls the Python code with `reticulate`.

## Files

| File | What it is |
|---|---|
| `model.py` | The model. Read it from top to bottom: every step is a small function with a docstring explaining it. |
| `example.py` | A worked example: runs the base case and carries out the incremental analysis (dominance, ICERs, threshold). |
| `overview.md` | Description of the model, its parameters and its calibration, shown in the Overview tab of THALASSA. |
| `thalassa_interface.R` | Describes the model to THALASSA and calls `model.py` through `reticulate`. |
| `environment.yml` | Conda environment THALASSA builds for the model (R, reticulate, Python, numpy, pandas). |
| `tests/test_model.py` | Sanity checks that any correct cohort model should pass. |

## Running it in Python

```bash
pip install -r requirements.txt
python example.py                    # base case and ICERs
python -m unittest discover tests    # sanity checks
```

To use the model in your own code:

```python
import model

results = model.simulate(model.STRATEGIES, **model.base_case())
print(results["summary"])            # cost (C) and QALYs (E) per strategy

# Any input can be changed by name.
cheaper = model.simulate(model.STRATEGIES, **(model.base_case() | {"cost_immuno": 2000}))
```

## Running it in R

```r
Sys.setenv(RETICULATE_PYTHON = "/path/to/python")   # a Python with numpy and pandas
source("thalassa_interface.R")

params <- get.parameters()
pars <- setNames(lapply(params, `[[`, "base.value"), sapply(params, `[[`, "name"))
run.simulation(c("chemotherapy", "immunotherapy", "biomarker_immunotherapy"), pars)$summary
```

Parameter names use dots in R (`cost.immuno`) and underscores in Python
(`cost_immuno`). The interface converts between them, and reads the base values
from `model.py`, so `model.py` is the only place where they are set.
