# k_star_demo.ipynb

Companion notebook for the associated manuscript. Runtime < 3 minutes;
fourteen working cells cross-check every universal-n combinatorial
identity used by the paper against the Lean4 machine-verified
axiom report.

## Run locally

```
pip install -r notebook/requirements.txt
jupyter notebook notebook/k_star_demo.ipynb
```

Then select `Run All`.

## Run on Binder (no install)

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/kootru-repo/kstar-verify/v1.0.0-submission-submission?filepath=notebook%2Fk_star_demo.ipynb)

## Run on Colab

1. Open https://colab.research.google.com/
2. File > Open Notebook > GitHub tab.
3. Enter `kootru-repo/kstar-verify` and select `notebook/k_star_demo.ipynb`.

## Regenerate the notebook

The notebook is built programmatically from [`build_notebook.py`](build_notebook.py)
so that diffs to it are reviewable as plain Python:

```
python notebook/build_notebook.py
```
