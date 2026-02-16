# matRockmag – Python Port

Python translation of the **matRockmag** MATLAB toolbox by Robert E. Kopp (2008).
Parses `.rmg` rock-magnetism data files and produces a full analysis dashboard,
statistics tables, curve fits, and simulations.

---

## File map

| Python file | Original MATLAB | Purpose |
|---|---|---|
| `run_rmg_analysis.py` | `RMGFileFullAnalysis.m` | **Entry point** |
| `rmg_import.py` | `RMGImport.m` | Parse `.rmg` file → `RmgData` object |
| `rmg_curves.py` | `RmgSIRMCurve`, `RmgARMCurve`, `RmgIRMBackfieldCurve`, `RmgRRMCurve`, `RmgLowrieFullerCurves`, `RmgFullerCurves`, `RmgExtractAFOfStep`, `RmgAFLevelFind`, `moving` | Curve extraction |
| `rmg_stats.py` | `RmgStats.m` | Rock-mag statistics |
| `rmg_plots.py` | All 8 plot functions + `RmgDataFullAnalysis` | Plotting + 3×3 dashboard |
| `sgg.py` | `SGG.m` | Skewed Generalized Gaussian PDF |
| `fit_sgg.py` | `fitSGG.m`, `fitSGGComps.m`, `fitSGGMulti.m` | SGG curve fitting |
| `rmg_curve_fits.py` | `RmgSIRMDerivativeCurveFits.m`, `RmgSIRMDerivativeCurveFitsSGG.m`, `RmgSIRMDerivativeCurveFitComps.m`, `RmgAFDerivativeCurveFits.m` | IRM/AF derivative fitting |
| `rmg_write_tables.py` | `RmgStatsWriteTable.m`, `RmgIRMFitsWriteTable.m` | Export `.asc` tables |
| `rmg_simulate.py` | `RmgIRMSimulate.m`, `RmgSimulateSwitchingFieldDistribution.m`, `RmgFractionalARM.m`, `RmgARMCurveSubtract.m`, `CODICAErrorSmoothing.m` | Simulations & specialised analysis |

---

## Requirements

```
numpy>=1.23
scipy>=1.9
matplotlib>=3.6
```

Install: `pip install -r requirements.txt`

---

## Quick start

```bash
python run_rmg_analysis.py            # interactive: prompts for file path(s)
python run_rmg_analysis.py sample.rmg # direct argument
```

The script shows the 3×3 dashboard and offers to save a PNG and stats `.asc` table.

---

## Scripting examples

```python
# ── Load and plot ─────────────────────────────────────────────────────────────
from rmg_import import rmg_import
from rmg_plots  import rmg_data_full_analysis
import matplotlib.pyplot as plt

data = rmg_import('sample.rmg')
fig  = rmg_data_full_analysis(data)
plt.savefig('output.png', dpi=150, bbox_inches='tight')

# ── Stats table ───────────────────────────────────────────────────────────────
from rmg_write_tables import rmg_stats_write_table
rmg_stats_write_table([data], 'rockmagstats-summary')  # → .asc

# ── SGG curve fitting (auto component selection) ─────────────────────────────
from rmg_curves      import rmg_sirm_curve
from rmg_curve_fits  import rmg_sirm_derivative_curve_fits_sgg

sirm = rmg_sirm_curve(data)[0]
fits = rmg_sirm_derivative_curve_fits_sgg(sirm)
print(f"IRM best fit: {fits['IRMfit']['bestfitSGGn']} SGG component(s)")

# ── Gaussian component fitting ────────────────────────────────────────────────
from rmg_curve_fits import rmg_sirm_derivative_curve_fits
fitted = rmg_sirm_derivative_curve_fits(sirm)
print(f"IRM Gaussians: {fitted.IRM.logderivbestfitgaussians}")

# ── IRM simulation ────────────────────────────────────────────────────────────
import numpy as np
from rmg_simulate import rmg_irm_simulate
B   = np.linspace(0, 0.3, 100)
Ban = np.array([0.05, 0.10, 0.20])
lw  = np.array([0.01, 0.02, 0.03])
result = rmg_irm_simulate(B, Ban, lw)

# ── CODICA error smoothing ────────────────────────────────────────────────────
from rmg_simulate import codica_error_smoothing
smoothed = codica_error_smoothing(fields, irm_values)['rescaled']
```

---

## Dashboard panels (3×3)

| | Col 1 | Col 2 | Col 3 |
|---|---|---|---|
| **Row 1** | IRM acq + AF demag | dIRM/dB derivative | ARM acquisition |
| **Row 2** | Lowrie-Fuller (ARM+IRM) | Lowrie-Fuller (ARM only) | Fuller plot |
| **Row 3** | RRM | Backfield IRM | Stats text box |

---

## Not ported

- MATLAB GUI (`matRockmag.fig`) — use scripts or Jupyter notebook instead
- `MPMSImport.m` — different instrument format
- `RmgStatDepthProfiles.m` — trivial to add (calls `rmg_stats` + matplotlib)
- `RmgLowrieFullerDerivativePlot.m` — add on request

---

## License

Original MATLAB © 2008 Robert E. Kopp, GNU GPL v.3. Python port 2025, same license.
