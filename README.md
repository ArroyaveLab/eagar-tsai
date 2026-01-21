# ET Model (Eagar–Tsai)

This project implements the Eagar–Tsai moving heat source model to estimate melt pool dimensions for a scanning laser/beam over a semi‑infinite solid. The temperature field is computed from a 1D integral and evaluated numerically (optionally using a compiled C integrand for speed). Melt pool dimensions are extracted from the liquidus isotherm.

Attribution: The C integrand implementation is based on a reformulation by Sasha Rubenchik (LLNL, 2015).

Reference: T. W. Eagar and N.-S. Tsai, “Temperature Fields Produced by Traveling Distributed Heat Sources,” Welding Journal (Research Supplement), December 1983, pp. 346‑s–354‑s (see `original_ET_paper.pdf`).

Refactor note: Brent Vela refactored this code on January 21, 2026.

## Inputs (per row)

Beam/process:
- `Velocity_m/s` (scan speed, v) [m/s]
- `Power` (laser power, P) [W]
- `Beam_diameter_m` (beam diameter, 2σ) [m]
- `Absorptivity` (A) [unitless]

Material (liquidus properties):
- `T_liquidus` [K]
- `thermal_cond_liq` (k) [W/(m·K)]
- `Density_kg/m3` (ρ) [kg/m^3]
- `Cp_J/kg` (cp) [J/(kg·K)]

Fixed/implicit:
- Ambient temperature `t0 = 300 K`
- Domain size and spatial resolution (defaults in `compute_melt_pool`)

## Outputs (per row)

- `melt_length` [m]
- `melt_width` [m]
- `melt_depth` [m]
- `melt_length_um` [µm]
- `melt_width_um` [µm]
- `melt_depth_um` [µm]
- `peakT` [K]
- `minT` [K]

## Equations (as implemented)

Thermal diffusivity:
```
alpha = k / (rho * cp)
```

Non‑dimensional parameter:
```
p = alpha / (v * sigma)
```

Prefactor:
```
Ts = (A*P) / (pi*(k/alpha)*sqrt(pi*alpha*v*(sigma^3)))
```

Temperature field at (x,y,z):
```
T = t0 + Ts * ∫_0^∞ f(t, x, y, z, p) dt
```

Integrand:
```
f(t,x,y,z,p) = 1 / ((4 p t + 1) * sqrt(t))
               * exp( -z^2/(4t) - ((y^2 + (x - t)^2)/(4 p t + 1)) )
```

## Assumptions

- Semi‑infinite solid; domain only used for numerical evaluation.
- Constant material properties at liquidus.
- No melt flow, vaporization, or latent heat effects.
- Gaussian heat source with constant absorptivity.
- Steady‑state moving source; integral is evaluated numerically.

## Platform Build Notes (eagar_tsai_integrand)

This project uses a small C helper (`eagar_tsai_integrand.c`) that can be compiled into a shared library to speed up the Eagar–Tsai integration. The Python code will use the compiled library if it is present, and fall back to the slow interpreted integrator if not.

### Windows

#### Option A: MSVC (Visual Studio Developer Command Prompt)

1. Open the **Developer Command Prompt for VS**.
2. From the project folder, run:
   ```bat
   build_windows_dll.bat
   ```

#### Option B: MinGW-w64 (gcc)

1. Install MinGW-w64 and ensure `gcc` is on your PATH.
2. From the project folder, run:
   ```bat
   build_windows_dll.bat
   ```

#### PowerShell alternative

If you prefer PowerShell, run:
```powershell
.\build_windows_dll.ps1
```

This produces `libeagar_tsai_integrand.dll` in the same directory as `et_melt_pool_script.py`.

### Linux

From the project folder:
```bash
gcc -O3 -fPIC -shared -o libeagar_tsai_integrand.so eagar_tsai_integrand.c
```

This produces `libeagar_tsai_integrand.so` in the same directory as `et_melt_pool_script.py`.

### macOS

From the project folder:
```bash
clang -O3 -fPIC -shared -o libeagar_tsai_integrand.dylib eagar_tsai_integrand.c
```

This produces `libeagar_tsai_integrand.dylib` in the same directory as `et_melt_pool_script.py`.

### Notes

- The Python code looks for these files by name in the working directory:
  - Windows: `libeagar_tsai_integrand.dll`
  - Linux: `libeagar_tsai_integrand.so`
  - macOS: `libeagar_tsai_integrand.dylib`
- The exported C symbol name that Python loads is `eagar_tsai_integrand`.
- If the compiled library is not found, the code will fall back to the interpreted integrator (much slower).
- If you change the C function name or file name, rebuild the shared library before running Python.


## Validation Against Old Predictions

If you have legacy ET outputs (e.g., `ET_width_old`, `ET_depth_old`) in the input sheet, you can validate the new run by:

1. Run the model to generate fresh `melt_width`, `melt_depth`, and `melt_length`.
2. Compare the new outputs to the old columns using summary error metrics and parity plots.

Example (Python):
```python
import pandas as pd
from et_melt_pool_script import compute_melt_pool

# Load input with old ET columns (renamed to *_old)
df = pd.read_excel("et_input_data_example.xlsx")

# Map material columns used by the model
# (Adjust names if your sheet uses different headers)
df["T_liquidus"] = df["PROP LT (K)"]
df["thermal_cond_liq"] = df["PROP LT THCD (W/(mK))"]
df["Density_kg/m3"] = df["PROP RT Density (kg/m3)"]
df["Cp_J/kg"] = df["PROP LT C (J/(kg K))"]
df["Beam_diameter_m"] = df["Beam Diam (m)"]

# Run ET model
out = compute_melt_pool(df, workers=1, chunk_size=50)

# Simple error metrics vs. legacy columns
for col_new, col_old in [("melt_width", "ET_width_old"), ("melt_depth", "ET_depth_old")]:
    if col_old in out.columns:
        diff = out[col_new] - out[col_old]
        mae = diff.abs().mean()
        rmse = (diff.pow(2).mean()) ** 0.5
        print(col_new, "MAE:", mae, "RMSE:", rmse)
```

Recommended checks:
- Parity plots for width/depth.
- MAE/RMSE summary and max absolute error.
- Stratify errors by process parameters (power, speed) to spot systematic drift.
