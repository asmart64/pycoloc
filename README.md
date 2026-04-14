# pycoloc

`pycoloc` is a desktop app for two-channel fluorescence colocalization analysis.
It computes Manders-style metrics with automatic Costes thresholding and a
Costes randomization significance test, and provides interactive plots and ROI
selection for exploratory analysis.

The application is implemented in Python with a Tkinter GUI and Matplotlib
figures.

## What the app does

The app helps answer these practical questions for two aligned channels (Ch1,
Ch2):

- How much of Ch1 overlaps Ch2, and vice versa?
- What thresholds should be used objectively (Costes auto-threshold)?
- Is the observed channel correlation stronger than a randomized null model?

It includes:

- Loading of two single-channel images (TIFF, PNG, JPEG, BMP; 8-bit or 16-bit)
- Synthetic demo datasets for quick testing of different overlap regimes
- Interactive rectangular ROI drawing on Ch1, mirrored on Ch2
- Costes threshold estimation by scanning threshold pairs along regression
- Manders and related overlap metrics at both auto and manual thresholds
- Costes randomization test (block shuffling, PSF-like block size)
- Histogram and scatter visual diagnostics

## Main features

### 1) Image loading and preprocessing

- Accepts TIFF via `tifffile` (preferred for microscopy formats)
- Falls back to Pillow for common image formats
- Handles common dimensionality cases by reducing to 2D grayscale

### 2) ROI-driven analysis

- Draw and edit ROI on Channel 1 with interactive handles
- ROI coordinates are mirrored to Channel 2 display
- Clear ROI to revert analysis to the full image

### 3) Costes automatic thresholding

Given flattened pixel vectors `c1` and `c2`:

1. Fit linear regression `Ch2 = a * Ch1 + b`
2. Scan candidate `T1` values from high to low
3. Pair each `T1` with `T2 = clip(a*T1 + b)`
4. Compute Pearson `r` on pixels below both thresholds
5. Select the first threshold where background correlation becomes
	 non-positive (`r <= 0`) with minimum foreground support constraints

This gives robust automatic thresholds while avoiding degenerate tiny-tail
foreground solutions.

### 4) Metrics reported

The results panel reports:

- Costes thresholds: `T1`, `T2`
- Manders-style coefficients: `M1`, `M2`, `tM1`, `tM2`
- Pearson correlation above thresholds
- Manders overlap coefficient `R` plus `k1`, `k2`
- Costes randomization outputs:
	- observed `r`
	- randomized mean and standard deviation
	- significance percentage
	- Monte Carlo p-value

### 5) Visualization

- Two intensity histograms (one per channel) with threshold markers
- Pixel scatter plot with threshold crosshairs and regression line
- Costes curve: background `r` versus `T1`, including zero-crossing context

### 6) Manual threshold mode

You can type `T1` and `T2` values to immediately recalculate all metrics and
refresh plots, which is useful for sensitivity checks and reporting.

## Demo presets

The built-in synthetic generator includes presets such as:

- High overlap
- Partial overlap
- Offset structures
- Low overlap
- Randomization demo
- Randomization control

These are intended for sanity checks and method demonstration.

## Dependencies

Install required packages:

```bash
pip install numpy scipy matplotlib tifffile Pillow
```

## Run

From the repository root:

```bash
python colocalization_app.py
```

## Typical workflow

1. Load Channel 1 and Channel 2 (or load a demo set)
2. Optionally draw/refine ROI on Channel 1
3. Click "Run Costes + Analyze"
4. Inspect thresholds, coefficients, and randomization significance
5. Optionally tweak manual thresholds and compare outcomes

## Notes

- Both channels must have the same spatial dimensions.
- For TIFF export and optimal TIFF loading, `tifffile` is recommended.
- Randomization block size is controlled by the PSF field in the UI.

## Repository layout

- `colocalization_app.py`: main GUI and analysis logic
- `requirements.txt`: Python dependencies

## License

This repository currently includes a `LICENSE` file at the root.
