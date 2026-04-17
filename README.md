# pycoloc

`pycoloc` is a desktop app for two-channel fluorescence colocalization analysis.
It computes Manders-style metrics with automatic Costes thresholding and a
Costes randomization significance test, and provides interactive plots and ROI
selection for exploratory analysis.

The application is implemented in Python with a Tkinter GUI and Matplotlib
figures. The UI is organized into three separate windows:

- **Analyzer Controls** — all buttons, options, and status display
- **Colocalization Analyzer** — analysis plots (histograms, Costes scatter, results)
- **Channel Images** — channel display with ROI drawing tools

## What the app does

The app helps answer these practical questions for two aligned channels (Ch1,
Ch2):

- How much of Ch1 overlaps Ch2, and vice versa?
- What thresholds should be used objectively (Costes auto-threshold)?
- Is the observed channel correlation stronger than a randomized null model?

It includes:

- Loading of two single-channel images (TIFF, PNG, JPEG, BMP; 8-bit or 16-bit)
- Synthetic demo datasets for quick testing of different overlap regimes
- Interactive rectangular or freehand (lasso) ROI drawing on Ch1, mirrored on Ch2
- Costes threshold estimation by scanning threshold pairs along regression
- Manders and related overlap metrics at both auto and manual thresholds
- Costes randomization test (block shuffling, PSF-like block size)
- Histogram and scatter visual diagnostics
- Auto-stretch display toggle for contrast enhancement (robust 2–98th percentile)
- Hot-pixel pre-processing (optional median-based despiking)
- PDF summary export

## Main features

### 1) Image loading and preprocessing

- Accepts TIFF via `tifffile` (preferred for microscopy formats)
- Falls back to Pillow for common image formats
- Handles common dimensionality cases by reducing to 2D grayscale
- Optional hot-pixel suppression using local median filter

### 2) ROI-driven analysis

- **Rectangle mode**: draw and interactively resize a rectangular ROI on Channel 1; mirrored on Channel 2
- **Lasso mode**: freehand polygon ROI for arbitrary region selection
- Clear ROI to revert analysis to the full image
- Costes randomization is disabled for lasso ROIs (non-rectangular mask incompatible with block shuffling)

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

- Two intensity histograms (one per channel) with threshold markers and vertical grids
- Pixel scatter plot with threshold crosshairs and regression line
- Costes curve: background `r` versus `T1`, including zero-crossing context

### 6) Manual threshold mode

Type `T1` and `T2` values to immediately recalculate all metrics and
refresh plots — useful for sensitivity checks and reporting. Editing the PSF
field automatically re-runs the Costes algorithm.

### 7) Auto-stretch display

Toggle "Auto-stretch display" to enhance visual contrast using a robust
2nd–98th percentile mapping. The analysis always uses original pixel values.

## Demo presets

The built-in synthetic generator includes presets:

| Preset | Description |
|---|---|
| High overlap | Strongly shared puncta with minimal independent signal |
| Partial overlap | Mixed shared, offset, and channel-specific structures |
| Offset structures | Similar objects with systematic spatial offsets |
| Low overlap | Mostly independent objects with weak diffuse correlation |
| Pure bleedthrough | Ch2 dominated by linear bleed-through from Ch1 |
| Randomization demo | High observed r, low shuffled r (strong shared structure) |
| Randomization control | Mostly independent channels; shuffled r ≈ observed r |

Demo sets can also be exported as uint16 TIFF files via "Export Demo TIFFs".

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
2. Optionally draw/refine a rectangular or lasso ROI on Channel 1
3. Click **Run Costes + Analyze**
4. Inspect thresholds, coefficients, and randomization significance in the Results tab
5. Optionally tweak manual thresholds to compare outcomes
6. Save a PDF summary report if needed

## Notes

- Both channels must have the same spatial dimensions.
- For TIFF export and optimal TIFF loading, `tifffile` is recommended.
- Randomization block size is controlled by the PSF field (in pixels).
- Lasso ROI does not support Costes randomization.
- All three windows must remain open for full functionality; closing any one exits the app.

## Repository layout

- `colocalization_app.py`: main GUI and analysis logic (single-file app)
- `requirements.txt`: Python dependencies

## License

This repository currently includes a `LICENSE` file at the root.
