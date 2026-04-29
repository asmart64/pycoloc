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
- ROI mask import/export for file-based region definitions
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
- **Mask file workflow**: load a binary mask image as ROI or save the current ROI to a new mask file
- Clear ROI to revert analysis to the full image
- Costes randomization supports rectangular, lasso, and file-loaded mask ROIs

### 3) Costes automatic thresholding

Given flattened pixel vectors `c1` and `c2`:

1. Fit linear regression `Ch2 = a * Ch1 + b`
2. Follow the Coloc2 stepping rule: step on `Ch1` when `-1 < a < 1`, otherwise step on `Ch2`
3. Round each candidate threshold pair to integer image levels, matching Coloc2 behavior
4. Compute Pearson `r` on the Costes background set `c1 < T1 OR c2 < T2`
5. Stop using the Coloc2-style SimpleStepper criteria when background `r` becomes very small, non-finite, or starts increasing again. An importante difference has been introduced with respect to the latter: the incremental ratio computed as r/r_previous is prone to statistical noise and that prompted unreasonable early stops. Now the search for T1 and T2 stops only if the relative increase of r with respect to the running minimum r value is larger then some percentage( e.g. 1%).

This aligns the automatic threshold search more closely with Coloc2 instead of the earlier simplified scan.

### Coloc2 alignment changes

The recent Costes updates were introduced specifically to match Coloc2 more closely:

- Threshold stepping now follows the active channel rule used by Coloc2 (`Ch1` for slopes between `-1` and `1`, otherwise `Ch2`)
- Background correlation is evaluated on the union mask `c1 < T1 OR c2 < T2`, matching Coloc2 `ThresholdMode.Below`
- Threshold candidates are rounded to integer image intensities before evaluation
- The Costes curve now stores and displays both members of the threshold pair, so the `Pearson vs T1` plot also exposes the paired `T2` values
- Costes randomization now works for lasso ROIs by shuffling blocks inside the lasso bounding box and evaluating Pearson only on pixels inside the lasso mask

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
- Costes curve: background `r` versus `T1`, with paired `T2` values shown on the top axis

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
## venv setup
cd to the app folder and create a virtual environment with a python >=3.12:
```bash
cd pycoloc-x.x.x
python3.12 -m venv ./
```
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
2. Optionally draw/refine a rectangular or lasso ROI on Channel 1, or load a mask file as ROI
3. Click **Run Costes + Analyze**
4. Inspect thresholds, coefficients, and randomization significance in the Results tab
5. Optionally tweak manual thresholds to compare outcomes
6. Save a PDF summary report if needed

## Notes

- Both channels must have the same spatial dimensions.
- For TIFF export and optimal TIFF loading, `tifffile` is recommended.
- Randomization block size is controlled by the PSF field (in pixels).
- For lasso and file-loaded mask ROIs, Costes randomization is evaluated on masked pixels only.
- All three windows must remain open for full functionality; closing any one exits the app.

## Repository layout

- `colocalization_app.py`: main GUI and analysis logic (single-file app)
- `requirements.txt`: Python dependencies

## License

This repository currently includes a `LICENSE` file at the root.
