#!/usr/bin/env python3
"""
Colocalization Analyzer
=======================
Estimates Mander's Overlap Coefficients (M1, M2, tM1, tM2) for two-channel
fluorescence images, using the Costes automatic thresholding algorithm.

Features
--------
- Load two single-channel images (TIFF, PNG, JPEG, BMP; 8-bit or 16-bit).
- Generate synthetic demo image pairs with multiple colocalization regimes
    for testing the analysis workflow.
- Draw a rectangular ROI on the Channel-1 display; analysis is restricted to
  that region (or the full image when no ROI is defined).
- Costes algorithm: iterates the threshold pair (T1, T2) down the linear
  regression line of Ch2 vs Ch1 until Pearson's r drops to ≤ 0.
- Intensity histograms with threshold markers.
- Scatter plot coloured by colocalization status and Pearson's-r-vs-threshold
  curve that illustrates how the Costes threshold is found.
- Manual threshold input boxes for live recalculation of the Manders coefficients.
- Save a PDF summary report with metrics and current plots.

Dependencies
------------
    pip install numpy scipy matplotlib tifffile Pillow
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from pathlib import Path
from datetime import datetime
import subprocess

import numpy as np
import warnings
warnings.filterwarnings("ignore")

import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Polygon, Rectangle
from matplotlib.path import Path as MplPath
from matplotlib.widgets import LassoSelector, RectangleSelector

from scipy import stats, ndimage

# ── optional image-loading libraries ──────────────────────────────────────────
try:
    import tifffile
    _HAS_TIFFFILE = True
except ImportError:
    _HAS_TIFFFILE = False

try:
    from PIL import Image
    _HAS_PIL = True
except ImportError:
    _HAS_PIL = False

# ── custom colormaps: black → colour ──────────────────────────────────────────
def _linear_cmap(r, g, b, name):
    cdict = {
        "red":   [(0.0, 0.0, 0.0), (1.0, r, r)],
        "green": [(0.0, 0.0, 0.0), (1.0, g, g)],
        "blue":  [(0.0, 0.0, 0.0), (1.0, b, b)],
    }
    return mcolors.LinearSegmentedColormap(name, cdict)

_CYAN    = _linear_cmap(0, 1, 1, "cyan_cm")
_MAGENTA = _linear_cmap(1, 0, 1, "magenta_cm")

# ── colour palette ─────────────────────────────────────────────────────────────
BG     = "#1e1e2e"
BG2    = "#2b2b3b"
FG     = "#cdd6f4"
YELLOW = "#f9e2af"
CYAN   = "#89dceb"
MAGENTA = "#f38ba8"
GREEN  = "#a6e3a1"
GRAY   = "#6c7086"


def _detect_app_version() -> str:
    """Return a short build identifier, preferring the current git commit."""
    try:
        repo_dir = Path(__file__).resolve().parent
        result = subprocess.run(
            ["git", "-C", str(repo_dir), "rev-parse", "--short", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
            timeout=1.0,
        )
        commit = result.stdout.strip()
        if commit:
            return f"git-{commit}"
    except Exception:
        pass
    return "dev"


APP_VERSION = _detect_app_version()


# ══════════════════════════════════════════════════════════════════════════════
class ColocalizationApp:
    """Main application window."""

    def __init__(self, root: tk.Tk):
        self.root = root
        self.root.title(f"Colocalization Analyzer — Mander's Coefficients ({APP_VERSION})")
        self.root.configure(bg=BG2)
        # Use a screen-aware startup size so the window never opens off-screen.
        self._set_initial_window_geometry()

        # application state
        self.ch1: np.ndarray | None = None
        self.ch2: np.ndarray | None = None
        self.roi: tuple | None = None          # (x1, y1, x2, y2) pixel coords
        self.roi_polygon: np.ndarray | None = None
        self.roi_mask: np.ndarray | None = None
        self._roi_kind: str | None = None
        self._roi_mode_var = tk.StringVar(value="Rectangle")
        self._roi_active = True
        self._demo_presets = {
            "High overlap": "Strongly shared puncta with mild independent signal.",
            "Partial overlap": "Mixed shared, offset, and channel-specific structures.",
            "Offset structures": "Similar objects with systematic spatial offsets.",
            "Low overlap": "Mostly independent objects with weak diffuse correlation.",
            "Pure bleedthrough": "Channel 2 is dominated by linear bleed-through from channel 1.",
            "Randomization demo": "Aligned mesoscale structure: high observed r, low shuffled r.",
            "Randomization control": "Mostly independent channels: shuffled r similar to observed r.",
        }
        self._last_costes_slope: float | None = None
        self._last_costes_intercept: float | None = None
        self._last_costes_curve_t1 = np.array([])
        self._last_costes_curve_t2 = np.array([])
        self._last_costes_curve_r = np.array([])
        self._roi_overlay_after_id = None
        self._controls_win: tk.Toplevel | None = None
        self._images_win: tk.Toplevel | None = None
        self._tooltip_win: tk.Toplevel | None = None
        self._tooltip_after_id = None
        self._despike_var = tk.BooleanVar(value=False)
        self._autostretch_var = tk.BooleanVar(value=False)
        self._histeq_var = tk.BooleanVar(value=False)
        self._clahe_var = tk.BooleanVar(value=False)
        self._orthreg_var = tk.BooleanVar(value=True)
        self._costes_debug_log_var = tk.BooleanVar(value=False)
        self._costes_stop_pct_var = tk.StringVar(value="1.0")

        self._build_ui()
        self.root.after(0, self._ensure_window_visible)
        self.root.after(0, self._position_controls_window)
        self.root.after(0, self._position_images_window)
        self.root.protocol("WM_DELETE_WINDOW", self._on_close)

    def _set_initial_window_geometry(self):
        """Choose a default size that fits the current display and center it."""
        screen_w = int(self.root.winfo_screenwidth())
        screen_h = int(self.root.winfo_screenheight())

        # Keep margins so title bar and edges remain visible on smaller screens.
        # Reduced width by ~40%, increased height to show all controls including manual thresholds.
        target_w = min(660, max(560, int((screen_w - 120) * 0.60)))
        target_h = min(900, max(640, int((screen_h - 160) * 0.90)))

        x = max(0, (screen_w - target_w) // 2)
        y = max(0, (screen_h - target_h) // 2)
        self.root.geometry(f"{target_w}x{target_h}+{x}+{y}")

        # Keep usability while allowing a substantially smaller main window.
        self.root.minsize(560, 640)

    def _ensure_window_visible(self):
        """Clamp the mapped window to the visible screen so toolbar is never clipped."""
        self.root.update_idletasks()

        screen_w = int(self.root.winfo_screenwidth())
        screen_h = int(self.root.winfo_screenheight())
        w = int(self.root.winfo_width())
        h = int(self.root.winfo_height())
        x = int(self.root.winfo_x())
        y = int(self.root.winfo_y())

        max_x = max(0, screen_w - w)
        max_y = max(0, screen_h - h)
        clamped_x = min(max(x, 0), max_x)
        clamped_y = min(max(y, 0), max_y)

        if clamped_x != x or clamped_y != y:
            self.root.geometry(f"{w}x{h}+{clamped_x}+{clamped_y}")

    def _position_controls_window(self):
        """Position the controls window near the main window and keep it visible."""
        if self._controls_win is None or not self._controls_win.winfo_exists():
            return

        self.root.update_idletasks()
        self._controls_win.update_idletasks()

        screen_w = int(self.root.winfo_screenwidth())
        screen_h = int(self.root.winfo_screenheight())

        main_x = int(self.root.winfo_x())
        main_y = int(self.root.winfo_y())
        main_w = int(self.root.winfo_width())

        cw_w = int(self._controls_win.winfo_width())
        cw_h = int(self._controls_win.winfo_height())

        x = main_x + max(0, (main_w - cw_w) // 2)
        y = max(0, main_y + 24)

        x = min(max(0, x), max(0, screen_w - cw_w))
        y = min(max(0, y), max(0, screen_h - cw_h))
        self._controls_win.geometry(f"+{x}+{y}")

    def _position_images_window(self):
        """Position the images window near the main window and keep it visible."""
        if self._images_win is None or not self._images_win.winfo_exists():
            return

        self.root.update_idletasks()
        self._images_win.update_idletasks()

        screen_w = int(self.root.winfo_screenwidth())
        screen_h = int(self.root.winfo_screenheight())

        main_x = int(self.root.winfo_x())
        main_y = int(self.root.winfo_y())
        main_w = int(self.root.winfo_width())

        iw_w = int(self._images_win.winfo_width())
        iw_h = int(self._images_win.winfo_height())

        x = main_x + main_w + 16
        y = max(0, main_y)

        if x + iw_w > screen_w:
            x = max(0, main_x - iw_w - 16)

        x = min(max(0, x), max(0, screen_w - iw_w))
        y = min(max(0, y), max(0, screen_h - iw_h))
        self._images_win.geometry(f"+{x}+{y}")

    # ── UI construction ────────────────────────────────────────────────────────

    def _build_ui(self):
        self._apply_theme()
        self._build_controls_window()
        self._build_images_window()

        # ── main window now hosts analysis panels only ────────────────────────
        right = ttk.Frame(self.root)
        right.pack(fill=tk.BOTH, expand=True, padx=6, pady=4)
        self._build_right_panel(right)

    def _build_images_window(self):
        """Build a dedicated resizable window for channel image display."""
        win = tk.Toplevel(self.root)
        self._images_win = win
        win.title("Channel Images")
        win.configure(bg=BG2)
        win.resizable(True, True)
        win.protocol("WM_DELETE_WINDOW", self._on_close)

        outer = ttk.LabelFrame(win, text="Images  (draw ROI on left image)", padding=4)
        outer.pack(fill=tk.BOTH, expand=True, padx=6, pady=6)
        self._build_image_panel(outer)

        win.update_idletasks()
        screen_w = int(win.winfo_screenwidth())
        screen_h = int(win.winfo_screenheight())
        w = min(max(640, int(win.winfo_reqwidth()) + 20), max(520, screen_w - 120))
        h = min(max(420, int(win.winfo_reqheight()) + 20), max(360, screen_h - 140))
        win.geometry(f"{w}x{h}")
        win.minsize(620, 400)

    def _build_controls_window(self):
        """Build a compact secondary window containing all menu controls."""
        win = tk.Toplevel(self.root)
        self._controls_win = win
        win.title("Analyzer Controls")
        win.configure(bg=BG2)
        win.resizable(False, False)
        win.transient(self.root)
        win.protocol("WM_DELETE_WINDOW", self._on_close)

        outer = ttk.Frame(win, padding=8)
        outer.pack(fill=tk.BOTH, expand=True)

        self._demo_preset_var = tk.StringVar(value="Partial overlap")
        self._status = tk.StringVar(value="Load two single-channel images to begin.")

        # Row 0
        btn_load_ch1 = ttk.Button(outer, text="Load Channel 1", command=self.load_ch1, width=16)
        btn_load_ch1.grid(row=0, column=0, padx=3, pady=2)
        self._attach_tooltip(btn_load_ch1, "Open grayscale image for channel 1")

        btn_load_ch2 = ttk.Button(outer, text="Load Channel 2", command=self.load_ch2, width=16)
        btn_load_ch2.grid(row=0, column=1, padx=3, pady=2)
        self._attach_tooltip(btn_load_ch2, "Open grayscale image for channel 2")

        self._demo_preset_box = ttk.Combobox(
            outer,
            textvariable=self._demo_preset_var,
            values=list(self._demo_presets.keys()),
            state="readonly",
            width=20,
        )
        self._demo_preset_box.grid(row=0, column=2, padx=3, pady=2)
        self._attach_tooltip(self._demo_preset_box, "Select a synthetic demo scenario")

        btn_demo_load = ttk.Button(outer, text="Load Demo Set", command=self.load_demo_set, width=14)
        btn_demo_load.grid(row=0, column=3, padx=3, pady=2)
        self._attach_tooltip(btn_demo_load, "Generate and load demo channels")

        btn_demo_export = ttk.Button(outer, text="Export Demo TIFFs", command=self.export_demo_set, width=16)
        btn_demo_export.grid(row=0, column=4, padx=3, pady=2)
        self._attach_tooltip(btn_demo_export, "Save demo channels as TIFF files")

        # Row 1
        ttk.Label(outer, text="ROI mode:").grid(row=1, column=0, padx=(3, 1), pady=2, sticky=tk.E)
        roi_mode_frame = ttk.Frame(outer)
        roi_mode_frame.grid(row=1, column=1, padx=3, pady=2, sticky=tk.W)
        rb_rect = ttk.Radiobutton(
            roi_mode_frame,
            text="Rectangle",
            value="Rectangle",
            variable=self._roi_mode_var,
            command=self._on_roi_mode_change,
        )
        rb_rect.pack(side=tk.LEFT, padx=(0, 8))
        self._attach_tooltip(rb_rect, "Use rectangular ROI with resize handles")

        rb_lasso = ttk.Radiobutton(
            roi_mode_frame,
            text="Lasso",
            value="Lasso",
            variable=self._roi_mode_var,
            command=self._on_roi_mode_change,
        )
        rb_lasso.pack(side=tk.LEFT)
        self._attach_tooltip(rb_lasso, "Use freehand ROI drawing")

        self._roi_btn = ttk.Button(
            outer,
            text="Cancel ROI" if self._roi_active else "Draw ROI",
            command=self.toggle_roi,
            width=14,
        )
        self._roi_btn.grid(row=1, column=2, padx=3, pady=2, sticky=tk.W)
        self._attach_tooltip(self._roi_btn, "Enable or disable ROI drawing")

        btn_clear_roi = ttk.Button(outer, text="Clear ROI", command=self.clear_roi, width=12)
        btn_clear_roi.grid(row=1, column=3, padx=3, pady=2)
        self._attach_tooltip(btn_clear_roi, "Remove ROI and use full image")

        btn_load_roi_mask = ttk.Button(outer, text="Load ROI Mask", command=self.load_roi_mask, width=14)
        btn_load_roi_mask.grid(row=1, column=4, padx=3, pady=2)
        self._attach_tooltip(btn_load_roi_mask, "Load a binary mask file and use it as ROI")

        btn_save_roi_mask = ttk.Button(outer, text="Save ROI Mask", command=self.save_roi_mask, width=14)
        btn_save_roi_mask.grid(row=1, column=5, padx=3, pady=2)
        self._attach_tooltip(btn_save_roi_mask, "Save the current ROI as a binary mask image")

        btn_run = ttk.Button(outer, text="Run Costes + Analyze", command=self.run_analysis, width=18)
        btn_run.configure(style="Analyze.TButton")
        btn_run.grid(row=1, column=6, padx=3, pady=2)
        self._attach_tooltip(btn_run, "Compute thresholds and colocalization metrics")

        # Row 2
        cb_despike = ttk.Checkbutton(
            outer,
            text="Preprocess hot pixels",
            variable=self._despike_var,
        )
        cb_despike.grid(row=2, column=0, columnspan=2, padx=3, pady=2, sticky=tk.W)
        self._attach_tooltip(cb_despike, "Replace isolated bright outliers with local median")

        cb_autostretch = ttk.Checkbutton(
            outer,
            text="Auto-stretch display",
            variable=self._autostretch_var,
            command=self._on_autostretch_toggle,
        )
        cb_autostretch.grid(row=2, column=2, padx=3, pady=2, sticky=tk.W)
        self._attach_tooltip(cb_autostretch, "Enhance display contrast using robust percentiles")

        cb_histeq = ttk.Checkbutton(
            outer,
            text="Histogram equalization display",
            variable=self._histeq_var,
            command=self._on_histeq_toggle,
        )
        cb_histeq.grid(row=2, column=3, padx=3, pady=2, sticky=tk.W)
        self._attach_tooltip(
            cb_histeq,
            "Aggressive display contrast enhancement using per-channel histogram equalization",
        )

        cb_clahe = ttk.Checkbutton(
            outer,
            text="CLAHE display",
            variable=self._clahe_var,
            command=self._on_clahe_toggle,
        )
        cb_clahe.grid(row=2, column=4, padx=3, pady=2, sticky=tk.W)
        self._attach_tooltip(
            cb_clahe,
            "Adaptive local histogram equalization (CLAHE) for uneven illumination",
        )

        # Row 3
        cb_orthreg = ttk.Checkbutton(
            outer,
            text="Orthogonal regression (Costes)",
            variable=self._orthreg_var,
        )
        cb_orthreg.grid(row=3, column=0, columnspan=2, padx=3, pady=2, sticky=tk.W)
        self._attach_tooltip(
            cb_orthreg,
            "Use orthogonal (total least squares) regression instead of OLS "
            "when fitting the Ch2 vs Ch1 line for Costes thresholding. "
            "Minimises perpendicular distances rather than vertical residuals.",
        )

        debug_log_frame = ttk.Frame(outer)
        debug_log_frame.grid(row=3, column=2, padx=3, pady=2, sticky=tk.W)

        cb_costes_debug_log = ttk.Checkbutton(
            debug_log_frame,
            text="Debug Costes log",
            variable=self._costes_debug_log_var,
        )
        cb_costes_debug_log.pack(side=tk.LEFT)
        self._attach_tooltip(
            cb_costes_debug_log,
            "When enabled, write timestamped Costes regression/step/stop traces "
            "to costes_threshold_debug.log for debugging only.",
        )

        btn_clear_costes_debug_log = ttk.Button(
            debug_log_frame,
            text="Clear",
            width=7,
            command=self._clear_costes_debug_log,
        )
        btn_clear_costes_debug_log.pack(side=tk.LEFT, padx=(6, 0))
        self._attach_tooltip(
            btn_clear_costes_debug_log,
            "Delete costes_threshold_debug.log.",
        )

        btn_pdf = ttk.Button(outer, text="Save PDF Summary", command=self.save_pdf_summary, width=18)
        btn_pdf.grid(row=3, column=3, padx=3, pady=2)
        self._attach_tooltip(btn_pdf, "Export current plots and metrics to PDF")

        btn_quit = ttk.Button(outer, text="Quit", command=self._on_close, width=10)
        btn_quit.grid(row=3, column=4, padx=3, pady=2)
        self._attach_tooltip(btn_quit, "Close the application")

        ttk.Separator(outer, orient=tk.HORIZONTAL).grid(row=4, column=0, columnspan=5, sticky=tk.EW, pady=(4, 3))
        ttk.Label(outer, textvariable=self._status, foreground=GRAY).grid(
            row=5, column=0, columnspan=5, sticky=tk.W, padx=2, pady=(0, 2)
        )

        for col in range(5):
            outer.columnconfigure(col, weight=1)

        win.update_idletasks()
        screen_w = int(win.winfo_screenwidth())
        screen_h = int(win.winfo_screenheight())
        w = min(max(760, int(win.winfo_reqwidth()) + 8), max(620, screen_w - 80))
        h = min(int(win.winfo_reqheight()) + 8, max(280, screen_h - 120))
        win.geometry(f"{w}x{h}")

    def _attach_tooltip(self, widget, text: str):
        """Attach a short hover tooltip to a Tk/ttk widget."""
        widget.bind("<Enter>", lambda e: self._schedule_tooltip(widget, text), add="+")
        widget.bind("<Leave>", lambda e: self._hide_tooltip(), add="+")
        widget.bind("<ButtonPress>", lambda e: self._hide_tooltip(), add="+")

    def _schedule_tooltip(self, widget, text: str):
        if self._tooltip_after_id is not None:
            try:
                self.root.after_cancel(self._tooltip_after_id)
            except tk.TclError:
                pass
        self._tooltip_after_id = self.root.after(
            350, lambda: self._show_tooltip(widget, text)
        )

    def _show_tooltip(self, widget, text: str):
        self._tooltip_after_id = None
        if not widget.winfo_exists():
            return

        self._hide_tooltip()
        tw = tk.Toplevel(self.root)
        self._tooltip_win = tw
        tw.wm_overrideredirect(True)
        tw.attributes("-topmost", True)

        label = ttk.Label(tw, text=text, background="#111122", foreground=FG, padding=(6, 3))
        label.pack()

        tw.update_idletasks()
        x = widget.winfo_rootx() + 10
        y = widget.winfo_rooty() + widget.winfo_height() + 6
        sw = self.root.winfo_screenwidth()
        sh = self.root.winfo_screenheight()
        ww = tw.winfo_width()
        wh = tw.winfo_height()
        x = min(max(0, x), max(0, sw - ww))
        y = min(max(0, y), max(0, sh - wh))
        tw.geometry(f"+{x}+{y}")

    def _hide_tooltip(self):
        if self._tooltip_after_id is not None:
            try:
                self.root.after_cancel(self._tooltip_after_id)
            except tk.TclError:
                pass
            self._tooltip_after_id = None
        if self._tooltip_win is not None and self._tooltip_win.winfo_exists():
            self._tooltip_win.destroy()
        self._tooltip_win = None

    def _on_close(self):
        """Close app safely by canceling scheduled callbacks first."""
        # Guard against re-entrant close events fired by multiple windows.
        if getattr(self, "_is_closing", False):
            return
        self._is_closing = True

        if self._tooltip_after_id is not None:
            try:
                self.root.after_cancel(self._tooltip_after_id)
            except tk.TclError:
                pass
            self._tooltip_after_id = None

        if self._roi_overlay_after_id is not None:
            try:
                self.root.after_cancel(self._roi_overlay_after_id)
            except tk.TclError:
                pass
            self._roi_overlay_after_id = None

        self._hide_tooltip()

        try:
            if self._controls_win is not None and self._controls_win.winfo_exists():
                self._controls_win.destroy()
        except tk.TclError:
            pass

        try:
            if self._images_win is not None and self._images_win.winfo_exists():
                self._images_win.destroy()
        except tk.TclError:
            pass

        # Ensure no Matplotlib/Tk resources keep the process alive.
        try:
            plt.close("all")
        except Exception:
            pass

        try:
            self.root.quit()
        except tk.TclError:
            pass

        try:
            self.root.destroy()
        except tk.TclError:
            pass

    def _apply_theme(self):
        style = ttk.Style()
        try:
            style.theme_use("clam")
        except tk.TclError:
            pass
        style.configure(".",       background=BG2, foreground=FG)
        style.configure("TFrame",  background=BG2)
        style.configure("TLabel",  background=BG2, foreground=FG)
        style.configure("TButton", background=BG,  foreground=FG, relief=tk.FLAT, padding=4)
        style.map("TButton", background=[("active", "#313244")])
        style.configure("TLabelframe",       background=BG2, foreground=FG)
        style.configure("TLabelframe.Label", background=BG2, foreground=CYAN)
        style.configure("TNotebook",         background=BG2)
        style.configure("TNotebook.Tab",     background=BG,  foreground=FG, padding=[8, 3])
        style.map("TNotebook.Tab",           background=[("selected", BG2)], foreground=[("selected", YELLOW)])
        style.configure("TScale",            background=BG2)
        style.configure("TSeparator",        background=GRAY)

        # Emphasized action button for starting analysis.
        style.configure(
            "Analyze.TButton",
            background="#2f7d4a",
            foreground="#f3fff3",
            padding=5,
            font=("TkDefaultFont", 9, "bold"),
        )
        style.map(
            "Analyze.TButton",
            background=[("active", "#3f965c")],
            foreground=[("active", "#ffffff")],
        )

    def _build_image_panel(self, parent):
        fig, axes = plt.subplots(1, 2, figsize=(7, 4), facecolor=BG2)
        self._img_fig   = fig
        self._img_axes  = axes
        for ax in axes:
            ax.set_facecolor("#111122")
            ax.tick_params(colors=FG)
        axes[0].set_title("Channel 1", color=CYAN,    fontsize=10)
        axes[1].set_title("Channel 2", color=MAGENTA, fontsize=10)
        fig.tight_layout(pad=1.2)

        canvas = FigureCanvasTkAgg(fig, parent)
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self._img_canvas = canvas

        # ROI selector on the Channel-1 axis
        rs_kwargs = dict(
            useblit=True, button=[1],
            minspanx=5, minspany=5,
            spancoords="data", interactive=True,
        )
        try:
            rs_kwargs["props"] = dict(edgecolor=YELLOW, facecolor=YELLOW, alpha=0.15)
            self._roi_sel = RectangleSelector(axes[0], self._on_roi_select, **rs_kwargs)
        except TypeError:
            # older matplotlib – no 'props' kwarg
            rs_kwargs.pop("props", None)
            self._roi_sel = RectangleSelector(axes[0], self._on_roi_select, **rs_kwargs)

        # Lasso selector for freehand ROI drawing.
        lasso_kwargs = dict(useblit=True, button=[1])
        try:
            lasso_kwargs["props"] = dict(color=YELLOW, alpha=0.9, linewidth=1.5)
            self._lasso_sel = LassoSelector(axes[0], self._on_lasso_select, **lasso_kwargs)
        except TypeError:
            lasso_kwargs.pop("props", None)
            lasso_kwargs["lineprops"] = dict(color=YELLOW, alpha=0.9, linewidth=1.5)
            self._lasso_sel = LassoSelector(axes[0], self._on_lasso_select, **lasso_kwargs)

        self._set_roi_selector_active()

    def _set_roi_selector_active(self):
        """Activate only the selector matching the currently selected ROI mode."""
        mode = self._roi_mode_var.get()
        rect_on = self._roi_active and mode == "Rectangle"
        lasso_on = self._roi_active and mode == "Lasso"
        self._roi_sel.set_active(rect_on)
        self._lasso_sel.set_active(lasso_on)

    def _on_roi_mode_change(self, _event=None):
        self._set_roi_selector_active()
        if self._roi_active:
            self._status.set(
                f"ROI mode set to {self._roi_mode_var.get()}. Draw on Channel-1 image."
            )

    def _build_right_panel(self, parent):
        nb = ttk.Notebook(parent)
        nb.pack(fill=tk.BOTH, expand=True)

        # Tab 1 – Intensity Histograms
        t_hist = ttk.Frame(nb)
        nb.add(t_hist, text="Intensity Histograms")
        hist_fig, hist_axes = plt.subplots(2, 1, figsize=(6, 5),
                                           facecolor=BG, tight_layout=True)
        self._hist_fig   = hist_fig
        self._hist_axes  = hist_axes
        self._hist_t1_line = None
        self._hist_t2_line = None
        self._hist_canvas = FigureCanvasTkAgg(hist_fig, t_hist)
        self._hist_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Tab 2 – Costes thresholding
        t_costes = ttk.Frame(nb)
        nb.add(t_costes, text="Costes Thresholding")
        costes_fig, costes_axes = plt.subplots(1, 2, figsize=(8, 4),
                                               facecolor=BG, tight_layout=True)
        self._costes_fig    = costes_fig
        self._costes_axes   = costes_axes
        self._costes_canvas = FigureCanvasTkAgg(costes_fig, t_costes)
        self._costes_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Tab 3 – Results
        t_res = ttk.Frame(nb)
        nb.add(t_res, text="Results")
        self._build_results_tab(t_res)

        # Manual threshold controls (entry-only for responsiveness)
        sf = ttk.LabelFrame(parent, text="Manual Threshold Input")
        sf.pack(fill=tk.X, padx=4, pady=4)
        self._build_sliders(sf)

    def _build_results_tab(self, parent):
        outer = ttk.Frame(parent, padding=14)
        outer.pack(fill=tk.BOTH, expand=True)

        ttk.Label(
            outer,
            text=f"App version: {APP_VERSION}",
            foreground=CYAN,
            font=("TkDefaultFont", 9, "bold"),
        ).pack(anchor=tk.E, pady=(0, 6))

        # formula reference
        ref = (
            "M1  = Σ ch1ᵢ [ch2ᵢ > T₂] / Σ ch1ᵢ\n"
            "M2  = Σ ch2ᵢ [ch1ᵢ > T₁] / Σ ch2ᵢ\n"
            "tM1 = Σ ch1ᵢ [ch1ᵢ > T₁ ∧ ch2ᵢ > T₂] / Σ ch1ᵢ\n"
            "tM2 = Σ ch2ᵢ [ch1ᵢ > T₁ ∧ ch2ᵢ > T₂] / Σ ch2ᵢ"
        )
        ttk.Label(outer, text=ref, foreground=GRAY,
                  font=("Courier New", 9), justify=tk.LEFT).pack(anchor=tk.W, pady=(0, 8))
        ttk.Separator(outer, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=4)

        grid = ttk.Frame(outer)
        grid.pack(fill=tk.BOTH, expand=True)

        rows = [
            ("Costes Threshold  T₁ (Ch1)",              "costes_t1"),
            ("Costes Threshold  T₂ (Ch2)",              "costes_t2"),
            None,
            ("Mander's M1   — Ch1 fraction colocalizing with Ch2", "m1"),
            ("Mander's M2   — Ch2 fraction colocalizing with Ch1", "m2"),
            None,
            ("Thresholded  tM1  (both channels above Costes T)",   "tm1"),
            ("Thresholded  tM2  (both channels above Costes T)",   "tm2"),
            None,
            ("Pearson's r  (pixels above both thresholds)",        "pearson"),
            ("Overlap coefficient R  (Manders 1993)",              "overlap"),
            ("k1   = Σ(ch1·ch2) / Σ(ch1²)",                       "k1"),
            ("k2   = Σ(ch1·ch2) / Σ(ch2²)",                       "k2"),
            None,
            ("Costes randomization r  (all ROI pixels)",          "costes_r_obs"),
            ("Randomized r mean ± SD",                             "costes_r_rand"),
            ("Costes randomization significance",                  "costes_sig"),
            ("Monte Carlo p  (random r ≥ observed r)",            "costes_p"),
        ]

        self._rv = {}
        row_idx = 0
        for item in rows:
            if item is None:
                ttk.Separator(grid, orient=tk.HORIZONTAL).grid(
                    row=row_idx, column=0, columnspan=2, sticky=tk.EW, pady=5)
                row_idx += 1
                continue
            label, key = item
            ttk.Label(grid, text=label + ":", anchor=tk.W).grid(
                row=row_idx, column=0, sticky=tk.W, pady=2)
            var = tk.StringVar(value="—")
            self._rv[key] = var
            ttk.Label(grid, textvariable=var,
                      font=("Courier New", 11, "bold"),
                      foreground=GREEN).grid(row=row_idx, column=1, sticky=tk.W, padx=12)
            row_idx += 1

        grid.columnconfigure(0, weight=1)

    def _build_sliders(self, parent):
        self._t1_max = 255.0
        self._t2_max = 255.0
        self._psf_var = tk.StringVar(value="3.0")
        self._bg1_var = tk.StringVar(value="0.0")
        self._bg2_var = tk.StringVar(value="0.0")
        self._t1_entry_var = tk.StringVar(value="0.0")
        self._t2_entry_var = tk.StringVar(value="0.0")

        ttk.Label(parent, text="T₁ (Ch1):").grid(row=0, column=0, padx=6, sticky=tk.W, pady=3)
        self._entry_t1 = ttk.Entry(parent, textvariable=self._t1_entry_var, width=9)
        self._entry_t1.grid(row=0, column=1, padx=4, sticky=tk.W)
        self._entry_t1.bind("<Return>",    lambda e: self._entry_update(1))
        self._entry_t1.bind("<FocusOut>",  lambda e: self._entry_update(1))

        ttk.Label(parent, text="T₂ (Ch2):").grid(row=1, column=0, padx=6, sticky=tk.W, pady=3)
        self._entry_t2 = ttk.Entry(parent, textvariable=self._t2_entry_var, width=9)
        self._entry_t2.grid(row=1, column=1, padx=4, sticky=tk.W)
        self._entry_t2.bind("<Return>",   lambda e: self._entry_update(2))
        self._entry_t2.bind("<FocusOut>", lambda e: self._entry_update(2))

        ttk.Button(parent, text="Apply", command=lambda: self._slider_update(draw_plots=True)).grid(
            row=0, column=2, rowspan=2, padx=8, sticky=tk.NS)

        ttk.Label(parent, text="PSF (px, Costes randomization block):").grid(
            row=2, column=0, padx=6, sticky=tk.W, pady=3)
        self._entry_psf = ttk.Entry(parent, textvariable=self._psf_var, width=9)
        self._entry_psf.grid(row=2, column=1, padx=4, sticky=tk.W)
        self._entry_psf.bind("<Return>", lambda e: self._on_psf_change())
        self._entry_psf.bind("<FocusOut>", lambda e: self._on_psf_change())

        ttk.Label(parent, text="Costes rise stop (% over min r):").grid(
            row=3, column=0, padx=6, sticky=tk.W, pady=3)
        self._entry_costes_stop_pct = ttk.Entry(parent, textvariable=self._costes_stop_pct_var, width=9)
        self._entry_costes_stop_pct.grid(row=3, column=1, padx=4, sticky=tk.W)
        self._entry_costes_stop_pct.bind("<Return>", lambda e: self._on_costes_stop_pct_change())
        self._entry_costes_stop_pct.bind("<FocusOut>", lambda e: self._on_costes_stop_pct_change())

        ttk.Label(parent, text="Background Ch1:").grid(
            row=4, column=0, padx=6, sticky=tk.W, pady=3)
        self._entry_bg1 = ttk.Entry(parent, textvariable=self._bg1_var, width=9)
        self._entry_bg1.grid(row=4, column=1, padx=4, sticky=tk.W)
        self._entry_bg1.bind("<Return>", lambda e: self._on_background_change())
        self._entry_bg1.bind("<FocusOut>", lambda e: self._on_background_change())

        ttk.Label(parent, text="Background Ch2:").grid(
            row=5, column=0, padx=6, sticky=tk.W, pady=3)
        self._entry_bg2 = ttk.Entry(parent, textvariable=self._bg2_var, width=9)
        self._entry_bg2.grid(row=5, column=1, padx=4, sticky=tk.W)
        self._entry_bg2.bind("<Return>", lambda e: self._on_background_change())
        self._entry_bg2.bind("<FocusOut>", lambda e: self._on_background_change())

        ttk.Button(parent, text="Apply Stop", command=self._on_costes_stop_pct_change).grid(
            row=3, column=2, padx=8, sticky=tk.NS)
        ttk.Button(parent, text="Apply BG", command=self._on_background_change).grid(
            row=4, column=2, rowspan=2, padx=8, sticky=tk.NS)

        ttk.Label(parent, text="(press Enter or Tab to apply typed value)",
                  foreground=GRAY, font=("TkDefaultFont", 8)).grid(
            row=6, column=0, columnspan=3, sticky=tk.W, padx=6, pady=(0, 2))

        parent.columnconfigure(1, weight=1)

    def _get_psf_block_size(self) -> int:
        """Get randomization block size in pixels from the PSF field."""
        try:
            psf_px = float(self._psf_var.get())
        except ValueError:
            psf_px = 3.0
        if not np.isfinite(psf_px):
            psf_px = 3.0
        psf_px = max(2.0, psf_px)
        self._psf_var.set(f"{psf_px:.1f}")
        return int(round(psf_px))

    def _get_costes_rise_stop_fraction(self) -> float:
        """Get the Costes rise-stop threshold as a fractional value."""
        try:
            pct = float(self._costes_stop_pct_var.get())
        except ValueError:
            pct = 1.0
        if not np.isfinite(pct):
            pct = 1.0
        pct = max(0.0, pct)
        self._costes_stop_pct_var.set(f"{pct:.2f}")
        return pct / 100.0

    def _on_psf_change(self):
        """Apply PSF edits and re-run analysis when images are loaded."""
        self._get_psf_block_size()  # sanitize and normalize field formatting

        if self.ch1 is None or self.ch2 is None:
            return

        self._status.set("PSF updated - rerunning analysis...")
        self.root.update_idletasks()
        self.run_analysis()

    def _on_costes_stop_pct_change(self):
        """Apply Costes rise-stop edits and re-run analysis when images are loaded."""
        stop_pct = 100.0 * self._get_costes_rise_stop_fraction()

        if self.ch1 is None or self.ch2 is None:
            self._status.set(f"Costes rise stop set to {stop_pct:.2f}% over min r.")
            return

        self._status.set(
            f"Costes rise stop set to {stop_pct:.2f}% over min r - rerunning analysis..."
        )
        self.root.update_idletasks()
        self.run_analysis()

    def _get_background_levels(self) -> tuple[float, float]:
        """Get per-channel constant background values for subtraction."""
        try:
            bg1 = float(self._bg1_var.get())
        except ValueError:
            bg1 = 0.0
        try:
            bg2 = float(self._bg2_var.get())
        except ValueError:
            bg2 = 0.0

        if not np.isfinite(bg1):
            bg1 = 0.0
        if not np.isfinite(bg2):
            bg2 = 0.0

        bg1 = max(0.0, bg1)
        bg2 = max(0.0, bg2)
        self._bg1_var.set(f"{bg1:.1f}")
        self._bg2_var.set(f"{bg2:.1f}")
        return bg1, bg2

    def _apply_background_subtraction(self, c1: np.ndarray, c2: np.ndarray):
        """Subtract constant per-channel background and clamp to non-negative."""
        bg1, bg2 = self._get_background_levels()
        if bg1 > 0.0:
            c1 = np.clip(c1 - bg1, 0.0, None)
        if bg2 > 0.0:
            c2 = np.clip(c2 - bg2, 0.0, None)
        return c1, c2

    def _on_background_change(self):
        """Apply background subtraction edits and refresh analysis if possible."""
        bg1, bg2 = self._get_background_levels()
        self._sync_slider_range()

        if self.ch1 is None or self.ch2 is None:
            self._status.set(
                f"Background subtraction set (Ch1={bg1:.1f}, Ch2={bg2:.1f})."
            )
            return

        self._status.set(
            f"Background subtraction set (Ch1={bg1:.1f}, Ch2={bg2:.1f}) - rerunning analysis..."
        )
        self.root.update_idletasks()
        self.run_analysis()

    def _entry_update(self, channel: int):
        """Parse the manually typed threshold value and apply it."""
        var    = self._t1_entry_var if channel == 1 else self._t2_entry_var
        try:
            val = float(var.get())
        except ValueError:
            var.set("0.0")
            return
        lo = 0.0
        hi = self._t1_max if channel == 1 else self._t2_max
        val = max(lo, min(hi, val))
        var.set(f"{val:.1f}")
        self._slider_update(draw_plots=True)

    # ── image loading ──────────────────────────────────────────────────────────

    def _load_file(self, title: str):
        path = filedialog.askopenfilename(
            title=title,
            filetypes=[
                ("Images", "*.tif *.tiff *.png *.jpg *.jpeg *.bmp"),
                ("TIFF",   "*.tif *.tiff"),
                ("All",    "*.*"),
            ],
        )
        if not path:
            return None

        return self._read_image_file(path)

    def _read_image_file(self, path: str):
        path = str(path)

        # ── try tifffile first (best for microscopy TIFFs) ──────────────────
        if _HAS_TIFFFILE and path.lower().endswith((".tif", ".tiff")):
            try:
                img = tifffile.imread(path)
                img = self._squeeze_to_2d(img)
                return img.astype(np.float32)
            except Exception:
                pass

        # ── fallback to Pillow ───────────────────────────────────────────────
        if _HAS_PIL:
            try:
                pil = Image.open(path)
                # keep 16-bit precision where possible
                if pil.mode in ("I;16", "I;16B"):
                    arr = np.frombuffer(pil.tobytes(), dtype=np.uint16
                                        ).reshape(pil.size[1], pil.size[0])
                elif pil.mode == "I":
                    arr = np.array(pil, dtype=np.int32)
                else:
                    arr = np.array(pil.convert("L"), dtype=np.uint8)
                return arr.astype(np.float32)
            except Exception as exc:
                messagebox.showerror("Load error", str(exc))
                return None

        messagebox.showerror(
            "Missing library",
            "Install 'tifffile' and/or 'Pillow':\n  pip install tifffile Pillow",
        )
        return None

    @staticmethod
    def _squeeze_to_2d(img: np.ndarray) -> np.ndarray:
        """Reduce tifffile output to a 2-D (Y, X) grayscale array."""
        if img.ndim == 2:
            return img
        if img.ndim == 3:
            # XYC (H, W, C)
            if img.shape[2] in (3, 4):
                return img[..., :3].mean(axis=2)
            # CXY (C, H, W)
            if img.shape[0] in (3, 4):
                return img[:3].mean(axis=0)
            # single extra dim – take first slice (Z or T)
            return img[0]
        # ≥ 4-D: take first slice along all leading dims
        while img.ndim > 2:
            img = img[0]
        return img

    def load_ch1(self):
        img = self._load_file("Open Channel 1")
        if img is not None:
            self._set_channels(ch1=img, ch2=self.ch2)
            self._status.set(
                f"Ch1 loaded — shape {img.shape}, "
                f"range [{img.min():.0f}, {img.max():.0f}]"
            )

    def load_ch2(self):
        img = self._load_file("Open Channel 2")
        if img is not None:
            self._set_channels(ch1=self.ch1, ch2=img)
            self._status.set(
                f"Ch2 loaded — shape {img.shape}, "
                f"range [{img.min():.0f}, {img.max():.0f}]"
            )

    def load_demo_set(self):
        preset = self._demo_preset_var.get()
        ch1, ch2 = self._generate_demo_images(preset=preset)
        self._set_channels(ch1=ch1, ch2=ch2)
        self._status.set(
            f"Demo set loaded — {preset}: {self._demo_presets[preset]}"
        )

    def export_demo_set(self):
        if not _HAS_TIFFFILE:
            messagebox.showerror(
                "Missing dependency",
                "TIFF export requires tifffile. Install with: pip install tifffile",
            )
            return

        preset = self._demo_preset_var.get()
        ch1, ch2 = self._generate_demo_images(preset=preset)
        preset_slug = preset.lower().replace(" ", "_")

        out_dir = Path(__file__).resolve().parent / "demo_exports"
        out_dir.mkdir(parents=True, exist_ok=True)

        ch1_path = out_dir / f"demo_{preset_slug}_ch1.tif"
        ch2_path = out_dir / f"demo_{preset_slug}_ch2.tif"

        # Export as uint16 for microscopy-style tooling compatibility.
        tifffile.imwrite(ch1_path, np.clip(ch1, 0, 65535).astype(np.uint16))
        tifffile.imwrite(ch2_path, np.clip(ch2, 0, 65535).astype(np.uint16))

        self._status.set(
            f"Demo TIFFs exported: {ch1_path.name}, {ch2_path.name}"
        )

    def _set_channels(self, ch1: np.ndarray | None, ch2: np.ndarray | None):
        self.ch1 = ch1
        self.ch2 = ch2
        self.roi = None
        self.roi_polygon = None
        self.roi_mask = None
        self._roi_kind = None
        if hasattr(self, "_roi_sel"):
            self._roi_sel.clear()
        self._refresh_images()
        self._sync_slider_range()

    def _mask_to_roi_bounds(self, mask: np.ndarray):
        ys, xs = np.where(mask)
        return (int(xs.min()), int(ys.min()), int(xs.max()) + 1, int(ys.max()) + 1)

    def _apply_mask_roi(self, mask: np.ndarray, source: str):
        if self.ch1 is None:
            messagebox.showwarning("No image", "Load Channel 1 before using an ROI mask.")
            return False

        if mask.shape != self.ch1.shape[:2]:
            messagebox.showerror(
                "Mask size mismatch",
                "The ROI mask must match the loaded image size.\n"
                f"Mask: {mask.shape}   Ch1: {self.ch1.shape[:2]}",
            )
            return False

        n_sel = int(mask.sum())
        if n_sel == 0:
            messagebox.showerror("Empty mask", "The selected ROI mask does not contain any non-zero pixels.")
            return False

        self.roi = self._mask_to_roi_bounds(mask)
        self.roi_polygon = None
        self.roi_mask = mask
        self._roi_kind = "mask"
        self._roi_active = False
        self._set_roi_selector_active()
        self._roi_sel.clear()
        self._roi_btn.config(text="Draw ROI")
        self._refresh_images()
        self._status.set(f"ROI mask loaded from {source} ({n_sel} pixels).")
        return True

    def _current_roi_mask_image(self):
        if self.ch1 is None:
            return None

        h, w = self.ch1.shape[:2]
        mask = np.zeros((h, w), dtype=bool)

        if self._roi_kind in {"lasso", "mask"} and self.roi_mask is not None:
            if self.roi_mask.shape != (h, w):
                return None
            return self.roi_mask.astype(bool, copy=True)

        if self._roi_kind == "rect" and self.roi is not None:
            x1, y1, x2, y2 = self.roi
            mask[y1:y2, x1:x2] = True
            return mask

        return None

    def load_roi_mask(self):
        if self.ch1 is None:
            messagebox.showwarning("No image", "Load Channel 1 before loading an ROI mask.")
            return

        path = filedialog.askopenfilename(
            title="Open ROI Mask",
            filetypes=[
                ("Images", "*.tif *.tiff *.png *.jpg *.jpeg *.bmp"),
                ("TIFF", "*.tif *.tiff"),
                ("All", "*.*"),
            ],
        )
        if not path:
            return

        img = self._read_image_file(path)
        if img is None:
            return

        mask = np.asarray(img) > 0
        self._apply_mask_roi(mask, Path(path).name)

    def save_roi_mask(self):
        if self.ch1 is None:
            messagebox.showwarning("No image", "Load Channel 1 before saving an ROI mask.")
            return

        mask = self._current_roi_mask_image()
        if mask is None or not np.any(mask):
            messagebox.showwarning("No ROI", "Define or load an ROI before saving a mask.")
            return

        path = filedialog.asksaveasfilename(
            title="Save ROI Mask",
            defaultextension=".png",
            filetypes=[
                ("PNG", "*.png"),
                ("TIFF", "*.tif *.tiff"),
                ("BMP", "*.bmp"),
                ("All", "*.*"),
            ],
        )
        if not path:
            return

        out = (mask.astype(np.uint8) * 255)
        suffix = Path(path).suffix.lower()
        try:
            if suffix in {".tif", ".tiff"}:
                if not _HAS_TIFFFILE:
                    messagebox.showerror(
                        "Missing dependency",
                        "TIFF export requires tifffile. Install with: pip install tifffile",
                    )
                    return
                tifffile.imwrite(path, out)
            else:
                if not _HAS_PIL:
                    messagebox.showerror(
                        "Missing dependency",
                        "PNG/BMP export requires Pillow. Install with: pip install Pillow",
                    )
                    return
                Image.fromarray(out, mode="L").save(path)
        except Exception as exc:
            messagebox.showerror("Save error", str(exc))
            return

        self._status.set(f"ROI mask saved to {Path(path).name}.")

    @staticmethod
    def _generate_demo_images(shape: tuple[int, int] = (256, 256),
                              seed: int = 7,
                              preset: str = "Partial overlap") -> tuple[np.ndarray, np.ndarray]:
        """
        Create a deterministic synthetic two-channel dataset for testing.

        Design rationale
        ----------------
        The Costes algorithm finds a threshold by scanning T1 downward and
        monitoring Pearson r of *below*-threshold pixels until r ≤ 0.  For
        this to work, the background must be truly decorrelated between the
        two channels.  Wide Gaussian spots (σ > 6 px, amplitude >> background)
        have tails that span many pixels and maintain positive correlation at
        ANY threshold below the spot peak, preventing the r-vs-T curve from
        ever crossing zero.

        We therefore use:
        - σ = 4–6 px (tight), amplitude 800–1400 DN above flat 100 DN base.
          The tail decays to background noise floor within ~10 px of the spot
          centre, leaving ~95 % of all pixels as genuinely uncorrelated background.
        - Background: independent per-channel Gaussian noise (no shared spatial
          structure).  This forces r(background) ≈ 0 at Any threshold that
          excludes the spot footprints, giving a clean Costes crossing.
        """
        seed_map = {
            "High overlap":     seed,
            "Partial overlap":  seed + 11,
            "Offset structures": seed + 23,
            "Low overlap":      seed + 37,
            "Pure bleedthrough": seed + 47,
            "Randomization demo": seed + 53,
            "Randomization control": seed + 67,
        }
        rng = np.random.default_rng(seed_map.get(preset, seed + 11))
        h, w = shape
        yy, xx = np.mgrid[0:h, 0:w]

        def gauss(cx, cy, sigma, amp):
            return amp * np.exp(-((xx - cx) ** 2 + (yy - cy) ** 2) / (2.0 * sigma ** 2))

        # Flat, fully INDEPENDENT background – ensures r(background) ≈ 0
        ch1 = 100.0 + rng.normal(0.0, 28.0, shape)
        ch2 =  95.0 + rng.normal(0.0, 28.0, shape)

        if preset == "High overlap":
            # High-overlap should have low diffuse baseline and dominant shared
            # puncta, otherwise whole-image denominators suppress tM1/tM2.
            ch1 = 10.0 + rng.normal(0.0, 3.0, shape)
            ch2 = 10.0 + rng.normal(0.0, 3.0, shape)

            # 8 tightly co-localised spots; minimal channel-specific signal
            for cx, cy, s, a1, a2 in [
                ( 54,  64, 4.8, 1300, 1260),
                ( 96, 120, 5.6, 1580, 1520),
                (148,  84, 4.8, 1220, 1180),
                (192, 152, 6.0, 1760, 1700),
                (122, 206, 5.2, 1380, 1330),
                ( 62, 188, 4.8, 1160, 1110),
                (206,  72, 4.4, 1260, 1210),
                (168, 208, 4.8, 1240, 1190),
            ]:
                ch1 += gauss(cx,     cy, s,        a1)
                ch2 += gauss(cx + 1, cy, s * 1.05, a2)
            # Small channel-specific spots (~20-25% of shared amplitude)
            for cx, cy, s, a in [(36, 34, 3.8, 300), (224, 44, 4.0, 320)]:
                ch1 += gauss(cx, cy, s, a)
            for cx, cy, s, a in [(42, 220, 3.8, 290), (210, 108, 4.0, 315)]:
                ch2 += gauss(cx, cy, s, a)

        elif preset == "Offset structures":
            # 5 pairs displaced by 10–14 px (≈ 2–3 σ) – large offset, low true overlap
            for cx, cy, s, a1, a2, dx, dy in [
                ( 58,  68, 5.0, 1300, 1250,  12,   9),
                (104, 124, 6.0, 1450, 1380, -13,   8),
                (150,  90, 5.5, 1200, 1150,  10, -11),
                (192, 156, 6.5, 1550, 1480, -11,  -9),
                (118, 204, 5.0, 1280, 1200,   9,  12),
            ]:
                ch1 += gauss(cx,      cy,      s,        a1)
                ch2 += gauss(cx + dx, cy + dy, s * 1.05, a2)
            # Unique spots
            for cx, cy, s, a in [(42, 34, 4.5, 750), (226, 114, 5.0, 700)]:
                ch1 += gauss(cx, cy, s, a)
            for cx, cy, s, a in [(218, 40, 4.5, 720), (38, 150, 5.0, 680)]:
                ch2 += gauss(cx, cy, s, a)

        elif preset == "Low overlap":
            # Low-overlap should still yield a moderate Costes threshold, not
            # an extreme value driven by a handful of very bright isolated
            # puncta. Use lower-amplitude unique spots plus a weak broad shared
            # haze so the r-vs-threshold curve crosses zero at a sensible level.
            ch1 = 22.0 + rng.normal(0.0, 4.0, shape)
            ch2 = 22.0 + rng.normal(0.0, 4.0, shape)

            diffuse = rng.normal(0.0, 1.0, shape)
            ch1 += 5.5 * diffuse
            ch2 += 5.5 * diffuse

            # Mostly independent spot populations with matched intensity budgets.
            low_overlap_spots = [
                (4.4, 340),
                (4.8, 410),
                (4.6, 370),
                (5.0, 440),
                (4.4, 390),
            ]
            for (cx, cy), (s, a) in zip(
                [(40, 44), (92, 126), (156, 84), (214, 152), (122, 212)],
                low_overlap_spots,
            ):
                ch1 += gauss(cx, cy, s, a)
            for (cx, cy), (s, a) in zip(
                [(210, 42), (156, 188), (78, 210), (44, 146), (206, 220)],
                low_overlap_spots,
            ):
                ch2 += gauss(cx, cy, s, a)
            # Weak shared component
            ch1 += gauss(128, 128, 8.5, 70)
            ch2 += gauss(128, 128, 8.5, 70)

        elif preset == "Pure bleedthrough":
            # Simulate optical bleed-through: Ch2 is mostly a scaled copy of Ch1
            # plus detector/background noise, with no channel-specific structures.
            ch1 = 28.0 + rng.normal(0.0, 5.5, shape)

            for cx, cy, s, a in [
                (42, 56, 4.6, 900),
                (88, 118, 5.2, 980),
                (136, 84, 4.8, 860),
                (182, 142, 5.8, 1040),
                (220, 92, 4.4, 820),
                (164, 204, 5.0, 920),
                (74, 196, 4.6, 840),
            ]:
                ch1 += gauss(cx, cy, s, a)

            # Add broad diffuse structure to make the bleed-through pattern obvious.
            mesoscale = (
                0.9 * np.sin(xx / 15.0)
                + 0.7 * np.cos(yy / 18.0)
                + 0.5 * np.sin((xx + yy) / 24.0)
            )
            mesoscale = (mesoscale - mesoscale.min()) / (mesoscale.max() - mesoscale.min() + 1e-12)
            ch1 += 180.0 * mesoscale

            bleed_frac = 0.38
            ch2 = (12.0 + rng.normal(0.0, 5.0, shape)) + (bleed_frac * ch1)

        elif preset == "Randomization demo":
            # Designed to make Costes randomization behavior obvious:
            # strong shared spatial structure yields high observed r, while
            # block shuffling Ch2 suppresses r toward random baseline.
            ch1 = 18.0 + rng.normal(0.0, 4.5, shape)
            ch2 = 18.0 + rng.normal(0.0, 4.5, shape)

            mesoscale = (
                np.sin(xx / 12.0)
                + 0.8 * np.cos(yy / 17.0)
                + 0.6 * np.sin((xx + yy) / 21.0)
            )
            mesoscale = (mesoscale - mesoscale.min()) / (mesoscale.max() - mesoscale.min() + 1e-12)
            ch1 += 240.0 * mesoscale
            ch2 += 228.0 * mesoscale

            for cx, cy, s, a in [
                (40, 52, 5.0, 620),
                (86, 94, 5.8, 700),
                (132, 70, 4.6, 560),
                (184, 118, 5.5, 680),
                (220, 168, 5.2, 640),
                (160, 208, 5.0, 600),
            ]:
                ch1 += gauss(cx, cy, s, a)
                ch2 += gauss(cx + 1, cy - 1, s * 1.03, a * 0.96)

            # Add channel-specific content so the dataset is not trivially identical.
            for cx, cy, s, a in [(26, 210, 4.0, 300), (222, 44, 4.4, 320)]:
                ch1 += gauss(cx, cy, s, a)
            for cx, cy, s, a in [(44, 24, 4.2, 290), (206, 222, 4.2, 310)]:
                ch2 += gauss(cx, cy, s, a)

            shared_noise = rng.normal(0.0, 1.0, shape)
            ch1 += 18.0 * shared_noise
            ch2 += 17.0 * shared_noise

        elif preset == "Randomization control":
            # Negative control for randomization: mostly independent channel
            # structures so observed r should be close to the randomized null.
            ch1 = 36.0 + rng.normal(0.0, 5.5, shape)
            ch2 = 36.0 + rng.normal(0.0, 5.5, shape)

            # Independent puncta populations with matched intensity budgets.
            for cx, cy, s, a in [
                (34, 42, 4.2, 210),
                (84, 122, 4.8, 250),
                (142, 86, 4.4, 225),
                (208, 156, 5.2, 270),
                (126, 214, 4.6, 240),
                (226, 56, 4.2, 220),
            ]:
                ch1 += gauss(cx, cy, s, a)

            for cx, cy, s, a in [
                (218, 40, 4.2, 210),
                (160, 206, 4.8, 250),
                (70, 214, 4.4, 225),
                (44, 150, 5.2, 270),
                (206, 222, 4.6, 240),
                (120, 34, 4.2, 220),
            ]:
                ch2 += gauss(cx, cy, s, a)

            # Only a very weak shared component to avoid exact zero-correlation.
            weak_shared = rng.normal(0.0, 1.0, shape)
            ch1 += 2.0 * weak_shared
            ch2 += 2.0 * weak_shared

        else:
            # Partial overlap (default):
            # Use matched baseline statistics and balanced structure budgets so
            # M1 and M2 are moderately close (while still not identical).
            ch1 = 36.0 + rng.normal(0.0, 5.0, shape)
            ch2 = 36.0 + rng.normal(0.0, 5.0, shape)

            # 3 fully shared + 2 slightly offset pairs + 4 channel-specific each
            for cx, cy, s, a1, a2 in [
                ( 98,  78, 4.8, 900, 880),
                (150, 170, 5.2, 860, 840),
                ( 64, 194, 4.6, 820, 800),
            ]:
                ch1 += gauss(cx, cy, s,        a1)
                ch2 += gauss(cx, cy, s * 1.05, a2)
            # Slightly offset pairs (offset 6–7 px ≈ 1.2 σ – partial pixel-level overlap)
            for cx, cy, s, a1, a2, dx, dy in [
                (200, 120, 4.8, 760, 740,  7,  6),
                (110, 206, 5.0, 730, 710, -6,  7),
            ]:
                ch1 += gauss(cx,      cy,      s, a1)
                ch2 += gauss(cx + dx, cy + dy, s, a2)
            # Channel-specific
            for cx, cy, s, a in [
                ( 38,  34, 4.2, 840),
                (142,  44, 4.6, 790),
                (226, 208, 4.4, 850),
                (208,  24, 4.4, 740),
            ]:
                ch1 += gauss(cx, cy, s, a)
            for cx, cy, s, a in [
                (224,  42, 4.2, 830),
                ( 44, 136, 4.8, 800),
                (160, 226, 4.4, 860),
                (186,  30, 4.6, 760),
            ]:
                ch2 += gauss(cx, cy, s, a)

            # Equalize global intensity budget to avoid built-in channel bias.
            ch2 *= ch1.sum() / max(ch2.sum(), 1e-12)

        # A handful of hot pixels to test threshold robustness (not too bright)
        hot_lo, hot_hi = (200, 450)
        if preset == "Randomization control":
            # Keep this control at a moderate dynamic range to avoid
            # high-threshold picks driven by very sparse bright outliers.
            hot_lo, hot_hi = (40, 120)

        hy = rng.integers(0, h, size=12)
        hx = rng.integers(0, w, size=12)
        ch1[hy, hx] += rng.uniform(hot_lo, hot_hi, size=12)
        hy = rng.integers(0, h, size=12)
        hx = rng.integers(0, w, size=12)
        ch2[hy, hx] += rng.uniform(hot_lo, hot_hi, size=12)

        ch1 = np.clip(ch1, 0, 4095).astype(np.float32)
        ch2 = np.clip(ch2, 0, 4095).astype(np.float32)
        return ch1, ch2

    def _sync_slider_range(self):
        bg1, bg2 = self._get_background_levels()
        if self.ch1 is not None:
            self._t1_max = max(0.0, float(self.ch1.max()) - bg1)
        if self.ch2 is not None:
            self._t2_max = max(0.0, float(self.ch2.max()) - bg2)

    @staticmethod
    def _autostretch_for_display(img: np.ndarray) -> np.ndarray:
        """Create a contrast-stretched copy for display only."""
        if img.size == 0:
            return img
        lo, hi = np.percentile(img, [2.0, 98.0])
        if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
            return img
        stretched = np.clip((img - lo) / (hi - lo), 0.0, 1.0)
        return stretched.astype(np.float32)

    @staticmethod
    def _histeq_for_display(img: np.ndarray) -> np.ndarray:
        """Create a histogram-equalized copy for display only."""
        if img.size == 0:
            return img

        flat = img.ravel().astype(np.float64)
        lo, hi = np.percentile(flat, [0.5, 99.5])
        if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
            return img

        clipped = np.clip(flat, lo, hi)
        hist, bin_edges = np.histogram(clipped, bins=512, range=(lo, hi))
        cdf = hist.cumsum().astype(np.float64)
        if cdf[-1] <= 0:
            return img
        cdf /= cdf[-1]

        vals = np.interp(clipped, bin_edges[:-1], cdf)
        return vals.reshape(img.shape).astype(np.float32)

    @staticmethod
    def _clahe_for_display(img: np.ndarray,
                           tile_size: int = 32,
                           clip_limit: float = 0.01,
                           bins: int = 256) -> np.ndarray:
        """Create a CLAHE-enhanced copy for display only."""
        if img.size == 0:
            return img

        arr = img.astype(np.float64)
        lo, hi = np.percentile(arr, [0.5, 99.5])
        if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
            return img

        arr = np.clip((arr - lo) / (hi - lo), 0.0, 1.0)
        h, w = arr.shape
        ny = int(np.ceil(h / tile_size))
        nx = int(np.ceil(w / tile_size))
        lut_map = np.zeros((ny, nx, bins), dtype=np.float64)

        for ty in range(ny):
            y0 = ty * tile_size
            y1 = min((ty + 1) * tile_size, h)
            for tx in range(nx):
                x0 = tx * tile_size
                x1 = min((tx + 1) * tile_size, w)
                tile = arr[y0:y1, x0:x1]
                hist, _ = np.histogram(tile, bins=bins, range=(0.0, 1.0))

                max_count = max(int(clip_limit * tile.size), 1)
                excess = np.maximum(hist - max_count, 0)
                clipped = hist - excess
                redistribute = int(excess.sum())
                if redistribute > 0:
                    q, r = divmod(redistribute, bins)
                    clipped += q
                    if r > 0:
                        clipped[:r] += 1

                cdf = clipped.cumsum().astype(np.float64)
                if cdf[-1] > 0:
                    cdf /= cdf[-1]
                lut_map[ty, tx] = cdf

        out = np.empty_like(arr, dtype=np.float64)
        ys = np.arange(h)
        xs = np.arange(w)
        for y in ys:
            gy = y / tile_size - 0.5
            y0 = int(np.floor(gy))
            y1 = y0 + 1
            wy = gy - y0
            y0 = min(max(y0, 0), ny - 1)
            y1 = min(max(y1, 0), ny - 1)

            row = arr[y]
            bidx = np.minimum((row * (bins - 1)).astype(np.int32), bins - 1)
            for x in xs:
                gx = x / tile_size - 0.5
                x0 = int(np.floor(gx))
                x1 = x0 + 1
                wx = gx - x0
                x0 = min(max(x0, 0), nx - 1)
                x1 = min(max(x1, 0), nx - 1)

                v00 = lut_map[y0, x0, bidx[x]]
                v10 = lut_map[y0, x1, bidx[x]]
                v01 = lut_map[y1, x0, bidx[x]]
                v11 = lut_map[y1, x1, bidx[x]]

                top = (1.0 - wx) * v00 + wx * v10
                bot = (1.0 - wx) * v01 + wx * v11
                out[y, x] = (1.0 - wy) * top + wy * bot

        return out.astype(np.float32)

    def _on_autostretch_toggle(self):
        if self._autostretch_var.get():
            self._histeq_var.set(False)
            self._clahe_var.set(False)
        self._refresh_images()
        mode = "ON" if self._autostretch_var.get() else "OFF"
        self._status.set(f"Auto-stretch display {mode}.")

    def _on_histeq_toggle(self):
        if self._histeq_var.get():
            self._autostretch_var.set(False)
            self._clahe_var.set(False)
        self._refresh_images()
        mode = "ON" if self._histeq_var.get() else "OFF"
        self._status.set(f"Histogram equalization display {mode}.")

    def _on_clahe_toggle(self):
        if self._clahe_var.get():
            self._autostretch_var.set(False)
            self._histeq_var.set(False)
        self._refresh_images()
        mode = "ON" if self._clahe_var.get() else "OFF"
        self._status.set(f"CLAHE display {mode}.")

    # ── image display ──────────────────────────────────────────────────────────

    def _refresh_images(self):
        ax0, ax1 = self._img_axes
        ax0.cla(); ax1.cla()
        ax0.set_facecolor("#111122"); ax1.set_facecolor("#111122")
        ax0.set_title("Channel 1", color=CYAN,    fontsize=10)
        ax1.set_title("Channel 2", color=MAGENTA, fontsize=10)

        if self.ch1 is not None:
            ch1_disp = self.ch1
            if self._clahe_var.get():
                ch1_disp = self._clahe_for_display(self.ch1)
                ax0.imshow(ch1_disp, cmap=_CYAN, origin="upper", vmin=0.0, vmax=1.0)
            elif self._histeq_var.get():
                ch1_disp = self._histeq_for_display(self.ch1)
                ax0.imshow(ch1_disp, cmap=_CYAN, origin="upper", vmin=0.0, vmax=1.0)
            elif self._autostretch_var.get():
                ch1_disp = self._autostretch_for_display(self.ch1)
                ax0.imshow(ch1_disp, cmap=_CYAN, origin="upper", vmin=0.0, vmax=1.0)
            else:
                ax0.imshow(ch1_disp, cmap=_CYAN, origin="upper", vmin=0, vmax=self.ch1.max())
        if self.ch2 is not None:
            ch2_disp = self.ch2
            if self._clahe_var.get():
                ch2_disp = self._clahe_for_display(self.ch2)
                ax1.imshow(ch2_disp, cmap=_MAGENTA, origin="upper", vmin=0.0, vmax=1.0)
            elif self._histeq_var.get():
                ch2_disp = self._histeq_for_display(self.ch2)
                ax1.imshow(ch2_disp, cmap=_MAGENTA, origin="upper", vmin=0.0, vmax=1.0)
            elif self._autostretch_var.get():
                ch2_disp = self._autostretch_for_display(self.ch2)
                ax1.imshow(ch2_disp, cmap=_MAGENTA, origin="upper", vmin=0.0, vmax=1.0)
            else:
                ax1.imshow(ch2_disp, cmap=_MAGENTA, origin="upper", vmin=0, vmax=self.ch2.max())

        if self.roi is not None:
            x1, y1, x2, y2 = self.roi
            # Channel 1 ROI is shown by the live RectangleSelector when active.
            if self._roi_kind == "rect" and not self._roi_active:
                ax0.add_patch(Rectangle(
                    (x1, y1), x2 - x1, y2 - y1,
                    linewidth=2, edgecolor=YELLOW, facecolor="none",
                ))
            if self._roi_kind == "rect":
                ax1.add_patch(Rectangle(
                    (x1, y1), x2 - x1, y2 - y1,
                    linewidth=2, edgecolor=YELLOW, facecolor="none",
                ))

        if self._roi_kind == "lasso" and self.roi_polygon is not None:
            poly = Polygon(self.roi_polygon, closed=True, fill=False,
                           edgecolor=YELLOW, linewidth=2)
            ax0.add_patch(poly)
            ax1.add_patch(Polygon(self.roi_polygon, closed=True, fill=False,
                                  edgecolor=YELLOW, linewidth=2))

        if self._roi_kind == "mask" and self.roi_mask is not None:
            self._add_mask_outline(ax0, self.roi_mask)
            self._add_mask_outline(ax1, self.roi_mask)

        for ax in self._img_axes:
            ax.axis("off")
        self._img_fig.tight_layout(pad=1.2)
        self._img_canvas.draw_idle()

    # ── ROI ────────────────────────────────────────────────────────────────────

    def toggle_roi(self):
        if self.ch1 is None:
            messagebox.showwarning("No image", "Load Channel 1 first.")
            return
        self._roi_active = not self._roi_active
        self._set_roi_selector_active()
        self._roi_btn.config(
            text="Cancel ROI" if self._roi_active else "Draw ROI"
        )
        if self._roi_active:
            self._status.set(
                f"{self._roi_mode_var.get()} ROI active: draw on the Channel-1 image."
            )

    def _on_roi_select(self, eclick, erelease):
        if eclick.xdata is None or erelease.xdata is None:
            return
        x1, x2 = sorted([eclick.xdata, erelease.xdata])
        y1, y2 = sorted([eclick.ydata, erelease.ydata])
        if self.ch1 is not None:
            h, w = self.ch1.shape[:2]
            x1 = max(0, int(x1)); x2 = min(w, int(x2))
            y1 = max(0, int(y1)); y2 = min(h, int(y2))
        if (x2 - x1) > 5 and (y2 - y1) > 5:
            self.roi = (x1, y1, x2, y2)
            self.roi_polygon = None
            self.roi_mask = None
            self._roi_kind = "rect"
            # Keep the selector active so the ROI remains editable via handles.
            self._roi_active = True
            self._set_roi_selector_active()
            self._roi_sel.extents = (x1, x2, y1, y2)
            self._roi_btn.config(text="Cancel ROI")
            self._schedule_secondary_roi_overlay()
            self._img_canvas.draw_idle()
        self._status.set(
            f"ROI: x=[{x1}, {x2}], y=[{y1}, {y2}]  —  Drag handles to refine or click 'Run Costes + Analyze'."
        )

    def _on_lasso_select(self, verts):
        if self.ch1 is None or self.ch2 is None:
            return
        if len(verts) < 3:
            return

        h, w = self.ch1.shape[:2]
        pts = np.asarray(verts, dtype=np.float64)
        path = MplPath(pts)

        yy, xx = np.mgrid[0:h, 0:w]
        pix = np.column_stack((xx.ravel(), yy.ravel()))
        mask = path.contains_points(pix).reshape(h, w)
        n_sel = int(mask.sum())
        if n_sel < 25:
            self._status.set("Lasso ROI too small. Draw a larger area.")
            return

        self.roi_polygon = pts
        self.roi_mask = mask
        self._roi_kind = "lasso"

        ys, xs = np.where(mask)
        self.roi = (int(xs.min()), int(ys.min()), int(xs.max()) + 1, int(ys.max()) + 1)

        self._roi_active = True
        self._set_roi_selector_active()
        self._roi_btn.config(text="Cancel ROI")
        self._refresh_images()
        self._status.set(
            f"Lasso ROI selected ({n_sel} pixels). Click 'Run Costes + Analyze'."
        )

    def clear_roi(self):
        self.roi = None
        self.roi_polygon = None
        self.roi_mask = None
        self._roi_kind = None
        self._roi_sel.clear()
        self._set_roi_selector_active()
        self._refresh_images()
        self._update_secondary_roi_overlay()
        self._status.set("ROI cleared — full image will be used.")

    def _schedule_secondary_roi_overlay(self):
        """Coalesce Ch2 ROI overlay updates to the next idle UI tick."""
        if self._roi_overlay_after_id is not None:
            try:
                self.root.after_cancel(self._roi_overlay_after_id)
            except tk.TclError:
                pass
        self._roi_overlay_after_id = self.root.after_idle(self._update_secondary_roi_overlay)

    def _update_secondary_roi_overlay(self):
        """Update the passive ROI rectangle on Channel 2 without resetting selector handles."""
        self._roi_overlay_after_id = None
        ax1 = self._img_axes[1]
        for p in list(ax1.patches):
            p.remove()
        for coll in list(ax1.collections):
            if getattr(coll, "_roi_mask_overlay", False):
                coll.remove()
        if self._roi_kind == "rect" and self.roi is not None:
            x1, y1, x2, y2 = self.roi
            ax1.add_patch(Rectangle(
                (x1, y1), x2 - x1, y2 - y1,
                linewidth=2, edgecolor=YELLOW, facecolor="none",
            ))
        elif self._roi_kind == "lasso" and self.roi_polygon is not None:
            ax1.add_patch(Polygon(self.roi_polygon, closed=True, fill=False,
                                  edgecolor=YELLOW, linewidth=2))
        elif self._roi_kind == "mask" and self.roi_mask is not None:
            self._add_mask_outline(ax1, self.roi_mask)
        self._img_canvas.draw_idle()

    def _add_mask_outline(self, ax, mask: np.ndarray):
        contour = ax.contour(
            mask.astype(float),
            levels=[0.5],
            colors=[YELLOW],
            linewidths=2,
            origin="upper",
        )
        for coll in contour.collections:
            coll._roi_mask_overlay = True

    # ── data extraction ────────────────────────────────────────────────────────

    @staticmethod
    def _despike_hot_pixels(img: np.ndarray) -> np.ndarray:
        """Suppress isolated bright outliers using a conservative 3x3 local test."""
        if img.size == 0:
            return img

        med = ndimage.median_filter(img, size=3, mode="reflect")
        diff = img - med
        local_mad = ndimage.median_filter(np.abs(diff), size=3, mode="reflect")

        p1, p99 = np.percentile(img, [1.0, 99.0])
        robust_range = max(float(p99 - p1), 1.0)
        threshold = (6.0 * local_mad) + (0.02 * robust_range)
        mask = diff > threshold

        if not np.any(mask):
            return img

        out = img.copy()
        out[mask] = med[mask]
        return out

    def _get_roi_arrays(self, preprocess: bool = False):
        if self.ch1 is None or self.ch2 is None:
            return None, None
        if self.ch1.shape[:2] != self.ch2.shape[:2]:
            messagebox.showerror(
                "Shape mismatch",
                "Both channels must have identical spatial dimensions.\n"
                f"Ch1: {self.ch1.shape}   Ch2: {self.ch2.shape}",
            )
            return None, None
        if self.roi and self._roi_kind == "rect":
            x1, y1, x2, y2 = self.roi
            c1 = self.ch1[y1:y2, x1:x2]
            c2 = self.ch2[y1:y2, x1:x2]
        else:
            c1 = self.ch1
            c2 = self.ch2

        c1, c2 = self._apply_background_subtraction(c1, c2)

        if preprocess and self._despike_var.get():
            c1 = self._despike_hot_pixels(c1)
            c2 = self._despike_hot_pixels(c2)

        return c1, c2

    def _get_pixels(self, preprocess: bool = False):
        if self.ch1 is None or self.ch2 is None:
            return None, None
        if self.ch1.shape[:2] != self.ch2.shape[:2]:
            messagebox.showerror(
                "Shape mismatch",
                "Both channels must have identical spatial dimensions.\n"
                f"Ch1: {self.ch1.shape}   Ch2: {self.ch2.shape}",
            )
            return None, None

        if self._roi_kind in {"lasso", "mask"} and self.roi_mask is not None:
            c1_img = self.ch1
            c2_img = self.ch2
            c1_img, c2_img = self._apply_background_subtraction(c1_img, c2_img)
            if preprocess and self._despike_var.get():
                c1_img = self._despike_hot_pixels(c1_img)
                c2_img = self._despike_hot_pixels(c2_img)
            c1 = c1_img[self.roi_mask]
            c2 = c2_img[self.roi_mask]
        else:
            c1_img, c2_img = self._get_roi_arrays(preprocess=preprocess)
            if c1_img is None:
                return None, None
            c1 = c1_img.ravel()
            c2 = c2_img.ravel()
        return c1, c2

    def _get_randomization_arrays(self, preprocess: bool = False):
        """Return 2D analysis arrays plus an optional ROI mask for randomization."""
        if self.ch1 is None or self.ch2 is None:
            return None, None, None
        if self.ch1.shape[:2] != self.ch2.shape[:2]:
            messagebox.showerror(
                "Shape mismatch",
                "Both channels must have identical spatial dimensions.\n"
                f"Ch1: {self.ch1.shape}   Ch2: {self.ch2.shape}",
            )
            return None, None, None

        if self._roi_kind in {"lasso", "mask"} and self.roi_mask is not None:
            c1_img = self.ch1
            c2_img = self.ch2
            roi_mask = self.roi_mask
        elif self.roi and self._roi_kind == "rect":
            x1, y1, x2, y2 = self.roi
            c1_img = self.ch1[y1:y2, x1:x2]
            c2_img = self.ch2[y1:y2, x1:x2]
            roi_mask = None
        else:
            c1_img = self.ch1
            c2_img = self.ch2
            roi_mask = None

        c1_img, c2_img = self._apply_background_subtraction(c1_img, c2_img)
        if preprocess and self._despike_var.get():
            c1_img = self._despike_hot_pixels(c1_img)
            c2_img = self._despike_hot_pixels(c2_img)

        return c1_img, c2_img, roi_mask

    # ── Costes algorithm ───────────────────────────────────────────────────────

    def _costes_debug_log_path(self) -> Path:
        return Path(__file__).resolve().parent / "costes_threshold_debug.log"

    def _costes_debug_log_write(self, message: str):
        """Append a timestamped line to the Costes debug log."""
        if not self._costes_debug_log_var.get():
            return
        ts = datetime.now().isoformat(timespec="milliseconds")
        try:
            with self._costes_debug_log_path().open("a", encoding="utf-8") as fh:
                fh.write(f"[{ts}] {message}\n")
        except Exception:
            # Logging must not interfere with analysis.
            pass

    def _clear_costes_debug_log(self):
        """Delete the Costes debug log file used for troubleshooting runs."""
        log_path = self._costes_debug_log_path()
        try:
            if log_path.exists():
                log_path.unlink()
                self._status.set("Costes debug log cleared.")
            else:
                self._status.set("Costes debug log is already empty.")
        except Exception as exc:
            self._status.set(f"Could not clear Costes debug log: {exc}")

    @staticmethod
    def _orthogonal_regression(c1: np.ndarray, c2: np.ndarray):
        """
        Orthogonal (total least squares) regression of Ch2 on Ch1.

        Minimises perpendicular distances to the fitted line rather than
        vertical residuals. Equivalent to the first principal component of
        the (c1, c2) scatter cloud, forced through the means.

        Returns (slope, intercept) in the same convention as linregress.
        """
        x = c1.astype(np.float64)
        y = c2.astype(np.float64)
        mx, my = x.mean(), y.mean()
        xc, yc = x - mx, y - my
        # 2×2 scatter matrix
        sxx = float((xc * xc).sum())
        sxy = float((xc * yc).sum())
        syy = float((yc * yc).sum())
        # Slope from first eigenvector of the scatter matrix
        b = (syy - sxx + np.sqrt((syy - sxx) ** 2 + 4.0 * sxy ** 2)) / (2.0 * sxy) if abs(sxy) > 1e-12 else (1.0 if syy >= sxx else 0.0)
        a = my - b * mx
        return float(b), float(a)

    @staticmethod
    def costes_threshold(c1: np.ndarray, c2: np.ndarray,
                         n_steps: int = 1024,
                         orthogonal: bool = False,
                         debug_hook=None,
                         rise_stop_fraction: float = 0.01):
        """
        Costes automatic threshold algorithm (Costes et al., 2004).

        Steps
        -----
        1. Fit regression  Ch2 = a·Ch1 + b  (OLS or orthogonal, see `orthogonal`).
          2. Starting from channel max, decrease the working threshold by 1
              intensity unit per step (Coloc2 SimpleStepper). If |a| < 1 the
              working threshold is Ch1, otherwise Ch2; the paired threshold is
              mapped via the regression line.
        3. At each step compute Pearson's r for pixels BELOW the threshold
           pair: c1 < T1 OR c2 < T2 (union / background region, per Costes 2004
           and Coloc2).  These are the putative background
           pixels; they should be uncorrelated (r ≤ 0) when the threshold
           sits above all real signal.
          4. Stop when r_background becomes very small (< 1e-4), becomes NaN,
              or rises above the running minimum r by at least the configured
              percentage.

        Parameters
        ----------
        orthogonal : bool
            When True use orthogonal (total least squares) regression instead
            of ordinary least squares.  This treats both channels symmetrically
            and minimises perpendicular rather than vertical residuals.

        Returns
        -------
        t1, t2 : optimal threshold values
        slope, intercept : regression coefficients
        curve_t1, curve_t2, curve_r : threshold-pair scan and Pearson r of
            below-threshold pixels at each step
        """
        def _emit(msg: str):
            if debug_hook is None:
                return
            try:
                debug_hook(msg)
            except Exception:
                pass

        if orthogonal:
            slope, intercept = ColocalizationApp._orthogonal_regression(c1, c2)
            reg_mode = "orthogonal"
        else:
            res = stats.linregress(c1, c2)
            slope, intercept = float(res.slope), float(res.intercept)
            reg_mode = "ols"

        _emit(
            "REGRESSION "
            f"mode={reg_mode} slope={slope:.8g} intercept={intercept:.8g} "
            f"n_pixels={len(c1)} rise_stop_fraction={rise_stop_fraction:.8g}"
        )

        max1, min1 = float(c1.max()), float(c1.min())
        max2, min2 = float(c2.max()), float(c2.min())

        curve_t1, curve_t2, curve_r = [], [], []
        best_r_any = float("inf")
        t1_opt = float(np.percentile(c1, 95))
        t2_opt = float(np.clip(slope * t1_opt + intercept, min2, max2))
        t1_opt_any = t1_opt
        t2_opt_any = t2_opt

        # Coloc2 behavior: step on Ch1 if -1 < slope < 1, otherwise on Ch2.
        step_on_ch1 = (-1.0 < slope < 1.0)
        work_t = max1 if step_on_ch1 else max2
        if step_on_ch1:
            _emit(
                "STEP_AXIS driver=ch1 driven_threshold=T1 paired_threshold=T2 "
                f"rule='-1 < slope < 1' slope={slope:.8g} start_work_t={work_t:.6g}"
            )
        else:
            _emit(
                "STEP_AXIS driver=ch2 driven_threshold=T2 paired_threshold=T1 "
                f"rule='slope <= -1 or slope >= 1' slope={slope:.8g} start_work_t={work_t:.6g}"
            )

        # Coloc2 SimpleStepper defaults.
        current_r = 1.0
        last_r = float("inf")
        finished = False
        step_idx = 0
        min_r_seen = float("inf")
        # Track the last step where r was still >= 0.0001 (above the noise floor).
        # When we stop for any reason, this is what we return as the threshold.
        last_ok_t1 = t1_opt
        last_ok_t2 = t2_opt
        last_ok_r = float("nan")

        while not finished:
            if step_on_ch1:
                t1_raw = work_t
                t2_raw = slope * work_t + intercept
            else:
                t2_raw = work_t
                if abs(slope) < 1e-12:
                    t1_raw = max1
                else:
                    t1_raw = (work_t - intercept) / slope

            # Coloc2 rounds thresholds to integer image levels.
            t1_cand = float(np.clip(np.floor(t1_raw + 0.5), min1, max1))
            t2_cand = float(np.clip(np.floor(t2_raw + 0.5), min2, max2))

            # ThresholdMode.Below in Coloc2: ch1 < T1 OR ch2 < T2.
            mask = (c1 < t1_cand) | (c2 < t2_cand)
            if mask.sum() >= 3:
                c1b, c2b = c1[mask], c2[mask]
                if c1b.std() >= 1e-9 and c2b.std() >= 1e-9:
                    r = float(stats.pearsonr(c1b, c2b)[0])
                else:
                    r = float("nan")
            else:
                r = float("nan")

            curve_t1.append(t1_cand)
            curve_t2.append(t2_cand)
            curve_r.append(r)

            _emit(
                "STEP "
                f"idx={step_idx} work_t={work_t:.6g} t1={t1_cand:.6g} t2={t2_cand:.6g} "
                f"mask_n={int(mask.sum())} r={r:.8g} prev_r={last_r:.8g} min_r={min_r_seen:.8g}"
            )

            if np.isfinite(r) and r < best_r_any:
                best_r_any = r
                t1_opt_any = t1_cand
                t2_opt_any = t2_cand

            if np.isfinite(r) and r < min_r_seen:
                min_r_seen = r

            # Keep a running record of the last step where r was still above
            # the noise floor (>= 0.0001).  This is the threshold we will
            # return regardless of which stop condition fires.
            if np.isfinite(r) and r >= 0.0001:
                last_ok_t1 = t1_cand
                last_ok_t2 = t2_cand
                last_ok_r = r

            # Emulate Coloc2 SimpleStepper update() and termination criteria.
            last_r, current_r = current_r, r
            next_work_t = work_t - 1.0
            stop_reasons = []
            if not np.isfinite(r):
                stop_reasons.append("non_finite_r")
            if next_work_t < 1.0:
                stop_reasons.append("next_work_t_below_1")
            if np.isfinite(r) and r < 0.0001:
                stop_reasons.append("r_below_0.0001")
            min_r_excess = rise_stop_fraction * abs(min_r_seen) if np.isfinite(min_r_seen) else float("inf")
            if np.isfinite(min_r_seen) and np.isfinite(r) and (r - min_r_seen) > min_r_excess:
                stop_reasons.append("r_above_min_r")
            finished = bool(stop_reasons)

            if finished:
                # Always return the last threshold where r was still >= 0.0001.
                t1_opt = last_ok_t1
                t2_opt = last_ok_t2
                _emit(
                    "STOP "
                    f"idx={step_idx} reasons={','.join(stop_reasons)} "
                    f"final_t1={t1_opt:.6g} final_t2={t2_opt:.6g} "
                    f"last_ok_r={last_ok_r:.8g} current_r={r:.8g} "
                    f"min_r={min_r_seen:.8g} next_work_t={next_work_t:.6g}"
                )

            work_t = next_work_t
            step_idx += 1

        # Fallback: if all evaluated r were non-finite, keep conservative defaults.
        if not np.isfinite(best_r_any):
            t1_opt, t2_opt = t1_opt_any, t2_opt_any
            _emit(
                "FALLBACK all_r_non_finite=true "
                f"fallback_t1={t1_opt:.6g} fallback_t2={t2_opt:.6g}"
            )

        _emit(
            "RESULT "
            f"t1={t1_opt:.6g} t2={t2_opt:.6g} best_r_any={best_r_any:.8g} "
            f"curve_len={len(curve_t1)}"
        )

        return (
            t1_opt, t2_opt,
            slope, intercept,
            np.asarray(curve_t1), np.asarray(curve_t2), np.asarray(curve_r),
        )

    # ── Mander's coefficients ──────────────────────────────────────────────────

    @staticmethod
    def manders(c1: np.ndarray, c2: np.ndarray,
                t1: float, t2: float) -> dict:
        """
        Compute Mander's Overlap Coefficients and related metrics.

        M1  : fraction of Ch1 signal that overlaps with Ch2 signal above T₂
        M2  : fraction of Ch2 signal that overlaps with Ch1 signal above T₁
        tM1 : thresholded M1 – only pixels where BOTH channels exceed threshold
        tM2 : thresholded M2
        overlap : Manders (1993) overlap coefficient R
        k1, k2  : partitioning coefficients
        pearson : Pearson's r restricted to pixels above both thresholds
        """
        s1, s2 = c1.sum(), c2.sum()
        mask2     = c2 > t2
        mask1     = c1 > t1
        mask_both = mask1 & mask2

        m1  = c1[mask2].sum() / s1     if s1 > 0 else 0.0
        m2  = c2[mask1].sum() / s2     if s2 > 0 else 0.0
        tm1 = c1[mask_both].sum() / s1 if s1 > 0 else 0.0
        tm2 = c2[mask_both].sum() / s2 if s2 > 0 else 0.0

        # Overlap R and k-coefficients (positive pixels only)
        pos  = (c1 > 0) & (c2 > 0)
        a, b = c1[pos].astype(np.float64), c2[pos].astype(np.float64)
        dR   = np.sqrt((a ** 2).sum() * (b ** 2).sum())
        d1   = (a ** 2).sum()
        d2   = (b ** 2).sum()
        ab   = (a * b).sum()
        overlap = ab / dR if dR > 0 else 0.0
        k1      = ab / d1 if d1 > 0 else 0.0
        k2      = ab / d2 if d2 > 0 else 0.0

        # Pearson r above thresholds
        c1b, c2b = c1[mask_both], c2[mask_both]
        if len(c1b) > 5 and c1b.std() > 0 and c2b.std() > 0:
            pearson = float(stats.pearsonr(c1b, c2b)[0])
        else:
            pearson = float("nan")

        return dict(m1=m1, m2=m2, tm1=tm1, tm2=tm2,
                    overlap=overlap, k1=k1, k2=k2, pearson=pearson)

    @staticmethod
    def costes_randomization(c1_img: np.ndarray, c2_img: np.ndarray,
                             n_iter: int = 100, block_size: int | None = None,
                             seed: int = 12345,
                             roi_mask: np.ndarray | None = None) -> dict:
        """
        Costes randomization significance test.

        Channel 2 is partitioned into square blocks and block order is
        randomized repeatedly. This preserves local intensity structure within
        each block while destroying cross-channel spatial correspondence.
        The block size should reflect the microscope PSF scale (in pixels), as
        done in Coloc2. When `roi_mask` is provided, the shuffle is performed
        within the ROI bounding box and Pearson's r is evaluated only on pixels
        selected by the mask.
        """
        h, w = c1_img.shape
        mask = None
        if roi_mask is not None:
            mask = np.asarray(roi_mask, dtype=bool)
            if mask.shape != c1_img.shape:
                raise ValueError("roi_mask must match the image shape for Costes randomization")
            if not np.any(mask):
                return dict(observed_r=float("nan"), random_mean=float("nan"),
                            random_std=float("nan"), significance=float("nan"),
                            p_value=float("nan"), block_size=0, n_iter=0)

            ys, xs = np.where(mask)
            y0, y1 = int(ys.min()), int(ys.max()) + 1
            x0, x1 = int(xs.min()), int(xs.max()) + 1
            c1_img = c1_img[y0:y1, x0:x1]
            c2_img = c2_img[y0:y1, x0:x1]
            mask = mask[y0:y1, x0:x1]
            h, w = c1_img.shape

        min_dim = min(h, w)
        if min_dim < 8:
            return dict(observed_r=float("nan"), random_mean=float("nan"),
                        random_std=float("nan"), significance=float("nan"),
                        p_value=float("nan"), block_size=0, n_iter=0)

        if block_size is None:
            if min_dim >= 128:
                block_size = 8
            elif min_dim >= 64:
                block_size = 4
            else:
                block_size = 2

        block_size = max(2, min(block_size, min_dim))
        crop_h = (h // block_size) * block_size
        crop_w = (w // block_size) * block_size
        if crop_h < block_size or crop_w < block_size:
            return dict(observed_r=float("nan"), random_mean=float("nan"),
                        random_std=float("nan"), significance=float("nan"),
                        p_value=float("nan"), block_size=block_size, n_iter=0)

        if mask is None:
            c1_crop = c1_img[:crop_h, :crop_w]
            c2_crop = c2_img[:crop_h, :crop_w]
            mask_crop = None
        else:
            best_score = -1
            best_offsets = (0, 0)
            max_y0 = h - crop_h
            max_x0 = w - crop_w
            for y_off in range(max_y0 + 1):
                for x_off in range(max_x0 + 1):
                    score = int(mask[y_off:y_off + crop_h, x_off:x_off + crop_w].sum())
                    if score > best_score:
                        best_score = score
                        best_offsets = (y_off, x_off)

            if best_score < 3:
                return dict(observed_r=float("nan"), random_mean=float("nan"),
                            random_std=float("nan"), significance=float("nan"),
                            p_value=float("nan"), block_size=block_size, n_iter=0)

            y_off, x_off = best_offsets
            c1_crop = c1_img[y_off:y_off + crop_h, x_off:x_off + crop_w]
            c2_crop = c2_img[y_off:y_off + crop_h, x_off:x_off + crop_w]
            mask_crop = mask[y_off:y_off + crop_h, x_off:x_off + crop_w]

        if mask_crop is None:
            c1_flat = c1_crop.ravel()
            c2_flat = c2_crop.ravel()
        else:
            c1_flat = c1_crop[mask_crop]
            c2_flat = c2_crop[mask_crop]

        if c1_flat.std() < 1e-9 or c2_flat.std() < 1e-9:
            return dict(observed_r=float("nan"), random_mean=float("nan"),
                        random_std=float("nan"), significance=float("nan"),
                        p_value=float("nan"), block_size=block_size, n_iter=0)

        observed_r = float(stats.pearsonr(c1_flat, c2_flat)[0])

        ny = crop_h // block_size
        nx = crop_w // block_size
        blocks = c2_crop.reshape(ny, block_size, nx, block_size).transpose(0, 2, 1, 3)
        blocks = blocks.reshape(ny * nx, block_size, block_size)

        rng = np.random.default_rng(seed)
        randomized_r = np.empty(n_iter, dtype=np.float64)
        for idx in range(n_iter):
            perm = rng.permutation(blocks.shape[0])
            shuffled = blocks[perm].reshape(ny, nx, block_size, block_size)
            shuffled = shuffled.transpose(0, 2, 1, 3).reshape(crop_h, crop_w)
            shuffled_flat = shuffled.ravel() if mask_crop is None else shuffled[mask_crop]
            randomized_r[idx] = float(stats.pearsonr(c1_flat, shuffled_flat)[0])

        random_mean = float(randomized_r.mean())
        random_std = float(randomized_r.std(ddof=1)) if n_iter > 1 else 0.0
        significance = 100.0 * float(np.mean(randomized_r < observed_r))
        p_value = (float(np.sum(randomized_r >= observed_r)) + 1.0) / (n_iter + 1.0)

        return dict(
            observed_r=observed_r,
            random_mean=random_mean,
            random_std=random_std,
            significance=significance,
            p_value=p_value,
            block_size=block_size,
            n_iter=n_iter,
        )

    # ── analysis entry point ───────────────────────────────────────────────────

    def run_analysis(self):
        use_preprocess = self._despike_var.get()

        if self._roi_kind == "lasso":
            c1, c2 = self._get_pixels(preprocess=use_preprocess)
            c1_img = None
            c2_img = None
        else:
            c1_img, c2_img = self._get_roi_arrays(preprocess=use_preprocess)
            if c1_img is None:
                messagebox.showwarning("Missing data", "Load both channel images first.")
                return
            c1 = c1_img.ravel()
            c2 = c2_img.ravel()

        if c1 is None:
            messagebox.showwarning("Missing data", "Load both channel images first.")
            return

        self._status.set("Running Costes algorithm…")
        self.root.update_idletasks()

        debug_hook = self._costes_debug_log_write if self._costes_debug_log_var.get() else None
        rise_stop_fraction = self._get_costes_rise_stop_fraction()
        if debug_hook is not None:
            debug_hook(
                "RUN_START "
                f"orthogonal={self._orthreg_var.get()} preprocess_hot_pixels={use_preprocess} "
                f"roi_kind={self._roi_kind if self._roi_kind else 'full'} n_pixels={len(c1)} "
                f"rise_stop_pct={100.0 * rise_stop_fraction:.4f}"
            )

        t1, t2, slope, intercept, curve_t1, curve_t2, curve_r = self.costes_threshold(
            c1, c2,
            orthogonal=self._orthreg_var.get(),
            debug_hook=debug_hook,
            rise_stop_fraction=rise_stop_fraction,
        )
        if debug_hook is not None:
            debug_hook(
                "RUN_END "
                f"t1={t1:.6g} t2={t2:.6g} slope={slope:.8g} intercept={intercept:.8g}"
            )
        self._last_costes_slope = slope
        self._last_costes_intercept = intercept
        self._last_costes_curve_t1 = curve_t1
        self._last_costes_curve_t2 = curve_t2
        self._last_costes_curve_r = curve_r
        res = self.manders(c1, c2, t1, t2)
        self._status.set("Running Costes randomization…")
        self.root.update_idletasks()
        rand_c1_img, rand_c2_img, rand_mask = self._get_randomization_arrays(
            preprocess=use_preprocess
        )
        rand = self.costes_randomization(
            rand_c1_img,
            rand_c2_img,
            block_size=self._get_psf_block_size(),
            roi_mask=rand_mask,
        )

        # update threshold entry fields
        self._t1_entry_var.set(f"{t1:.1f}")
        self._t2_entry_var.set(f"{t2:.1f}")

        # update results
        self._rv["costes_t1"].set(f"{t1:.2f}")
        self._rv["costes_t2"].set(f"{t2:.2f}")
        for k in ("m1", "m2", "tm1", "tm2", "overlap", "k1", "k2"):
            self._rv[k].set(f"{res[k]:.4f}")
        p = res["pearson"]
        self._rv["pearson"].set(f"{p:.4f}" if not np.isnan(p) else "N/A")
        obs_r = rand["observed_r"]
        if np.isnan(obs_r):
            self._rv["costes_r_obs"].set("N/A")
            self._rv["costes_r_rand"].set("N/A")
            self._rv["costes_sig"].set("N/A")
            self._rv["costes_p"].set("N/A")
        else:
            self._rv["costes_r_obs"].set(
                f"r={obs_r:.4f}  (block={rand['block_size']}, n={rand['n_iter']})"
            )
            self._rv["costes_r_rand"].set(
                f"{rand['random_mean']:.4f} ± {rand['random_std']:.4f}"
            )
            self._rv["costes_sig"].set(f"{rand['significance']:.1f} %")
            self._rv["costes_p"].set(f"{rand['p_value']:.4f}")

        # plots
        self._draw_histograms(c1, c2, t1, t2)
        self._draw_costes(c1, c2, t1, t2, slope, intercept, curve_t1, curve_t2, curve_r)

        self._status.set(
            f"Done  |  T₁={t1:.1f}  T₂={t2:.1f}  |  "
            f"M1={res['m1']:.3f}  M2={res['m2']:.3f}  "
            f"tM1={res['tm1']:.3f}  tM2={res['tm2']:.3f}  |  "
            f"Costes sig={rand['significance']:.1f}%"
            + ("  |  Orthogonal reg" if self._orthreg_var.get() else "")
            + ("  |  hot-pixel preprocess ON" if use_preprocess else "")
            + (
                f"  |  BG sub Ch1={self._bg1_var.get()}, Ch2={self._bg2_var.get()}"
                if (float(self._bg1_var.get()) > 0.0 or float(self._bg2_var.get()) > 0.0)
                else ""
            )
        )

    def _slider_update(self, _=None, draw_plots: bool = True):
        try:
            t1 = float(self._t1_entry_var.get())
        except ValueError:
            t1 = 0.0
        try:
            t2 = float(self._t2_entry_var.get())
        except ValueError:
            t2 = 0.0
        t1 = max(0.0, min(self._t1_max, t1))
        t2 = max(0.0, min(self._t2_max, t2))
        self._t1_entry_var.set(f"{t1:.1f}")
        self._t2_entry_var.set(f"{t2:.1f}")
        self._rv["costes_t1"].set(f"{t1:.2f}")
        self._rv["costes_t2"].set(f"{t2:.2f}")
        c1, c2 = self._get_pixels(preprocess=self._despike_var.get())
        if c1 is None:
            return
        res = self.manders(c1, c2, t1, t2)
        for k in ("m1", "m2", "tm1", "tm2", "overlap", "k1", "k2"):
            self._rv[k].set(f"{res[k]:.4f}")
        p = res["pearson"]
        self._rv["pearson"].set(f"{p:.4f}" if not np.isnan(p) else "N/A")
        if draw_plots:
            self._draw_histograms(c1, c2, t1, t2)
            if self._last_costes_slope is not None and self._last_costes_intercept is not None:
                self._draw_costes(
                    c1, c2, t1, t2,
                    self._last_costes_slope,
                    self._last_costes_intercept,
                    self._last_costes_curve_t1,
                    self._last_costes_curve_t2,
                    self._last_costes_curve_r,
                )

        self._status.set(
            f"Manual thresholds applied  |  T₁={t1:.1f}  T₂={t2:.1f}  |  "
            f"M1={res['m1']:.3f}  M2={res['m2']:.3f}  "
            f"tM1={res['tm1']:.3f}  tM2={res['tm2']:.3f}"
            + (
                f"  |  BG sub Ch1={self._bg1_var.get()}, Ch2={self._bg2_var.get()}"
                if (float(self._bg1_var.get()) > 0.0 or float(self._bg2_var.get()) > 0.0)
                else ""
            )
        )

    # ── plots ──────────────────────────────────────────────────────────────────

    def _style_ax(self, ax, title="", xlabel="", ylabel="", title_color=FG):
        ax.set_facecolor("#11111e")
        ax.tick_params(colors=FG, labelsize=8)
        for sp in ax.spines.values():
            sp.set_edgecolor(GRAY)
        ax.set_title(title,  color=title_color, fontsize=10)
        ax.set_xlabel(xlabel, color=FG,          fontsize=9)
        ax.set_ylabel(ylabel, color=FG,          fontsize=9)

    def _draw_histograms(self, c1: np.ndarray, c2: np.ndarray,
                         t1: float, t2: float):
        ax1, ax2 = self._hist_axes
        ax1.cla(); ax2.cla()

        def _hist_one(ax, data, color, title, thresh):
            n_bins = min(512, max(64, int(data.max() - data.min()) + 1))
            # all pixels
            ax.hist(data, bins=n_bins, color=color, alpha=0.70,
                    log=True, linewidth=0, label="All pixels")
            # pixels above threshold (highlight)
            above = data[data >= thresh]
            pct   = 100 * len(above) / max(len(data), 1)
            if len(above):
                ax.hist(above, bins=n_bins, color="white", alpha=0.30,
                        log=True, linewidth=0,
                        label=f"Above threshold  ({pct:.1f} %)")
            threshold_line = ax.axvline(
                thresh, color=YELLOW, linewidth=2, linestyle="--",
                label=f"T = {thresh:.1f}"
            )
            self._style_ax(ax, title=title,
                           xlabel="Intensity", ylabel="Pixel count (log)",
                           title_color=color)
            ax.minorticks_on()
            ax.grid(axis="x", which="major", color=GRAY, alpha=0.26, linestyle="-", linewidth=0.6)
            ax.grid(axis="x", which="minor", color=GRAY, alpha=0.14, linestyle="-", linewidth=0.35)
            ax.legend(fontsize=7, facecolor="#111122",
                      labelcolor="white", framealpha=0.9)
            return threshold_line

        self._hist_t1_line = _hist_one(ax1, c1, CYAN,    "Channel 1  Intensity Histogram", t1)
        self._hist_t2_line = _hist_one(ax2, c2, MAGENTA, "Channel 2  Intensity Histogram", t2)

        self._hist_fig.tight_layout()
        self._hist_canvas.draw_idle()

    def _update_hist_threshold_lines(self, t1: float, t2: float):
        """Move only the threshold guide lines without recomputing plots."""
        if self._hist_t1_line is None or self._hist_t2_line is None:
            return
        self._hist_t1_line.set_xdata([t1, t1])
        self._hist_t2_line.set_xdata([t2, t2])
        self._hist_t1_line.set_label(f"T = {t1:.1f}")
        self._hist_t2_line.set_label(f"T = {t2:.1f}")
        self._hist_canvas.draw_idle()

    def _draw_costes(self, c1: np.ndarray, c2: np.ndarray,
                     t1: float, t2: float,
                     slope: float, intercept: float,
                     curve_t1: np.ndarray, curve_t2: np.ndarray,
                     curve_r: np.ndarray):

        ax_sc, ax_r = self._costes_axes
        ax_sc.cla(); ax_r.cla()

        # ── scatter plot ───────────────────────────────────────────────────────
        n = len(c1)
        rng = np.random.default_rng(42)
        if n > 30_000:
            idx = rng.choice(n, 30_000, replace=False)
            sc1, sc2 = c1[idx], c2[idx]
        else:
            sc1, sc2 = c1, c2

        above = (sc1 >= t1) & (sc2 >= t2)
        ax_sc.scatter(sc1[~above], sc2[~above], c="#2d4a5e", s=1,
                      alpha=0.4, rasterized=True, label="Non-colocalizing")
        ax_sc.scatter(sc1[above],  sc2[above],  c="gold",    s=2,
                      alpha=0.8, rasterized=True, label="Colocalizing")

        # regression line
        xr = np.array([c1.min(), c1.max()])
        reg_label = "Orthogonal" if self._orthreg_var.get() else "OLS"
        ax_sc.plot(xr, slope * xr + intercept, color="tomato",
                   linewidth=1.5, label=f"{reg_label} regression  a={slope:.3f}  b={intercept:.3f}")

        # threshold crosshairs
        ax_sc.axvline(t1, color=CYAN,    linewidth=1.5, linestyle="--",
                      label=f"T₁ = {t1:.1f}")
        ax_sc.axhline(t2, color=MAGENTA, linewidth=1.5, linestyle="--",
                      label=f"T₂ = {t2:.1f}")

        # shade colocalization quadrant
        ax_sc.axvspan(t1, c1.max(), alpha=0.06, color=CYAN)
        ax_sc.axhspan(t2, c2.max(), alpha=0.06, color=MAGENTA)

        self._style_ax(ax_sc,
                       title="Pixel Scatter  (Costes thresholds)",
                       xlabel="Channel 1 Intensity",
                       ylabel="Channel 2 Intensity")
        ax_sc.legend(fontsize=7, facecolor="#111122",
                     labelcolor="white", framealpha=0.9, markerscale=4)

        # ── Costes r-vs-threshold curve ────────────────────────────────────────
        if len(curve_t1) > 1:
            curve_r_arr = np.asarray(curve_r)
            ax_r.plot(curve_t1, curve_r_arr, color="white",
                      linewidth=1.8, label="r(background)  vs  T₁")
            # r > 0: background pixels still correlate → threshold too high
            ax_r.fill_between(curve_t1, curve_r_arr, 0,
                               where=(curve_r_arr > 0),
                               color="steelblue", alpha=0.30,
                               label="r > 0  (signal leaks into background)")
            # r ≤ 0: background is uncorrelated → threshold is at or above signal
            ax_r.fill_between(curve_t1, curve_r_arr, 0,
                               where=(curve_r_arr <= 0),
                               color="#a6e3a1", alpha=0.35,
                               label="r ≤ 0  (background uncorrelated ✓)")
            ax_r.axhline(0, color="red", linewidth=1.2, linestyle=":")
            ax_r.axvline(t1, color=YELLOW, linewidth=2, linestyle="-.",
                         label=f"Costes thresholds = ({t1:.1f}, {t2:.1f})")
            self._style_ax(ax_r,
                           title="Costes:  r of BELOW-threshold pixels vs T₁\n"
                                 "(paired T₂ shown on top axis)",
                           xlabel="Ch1 Threshold  (decreasing →)",
                           ylabel="Pearson's r  (background pixels)")

            secax = ax_r.secondary_xaxis("top")
            sample_idx = np.linspace(0, len(curve_t1) - 1, min(6, len(curve_t1)), dtype=int)
            sample_idx = np.unique(sample_idx)
            tick_pairs = []
            for idx in sample_idx:
                pos = float(curve_t1[idx])
                if any(abs(pos - existing_pos) < 1e-9 for existing_pos, _ in tick_pairs):
                    continue
                tick_pairs.append((pos, f"{float(curve_t2[idx]):.1f}"))
            if tick_pairs:
                secax.set_xticks([pos for pos, _ in tick_pairs])
                secax.set_xticklabels([label for _, label in tick_pairs])
            secax.set_xlabel("Paired Ch2 Threshold  (T₂)")
            secax.tick_params(axis="x", colors=FG, labelsize=8)
            secax.xaxis.label.set_color(FG)

            ax_r.legend(fontsize=7, facecolor="#111122",
                        labelcolor="white", framealpha=0.9)
            ax_r.invert_xaxis()   # show threshold decreasing left to right
        else:
            self._style_ax(ax_r, title="Costes curve  (insufficient data)",
                           xlabel="Ch1 Threshold", ylabel="Pearson's r")

        self._costes_fig.tight_layout()
        self._costes_canvas.draw_idle()

    # ── PDF export ────────────────────────────────────────────────────────────

    def save_pdf_summary(self):
        """Export a PDF report with current metrics and visualizations."""
        if self.ch1 is None or self.ch2 is None:
            messagebox.showwarning("Missing data", "Load both channel images first.")
            return

        default_name = f"colocalization_summary_{datetime.now():%Y%m%d_%H%M%S}.pdf"
        out_path = filedialog.asksaveasfilename(
            title="Save PDF Summary",
            defaultextension=".pdf",
            initialfile=default_name,
            filetypes=[("PDF", "*.pdf")],
        )
        if not out_path:
            return

        if not out_path.lower().endswith(".pdf"):
            out_path += ".pdf"

        if self._roi_kind == "lasso" and self.roi_mask is not None and self.roi is not None:
            x1, y1, x2, y2 = self.roi
            roi_text = (
                f"Lasso ({int(self.roi_mask.sum())} px), bbox: x=[{x1}, {x2}], y=[{y1}, {y2}]"
            )
        elif self.roi is None:
            roi_text = "Full image"
        else:
            x1, y1, x2, y2 = self.roi
            roi_text = f"x=[{x1}, {x2}], y=[{y1}, {y2}]"

        has_analysis = self._rv["costes_t1"].get() != "—"
        metrics = [
            ("Costes T1", self._rv["costes_t1"].get()),
            ("Costes T2", self._rv["costes_t2"].get()),
            ("M1", self._rv["m1"].get()),
            ("M2", self._rv["m2"].get()),
            ("tM1", self._rv["tm1"].get()),
            ("tM2", self._rv["tm2"].get()),
            ("Pearson r", self._rv["pearson"].get()),
            ("Overlap R", self._rv["overlap"].get()),
            ("k1", self._rv["k1"].get()),
            ("k2", self._rv["k2"].get()),
            ("Costes observed r", self._rv["costes_r_obs"].get()),
            ("Randomized r mean ± SD", self._rv["costes_r_rand"].get()),
            ("Costes significance", self._rv["costes_sig"].get()),
            ("Monte Carlo p", self._rv["costes_p"].get()),
        ]

        try:
            with PdfPages(out_path) as pdf:
                # First page: human-readable metrics table and run metadata.
                summary_fig = plt.figure(figsize=(8.27, 11.69), facecolor="white")
                ax = summary_fig.add_axes([0.06, 0.06, 0.88, 0.88])
                ax.axis("off")

                lines = [
                    "Colocalization Analyzer - PDF Summary",
                    "",
                    f"App version: {APP_VERSION}",
                    f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
                    f"Image shape: {self.ch1.shape}",
                    f"ROI: {roi_text}",
                    f"Costes regression: {'Orthogonal (TLS)' if self._orthreg_var.get() else 'Ordinary least squares'}",
                    f"Background subtraction: Ch1={self._bg1_var.get()}  Ch2={self._bg2_var.get()}",
                    "",
                ]

                if has_analysis:
                    lines.append("Results")
                    lines.append("-" * 72)
                    for label, value in metrics:
                        lines.append(f"{label:<28} {value}")
                else:
                    lines.append("No analysis results available yet.")
                    lines.append("Run 'Costes + Analyze' first to populate metrics.")

                ax.text(
                    0.0, 1.0,
                    "\n".join(lines),
                    va="top", ha="left",
                    fontsize=10,
                    family="monospace",
                    color="black",
                )
                pdf.savefig(summary_fig)
                plt.close(summary_fig)

                # Following pages: current UI figures exactly as shown.
                pdf.savefig(self._img_fig)
                pdf.savefig(self._hist_fig)
                pdf.savefig(self._costes_fig)

            self._status.set(f"PDF summary saved: {Path(out_path).name}")
            messagebox.showinfo("Export complete", f"Saved PDF summary:\n{out_path}")
        except Exception as exc:
            messagebox.showerror("Export error", str(exc))


# ══════════════════════════════════════════════════════════════════════════════
def main():
    root = tk.Tk()
    try:
        root.tk.call("tk", "scaling", 1.5)
    except tk.TclError:
        pass
    ColocalizationApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()
