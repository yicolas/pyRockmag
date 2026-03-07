"""
pyrockmag_gui.py
================
Tkinter GUI front-end for the pyrockmag package.

Covers every module in the package:
  rmg_import          – load .rmg files
  mpms_import         – load MPMS files
  rmg_batch_plotter   – batch plot engine (IRM, dIRM, ARM, LowrieFuller,
                        AF, Fuller, RRM, Backfield, StatBox)
  rmg_plots           – full 3×3 dashboard
  rmg_forc            – FORC diagrams + FORCinel workflow
  rmg_sam             – generate .sam header files
  rmg_curve_fits      – IRM/AF derivative Gaussian / SGG curve fits
  rmg_simulate        – IRM simulation, CODICA smoothing
  rmg_stats           – coercivity statistics
  rmg_write_tables    – export stats and IRM-fits tables
  fit_sgg / sgg       – SGG fitting (called via rmg_curve_fits)

Original MATLAB code © 2008 Robert E. Kopp, GNU GPL v3.
Python port 2025, same licence.

Usage
-----
    python pyrockmag_gui.py
"""

import os
import glob
import threading
import tkinter as tk
from tkinter import ttk, filedialog, messagebox


# ---------------------------------------------------------------------------
# Lazy pyrockmag imports — GUI opens even if deps not yet installed
# ---------------------------------------------------------------------------
def _get(module, *names):
    """Import module and return requested callables, or raise ImportError."""
    import importlib
    m = importlib.import_module(module)
    return tuple(getattr(m, n) for n in names)


# macOS requires ALL matplotlib GUI calls (subplots, show, etc.) to happen on
# the main OS thread.  We use the MacOSX backend when available (native), fall
# back to TkAgg.  plt.ion() puts matplotlib in interactive mode so plt.show()
# is non-blocking and figures appear immediately without a separate event loop.
import matplotlib
try:
    matplotlib.use("MacOSX")    # native macOS backend, always on main thread
except Exception:
    matplotlib.use("TkAgg")
import matplotlib.pyplot as plt  # noqa: E402
plt.ion()                        # interactive mode: show() is non-blocking

ROUTINE_TYPES = ["IRM", "dIRM", "ARM", "LowrieFuller", "AF",
                 "Fuller", "RRM", "Backfield", "StatBox"]

FORC_OPTIONS = [
    "Hysteresis (all curves)",
    "FORC diagram – hot_r",
    "FORC diagram – plasma_r",
    "FORC diagram – inferno_r",
    "FORC contour plots (all colormaps)",
    "All plots (hysteresis + scatter + contour)",
    "FORCinel v3 workflow",
]

CURVE_FIT_OPTIONS = [
    "Gaussian (auto-select n)",
    "SGG multi-component",
    "Fixed-component amplitudes",
    "AF derivative Gaussians",
]

SIMULATE_OPTIONS = [
    "IRM simulation",
    "Switching field distribution",
    "Fractional ARM",
    "ARM curve subtraction",
    "CODICA error smoothing",
]


# ===========================================================================
# Main GUI
# ===========================================================================
class PyRockmagGUI:
    def __init__(self, root: tk.Tk) -> None:
        self.root = root
        self.root.title("pyrockmag")
        self.root.resizable(True, True)

        self._rmg_data: list = []      # loaded RmgData objects
        self._mpms_data: list = []     # loaded MPMS objects
        self._file_paths: list = []    # absolute paths parallel to lb_files

        # Tk variables
        self.var_pathname     = tk.StringVar(value=os.getcwd() + os.sep)
        self.var_file_prefix  = tk.StringVar()
        self.var_af_level     = tk.StringVar()
        self.var_multisamples = tk.BooleanVar()
        self.var_subplots     = tk.BooleanVar()
        self.var_autosave_eps = tk.BooleanVar()
        self.var_autosave_fig = tk.BooleanVar()
        self.var_rmgstats     = tk.BooleanVar()
        self.var_irm_fits     = tk.BooleanVar()
        self.var_status       = tk.StringVar(value="Ready.")

        self._build_ui()
        self._refresh_file_list()

    # -----------------------------------------------------------------------
    # UI layout
    # -----------------------------------------------------------------------
    def _build_ui(self) -> None:
        PAD = dict(padx=5, pady=4)

        # ── Notebook (tabs) ───────────────────────────────────────────
        self.nb = ttk.Notebook(self.root)
        self.nb.pack(fill=tk.BOTH, expand=True, padx=6, pady=6)

        self._tab_main   = ttk.Frame(self.nb)
        self._tab_forc   = ttk.Frame(self.nb)
        self._tab_fits   = ttk.Frame(self.nb)
        self._tab_sam    = ttk.Frame(self.nb)
        self._tab_sim    = ttk.Frame(self.nb)
        self._tab_mpms   = ttk.Frame(self.nb)

        self.nb.add(self._tab_main,  text="  Main  ")
        self.nb.add(self._tab_forc,  text="  FORC  ")
        self.nb.add(self._tab_fits,  text="  Curve Fits  ")
        self.nb.add(self._tab_sam,   text="  SAM Files  ")
        self.nb.add(self._tab_sim,   text="  Simulate  ")
        self.nb.add(self._tab_mpms,  text="  MPMS  ")

        self._build_main_tab(PAD)
        self._build_forc_tab(PAD)
        self._build_fits_tab(PAD)
        self._build_sam_tab(PAD)
        self._build_sim_tab(PAD)
        self._build_mpms_tab(PAD)

        # ── Status bar ────────────────────────────────────────────────
        ttk.Label(self.root, textvariable=self.var_status,
                  relief=tk.SUNKEN, anchor=tk.W).pack(
            side=tk.BOTTOM, fill=tk.X, padx=6, pady=(0, 4))

    # -----------------------------------------------------------------------
    # TAB 1 – Main (load + batch plot)
    # -----------------------------------------------------------------------
    def _build_main_tab(self, PAD):
        tab = self._tab_main

        # Path row
        frm_path = ttk.LabelFrame(tab, text="Data Path")
        frm_path.grid(row=0, column=0, columnspan=3, sticky="ew", **PAD)

        self.ent_pathname = ttk.Entry(
            frm_path, textvariable=self.var_pathname, width=55)
        self.ent_pathname.grid(row=0, column=0, padx=4, pady=4, sticky="ew")
        self.ent_pathname.bind("<Return>", lambda _: self._on_pathname_return())
        ttk.Button(frm_path, text="Browse…",
                   command=self._browse_path).grid(row=0, column=1, padx=4)
        frm_path.columnconfigure(0, weight=1)

        # Available files
        frm_files = ttk.LabelFrame(tab, text="Available Files (.rmg / .frc / .cor)")
        frm_files.grid(row=1, column=0, sticky="nsew", **PAD)

        self.lb_files, sb1 = _listbox(frm_files, height=10, width=30)
        self.lb_files.grid(row=0, column=0, sticky="nsew", padx=(4,0), pady=4)
        sb1.grid(row=0, column=1, sticky="ns", pady=4)

        bf = ttk.Frame(frm_files)
        bf.grid(row=1, column=0, columnspan=2, sticky="w", padx=4, pady=(0,4))
        ttk.Button(bf, text="Load Selected",
                   command=self._load_rmg).pack(side=tk.LEFT, padx=2)
        ttk.Button(bf, text="Select All",
                   command=lambda: self.lb_files.select_set(0, tk.END)).pack(
                       side=tk.LEFT, padx=2)
        frm_files.rowconfigure(0, weight=1)
        frm_files.columnconfigure(0, weight=1)

        # Loaded data
        frm_loaded = ttk.LabelFrame(tab, text="Loaded Data")
        frm_loaded.grid(row=2, column=0, sticky="nsew", **PAD)

        self.lb_rmgdata, sb2 = _listbox(frm_loaded, height=8, width=30)
        self.lb_rmgdata.grid(row=0, column=0, sticky="nsew", padx=(4,0), pady=4)
        sb2.grid(row=0, column=1, sticky="ns", pady=4)

        bl = ttk.Frame(frm_loaded)
        bl.grid(row=1, column=0, columnspan=2, sticky="w", padx=4, pady=(0,4))
        ttk.Button(bl, text="Select All",
                   command=lambda: self.lb_rmgdata.select_set(0, tk.END)).pack(
                       side=tk.LEFT, padx=2)
        ttk.Button(bl, text="Clear All",
                   command=self._clear_data).pack(side=tk.LEFT, padx=2)
        frm_loaded.rowconfigure(0, weight=1)
        frm_loaded.columnconfigure(0, weight=1)

        # Routines
        frm_routines = ttk.LabelFrame(tab, text="Analysis Routines")
        frm_routines.grid(row=1, column=1, sticky="nsew", **PAD)

        self.lb_routines, sb3 = _listbox(frm_routines, height=10, width=16)
        self.lb_routines.grid(row=0, column=0, sticky="nsew", padx=(4,0), pady=4)
        sb3.grid(row=0, column=1, sticky="ns", pady=4)
        for r in ROUTINE_TYPES:
            self.lb_routines.insert(tk.END, r)

        br = ttk.Frame(frm_routines)
        br.grid(row=1, column=0, columnspan=2, sticky="w", padx=4, pady=(0,4))
        ttk.Button(br, text="Select All",
                   command=lambda: self.lb_routines.select_set(0, tk.END)).pack(
                       side=tk.LEFT, padx=2)
        ttk.Button(br, text="Clear",
                   command=lambda: self.lb_routines.select_clear(0, tk.END)).pack(
                       side=tk.LEFT, padx=2)
        frm_routines.rowconfigure(0, weight=1)
        frm_routines.columnconfigure(0, weight=1)

        # Options
        frm_opts = ttk.LabelFrame(tab, text="Options")
        frm_opts.grid(row=2, column=1, sticky="nsew", **PAD)

        checks = [
            ("Multi-samples (overlay)",   self.var_multisamples),
            ("Subplots (grid layout)",     self.var_subplots),
            ("Autosave EPS",               self.var_autosave_eps),
            ("Autosave FIG / PNG",         self.var_autosave_fig),
            ("Export RmgStats Table",      self.var_rmgstats),
            ("Export IRM Fits Table",      self.var_irm_fits),
        ]
        for i, (lbl, var) in enumerate(checks):
            ttk.Checkbutton(frm_opts, text=lbl, variable=var).grid(
                row=i, column=0, sticky="w", padx=8, pady=2)

        ttk.Separator(frm_opts, orient=tk.HORIZONTAL).grid(
            row=len(checks), column=0, sticky="ew", padx=6, pady=4)

        r0 = len(checks) + 1
        ttk.Label(frm_opts, text="File Prefix:").grid(
            row=r0, column=0, sticky="w", padx=8)
        ttk.Entry(frm_opts, textvariable=self.var_file_prefix, width=16).grid(
            row=r0+1, column=0, sticky="ew", padx=8, pady=2)
        ttk.Label(frm_opts, text="AF Level (mT):").grid(
            row=r0+2, column=0, sticky="w", padx=8)
        ttk.Entry(frm_opts, textvariable=self.var_af_level, width=16).grid(
            row=r0+3, column=0, sticky="ew", padx=8, pady=(2,6))
        frm_opts.columnconfigure(0, weight=1)

        # Run buttons
        frm_run = ttk.LabelFrame(tab, text="Execute")
        frm_run.grid(row=1, column=2, rowspan=2, sticky="nsew", **PAD)

        self.btn_run = ttk.Button(frm_run, text="▶  Run Batch",
                                   command=self._run_batch, width=16)
        self.btn_run.grid(row=0, column=0, padx=8, pady=6, sticky="ew")

        ttk.Button(frm_run, text="▶  Full Dashboard",
                   command=self._run_dashboard, width=16).grid(
                       row=1, column=0, padx=8, pady=4, sticky="ew")

        ttk.Button(frm_run, text="▶  Inspect Sample",
                   command=self._inspect_sample, width=16).grid(
                       row=2, column=0, padx=8, pady=4, sticky="ew")

        ttk.Button(frm_run, text="▶  Show Coercivities",
                   command=self._show_coercivities, width=16).grid(
                       row=3, column=0, padx=8, pady=4, sticky="ew")

        ttk.Separator(frm_run, orient=tk.HORIZONTAL).grid(
            row=4, column=0, sticky="ew", padx=8, pady=4)

        ttk.Button(frm_run, text="Export Stats Table",
                   command=self._export_stats, width=16).grid(
                       row=5, column=0, padx=8, pady=4, sticky="ew")

        ttk.Button(frm_run, text="Export IRM Fits",
                   command=self._export_irm_fits, width=16).grid(
                       row=6, column=0, padx=8, pady=4, sticky="ew")

        ttk.Separator(frm_run, orient=tk.HORIZONTAL).grid(
            row=7, column=0, sticky="ew", padx=8, pady=4)

        self.lbl_busy = ttk.Label(frm_run, text="",
                                   foreground="red",
                                   font=("TkDefaultFont", 9, "bold"))
        self.lbl_busy.grid(row=8, column=0, padx=8, pady=4)

        frm_run.columnconfigure(0, weight=1)

        tab.columnconfigure(0, weight=2)
        tab.columnconfigure(1, weight=1)
        tab.columnconfigure(2, weight=1)
        tab.rowconfigure(1, weight=1)
        tab.rowconfigure(2, weight=1)

    # -----------------------------------------------------------------------
    # TAB 2 – FORC
    # -----------------------------------------------------------------------
    def _build_forc_tab(self, PAD):
        tab = self._tab_forc

        ttk.Label(tab, text="Select loaded samples, then choose a FORC plot option.",
                  foreground="gray").grid(row=0, column=0, columnspan=2,
                                          sticky="w", padx=8, pady=(8,4))

        # Sample selector (mirrors lb_rmgdata)
        frm_s = ttk.LabelFrame(tab, text="Select Samples")
        frm_s.grid(row=1, column=0, sticky="nsew", **PAD)
        self.lb_forc_data, sb = _listbox(frm_s, height=12, width=28)
        self.lb_forc_data.grid(row=0, column=0, sticky="nsew", padx=(4,0), pady=4)
        sb.grid(row=0, column=1, sticky="ns", pady=4)
        frm_s.rowconfigure(0, weight=1); frm_s.columnconfigure(0, weight=1)

        # Options
        frm_o = ttk.LabelFrame(tab, text="Plot Options")
        frm_o.grid(row=1, column=1, sticky="nsew", **PAD)

        self.lb_forc_opts, sb2 = _listbox(frm_o, height=8, width=30)
        self.lb_forc_opts.grid(row=0, column=0, sticky="nsew", padx=(4,0), pady=4)
        sb2.grid(row=0, column=1, sticky="ns", pady=4)
        for opt in FORC_OPTIONS:
            self.lb_forc_opts.insert(tk.END, opt)
        self.lb_forc_opts.selection_set(5)  # default: All plots
        frm_o.rowconfigure(0, weight=1); frm_o.columnconfigure(0, weight=1)

        # Smoothing factor
        frm_sf = ttk.Frame(frm_o)
        frm_sf.grid(row=1, column=0, columnspan=2, sticky="ew", padx=4, pady=4)
        ttk.Label(frm_sf, text="Smoothing factor (FORCinel):").pack(side=tk.LEFT, padx=4)
        self.var_forc_smooth = tk.StringVar(value="3")
        ttk.Entry(frm_sf, textvariable=self.var_forc_smooth, width=5).pack(
            side=tk.LEFT, padx=4)

        ttk.Button(frm_o, text="▶  Plot FORC",
                   command=self._run_forc).grid(
                       row=2, column=0, padx=8, pady=6, sticky="ew")

        # Generate measurement script
        frm_g = ttk.LabelFrame(tab, text="Generate FORC Measurement Script")
        frm_g.grid(row=2, column=0, columnspan=2, sticky="ew", **PAD)

        script_fields = [
            ("Ha max (mT):", "var_forc_hamax", "500"),
            ("Ha min (mT):", "var_forc_hamin", "5"),
            ("Steps:",       "var_forc_steps", "100"),
        ]
        for col, (lbl, var, default) in enumerate(script_fields):
            ttk.Label(frm_g, text=lbl).grid(row=0, column=col*2, padx=6, pady=4)
            setattr(self, var, tk.StringVar(value=default))
            ttk.Entry(frm_g, textvariable=getattr(self, var), width=8).grid(
                row=0, column=col*2+1, padx=4, pady=4)
        ttk.Button(frm_g, text="Generate Script",
                   command=self._generate_forc_script).grid(
                       row=0, column=len(script_fields)*2, padx=8, pady=4)

        tab.columnconfigure(0, weight=1)
        tab.columnconfigure(1, weight=1)
        tab.rowconfigure(1, weight=1)

    # -----------------------------------------------------------------------
    # TAB 3 – Curve Fits
    # -----------------------------------------------------------------------
    def _build_fits_tab(self, PAD):
        tab = self._tab_fits

        ttk.Label(tab,
                  text="Fit Gaussian or SGG components to IRM/AF derivative curves.",
                  foreground="gray").grid(row=0, column=0, columnspan=2,
                                          sticky="w", padx=8, pady=(8,4))

        frm_s = ttk.LabelFrame(tab, text="Select Samples")
        frm_s.grid(row=1, column=0, sticky="nsew", **PAD)
        self.lb_fits_data, sb = _listbox(frm_s, height=10, width=28)
        self.lb_fits_data.grid(row=0, column=0, sticky="nsew", padx=(4,0), pady=4)
        sb.grid(row=0, column=1, sticky="ns", pady=4)
        frm_s.rowconfigure(0, weight=1); frm_s.columnconfigure(0, weight=1)

        frm_o = ttk.LabelFrame(tab, text="Fit Method")
        frm_o.grid(row=1, column=1, sticky="nsew", **PAD)

        self.lb_fit_opts, sb2 = _listbox(frm_o, height=5, width=30)
        self.lb_fit_opts.grid(row=0, column=0, sticky="nsew", padx=(4,0), pady=4)
        sb2.grid(row=0, column=1, sticky="ns", pady=4)
        for opt in CURVE_FIT_OPTIONS:
            self.lb_fit_opts.insert(tk.END, opt)
        self.lb_fit_opts.selection_set(0)
        frm_o.rowconfigure(0, weight=1); frm_o.columnconfigure(0, weight=1)

        # Fixed-component means/stds entry (for option 3)
        frm_fixed = ttk.LabelFrame(tab,
            text="Fixed Component Means & Stds (log10 mT, comma-separated)")
        frm_fixed.grid(row=2, column=0, columnspan=2, sticky="ew", **PAD)

        ttk.Label(frm_fixed, text="Means:").grid(row=0, column=0, padx=6, pady=4)
        self.var_fit_means = tk.StringVar(value="1.0, 1.5, 2.0")
        ttk.Entry(frm_fixed, textvariable=self.var_fit_means, width=30).grid(
            row=0, column=1, padx=4, pady=4, sticky="ew")

        ttk.Label(frm_fixed, text="Stds:").grid(row=1, column=0, padx=6, pady=4)
        self.var_fit_stds = tk.StringVar(value="0.2, 0.2, 0.2")
        ttk.Entry(frm_fixed, textvariable=self.var_fit_stds, width=30).grid(
            row=1, column=1, padx=4, pady=4, sticky="ew")
        frm_fixed.columnconfigure(1, weight=1)

        self.var_fits_plot  = tk.BooleanVar(value=True)
        self.var_fits_save  = tk.BooleanVar(value=False)
        ttk.Checkbutton(frm_o, text="Show plot",
                        variable=self.var_fits_plot).grid(
                            row=1, column=0, sticky="w", padx=8)
        ttk.Checkbutton(frm_o, text="Export fits table",
                        variable=self.var_fits_save).grid(
                            row=2, column=0, sticky="w", padx=8)
        ttk.Button(frm_o, text="▶  Run Fits",
                   command=self._run_curve_fits).grid(
                       row=3, column=0, padx=8, pady=6, sticky="ew")

        tab.columnconfigure(0, weight=1)
        tab.columnconfigure(1, weight=1)
        tab.rowconfigure(1, weight=1)

    # -----------------------------------------------------------------------
    # TAB 4 – SAM Files
    # -----------------------------------------------------------------------
    def _build_sam_tab(self, PAD):
        tab = self._tab_sam

        ttk.Label(tab,
                  text="Generate CIT-format .sam header files for paleomagnetic analysis.",
                  foreground="gray").grid(row=0, column=0, columnspan=2,
                                          sticky="w", padx=8, pady=(8,4))

        frm_mode = ttk.LabelFrame(tab, text="Mode")
        frm_mode.grid(row=1, column=0, sticky="nsew", **PAD)

        self.var_sam_mode = tk.StringVar(value="blank")
        ttk.Radiobutton(frm_mode, text="Generate blank template",
                        variable=self.var_sam_mode,
                        value="blank").grid(row=0, column=0, sticky="w", padx=8, pady=4)
        ttk.Radiobutton(frm_mode, text="Interactive (enter orientations)",
                        variable=self.var_sam_mode,
                        value="interactive").grid(row=1, column=0, sticky="w", padx=8, pady=4)

        frm_fields = ttk.LabelFrame(tab, text="Blank Template Parameters")
        frm_fields.grid(row=2, column=0, sticky="ew", **PAD)

        sam_params = [
            ("Site ID (≤8 chars):",   "var_sam_site",    "SITE01"),
            ("Number of samples:",    "var_sam_nsamples", "10"),
            ("Sample extension:",     "var_sam_ext",      "a"),
            ("Output directory:",     "var_sam_outdir",   "."),
        ]
        for i, (lbl, var, default) in enumerate(sam_params):
            ttk.Label(frm_fields, text=lbl).grid(row=i, column=0, sticky="w", padx=8, pady=3)
            setattr(self, var, tk.StringVar(value=default))
            ttk.Entry(frm_fields, textvariable=getattr(self, var), width=20).grid(
                row=i, column=1, sticky="ew", padx=4, pady=3)
        frm_fields.columnconfigure(1, weight=1)

        frm_opts = ttk.LabelFrame(tab, text="Options")
        frm_opts.grid(row=3, column=0, sticky="ew", **PAD)
        self.var_sam_core_strike = tk.BooleanVar(value=True)
        self.var_sam_validate    = tk.BooleanVar(value=True)
        ttk.Checkbutton(frm_opts, text="Core strike at 90°",
                        variable=self.var_sam_core_strike).grid(
                            row=0, column=0, sticky="w", padx=8, pady=3)
        ttk.Checkbutton(frm_opts, text="Validate 8.3 filenames",
                        variable=self.var_sam_validate).grid(
                            row=1, column=0, sticky="w", padx=8, pady=3)

        ttk.Button(tab, text="▶  Generate SAM Files",
                   command=self._run_sam).grid(
                       row=4, column=0, padx=8, pady=8, sticky="w")

        tab.columnconfigure(0, weight=1)

    # -----------------------------------------------------------------------
    # TAB 5 – Simulate
    # -----------------------------------------------------------------------
    def _build_sim_tab(self, PAD):
        tab = self._tab_sim

        ttk.Label(tab,
                  text="IRM simulation, switching field distribution, CODICA smoothing.",
                  foreground="gray").grid(row=0, column=0, columnspan=2,
                                          sticky="w", padx=8, pady=(8,4))

        frm_s = ttk.LabelFrame(tab, text="Select Samples (for ARM / CODICA)")
        frm_s.grid(row=1, column=0, sticky="nsew", **PAD)
        self.lb_sim_data, sb = _listbox(frm_s, height=10, width=28)
        self.lb_sim_data.grid(row=0, column=0, sticky="nsew", padx=(4,0), pady=4)
        sb.grid(row=0, column=1, sticky="ns", pady=4)
        frm_s.rowconfigure(0, weight=1); frm_s.columnconfigure(0, weight=1)

        frm_o = ttk.LabelFrame(tab, text="Simulation Type")
        frm_o.grid(row=1, column=1, sticky="nsew", **PAD)

        self.lb_sim_opts, sb2 = _listbox(frm_o, height=6, width=30)
        self.lb_sim_opts.grid(row=0, column=0, sticky="nsew", padx=(4,0), pady=4)
        sb2.grid(row=0, column=1, sticky="ns", pady=4)
        for opt in SIMULATE_OPTIONS:
            self.lb_sim_opts.insert(tk.END, opt)
        self.lb_sim_opts.selection_set(0)
        frm_o.rowconfigure(0, weight=1); frm_o.columnconfigure(0, weight=1)

        # IRM simulation parameters
        frm_p = ttk.LabelFrame(tab, text="IRM Simulation Parameters")
        frm_p.grid(row=2, column=0, columnspan=2, sticky="ew", **PAD)

        sim_params = [
            ("Ban values (T, comma-sep):", "var_sim_ban",  "0.01, 0.05, 0.1"),
            ("Linewidth values (T):",       "var_sim_lw",   "0.005, 0.01"),
            ("B range (T), B_max:",         "var_sim_bmax", "0.3"),
            ("CODICA smooth span:",         "var_sim_span", "0"),
        ]
        for i, (lbl, var, default) in enumerate(sim_params):
            ttk.Label(frm_p, text=lbl).grid(row=i, column=0, sticky="w", padx=8, pady=3)
            setattr(self, var, tk.StringVar(value=default))
            ttk.Entry(frm_p, textvariable=getattr(self, var), width=24).grid(
                row=i, column=1, sticky="ew", padx=4, pady=3)
        frm_p.columnconfigure(1, weight=1)

        ttk.Button(tab, text="▶  Run Simulation",
                   command=self._run_simulate).grid(
                       row=3, column=0, padx=8, pady=8, sticky="w")

        tab.columnconfigure(0, weight=1)
        tab.columnconfigure(1, weight=1)
        tab.rowconfigure(1, weight=1)

    # -----------------------------------------------------------------------
    # TAB 6 – MPMS
    # -----------------------------------------------------------------------
    def _build_mpms_tab(self, PAD):
        tab = self._tab_mpms

        ttk.Label(tab,
                  text="Load and plot MPMS (SQUID magnetometer) data files.",
                  foreground="gray").grid(row=0, column=0, columnspan=2,
                                          sticky="w", padx=8, pady=(8,4))

        frm_load = ttk.LabelFrame(tab, text="Load MPMS File")
        frm_load.grid(row=1, column=0, sticky="ew", **PAD)

        self.var_mpms_path = tk.StringVar(value="")
        ttk.Entry(frm_load, textvariable=self.var_mpms_path, width=45).grid(
            row=0, column=0, padx=4, pady=4, sticky="ew")
        ttk.Button(frm_load, text="Browse…",
                   command=self._browse_mpms).grid(row=0, column=1, padx=4)
        ttk.Button(frm_load, text="Load",
                   command=self._load_mpms).grid(row=0, column=2, padx=4)
        frm_load.columnconfigure(0, weight=1)

        frm_loaded = ttk.LabelFrame(tab, text="Loaded MPMS Data")
        frm_loaded.grid(row=2, column=0, sticky="nsew", **PAD)
        self.lb_mpms, sb = _listbox(frm_loaded, height=8, width=40)
        self.lb_mpms.grid(row=0, column=0, sticky="nsew", padx=(4,0), pady=4)
        sb.grid(row=0, column=1, sticky="ns", pady=4)
        frm_loaded.rowconfigure(0, weight=1)
        frm_loaded.columnconfigure(0, weight=1)

        frm_plot = ttk.LabelFrame(tab, text="Plot")
        frm_plot.grid(row=3, column=0, sticky="ew", **PAD)
        ttk.Button(frm_plot, text="▶  Plot Selected MPMS Data",
                   command=self._plot_mpms).pack(padx=8, pady=6)

        tab.columnconfigure(0, weight=1)
        tab.rowconfigure(2, weight=1)

    # -----------------------------------------------------------------------
    # Helpers
    # -----------------------------------------------------------------------
    def _set_busy(self, busy: bool) -> None:
        state = tk.DISABLED if busy else tk.NORMAL
        self.btn_run.config(state=state)
        self.lbl_busy.config(text="BUSY…" if busy else "")
        self.root.update_idletasks()

    def _set_status(self, msg: str) -> None:
        self.var_status.set(msg)
        self.root.update_idletasks()

    def _run_on_main(self, plot_fn):
        """
        Run plot_fn() synchronously on the main thread.
        With plt.ion() active, figures appear immediately without plt.show().
        We call plt.pause(0.05) to flush the event queue so windows render.
        """
        try:
            plot_fn()
            plt.pause(0.05)   # flush GUI events so figure windows appear
        except Exception as exc:
            messagebox.showerror("Plot error", str(exc))
            self._set_status("Plot error — see dialog.")
        finally:
            self._set_busy(False)

    def _run_in_thread(self, compute_fn, plot_fn=None):
        """
        Two-phase execution pattern for matplotlib / macOS thread safety:
          1. compute_fn() — pure data work — runs in a background thread.
          2. plot_fn()    — all matplotlib calls — is scheduled on the main
             thread via root.after(0, ...) after compute completes.

        If plot_fn is None, compute_fn must contain no matplotlib calls.
        """
        def _worker():
            try:
                compute_fn()
                if plot_fn is not None:
                    self.root.after(0, lambda: self._run_on_main(plot_fn))
                else:
                    self.root.after(0, lambda: self._set_busy(False))
            except Exception as exc:
                self.root.after(0, lambda e=exc: (
                    messagebox.showerror("Error", str(e)),
                    self._set_status("Error — see dialog.")
                ))
                self.root.after(0, lambda: self._set_busy(False))

        self._set_busy(True)
        threading.Thread(target=_worker, daemon=True).start()

    def _selected_rmg(self, listbox=None):
        """Return list of RmgData objects selected in listbox (default: lb_rmgdata)."""
        lb = listbox or self.lb_rmgdata
        return [self._rmg_data[i] for i in lb.curselection()]

    def _require_selection(self, listbox=None, minimum=1) -> bool:
        sel = self._selected_rmg(listbox)
        if len(sel) < minimum:
            messagebox.showwarning("No selection",
                                   "Select at least one loaded sample first.")
            return False
        return True

    # -----------------------------------------------------------------------
    # Path / file list
    # -----------------------------------------------------------------------
    def _refresh_file_list(self) -> None:
        """
        Recursively walk the data path and populate the file list with all
        .rmg / .frc / .cor files (including .frc.rmg and .cor.rmg variants).
        The listbox displays paths relative to the root so subfolder structure
        is visible.  Full absolute paths are stored in self._file_paths so
        _load_rmg can use them directly.
        """
        pathname = self.var_pathname.get().rstrip(os.sep)
        self.lb_files.delete(0, tk.END)
        self._file_paths = []          # parallel list of absolute paths
        try:
            extensions = {".rmg", ".frc", ".cor"}
            compound   = {".frc.rmg", ".cor.rmg"}

            seen_abs = set()
            entries  = []             # (rel_display, abs_path)

            for root, dirs, files in os.walk(pathname):
                dirs.sort()           # deterministic traversal order
                for fname in sorted(files):
                    abs_path = os.path.join(root, fname)
                    if abs_path in seen_abs:
                        continue
                    flower = fname.lower()
                    # match compound extensions first, then single
                    if any(flower.endswith(c) for c in compound) or                        any(flower.endswith(e) for e in extensions):
                        seen_abs.add(abs_path)
                        rel = os.path.relpath(abs_path, pathname)
                        entries.append((rel, abs_path))

            entries.sort(key=lambda x: x[0].lower())
            for rel, abs_path in entries:
                self.lb_files.insert(tk.END, rel)
                self._file_paths.append(abs_path)

            self._set_status(
                f"{len(entries)} file(s) (.rmg/.frc/.cor) found under {pathname}"
            )
        except Exception as exc:
            self._set_status(f"Path error: {exc}")

    def _on_pathname_return(self) -> None:
        try:
            os.chdir(self.var_pathname.get())
        except Exception:
            pass
        self._refresh_file_list()

    def _browse_path(self) -> None:
        d = filedialog.askdirectory(title="Select data folder")
        if d:
            self.var_pathname.set(d + os.sep)
            try:
                os.chdir(d)
            except Exception:
                pass
            self._refresh_file_list()

    # -----------------------------------------------------------------------
    # Load / clear data
    # -----------------------------------------------------------------------
    def _load_rmg(self) -> None:
        sel = self.lb_files.curselection()
        if not sel:
            messagebox.showwarning("No selection",
                                   "Select files to load (.rmg, .frc, .cor).")
            return
        try:
            rmg_import, = _get("rmg_import", "rmg_import")
        except ImportError as exc:
            messagebox.showerror("Import error", str(exc))
            return

        self._set_busy(True)
        existing = {d.samplename for d in self._rmg_data}
        n = 0
        for i in sel:
            fname = self.lb_files.get(i)          # relative display path
            # Use stored absolute path when available (recursive scan)
            if hasattr(self, "_file_paths") and i < len(self._file_paths):
                fpath = self._file_paths[i]
            else:
                fpath = os.path.join(self.var_pathname.get(), fname)
            # samplename: strip all known compound + single extensions
            sname = fname
            for ext in (".frc.rmg", ".cor.rmg", ".rmg", ".frc", ".cor"):
                if sname.lower().endswith(ext):
                    sname = sname[: -len(ext)]
                    break
            sname = os.path.basename(sname)       # drop subfolder prefix
            if sname in existing:
                continue
            try:
                data = rmg_import(fpath)
                self._rmg_data.append(data)
                existing.add(sname)
                n += 1
            except Exception as exc:
                messagebox.showwarning("Load error",
                                       f"{fname}:\n{exc}")

        self._refresh_all_sample_lists()
        self._set_busy(False)
        self._set_status(
            f"Loaded {n} sample(s). {len(self._rmg_data)} total.")

    def _clear_data(self) -> None:
        self._rmg_data.clear()
        self._refresh_all_sample_lists()
        self._set_status("Data cleared.")

    def _refresh_all_sample_lists(self) -> None:
        """Keep all per-tab sample listboxes in sync with _rmg_data."""
        names = [d.samplename for d in self._rmg_data]
        for lb in (self.lb_rmgdata, self.lb_forc_data,
                   self.lb_fits_data, self.lb_sim_data):
            cur = list(lb.curselection())
            lb.delete(0, tk.END)
            for n in names:
                lb.insert(tk.END, n)
            for i in cur:
                if i < len(names):
                    lb.selection_set(i)

    # -----------------------------------------------------------------------
    # TAB 1 actions
    # -----------------------------------------------------------------------
    def _run_batch(self) -> None:
        datasets  = self._selected_rmg()
        routines  = [ROUTINE_TYPES[i] for i in self.lb_routines.curselection()]
        if not datasets:
            messagebox.showwarning("No data", "Select loaded samples first.")
            return
        if not routines:
            messagebox.showwarning("No routines", "Select at least one routine.")
            return
        try:
            rmg_batch_plotter, = _get("rmg_batch_plotter", "rmg_batch_plotter")
            rmg_stats_write_table, = _get("rmg_write_tables", "rmg_stats_write_table")
            rmg_irm_fits_write_table, = _get("rmg_write_tables", "rmg_irm_fits_write_table")
        except ImportError as exc:
            messagebox.showerror("Import error", str(exc))
            return

        prefix = self.var_file_prefix.get().strip()
        kwargs = dict(
            multisample  = self.var_multisamples.get(),
            subplots     = self.var_subplots.get(),
            autosave_eps = self.var_autosave_eps.get(),
            autosave_fig = self.var_autosave_fig.get(),
            file_prefix  = prefix,
        )
        af_str = self.var_af_level.get().strip()
        if af_str:
            try:
                kwargs["af_level"] = float(af_str)
            except ValueError:
                pass

        do_stats = self.var_rmgstats.get()
        do_fits  = self.var_irm_fits.get()

        def _plot():
            rmg_batch_plotter(datasets, routines, **kwargs)
            if do_stats:
                rmg_stats_write_table(datasets, prefix + "rockmagstats-summary")
            if do_fits:
                rmg_irm_fits_write_table(datasets, prefix + "irmfits-summary")
            self._set_status(
                f"Done — {len(datasets)} sample(s), {', '.join(routines)}")

        self._set_busy(True)
        self.root.after(0, lambda: self._run_on_main(_plot))

    def _run_dashboard(self) -> None:
        datasets = self._selected_rmg()
        if not datasets:
            messagebox.showwarning("No data", "Select loaded samples first.")
            return
        try:
            rmg_data_full_analysis, = _get("rmg_plots", "rmg_data_full_analysis")
        except ImportError as exc:
            messagebox.showerror("Import error", str(exc))
            return

        def _plot():
            rmg_data_full_analysis(datasets)
            self._set_status("Dashboard complete.")

        self._set_busy(True)
        self.root.after(0, lambda: self._run_on_main(_plot))

    def _inspect_sample(self) -> None:
        datasets = self._selected_rmg()
        if not datasets:
            messagebox.showwarning("No data", "Select one sample to inspect.")
            return
        try:
            rmg_stats, = _get("rmg_stats", "rmg_stats")
        except ImportError as exc:
            messagebox.showerror("Import error", str(exc))
            return

        stats = rmg_stats(datasets[:1])
        if not stats:
            messagebox.showinfo("Stats", "No statistics computed.")
            return
        s = stats[0]
        lines = [f"Sample: {s['sample']}", ""]
        for k, v in s.items():
            if k in ("sample", "units"):
                continue
            unit = s.get("units", {}).get(k, "")
            try:
                lines.append(f"  {k:30s} = {v:.5g}  {unit}")
            except Exception:
                lines.append(f"  {k:30s} = {v}")
        messagebox.showinfo(f"Stats — {s['sample']}", "\n".join(lines))
        self._set_status("Inspect done.")

    def _show_coercivities(self) -> None:
        if not self._rmg_data:
            messagebox.showwarning("No data", "Load samples first.")
            return
        try:
            rmg_stats, = _get("rmg_stats", "rmg_stats")
        except ImportError as exc:
            messagebox.showerror("Import error", str(exc))
            return

        stats = rmg_stats(self._rmg_data)
        lines = [f"{'Sample':<25} {'Hcr':>10} {'MDF_IRM':>10} {'MDF_ARM':>10}"]
        lines.append("-" * 60)
        for s in stats:
            hcr = s.get("Hcr", float("nan"))
            mdf = s.get("MDFofIRM", float("nan"))
            mda = s.get("MDFofARM", float("nan"))
            try:
                lines.append(f"{s['sample']:<25} {hcr:>10.4g} {mdf:>10.4g} {mda:>10.4g}")
            except Exception:
                lines.append(s["sample"])
        messagebox.showinfo("Coercivities (mT)", "\n".join(lines))
        self._set_status("Coercivities shown.")

    def _export_stats(self) -> None:
        datasets = self._selected_rmg()
        if not datasets:
            messagebox.showwarning("No data", "Select samples to export.")
            return
        try:
            rmg_stats_write_table, = _get("rmg_write_tables", "rmg_stats_write_table")
        except ImportError as exc:
            messagebox.showerror("Import error", str(exc))
            return

        prefix = self.var_file_prefix.get().strip()
        stem = prefix + "rockmagstats-summary"

        def _work():
            rmg_stats_write_table(datasets, stem)
            self.root.after(0, lambda: self._set_status(
                f"Stats table written → {stem}.asc"))

        self._run_in_thread(_work)

    def _export_irm_fits(self) -> None:
        datasets = self._selected_rmg()
        if not datasets:
            messagebox.showwarning("No data", "Select samples to export.")
            return
        try:
            rmg_irm_fits_write_table, = _get("rmg_write_tables",
                                              "rmg_irm_fits_write_table")
        except ImportError as exc:
            messagebox.showerror("Import error", str(exc))
            return

        prefix = self.var_file_prefix.get().strip()
        stem = prefix + "irmfits-summary"

        def _work():
            rmg_irm_fits_write_table(datasets, stem)
            self.root.after(0, lambda: self._set_status(
                f"IRM fits table written → {stem}.asc"))

        self._run_in_thread(_work)

    # -----------------------------------------------------------------------
    # TAB 2 – FORC actions
    # -----------------------------------------------------------------------
    def _run_forc(self) -> None:
        sel_idx  = self.lb_forc_data.curselection()
        opt_idx  = self.lb_forc_opts.curselection()
        if not sel_idx:
            messagebox.showwarning("No selection", "Select samples.")
            return
        if not opt_idx:
            messagebox.showwarning("No option", "Select a FORC plot option.")
            return

        datasets = [self._rmg_data[i] for i in sel_idx]
        option   = opt_idx[0] + 1  # 1-indexed to match run_pyrockmag choices

        try:
            smooth = int(self.var_forc_smooth.get().strip() or "3")
        except ValueError:
            smooth = 3

        try:
            (rmg_extract_forc_data,
             rmg_plot_forc_curves,
             rmg_forc_diagram,
             process_forc_forcinel_workflow) = _get(
                "rmg_forc",
                "rmg_extract_forc_data",
                "rmg_plot_forc_curves",
                "rmg_forc_diagram",
                "process_forc_forcinel_workflow",
            )
        except ImportError as exc:
            messagebox.showerror("Import error", str(exc))
            return

        def _plot():
            for d in datasets:
                if option == 1:
                    rmg_plot_forc_curves(d)
                elif option in (2, 3, 4):
                    cmaps = ["hot_r", "plasma_r", "inferno_r"]
                    rmg_forc_diagram(d, cmap=cmaps[option-2],
                                     smoothing_factor=smooth)
                elif option == 5:
                    for c in ["hot_r", "plasma_r", "inferno_r"]:
                        rmg_forc_diagram(d, cmap=c,
                                         smoothing_factor=smooth,
                                         contour=True)
                elif option == 6:
                    rmg_plot_forc_curves(d)
                    for c in ["hot_r", "plasma_r", "inferno_r"]:
                        rmg_forc_diagram(d, cmap=c, smoothing_factor=smooth)
                        rmg_forc_diagram(d, cmap=c, smoothing_factor=smooth,
                                         contour=True)
                elif option == 7:
                    process_forc_forcinel_workflow(d,
                                                   smoothing_factor=smooth)
            self._set_status("FORC plot done.")

        self._set_busy(True)
        self.root.after(0, lambda: self._run_on_main(_plot))

    def _generate_forc_script(self) -> None:
        try:
            generate_forc_script, = _get("rmg_forc", "generate_forc_script")
        except ImportError as exc:
            messagebox.showerror("Import error", str(exc))
            return
        try:
            hamax  = float(self.var_forc_hamax.get())
            hamin  = float(self.var_forc_hamin.get())
            steps  = int(self.var_forc_steps.get())
        except ValueError:
            messagebox.showerror("Input error", "Check Ha max, Ha min, Steps.")
            return

        path = filedialog.asksaveasfilename(
            title="Save FORC script",
            defaultextension=".rmg",
            filetypes=[("RMG script", "*.rmg"), ("All files", "*.*")])
        if not path:
            return

        def _work():
            generate_forc_script(path, ha_max=hamax, ha_min=hamin, steps=steps)
            self.root.after(0, lambda: self._set_status(
                f"FORC script written → {os.path.basename(path)}"))

        self._run_in_thread(_work)

    # -----------------------------------------------------------------------
    # TAB 3 – Curve Fits actions
    # -----------------------------------------------------------------------
    def _run_curve_fits(self) -> None:
        sel_idx = self.lb_fits_data.curselection()
        opt_idx = self.lb_fit_opts.curselection()
        if not sel_idx:
            messagebox.showwarning("No selection", "Select samples.")
            return
        if not opt_idx:
            messagebox.showwarning("No option", "Select a fit method.")
            return

        datasets = [self._rmg_data[i] for i in sel_idx]
        method   = opt_idx[0]

        try:
            (rmg_sirm_derivative_curve_fits,
             rmg_sirm_derivative_curve_fits_sgg,
             rmg_sirm_derivative_curve_fit_comps,
             rmg_af_derivative_curve_fits) = _get(
                "rmg_curve_fits",
                "rmg_sirm_derivative_curve_fits",
                "rmg_sirm_derivative_curve_fits_sgg",
                "rmg_sirm_derivative_curve_fit_comps",
                "rmg_af_derivative_curve_fits",
            )
            rmg_sirm_curve, = _get("rmg_curves", "rmg_sirm_curve")
        except ImportError as exc:
            messagebox.showerror("Import error", str(exc))
            return

        def _parse_list(s):
            return [float(x.strip()) for x in s.split(",") if x.strip()]

        do_plot = self.var_fits_plot.get()
        do_save = self.var_fits_save.get()
        prefix  = self.var_file_prefix.get().strip()

        fit_results = []   # populated by compute, consumed by plot

        def _compute():
            for d in datasets:
                curves = rmg_sirm_curve(d)
                best = next((c for c in curves if c.doesExist), None)
                if best is None:
                    continue
                if method == 0:
                    result = rmg_sirm_derivative_curve_fits(best)
                elif method == 1:
                    result = rmg_sirm_derivative_curve_fits_sgg(best)
                elif method == 2:
                    means = _parse_list(self.var_fit_means.get())
                    stds  = _parse_list(self.var_fit_stds.get())
                    result = rmg_sirm_derivative_curve_fit_comps(best, means, stds)
                elif method == 3:
                    result = rmg_af_derivative_curve_fits(best)
                else:
                    result = None
                fit_results.append((d, best, result))
                if do_save:
                    try:
                        rmg_irm_fits_write_table, = _get(
                            "rmg_write_tables", "rmg_irm_fits_write_table")
                        rmg_irm_fits_write_table([d], prefix + d.samplename + "-fits")
                    except Exception:
                        pass

        def _plot():
            if do_plot:
                for d, best, result in fit_results:
                    fig, ax = plt.subplots()
                    ax.set_title(f"Curve fit — {d.samplename}")
                    ax.set_xlabel("log10 B (mT)")
                    try:
                        xv = best.IRM.logDerivFields
                        yv = best.IRM.logderivSmooth
                        ax.plot(xv, yv, 'k.', label="data")
                        if method in (0, 1) and hasattr(result, "IRM"):
                            bf = result.IRM.logderivbestfit
                            ax.plot(xv, bf["fit"](xv), "r-", label="best fit")
                    except Exception:
                        pass
                    ax.legend()
            self._set_status("Curve fits done.")

        self._run_in_thread(_compute, _plot)

    # -----------------------------------------------------------------------
    # TAB 4 – SAM actions
    # -----------------------------------------------------------------------
    def _run_sam(self) -> None:
        mode = self.var_sam_mode.get()
        try:
            if mode == "blank":
                generate_blank_sam_template, = _get(
                    "rmg_sam", "generate_blank_sam_template")
            else:
                create_sam_interactive, = _get("rmg_sam", "create_sam_interactive")
        except ImportError as exc:
            messagebox.showerror("Import error", str(exc))
            return

        if mode == "blank":
            site_id = self.var_sam_site.get().strip()
            try:
                n_samples = int(self.var_sam_nsamples.get())
            except ValueError:
                n_samples = 10
            out_dir = self.var_sam_outdir.get().strip() or "."

            def _work():
                generate_blank_sam_template(site_id, n_samples, out_dir)
                self.root.after(0, lambda: self._set_status(
                    f"SAM template written for site {site_id} in {out_dir}"))

            self._run_in_thread(_work)
        else:
            # Interactive mode — must run in main thread (stdin required)
            self._set_busy(False)
            try:
                create_sam_interactive()
                self._set_status("Interactive SAM session complete.")
            except Exception as exc:
                messagebox.showerror("SAM error", str(exc))

    # -----------------------------------------------------------------------
    # TAB 5 – Simulate actions
    # -----------------------------------------------------------------------
    def _run_simulate(self) -> None:
        opt_idx = self.lb_sim_opts.curselection()
        if not opt_idx:
            messagebox.showwarning("No option", "Select a simulation type.")
            return
        method = opt_idx[0]

        try:
            (rmg_irm_simulate,
             rmg_simulate_switching_field_dist,
             rmg_fractional_arm,
             rmg_arm_curve_subtract,
             codica_error_smoothing) = _get(
                "rmg_simulate",
                "rmg_irm_simulate",
                "rmg_simulate_switching_field_dist",
                "rmg_fractional_arm",
                "rmg_arm_curve_subtract",
                "codica_error_smoothing",
            )
            import numpy as np
        except ImportError as exc:
            messagebox.showerror("Import error", str(exc))
            return

        def _parse(s):
            return np.array([float(x.strip()) for x in s.split(",") if x.strip()])

        sim_data = {}   # shared between compute and plot closures

        def _compute():
            if method in (0, 1):
                try:
                    ban_vals = _parse(self.var_sim_ban.get())
                    lw_vals  = _parse(self.var_sim_lw.get())
                    bmax     = float(self.var_sim_bmax.get() or "0.3")
                    b_vals   = np.linspace(0, bmax, 200)
                except Exception as exc:
                    self.root.after(0, lambda e=exc: messagebox.showerror(
                        "Parameter error", str(e)))
                    return
                if method == 0:
                    result = rmg_irm_simulate(b_vals, ban_vals, lw_vals)
                else:
                    result = rmg_simulate_switching_field_dist(
                        b_vals, ban_vals, lw_vals)
                sim_data["type"] = "sim"
                sim_data["b_vals"] = b_vals
                sim_data["result"] = result

            elif method in (2, 3):
                sel_datasets = [self._rmg_data[i]
                                 for i in self.lb_sim_data.curselection()]
                sim_data["type"] = "arm"
                sim_data["items"] = []
                for d in sel_datasets:
                    if method == 2:
                        fa = rmg_fractional_arm(d)
                        if fa.get("doesExist"):
                            sim_data["items"].append((d.samplename, fa))

            elif method == 4:
                sel_datasets = [self._rmg_data[i]
                                 for i in self.lb_sim_data.curselection()]
                from rmg_curves import rmg_sirm_curve
                try:
                    span = int(self.var_sim_span.get() or "0")
                except ValueError:
                    span = 0
                sim_data["type"] = "codica"
                sim_data["items"] = []
                for d in sel_datasets:
                    curves = rmg_sirm_curve(d)
                    best = next((c for c in curves if c.doesExist), None)
                    if best is None:
                        continue
                    x = best.IRM.treatmentDCFields
                    y = best.IRM.fracmags
                    res = codica_error_smoothing(x, y, smooth_span=span)
                    sim_data["items"].append((d.samplename, x, y, res))

        def _plot():
            t = sim_data.get("type")
            if t == "sim":
                b_vals = sim_data["b_vals"]
                result = sim_data["result"]
                fig, ax = plt.subplots()
                ax.plot(b_vals * 1000, result.sum(axis=(1, 2)))
                ax.set_xlabel("B (mT)")
                ax.set_ylabel("Absorption (arb)")
                ax.set_title(SIMULATE_OPTIONS[method])
            elif t == "arm":
                for name, fa in sim_data.get("items", []):
                    fig, ax = plt.subplots()
                    ax.plot(fa["treatmentDCField"] * 1000, fa["fracmags"])
                    ax.set_xlabel("B_DC (mT)")
                    ax.set_ylabel("Fractional ARM")
                    ax.set_title(f"Fractional ARM — {name}")
            elif t == "codica":
                for name, x, y, res in sim_data.get("items", []):
                    fig, ax = plt.subplots()
                    ax.plot(x * 1000, y, 'k.', label="original")
                    ax.plot(x * 1000, res["rescaled"], 'r-', label="CODICA smooth")
                    ax.set_xscale("log")
                    ax.set_xlabel("B (mT)")
                    ax.set_ylabel("fIRM")
                    ax.set_title(f"CODICA — {name}")
                    ax.legend()
            self._set_status("Simulation done.")

        self._run_in_thread(_compute, _plot)

    # -----------------------------------------------------------------------
    # TAB 6 – MPMS actions
    # -----------------------------------------------------------------------
    def _browse_mpms(self) -> None:
        path = filedialog.askopenfilename(
            title="Select MPMS file",
            filetypes=[("MPMS / dat files", "*.dat *.rso *.sqd *.txt"),
                       ("All files", "*.*")])
        if path:
            self.var_mpms_path.set(path)

    def _load_mpms(self) -> None:
        path = self.var_mpms_path.get().strip()
        if not path:
            messagebox.showwarning("No file", "Enter or browse for an MPMS file.")
            return
        try:
            mpms_import, = _get("mpms_import", "mpms_import")
        except ImportError as exc:
            messagebox.showerror("Import error",
                                 f"mpms_import.py not found:\n{exc}")
            return
        try:
            data = mpms_import(path)
            self._mpms_data.append(data)
            self.lb_mpms.insert(tk.END, os.path.basename(path))
            self._set_status(f"MPMS loaded: {os.path.basename(path)}")
        except Exception as exc:
            messagebox.showerror("Load error", str(exc))

    def _plot_mpms(self) -> None:
        sel = self.lb_mpms.curselection()
        if not sel:
            messagebox.showwarning("No selection", "Select an MPMS dataset.")
            return
        def _plot():
            for i in sel:
                d = self._mpms_data[i]
                fig, ax = plt.subplots()
                try:
                    ax.plot(d.field, d.moment)
                    ax.set_xlabel("Field (Oe)")
                    ax.set_ylabel("Moment (emu)")
                    ax.set_title(str(d))
                except AttributeError:
                    ax.text(0.5, 0.5, str(d), transform=ax.transAxes,
                            ha="center", va="center")
            self._set_status("MPMS plot done.")

        self._set_busy(True)
        self.root.after(0, lambda: self._run_on_main(_plot))


# ===========================================================================
# Widget helpers
# ===========================================================================
def _listbox(parent, height: int, width: int):
    lb = tk.Listbox(parent, selectmode=tk.EXTENDED,
                    height=height, width=width, exportselection=False)
    sb = ttk.Scrollbar(parent, orient=tk.VERTICAL, command=lb.yview)
    lb.config(yscrollcommand=sb.set)
    return lb, sb


# ===========================================================================
# Entry point
# ===========================================================================
def main() -> None:
    root = tk.Tk()
    PyRockmagGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
