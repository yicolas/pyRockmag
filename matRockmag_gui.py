"""
matRockmag GUI - Python/Tkinter port of the MATLAB GUIDE matRockmag.m interface
Original: Robert E. Kopp, licensed under GNU GPL v3
Python port: 2026

Requirements: Python 3.x (tkinter is included in standard library)
Usage: python matRockmag_gui.py
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import os
import glob


ROUTINE_TYPES = ["IRM", "dIRM", "ARM", "LowrieFuller", "AF", "Fuller", "RRM", "Backfield", "StatBox"]

# Global equivalent of MATLAB's "global RmgData"
RmgData = []


class MatRockmagGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("matRockmag")
        self.root.resizable(False, False)

        # --- State variables ---
        self.pathname_var    = tk.StringVar(value=os.getcwd() + os.sep)
        self.fileprefix_var  = tk.StringVar()
        self.aflevel_var     = tk.StringVar()
        self.chk_multisamples = tk.BooleanVar()
        self.chk_subplots     = tk.BooleanVar()
        self.chk_autosave_eps = tk.BooleanVar()
        self.chk_autosave_fig = tk.BooleanVar()
        self.chk_rmgstats     = tk.BooleanVar()

        self._build_ui()
        self._refresh_file_list()

    # ------------------------------------------------------------------
    # UI construction
    # ------------------------------------------------------------------

    def _build_ui(self):
        pad = dict(padx=6, pady=4)

        # ── Row 0: Path bar ───────────────────────────────────────────
        path_frame = ttk.LabelFrame(self.root, text="Data Path")
        path_frame.grid(row=0, column=0, columnspan=2, sticky="ew", **pad)

        self.txt_pathname = ttk.Entry(path_frame, textvariable=self.pathname_var, width=55)
        self.txt_pathname.grid(row=0, column=0, padx=4, pady=4, sticky="ew")
        self.txt_pathname.bind("<Return>", lambda e: self._on_pathname_return())

        ttk.Button(path_frame, text="Browse…", command=self._browse_path).grid(
            row=0, column=1, padx=4, pady=4)

        path_frame.columnconfigure(0, weight=1)

        # ── Row 1 left: Available files ───────────────────────────────
        files_frame = ttk.LabelFrame(self.root, text="Available .rmg Files")
        files_frame.grid(row=1, column=0, sticky="nsew", **pad)

        self.list_files = tk.Listbox(files_frame, selectmode=tk.EXTENDED,
                                     width=30, height=10, exportselection=False)
        sb_files = ttk.Scrollbar(files_frame, orient=tk.VERTICAL,
                                 command=self.list_files.yview)
        self.list_files.config(yscrollcommand=sb_files.set)
        self.list_files.grid(row=0, column=0, sticky="nsew", padx=(4, 0), pady=4)
        sb_files.grid(row=0, column=1, sticky="ns", pady=4)

        btn_frame_files = ttk.Frame(files_frame)
        btn_frame_files.grid(row=1, column=0, columnspan=2, sticky="ew", padx=4, pady=(0, 4))
        ttk.Button(btn_frame_files, text="Load Selected", command=self._load_files).pack(side=tk.LEFT, padx=2)

        files_frame.rowconfigure(0, weight=1)
        files_frame.columnconfigure(0, weight=1)

        # ── Row 1 right: Routines ─────────────────────────────────────
        routines_frame = ttk.LabelFrame(self.root, text="Analysis Routines")
        routines_frame.grid(row=1, column=1, sticky="nsew", **pad)

        self.list_routines = tk.Listbox(routines_frame, selectmode=tk.EXTENDED,
                                        width=20, height=10, exportselection=False)
        sb_routines = ttk.Scrollbar(routines_frame, orient=tk.VERTICAL,
                                    command=self.list_routines.yview)
        self.list_routines.config(yscrollcommand=sb_routines.set)
        self.list_routines.grid(row=0, column=0, sticky="nsew", padx=(4, 0), pady=4)
        sb_routines.grid(row=0, column=1, sticky="ns", pady=4)

        for r in ROUTINE_TYPES:
            self.list_routines.insert(tk.END, r)

        routines_frame.rowconfigure(0, weight=1)
        routines_frame.columnconfigure(0, weight=1)

        # ── Row 2 left: Loaded data ───────────────────────────────────
        loaded_frame = ttk.LabelFrame(self.root, text="Loaded Data (RmgData)")
        loaded_frame.grid(row=2, column=0, sticky="nsew", **pad)

        self.list_rmgdata = tk.Listbox(loaded_frame, selectmode=tk.EXTENDED,
                                       width=30, height=8, exportselection=False)
        sb_loaded = ttk.Scrollbar(loaded_frame, orient=tk.VERTICAL,
                                  command=self.list_rmgdata.yview)
        self.list_rmgdata.config(yscrollcommand=sb_loaded.set)
        self.list_rmgdata.grid(row=0, column=0, sticky="nsew", padx=(4, 0), pady=4)
        sb_loaded.grid(row=0, column=1, sticky="ns", pady=4)

        btn_frame_loaded = ttk.Frame(loaded_frame)
        btn_frame_loaded.grid(row=1, column=0, columnspan=2, sticky="ew", padx=4, pady=(0, 4))
        ttk.Button(btn_frame_loaded, text="Clear All", command=self._clear_data).pack(side=tk.LEFT, padx=2)

        loaded_frame.rowconfigure(0, weight=1)
        loaded_frame.columnconfigure(0, weight=1)

        # ── Row 2 right: Options ──────────────────────────────────────
        options_frame = ttk.LabelFrame(self.root, text="Options")
        options_frame.grid(row=2, column=1, sticky="nsew", **pad)

        ttk.Checkbutton(options_frame, text="Multi-samples",
                        variable=self.chk_multisamples).grid(row=0, column=0, sticky="w", padx=6, pady=2)
        ttk.Checkbutton(options_frame, text="Subplots",
                        variable=self.chk_subplots).grid(row=1, column=0, sticky="w", padx=6, pady=2)
        ttk.Checkbutton(options_frame, text="Autosave EPS",
                        variable=self.chk_autosave_eps).grid(row=2, column=0, sticky="w", padx=6, pady=2)
        ttk.Checkbutton(options_frame, text="Autosave FIG",
                        variable=self.chk_autosave_fig).grid(row=3, column=0, sticky="w", padx=6, pady=2)
        ttk.Checkbutton(options_frame, text="Write RmgStats Table",
                        variable=self.chk_rmgstats).grid(row=4, column=0, sticky="w", padx=6, pady=2)

        ttk.Separator(options_frame, orient=tk.HORIZONTAL).grid(
            row=5, column=0, columnspan=2, sticky="ew", padx=6, pady=4)

        ttk.Label(options_frame, text="File Prefix:").grid(row=6, column=0, sticky="w", padx=6)
        ttk.Entry(options_frame, textvariable=self.fileprefix_var, width=18).grid(
            row=7, column=0, sticky="ew", padx=6, pady=2)

        ttk.Label(options_frame, text="AF Level (mT):").grid(row=8, column=0, sticky="w", padx=6)
        ttk.Entry(options_frame, textvariable=self.aflevel_var, width=18).grid(
            row=9, column=0, sticky="ew", padx=6, pady=2)

        options_frame.columnconfigure(0, weight=1)

        # ── Row 3: Run bar ────────────────────────────────────────────
        run_frame = ttk.Frame(self.root)
        run_frame.grid(row=3, column=0, columnspan=2, sticky="ew", **pad)

        self.btn_run = ttk.Button(run_frame, text="Run", command=self._run, width=12)
        self.btn_run.pack(side=tk.LEFT, padx=4)

        self.lbl_busy = ttk.Label(run_frame, text="", foreground="red", font=("TkDefaultFont", 9, "bold"))
        self.lbl_busy.pack(side=tk.LEFT, padx=8)

        # ── Status bar ────────────────────────────────────────────────
        self.status_var = tk.StringVar(value="Ready.")
        ttk.Label(self.root, textvariable=self.status_var, relief=tk.SUNKEN,
                  anchor=tk.W).grid(row=4, column=0, columnspan=2, sticky="ew", padx=6, pady=(0, 4))

        self.root.columnconfigure(0, weight=1)
        self.root.columnconfigure(1, weight=1)

    # ------------------------------------------------------------------
    # Callbacks
    # ------------------------------------------------------------------

    def _refresh_file_list(self):
        """Scan current pathname for .rmg files and populate listFiles."""
        pathname = self.pathname_var.get()
        self.list_files.delete(0, tk.END)
        try:
            files = sorted(glob.glob(os.path.join(pathname, "*.rmg")))
            for f in files:
                self.list_files.insert(tk.END, os.path.basename(f))
            self.status_var.set(f"Found {len(files)} .rmg file(s) in {pathname}")
        except Exception as e:
            self.status_var.set(f"Error scanning path: {e}")

    def _on_pathname_return(self):
        pathname = self.pathname_var.get()
        try:
            os.chdir(pathname)
        except Exception:
            pass
        self._refresh_file_list()

    def _browse_path(self):
        chosen = filedialog.askdirectory(title="Select data folder")
        if chosen:
            path = chosen + os.sep
            self.pathname_var.set(path)
            try:
                os.chdir(path)
            except Exception:
                pass
            self._refresh_file_list()

    def _load_files(self):
        """Load selected files into RmgData (calls RMGImport equivalent)."""
        global RmgData

        sel_indices = self.list_files.curselection()
        if not sel_indices:
            messagebox.showwarning("No selection", "Please select one or more .rmg files to load.")
            return

        self._set_busy(True)
        pathname = self.pathname_var.get()

        for i in sel_indices:
            filename = self.list_files.get(i)
            filepath = os.path.join(pathname, filename)
            samplename = os.path.splitext(filename)[0]

            # Avoid duplicates
            existing_names = [d["samplename"] for d in RmgData]
            if samplename in existing_names:
                continue

            # In a full port this would call RMGImport(filepath).
            # Here we store a placeholder dict matching the expected structure.
            RmgData.append({
                "samplename": samplename,
                "filepath": filepath,
                "loaded": True,
            })

        self._refresh_rmgdata_list()
        self._set_busy(False)
        self.status_var.set(f"Loaded. {len(RmgData)} sample(s) in RmgData.")

    def _clear_data(self):
        global RmgData
        RmgData = []
        self.list_rmgdata.delete(0, tk.END)
        self.status_var.set("RmgData cleared.")

    def _refresh_rmgdata_list(self):
        self.list_rmgdata.delete(0, tk.END)
        for d in RmgData:
            self.list_rmgdata.insert(tk.END, d["samplename"])

    def _run(self):
        """Equivalent of btnRun_Callback — calls RmgBatchPlotter with selected options."""
        global RmgData

        sel_data = list(self.list_rmgdata.curselection())
        sel_routines = list(self.list_routines.curselection())

        if not RmgData:
            messagebox.showwarning("No data", "No data loaded. Please load .rmg files first.")
            return
        if not sel_data:
            messagebox.showwarning("No selection", "Please select one or more loaded datasets.")
            return
        if not sel_routines:
            messagebox.showwarning("No routines", "Please select at least one analysis routine.")
            return

        self._set_busy(True)

        selected_datasets = [RmgData[i] for i in sel_data]
        selected_routine_names = [ROUTINE_TYPES[i] for i in sel_routines]

        # Build options list (mirrors MATLAB cell array construction)
        options = []
        if self.chk_multisamples.get():
            options.append("Multisample")
        if self.chk_subplots.get():
            options.append("Subplots")
        if self.chk_autosave_eps.get():
            options.append("AutosaveEPS")
        if self.chk_autosave_fig.get():
            options.append("AutosaveFIG")

        fileprefix = self.fileprefix_var.get().strip()
        if fileprefix:
            options += ["FilePrefix", fileprefix]

        # ── In a full port, call the Python equivalent of RmgBatchPlotter here ──
        # RmgBatchPlotter(selected_datasets, selected_routine_names, *options)
        print("=" * 60)
        print("RmgBatchPlotter called with:")
        print(f"  Datasets  : {[d['samplename'] for d in selected_datasets]}")
        print(f"  Routines  : {selected_routine_names}")
        print(f"  Options   : {options}")

        if self.chk_rmgstats.get():
            # RmgStatsWriteTable(selected_datasets, fileprefix + "rockmagstats-summary")
            print(f"  RmgStats  : writing '{fileprefix}rockmagstats-summary'")
        print("=" * 60)

        self._set_busy(False)
        self.status_var.set(
            f"Run complete — {len(selected_datasets)} sample(s), "
            f"routines: {', '.join(selected_routine_names)}"
        )

    def _set_busy(self, busy: bool):
        if busy:
            self.lbl_busy.config(text="BUSY…")
            self.btn_run.config(state=tk.DISABLED)
            self.root.update_idletasks()
        else:
            self.lbl_busy.config(text="")
            self.btn_run.config(state=tk.NORMAL)


# ------------------------------------------------------------------
# Entry point
# ------------------------------------------------------------------

def main():
    root = tk.Tk()
    app = MatRockmagGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
