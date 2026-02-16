#!/usr/bin/env python3
"""
run_rmg_analysis.py
-------------------
Command-line equivalent of the matRockmag MATLAB GUI (matRockmag.m).

The original GUI provided:
  • A file browser to load one or more .rmg files from a directory
  • A loaded-samples list (accumulate multiple files)
  • A routine selector (IRM, dIRM, ARM, LowrieFuller, AF, Fuller, RRM,
      Backfield, StatBox – any combination)
  • Options: Multisample overlay, Subplots grid, AutosaveEPS/PNG,
      custom file prefix
  • A "Run" button  →  calls RmgBatchPlotter
  • A stats checkbox  →  calls RmgStatsWriteTable

This script replicates all of that interactively at the command line.

Usage
-----
    python run_rmg_analysis.py                  # interactive menu
    python run_rmg_analysis.py sample.rmg       # quick full analysis
    python run_rmg_analysis.py a.rmg b.rmg      # load multiple, then menu

Additional modules for scripting
---------------------------------
    from rmg_batch_plotter import rmg_batch_plotter, ALL_ROUTINES
    from rmg_write_tables  import rmg_stats_write_table
    from rmg_curve_fits    import rmg_sirm_derivative_curve_fits_sgg
    from fit_sgg           import fit_sgg, fit_sgg_multi
    from rmg_simulate      import rmg_irm_simulate, codica_error_smoothing
"""

import sys
import os
import textwrap

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from rmg_import        import rmg_import
from rmg_batch_plotter import rmg_batch_plotter, ALL_ROUTINES
from rmg_plots         import rmg_data_full_analysis
from rmg_write_tables  import rmg_stats_write_table

try:
    import matplotlib
    matplotlib.use('TkAgg')
    import matplotlib.pyplot as plt
    MATPLOTLIB_OK = True
except Exception:
    MATPLOTLIB_OK = False


BANNER = textwrap.dedent("""
    ╔══════════════════════════════════════════════════════╗
    ║          matRockmag  –  Python Port                  ║
    ║  Rock-magnetic analysis of .rmg data files           ║
    ╚══════════════════════════════════════════════════════╝
""")


# ─────────────────────────────────────────────────────────────────────────────
#  File loading helpers  (mirrors btnLoadPath + btnLoad in matRockmag.m)
# ─────────────────────────────────────────────────────────────────────────────

def list_rmg_files(directory):
    """Recursively find all .rmg files under directory. Returns absolute paths."""
    found = []
    try:
        for root, dirs, files in os.walk(directory):
            dirs.sort()
            for f in sorted(files):
                if f.lower().endswith('.rmg'):
                    found.append(os.path.join(root, f))
    except Exception:
        pass
    return found


def browse_and_load(loaded):
    """
    Interactive file browser.
    Recursively finds .rmg files under a directory tree; user picks which to load.
    Appends to the `loaded` list.  Mirrors btnLoadPath + btnLoad.
    """
    print("\n── Load files ──────────────────────────────────────────────")
    raw_dir = input("  Enter directory path (or press Enter for current dir): ").strip()
    directory = os.path.expandvars(os.path.expanduser(raw_dir)) if raw_dir else os.getcwd()

    paths = list_rmg_files(directory)
    if not paths:
        print(f"  No .rmg files found under: {directory}")
        return

    print(f"\n  Found {len(paths)} .rmg file(s) under {directory}:")
    for i, p in enumerate(paths):
        rel = os.path.relpath(p, directory)
        print(f"    [{i+1:3d}] {rel}")

    raw = input("\n  Enter numbers to load (e.g. 1 3 4), or 'all': ").strip()
    if raw.lower() == 'all':
        indices = list(range(len(paths)))
    else:
        indices = []
        for tok in raw.split():
            try:
                n = int(tok) - 1
                if 0 <= n < len(paths):
                    indices.append(n)
            except ValueError:
                pass

    for idx in indices:
        path = paths[idx]
        try:
            data = rmg_import(path)
            loaded.append(data)
            rel = os.path.relpath(path, directory)
            mass_str = f', mass={data.mass*1e3:.4g} g' if data.mass == data.mass else ''
            print(f"  ✓ Loaded: {data.samplename}  ({len(data.steptypes)} steps{mass_str})  [{rel}]")
        except Exception as exc:
            print(f"  ✗ Error loading {os.path.basename(path)}: {exc}")


def load_direct(paths, loaded):
    """Load files given as command-line arguments."""
    for p in paths:
        path = os.path.expandvars(os.path.expanduser(p))
        if not os.path.isfile(path):
            print(f"  ✗ File not found: {path}")
            continue
        try:
            data = rmg_import(path)
            loaded.append(data)
            print(f"  ✓ Loaded: {data.samplename}  ({len(data.steptypes)} steps)")
        except Exception as exc:
            print(f"  ✗ Error: {exc}")


# ─────────────────────────────────────────────────────────────────────────────
#  Routine / options selection  (mirrors listRoutine + checkboxes in GUI)
# ─────────────────────────────────────────────────────────────────────────────

def select_routines():
    """Prompt user to choose which plot routines to run."""
    print("\n── Select routines ─────────────────────────────────────────")
    for i, r in enumerate(ALL_ROUTINES):
        print(f"    [{i+1}] {r}")
    raw = input("  Enter numbers (e.g. 1 2 3), or 'all' [default: all]: ").strip()
    if not raw or raw.lower() == 'all':
        return list(ALL_ROUTINES)
    selected = []
    for tok in raw.split():
        try:
            n = int(tok) - 1
            if 0 <= n < len(ALL_ROUTINES):
                selected.append(ALL_ROUTINES[n])
        except ValueError:
            pass
    return selected if selected else list(ALL_ROUTINES)


def select_options():
    """Prompt for the four option flags from the GUI."""
    print("\n── Options ─────────────────────────────────────────────────")

    def yn(prompt, default='n'):
        raw = input(f"  {prompt} [y/N]: ").strip().lower()
        return raw in ('y', 'yes')

    multisample   = yn("Multisample (overlay all loaded samples on each plot)?")
    subplots      = yn("Subplots (show all routines in one figure)?")
    autosave_png  = yn("Autosave PNG?")
    autosave_eps  = yn("Autosave EPS?")
    stats_table   = yn("Export statistics table (.asc)?")

    raw_prefix = input("  File prefix (e.g. 'run1-', or Enter for none): ").strip()

    formats = []
    if autosave_png:
        formats.append('png')
    if autosave_eps:
        formats.append('eps')

    return {
        'multisample':      multisample,
        'subplots':         subplots,
        'autosave_formats': formats,
        'stats_table':      stats_table,
        'file_prefix':      raw_prefix,
    }


# ─────────────────────────────────────────────────────────────────────────────
#  Main menu loop  (mirrors the GUI window)
# ─────────────────────────────────────────────────────────────────────────────

def show_loaded(loaded):
    if not loaded:
        print("  (no samples loaded)")
    else:
        for i, d in enumerate(loaded):
            print(f"    [{i+1}] {d.samplename}")


def main():
    print(BANNER)

    loaded = []   # global RmgData equivalent

    # ── if files were passed on the command line, load them straight away ──
    if len(sys.argv) > 1:
        print("Loading command-line files…")
        load_direct(sys.argv[1:], loaded)

    # ── interactive menu loop ──────────────────────────────────────────────
    while True:
        print("\n── Loaded samples ──────────────────────────────────────────")
        show_loaded(loaded)

        print("\n── Actions ─────────────────────────────────────────────────")
        print("    [L] Load more files from a directory")
        print("    [R] Run selected routines (batch plotter)")
        print("    [F] Full 3×3 analysis dashboard (all routines, subplots)")
        print("    [S] Export statistics table only")
        print("    [C] Clear all loaded samples")
        print("    [Q] Quit")

        choice = input("\n  Choice: ").strip().upper()

        if choice == 'L':
            browse_and_load(loaded)

        elif choice == 'R':
            if not loaded:
                print("  No samples loaded.")
                continue
            routines = select_routines()
            opts     = select_options()

            # optionally let user pick a subset of loaded samples
            if len(loaded) > 1:
                print("\n── Select samples to plot ───────────────────────────────")
                show_loaded(loaded)
                raw = input("  Enter numbers (or 'all') [default: all]: ").strip()
                if raw and raw.lower() != 'all':
                    sel = []
                    for tok in raw.split():
                        try:
                            n = int(tok) - 1
                            if 0 <= n < len(loaded):
                                sel.append(loaded[n])
                        except ValueError:
                            pass
                    data_to_plot = sel if sel else loaded
                else:
                    data_to_plot = loaded
            else:
                data_to_plot = loaded

            print("\nRunning batch plotter…")
            figs = rmg_batch_plotter(
                data_to_plot,
                routines,
                multisample      = opts['multisample'],
                subplots         = opts['subplots'],
                autosave_formats = opts['autosave_formats'],
                file_prefix      = opts['file_prefix'],
            )

            if opts['stats_table']:
                stem = opts['file_prefix'] + 'rockmagstats-summary'
                rmg_stats_write_table(data_to_plot, stem)

            if figs and MATPLOTLIB_OK:
                print(f"\n{len(figs)} figure(s) ready. Showing (close all windows to continue)…")
                plt.show()

        elif choice == 'F':
            if not loaded:
                print("  No samples loaded.")
                continue
            print("\nGenerating full 3×3 dashboard…")
            fig = rmg_data_full_analysis(loaded)
            raw = input("Save figure as PNG? [y/N] ").strip().lower()
            if raw in ('y', 'yes'):
                name = (loaded[0].samplename if len(loaded) == 1 else 'multisamples')
                fig.savefig(f'{name}_rockmag.png', dpi=150, bbox_inches='tight')
                print(f"  Saved → {name}_rockmag.png")
            if MATPLOTLIB_OK:
                plt.show()

        elif choice == 'S':
            if not loaded:
                print("  No samples loaded.")
                continue
            raw_prefix = input("  File prefix (or Enter for none): ").strip()
            stem = raw_prefix + 'rockmagstats-summary'
            rmg_stats_write_table(loaded, stem)

        elif choice == 'C':
            loaded.clear()
            print("  All samples cleared.")

        elif choice == 'Q':
            print("Goodbye.")
            break

        else:
            print("  Unrecognised choice.")


if __name__ == '__main__':
    main()
