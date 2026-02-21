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
from rmg_stats         import rmg_stats
from rmg_forc          import (rmg_extract_forc_data, rmg_extract_forc_data_complete, prompt_forc_axis_limits,
                                rmg_plot_forc_curves, rmg_forc_diagram, 
                                generate_forc_script, convert_rmg_to_forc,
                                plot_forc_diagram_standard, plot_forc_hysteresis_complete,
                                process_forc_forcinel_workflow)
from rmg_sam           import (SAMFile, SampleOrientation, generate_blank_sam_template,
                                create_sam_interactive)

try:
    import matplotlib
    matplotlib.use('TkAgg')
    import matplotlib.pyplot as plt
    MATPLOTLIB_OK = True
except Exception:
    MATPLOTLIB_OK = False


BANNER = """
===============================================================
          ░▒▓███████▓▒░░▒▓█▓▒░░▒▓█▓▒░ 
          ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░ 
          ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░ 
          ░▒▓███████▓▒░ ░▒▓██████▓▒░  
          ░▒▓█▓▒░         ░▒▓█▓▒░     
          ░▒▓█▓▒░         ░▒▓█▓▒░     
          ░▒▓█▓▒░         ░▒▓█▓▒░     

     ░▒▓███████▓▒░ ░▒▓██████▓▒░  ░▒▓██████▓▒░░▒▓█▓▒░░▒▓█▓▒░ 
     ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░ 
     ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░ 
     ░▒▓███████▓▒░ ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓███████▓▒░  
     ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░ 
     ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░ 
     ░▒▓█▓▒░░▒▓█▓▒░░▒▓██████▓▒░ ░▒▓██████▓▒░░▒▓█▓▒░░▒▓█▓▒░ 

     ░▒▓██████████████▓▒░ ░▒▓██████▓▒░ ░▒▓██████▓▒░░▒▓█▓▒░ 
     ░▒▓█▓▒░░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░ 
     ░▒▓█▓▒░░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░             
     ░▒▓█▓▒░░▒▓█▓▒░░▒▓█▓▒░▒▓████████▓▒░▒▓█▓▒▒▓███▓▒░▒▓█▓▒░ 
     ░▒▓█▓▒░░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░       
     ░▒▓█▓▒░░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░ 
     ░▒▓█▓▒░░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░░▒▓██████▓▒░░▒▓█▓▒░ 

                    v1.0 - CIT Format Suite
        
[ANALYSIS] Hysteresis • Coercivity • IRM • ARM • Backfield
[FORC] Diagrams • FORCinel v3.0.8 • Multiple Colormaps
[PALEOMAG] .SAM Generation • Sun Compass • IGRF • 8.3 Format
[UTILITIES] Batch Plotting • Lowrie-Fuller • Statistics Export
===============================================================
"""


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
#  Sample management
# ─────────────────────────────────────────────────────────────────────────────

def remove_samples(loaded):
    """Let user remove samples from the loaded list."""
    if not loaded:
        print("  No samples to remove.")
        return
    
    print("\n── Remove samples ──────────────────────────────────────────")
    show_loaded(loaded)
    raw = input("\n  Enter numbers to remove (e.g. 2 5 7): ").strip()
    
    to_remove = []
    for tok in raw.split():
        try:
            n = int(tok) - 1
            if 0 <= n < len(loaded):
                to_remove.append(n)
        except ValueError:
            pass
    
    if not to_remove:
        print("  No samples removed.")
        return
    
    # Remove in reverse order to preserve indices
    for idx in sorted(to_remove, reverse=True):
        removed = loaded.pop(idx)
        print(f"  ✗ Removed: {removed.samplename}")


def group_samples(loaded):
    """
    Let user organize samples into groups for plotting.
    Returns: list of lists [[group1_samples], [group2_samples], ...]
    """
    print("\n── Group samples for plotting ──────────────────────────────")
    print("  Samples in the same group will be overlaid on the same axes.")
    print("  Each group gets its own figure set.")
    show_loaded(loaded)
    
    groups = []
    while True:
        print(f"\n  Current groups: {len(groups)}")
        for g_idx, grp in enumerate(groups, start=1):
            names = ', '.join(d.samplename for d in grp)
            print(f"    Group {g_idx}: {names}")
        
        print("\n  Options:")
        print("  + [N] Create new group")
        print("  + [D] Done (proceed with current groups)")
        print("  + [A] All together (one group with all samples)")
        print("  + [S] All separate (one sample per group)")
        
        choice = input("  Choice: ").strip().upper()
        
        if choice == 'D':
            break
        elif choice == 'A':
            return [loaded]
        elif choice == 'S':
            return [[d] for d in loaded]
        elif choice == 'N':
            raw = input(f"  Enter sample numbers for this group (e.g. 1 3 4): ").strip()
            indices = []
            for tok in raw.split():
                try:
                    n = int(tok) - 1
                    if 0 <= n < len(loaded):
                        indices.append(n)
                except ValueError:
                    pass
            if indices:
                grp = [loaded[i] for i in indices]
                groups.append(grp)
                names = ', '.join(d.samplename for d in grp)
                print(f"  ✓ Group {len(groups)} created: {names}")
            else:
                print("  No valid samples selected.")
    
    return groups if groups else [loaded]


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

    # NOTE: multisample is now handled via grouping, not a global flag
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
        'subplots':         subplots,
        'autosave_formats': formats,
        'stats_table':      stats_table,
        'file_prefix':      raw_prefix,
    }


# ─────────────────────────────────────────────────────────────────────────────
#  Info / inspection tools
# ─────────────────────────────────────────────────────────────────────────────

def inspect_sample(loaded):
    """Show detailed stats for one sample."""
    if not loaded:
        print("  No samples loaded.")
        return
    
    print("\n── Inspect sample ──────────────────────────────────────────")
    show_loaded(loaded)
    raw = input("\n  Enter sample number to inspect: ").strip()
    try:
        n = int(raw) - 1
        if 0 <= n < len(loaded):
            d = loaded[n]
            stats = rmg_stats([d])[0]
            
            print(f"\n  Sample: {stats['sample']}")
            print(f"  Steps:  {len(d.steptypes)}")
            print(f"  Blocks: {', '.join(d.BlockType)}")
            print()
            
            def show(key, label, unit=''):
                v = stats.get(key)
                if v is not None and v == v:  # not None and not NaN
                    u = stats['units'].get(key, unit)
                    print(f"    {label:30s} = {v:.4g} {u}".strip())
            
            show('susceptibility', 'χ')
            show('sIRM', 'sIRM')
            show('sIRMperkg', 'sIRM/kg')
            show('Hcr', 'Hcr (coercivity of remanence)')
            show('CisowskiR', 'R (Cisowski)')
            show('MDFofIRM', 'MDF of IRM')
            show('DPofIRM', 'DP of IRM')
            show('ARMsusceptibilityToIRM', 'k_ARM/IRM')
            show('ARMtoIRMat100uT', 'ARM/IRM at 0.1 mT')
            show('MDFofARM', 'MDF of ARM')
        else:
            print("  Invalid sample number.")
    except ValueError:
        print("  Invalid input.")


def show_coercivities(loaded):
    """Print coercivity values (Hcr, MDF_IRM, MDF_ARM) for all loaded samples."""
    if not loaded:
        print("  No samples loaded.")
        return
    
    print("\n── Coercivity values ───────────────────────────────────────")
    
    # Suppress warnings during stats computation
    import warnings
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=RuntimeWarning)
        stats_list = rmg_stats(loaded)
    
    # Filter out samples with all NaN coercivity values
    valid_stats = []
    for s in stats_list:
        hcr = s.get('Hcr')
        mdf_irm = s.get('MDFofIRM')
        mdf_arm = s.get('MDFofARM')
        # Keep sample if at least one value is not NaN
        if any(v is not None and v == v for v in [hcr, mdf_irm, mdf_arm]):
            valid_stats.append(s)
    
    if not valid_stats:
        print("  No samples with coercivity data.")
        return
    
    # Ask for filtering/sorting options
    print(f"\n  Found {len(valid_stats)} sample(s) with coercivity data.")
    print("  Options:")
    print("  + [A] Show all samples")
    print("  + [R] Filter by Hcr range")
    print("  + [S] Sort by parameter")
    
    choice = input("  Choice [A]: ").strip().upper() or 'A'
    
    filtered_stats = valid_stats
    
    if choice == 'R':
        raw_min = input("  Minimum Hcr (mT) [0]: ").strip()
        raw_max = input("  Maximum Hcr (mT) [10000]: ").strip()
        try:
            hcr_min = float(raw_min) if raw_min else 0.0           # already in mT
            hcr_max = float(raw_max) if raw_max else 10000.0       # already in mT
            filtered_stats = []
            for s in valid_stats:
                hcr = s.get('Hcr')
                if hcr is not None and hcr == hcr:  # not None and not NaN
                    if hcr_min <= hcr <= hcr_max:
                        filtered_stats.append(s)
            print(f"  → {len(filtered_stats)} sample(s) in range [{hcr_min:.0f}, {hcr_max:.0f}] mT")
        except ValueError:
            print("  Invalid range, showing all samples.")
            filtered_stats = valid_stats
    
    elif choice == 'S':
        print("  Sort by:")
        print("  + [1] Hcr")
        print("  + [2] MDF_IRM")
        print("  + [3] MDF_ARM")
        print("  + [4] Sample name")
        sort_choice = input("  Choice [1]: ").strip() or '1'
        
        if sort_choice == '1':
            filtered_stats.sort(key=lambda s: s.get('Hcr') if s.get('Hcr') == s.get('Hcr') else float('inf'))
        elif sort_choice == '2':
            filtered_stats.sort(key=lambda s: s.get('MDFofIRM') if s.get('MDFofIRM') == s.get('MDFofIRM') else float('inf'))
        elif sort_choice == '3':
            filtered_stats.sort(key=lambda s: s.get('MDFofARM') if s.get('MDFofARM') == s.get('MDFofARM') else float('inf'))
        elif sort_choice == '4':
            filtered_stats.sort(key=lambda s: s['sample'])
    
    if not filtered_stats:
        print("  No samples match criteria.")
        return
    
    # Print table header
    print(f"\n  {'Sample':<30s}  {'Hcr (mT)':>12s}  {'MDF_IRM (mT)':>12s}  {'MDF_ARM (mT)':>12s}")
    print("  " + "─" * 72)
    
    def fmt(v):
        if v is None or v != v:  # None or NaN
            return "—".rjust(12)
        return f"{v:12.2f}"
    
    for s in filtered_stats:
        name = s['sample'][:30]
        hcr = s.get('Hcr')
        mdf_irm = s.get('MDFofIRM')
        mdf_arm = s.get('MDFofARM')
        print(f"  {name:<30s}  {fmt(hcr)}  {fmt(mdf_irm)}  {fmt(mdf_arm)}")


def export_all_stats(loaded):
    """Export statistics table for all loaded samples."""
    if not loaded:
        print("  No samples loaded.")
        return
    
    print("\n── Export all statistics ───────────────────────────────────")
    raw_prefix = input("  File prefix (or Enter for 'all_samples'): ").strip()
    stem = raw_prefix if raw_prefix else 'all_samples'
    
    rmg_stats_write_table(loaded, stem)
    print(f"  ✓ Exported {len(loaded)} sample(s) to {stem}.asc")


# ─────────────────────────────────────────────────────────────────────────────
#  Main menu loop  (mirrors the GUI window)
# ─────────────────────────────────────────────────────────────────────────────










def select_scatter_brightness():
    """
    Prompt user to select scatter plot brightness mode.
    
    Returns
    -------
    dict : Settings for point size and alpha
    """
    print("\n▼══════════════════════════════════════════════════════════▼")
    print("           Scatter Plot Brightness")
    print("▲══════════════════════════════════════════════════════════▲")
    print("\nBrightness mode (for scatter plots only):")
    print("  + [1] Standard - Subtle, transparent points")
    print("  + [2] Bright - Larger, more opaque points")
    print("  + [3] Maximum - Largest, fully opaque (vibrant)")
    
    choice = input("\nChoice [2]: ").strip() or "2"
    
    settings_map = {
        '1': {'size': 1, 'alpha': 0.6, 'desc': 'Standard'},
        '2': {'size': 3, 'alpha': 0.9, 'desc': 'Bright'},
        '3': {'size': 5, 'alpha': 1.0, 'desc': 'Maximum'}
    }
    
    settings = settings_map.get(choice, settings_map['2'])
    
    print(f"✓ Brightness: {settings['desc']} (size={settings['size']}, alpha={settings['alpha']})")
    
    return settings

def select_forc_grid_resolution():
    """
    Prompt user to select FORC grid resolution.
    
    Higher resolution = more data points = smoother contours in cropped regions.
    
    Returns
    -------
    int : Grid size (number of points per axis)
    """
    print("\n▼══════════════════════════════════════════════════════════▼")
    print("           FORC Data Resolution")
    print("▲══════════════════════════════════════════════════════════▲")
    print("\nGrid resolution (affects data smoothness when cropped):")
    print("  + [1] Standard (100×100) - Fast, good for full range")
    print("  + [2] Fine (200×200) - Better for moderate cropping")
    print("  + [3] Very Fine (300×300) - Best for heavy cropping")
    print("  + [4] Ultra Fine (500×500) - Maximum detail, slower")
    print("\nNote: Higher resolution = smoother contours but slower processing")
    
    choice = input("\nChoice [2]: ").strip() or "2"
    
    grid_map = {
        '1': 100,
        '2': 200,
        '3': 300,
        '4': 500
    }
    
    grid_size = grid_map.get(choice, 200)
    
    print(f"✓ Grid resolution: {grid_size}×{grid_size} = {grid_size**2:,} points")
    
    return grid_size

def set_plot_dpi():
    """
    Set matplotlib DPI for high-resolution plots.
    
    Returns
    -------
    int : DPI value
    """
    import matplotlib as mpl
    
    print("\n▼══════════════════════════════════════════════════════════▼")
    print("           Plot Resolution")
    print("▲══════════════════════════════════════════════════════════▲")
    print("\nResolution options:")
    print("  + [1] Screen (100 DPI) - Fast, lower quality")
    print("  + [2] Standard (150 DPI) - Good balance")
    print("  + [3] High (300 DPI) - Publication quality")
    print("  + [4] Ultra (600 DPI) - Maximum detail")
    
    choice = input("\nChoice [3]: ").strip() or "3"
    
    dpi_map = {
        '1': 100,
        '2': 150,
        '3': 300,
        '4': 600
    }
    
    dpi = dpi_map.get(choice, 300)
    
    # Set matplotlib default DPI
    mpl.rcParams['figure.dpi'] = dpi
    mpl.rcParams['savefig.dpi'] = dpi
    
    print(f"✓ Resolution set to {dpi} DPI")
    
    return dpi

def select_forc_colormap():
    """
    Prompt user to select a colormap for FORC diagrams.
    
    Returns
    -------
    str : Selected colormap name with _r suffix
    """
    print("\n▼══════════════════════════════════════════════════════════▼")
    print("           Colormap Selection")
    print("▲══════════════════════════════════════════════════════════▲")
    print("\nAvailable colormaps:")
    print("  + [1] viridis   - Green/blue/purple (good for weak signals)")
    print("  + [2] inferno   - Yellow/orange/purple (high contrast)")
    print("  + [3] rainbow   - Full spectrum (colorful)")
    print("  + [4] spring    - Magenta/yellow (vibrant)")
    print("  + [5] gray      - Grayscale (publication B&W)")
    print("  + [6] hot       - White/red/black (classic)")
    
    choice = input("\nChoice [2]: ").strip() or "2"
    
    colormap_map = {
        '1': 'viridis_r',
        '2': 'inferno_r',
        '3': 'rainbow',
        '4': 'spring',
        '5': 'gray_r',
        '6': 'hot_r'
    }
    
    colormap = colormap_map.get(choice, 'inferno_r')
    colormap_name = {
        'viridis_r': 'viridis',
        'inferno_r': 'inferno',
        'rainbow': 'rainbow',
        'spring': 'spring',
        'gray_r': 'gray',
        'hot_r': 'hot'
    }.get(colormap, 'inferno')
    
    print(f"✓ Selected: {colormap_name}")
    
    return colormap

def plot_forc_data(loaded):
    """Plot FORC curves and diagram for selected sample(s)."""
    if not loaded:
        print("  No samples loaded.")
        return
    
    print("\n── FORC Analysis ───────────────────────────────────────────")
    show_loaded(loaded)
    
    raw = input("\n  Enter sample number(s) to plot (e.g. 1 3 4): ").strip()
    indices = []
    for tok in raw.split():
        try:
            n = int(tok) - 1
            if 0 <= n < len(loaded):
                indices.append(n)
        except ValueError:
            pass
    
    if not indices:
        print("  No valid samples selected.")
        return
    
    data_list = [loaded[i] for i in indices]
    
    # Check which samples have FORC data
    forc_samples = []
    for d in data_list:
        forc = rmg_extract_forc_data(d)
        if forc.exists:
            forc_samples.append(d)
            print(f"  ✓ {d.samplename}: {len(forc.curves)} FORC curves")
        else:
            print(f"  ✗ {d.samplename}: No FORC data")
    
    if not forc_samples:
        print("\n  No FORC data found in selected samples.")
        return
    
    # Plot options
    print("\n  FORC Plot Options:")
    print("  + [1] Complete hysteresis (all curves with reversals)")
    print("  + [2] FORC diagram - hot_r colormap (classic red/orange)")
    print("  + [3] FORC diagram - plasma_r colormap (yellow/pink/purple)")
    print("  + [4] FORC diagram - inferno_r colormap (high contrast)")
    print("  + [5] FORC diagrams - contour plots (hot_r + plasma_r + inferno_r)")
    print("  + [6] All plots (hysteresis + 3 scatter + 3 contour = 7 plots)")
    print("  + [7] FORCinel v3 workflow (export .frc + processed)")
    
    choice = input("  Choice [6]: ").strip() or '6'
    
    # Extract FORC data with complete reversal points
    forc_complete_list = []
    for d in forc_samples:
        forc_complete = rmg_extract_forc_data_complete(d)
        if forc_complete.exists:
            forc_complete_list.append((d, forc_complete))
    
    if not forc_complete_list:
        print("\n  No FORC data with reversal points found.")
        return
    
    # Option 1: Complete hysteresis
    if choice == '1':
        print("\n  Plotting complete FORC hysteresis...")
        for d, forc_complete in forc_complete_list:
            fig, ax = plt.subplots(figsize=(12, 7))
            plot_forc_hysteresis_complete(forc_complete, ax=ax)
            ax.set_title(f'FORC Hysteresis: {d.samplename}')
            plt.tight_layout()
    
    # Option 2: Hot colormap
    elif choice == '2':
        print("\n  Plotting FORC scatter diagrams...")
        
        # Set grid resolution
        grid_size = select_forc_grid_resolution()
        
        # Set plot resolution
        dpi = set_plot_dpi()
        
        # Select colormap
        colormap = select_forc_colormap()
        
        for d, forc_complete in forc_complete_list:
            print(f"\n  Processing {d.samplename}...")
            
            result = process_forc_forcinel_workflow(d, smoothing_factor=3, grid_size=grid_size)
            
            # Prompt for axis limits
            limits = prompt_forc_axis_limits(result)
            bc_lim = limits['bc_lim'] if limits else None
            bu_lim = limits['bu_lim'] if limits else None
            
            # Generate plot
            fig, ax = plt.subplots(figsize=(12, 10))
            plot_forc_diagram_standard(result, ax=ax, show_points=True, colormap=colormap,
                                      bc_lim=bc_lim, bu_lim=bu_lim)
            ax.set_title(f'{d.samplename} - Scatter')
            plt.tight_layout()
    
    # Option 3: Scatter diagram with user choice
    elif choice == '3':
        print("\n  Plotting FORC scatter diagrams...")
        
        # Set grid resolution
        grid_size = select_forc_grid_resolution()
        
        # Set plot resolution
        dpi = set_plot_dpi()
        
        # Select colormap
        colormap = select_forc_colormap()
        
        for d, forc_complete in forc_complete_list:
            print(f"\n  Processing {d.samplename}...")
            
            result = process_forc_forcinel_workflow(d, smoothing_factor=3, grid_size=grid_size)
            
            # Prompt for axis limits
            limits = prompt_forc_axis_limits(result)
            bc_lim = limits['bc_lim'] if limits else None
            bu_lim = limits['bu_lim'] if limits else None
            
            # Generate plot
            fig, ax = plt.subplots(figsize=(12, 10))
            plot_forc_diagram_standard(result, ax=ax, show_points=True, colormap=colormap,
                                      bc_lim=bc_lim, bu_lim=bu_lim)
            ax.set_title(f'{d.samplename} - Scatter')
            plt.tight_layout()
    
    # Option 4: Scatter diagram with user choice
    elif choice == '4':
        print("\n  Plotting FORC scatter diagrams...")
        
        # Set grid resolution
        grid_size = select_forc_grid_resolution()
        
        # Set plot resolution
        dpi = set_plot_dpi()
        
        # Select colormap
        colormap = select_forc_colormap()
        
        for d, forc_complete in forc_complete_list:
            print(f"\n  Processing {d.samplename}...")
            
            result = process_forc_forcinel_workflow(d, smoothing_factor=3, grid_size=grid_size)
            
            # Prompt for axis limits
            limits = prompt_forc_axis_limits(result)
            bc_lim = limits['bc_lim'] if limits else None
            bu_lim = limits['bu_lim'] if limits else None
            
            # Generate plot
            fig, ax = plt.subplots(figsize=(12, 10))
            plot_forc_diagram_standard(result, ax=ax, show_points=True, colormap=colormap,
                                      bc_lim=bc_lim, bu_lim=bu_lim)
            ax.set_title(f'{d.samplename} - Scatter')
            plt.tight_layout()
    
    # Option 5: Contour plots (all three colormaps)
    elif choice == '5':
        print("\n  Plotting FORC contour diagrams...")
        
        # Set grid resolution
        grid_size = select_forc_grid_resolution()
        
        # Set plot resolution
        dpi = set_plot_dpi()
        
        # Select colormap
        colormap = select_forc_colormap()
        
        for d, forc_complete in forc_complete_list:
            print(f"\n  Processing {d.samplename}...")
            
            # Process FORC data first
            result = process_forc_forcinel_workflow(d, smoothing_factor=3, grid_size=grid_size)
            
            # Prompt for axis limits BEFORE plotting
            limits = prompt_forc_axis_limits(result)
            
            # Extract limits if provided
            bc_lim = limits['bc_lim'] if limits else None
            bu_lim = limits['bu_lim'] if limits else None
            
            # Generate single contour plot with chosen colormap
            fig, ax = plt.subplots(figsize=(12, 10))
            
            plot_forc_diagram_standard(result, ax=ax, show_points=False, colormap=colormap,
                                      bc_lim=bc_lim, bu_lim=bu_lim)
            ax.set_title(f'{d.samplename} - Contour')
            
            plt.tight_layout()
    
    elif choice == '6':
        print("\n  Creating all FORC plots...")
        
        # Set grid resolution
        grid_size = select_forc_grid_resolution()
        
        # Set plot resolution
        dpi = set_plot_dpi()
        
        for d, forc_complete in forc_complete_list:
            print(f"\n  Processing {d.samplename}...")
            
            # Hysteresis
            print("    - Complete hysteresis")
            fig1, ax1 = plt.subplots(figsize=(12, 7))
            plot_forc_hysteresis_complete(forc_complete, ax=ax1)
            ax1.set_title(f'FORC Hysteresis: {d.samplename}')
            plt.tight_layout()
            
            # Process once for all diagrams
            result = process_forc_forcinel_workflow(d, smoothing_factor=3, grid_size=grid_size)
            
            # Prompt for axis limits
            limits = prompt_forc_axis_limits(result)
            bc_lim = limits['bc_lim'] if limits else None
            bu_lim = limits['bu_lim'] if limits else None
            
            # Hot colormap
            print("    - FORC diagram (hot)")
            fig2, ax2 = plt.subplots(figsize=(8, 9))
            plot_forc_diagram_standard(result, ax=ax2, show_points=True, colormap='hot_r', bc_lim=bc_lim, bu_lim=bu_lim)
            ax2.set_title(f'FORC (hot): {d.samplename}')
            plt.tight_layout()
            
            # Plasma colormap
            print("    - FORC diagram (plasma)")
            fig3, ax3 = plt.subplots(figsize=(8, 9))
            plot_forc_diagram_standard(result, ax=ax3, show_points=True, colormap='plasma_r', bc_lim=bc_lim, bu_lim=bu_lim)
            ax3.set_title(f'FORC (plasma): {d.samplename}')
            plt.tight_layout()
            
            # Inferno colormap
            print("    - FORC diagram (inferno)")
            fig4, ax4 = plt.subplots(figsize=(8, 9))
            plot_forc_diagram_standard(result, ax=ax4, show_points=True, colormap='inferno_r', bc_lim=bc_lim, bu_lim=bu_lim)
            ax4.set_title(f'FORC (inferno): {d.samplename}')
            plt.tight_layout()
            
            # Contour plots (all three colormaps)
            print("    - FORC diagram (contour hot_r)")
            fig5, ax5 = plt.subplots(figsize=(8, 9))
            plot_forc_diagram_standard(result, ax=ax5, show_points=False, colormap='hot_r', bc_lim=bc_lim, bu_lim=bu_lim)
            ax5.set_title(f'FORC Contour (hot_r): {d.samplename}')
            plt.tight_layout()
            
            print("    - FORC diagram (contour plasma_r)")
            fig6, ax6 = plt.subplots(figsize=(8, 9))
            plot_forc_diagram_standard(result, ax=ax6, show_points=False, colormap='plasma_r', bc_lim=bc_lim, bu_lim=bu_lim)
            ax6.set_title(f'FORC Contour (plasma_r): {d.samplename}')
            plt.tight_layout()
            
            print("    - FORC diagram (contour inferno_r)")
            fig7, ax7 = plt.subplots(figsize=(8, 9))
            plot_forc_diagram_standard(result, ax=ax7, show_points=False, colormap='inferno_r', bc_lim=bc_lim, bu_lim=bu_lim)
            ax7.set_title(f'FORC Contour (inferno_r): {d.samplename}')
            plt.tight_layout()
    
    # Option 7: FORCinel workflow
    elif choice == '7':
        # FORCinel v3 workflow
        print("\n  FORCinel v3 Workflow")
        print("  " + "="*50)
        
        # Set grid resolution
        grid_size = select_forc_grid_resolution()
        
        # Set plot resolution
        dpi = set_plot_dpi()
        
        # Get parameters
        try:
            sf_raw = input("  Smoothing factor (2-4, default 3): ").strip()
            smoothing_factor = int(sf_raw) if sf_raw else 3
            
            mult_raw = input("  Moment multiplier (default 1e6 for typical samples): ").strip()
            multiplier = float(mult_raw) if mult_raw else 1e6
            
            output_dir = input("  Output directory [current dir]: ").strip() or "."
        except ValueError:
            print("  Invalid input, using defaults.")
            smoothing_factor = 3
            multiplier = 1e6
            output_dir = "."
        
        for d, forc_complete in forc_complete_list:
            print(f"\n  Processing {d.samplename}...")
            
            # Process FORC data
            try:
                result = process_forc_forcinel_workflow(d, smoothing_factor=smoothing_factor, grid_size=grid_size)
                
                # Prompt for axis limits
                limits = prompt_forc_axis_limits(result)
                bc_lim = limits['bc_lim'] if limits else None
                bu_lim = limits['bu_lim'] if limits else None
                
                # Plot FORCinel-style diagram
                fig, ax = plt.subplots(figsize=(12, 10))
                plot_forc_diagram_standard(result, ax=ax, show_points=True, colormap='hot_r',
                                          bc_lim=bc_lim, bu_lim=bu_lim)
                ax.set_title(f'FORC Diagram: {d.samplename} (SF={smoothing_factor})')
                plt.show()
                
                # Save plot
                plot_path = os.path.join(output_dir, f"{d.samplename}_forc_SF{smoothing_factor}.png")
                fig.savefig(plot_path, dpi=150, bbox_inches='tight')
                print(f"    ✓ Saved plot: {plot_path}")
                
                # Export .frc file
                print(f"    Note: .frc export requires original file path")
                print(f"          Use 'convert_rmg_to_forc()' function directly")
                
            except Exception as e:
                print(f"    ✗ Error processing {d.samplename}: {e}")
    
    if MATPLOTLIB_OK and choice in ('1', '2', '3', '4', '5', '6'):
        plt.show()


def generate_forc_measurement():
    """Generate a .rmg template file for FORC measurement."""
    print("\n── Generate FORC Measurement Script ───────────────────────")
    print("\nThis creates a .rmg template file that the magnetometer can load")
    print("to run a FORC measurement sequence.\n")
    
    try:
        start = float(input("  Start field (mT) [default: 0]: ").strip() or "0")
        stop = float(input("  Stop field (mT) [default: 500]: ").strip() or "500")
        
        # Ask for linear or exponential spacing
        print("\n  Field Spacing:")
        print("  + [1] Linear (constant step size)")
        print("  + [2] Exponential (growing step size)")
        spacing_choice = input("  Choice [1]: ").strip() or "1"
        exponential = (spacing_choice == '2')
        
        if exponential:
            step = float(input("  Minimum step size (mT) [default: 5]: ").strip() or "5")
            exp_base = float(input("  Exponential base (typically 1.2-1.5) [default: 1.3]: ").strip() or "1.3")
        else:
            step = float(input("  Step size (mT) [default: 5]: ").strip() or "5")
            exp_base = 1.3  # Not used for linear
            
    except ValueError:
        print("  Invalid input.")
        return
    
    # Ask for saturation or standard FORC
    print("\n  FORC Type:")
    print("  + [1] Standard FORC (measure up to reversal field)")
    print("  + [2] Saturation FORC (measure full curve at each reversal)")
    choice = input("  Choice [1]: ").strip() or "1"
    saturation = (choice == '2')
    
    # Output path
    output_dir = input("\n  Output directory [current dir]: ").strip() or "."
    
    # Generate descriptive filename with all settings
    if exponential:
        # Exponential: include min step and base
        if saturation:
            filename = f"FORC_{start:.1f}-{stop:.1f}mT_exp{exp_base:.2f}_min{step:.1f}mT_SATURATION.rmg"
        else:
            filename = f"FORC_{start:.1f}-{stop:.1f}mT_exp{exp_base:.2f}_min{step:.1f}mT_standard.rmg"
    else:
        # Linear: include step size
        if saturation:
            filename = f"FORC_{start:.1f}-{stop:.1f}mT_linear_{step:.1f}mT_SATURATION.rmg"
        else:
            filename = f"FORC_{start:.1f}-{stop:.1f}mT_linear_{step:.1f}mT_standard.rmg"
    
    output_path = os.path.join(output_dir, filename)
    
    try:
        result = generate_forc_script(start, stop, step, output_path, saturation, 
                                      exponential, exp_base)
        print(f"\n  ✓ Generated FORC script: {result}")
        
        # Calculate number of measurements
        num_curves = int((stop - start) / step) + 1
        if saturation:
            points_per_curve = int((stop - start) / step) + 1
        else:
            points_per_curve = (num_curves + 1) // 2  # average
        total_points = num_curves * points_per_curve
        
        print(f"  ✓ {num_curves} reversal fields")
        print(f"  ✓ ~{total_points} total measurement points")
    except Exception as e:
        print(f"  ✗ Error: {e}")


def export_forc_to_generic(loaded):
    """Convert loaded FORC data to generic .frc format."""
    if not loaded:
        print("  No samples loaded.")
        return
    
    print("\n── Export FORC to Generic Format ──────────────────────────")
    show_loaded(loaded)
    
    raw = input("\n  Enter sample number to export: ").strip()
    try:
        n = int(raw) - 1
        if 0 <= n < len(loaded):
            d = loaded[n]
        else:
            print("  Invalid sample number.")
            return
    except ValueError:
        print("  Invalid input.")
        return
    
    # Check for FORC data
    forc = rmg_extract_forc_data(d)
    if not forc.exists:
        print(f"  ✗ No FORC data found in {d.samplename}")
        return
    
    # Get multiplier
    mult_raw = input("  Multiplier for moment values [1.0]: ").strip()
    multiplier = float(mult_raw) if mult_raw else 1.0
    
    # Output path
    output_dir = input("  Output directory [current dir]: ").strip() or "."
    output_name = f"{d.samplename}.frc"
    output_path = os.path.join(output_dir, output_name)
    
    try:
        result_path, data_array = convert_rmg_to_forc(
            d.samplename,  # This won't work, need actual path
            output_path,
            multiplier
        )
        print(f"\n  ✓ Exported {len(forc.curves)} curves to {result_path}")
        print(f"  ✓ Total {len(data_array)} data points")
    except Exception as e:
        print(f"  ✗ Error: {e}")

def show_loaded(loaded):
    if not loaded:
        print("  (no samples loaded)")
    else:
        for i, d in enumerate(loaded):
            print(f"    [{i+1}] {d.samplename}")




def create_sam_files():
    """Generate .sam header files for paleomagnetic analysis."""
    print("\n" + "="*70)
    print("SAM HEADER FILE GENERATION")
    print("="*70)
    print("\nOptions:")
    print("  [1] Create blank SAM template (no calculations)")
    print("  [2] Interactive SAM creation (with sun compass & IGRF)")
    print("  [3] Return to main menu")
    
    choice = input("\nChoice [1]: ").strip() or '1'
    
    if choice == '1':
        # Blank template with detailed prompts
        generate_blank_sam_template(interactive=True)
    
    elif choice == '2':
        # Interactive with calculations
        create_sam_interactive()
    
    elif choice == '3':
        return
    else:
        print("Invalid choice")

def main():
    print(BANNER)

    loaded = []   # global RmgData equivalent

    # ── if files were passed on the command line, load them straight away ──
    if len(sys.argv) > 1:
        print("Loading command-line files…")
        load_direct(sys.argv[1:], loaded)

    # ── interactive menu loop ──────────────────────────────────────────────
    while True:
        print("\n▼═══════════════════════════════════════════════════════════▼")
        show_loaded(loaded)

        print("\n▼═══════════════════════════════════════════════════════════▼")
        print("  + [L] Load more files from a directory")
        print("  + [X] Remove samples from loaded list")
        print("  + [R] Run selected routines (batch plotter)")
        print("  + [F] 8-Panel Rock-Mag Dashboard (IRM, ARM, Backfield, etc.)")
        print("  + [I] Inspect one sample (show detailed stats)")
        print("  + [H] Show coercivity values (Hcr, MDF) for all samples")
        print("  + [E] Export all statistics to table")
        print("  + [P] Plot FORC data (First Order Reversal Curves)")
        print("  + [G] Generate FORC measurement script")
        print("  + [S] Generate SAM header files (paleomag)")
        print("  + [C] Clear all loaded samples")
        print("  + [Q] Quit")

        choice = input("\n  Choice: ").strip().upper()

        if choice == 'L':
            browse_and_load(loaded)

        elif choice == 'X':
            remove_samples(loaded)

        elif choice == 'R':
            if not loaded:
                print("  No samples loaded.")
                continue
            
            routines = select_routines()
            opts     = select_options()
            groups   = group_samples(loaded)
            
            print(f"\nRunning batch plotter for {len(groups)} group(s)…")
            all_figs = []
            
            for g_idx, grp in enumerate(groups, start=1):
                multisample = len(grp) > 1
                grp_name = f"group{g_idx}" if len(groups) > 1 else ""
                prefix = opts['file_prefix'] + grp_name + ("-" if grp_name else "")
                
                figs = rmg_batch_plotter(
                    grp,
                    routines,
                    multisample      = multisample,
                    subplots         = opts['subplots'],
                    autosave_formats = opts['autosave_formats'],
                    file_prefix      = prefix,
                )
                all_figs.extend(figs)
                
                if opts['stats_table']:
                    stem = prefix + 'rockmagstats-summary'
                    rmg_stats_write_table(grp, stem)
            
            if all_figs and MATPLOTLIB_OK:
                print(f"\n{len(all_figs)} figure(s) ready. Showing (close all windows to continue)…")
                plt.show()

        elif choice == 'F':
            if not loaded:
                print("  No samples loaded.")
                continue
            
            groups = group_samples(loaded)
            
            print(f"\nGenerating 8-panel dashboard for {len(groups)} group(s)…")
            all_figs = []
            
            for g_idx, grp in enumerate(groups, start=1):
                fig = rmg_data_full_analysis(grp)
                all_figs.append(fig)
                
                raw = input(f"  Save group {g_idx} dashboard as PNG? [y/N] ").strip().lower()
                if raw in ('y', 'yes'):
                    if len(grp) == 1:
                        name = grp[0].samplename
                    else:
                        name = f'group{g_idx}_multisamples'
                    fig.savefig(f'{name}_rockmag.png', dpi=150, bbox_inches='tight')
                    print(f"    ✓ Saved → {name}_rockmag.png")
            
            if all_figs and MATPLOTLIB_OK:
                plt.show()

        elif choice == 'I':
            inspect_sample(loaded)

        elif choice == 'H':
            show_coercivities(loaded)

        elif choice == 'E':
            export_all_stats(loaded)

        elif choice == 'P':
            plot_forc_data(loaded)

        
        elif choice == 'S':
            create_sam_files()
        
        elif choice == 'G':
            generate_forc_measurement()

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
