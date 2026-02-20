"""
rmg_forc.py
-----------
FORC (First Order Reversal Curve) analysis tools for rock magnetism.

FORC measurements involve:
1. Saturating the sample in a positive field (Ha_max)
2. Reversing to a reversal field (Hr)
3. Sweeping back up through increasing fields (Ha)
4. Measuring magnetization M(Hr, Ha) at each point

This creates a 2D dataset that reveals magnetic domain interactions
and coercivity distributions.

Functions:
- rmg_extract_forc_data: Extract FORC curves from .rmg file
- rmg_plot_forc_curves: Plot raw FORC curves
- rmg_forc_diagram: Calculate and plot FORC distribution
- generate_forc_script: Generate .rmg template for FORC measurement
- convert_rmg_to_forc: Export FORC data to generic .frc format

Copyright (C) 2025 Python port, GNU GPL v.3
"""

import numpy as np
import os
from typing import List, Tuple, Optional


def rmg_extract_forc_data(data):
    """
    Extract FORC curves from an RmgData object.
    
    FORC measurements use step type 'FORCz' or 'IRMz' with fields in data.levels.
    Pattern: negative reversal field followed by positive sweep fields.
    
    Parameters
    ----------
    data : RmgData
    
    Returns
    -------
    FORCData object with curves list
    """
    class FORCData:
        pass
    
    y = FORCData()
    y.exists = False
    y.curves = []
    
    # Find FORC steps
    forc_steps = []
    for i, stype in enumerate(data.steptypes):
        if stype.lower() in ('forcz', 'irmz', 'irm'):
            forc_steps.append(i)
    
    if len(forc_steps) == 0:
        return y
    
    # CRITICAL: Fields are in data.levels (raw Oe), not treatmentDCFields
    OE_TO_T = 1e-4
    
    current_curve = None
    skip_first_positive = True
    
    for i in forc_steps:
        field_oe = data.levels[i]
        field_T = field_oe * OE_TO_T
        mz = data.mvector[2, i]
        
        # Skip initial saturation pulse
        if skip_first_positive and field_oe > 0:
            skip_first_positive = False
            continue
        
        # Negative = new curve (reversal)
        if field_oe < 0:
            if current_curve is not None and len(current_curve['Ha']) > 0:
                y.curves.append(current_curve)
            
            current_curve = {
                'Hr': field_T,
                'Ha': [],
                'M': []
            }
        
        # Positive/zero = measurement point
        elif current_curve is not None:
            current_curve['Ha'].append(field_T)
            current_curve['M'].append(mz)
    
    if current_curve is not None and len(current_curve['Ha']) > 0:
        y.curves.append(current_curve)
    
    for curve in y.curves:
        curve['Ha'] = np.array(curve['Ha'])
        curve['M'] = np.array(curve['M'])
    
    y.exists = len(y.curves) > 0
    return y



def rmg_plot_forc_curves(data_list, ax=None):
    """
    Plot raw FORC curves (M vs Ha for each Hr).
    
    Parameters
    ----------
    data_list : list of RmgData or single RmgData
    ax : matplotlib Axes, optional
    """
    import matplotlib.pyplot as plt
    
    if not isinstance(data_list, list):
        data_list = [data_list]
    
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))
    
    colors = plt.cm.viridis(np.linspace(0, 1, 20))
    
    for data in data_list:
        forc = rmg_extract_forc_data(data)
        
        if not forc.exists:
            continue
        
        for i, curve in enumerate(forc.curves):
            color = colors[i % len(colors)]
            label = f"Hr = {curve['Hr']*1000:.1f} mT" if i < 10 else None
            ax.plot(curve['Ha'] * 1000, curve['M'] * 1e9, 
                   'o-', color=color, markersize=3, linewidth=1, label=label)
    
    if ax.lines:
        ax.set_xlabel('Applied Field Ha (mT)')
        ax.set_ylabel('Magnetization M (nAm²)')
        ax.set_title('FORC Curves')
        ax.grid(True, alpha=0.3)
        if len(forc.curves) <= 10:
            ax.legend(fontsize='small', loc='best')
    else:
        ax.text(0.5, 0.5, 'No FORC data', 
               ha='center', va='center', transform=ax.transAxes)
        ax.set_xticks([])
        ax.set_yticks([])
    
    return ax


def rmg_forc_diagram(data, smoothing_factor=3, ax=None):
    """
    Calculate and plot FORC distribution diagram.
    
    The FORC distribution ρ(Hr, Hc) is the mixed second derivative:
        ρ = -0.5 * ∂²M / ∂Hr∂Ha
    
    Transformed to (Hc, Hu) coordinates:
        Hc = (Ha - Hr) / 2  (coercivity)
        Hu = (Ha + Hr) / 2  (interaction field)
    
    Parameters
    ----------
    data : RmgData
    smoothing_factor : int
        Smoothing window size
    ax : matplotlib Axes, optional
    
    Returns
    -------
    dict with keys:
        'Hc': ndarray - coercivity axis
        'Hu': ndarray - interaction field axis  
        'rho': 2D ndarray - FORC distribution
    """
    import matplotlib.pyplot as plt
    from scipy.ndimage import gaussian_filter
    from scipy.interpolate import griddata
    
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))
    
    forc = rmg_extract_forc_data(data)
    
    if not forc.exists or len(forc.curves) < 3:
        ax.text(0.5, 0.5, 'Insufficient FORC data\n(need ≥3 curves)', 
               ha='center', va='center', transform=ax.transAxes)
        ax.set_xticks([])
        ax.set_yticks([])
        return None
    
    # Extract all (Hr, Ha, M) points
    Hr_all = []
    Ha_all = []
    M_all = []
    
    for curve in forc.curves:
        n = len(curve['Ha'])
        Hr_all.extend([curve['Hr']] * n)
        Ha_all.extend(curve['Ha'])
        M_all.extend(curve['M'])
    
    Hr_all = np.array(Hr_all)
    Ha_all = np.array(Ha_all)
    M_all = np.array(M_all)
    
    # Create regular grid for interpolation
    Hr_min, Hr_max = Hr_all.min(), Hr_all.max()
    Ha_min, Ha_max = Ha_all.min(), Ha_all.max()
    
    n_grid = 100
    Hr_grid = np.linspace(Hr_min, Hr_max, n_grid)
    Ha_grid = np.linspace(Ha_min, Ha_max, n_grid)
    Hr_mesh, Ha_mesh = np.meshgrid(Hr_grid, Ha_grid)
    
    # Interpolate M onto regular grid
    M_grid = griddata((Hr_all, Ha_all), M_all, (Hr_mesh, Ha_mesh), method='cubic')
    
    # Calculate mixed second derivative: -0.5 * d²M/dHr/dHa
    # Use central differences
    dHr = Hr_grid[1] - Hr_grid[0]
    dHa = Ha_grid[1] - Ha_grid[0]
    
    # First derivative with respect to Hr
    dM_dHr = np.gradient(M_grid, dHr, axis=1)
    
    # Second derivative with respect to Ha
    d2M_dHrdHa = np.gradient(dM_dHr, dHa, axis=0)
    
    # FORC distribution
    rho_raw = -0.5 * d2M_dHrdHa
    
    # Multiply by -1 to flip sign for visualization
    # This makes strong FORC signals (large magnitude) appear as high values (yellow)
    # and weak signals appear as low values (purple/black)
    rho = -1.0 * rho_raw
    
    # Apply smoothing
    if smoothing_factor > 0:
        rho = gaussian_filter(rho, sigma=smoothing_factor)
    
    # Transform to (Hc, Hu) coordinates
    Hc_mesh = (Ha_mesh - Hr_mesh) / 2  # coercivity
    Hu_mesh = (Ha_mesh + Hr_mesh) / 2  # interaction field
    
    # Plot FORC diagram
    # Mask region where Ha < Hr (unphysical)
    mask = Ha_mesh < Hr_mesh
    rho_masked = np.ma.masked_where(mask, rho)
    
    levels = np.linspace(rho[~mask].min(), rho[~mask].max(), 20)
    
    contour = ax.contourf(Hc_mesh * 1000, Hu_mesh * 1000, rho_masked, 
                          levels=levels, cmap='RdBu_r', extend='both')
    
    ax.contour(Hc_mesh * 1000, Hu_mesh * 1000, rho_masked,
              levels=levels, colors='k', linewidths=0.5, alpha=0.3)
    
    ax.set_xlabel('Coercivity Hc (mT)')
    ax.set_ylabel('Interaction Field Hu (mT)')
    ax.set_title(f'FORC Diagram: {data.samplename}')
    ax.set_aspect('equal')
    
    plt.colorbar(contour, ax=ax, label='FORC distribution ρ')
    
    return {
        'Hc': Hc_mesh * 1000,
        'Hu': Hu_mesh * 1000,
        'rho': rho
    }


def generate_forc_script(start_mT: float, stop_mT: float, step_mT: float,
                        output_path: str, saturation: bool = False):
    """
    Generate a .rmg template file for FORC measurement.
    
    Creates a script that the magnetometer can load to run FORC measurements.
    
    Parameters
    ----------
    start_mT : float
        Starting field (mT)
    stop_mT : float
        Maximum field (mT)
    step_mT : float
        Field step size (mT)
    output_path : str
        Path to save .rmg file
    saturation : bool
        If True, measure full saturation curve at each reversal
        If False, only measure up to the reversal field (standard FORC)
    
    Returns
    -------
    str : path to created file
    """
    # Header line (column names)
    header = "Level,Bias Field (G),Spin Speed (rps),Hold Time (s),Mz (emu),Std. Dev. Z,Mz/Vol,Moment Susceptibility (emu/Oe),Mx (emu),Std. Dev. X,My (emu),Std. Dev. Y\r\n"
    
    # Template line (placeholder values, instrument fills in real data)
    dummy_vals = " , ".join(["1.00000000000001E-09"] * 12)
    
    # Start with AFmax baseline
    output = header
    output += f"AFmax, 1250 , 0 , 0 , 0 , {dummy_vals}\r\n"
    
    # NRM baseline
    output += f"NRM, 0 , 0 , 0 , 0 , {dummy_vals}\r\n"
    
    # Generate FORC curves
    # Field values are in Oersted (multiply mT by 10)
    curve_max_vals = np.arange(start_mT, stop_mT + step_mT, step_mT)
    
    for curve_max in curve_max_vals:
        # Start each curve with reversal field (negative)
        output += f"FORCz, {int(curve_max * -10)} , 0 , 0 , 0 , {dummy_vals}\r\n"
        
        # Measurement sweep (positive fields)
        if saturation:
            # Go all the way to stop field (saturation FORC)
            field_vals = np.arange(start_mT, stop_mT + step_mT, step_mT)
        else:
            # Go only up to reversal field (standard FORC)
            field_vals = np.arange(start_mT, curve_max + step_mT, step_mT)
        
        for field in field_vals:
            output += f"FORCz, {int(field * 10)} , 0 , 0 , 0 , {dummy_vals}\r\n"
        
        # Add AFmax at end of each curve for saturation FORC
        if saturation:
            output += f"AFmax, 1250 , 0 , 0 , 0 , {dummy_vals}\r\n"
    
    # Final AFmax
    if not saturation:
        output += f"AFmax, 1250 , 0 , 0 , 0 , {dummy_vals}\r\n"
    
    # Write to file
    with open(output_path, 'w', newline='') as f:
        f.write(output)
    
    return output_path


def convert_rmg_to_forc(rmg_path: str, output_path: Optional[str] = None, 
                       multiplier: float = 1.0):
    """
    Convert .rmg FORC file to generic .frc format (FORCinel v3-compatible).
    
    FORCinel v3.0.6 requirements:
    - Generic FORC format (.frc extension)
    - Six digits precision (no scientific notation)
    - Fields in Tesla, moments scaled appropriately
    - Blank lines separate curves
    - Ends with "END"
    
    Generic FORC format:
        Ha,M
        Ha,M
        ...
        [blank line]  ← separates curves
        Ha,M
        ...
        END
    
    Parameters
    ----------
    rmg_path : str
        Path to .rmg file with FORC data
    output_path : str, optional
        Output .frc file path (default: same name with .frc extension)
    multiplier : float
        Multiply moment values by this factor.
        FORCinel requires 6-digit precision without scientific notation.
        For moments in 10⁻⁶ to 10⁻⁴ emu range: use 1e6
        For stronger samples: adjust accordingly
        For reversed polarity: use negative multiplier
    
    Returns
    -------
    tuple: (output_path, array of [Ha, M] data)
    """
    from rmg_import import rmg_import
    
    # Load the data
    data = rmg_import(rmg_path)
    
    # Extract FORC curves
    forc = rmg_extract_forc_data(data)
    
    if not forc.exists:
        raise ValueError("No FORC data found in file")
    
    # Build output path
    if output_path is None:
        base = os.path.splitext(rmg_path)[0]
        output_path = base + '.frc'
    
    # Collect all data points
    output_lines = []
    forc_array = []
    
    for curve in forc.curves:
        for ha, m in zip(curve['Ha'], curve['M']):
            ha_T = ha
            m_scaled = m * multiplier
            output_lines.append(f"{ha_T:.6f},{m_scaled:.6f}")
            forc_array.append([ha_T, m_scaled])
        # Blank line between curves
        output_lines.append("")
    
    output_lines.append("END")
    
    # Write file
    with open(output_path, 'w') as f:
        f.write('\r'.join(output_lines))
    
    return output_path, np.array(forc_array)


def process_forc_forcinel_workflow(data, smoothing_factor: int = 3, grid_size: int = 100,
                                   output_dir: str = None,
                                   samplename: str = None) -> dict:
    """
    Complete FORCinel v3-compatible FORC processing workflow.
    
    Implements the workflow from "How to analyze remFORC data v3":
    1. Extract FORC curves from data
    2. Interpolate to regular grid
    3. Apply smoothing (FORCinel default: SF=3 or 4)
    4. Calculate FORC distribution ρ = -0.5 × ∂²M/∂Hr∂Ha
    5. Transform to (Hc, Hu) coordinates
    
    Parameters
    ----------
    data : RmgData
        Loaded .rmg data
    smoothing_factor : int
        FORCinel smoothing parameter (default: 3)
        FORCinel recommends 2-4, with 2 for noisy data
    grid_size : int
        Grid resolution (default: 100)
        Higher values (200, 300) give finer detail but slower processing
    output_dir : str, optional
        Directory for output files
    samplename : str, optional
        Sample name for output files
    
    Returns
    -------
    dict with processed FORC data:
        'forc_data': extracted curves
        'Hr_grid', 'Ha_grid': field grids (Tesla)
        'M_grid': magnetization grid (Am²)
        'Hc', 'Hu': transformed coordinates (Tesla)
        'rho': FORC distribution
        'valid_mask': region where Ha >= Hr
    """
    from scipy.interpolate import griddata
    from scipy.ndimage import uniform_filter
    
    # Extract FORC curves
    forc = rmg_extract_forc_data(data)
    
    if not forc.exists:
        raise ValueError("No FORC data found")
    
    # Collect all data points
    Hr_all = []
    Ha_all = []
    M_all = []
    
    for curve in forc.curves:
        n = len(curve['Ha'])
        Hr_all.extend([curve['Hr']] * n)
        Ha_all.extend(curve['Ha'])
        M_all.extend(curve['M'])
    
    Hr_all = np.array(Hr_all)
    Ha_all = np.array(Ha_all)
    M_all = np.array(M_all)
    
    # Create regular grid (default 100x100, user can increase)
    # grid_size is now a parameter
    Hr_min, Hr_max = Hr_all.min(), Hr_all.max()
    Ha_min, Ha_max = Ha_all.min(), Ha_all.max()
    
    Hr_1d = np.linspace(Hr_min, Hr_max, grid_size)
    Ha_1d = np.linspace(Ha_min, Ha_max, grid_size)
    Hr_grid, Ha_grid = np.meshgrid(Hr_1d, Ha_1d, indexing='ij')
    
    # Interpolate magnetization onto grid
    M_grid = griddata(
        (Hr_all, Ha_all),
        M_all,
        (Hr_grid, Ha_grid),
        method='cubic',
        fill_value=np.nan
    )
    
    # Apply smoothing (FORCinel approach)
    if smoothing_factor > 0:
        M_smooth = uniform_filter(M_grid, size=smoothing_factor, mode='nearest')
    else:
        M_smooth = M_grid.copy()
    
    # Calculate FORC distribution: ρ = -0.5 × ∂²M/∂Hr∂Ha
    dHr = Hr_grid[1, 0] - Hr_grid[0, 0]
    dHa = Ha_grid[0, 1] - Ha_grid[0, 0]
    
    # Mixed second derivative
    dM_dHr = np.gradient(M_smooth, dHr, axis=0)
    d2M_dHrdHa = np.gradient(dM_dHr, dHa, axis=1)
    rho_raw = -0.5 * d2M_dHrdHa
    
    # Multiply by -1 to flip sign for visualization
    # Strong FORC signals (large magnitude) → high values (yellow)
    # Weak signals (small magnitude) → low values (purple/black)
    rho = -1.0 * rho_raw
    
    # Transform to (Hc, Hu) coordinates
    Hc = (Ha_grid - Hr_grid) / 2.0
    Hu = (Ha_grid + Hr_grid) / 2.0
    
    # Mask unphysical region (Ha < Hr)
    valid_mask = Ha_grid >= Hr_grid
    
    return {
        'forc_data': forc,
        'Hr_grid': Hr_grid,
        'Ha_grid': Ha_grid,
        'M_grid': M_grid,
        'M_smooth': M_smooth,
        'Hc': Hc,
        'Hu': Hu,
        'rho': rho,
        'valid_mask': valid_mask,
        'smoothing_factor': smoothing_factor
    }


def plot_forc_forcinel_style(result: dict, ax=None, colormap='RdBu_r', 
                             levels: int = 20, show_colorbar: bool = True):
    """
    Plot FORC diagram in FORCinel v3 style.
    
    Displays bottom half of FORC diamond (Hu <= 0) as per FORCinel.
    
    Parameters
    ----------
    result : dict
        Output from process_forc_forcinel_workflow
    ax : matplotlib Axes
    colormap : str
        Colormap (FORCinel uses RdBu_r style)
    levels : int
        Number of contour levels
    show_colorbar : bool
        Display colorbar
    """
    import matplotlib.pyplot as plt
    
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 7))
    
    # Convert to mT for plotting
    Hc = result['Hc'] * 1000
    Hu = result['Hu'] * 1000
    rho = result['rho']
    mask = result['valid_mask']
    
    # Mask invalid region and focus on lower half (Hu <= 0)
    rho_masked = np.ma.masked_where(~mask | (Hu > 0), rho)
    
    # Determine contour levels from valid data
    valid_rho = rho[mask & (Hu <= 0)]
    if np.any(np.isfinite(valid_rho)):
        vmin, vmax = np.percentile(valid_rho[np.isfinite(valid_rho)], [2, 98])
    else:
        vmin, vmax = -1, 1
    
    # Plot filled contours
    contour = ax.contourf(Hc, Hu, rho_masked, 
                          levels=levels, cmap=colormap, 
                          vmin=vmin, vmax=vmax, extend='both')
    
    # Add contour lines
    ax.contour(Hc, Hu, rho_masked, 
              levels=levels, colors='k', linewidths=0.5, alpha=0.3)
    
    ax.set_xlabel('Coercivity Hc (mT)')
    ax.set_ylabel('Interaction Field Hu (mT)')
    ax.set_title('FORC Diagram (FORCinel-style)')
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.2, linestyle='--')
    ax.axhline(0, color='k', linewidth=0.8, alpha=0.5)
    ax.axvline(0, color='k', linewidth=0.8, alpha=0.5)
    
    if show_colorbar:
        plt.colorbar(contour, ax=ax, label='FORC Distribution ρ')
    
    return ax
"""
Fixed FORC plotting to match standard FORC diagram output.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize






def prompt_forc_axis_limits(result):
    """
    Prompt user for custom Bc and Bu axis limits BEFORE plotting.
    
    Shows the data range and asks if user wants to set custom limits.
    
    Parameters
    ----------
    result : dict
        Result from process_forc_forcinel_workflow containing Bc and Bu data
    
    Returns
    -------
    dict or None : Dictionary with 'bc_lim' and 'bu_lim' tuples, or None for auto
    """
    import numpy as np
    
    # Extract Bc and Bu from result
    # The result uses 'Hc' and 'Hu' as keys
    bc_data = result['Hc']
    bu_data = result['Hu']
    
    # Calculate data ranges (data is in Tesla, convert to mT for display)
    bc_min, bc_max = np.nanmin(bc_data) * 1000, np.nanmax(bc_data) * 1000
    bu_min, bu_max = np.nanmin(bu_data) * 1000, np.nanmax(bu_data) * 1000
    
    print("\n▼══════════════════════════════════════════════════════════▼")
    print("           FORC Axis Range Settings")
    print("▲══════════════════════════════════════════════════════════▲")
    print(f"\nData range:")
    print(f"  Bc (coercivity):      {bc_min:.1f} to {bc_max:.1f} mT")
    print(f"  Bu (interaction):     {bu_min:.1f} to {bu_max:.1f} mT")
    print("\nNote: Enter limits in mT (milliTesla)")
    
    # Ask if user wants custom limits
    print("\nOptions:")
    print("  [1] Auto-scale (use full data range)")
    print("  [2] Set custom Bc and Bu limits")
    
    choice = input("\nChoice [1]: ").strip() or "1"
    
    if choice != "2":
        print("✓ Using auto-scale")
        return None
    
    # Get custom limits
    print("\n--- Set Custom Limits ---")
    print("(Press Enter to keep default)")
    
    try:
        print("\nBc Axis (Coercivity):")
        bc_min_input = input(f"  Min (mT) [{bc_min:.1f}]: ").strip()
        bc_min_new = float(bc_min_input) if bc_min_input else bc_min
        
        bc_max_input = input(f"  Max (mT) [{bc_max:.1f}]: ").strip()
        bc_max_new = float(bc_max_input) if bc_max_input else bc_max
        
        print("\nBu Axis (Interaction Field):")
        bu_min_input = input(f"  Min (mT) [{bu_min:.1f}]: ").strip()
        bu_min_new = float(bu_min_input) if bu_min_input else bu_min
        
        bu_max_input = input(f"  Max (mT) [{bu_max:.1f}]: ").strip()
        bu_max_new = float(bu_max_input) if bu_max_input else bu_max
        
        print(f"\n✓ Custom limits: Bc [{bc_min_new:.1f}, {bc_max_new:.1f}] mT, Bu [{bu_min_new:.1f}, {bu_max_new:.1f}] mT")
        
        # Convert from mT back to Tesla for internal use (plot function will convert to mT again)
        # NO - keep in mT since plot_forc_diagram_standard expects mT after *1000 conversion
        return {
            'bc_lim': (bc_min_new, bc_max_new),
            'bu_lim': (bu_min_new, bu_max_new)
        }
        
    except ValueError:
        print("✗ Invalid input, using auto-scale")
        return None

def adjust_forc_axis_limits(result, current_bc_lim=None, current_bu_lim=None):
    """
    Interactive axis limit adjustment for FORC diagrams.
    
    Prompts user to set custom Bc and Bu axis limits after viewing initial plot.
    
    Parameters
    ----------
    result : dict
        Result from process_forc_forcinel_workflow containing Bc and Bu data
    current_bc_lim : tuple, optional
        Current Bc limits (min, max) in mT
    current_bu_lim : tuple, optional
        Current Bu limits (min, max) in mT
    
    Returns
    -------
    dict : Dictionary with 'bc_lim' and 'bu_lim' tuples, or None if no adjustment
    """
    import numpy as np
    
    # Extract Bc and Bu from result
    # The result uses 'Hc' and 'Hu' as keys
    bc_data = result['Hc']
    bu_data = result['Hu']
    
    # Calculate current data ranges
    bc_min, bc_max = np.nanmin(bc_data), np.nanmax(bc_data)
    bu_min, bu_max = np.nanmin(bu_data), np.nanmax(bu_data)
    
    # Show current ranges
    print("\n▼══════════════════════════════════════════════════════════▼")
    print("           FORC Axis Limit Adjustment")
    print("▲══════════════════════════════════════════════════════════▲")
    
    if current_bc_lim:
        print(f"\nCurrent Bc range: {current_bc_lim[0]:.1f} - {current_bc_lim[1]:.1f} mT")
    else:
        print(f"\nCurrent Bc range (auto): {bc_min:.1f} - {bc_max:.1f} mT")
    
    if current_bu_lim:
        print(f"Current Bu range: {current_bu_lim[0]:.1f} - {current_bu_lim[1]:.1f} mT")
    else:
        print(f"Current Bu range (auto): {bu_min:.1f} - {bu_max:.1f} mT")
    
    # Ask if user wants to adjust
    adjust = input("\n+ Adjust axis limits? (y/n) [n]: ").strip().lower()
    
    if adjust != 'y':
        return None
    
    # Get new Bc limits
    print("\n--- Bc Axis (Coercivity) ---")
    try:
        bc_min_input = input(f"Bc min (mT) [{bc_min:.1f}]: ").strip()
        bc_min_new = float(bc_min_input) if bc_min_input else bc_min
        
        bc_max_input = input(f"Bc max (mT) [{bc_max:.1f}]: ").strip()
        bc_max_new = float(bc_max_input) if bc_max_input else bc_max
    except ValueError:
        print("Invalid input, using auto-scale for Bc")
        bc_min_new, bc_max_new = bc_min, bc_max
    
    # Get new Bu limits
    print("\n--- Bu Axis (Interaction Field) ---")
    try:
        bu_min_input = input(f"Bu min (mT) [{bu_min:.1f}]: ").strip()
        bu_min_new = float(bu_min_input) if bu_min_input else bu_min
        
        bu_max_input = input(f"Bu max (mT) [{bu_max:.1f}]: ").strip()
        bu_max_new = float(bu_max_input) if bu_max_input else bu_max
    except ValueError:
        print("Invalid input, using auto-scale for Bu")
        bu_min_new, bu_max_new = bu_min, bu_max
    
    print(f"\n✓ New limits: Bc [{bc_min_new:.1f}, {bc_max_new:.1f}] mT, Bu [{bu_min_new:.1f}, {bu_max_new:.1f}] mT")
    
    return {
        'bc_lim': (bc_min_new, bc_max_new),
        'bu_lim': (bu_min_new, bu_max_new)
    }

def plot_forc_diagram_standard(result, ax=None, show_points=True, 
                               colormap='hot', vmin=None, vmax=None,
                               bc_lim=None, bu_lim=None):
    """
    Plot FORC diagram in standard format matching published figures.
    
    Shows scatter plot of (Bc, Bu) points with FORC distribution values.
    Bc = coercivity = (Ha - Hr) / 2
    Bu = interaction field = (Ha + Hr) / 2
    
    Parameters
    ----------
    result : dict from process_forc_forcinel_workflow
    ax : matplotlib Axes
    show_points : bool
        If True, show scatter plot of data points
        If False, show contour plot
    colormap : str
        Best options (use '_r' for reversed):
        - 'inferno_r' : High contrast (yellow=high, purple=low) [RECOMMENDED]
        - 'plasma_r'  : Modern (yellow=high, purple=low)
        - 'hot_r'     : Classic (white=high, red=low)
        - 'RdBu_r'    : Diverging (red=pos, blue=neg)
    vmin, vmax : float
        Color scale limits
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 9))
    
    # Use Bc/Bu notation (standard in FORC literature)
    Bc = result['Hc']  # Already calculated as (Ha - Hr)/2
    Bu = result['Hu']  # Already calculated as (Ha + Hr)/2
    rho = result['rho']
    mask = result['valid_mask']
    
    # Flatten for scatter plot
    Bc_flat = Bc[mask].flatten()
    Bu_flat = Bu[mask].flatten()
    rho_flat = rho[mask].flatten()
    
    # Remove NaN
    finite = np.isfinite(rho_flat)
    Bc_flat = Bc_flat[finite]
    Bu_flat = Bu_flat[finite]
    rho_flat = rho_flat[finite]
    
    # Determine color scale
    if vmin is None or vmax is None:
        # Use tighter percentiles for better color contrast
        # 5-95 instead of 1-99 makes colors more vibrant
        vmin_auto = np.percentile(rho_flat, 5)
        vmax_auto = np.percentile(rho_flat, 95)
        if vmin is None:
            vmin = vmin_auto
        if vmax is None:
            vmax = vmax_auto
    
    if show_points:
        # Scatter plot - BRIGHT MODE for vibrant colors
        # Larger points (s=2-3) for better coverage
        # Higher alpha (0.8-1.0) for more saturated colors
        # Tighter vmin/vmax (5-95 percentile) for better contrast
        scatter = ax.scatter(Bc_flat, Bu_flat, c=rho_flat, 
                            cmap=colormap, s=3, alpha=0.9,
                            vmin=vmin, vmax=vmax)
    else:
        # Contour plot with FIXED color scale
        # CRITICAL FIX: Clip data to vmin/vmax to prevent interpolation artifacts
        rho_clipped = np.clip(rho, vmin, vmax)
        rho_clipped_masked = np.ma.masked_where(~mask, rho_clipped)
        
        # Create levels within the percentile range
        levels = np.linspace(vmin, vmax, 20)
        
        scatter = ax.contourf(Bc, Bu, rho_clipped_masked, 
                             levels=levels, cmap=colormap,
                             vmin=vmin, vmax=vmax, extend='neither')
        ax.contour(Bc, Bu, rho_clipped_masked, 
                  levels=levels, colors='k', linewidths=0.3, alpha=0.2)
    
    # Convert to mT for display (data is in Tesla)
    # Get current data limits
    current_xlim = ax.get_xlim()
    current_ylim = ax.get_ylim()
    
    # Convert Tesla to mT (multiply by 1000)
    # This affects both the data and the axis limits
    ax.set_xlim(current_xlim[0] * 1000, current_xlim[1] * 1000)
    ax.set_ylim(current_ylim[0] * 1000, current_ylim[1] * 1000)
    
    # Re-plot with mT units
    ax.clear()
    Bc_mT = Bc * 1000
    Bu_mT = Bu * 1000
    Bc_flat_mT = Bc_flat * 1000
    Bu_flat_mT = Bu_flat * 1000
    
    if show_points:
        # Scatter plot - BRIGHT MODE for vibrant colors
        scatter = ax.scatter(Bc_flat_mT, Bu_flat_mT, c=rho_flat, 
                            cmap=colormap, s=3, alpha=0.9,
                            vmin=vmin, vmax=vmax)
    else:
        rho_clipped = np.clip(rho, vmin, vmax)
        rho_clipped_masked = np.ma.masked_where(~mask, rho_clipped)
        levels = np.linspace(vmin, vmax, 20)
        scatter = ax.contourf(Bc_mT, Bu_mT, rho_clipped_masked, 
                             levels=levels, cmap=colormap,
                             vmin=vmin, vmax=vmax, extend='neither')
        ax.contour(Bc_mT, Bu_mT, rho_clipped_masked, 
                  levels=levels, colors='k', linewidths=0.3, alpha=0.2)
    
    # Formatting
    ax.set_xlabel('$B_c$ (mT)', fontsize=12)
    ax.set_ylabel('$B_u$ (mT)', fontsize=12)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
    ax.axhline(0, color='k', linewidth=0.8, linestyle='--', alpha=0.5)
    ax.axvline(0, color='k', linewidth=0.8, alpha=0.3)
    
    # Apply user-specified axis limits (in mT)
    if bc_lim is not None:
        ax.set_xlim(bc_lim)
    if bu_lim is not None:
        ax.set_ylim(bu_lim)
    
    # Colorbar
    cbar = plt.colorbar(scatter, ax=ax, label='Am²/T²')
    cbar.formatter.set_powerlimits((-3, 3))
    cbar.update_ticks()
    
    return ax


def plot_forc_raw_curves(forc_data, ax=None, max_curves=None):
    """
    Plot raw FORC curves (hysteresis plot).
    
    Matches Image 1: all reversal curves overlaid.
    
    Parameters
    ----------
    forc_data : output from rmg_extract_forc_data
    ax : matplotlib Axes
    max_curves : int
        Maximum number of curves to plot (None = all)
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 7))
    
    curves_to_plot = forc_data.curves
    if max_curves is not None:
        curves_to_plot = curves_to_plot[:max_curves]
    
    for i, curve in enumerate(curves_to_plot):
        Ha = curve['Ha']
        M = curve['M']
        ax.plot(Ha, M, 'b-', linewidth=0.5, alpha=0.7)
    
    ax.set_xlabel('Field (T)', fontsize=12)
    ax.set_ylabel('Remanence (Am^2)', fontsize=12)
    ax.set_title('Direct moment vs. Field - FORC')
    ax.grid(True, alpha=0.3)
    ax.axhline(0, color='k', linewidth=0.8, alpha=0.5)
    ax.axvline(0, color='k', linewidth=0.8, alpha=0.5)
    
    return ax
"""
Complete FORC functions with full hysteresis data.
"""

import numpy as np


def rmg_extract_forc_data_complete(data):
    """
    Extract FORC curves INCLUDING reversal points for complete hysteresis.
    
    Each curve contains:
    - Hr: reversal field (negative)
    - Mr: magnetization at reversal
    - Ha: applied fields (sweep back positive)
    - M: magnetization during sweep
    
    Returns full hysteresis data for plotting.
    """
    class FORCData:
        pass
    
    y = FORCData()
    y.exists = False
    y.curves = []
    
    # Find FORC steps
    forc_steps = []
    for i, stype in enumerate(data.steptypes):
        if stype.lower() in ('forcz', 'irmz', 'irm'):
            forc_steps.append(i)
    
    if len(forc_steps) == 0:
        return y
    
    OE_TO_T = 1e-4
    
    current_curve = None
    skip_first_positive = True
    
    for i in forc_steps:
        field_oe = data.levels[i]
        field_T = field_oe * OE_TO_T
        mz = data.mvector[2, i]
        
        # Skip initial saturation
        if skip_first_positive and field_oe > 0:
            skip_first_positive = False
            continue
        
        # Negative = new curve (reversal)
        if field_oe < 0:
            # Save previous curve
            if current_curve is not None and len(current_curve['Ha']) > 0:
                y.curves.append(current_curve)
            
            # Start new curve WITH the reversal point
            current_curve = {
                'Hr': field_T,      # reversal field
                'Mr': mz,           # magnetization AT reversal
                'Ha': [],           # applied fields (sweep)
                'M': []             # magnetization during sweep
            }
        
        # Positive/zero = measurement point
        elif current_curve is not None:
            current_curve['Ha'].append(field_T)
            current_curve['M'].append(mz)
    
    # Add last curve
    if current_curve is not None and len(current_curve['Ha']) > 0:
        y.curves.append(current_curve)
    
    # Convert to arrays
    for curve in y.curves:
        curve['Ha'] = np.array(curve['Ha'])
        curve['M'] = np.array(curve['M'])
    
    y.exists = len(y.curves) > 0
    return y


def plot_forc_hysteresis_complete(forc_data, ax=None):
    """
    Plot complete FORC hysteresis showing ALL reversal curves.
    
    This matches Image 1 - shows the full S-curve with all measurements.
    Each curve includes:
    - The reversal point (Hr, Mr)
    - The sweep back to positive (Ha, M)
    
    Parameters
    ----------
    forc_data : output from rmg_extract_forc_data_complete
    ax : matplotlib Axes
    """
    import matplotlib.pyplot as plt
    
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 7))
    
    # Collect ALL data points for the complete hysteresis plot
    all_fields = []
    all_moments = []
    
    for curve in forc_data.curves:
        # Add reversal point
        all_fields.append(curve['Hr'])
        all_moments.append(curve['Mr'])
        
        # Add sweep points
        all_fields.extend(curve['Ha'])
        all_moments.extend(curve['M'])
    
    all_fields = np.array(all_fields)
    all_moments = np.array(all_moments)
    
    # Plot as scatter to show all individual measurements
    ax.plot(all_fields, all_moments, 'b.', markersize=2, alpha=0.6)
    
    # Also plot individual curves to show the reversal pattern
    for curve in forc_data.curves:
        # Combine reversal + sweep for this curve
        curve_fields = np.concatenate([[curve['Hr']], curve['Ha']])
        curve_moments = np.concatenate([[curve['Mr']], curve['M']])
        ax.plot(curve_fields, curve_moments, 'b-', linewidth=0.3, alpha=0.5)
    
    ax.set_xlabel('Field (T)', fontsize=12)
    ax.set_ylabel('Remanence (Am²)', fontsize=12)
    ax.set_title('Direct moment vs. Field - FORC')
    ax.grid(True, alpha=0.3)
    ax.axhline(0, color='k', linewidth=0.8, alpha=0.5)
    ax.axvline(0, color='k', linewidth=0.8, alpha=0.5)
    
    # Add text box with info
    num_curves = len(forc_data.curves)
    total_points = len(all_fields)
    textstr = f'Curves: {num_curves}\nPoints: {total_points}'
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax.text(0.02, 0.98, textstr, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', bbox=props)
    
    return ax
