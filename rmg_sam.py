"""
rmg_sam.py
----------
Generate .sam header files and individual sample files for paleomagnetic analysis.

This module creates files compatible with the RAPID 2G magnetometer software
using the CIT (Caltech) format. Supports both manual entry and automatic
calculation of orientations using sun compass and IGRF magnetic field data.

Features:
- Generate blank .sam templates
- Sun compass orientation calculation
- IGRF-13 magnetic declination calculation
- Validation of sun vs magnetic compass data
- Batch processing of multiple samples
- Export to .sam header + individual sample files

Based on the SAM_Header workflow by the Swanson-Hysell Group.

Copyright (C) 2025 Python port, GNU GPL v.3
"""

import os
import sys
import numpy as np
from datetime import datetime
from typing import Dict, List, Optional, Tuple


class SampleOrientation:
    """Container for sample orientation data."""
    
    def __init__(self, sample_name: str):
        self.sample_name = sample_name
        self.comment = ""
        self.strat_level = 0.0
        
        # Core orientation
        self.core_strike = None  # Will be calculated or entered
        self.core_dip = None
        self.magnetic_core_strike = None  # From magnetic compass
        
        # Bedding
        self.bedding_strike = 0.0
        self.bedding_dip = 0.0
        self.correct_bedding = True
        
        # Physical properties
        self.mass = 1.0
        self.volume = 1.0
        
        # Sun compass data (optional)
        self.shadow_angle = None
        self.gmt_offset = None
        self.datetime = None  # datetime object
        
        # Calculated fields
        self.sun_core_strike = None
        self.igrf_dec = None
        self.calculated_mag_dec = None
        self.corrected_bedding_strike = None


class SAMFile:
    """Generate .sam header files and sample files."""
    
    def __init__(self, site_id: str, site_name: str = "", 
                 lat: float = 0.0, lon: float = 0.0, elevation: float = 0.0):
        """
        Initialize SAM file generator.
        
        Parameters
        ----------
        site_id : str
            Site identifier (max 8 characters for 8.3 filename format)
        site_name : str
            Descriptive site name
        lat : float
            Latitude in decimal degrees North
        lon : float
            Longitude in decimal degrees East
        elevation : float
            Elevation in meters above sea level
        """
        self.site_id = site_id[:8]  # 8.3 format: 8 char max
        self.site_name = site_name or site_id
        self.lat = lat
        self.lon = lon
        self.elevation = elevation
        self.samples: List[SampleOrientation] = []
    
    def add_sample(self, sample: SampleOrientation):
        """Add a sample to this site."""
        self.samples.append(sample)
    
    def calculate_sun_orientations(self, verbose: bool = True):
        """
        Calculate sun compass orientations for all samples with sun data.
        
        Uses solar position algorithm to determine true azimuth from shadow angle.
        Requires: shadow_angle, datetime, GMT offset, site lat/lon.
        """
        from rmg_sam_utilities import sundec
        
        for sample in self.samples:
            if sample.shadow_angle is not None and sample.datetime is not None:
                try:
                    sundata = {
                        'shadow_angle': sample.shadow_angle,
                        'lat': self.lat,
                        'lon': self.lon,
                        'delta_u': sample.gmt_offset or 0,
                        'date': sample.datetime.strftime('%Y:%m:%d:%H:%M')
                    }
                    sample.sun_core_strike = sundec(sundata)
                    
                    if verbose:
                        print(f"  {sample.sample_name}: Sun compass → {sample.sun_core_strike:.1f}°")
                        
                except Exception as e:
                    if verbose:
                        print(f"  {sample.sample_name}: Sun compass calculation failed - {e}")
    
    def calculate_igrf(self, verbose: bool = True):
        """
        Calculate IGRF magnetic declination for all samples.
        
        Uses IGRF-13 spherical harmonic model for accurate magnetic field.
        Requires: datetime, site lat/lon/elevation.
        """
        from rmg_sam_utilities import igrf, to_year_fraction
        
        for sample in self.samples:
            if sample.datetime is not None:
                try:
                    date_frac = to_year_fraction(sample.datetime)
                    alt_km = self.elevation / 1000.0
                    
                    igrf_result = igrf([date_frac, alt_km, self.lat, self.lon])
                    dec = igrf_result[0]
                    
                    # Convert to -180 to +180 range
                    if dec > 180:
                        dec = dec - 360
                    
                    sample.igrf_dec = dec
                    
                    if verbose:
                        print(f"  {sample.sample_name}: IGRF declination = {dec:+.1f}°")
                        
                except Exception as e:
                    if verbose:
                        print(f"  {sample.sample_name}: IGRF calculation failed - {e}")
    
    def validate_orientations(self, verbose: bool = True):
        """
        Compare sun compass vs magnetic compass orientations.
        
        Calculates local magnetic declination from the difference and warns
        if it differs from IGRF by more than 5°.
        """
        for sample in self.samples:
            if (sample.sun_core_strike is not None and 
                sample.magnetic_core_strike is not None):
                
                # Calculate local declination
                calc_dec = sample.sun_core_strike - sample.magnetic_core_strike
                
                # Normalize to -180 to +180
                if calc_dec > 180:
                    calc_dec -= 360
                elif calc_dec < -180:
                    calc_dec += 360
                
                sample.calculated_mag_dec = calc_dec
                
                if verbose:
                    print(f"  {sample.sample_name}: Calculated mag dec = {calc_dec:+.1f}°")
                    
                    # Compare with IGRF
                    if sample.igrf_dec is not None:
                        diff = abs(sample.igrf_dec - calc_dec)
                        if diff > 5:
                            print(f"    ⚠ WARNING: IGRF ({sample.igrf_dec:+.1f}°) vs "
                                  f"calculated ({calc_dec:+.1f}°) differ by {diff:.1f}°")
    
    def determine_final_orientations(self, prefer_sun: bool = True):
        """
        Determine final core orientations to use.
        
        Priority:
        1. Sun compass (most accurate) if prefer_sun=True
        2. Magnetic compass + IGRF correction
        
        Parameters
        ----------
        prefer_sun : bool
            If True, use sun compass when available. Otherwise use magnetic.
        """
        for sample in self.samples:
            # Determine core strike
            if prefer_sun and sample.sun_core_strike is not None:
                sample.core_strike = sample.sun_core_strike
                sample.comment = "sun compass orientation"
            elif sample.magnetic_core_strike is not None:
                if sample.igrf_dec is not None:
                    # Correct magnetic compass with IGRF
                    sample.core_strike = sample.magnetic_core_strike + sample.igrf_dec
                    sample.comment = "mag compass orientation (IGRF corrected)"
                else:
                    sample.core_strike = sample.magnetic_core_strike
                    sample.comment = "mag compass orientation (uncorrected)"
            else:
                # No orientation data available
                sample.core_strike = 0.0
                sample.comment = "no orientation data"
            
            # Normalize to 0-360
            if sample.core_strike is not None:
                sample.core_strike = sample.core_strike % 360
            
            # Correct bedding strike if requested
            if sample.correct_bedding and sample.igrf_dec is not None:
                sample.corrected_bedding_strike = (sample.bedding_strike + sample.igrf_dec) % 360
            else:
                sample.corrected_bedding_strike = sample.bedding_strike
    
    def write_sam_file(self, output_path: str):
        """
        Write .sam header file.
        
        Format (CIT):
        Line 1: CIT
        Line 2: Site name (comment)
        Line 3: lat lon mag_dec
        Line 4+: sample filenames
        
        Parameters
        ----------
        output_path : str
            Path for output .sam file
        """
        lines = []
        
        # Line 1: Format identifier
        lines.append("CIT")
        
        # Line 2: Site name
        lines.append(self.site_name)
        
        # Line 3: Location info
        # Format: " 36.2 245.3   0.0"
        # Mag dec set to 0.0 because orientations are pre-corrected
        lat_str = f"{self.lat:5.1f}"
        lon_str = f"{self.lon % 360:05.1f}"  # Force 0-360 range
        lines.append(f" {lat_str} {lon_str}   0.0")
        
        # Line 4+: Sample names
        for sample in self.samples:
            # Check if sample_name already includes site_id
            if sample.sample_name.startswith(self.site_id):
                sample_id = sample.sample_name
            else:
                sample_id = self.site_id + str(sample.sample_name)
            lines.append(sample_id)
        
        # Write file with Windows line endings (RAPID software requirement)
        with open(output_path, 'w') as f:
            f.write('\r\n'.join(lines) + '\r\n')
        
        print(f"\n✓ Created .sam file: {output_path}")
    
    def write_sample_files(self, output_dir: str):
        """
        Write individual sample files.
        
        Format (CIT):
        Line 1: site_id sample_name comment
        Line 2: strat_level core_strike core_dip bed_strike bed_dip mass
        Line 3+: Measurement data (added later by magnetometer)
        
        Parameters
        ----------
        output_dir : str
            Directory for output files
        """
        for sample in self.samples:
            # Check if sample_name already includes site_id
            if sample.sample_name.startswith(self.site_id):
                sample_id = sample.sample_name
            else:
                sample_id = self.site_id + str(sample.sample_name)
            filepath = os.path.join(output_dir, sample_id)
            
            lines = []
            
            # Line 1: sample ID and comment
            # Format: "site 1.0A Sample comment"
            # If sample_name already includes site_id, just use site_id + space + full_name
            # Otherwise use site_id + space + sample_name
            if sample.sample_name.startswith(self.site_id):
                # Sample name is full (e.g., "ObsFORC1")
                # Split to get just the number part
                name_only = sample.sample_name[len(self.site_id):]
                line1 = f"{self.site_id} {name_only}"
            else:
                line1 = f"{self.site_id} {sample.sample_name}"
            if sample.comment:
                line1 += f" {sample.comment}"
            lines.append(line1)
            
            # Line 2: Orientation data
            # Format: " 0.0    113.0  27.0  201.0  46.0   1.0"
            #         strat  core_s c_dip bed_s  bed_d  mass
            bed_strike = sample.corrected_bedding_strike or sample.bedding_strike
            
            line2 = (f" {sample.strat_level:6.1f}"
                    f"    {sample.core_strike or 0:5.1f}"
                    f" {sample.core_dip or 0:5.1f}"
                    f" {bed_strike:5.1f}"
                    f" {sample.bedding_dip:5.1f}"
                    f"   {sample.mass:.1f}")
            lines.append(line2)
            
            # Write file
            with open(filepath, 'w') as f:
                f.write('\r\n'.join(lines) + '\r\n')
        
        print(f"✓ Created {len(self.samples)} sample files in {output_dir}")


def generate_blank_sam_template(site_id: str = None, num_samples: int = None,
                                output_dir: str = None, interactive: bool = True):
    """
    Generate a blank .sam file template with empty sample files.
    
    Creates a template that can be manually filled in for sites where
    no automated calculations are needed.
    
    Parameters
    ----------
    site_id : str
        Site identifier (max 8 characters)
    num_samples : int
        Number of blank sample entries to create
    output_dir : str
        Output directory
    interactive : bool
        If True, prompt user for all details
    
    Returns
    -------
    str : Path to created .sam file
    """
    import os
    from rmg_sam_naming import (get_sample_naming_interactive, 
                                generate_sample_names, 
                                preview_sample_names)
    
    if interactive:
        print("="*70)
        print("BLANK SAM TEMPLATE GENERATOR")
        print("="*70)
        
        # Site ID (max 8 chars for 8.3 format)
        site_id = input("\nSite/Folder name (max 8 chars): ").strip()[:8]
        if not site_id:
            print("Error: Site name required")
            return None
        
        if len(site_id) > 8:
            print(f"Warning: Site name truncated to 8 characters: '{site_id}'")
        
        # Site name/description
        site_name = input("Site description [blank]: ").strip() or f"{site_id} Site"
        
        # Location - ask if user wants to set coordinates
        set_location = input("\nSet site coordinates? (y/n) [n]: ").strip().lower() == 'y'
        if set_location:
            try:
                lat = float(input("Latitude (decimal degrees N): ").strip())
                lon = float(input("Longitude (decimal degrees E): ").strip())
                elev = float(input("Elevation (meters) [0]: ").strip() or "0")
            except ValueError:
                print("Invalid coordinates, using 0,0")
                lat, lon, elev = 0.0, 0.0, 0.0
        else:
            lat, lon, elev = 0.0, 0.0, 0.0
        
        # Number of samples
        try:
            num_samples = int(input("\nNumber of samples to create: ").strip())
        except ValueError:
            print("Invalid number, using 10")
            num_samples = 10
        
        # Get sample naming preferences
        naming_prefs = get_sample_naming_interactive()
        
        # Preview sample names and validate
        while True:
            if preview_sample_names(site_id, num_samples, naming_prefs):
                confirm = input("\nAccept these sample names? (y/n) [y]: ").strip()
                if confirm.lower() != 'n':
                    break
            
            print("\nWould you like to:")
            print("  [1] Change naming preferences")
            print("  [2] Use shorter site name")
            print("  [3] Continue anyway")
            choice = input("Choice [1]: ").strip() or "1"
            
            if choice == "1":
                naming_prefs = get_sample_naming_interactive()
            elif choice == "2":
                site_id = input("New site name (max 8 chars): ").strip()[:8]
            else:
                break
        
        # Default values for all samples
        print("\n--- Default Values for All Samples ---")
        
        # Core plate strike
        use_default_strike = input("Set core plate strike to 90° for all samples? (y/n) [y]: ").strip()
        if use_default_strike.lower() != 'n':
            default_core_strike = 90.0
        else:
            try:
                default_core_strike = float(input("Core plate strike (°) [0]: ").strip() or "0")
            except ValueError:
                default_core_strike = 0.0
        
        # Core dip
        try:
            default_core_dip = float(input("Core plate dip (°) [0]: ").strip() or "0")
        except ValueError:
            default_core_dip = 0.0
        
        # Bedding
        try:
            default_bedding_strike = float(input("Bedding strike (°) [0]: ").strip() or "0")
            default_bedding_dip = float(input("Bedding dip (°) [0]: ").strip() or "0")
        except ValueError:
            default_bedding_strike, default_bedding_dip = 0.0, 0.0
        
        # Mass
        try:
            default_mass = float(input("Sample mass (g) [1.0]: ").strip() or "1.0")
        except ValueError:
            default_mass = 1.0
        
        # Comment
        default_comment = input("Sample description/comment [blank template]: ").strip() or "blank template"
        
        # Output location
        output_dir = input("\nSave location [current directory]: ").strip() or "."
        
    else:
        # Non-interactive mode - use provided values or defaults
        if site_id is None:
            raise ValueError("site_id required in non-interactive mode")
        site_name = f"{site_id} Site"
        lat, lon, elev = 0.0, 0.0, 0.0
        num_samples = num_samples or 10
        naming_prefs = {
            'use_integers': True,
            'sequential': True,
            'start_number': 1,
            'extension': ''
        }
        default_core_strike = 90.0
        default_core_dip = 0.0
        default_bedding_strike = 0.0
        default_bedding_dip = 0.0
        default_mass = 1.0
        default_comment = "blank template"
        output_dir = output_dir or "."
    
    # Generate sample names
    sample_list = generate_sample_names(
        site_id=site_id,
        num_samples=num_samples,
        **naming_prefs
    )
    
    # Create SAM file
    sam = SAMFile(site_id=site_id, site_name=site_name, 
                  lat=lat, lon=lon, elevation=elev)
    
    # Add samples
    for full_name, is_valid, warning in sample_list:
        if not is_valid:
            print(f"Warning: {full_name} - {warning}")
        
        sample = SampleOrientation(full_name)
        sample.core_strike = default_core_strike
        sample.core_dip = default_core_dip
        sample.bedding_strike = default_bedding_strike
        sample.bedding_dip = default_bedding_dip
        sample.mass = default_mass
        sample.comment = default_comment
        sam.add_sample(sample)
    
    # Create output directory if needed
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Write files
    sam_path = os.path.join(output_dir, f"{site_id}.sam")
    sam.write_sam_file(sam_path)
    sam.write_sample_files(output_dir)
    
    print(f"\n{'='*70}")
    print("BLANK TEMPLATE CREATED")
    print(f"{'='*70}")
    print(f"✓ Site: {site_id}")
    print(f"✓ Samples: {num_samples}")
    print(f"✓ Location: {output_dir}")
    print(f"✓ .sam file: {sam_path}")
    if lat != 0 or lon != 0:
        print(f"✓ Coordinates: {lat:.4f}°N, {lon:.4f}°E")
    else:
        print(f"✓ Coordinates: Not set (0, 0)")
    
    print(f"\nDefault values:")
    print(f"  Core strike: {default_core_strike}°")
    print(f"  Core dip: {default_core_dip}°")
    print(f"  Bedding: {default_bedding_strike}°/{default_bedding_dip}°")
    print(f"  Mass: {default_mass} g")
    print(f"  Comment: {default_comment}")
    
    # Show naming info
    if interactive:
        ext_info = f" with .{naming_prefs['extension']}" if naming_prefs['extension'] else " (no extension)"
        pattern = "Sequential" if naming_prefs['sequential'] else "Sub-samples"
        num_type = "integers" if naming_prefs['use_integers'] else "decimals"
        print(f"\nSample naming:")
        print(f"  Pattern: {pattern} using {num_type}{ext_info}")
    
    print(f"\nEdit individual sample files to modify values.")
    
    return sam_path


def create_sam_interactive():
    """
    Interactive SAM file creation with user prompts.
    
    Guides the user through creating a complete SAM file with
    optional sun compass and IGRF calculations.
    """
    print("="*70)
    print("INTERACTIVE SAM FILE GENERATOR")
    print("="*70)
    
    # Site information
    print("\n--- Site Information ---")
    site_id = input("Site ID (max 4 chars): ").strip()[:4]
    site_name = input("Site name [optional]: ").strip() or site_id
    
    try:
        lat = float(input("Latitude (decimal degrees N): ").strip())
        lon = float(input("Longitude (decimal degrees E): ").strip())
        elev = float(input("Elevation (meters) [0]: ").strip() or "0")
    except ValueError:
        print("Invalid number, using defaults")
        lat, lon, elev = 0.0, 0.0, 0.0
    
    sam = SAMFile(site_id, site_name, lat, lon, elev)
    
    # Sample information
    print("\n--- Sample Information ---")
    try:
        num_samples = int(input("Number of samples: ").strip())
    except ValueError:
        print("Invalid number, using 1 sample")
        num_samples = 1
    
    use_sun = input("Use sun compass data? (y/n) [n]: ").strip().lower() == 'y'
    use_mag = input("Use magnetic compass data? (y/n) [n]: ").strip().lower() == 'y'
    
    for i in range(num_samples):
        print(f"\n--- Sample {i+1} ---")
        sample_name = input(f"Sample name [{i+1}.0a]: ").strip() or f"{i+1}.0a"
        sample = SampleOrientation(sample_name)
        
        # Core orientation
        if use_mag:
            try:
                sample.magnetic_core_strike = float(input("Magnetic core strike (°): ").strip())
                sample.core_dip = float(input("Core dip (°): ").strip())
            except ValueError:
                print("Invalid number, skipping")
        
        # Bedding
        try:
            sample.bedding_strike = float(input("Bedding strike (°) [0]: ").strip() or "0")
            sample.bedding_dip = float(input("Bedding dip (°) [0]: ").strip() or "0")
        except ValueError:
            pass
        
        # Sun compass
        if use_sun:
            try:
                sample.shadow_angle = float(input("Shadow angle (°): ").strip())
                year = int(input("Year (YYYY): ").strip())
                month = int(input("Month (1-12): ").strip())
                day = int(input("Day: ").strip())
                hour = int(input("Hour (0-23): ").strip())
                minute = int(input("Minute: ").strip() or "0")
                sample.datetime = datetime(year, month, day, hour, minute)
                sample.gmt_offset = float(input("GMT offset (hrs to subtract): ").strip())
            except ValueError:
                print("Invalid date/time, skipping sun compass")
        
        # Physical
        try:
            sample.mass = float(input("Mass (g) [1.0]: ").strip() or "1.0")
        except ValueError:
            pass
        
        sam.add_sample(sample)
    
    # Calculations
    print("\n--- Processing ---")
    if use_sun:
        print("Calculating sun compass orientations...")
        sam.calculate_sun_orientations()
    
    print("Calculating IGRF magnetic declination...")
    sam.calculate_igrf()
    
    if use_sun and use_mag:
        print("Validating orientations...")
        sam.validate_orientations()
    
    print("Determining final orientations...")
    sam.determine_final_orientations(prefer_sun=use_sun)
    
    # Output
    print("\n--- Output ---")
    output_dir = input("Output directory [.]: ").strip() or "."
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    sam_path = os.path.join(output_dir, f"{site_id}.sam")
    sam.write_sam_file(sam_path)
    sam.write_sample_files(output_dir)
    
    print(f"\n{'='*70}")
    print("SAM FILE GENERATION COMPLETE")
    print(f"{'='*70}")
    print(f"Site: {site_id}")
    print(f"Samples: {len(sam.samples)}")
    print(f"Output: {output_dir}")


if __name__ == "__main__":
    if len(sys.argv) > 1:
        if sys.argv[1] == '--blank':
            site_id = input("Site ID: ").strip()
            num = int(input("Number of samples [10]: ").strip() or "10")
            output = input("Output directory [.]: ").strip() or "."
            generate_blank_sam_template(site_id, num, output)
        else:
            print("Usage: python rmg_sam.py [--blank]")
            print("  --blank : Generate blank template")
            print("  (no args) : Interactive mode")
    else:
        create_sam_interactive()
