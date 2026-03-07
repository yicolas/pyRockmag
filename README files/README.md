# pyRockmagCIT v1.0

**Comprehensive Rock Magnetic & Paleomagnetic Analysis Suite**

Python toolkit for rock magnetic data analysis, FORC diagram processing, and paleomagnetic sample orientation management. Designed for the CIT (Caltech) data format used with 2G Enterprises magnetometers.

---

## Overview

**pyRockmagCIT** is a complete Python port and expansion of the MATLAB `matRockmag` toolbox, with additional features for FORC analysis and SAM file generation. It provides a unified workflow from field data collection through laboratory analysis to publication-ready figures.

### Key Features

#### 🪨 **Rock Magnetic Analysis**
- Import and parse .rmg files from 2G magnetometers
- Statistical analysis of coercivity parameters (Hcr, Hc, MDF, SIRM)
- Automated plotting routines (hysteresis, IRM, ARM, backfield)
- Batch processing with grouped analyses
- Export statistics to tables

#### 📊 **FORC Diagram Processing**
- Complete FORCinel v3.0.8 workflow integration
- Multiple visualization styles (scatter, contour, multiple colormaps)
- Smoothing and interpolation algorithms
- FORC distribution calculations with coordinate transforms
- Export to generic .frc format
- Complete hysteresis plotting with reversal points

#### 🧭 **Paleomagnetic Sample Orientation**
- Generate .sam header files (CIT format)
- Sun compass orientation calculations
- IGRF-13 magnetic declination modeling
- Automatic validation (sun vs magnetic compass)
- 8.3 filename format compliance
- Blank template generation for field work

#### ⚙️ **Protocol Generation**
- Create FORC measurement scripts for 2G magnetometers
- Custom field ranges and step sizes
- Multiple FORC protocols (standard, exponential, 3-axis)

---

## Installation

### Requirements
- Python 3.8+
- NumPy
- SciPy
- Matplotlib
- Pandas

### Setup

```bash
# Clone or download the repository
cd pyRockmagCIT

# Install dependencies
pip install numpy scipy matplotlib pandas

# Run the application
python run_pyrockmag.py
```

---

## Quick Start

### Interactive Mode

```bash
python run_pyrockmag.py
```

**Main Menu:**
```
[L] Load .rmg files from directory
[R] Run batch analysis routines
[F] Full 3×3 dashboard analysis
[P] Plot FORC data
[G] Generate FORC protocol
[S] Generate SAM header files
[E] Export statistics
[Q] Quit
```

### Command Line

```bash
# Load files directly
python run_pyrockmag.py sample1.rmg sample2.rmg

# From directory
python run_pyrockmag.py /path/to/data/*.rmg
```

---

## Capabilities

### 1. Rock Magnetic Data Analysis

**Import .rmg Files:**
- Automatic parsing of 2G magnetometer data
- Support for multiple measurement types
- Replicate handling
- Quality control checks

**Statistical Analysis:**
- Coercivity of remanence (Hcr)
- Coercive force (Hc)
- Median destructive field (MDF)
- Saturation remanence (SIRM)
- S-ratios and ARM parameters

**Visualization:**
- Hysteresis loops
- IRM acquisition curves
- ARM curves
- Backfield demagnetization
- Lowrie-Fuller tests
- 

**Batch Processing:**
- Process multiple samples simultaneously
- Group by site, lithology, or custom criteria
- Generate comparison plots
- Export summary tables

### 2. FORC Diagram Analysis

**Data Processing:**
- Extract FORC curves from .rmg files
- Interpolate to regular grids (100×100)
- Apply smoothing (FORCinel algorithm)
- Calculate FORC distribution: ρ = -0.5 × ∂²M/∂Hr∂Ha
- Transform to (Bc, Bu) coercivity coordinates

**Visualization Options:**
- **Scatter plots**: hot_r, plasma_r, inferno_r colormaps
- **Contour plots**: Smooth interpolated surfaces
- **Hysteresis plots**: Complete curves with reversal points
- **Multiple outputs**: Generate all visualization types

**Export Formats:**
- Generic .frc format (FORCinel v3.0.6+ compatible)
- PNG/SVG publication-quality figures
- CSV data tables

**Key Features:**
- Positive colorbar values (strong signal = yellow)
- Clipped data to prevent interpolation artifacts
- Customizable smoothing factors (2-4 recommended)
- Full workflow automation

### 3. SAM Header File Generation

**Blank Templates:**
- Quick generation for field work
- User-defined defaults for all parameters
- 8.3 filename format compliance
- Core strike default: 90° (standard)
- Optional coordinate entry

**Interactive Creation with Calculations:**
- Sun compass orientation (solar position algorithm)
- IGRF-13 magnetic declination
- Altitude corrections
- Secular variation included
- Automatic validation (sun vs magnetic)

**Sample Naming:**
- Flexible numbering schemes:
  - Sequential (1, 2, 3, ...)
  - Sub-samples (1.1, 1.2, 1.3, ...)
  - Integers or decimals
- Custom extensions (.FRC, .0a, etc.)
- 8-character site names
- Full 8.3 format validation

**Coordinate Conventions:**
- Core plate strike/dip (perpendicular to core axis)
- Right-hand rule for bedding
- Automatic IGRF corrections
- Block sample compass support (clockwise → negate angle)

### 4. Protocol Generation

**FORC Measurement Scripts:**
- Standard FORC (linear field steps)
- Exponential FORC (log-spaced fields)
- 3-axis FORC
- Custom field ranges
- Configurable step sizes

**Output:**
- .rmg script format
- Copy-paste ready for 2G software
- Optimized measurement sequences

---

## Workflow Examples

### Rock Magnetic Analysis

```bash
python run_pyrockmag.py

[L] Load files → select directory
[R] Run routines → choose plots (hysteresis, IRM, ARM)
[E] Export statistics → save table

# Output: Plots + statistics.csv
```

### FORC Analysis

```bash
python run_pyrockmag.py

[L] Load FORC files
[P] Plot FORC → [6] All plots
# Generates: hysteresis + 3 scatter + 3 contour plots

[P] Plot FORC → [7] FORCinel workflow
# Smoothing factor: 3
# Exports: .frc file + processed diagrams
```

### Field Data Preparation

```bash
python run_pyrockmag.py

[S] Generate SAM → [1] Blank template
Site: ObsFORC (8 chars max)
Coordinates: n (skip for blank)
Samples: 20
Naming: [1] Integers, [1] Sequential
Extension: FRC
Core strike 90°: y
Mass: 10.5 g

# Output: ObsFORC.sam + 20 blank sample files
# Take to field, fill in orientations
```

### Lab Orientation Processing

```bash
python run_pyrockmag.py

[S] Generate SAM → [2] Interactive
# Enter all field data
# Sun compass + magnetic compass
# System calculates orientations
# Validates (warns if >5° difference)

# Output: .sam + samples with calculated orientations
# Ready for magnetometer!
```

---

## File Formats

### Input Formats

**.rmg Files:**
- 2G Enterprises magnetometer format
- ASCII text with structured headers
- Multiple measurement types supported
- Replicate measurements handled

### Output Formats

**.sam Files (CIT):**
```
CIT
Site Name
 lat lon 0.0
sample1
sample2
...
```

**Sample Files:**
```
site 1 comment
 strat  strike  dip  bed_s  bed_d  mass
```

**.frc Files (Generic FORC):**
```
Ha,M
Ha,M
...
[blank line]
Ha,M
...
END
```

---

## Module Reference

### Core Modules

**`rmg_import.py`**
- Parse .rmg files
- Extract measurement data
- Handle replicates

**`rmg_stats.py`**
- Coercivity calculations
- Statistical analysis
- Parameter extraction

**`rmg_plots.py`**
- All plotting routines
- Batch plot generation
- Publication-quality figures

**`rmg_forc.py`**
- FORC data extraction
- Distribution calculations
- Multiple visualization styles
- FORCinel workflow

**`rmg_sam.py`**
- SAM file generation
- Sample orientation calculations
- Interactive and programmatic modes

**`rmg_sam_utilities.py`**
- Sun compass calculations
- IGRF magnetic field modeling
- Coordinate transformations

**`rmg_sam_naming.py`**
- 8.3 filename validation
- Sample name generation
- Format compliance checking

**`run_pyrockmag.py`**
- Main interface
- Menu system
- Workflow orchestration

---

## Documentation

### Guides
- **SAM_Generation_Guide.md** - Complete SAM file tutorial
- **SAM_Blank_Template_Guide.md** - Field work preparation
- **SAM_File_Generation_Comparison.md** - Excel vs Python
- **FORC_Analysis_Guide.md** - FORC processing details

### References
- **SAM_file_format.pdf** - CIT format specification
- **How_to_analyze_remFORC_data_v3.docx** - FORCinel workflow

---

## Citations

### Original Software

**matRockmag (MATLAB):**
- Original MATLAB toolbox for rock magnetic analysis
- Developed at Caltech
- CIT data format standard

**FORCinel:**
- Harrison, R.J., and Feinberg, J.M. (2008). FORCinel: An improved algorithm for calculating first-order reversal curve distributions using locally weighted regression smoothing. *Geochemistry, Geophysics, Geosystems*, 9, Q05016.

**SAM_Header (Python):**
- Swanson-Hysell Group, UC Berkeley
- Original Python implementation of SAM file generation
- Sun compass and IGRF calculations

**IGRF Model:**
- Alken, P., et al. (2021). International Geomagnetic Reference Field: the thirteenth generation. *Earth, Planets and Space*, 73, 49.

### This Software

**pyRockmagCIT v1.0:**
```
Anderson, N. (2025). pyRockmagCIT: Comprehensive rock magnetic and 
paleomagnetic analysis suite. Python software. 
https://github.com/yicolas/pyRockmagCIT
```

**Key Contributions:**
- Complete Python port of matRockmag
- FORCinel v3.0.8 workflow integration
- Enhanced FORC visualization with multiple colormaps
- 8.3 filename format compliance for SAM files
- Unified interface for rock magnetic and paleomagnetic workflows

---

## Changelog

### v1.0.0 (2025-02-20)
- **Initial release**
- Complete port of matRockmag MATLAB toolbox
- FORCinel v3.0.8 workflow integration
- SAM header file generation with sun compass support
- IGRF-13 magnetic declination calculations
- 8.3 filename format support
- Interactive and batch processing modes
- Publication-quality plotting with multiple colormaps
- Comprehensive documentation

---

## Contributing

Contributions welcome! Areas for development:

- Additional plot types
- More statistical methods
- Web interface
- Database integration
- Additional file format support
- Extended IGRF models

---

## License

GNU GPL v3.0

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

---

## Support

For questions, bug reports, or feature requests:
- GitHub Issues: yicolas/issues
- Email: yick@duck.com

---

## Acknowledgments

- **Caltech Paleomagnetics Lab** - matRockmag MATLAB toolbox
- **Swanson-Hysell Group** - SAM_Header Python tools
- **Richard Harrison & Josh Feinberg** - FORCinel algorithm
- **IGRF Working Group** - Magnetic field modeling

---

**pyRockmagCIT** - From field to publication, one toolkit.
