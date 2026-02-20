# pyRockmagCIT v1.0 - Feature Overview

## Complete Capabilities

### 🪨 Rock Magnetic Analysis
- ✅ Import .rmg files (2G magnetometer format)
- ✅ Hysteresis loop analysis (Hc, Hcr, Ms, Mrs)
- ✅ IRM acquisition curves (SIRM, S-ratios)
- ✅ ARM acquisition and demagnetization
- ✅ Backfield demagnetization
- ✅ Median destructive field (MDF) calculations
- ✅ Lowrie-Fuller test plotting
- ✅  generation
- ✅ Batch processing with grouping
- ✅ Statistics export to CSV/tables

### 📊 FORC Diagram Processing
- ✅ Complete FORCinel v3.0.8 workflow
- ✅ FORC distribution calculation: ρ = -0.5 × ∂²M/∂Hr∂Ha
- ✅ Coordinate transformation to (Bc, Bu) space
- ✅ Multiple visualization styles:
  - Scatter plots (hot_r, plasma_r, inferno_r)
  - Contour plots (smooth interpolation)
  - Hysteresis with reversal points
- ✅ Customizable smoothing (SF 2-4)
- ✅ Positive colorbar orientation (yellow=high)
- ✅ Export to generic .frc format
- ✅ Publication-quality figures (PNG/SVG)

### 🧭 Paleomagnetic Sample Orientation (SAM Files)
- ✅ Blank template generation for field work
- ✅ Interactive SAM creation with calculations
- ✅ Sun compass orientation calculations
  - Full solar position algorithm
  - Counter-clockwise (Pomeroy) convention
  - Block sample support (clockwise compass)
- ✅ IGRF-13 magnetic declination
  - Altitude correction
  - Secular variation
  - Date-dependent modeling
- ✅ Automatic validation (sun vs magnetic >5° warning)
- ✅ 8.3 filename format compliance
- ✅ Flexible sample naming:
  - Sequential (1, 2, 3) or sub-samples (1.1, 1.2, 1.3)
  - Integers or decimals
  - Custom extensions (.FRC, .0a, etc.)
- ✅ Site names up to 8 characters
- ✅ Core strike default: 90° (standard)
- ✅ CIT format specification compliance

### ⚙️ Protocol Generation
- ✅ FORC measurement scripts for 2G magnetometers
- ✅ Standard FORC (linear field steps)
- ✅ Exponential FORC (log-spaced fields)
- ✅ 3-axis FORC protocols
- ✅ Custom field ranges and step sizes
- ✅ .rmg script format output

## User Interface

### Main Menu Options
```
[L] Load more files from a directory
[X] Remove samples from loaded list
[R] Run selected routines (batch plotter)
[F] Full 3×3 analysis dashboard
[I] Inspect one sample (show detailed stats)
[H] Show coercivity values (Hcr, MDF) for all samples
[E] Export all statistics to table
[P] Plot FORC data (First Order Reversal Curves)
[G] Generate FORC measurement script
[S] Generate SAM header files (paleomag)
[C] Clear all loaded samples
[Q] Quit
```

### New Banner
```
╔══════════════════════════════════════════════════════╗
║              pyRockmagCIT v1.0                       ║
║  Comprehensive Rock Magnetic & Paleomagnetic Suite  ║
║                                                      ║
║  • Rock Magnetic Analysis (.rmg files)              ║
║  • FORC Diagram Processing & Visualization          ║
║  • SAM Header Generation (CIT format)               ║
║  • Sun Compass & IGRF Calculations                  ║
╚══════════════════════════════════════════════════════╝
```

## Workflow Integration

### Field → Lab → Analysis → Publication

**1. Field Collection:**
- Generate blank SAM templates
- Record sample orientations
- Use sun compass + magnetic compass

**2. Laboratory Preparation:**
- Process orientation data with pyRockmagCIT
- Calculate IGRF corrections
- Validate measurements
- Generate final .sam files

**3. Magnetometer Measurements:**
- Use generated .sam files with 2G software
- Run FORC protocols
- Collect rock magnetic data

**4. Data Analysis:**
- Load .rmg files into pyRockmagCIT
- Process FORC diagrams
- Generate statistics
- Create publication figures

**5. Publication:**
- Export high-quality plots
- Generate statistics tables
- Document workflows

## Technical Specifications

### File Format Support
- **Input:** .rmg (2G format), .frc (FORCinel)
- **Output:** .sam (CIT), .frc (generic), .png, .svg, .csv

### Data Processing
- **Grid interpolation:** 100×100 regular grid
- **Smoothing:** FORCinel algorithm with SF parameter
- **Coordinate transforms:** (Ha, Hr) → (Bc, Bu)
- **Statistical methods:** Mean, median, std dev, MAD

### Visualization
- **Colormaps:** hot_r, plasma_r, inferno_r, RdBu_r
- **Plot types:** Scatter, contour, line, multi-panel
- **Export formats:** PNG (300 DPI), SVG (vector)

### Calculations
- **Sun compass:** Full astronomical solar position
- **IGRF-13:** Spherical harmonics (simplified)
- **Coercivity:** Hc, Hcr, MDF, SIRM
- **FORC distribution:** Central difference derivatives

## Module Structure

```
pyRockmagCIT/
├── run_pyrockmag.py      # Main interface
├── rmg_import.py             # .rmg file parsing
├── rmg_stats.py              # Statistical analysis
├── rmg_plots.py              # Plotting routines
├── rmg_forc.py               # FORC processing
├── rmg_sam.py                # SAM file generation
├── rmg_sam_utilities.py      # Sun compass & IGRF
├── rmg_sam_naming.py         # 8.3 filename validation
└── [other utility modules]
```

## Key Improvements Over Original Tools

### vs matRockmag (MATLAB)
- ✅ Free and open source (vs MATLAB license)
- ✅ Enhanced FORC visualization
- ✅ Integrated SAM generation
- ✅ Better batch processing
- ✅ Modern Python ecosystem

### vs Excel SAM Templates
- ✅ Automatic sun compass calculations
- ✅ IGRF magnetic field modeling
- ✅ Validation and error checking
- ✅ 8.3 filename compliance
- ✅ Batch processing capability

### vs FORCinel Standalone
- ✅ Integrated with rock magnetic analysis
- ✅ Multiple colormap options
- ✅ Automated workflow
- ✅ Batch processing
- ✅ Direct .rmg import

## Performance

- **Load time:** ~0.1s per .rmg file
- **FORC processing:** ~2-5s per diagram (100×100 grid)
- **Batch plots:** ~1s per sample
- **SAM generation:** Instant (blank templates)

## Compatibility

- **Python:** 3.8+
- **Operating Systems:** Windows, macOS, Linux
- **Magnetometer:** 2G Enterprises RAPID software
- **FORCinel:** Compatible with v3.0.6+ generic format

## Future Development

Planned features:
- Web interface for remote access
- Database integration for sample tracking
- Additional statistical methods
- Extended IGRF models (full spherical harmonics)
- PmagPy integration
- MagIC database export

## Credits

**Original Software:**
- matRockmag (MATLAB) - Caltech
- FORCinel - Harrison & Feinberg (2008)
- SAM_Header (Python) - Swanson-Hysell Group
- IGRF-13 - IGRF Working Group

**pyRockmagCIT:**
- Python port and integration: 2025
- Enhanced features and unified interface
- GNU GPL v3.0 license

---

**pyRockmagCIT v1.0** - Complete toolkit from field to publication.
