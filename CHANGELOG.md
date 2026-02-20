# Changelog

All notable changes to pyRockmagCIT will be documented in this file.

## [1.0.0] - 2025-02-20

### Added
- **Initial Release** of pyRockmagCIT
- Complete Python port of matRockmag MATLAB toolbox
- Rock magnetic analysis (.rmg file support)
  - Hysteresis loop analysis (Hc, Hcr, Ms, Mrs)
  - IRM acquisition and S-ratios
  - ARM curves
  - Backfield demagnetization
  - MDF calculations
  - Batch processing with grouping
  - Statistics export to CSV
- FORC diagram processing
  - FORCinel v3.0.8 workflow integration
  - Multiple colormaps (hot_r, plasma_r, inferno_r)
  - Scatter and contour plot options
  - Positive colorbar orientation (yellow=high)
  - Generic .frc format export
  - Smoothing factor customization
- SAM header file generation
  - Blank template generation for field work
  - Interactive creation with calculations
  - Sun compass orientation (full solar position algorithm)
  - IGRF-13 magnetic declination modeling
  - Automatic validation (sun vs magnetic compass)
  - 8.3 filename format compliance
  - Flexible sample naming schemes
  - Site names up to 8 characters
- FORC protocol generation
  - Standard FORC (linear fields)
  - Exponential FORC (log-spaced)
  - 3-axis FORC
  - Custom field ranges
- Interactive menu system
- Command-line file loading
- Comprehensive documentation

### Technical Details
- Python 3.8+ compatibility
- NumPy/SciPy/Matplotlib/Pandas dependencies
- Cross-platform (Windows, macOS, Linux)
- GNU GPL v3.0 license

### Known Issues
- IGRF implementation is simplified (full spherical harmonics planned for v1.1)
- Sun compass calculations assume standard Pomeroy convention (document block sample conversion)

---

## Future Releases (Planned)

### [1.1.0] - TBD
- Full IGRF-13 spherical harmonic implementation
- Web interface for remote access
- Database integration
- Additional statistical methods
- PmagPy integration
- MagIC database export

### [1.2.0] - TBD
- Real-time magnetometer data acquisition
- Live FORC plotting
- Machine learning mineral identification
- Automated report generation

---

## Version History

- **v1.0.0** (2025-02-20): Initial release

---

For detailed changes, see commit history at:
https://github.com/yicolas/pyRockmagCIT/commits
