# pyRockmagCIT v1.0 - Complete File List

## Files to Upload to GitHub Repository

### Main Application
- `run_pyrockmag.py` - **Main application entry point** (renamed from run_rmg_analysis.py)

### Core Modules
- `rmg_import.py` - Data import and .rmg file parsing
- `rmg_stats.py` - Statistical analysis and coercivity calculations
- `rmg_plots.py` - Plotting routines for all rock magnetic data
- `rmg_forc.py` - FORC diagram processing and visualization
- `rmg_sam.py` - SAM header file generation
- `rmg_sam_utilities.py` - Sun compass and IGRF calculations
- `rmg_sam_naming.py` - 8.3 filename validation utilities

### Additional Utility Modules (from original matRockmag)
- `RmgARMCurveSubtract.m` → Python equivalent in rmg_plots.py
- `RmgExtractAFOfStep.m` → Python equivalent in rmg_import.py
- `fitSGG.m` → Python equivalent in rmg_forc.py
- `fitSGGComps.m` → Python equivalent in rmg_forc.py
- `SGG.m` → Python equivalent in rmg_forc.py
- `RmgStatsWriteTable.m` → Python equivalent in rmg_stats.py
- `RmgStats.m` → Python equivalent in rmg_stats.py
- `RmgStatDepthProfiles.m` → Python equivalent in rmg_plots.py
- `RmgStatBox.m` → Python equivalent in rmg_plots.py
- (And all other modules from the matRockmag conversion)

### Documentation
- `README.md` - **Main documentation** (comprehensive overview)
- `FEATURES.md` - Detailed feature list and capabilities
- `CHANGELOG.md` - Version history and release notes
- `INSTALL.md` - Installation guide
- `LICENSE` - GNU GPL v3.0 license

### SAM/Paleomag Documentation
- `SAM_Generation_Guide.md` - Complete SAM file generation tutorial
- `SAM_Blank_Template_Guide.md` - Field work template guide
- `SAM_File_Generation_Comparison.md` - Excel vs Python comparison

### Reference Materials
- `SAM_file_format.pdf` - CIT format specification (included in uploads)
- `Instructions_for_block_sample_orientations.md` - Block sample notes

### Configuration
- `requirements.txt` - Python package dependencies
- `.gitignore` - Git ignore patterns (recommended)

---

## Recommended .gitignore

Create a `.gitignore` file with:

```
# Python
__pycache__/
*.py[cod]
*$py.class
*.so
.Python
env/
venv/
ENV/
*.egg-info/
dist/
build/

# Data files (don't commit user data)
*.rmg
*.frc
*.sam
*.csv
*.dat

# Output directories
output/
outputs/
figures/
plots/
temp/
tmp/

# OS files
.DS_Store
Thumbs.db
desktop.ini

# IDE
.vscode/
.idea/
*.swp
*.swo
*~

# Testing
.pytest_cache/
.coverage
htmlcov/

# Documentation builds
docs/_build/
```

---

## Repository Structure

Recommended structure for GitHub:

```
pyRockmagCIT/
├── run_pyrockmag.py              # Main entry point
├── rmg_*.py                      # All core modules
├── README.md                     # Main documentation
├── FEATURES.md                   # Feature overview
├── CHANGELOG.md                  # Version history
├── INSTALL.md                    # Installation guide
├── LICENSE                       # GPL v3.0
├── requirements.txt              # Dependencies
├── .gitignore                    # Git ignore rules
│
├── docs/                         # Documentation folder
│   ├── SAM_Generation_Guide.md
│   ├── SAM_Blank_Template_Guide.md
│   ├── SAM_File_Generation_Comparison.md
│   ├── Instructions_for_block_sample_orientations.md
│   └── SAM_file_format.pdf
│
├── examples/                     # Example data (optional)
│   ├── sample_data/
│   │   ├── example.rmg
│   │   └── example.frc
│   └── tutorials/
│       └── getting_started.md
│
└── tests/                        # Unit tests (optional, future)
    ├── test_import.py
    ├── test_stats.py
    └── test_forc.py
```

---

## GitHub Repository Setup

### 1. Create Repository
```bash
# On GitHub website:
# - Go to https://github.com/yicolas
# - Click "New repository"
# - Name: pyRockmagCIT
# - Description: Comprehensive Rock Magnetic & Paleomagnetic Analysis Suite
# - Public or Private (your choice)
# - Add README: NO (you already have one)
# - Add .gitignore: Python
# - Choose license: GPL-3.0
```

### 2. Initialize Local Repository
```bash
cd /path/to/pyRockmagCIT
git init
git add .
git commit -m "Initial commit - pyRockmagCIT v1.0"
```

### 3. Connect to GitHub
```bash
git remote add origin https://github.com/yicolas/pyRockmagCIT.git
git branch -M main
git push -u origin main
```

### 4. Add Topics (GitHub website)
Suggested topics for discoverability:
- `rock-magnetism`
- `paleomagnetism`
- `forc-diagrams`
- `magnetometer`
- `python`
- `geophysics`
- `earth-science`
- `2g-enterprises`
- `cit-format`

---

## File Checklist

Upload all files from `/mnt/user-data/outputs/rmg_python/`:

### Essential Files (Must Have)
- [ ] run_pyrockmag.py
- [ ] rmg_import.py
- [ ] rmg_stats.py
- [ ] rmg_plots.py
- [ ] rmg_forc.py
- [ ] rmg_sam.py
- [ ] rmg_sam_utilities.py
- [ ] rmg_sam_naming.py
- [ ] README.md
- [ ] LICENSE
- [ ] requirements.txt

### Documentation (Highly Recommended)
- [ ] FEATURES.md
- [ ] CHANGELOG.md
- [ ] INSTALL.md
- [ ] SAM_Generation_Guide.md
- [ ] SAM_Blank_Template_Guide.md
- [ ] SAM_File_Generation_Comparison.md

### Optional (But Useful)
- [ ] .gitignore
- [ ] Examples folder with sample data
- [ ] Tests folder (for future development)

---

## Release Notes for v1.0.0

When creating the first release on GitHub:

**Tag:** `v1.0.0`
**Title:** pyRockmagCIT v1.0.0 - Initial Release
**Description:**
```
# pyRockmagCIT v1.0.0 - Initial Release

Complete rock magnetic and paleomagnetic analysis suite in Python.

## Features
- 🪨 Rock magnetic analysis (.rmg files)
- 📊 FORC diagram processing (FORCinel v3.0.8)
- 🧭 SAM header file generation (CIT format)
- ⚙️ FORC protocol generation

## Installation
```bash
pip install numpy scipy matplotlib pandas
git clone https://github.com/yicolas/pyRockmagCIT.git
cd pyRockmagCIT
python run_pyrockmag.py
```

See INSTALL.md for detailed instructions.

## Documentation
- README.md - Overview
- FEATURES.md - Complete feature list
- SAM_Generation_Guide.md - Paleomag workflow

## License
GNU GPL v3.0
```

---

## Contact Information

- **GitHub:** https://github.com/yicolas/pyRockmagCIT
- **Issues:** https://github.com/yicolas/pyRockmagCIT/issues
- **Email:** yick@duck.com

---

All files are ready in `/mnt/user-data/outputs/rmg_python/` and `/mnt/user-data/outputs/`
