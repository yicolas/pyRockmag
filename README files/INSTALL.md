# Installation Guide - pyRockmagCIT v1.0

## Quick Install

### Requirements
- Python 3.8 or higher
- pip (Python package installer)

### Install Dependencies

```bash
pip install numpy scipy matplotlib pandas
```

### Download pyRockmagCIT

```bash
git clone https://github.com/yicolas/pyRockmagCIT.git
cd pyRockmagCIT
```

### Run the Application

```bash
python run_pyrockmag.py
```

---

## Detailed Installation

### Step 1: Install Python

**Windows:**
1. Download Python from https://www.python.org/downloads/
2. Run installer
3. **CHECK** "Add Python to PATH"
4. Click "Install Now"

**macOS:**
```bash
# Using Homebrew
brew install python3
```

**Linux (Ubuntu/Debian):**
```bash
sudo apt update
sudo apt install python3 python3-pip
```

### Step 2: Verify Installation

```bash
python --version
# Should show: Python 3.8.x or higher

pip --version
# Should show pip version
```

### Step 3: Install Required Packages

**Method 1: Using pip directly**
```bash
pip install numpy scipy matplotlib pandas
```

**Method 2: Using requirements.txt** (if provided)
```bash
pip install -r requirements.txt
```

**Method 3: Using conda** (if you have Anaconda/Miniconda)
```bash
conda install numpy scipy matplotlib pandas
```

### Step 4: Download pyRockmagCIT

**Option A: Using git** (recommended)
```bash
git clone https://github.com/yicolas/pyRockmagCIT.git
cd pyRockmagCIT
```

**Option B: Download ZIP**
1. Go to https://github.com/yicolas/pyRockmagCIT
2. Click "Code" → "Download ZIP"
3. Extract to desired location
4. Open terminal in extracted folder

### Step 5: Test Installation

```bash
python run_pyrockmag.py --help
# Should show usage information
```

Or simply run:
```bash
python run_pyrockmag.py
```

You should see the banner:
```
╔══════════════════════════════════════════════════════╗
║              pyRockmagCIT v1.0                       ║
║  Comprehensive Rock Magnetic & Paleomagnetic Suite  ║
╚══════════════════════════════════════════════════════╝
```

---

## Platform-Specific Notes

### Windows

**Using Anaconda** (recommended for Windows):
1. Download Anaconda from https://www.anaconda.com/download
2. Install Anaconda
3. Open Anaconda Prompt
4. Follow steps above

**Path Issues:**
If `python` command not found, try:
```bash
py run_pyrockmag.py
```

### macOS

**Xcode Command Line Tools** (needed for some packages):
```bash
xcode-select --install
```

### Linux

**Additional system dependencies** (Ubuntu/Debian):
```bash
sudo apt install python3-tk  # For matplotlib GUI
```

---

## Virtual Environment (Optional but Recommended)

Using a virtual environment keeps dependencies isolated:

**Create virtual environment:**
```bash
python -m venv pyrockmag_env
```

**Activate:**
- Windows: `pyrockmag_env\Scripts\activate`
- macOS/Linux: `source pyrockmag_env/bin/activate`

**Install packages:**
```bash
pip install numpy scipy matplotlib pandas
```

**Run application:**
```bash
python run_pyrockmag.py
```

**Deactivate when done:**
```bash
deactivate
```

---

## Troubleshooting

### "No module named numpy"
```bash
pip install numpy
```

### "matplotlib backend error"
**On Linux:**
```bash
sudo apt install python3-tk
```

**Set backend manually in code:**
Add to beginning of script:
```python
import matplotlib
matplotlib.use('TkAgg')  # or 'Qt5Agg'
```

### "Permission denied"
**On macOS/Linux:**
```bash
python3 run_pyrockmag.py
```

Or use pip with --user:
```bash
pip install --user numpy scipy matplotlib pandas
```

### "pip not found"
**Windows:**
```bash
python -m pip install numpy scipy matplotlib pandas
```

**macOS/Linux:**
```bash
python3 -m pip install numpy scipy matplotlib pandas
```

### Import errors with project modules
Make sure you're running from the pyRockmagCIT directory:
```bash
cd /path/to/pyRockmagCIT
python run_pyrockmag.py
```

---

## Updating

### Update from git:
```bash
cd pyRockmagCIT
git pull
```

### Update dependencies:
```bash
pip install --upgrade numpy scipy matplotlib pandas
```

---

## Uninstalling

### Remove repository:
```bash
rm -rf pyRockmagCIT  # macOS/Linux
rmdir /s pyRockmagCIT  # Windows
```

### Remove packages (if not needed for other projects):
```bash
pip uninstall numpy scipy matplotlib pandas
```

---

## Development Installation

For contributing to pyRockmagCIT:

```bash
git clone https://github.com/yicolas/pyRockmagCIT.git
cd pyRockmagCIT

# Create virtual environment
python -m venv venv
source venv/bin/activate  # or venv\Scripts\activate on Windows

# Install in editable mode with dev dependencies
pip install -e .
pip install pytest black flake8  # dev tools

# Run tests
pytest tests/

# Format code
black .
```

---

## System Requirements

### Minimum:
- CPU: 1 GHz processor
- RAM: 2 GB
- Disk: 500 MB free space
- Display: 1024×768

### Recommended:
- CPU: 2+ GHz multi-core
- RAM: 4+ GB
- Disk: 1+ GB free space
- Display: 1920×1080 or higher

---

## File Structure After Installation

```
pyRockmagCIT/
├── run_pyrockmag.py          # Main application
├── rmg_import.py              # Data import
├── rmg_stats.py               # Statistics
├── rmg_plots.py               # Plotting
├── rmg_forc.py                # FORC processing
├── rmg_sam.py                 # SAM generation
├── rmg_sam_utilities.py       # Sun compass/IGRF
├── rmg_sam_naming.py          # Filename validation
├── [other modules...]
├── README.md                  # Documentation
├── FEATURES.md
├── CHANGELOG.md
└── INSTALL.md                 # This file
```

---

## Getting Help

- Documentation: See README.md
- Issues: https://github.com/yicolas/pyRockmagCIT/issues
- Email: yick@duck.com

---

## Next Steps

After installation, see:
- **README.md** for feature overview
- **FEATURES.md** for detailed capabilities
- **SAM_Generation_Guide.md** for SAM file workflow
- Run `python run_pyrockmag.py` and explore the menu!

---

**pyRockmagCIT v1.0** - Ready to analyze!
