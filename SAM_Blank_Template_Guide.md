# Blank SAM Template Generation - Detailed Guide

## Overview

The blank template generator creates empty `.sam` header files and sample files with user-defined default values. Perfect for field work where you'll fill in orientation data later.

---

## Quick Start

```bash
python run_rmg_analysis.py
[S] Generate SAM header files
[1] Create blank SAM template
```

Follow the prompts - the system will ask for all necessary information.

---

## Interactive Prompts

### 1. Site ID
```
Site ID (max 4 chars):
```
- **Maximum 4 characters** (CIT format requirement)
- Examples: `GB20`, `NWF`, `TEST`
- This becomes the filename prefix

### 2. Site Description
```
Site description [blank]:
```
- Optional descriptive text
- Appears in .sam header file (line 2)
- Press Enter to use default: "{Site ID} Site"
- Example: "Gooseberry Basalts Flow 20"

### 3. Site Coordinates
```
Set site coordinates? (y/n) [n]:
```

**If NO (default):**
- Lat/Lon/Elevation set to 0
- No IGRF calculations
- Perfect for blank templates

**If YES:**
```
Latitude (decimal degrees N):
Longitude (decimal degrees E):
Elevation (meters) [0]:
```
- Enter decimal degrees (e.g., 47.8234)
- Elevation in meters above sea level
- Used for reference only (no calculations for blank templates)

### 4. Number of Samples
```
Number of samples to create:
```
- Enter the number of cores/samples
- Creates one file per sample
- Typical: 10-25 samples per site

### 5. Sample Naming
```
Sample naming (examples: 1.0a, 2.0a, ...)
Sample name prefix [blank = just numbers]:
```

**If blank (press Enter):**
- Samples named: `1.0a`, `2.0a`, `3.0a`, ...

**If you enter a prefix (e.g., "A"):**
- Samples named: `A1.0a`, `A2.0a`, `A3.0a`, ...
- Useful for distinguishing different groups

### 6. Core Plate Strike (DEFAULT VALUES FOR ALL SAMPLES)
```
Set core plate strike to 90° for all samples? (y/n) [y]:
```

**Why 90°?**
- Core plate strike = perpendicular to core axis
- **90° is the standard default** for most drilling
- Means core axis is horizontal (0° plunge)

**If YES (default - recommended):**
- All samples get 90° core plate strike

**If NO:**
```
Core plate strike (°) [0]:
```
- Enter custom value
- Applied to ALL samples

### 7. Core Plate Dip
```
Core plate dip (°) [0]:
```
- Default: 0° (horizontal cores)
- Dip of plane perpendicular to core axis
- Applied to all samples

### 8. Bedding Orientation
```
Bedding strike (°) [0]:
Bedding dip (°) [0]:
```
- Right-hand rule convention
- Default: 0°/0° (horizontal beds)
- Applied to all samples
- Edit individual files later if beds vary

### 9. Sample Mass
```
Sample mass (g) [1.0]:
```
- Enter typical sample mass in grams
- Default: 1.0 g
- Applied to all samples
- Edit individual files for actual masses later

### 10. Sample Description
```
Sample description/comment [blank template]:
```
- Descriptive text for all samples
- Examples:
  - "blank template"
  - "basalt cores"
  - "lake sediment"
  - "syn 2024-07-15"

### 11. Save Location
```
Save location [current directory]:
```
- Path where files will be created
- Press Enter for current directory
- Directory created automatically if needed
- Examples:
  - `.` (current directory)
  - `./field_data`
  - `/path/to/paleomag/GB20`

---

## Output Summary

After entering all information, you'll see:

```
======================================================================
BLANK TEMPLATE CREATED
======================================================================
✓ Site: GB20
✓ Samples: 10
✓ Location: ./field_data
✓ .sam file: ./field_data/GB20.sam
✓ Coordinates: Not set (0, 0)

Default values:
  Core strike: 90.0°
  Core dip: 0.0°
  Bedding: 0.0°/0.0°
  Mass: 1.0 g
  Comment: blank template

Edit individual sample files to modify values.
```

---

## Files Created

### .sam Header File
```
CIT
GB20 Site
   0.0 000.0   0.0
GB201.0a
GB202.0a
...
```

### Sample Files (one per core)
```
GB20 1.0a blank template
    0.0     90.0   0.0   0.0   0.0   1.0
```

Format: `strat_level  core_strike  core_dip  bed_strike  bed_dip  mass`

---

## Example Workflow

### Field Scenario

**Day 1 - Before Field Work:**
```bash
# Create blank template on laptop
python run_rmg_analysis.py → [S] → [1]

Site ID: NWF
Site description: North West Flow 2024
Set coordinates: n
Number of samples: 20
Sample prefix: [blank]
Core strike 90°: y
Core dip: 0
Bedding strike: 180
Bedding dip: 25
Sample mass: 10.5
Description: basalt core
Save location: ./NWF_2024

# Result: NWF.sam + 20 sample files
```

**Day 2 - In the Field:**
- Print the sample files or bring laptop
- Fill in actual orientation data as you core
- Edit files with text editor to update:
  - Core strike (from compass)
  - Individual sample masses
  - Notes/comments

**Day 3 - Back in Lab:**
- Transfer updated files to magnetometer computer
- Place in RAPID software folder
- Ready for measurement!

### Lab Scenario

**Preparing for Measurement:**
```bash
# Quick blank template
python run_rmg_analysis.py → [S] → [1]

Site ID: TEST
[Use all defaults, press Enter repeatedly]
Save location: ./RAPID/TEST

# Edit files manually with measured values
# Then run RAPID software
```

---

## Editing Sample Files

Sample files are plain text - edit with any editor:

```bash
# Before (blank template)
GB20 1.0a blank template
    0.0     90.0   0.0   0.0   0.0   1.0

# After (filled in)
GB20 1.0a measured 2024-07-15 14:30
    5.2    185.3  35.2  178.0  23.5  12.3
```

Changes:
- Comment: Added measurement date
- Strat level: 5.2 m
- Core strike: 185.3° (from compass)
- Core dip: 35.2° (from inclinometer)
- Bedding: 178.0°/23.5° (from field measurements)
- Mass: 12.3 g (actual measured mass)

---

## Tips and Best Practices

### 1. Coordinate Entry
- **For blank templates**: Skip coordinates (just press 'n')
- **For reference**: Enter if you want location recorded
- **Remember**: No IGRF calculations for blank templates

### 2. Core Strike Default
- **90° is standard** - use unless you have specific reason
- Assumes horizontal drilling (0° plunge)
- Can be edited per sample later if needed

### 3. Sample Naming
- **Keep it simple**: Usually just numbers work fine
- **Use prefix** if you have multiple groups (A-cores, B-cores)
- **Format**: Always ends in `.0a` (CIT convention)

### 4. Default Values
- Set to **most common** values for your site
- Individual files can be edited later
- Saves time if most samples are similar

### 5. File Organization
- **Create separate directories** for each site
- **Match folder name to Site ID** for RAPID software
- Example structure:
  ```
  paleomag_2024/
    GB20/
      GB20.sam
      GB201.0a
      GB202.0a
      ...
    NWF/
      NWF.sam
      NWF1.0a
      ...
  ```

### 6. Field Work Efficiency
- **Print sample files** or bring laptop
- Fill in values as you core
- Or write in field notebook and type later
- **Check file format** before going to magnetometer

### 7. Quality Control
- **Verify file names** match .sam header
- **Check units**: Degrees for angles, grams for mass
- **Test one file** in RAPID before processing all samples

---

## Troubleshooting

### "Site ID too long"
- Limit: 4 characters maximum
- Solution: Use abbreviations (e.g., "NWFB" not "NorthWestFlowBasalt")

### "Invalid number"
- Check that you entered numbers, not text
- Use decimals correctly: `10.5` not `10,5`
- Press Enter to use defaults

### "Directory not found"
- Directory is created automatically
- Check path permissions
- Use `.` for current directory if unsure

### "RAPID can't read files"
- Verify folder name = Site ID
- Check file has Windows line endings (\r\n)
- Open file in text editor - should show properly
- Ensure no special characters in filenames

---

## Advanced: Programmatic Use

```python
from rmg_sam import generate_blank_sam_template

# Non-interactive mode
sam_path = generate_blank_sam_template(
    site_id='GB20',
    num_samples=25,
    output_dir='./field_data',
    interactive=False  # Skip prompts, use defaults
)

# Result: 90° core strike, 0° dip, 1.0g mass, "blank template" comment
```

---

## Comparison: Blank vs Interactive

| Feature | Blank Template | Interactive Creation |
|---------|----------------|---------------------|
| Prompts | Detailed (10 questions) | Detailed (15+ questions) |
| Coordinates | Optional (usually skip) | Required |
| IGRF calculation | **Never** | **Always** |
| Sun compass | No | Yes (optional) |
| Validation | No | Yes (sun vs mag) |
| Use case | Field prep | Lab analysis |
| Time | 2 minutes | 5-10 minutes |

**Choose Blank Template when:**
- Preparing files before field work
- Will measure orientations later
- Don't need calculations
- Want simple, fast setup

**Choose Interactive when:**
- Have all orientation data
- Want automatic calculations
- Need sun compass processing
- Want IGRF corrections

---

## Next Steps

After creating blank template:

1. **Field Work:**
   - Bring files/printouts to field
   - Record measurements
   - Update files

2. **Lab Prep:**
   - Edit files with actual values
   - Verify all data entered
   - Check file format

3. **Measurement:**
   - Copy to RAPID folder
   - Open in RAPID software
   - Run magnetometer

4. **Analysis:**
   - Import measurement data back to rmg_python
   - Generate plots and statistics
   - Publish results!

