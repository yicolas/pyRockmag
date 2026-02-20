# FORC Scatter Plot - Brighter Colors

## Changes Made ✨

Your scatter plots now have **brighter, more vibrant colors** by default!

### **What Changed:**

| Parameter | Old Value | New Value | Effect |
|-----------|-----------|-----------|--------|
| **Point size** | s=1 | s=3 | 3× larger points = better coverage |
| **Opacity** | alpha=0.6 | alpha=0.9 | 90% opaque = more saturated colors |
| **Color range** | 1-99 percentile | 5-95 percentile | Better contrast |

---

## Before & After

### **Before (Muted):**
```
Point size: 1 pixel (tiny, sparse)
Transparency: 40% (washed out)
Color range: Very wide (low contrast)
Result: Pale, faded appearance
```

### **After (Vibrant):**
```
Point size: 3 pixels (visible, good coverage)
Transparency: 10% (saturated, opaque)
Color range: Tighter (high contrast)
Result: Bright, punchy colors!
```

---

## Why Scatter Was Muted

### **1. Tiny Points (s=1)**
- 1-pixel points are barely visible
- Lots of gaps between points
- Colors don't "fill in" the space
- Overall pale appearance

### **2. High Transparency (alpha=0.6)**
- 40% transparent = colors blend with white background
- Reduces color saturation
- Makes everything look washed out
- Hard to see weak signals

### **3. Wide Color Range (1-99 percentile)**
- Includes extreme outliers
- Compresses most data into middle of colormap
- Reduces contrast
- Colors look similar

---

## New Settings Explained

### **Point Size: s=3**

**Impact:**
- Each point is 3×3 pixels instead of 1×1
- **9× more coverage** per point
- Points overlap → continuous color
- No more gaps/white spaces

**When to adjust:**
- Very dense data? → Try s=2
- Sparse data? → Try s=4 or s=5
- High-res plots (600 DPI)? → Try s=5

### **Alpha: 0.9 (90% opaque)**

**Impact:**
- Colors are 90% saturated
- Only 10% transparency
- Much more vibrant appearance
- Still slight blending at overlap

**Why not 1.0 (fully opaque)?**
- 0.9 allows slight smoothing where points overlap
- Prevents harsh boundaries
- Maintains some depth perception

### **Color Range: 5-95 Percentile**

**Impact:**
- Ignores bottom 5% and top 5%
- Focuses colormap on main data distribution
- **Better contrast** in regions that matter
- Outliers don't dominate the scale

**Example:**
- Old: vmin=0.001, vmax=3.0 → Most data at 0.1-0.5 uses only 15% of colormap
- New: vmin=0.05, vmax=0.8 → Main data uses 80% of colormap = vibrant!

---

## Colormap Recommendations for Brightness

Now that scatter is brighter, different colormaps really shine:

### **1. Inferno (Default) ⭐**
```
Best for: Most FORC data
Why: High contrast, perceptually uniform
Colors: Purple → Red → Yellow
Bright mode: Stunning!
```

### **2. Viridis**
```
Best for: Weak signals, subtle features
Why: Colorblind-friendly, professional
Colors: Purple → Green → Yellow
Bright mode: Clear and scientific
```

### **3. Rainbow**
```
Best for: Presentations, maximum color
Why: Full spectrum, eye-catching
Colors: Blue → Green → Yellow → Red
Bright mode: Super vibrant!
```

### **4. Hot**
```
Best for: Classic FORC look
Why: Familiar, high dynamic range
Colors: Black → Red → Yellow → White
Bright mode: Dramatic contrast
```

### **5. Spring**
```
Best for: Fun, colorful displays
Why: Bright magentas and yellows
Colors: Magenta → Yellow
Bright mode: Electric!
```

### **6. Gray**
```
Best for: Publications (B&W)
Why: Print-friendly, clear structure
Colors: Black → White
Bright mode: Strong contrast
```

---

## Impact on Different Data Types

### **Weak FORC Signals:**
- **Old:** Barely visible, very pale
- **New:** Clear, distinguishable ✓

### **Strong Signals:**
- **Old:** Washed out even when strong
- **New:** Punchy, vibrant ✓

### **Sparse Data:**
- **Old:** Lots of white gaps
- **New:** Better coverage, smoother ✓

### **Dense Data:**
- **Old:** Muddy due to overlapping transparent points
- **New:** Rich, saturated colors ✓

---

## Technical Details

### **Point Rendering:**

At **300 DPI with s=3**:
- Each point covers ~9 pixels
- At 200×200 grid cropped to 0-100 mT:
  - ~40 points across Bc
  - Each point = ~27 screen pixels wide
  - Good overlap = continuous color

### **Alpha Blending Math:**

**Old (alpha=0.6):**
- 1 point: 60% color
- 2 overlapping: 84% color (60% + 40%×60%)
- 3 overlapping: 94% color

**New (alpha=0.9):**
- 1 point: 90% color (much brighter!)
- 2 overlapping: 99% color (nearly opaque)
- 3+ overlapping: 100% color (fully saturated)

### **Percentile Impact:**

Example FORC data distribution:
```
1st percentile:  0.001 (noise floor)
5th percentile:  0.020 (real signal starts)
50th percentile: 0.150 (median)
95th percentile: 0.650 (strong signals)
99th percentile: 2.800 (outliers/artifacts)
```

**Old scale (1-99%):**
- Range: 0.001 to 2.800 (factor of 2800×)
- Main signal (0.02-0.65) uses only 23% of colormap

**New scale (5-95%):**
- Range: 0.020 to 0.650 (factor of 32×)
- Main signal uses 100% of colormap = vibrant!

---

## Examples

### **Example 1: Magnetite Sample**

**Old settings:**
```
Bc 0-100 mT region: Pale yellow-orange
Hard to distinguish features
Low contrast between SD and MD
```

**New settings:**
```
Bc 0-100 mT region: Vibrant yellow-red-purple
Clear SD peak at ~40 mT (bright yellow!)
MD region shows as rich orange
Excellent contrast ✓
```

### **Example 2: Weak Volcanic Glass**

**Old settings:**
```
Signal barely visible
Looked like noise
Very pale colors throughout
```

**New settings:**
```
Signal clearly visible!
Distinct magnetic components
Good color separation
Easy to analyze ✓
```

---

## If You Want Even Brighter

Want to go **maximum vibrant**? Edit `rmg_forc.py`:

### **Option 1: Larger Points**
```python
# Line ~80 in plot_forc_diagram_standard
scatter = ax.scatter(..., s=5, alpha=0.9, ...)  # Change s=3 to s=5
```

### **Option 2: Fully Opaque**
```python
scatter = ax.scatter(..., s=3, alpha=1.0, ...)  # Change 0.9 to 1.0
```

### **Option 3: Tighter Color Range**
```python
# Line ~60
vmin_auto = np.percentile(rho_flat, 10)  # Change 5 to 10
vmax_auto = np.percentile(rho_flat, 90)  # Change 95 to 90
```

---

## Comparison with Contour Plots

| Feature | Scatter (Bright) | Contour |
|---------|-----------------|---------|
| **Color saturation** | Very high | Medium |
| **Detail level** | Individual points | Smooth regions |
| **Best for** | Exploring data | Publication figures |
| **Brightness** | Maximum | Controlled |
| **File size** | Smaller | Larger |

**Both are now optimized!**
- Scatter: Bright, vibrant, exploratory
- Contour: Smooth, professional, polished

---

## Tips for Best Results

### **1. Match Colormap to Brightness**
- Bright scatter + inferno = Perfect! ⭐
- Bright scatter + rainbow = Very colorful!
- Bright scatter + gray = High contrast B&W

### **2. Use with Higher Grid Resolution**
- 200×200 grid + bright scatter = Great
- 300×300 grid + bright scatter = Excellent coverage

### **3. Cropped Plots Benefit Most**
- Full range: Noticeable improvement
- Cropped 0-100 mT: **Dramatic improvement!** ✓

### **4. Combine with High DPI**
- 300 DPI + bright scatter = Publication quality
- 600 DPI + bright scatter = Perfect for posters

---

**Your FORC scatter plots are now vibrant and publication-ready!** 🌈✨
