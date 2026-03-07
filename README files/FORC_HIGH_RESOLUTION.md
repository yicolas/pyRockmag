# FORC High-Resolution Plots

## New Feature: Adjustable Plot Resolution ✨

Control the quality and detail of your FORC diagrams with selectable DPI (dots per inch) settings.

## Quick Reference

### **Resolution Options:**

```
▼══════════════════════════════════════════════════════════▼
           Plot Resolution
▲══════════════════════════════════════════════════════════▲

Resolution options:
  + [1] Screen (100 DPI) - Fast, lower quality
  + [2] Standard (150 DPI) - Good balance
  + [3] High (300 DPI) - Publication quality
  + [4] Ultra (600 DPI) - Maximum detail

Choice [3]: 
```

### **What Changed:**

**Figure Size:**
- Old: 9" × 7" 
- New: 12" × 10" (33% larger)

**Resolution Options:**
- Screen: 100 DPI → 1200 × 1000 pixels
- Standard: 150 DPI → 1800 × 1500 pixels  
- High: 300 DPI → 3600 × 3000 pixels ⭐ (Default)
- Ultra: 600 DPI → 7200 × 6000 pixels

---

## Use Cases

### **1. Screen Viewing (100 DPI)**
- Quick previews
- Fast iteration
- Lower memory usage
- **File size:** ~200 KB

### **2. Standard Quality (150 DPI)**
- General use
- Presentations (PowerPoint/Keynote)
- Email/web sharing
- **File size:** ~400 KB

### **3. Publication Quality (300 DPI)** ⭐ Recommended
- Journal submissions
- High-quality figures
- Print publications
- **File size:** ~1.5 MB
- **Most journals require 300 DPI minimum**

### **4. Ultra High (600 DPI)**
- Large format posters
- Maximum zoom capability
- Archival quality
- **File size:** ~6 MB
- **Use when fine details are critical**

---

## Resolution Comparison

### Pixel Counts at 12" × 10" Size:

| DPI | Total Pixels | Best For |
|-----|-------------|----------|
| 100 | 1.2 megapixels | Screen preview |
| 150 | 2.7 megapixels | Presentations |
| 300 | 10.8 megapixels | **Publications** ⭐ |
| 600 | 43.2 megapixels | Posters/archival |

### Cropped Plot Detail:

At **Bc: 0-100 mT** zoom:

**100 DPI:** ~360 pixels across Bc axis
- Adequate for viewing
- Limited zoom capability

**300 DPI:** ~1080 pixels across Bc axis ⭐
- Excellent detail
- Clean smooth contours
- Publication ready

**600 DPI:** ~2160 pixels across Bc axis
- Maximum detail
- Perfect for posters
- Overkill for most uses

---

## Workflow Examples

### Example 1: Quick Analysis
```bash
[5] FORC contour diagram

Resolution: [1] Screen (100 DPI)
Colormap: [2] inferno
Bc max: 100 mT

→ Fast generation
→ Good enough to see features
→ Small file for quick sharing
```

### Example 2: Paper Submission
```bash
[5] FORC contour diagram

Resolution: [3] High (300 DPI)  ← Default, press Enter
Colormap: [5] gray  ← B&W for print
Bc: 0-80 mT, Bu: ±50 mT

→ Publication quality
→ Meets journal requirements
→ Perfect for manuscripts
```

### Example 3: Conference Poster
```bash
[5] FORC contour diagram

Resolution: [4] Ultra (600 DPI)
Colormap: [2] inferno
Bc: 0-100 mT

→ Large format ready
→ Looks sharp at poster size
→ Can be viewed up close
```

---

## Technical Details

### How It Works

1. **DPI Selection** sets `matplotlib.rcParams['figure.dpi']`
2. **Figure Size** increased to 12" × 10" 
3. **Pixel Count** = Width(in) × Height(in) × DPI²
4. **File Size** scales roughly with pixel count

### Memory Usage

| DPI | RAM Usage | Generation Time |
|-----|-----------|-----------------|
| 100 | ~10 MB | 1-2 seconds |
| 150 | ~20 MB | 2-3 seconds |
| 300 | ~80 MB | 4-6 seconds |
| 600 | ~300 MB | 10-15 seconds |

### Recommended Settings by Purpose

**Exploratory Analysis:**
- DPI: 100-150
- Focus: Speed

**Presentations:**
- DPI: 150-300
- Focus: Clarity on screen

**Journal Papers:**
- DPI: 300 (minimum)
- Focus: Print quality
- Check journal requirements!

**Posters/Archival:**
- DPI: 600
- Focus: Maximum detail

---

## Tips & Best Practices

### 1. **Start with Default (300 DPI)**
Press Enter to use 300 DPI - suitable for 95% of uses.

### 2. **Use 600 DPI Sparingly**
Only when you need:
- Large format printing
- Extreme zoom capability
- Final archival version

### 3. **Match to Output**
- Screen only? → 100-150 DPI
- PowerPoint? → 150 DPI
- Paper? → 300 DPI
- Poster? → 600 DPI

### 4. **Cropped Plots Need Higher DPI**
When zooming to Bc 0-50 mT:
- Same pixels cover smaller range
- Higher DPI maintains smoothness
- 300 DPI recommended minimum

### 5. **File Management**
600 DPI files can be large (>5 MB):
- Use for final versions only
- Keep 300 DPI for working copies
- Compress for email if needed

---

## Before & After

### Low Resolution (100 DPI):
```
Pixels across Bc axis: ~360
Contour lines: Slightly jagged
Zoom capability: Limited
File size: 200 KB
```

### High Resolution (300 DPI):
```
Pixels across Bc axis: ~1080
Contour lines: Smooth and clean
Zoom capability: Excellent
File size: 1.5 MB
Publication ready: ✓
```

### Ultra Resolution (600 DPI):
```
Pixels across Bc axis: ~2160
Contour lines: Perfect
Zoom capability: Maximum
File size: 6 MB
Poster ready: ✓
```

---

## Common Questions

**Q: Why is 300 DPI the default?**
A: It meets publication requirements for most journals while keeping reasonable file sizes.

**Q: Can I change DPI after generating?**
A: No, you must regenerate. DPI affects how the plot is rendered, not just saved.

**Q: Will 600 DPI look better on screen?**
A: Not significantly. Screens are typically 72-150 DPI, so 300 DPI is already oversampled for viewing.

**Q: Which DPI for a presentation?**
A: 150 DPI is perfect. Projectors rarely exceed 100 DPI, so 150 DPI provides margin.

**Q: File too large?**
A: Use PNG with compression or reduce DPI. Or use vector format (SVG) for infinite zoom.

---

## Comparison Chart

| Purpose | DPI | Figure Size | Pixels | File Size |
|---------|-----|-------------|--------|-----------|
| Quick preview | 100 | 12×10" | 1200×1000 | ~200 KB |
| Presentation | 150 | 12×10" | 1800×1500 | ~400 KB |
| **Publication** | **300** | **12×10"** | **3600×3000** | **~1.5 MB** |
| Poster | 600 | 12×10" | 7200×6000 | ~6 MB |

---

**pyRockmagCIT v1.0** - Now with publication-quality high-resolution FORC diagrams!
