# FORC Diagram - Interactive Axis Limit Control

## New Feature: Bc and Bu Axis Adjustment ✨

**Focus on the critical zone of your FORC diagrams!**

### What It Does

Allows you to interactively zoom into specific regions of FORC diagrams by setting custom Bc (coercivity) and Bu (interaction field) axis limits after viewing the initial auto-scaled plot.

### Why This Matters

**Problem:** Full-range FORC plots often show mostly empty space, making it hard to see details in the magnetically active region.

**Solution:** View the full plot first, identify the interesting region, then zoom in!

**Example:**
- Full range: Bc 0-300 mT → mostly white space
- Zoomed: Bc 0-100 mT → clear magnetic interactions visible

---

## How to Use

### Workflow

```bash
python run_pyrockmag.py
[P] Plot FORC data
[7] FORCinel workflow (complete processing)

# 1. Initial plot displays (auto-scaled)
#    You see the full data range

# 2. Prompt appears:
▼══════════════════════════════════════════════════════════▼
           FORC Axis Limit Adjustment
▲══════════════════════════════════════════════════════════▲

Current Bc range (auto): 0.0 - 287.3 mT
Current Bu range (auto): -45.2 - 45.2 mT

+ Adjust axis limits? (y/n) [n]: y

# 3. Set new Bc limits
--- Bc Axis (Coercivity) ---
Bc min (mT) [0.0]: 0
Bc max (mT) [287.3]: 100

# 4. Set new Bu limits
--- Bu Axis (Interaction Field) ---
Bu min (mT) [-45.2]: -30
Bu max (mT) [45.2]: 30

✓ New limits: Bc [0.0, 100.0] mT, Bu [-30.0, 30.0] mT

# 5. Plot regenerates with new limits
#    Much clearer view of magnetic interactions!

# 6. Adjust again if needed
+ Adjust limits again? (y/n) [n]: n

# 7. Final plot saved
✓ Saved plot: sample_forc_SF3.png
```

---

## Use Cases

### 1. **Fine-Grained Sediments** (Low Coercivity)
```
Initial range: Bc 0-300 mT
Most data < 50 mT
Solution: Set Bc max = 50 mT
Result: 6x magnification of interesting region
```

### 2. **Magnetite Samples** (Moderate Coercivity)
```
Initial range: Bc 0-200 mT
Peak at 20-80 mT
Solution: Bc min = 10, Bc max = 100 mT
Result: Focus on magnetite interaction zone
```

### 3. **High Coercivity Samples** (Hematite, Goethite)
```
Initial range: Bc 0-500 mT
Data 100-400 mT
Solution: Bc min = 80, Bc max = 450 mT
Result: Exclude empty low-field region
```

### 4. **Strong Interactions**
```
Initial Bu range: -100 to +100 mT
Most activity ±30 mT
Solution: Bu min = -40, Bu max = 40 mT
Result: Better contrast in interaction field
```

---

## Examples

### Example 1: Magnetite FORC Diagram

**Initial Plot (auto):**
```
Bc range: 0 - 285 mT
Bu range: -48 - 48 mT
```
Most magnetic signal concentrated in Bc = 10-70 mT region.

**After Adjustment:**
```
+ Adjust axis limits? y
Bc min: 0
Bc max: 80
Bu min: -35
Bu max: 35
```

**Result:** Clear visualization of magnetite SD/PSD transition at Bc ~40 mT.

---

### Example 2: Volcanic Glass (Weak Signal)

**Initial Plot:**
```
Bc range: 0 - 150 mT
Bu range: -20 - 20 mT
```
Weak FORC signal barely visible.

**After Adjustment:**
```
Bc max: 60     # Exclude empty high-field region
Bu min: -10    # Tighten Bu range
Bu max: 10
```

**Result:** Enhanced visibility of weak magnetic interactions.

---

### Example 3: Iteration to Perfect View

**First Try:**
```
Bc: 0-100 mT
```
Still too wide.

**Second Adjustment:**
```
+ Adjust limits again? y
Bc: 20-80 mT    # Refined range
Bu: -25-25 mT   # Tightened
```

**Result:** Publication-quality focused view!

---

## Technical Details

### Implementation

- **When:** After initial plot generation, before saving
- **How:** Interactive prompts for Bc and Bu min/max
- **Where:** FORCinel workflow ([P] → [7])
- **Repeat:** Can adjust multiple times until satisfied

### Axis Definitions

**Bc (Coercivity):**
- Bc = (Ha - Hr) / 2
- Physical meaning: Average coercive field
- Typical range: 0-500 mT (sample dependent)
- Focus zone often: 10-100 mT for magnetite

**Bu (Interaction Field):**
- Bu = (Ha + Hr) / 2  
- Physical meaning: Magnetic interaction strength
- Typical range: ±50 mT
- Zero = non-interacting particles

### Default Behavior

- **No input:** Press Enter to keep current limits
- **Auto-scale:** Don't adjust (n) to use full range
- **Invalid input:** Falls back to auto-scale

---

## Tips & Best Practices

### 1. **Always View Full Plot First**
- See the complete data distribution
- Identify where the signal is concentrated
- Then zoom in strategically

### 2. **Leave Some Margin**
- Don't crop exactly at data edges
- Add ~10% margin for context
- Example: Data ends at 95 mT → set max to 105 mT

### 3. **Symmetric Bu Range** (Usually)
- Magnetic interactions often symmetric
- If data is -35 to +38 mT → use -40 to +40 mT
- Exception: Samples with strong bias field

### 4. **Iterate If Needed**
- First zoom: Broad (cut obvious empty space)
- Second zoom: Refined (optimal view)
- Third zoom: Publication perfection

### 5. **Different Samples, Different Ranges**
- Magnetite: Often Bc 0-80 mT
- Hematite: Often Bc 100-400 mT
- Goethite: Often Bc 200-500 mT
- Adjust per sample characteristics

### 6. **Save Original Too**
- First save: Auto-scaled (full view)
- Second save: Zoomed (focused view)
- Both useful for different purposes

---

## Comparison

### Before (Auto-Scale Only)

**Limitations:**
- ✗ Empty space dominates plot
- ✗ Hard to see fine details
- ✗ Poor use of color scale
- ✗ Not publication-ready

**Example:**
```
Magnetite sample: Bc 0-300 mT
Actual signal: Bc 15-65 mT (only 17% of plot!)
83% of plot is white space
```

### After (Interactive Zoom)

**Benefits:**
- ✓ Focused on magnetic signal
- ✓ Clear detail resolution
- ✓ Optimal color scale usage  
- ✓ Publication-quality

**Example:**
```
Zoomed to: Bc 0-80 mT
Signal occupies: ~60% of plot
Details clearly visible
Perfect for papers!
```

---

## Workflow Integration

### Research Paper Figures

```bash
# 1. Generate auto-scaled for overview
#    Shows full magnetic range

# 2. Zoom to critical zone
#    Bc: 0-100 mT (magnetite focus)
#    Bu: ±30 mT (interaction detail)

# 3. Save both versions
#    - Full: Supplementary material
#    - Zoomed: Main figure

# 4. Perfect for publication!
```

### Comparative Studies

```bash
# Process multiple samples
# Each gets custom zoom based on its characteristics

Sample A (magnetite): Bc 0-80 mT
Sample B (hematite): Bc 150-350 mT  
Sample C (mixture): Bc 0-200 mT

# All optimized individually!
```

---

## Future Enhancements

Planned for v1.1:
- [ ] Save zoom presets
- [ ] Batch apply same zoom to multiple samples
- [ ] Automatic optimal zoom detection
- [ ] Export zoom parameters to log file
- [ ] Compare before/after side-by-side

---

## Troubleshooting

**Q: Zoomed plot looks empty?**
A: Limits too narrow. Increase Bc/Bu range.

**Q: Still see white space?**
A: Try narrower limits, focus more tightly.

**Q: Lost some data?**
A: Check if your limits excluded real signal. View full plot again.

**Q: Can't see interactions?**
A: Tighten Bu range to enhance weak interactions.

**Q: Want to reset?**
A: Just press 'n' and re-run to get auto-scale again.

---

## Command Reference

**During FORC workflow:**
```
+ Adjust axis limits? (y/n) [n]
  → y: Set custom limits
  → n: Keep current (auto or previous)

+ Adjust limits again? (y/n) [n]
  → y: Refine further
  → n: Finalize and save
```

**Input format:**
```
Bc min (mT) [default]: <number or Enter>
Bc max (mT) [default]: <number or Enter>
Bu min (mT) [default]: <number or Enter>
Bu max (mT) [default]: <number or Enter>
```

---

**pyRockmagCIT v1.0** - Now with interactive FORC axis control for perfect diagrams!
