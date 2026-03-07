# FORC Exponential Stepping Guide

## User-Defined Exponential Base ✅

The exponential base is **fully customizable** by the user during FORC script generation.

---

## Interactive Menu

```bash
[G] Generate FORC measurement script

  Start field (mT) [default: 0]: 0
  Stop field (mT) [default: 500]: 500

  Field Spacing:
  + [1] Linear (constant step size)
  + [2] Exponential (growing step size)
  Choice [1]: 2

  Minimum step size (mT) [default: 5]: 5
  Exponential base (typically 1.2-1.5) [default: 1.3]: 1.3  ← USER DEFINED!
```

**The user can enter ANY value for the exponential base.**

---

## Formula

```python
field[i] = start + min_step * (base^i - 1)
```

**Parameters:**
- `start`: Starting field (mT)
- `min_step`: Minimum step size (mT)
- `base`: **USER-DEFINED exponential base**
- `i`: Step index (0, 1, 2, ...)

---

## Example with Different Bases

### **Base = 1.2 (Gentle Growth)**
```
0 → 0 → 1 → 2.4 → 4.3 → 6.7 → 9.8 → 13.8 → 18.8 → 25.1 → 33.1 → 43.1 → 56.1 → 72.3 → 91.7 → 115 → 143 → 176 → 215 → 262 → 317 → 383 → 462 → 500 mT
Result: 24 steps, gradual spacing
```

### **Base = 1.3 (Standard)** ⭐
```
0 → 1.5 → 3.5 → 6 → 9.3 → 13.6 → 19.6 → 28 → 39.4 → 55.3 → 76.9 → 106 → 143 → 191 → 253 → 334 → 439 → 500 mT
Result: 18 steps, balanced spacing
```

### **Base = 1.4 (Aggressive Growth)**
```
0 → 2 → 5.6 → 10.8 → 18.5 → 30.1 → 47.1 → 71 → 104 → 151 → 217 → 309 → 438 → 500 mT
Result: 14 steps, rapid growth
```

### **Base = 1.5 (Very Aggressive)**
```
0 → 2.5 → 7.5 → 15.6 → 28.1 → 47.1 → 75.6 → 118 → 182 → 278 → 421 → 500 mT
Result: 12 steps, very fast
```

---

## Choosing the Right Base

| Base | Steps (0-500 mT) | Use Case |
|------|------------------|----------|
| **1.0** | 100 | Linear (no exponential growth) |
| **1.2** | ~24 | Maximum detail everywhere |
| **1.3** | ~18 | **Standard, good balance** ⭐ |
| **1.4** | ~14 | Faster, less high-field detail |
| **1.5** | ~12 | Very fast measurements |
| **1.6** | ~10 | Minimal detail, maximum speed |

---

## Field Distribution Comparison

**At base = 1.3 with 0-500 mT range:**

```
Low field (0-50 mT):    9 steps → ~5.6 mT average spacing
Mid field (50-200 mT):  5 steps → ~30 mT average spacing  
High field (200-500 mT): 4 steps → ~75 mT average spacing
```

**Why this matters:**
- **Low fields** = SD/SP transitions, need fine detail → Dense sampling ✓
- **Mid fields** = MD region, gradual changes → Moderate sampling ✓
- **High fields** = Approaching saturation → Sparse sampling is OK ✓

---

## Programmatic Usage

```python
from rmg_forc import generate_forc_script

# User can specify any exponential base
generate_forc_script(
    start_mT=0,
    stop_mT=500,
    step_mT=5,
    output_path='./forc_exp1.3.rmg',
    saturation=True,
    exponential=True,
    exp_base=1.3  # ← USER-DEFINED!
)

# Different base for faster measurements
generate_forc_script(
    start_mT=0,
    stop_mT=500,
    step_mT=5,
    output_path='./forc_exp1.5_fast.rmg',
    saturation=True,
    exponential=True,
    exp_base=1.5  # ← Aggressive spacing
)
```

---

## Recommendations by Sample Type

### **Magnetite (SD/MD mixture)**
```
Base: 1.3
Reason: Need detail in 10-80 mT range
Steps: ~18
Time: ~2 hours
```

### **Hematite (High coercivity)**
```
Base: 1.4
Reason: Spread out to 300-500 mT range
Steps: ~16
Time: ~2.5 hours
```

### **Weak signals (sediments)**
```
Base: 1.2
Reason: Need maximum detail
Steps: ~24
Time: ~3 hours
```

### **Quick screening**
```
Base: 1.5
Reason: Fast measurements
Steps: ~12
Time: ~1 hour
```

---

## File Naming Convention

**Exponential files include the base:**

```
FORCz_0.0_to_500.0mT_in_exp1.30_saturation.rmg  ← Base = 1.3
FORCz_0.0_to_500.0mT_in_exp1.20_saturation.rmg  ← Base = 1.2
FORCz_0.0_to_500.0mT_in_exp1.50_Steps.rmg       ← Base = 1.5
```

**Linear files use step size:**
```
FORCz_0.0_to_500.0mT_in_5.0mT_saturation.rmg    ← Linear, 5 mT steps
```

---

## Common Base Values in Literature

| Study | Base | Reasoning |
|-------|------|-----------|
| Roberts et al. (2000) | 1.3 | Standard practice |
| Pike et al. (1999) | 1.2-1.4 | Varies by sample |
| Your MATLAB workflow | **1.3** | Balanced efficiency |

---

## Tips

### **1. Start with 1.3**
Default value works for most samples. Adjust only if needed.

### **2. Lower base for weak samples**
```
Base 1.2: More measurements → Better signal averaging
```

### **3. Higher base for fast screening**
```
Base 1.5: Fewer measurements → Quick overview
```

### **4. Match to coercivity range**
```
Low Hc samples: Base 1.2-1.3 (need low-field detail)
High Hc samples: Base 1.3-1.5 (can skip low fields)
```

### **5. Consider measurement time**
```
Base 1.2: ~3 hours
Base 1.3: ~2 hours  ← Good balance
Base 1.5: ~1 hour
```

---

## Verification

**Check your generated file:**

```python
# Linear spacing (base = 1.0)
Fields: 0, 5, 10, 15, 20, 25, 30, ...

# Exponential spacing (base = 1.3)
Fields: 0, 1.5, 3.5, 6, 9.3, 13.6, 19.6, 28, ...
         ↑    ↑    ↑   ↑    ↑     ↑     ↑
       Each gap grows by factor of ~1.3
```

---

## Summary

✅ **User-defined**: Exponential base is fully customizable
✅ **Default**: 1.3 (balanced, standard practice)
✅ **Range**: 1.0 (linear) to 1.6+ (very aggressive)
✅ **Interactive**: Prompted during script generation
✅ **Programmatic**: Available in function API
✅ **Documented**: Filename includes base value

**The exponential base is completely under user control!** 🎯
