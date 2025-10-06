# Exercise 10 – Moving Average Filters
tags: DSP, FIR, IIR, Moving Average, Frequency Response, Group Delay  
Author: Juan Rodriguez Esteban

## What this repository contains
- «Exercise10.m»: MATLAB script that implements 5‑point and 9‑point **moving‑average FIR** filters, computes magnitude and dB responses, shows pole–zero locations, and plots group delay. It also includes a safe «zero‑phase» fallback (forward–backward filtering) when «filtfilt» is not available.
- The code is **normalized** (DC gain = 1) for fair comparisons. It also shows unnormalized responses to visualise scaling only.

## How to run
1. Open «Exercise10.m» in MATLAB.
2. Run the script. No toolboxes are required.
3. Figures produced:
   - «Time‑Domain Responses»: input, causal outputs, and optional zero‑phase outputs.
   - «Frequency Responses»: linear magnitude, dB magnitude, and optional unnormalized comparison.
   - «Pole‑Zero Plots»: zeros on the unit circle and poles at the origin (FIR ⇒ «a = 1»).
   - «Group Delay»: measured delay vs the theoretical constant \((M-1)/2\) for linear‑phase FIRs.

## Results (at a glance)
- **Smoothing strength**: MA‑9 smooths more than MA‑5 (wider time window ⇒ stronger attenuation of high‑frequency components).
- **Transition bandwidth**: MA‑9 has a **narrower** main lobe in frequency, giving better separation of low vs high frequencies.
- **Delay**: causal moving averages are linear‑phase with constant group delay \((M-1)/2\). Thus, MA‑5 ≈ 2 samples; MA‑9 ≈ 4 samples. Using the zero‑phase path removes this delay for offline work.
- **Zeros**: an M‑point moving average places M‑1 zeros uniformly on the unit circle (excluding \(z=1\)), which creates deep notches at harmonic frequencies of \(\omega=2\pi/M\).

### Quick comparison: MA‑5 vs MA‑9
| Property | MA‑5 | MA‑9 | Practical effect |
|---|---:|---:|---|
| Window length \(M\) | 5 | 9 | Longer window averages more ⇒ smoother output |
| DC gain (normalized) | 1 | 1 | Same passband level at \(\omega=0\) |
| Group delay (samples) | \((5-1)/2 = 2\) | \((9-1)/2 = 4\) | Larger inherent latency for causal processing |
| Stop‑band attenuation | Lower | Higher | MA‑9 suppresses highs better (but with more ringing) |
| Transition width | Wider | Narrower | MA‑9 separates low/high frequencies better |

## Interpreting the figures
- «Time‑Domain»: MA‑9 looks smoother but lags more; «zero‑phase» curves align with input features while retaining smoothing.
- «Frequency Responses»: normalized curves start at 1 (DC) and roll off; MA‑9 shows a tighter low‑pass shape than MA‑5.
- «Pole‑Zero»: only zeros on the unit circle (and a pole at the origin due to FIR form) ⇒ unconditional stability and linear phase for symmetric coefficients.
- «Group Delay»: flat curves near \((M-1)/2\) confirm linear‑phase FIR behaviour; the dashed/aux lines indicate the theoretical values.

## FIR vs IIR – concise explanation
**FIR (Finite Impulse Response)** filters:
- Use **only current and past inputs** (non‑recursive), coefficient vector «b», and typically set «a = 1». Always **BIBO‑stable**. Linear‑phase designs (e.g., moving average, windowed‑sinc) are straightforward, and delay is predictable \((M-1)/2\).

![Figure 1 Time-Domain Responses.bmp](https://github.com/user-attachments/files/22728392/Figure.1.Time-Domain.Responses.bmp)


- Output equals **convolution** of the input with the impulse response; «filter(b,1,x)» and «conv(x,b,'same')» match in practice.

![Figure 2 Frequency Responses.bmp](https://github.com/user-attachments/files/22728397/Figure.2.Frequency.Responses.bmp)


**IIR (Infinite Impulse Response)** filters:
- Use **feedback** (depend on past outputs) with numerator «b» and denominator «a». They can achieve a given response with **lower order** (fewer taps) but require **stability** checks (poles must be inside the unit circle). Butterworth/Chebyshev/Biquad designs are typical.

  ![zero plot.bmp](https://github.com/user-attachments/files/22728404/zero.plot.bmp)

  
- IIRs usually have **non‑linear phase**, so they distort waveform timing; however, offline «filtfilt» can eliminate phase distortion by forward–backward filtering.

![Figure 4 Group Delay.bmp](https://github.com/user-attachments/files/22728408/Figure.4.Group.Delay.bmp)


**When to choose which**
- Choose **FIR** for linear‑phase requirements, robust stability, or arbitrary magnitude shapes (via design methods like windowed «fir1»). fileciteturn2file8
- Choose **IIR** when efficiency (sharp transition with few coefficients) matters and moderate phase distortion is acceptable—or when «filtfilt» is available for offline use. A side‑by‑side demo shows their different delay characteristics and smoothness.


## Reproducibility notes
- The script intentionally avoids «freqz», «tf2zpk», «grpdelay», and «zplane» dependencies by bundling safe fallbacks. If you do have the Signal Processing Toolbox or «signal» package, the wrappers will call them automatically.
- For offline «zero‑phase» runs without «filtfilt», the script uses a simple forward–backward implementation with edge reflection.

---

© 2025 Juan Rodriguez Esteban
