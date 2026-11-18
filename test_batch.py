"""Test der Serienauswertung: 12 Bilder mit Rauschen, dazu schlechte Bilder und ein Ausreisser."""
import numpy as np

from pdi_batch import analyze_series
from pdi_core import Params, to_gray
from pdi_synth import make_interferogram

BASE = {4: 0.30, 5: 0.15, 8: 0.20, 11: 0.35, 13: 0.05}
rng = np.random.default_rng(1)
frames, kinds = [], []
for s in range(12):
    # kleine zufaellige Aenderung von Kipp -> wie Vibration; Wahrheit bleibt gleich
    fr = (22 + rng.normal(0, 0.6), 6 + rng.normal(0, 0.6))
    img, truth = make_interferogram(coeffs=BASE, fringes=fr, noise=0.02, seed=s)
    frames.append(to_gray(img)); kinds.append("ok")
img, _ = make_interferogram(coeffs=BASE, contrast=0.04, noise=0.04, seed=50)         # fast ohne Streifen
frames.append(to_gray(img)); kinds.append("kontrastlos")
img, _ = make_interferogram(coeffs={**BASE, 8: 0.9, 5: -0.5}, seed=51)               # anderer Spiegelzustand (Luftschlieren)
frames.append(to_gray(img)); kinds.append("Ausreisser")
img, _ = make_interferogram(coeffs=BASE, fringes=(2.0, 0.5), seed=52)               # Traeger fehlt
frames.append(to_gray(img)); kinds.append("ohne Traeger")

ct = np.zeros(36)
for k, v in BASE.items():
    ct[k - 1] = v

p = Params()
res = analyze_series(frames, p)
print(f"verwendet {res.n_used} von {res.n_total}; Wiederholbarkeit {res.repeatability_rms:.4f} Wellen")
for i, inf in enumerate(res.frame_info):
    print(f"  {i:2d} {kinds[i]:12s} used={inf.get('used')} valid={inf.get('valid', float('nan')):.2f} resid={inf.get('resid', float('nan')):.3f}")
c = res.coeffs
s = np.sign(np.dot(c[4:], ct[4:]))
print("Fehler Serien-Median (j>=5): %.4f Wellen" % np.sqrt(np.sum((s * c[4:] - ct[4:]) ** 2)))
print("Kennzahlen:", {k: round(v, 3) if isinstance(v, float) else v for k, v in res.metrics.items()})
for w in res.warnings:
    print("!", w)

# mit Sollform: Vorzeichen automatisch
p2 = Params(mirror_d_mm=355.6, mirror_f_mm=1411)
res2 = analyze_series(frames[:12], p2)
print("mit Sollform: alpha %.3f +- %s, verwendet %d" % (res2.alpha, None if res2.alpha_std is None else round(res2.alpha_std, 3), res2.n_used))
