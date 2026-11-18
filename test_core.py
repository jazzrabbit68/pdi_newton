"""Testet pdi_core gegen kuenstliche Interferogramme mit bekannter Wellenfront."""
import sys

import numpy as np

from pdi_core import Params, analyze, to_gray, zernike_name
from pdi_synth import make_interferogram

BASE = {4: 0.30, 5: 0.15, 8: 0.20, 11: 0.35, 13: 0.05}

SCEN = {
    "A normal": dict(coeffs=BASE),
    "B schwache Zone + Abfall": dict(coeffs=BASE, weak_zones=[(-0.45, 0.25, 0.22, 0.08)], falloff=0.7),
    "C uebersteuert (gain 1.7)": dict(coeffs=BASE, gain=1.7),
    "D wenig Kontrast + Rauschen": dict(coeffs=BASE, contrast=0.25, noise=0.03),
    "E elliptisch + Halo": dict(coeffs=BASE, ellipticity=0.06, halo=0.25),
    "F grosse Fehler (Zonen)": dict(coeffs={4: 0.3, 5: 0.3, 8: 0.5, 11: 0.8, 12: 0.15, 13: 0.15}),
    "G kaum Traeger (Bullseye)": dict(coeffs=BASE, fringes=(1.5, 0.5)),
}


def run(name, kw, verbose=True):
    img, truth = make_interferogram(**kw)
    gray = to_gray(img)
    r = analyze(gray, Params())
    j = len(r.coeffs)
    ct = np.zeros(j)
    for k, v in truth["coeffs"].items():
        ct[k - 1] = v
    ce = r.coeffs.copy()
    # Vorzeichen ist nicht bestimmbar -> anpassen
    s = np.sign(np.dot(ce[3:], ct[3:])) or 1.0
    ce *= s
    d = (ce - ct)[3:]
    rms_err = np.sqrt(np.sum(d ** 2))
    rms_true = np.sqrt(np.sum(ct[3:] ** 2))
    print(f"{name:32s} Streifen {r.n_fringes:5.1f} gültig {100 * r.metrics['valid_fraction']:3.0f}%  "
          f"RMS wahr {rms_true:.3f}  Fehler {rms_err:.4f} Wellen  (Vorz. {'+' if s > 0 else '-'})")
    if verbose:
        for w in r.warnings:
            print("    ! ", w)
    return r, truth, img


if __name__ == "__main__":
    for n, kw in SCEN.items():
        run(n, kw)
