"""pdi_compare.py - Soll/Ist-Ansicht: aus dem Auswerteergebnis Streifenbilder rechnen.

Soll  = Streifenbild der Sollform (Parabel bzw. ideale Flaeche) mit denselben Kipp-/Defokus-Anteilen wie die Aufnahme
Ist   = Streifenbild der gemessenen Wellenfront; Abweichung von der Sollform optional um 'gain' ueberhoeht,
        damit man kleine Zonenfehler auch mit dem Auge sieht.
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi

from pdi_core import Result, zernike_basis


def _phase_models(r: Result, gain: float):
    """Gibt (theta, phi_soll, phi_ist, inside) zurueck; Phasen in Radiant auf dem Bildraster."""
    h, w = r.gray.shape
    yy, xx = np.mgrid[0:h, 0:w].astype(float)
    fx0, fy0 = r.carrier
    theta = 2 * np.pi * (fx0 * (xx - r.pupil.cx) + fy0 * (yy - r.pupil.cy))
    inside = r.rho <= 1.0
    th = np.arctan2(r.y, r.x)
    nterms = len(r.c_raw)
    B = zernike_basis(nterms, r.rho[inside], th[inside])
    c_dev = np.array(r.coeffs, float)
    c_dev[:4] = 0.0                              # nur die echten Formfehler (j >= 5)
    c_soll = r.c_raw - c_dev
    c_ist = r.c_raw + (gain - 1.0) * c_dev
    phi_s = np.zeros((h, w))
    phi_i = np.zeros((h, w))
    phi_s[inside] = 2 * np.pi * (B @ c_soll)
    phi_i[inside] = 2 * np.pi * (B @ c_ist)
    return theta, phi_s, phi_i, inside


def _hp(g, inside):
    g = np.asarray(g, float)
    d = g - ndi.gaussian_filter(g, 12)
    s = np.std(d[inside])
    return d / s if s > 0 else d


def render(r: Result, gain: float = 1.0):
    """Liefert dict: soll, ist (Streifenbilder 0..1, NaN ausserhalb), real (Original), corr (Korrelation
    Modell-Ist mit Original), crop (Zeilen-/Spaltenbereich der Pupille), sign (+1/-1 verwendet)."""
    theta, phi_s, phi_i, inside = _phase_models(r, 1.0)     # gain 1 fuer Vorzeichenwahl
    real = _hp(r.gray, inside)
    best = None
    for sgn in (1.0, -1.0):
        m = np.cos(theta + sgn * phi_i)
        c = float(np.corrcoef(real[inside], _hp(0.5 + 0.5 * m, inside)[inside])[0, 1])
        if best is None or c > best[0]:
            best = (c, sgn)
    corr, sgn = best
    _, phi_s, phi_i, _ = _phase_models(r, gain)
    soll = 0.5 + 0.5 * np.cos(theta + sgn * phi_s)
    ist = 0.5 + 0.5 * np.cos(theta + sgn * phi_i)
    soll = np.where(inside, soll, np.nan)
    ist = np.where(inside, ist, np.nan)
    ys, xs = np.where(inside)
    crop = (slice(ys.min(), ys.max() + 1), slice(xs.min(), xs.max() + 1))
    return {"soll": soll, "ist": ist, "real": np.where(inside, r.gray, np.nan), "corr": corr,
            "crop": crop, "sign": sgn, "inside": inside}


def deviation_map_waves(r: Result, nd: int = 201):
    """Abweichung von der Sollform (Wellen) auf Pupillenraster (ohne Piston/Kipp/Defokus)."""
    from pdi_batch import map_from_coeffs
    from pdi_core import Params
    return map_from_coeffs(r.coeffs, Params(remove_defocus=True), n=nd)
