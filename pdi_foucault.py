"""pdi_foucault.py - Bruecke Foucault-Messung <-> PDI-Werkzeug.

Nutzt kern_foucault.py aus dem Spiegelrechner (eine einzige Quelle der Wahrheit, keine Kopie).
Ordner mit kern_foucault.py: Umgebungsvariable PDI_FOUCAULT_DIR, Argument `folder`, oder ein Ordner
'newton_spiegel*' neben diesem Programm.

Vorzeichen/Konvention (numerisch geprueft): Foucault-W ist die Abweichung von der Parabel in Wellen
(Kugel: W = -T*rho^4, T = D^4/(512 f^3 lambda)). Das PDI-Werkzeug mit Sollform liefert dieselbe Groesse
(gemessen minus Parabel). Gemessene PDI-Wellenfront gegen die Messkugel = T*rho^4 + W_Foucault.
"""
from __future__ import annotations

import glob
import importlib
import os
import sys

import numpy as np

LAM_FOUCAULT_NM = 550.0      # wie im Spiegelrechner (fest)


def load_kern(folder: str | None = None):
    cands = [folder, os.environ.get("PDI_FOUCAULT_DIR")]
    here = os.path.dirname(os.path.abspath(__file__))
    for pat in ("newton_spiegel*", "../newton_spiegel*", "../*/newton_spiegel*"):
        cands += sorted(glob.glob(os.path.join(here, pat)))
    for c in cands:
        if c and os.path.isfile(os.path.join(c, "kern_foucault.py")):
            if c not in sys.path:
                sys.path.insert(0, c)
            return importlib.import_module("kern_foucault")
    raise ImportError("kern_foucault.py nicht gefunden - Ordner des Spiegelrechners angeben "
                      "(Argument folder oder Umgebungsvariable PDI_FOUCAULT_DIR).")


def target_waves(rho, D_mm, f_mm, lam_nm, k=-1.0):
    """Sollwellenfront Parabel gegen Messkugel (Wellen, Doppelweg) - identisch zu pdi_core."""
    r, R = D_mm / 2.0, 2.0 * f_mm
    return (-k) * (np.asarray(rho) * r) ** 4 / (4.0 * R ** 3 * lam_nm * 1e-6)


def foucault_profile(D_mm, f_mm, df_mm, n_zonen=6, rand_geteilt=False, folder=None):
    """Foucault-Zonenmessung -> Abweichung von der Parabel W_c(rho) in Wellen (550 nm),
    ohne Piston und mit bestem Defokus (wie im Spiegelrechner). df_mm: Messwerte je Zone (Zone 1 = 0)."""
    kf = load_kern(folder)
    r = kf.foucault_zonen(D_mm, n_zonen, rand_geteilt)
    g = kf.foucault_zonengrenzen(D_mm, n_zonen, rand_geteilt)
    df = np.asarray(df_mm, float)
    if len(df) != len(r):
        raise ValueError(f"{len(r)} Zonenwerte erwartet, {len(df)} bekommen.")
    W = kf.foucault_wellenfront(r, df, D_mm, f_mm, LAM_FOUCAULT_NM, best_fit_roc=False, r_grenz=g)
    kz = kf.foucault_kennzahlen(W, r, D_mm, LAM_FOUCAULT_NM, use_direct=False, grad_modus="fest", r_grenz=g)
    return {"rho": np.asarray(kz["rho_f"]), "W_c": np.asarray(kz["W_c"]), "zone_r_mm": r,
            "zone_rho": np.asarray(kz["rho_zonen"]), "W_c_zonen": np.asarray(kz["W_c_zonen"]),
            "rms": kz["W_rms"], "ptv": kz["W_ptv"], "lam_nm": LAM_FOUCAULT_NM}


def zone_radii(D_mm, n_zonen=6, rand_geteilt=False, folder=None):
    """Zonenmittelradien (mm) wie im Spiegelrechner."""
    return np.asarray(load_kern(folder).foucault_zonen(D_mm, n_zonen, rand_geteilt), float)


def parabola_readings(D_mm, f_mm, n_zonen=6, rand_geteilt=False, folder=None):
    """df-Werte (mm, Zone 1 = 0) einer perfekten Parabel: W(df) ist linear in df, daher per Ausgleichsrechnung."""
    kf = load_kern(folder)
    r = kf.foucault_zonen(D_mm, n_zonen, rand_geteilt)
    g = kf.foucault_zonengrenzen(D_mm, n_zonen, rand_geteilt)
    n_zonen = len(r)            # bei geteilter Randzone eine Zone mehr

    def W(df):
        return kf.foucault_wellenfront(r, df, D_mm, f_mm, LAM_FOUCAULT_NM, best_fit_roc=False, r_grenz=g)
    w0 = W(np.zeros(n_zonen))
    A = np.array([W(np.eye(n_zonen)[k]) - w0 for k in range(1, n_zonen)]).T
    sol, *_ = np.linalg.lstsq(A[1:], -w0[1:], rcond=None)
    return np.concatenate([[0.0], sol])


def foucault_profile_auto(D_mm, f_mm, df_mm, folder=None):
    """Wie foucault_profile; erkennt die geteilte Randzone an der Wertezahl: der Spiegelrechner bietet 4, 6 oder 8 Zonen,
    mit geteilter Randzone also 5, 7 oder 9 Werte."""
    n = len(df_mm)
    if n in (5, 7, 9):
        return foucault_profile(D_mm, f_mm, df_mm, n - 1, True, folder)
    return foucault_profile(D_mm, f_mm, df_mm, n, False, folder)


def sphere_readings(n_zonen=6):
    """Foucault-Messwerte einer perfekten Kugel (alle df = 0)."""
    return np.zeros(n_zonen)


def dummy_wavefront_fn(fp, D_mm, f_mm, lam_pdi_nm, refocus=True, k=-1.0):
    """Radiale Funktion rho -> Wellen (bei lam_pdi) fuer ein Dummy-PDI-Bild:
    Parabel gegen Messkugel + Foucault-Abweichung. refocus: Defokus T*1.5*rho^2 abziehen
    (Plaettchen leicht verschoben wie in der Praxis -> weniger dichte Randstreifen)."""
    T = float(target_waves(1.0, D_mm, f_mm, lam_pdi_nm, k))
    scale = fp["lam_nm"] / lam_pdi_nm

    def fn(rho):
        w = T * rho ** 4 + np.interp(rho, fp["rho"], fp["W_c"]) * scale
        if refocus:
            w = w - 1.5 * T * rho ** 2
        return w
    return fn


def zone_profile(coeffs, params, nbins=24):
    """Datenbasiertes Zonenprofil (Mittel ueber den Umfang) aus PDI-Zernike-Koeffizienten.
    Wellenfront in Wellen; Piston/Kipp/Defokus entfernt (Defokus wie params.remove_defocus)."""
    from pdi_batch import map_from_coeffs
    m = map_from_coeffs(coeffs, params, remove_defocus=True)
    n = m.shape[0]
    g = np.linspace(-1, 1, n)
    X, Y = np.meshgrid(g, g)
    rho = np.hypot(X, Y)
    edges = np.linspace(0, 1, nbins + 1)
    cen, prof = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        sel = (rho >= a) & (rho < b) & np.isfinite(m)
        if sel.any():
            cen.append(0.5 * (a + b))
            prof.append(float(np.mean(m[sel])))
    return np.array(cen), np.array(prof)


def compare_profiles(cen, prof_pdi_waves, fp, lam_pdi_nm):
    """PDI-Zonenprofil gegen Foucault (beide in Oberflaechen-nm, Differenz RMS/PV)."""
    wf = np.interp(cen, fp["rho"], fp["W_c"]) * fp["lam_nm"] / 2.0      # Oberflaeche nm
    wp = np.asarray(prof_pdi_waves) * lam_pdi_nm / 2.0
    wf = wf - np.mean(wf)
    wp = wp - np.mean(wp)
    d = wp - wf
    return {"foucault_nm": wf, "pdi_nm": wp, "diff_nm": d,
            "diff_rms_nm": float(np.sqrt(np.mean(d ** 2))), "diff_pv_nm": float(np.ptp(d))}


def save_profile_csv(path, cen, prof_pdi_waves, lam_pdi_nm, cmp=None):
    with open(path, "w", encoding="utf-8") as f:
        head = "rho;PDI_wellenfront_lambda;PDI_oberflaeche_nm"
        if cmp is not None:
            head += ";Foucault_oberflaeche_nm;Differenz_nm"
        f.write(head + "\n")
        for i, (c, p) in enumerate(zip(cen, prof_pdi_waves)):
            row = f"{c:.4f};{p:.5f};{p * lam_pdi_nm / 2:.2f}"
            if cmp is not None:
                row += f";{cmp['foucault_nm'][i]:.2f};{cmp['diff_nm'][i]:.2f}"
            f.write(row.replace(".", ",") + "\n")
