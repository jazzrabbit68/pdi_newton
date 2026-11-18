"""pdi_batch.py - Serienauswertung (Webcam-Video / mehrere Bilder) mit Qualitaetsbewertung und robuster Mittelung."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from pdi_core import Params, Pupil, Result, analyze, find_pupil, zernike_basis


def map_from_coeffs(c: np.ndarray, p: Params, n: int = 201, remove_defocus: bool | None = None):
    """Wellenfrontkarte (Wellen) auf einem Einheitskreis-Raster; Piston/Kipp (und optional Defokus) entfernt."""
    rd = p.remove_defocus if remove_defocus is None else remove_defocus
    g = np.linspace(-1, 1, n)
    X, Y = np.meshgrid(g, g)
    rho = np.hypot(X, Y)
    th = np.arctan2(Y, X)
    inside = rho <= 1.0
    ce = np.array(c, float)
    ce[:3] = 0.0
    if rd:
        ce[3] = 0.0
    m = np.full(rho.shape, np.nan)
    m[inside] = zernike_basis(len(ce), rho[inside], th[inside]) @ ce
    return m


def summarize(c: np.ndarray, p: Params) -> dict:
    """Kennzahlen aus einem Koeffizientenvektor (RMS-Wellen, orthonormal)."""
    ce = np.array(c, float)
    ce[:3] = 0.0
    if p.remove_defocus:
        ce[3] = 0.0
    rms = float(np.sqrt(np.sum(ce ** 2)))
    m = map_from_coeffs(c, p)
    pv = float(np.nanmax(m) - np.nanmin(m))
    lam = p.wavelength_nm
    return {
        "rms_wf_waves": rms, "pv_wf_waves": pv,
        "rms_surface_nm": rms * lam / 2.0, "pv_surface_nm": pv * lam / 2.0,
        "strehl": float(np.exp(-(2 * np.pi * rms) ** 2)),
        "defocus_removed": p.remove_defocus,
    }


@dataclass
class SeriesResult:
    coeffs: np.ndarray                 # robuster Median (wie Result.coeffs)
    coeff_std: np.ndarray              # Streuung der verwendeten Bilder
    n_total: int
    n_used: int
    used: np.ndarray                   # bool pro Bild
    frame_info: list = field(default_factory=list)   # dict pro Bild
    best: Result | None = None         # bestes Einzelbild (fuer Anzeige)
    pupil: Pupil | None = None
    metrics: dict = field(default_factory=dict)
    alpha: float | None = None
    alpha_std: float | None = None
    repeatability_rms: float = float("nan")   # typische Abweichung eines Bildes vom Median (Wellen)
    warnings: list = field(default_factory=list)


def analyze_series(grays, params: Params | None = None, pupil: Pupil | None = None, progress=None) -> SeriesResult:
    p = params or Params()
    grays = [np.asarray(g, float) for g in grays]
    n = len(grays)
    if n == 0:
        raise ValueError("Keine Bilder.")
    warns: list[str] = []
    if pupil is None:
        pupil = find_pupil(np.mean(grays, axis=0) if n > 1 else grays[0])

    results: list[Result | None] = []
    info: list[dict] = []
    for i, g in enumerate(grays):
        try:
            r = analyze(g, p, pupil)
        except Exception as e:      # ein defektes Bild darf die Serie nicht stoppen
            results.append(None)
            info.append({"index": i, "error": str(e), "ok": False})
            if progress:
                progress(i + 1, n)
            continue
        results.append(r)
        info.append({"index": i, "ok": True, "valid": r.metrics["valid_fraction"],
                     "resid": r.metrics["fit_residual_rms_waves"], "alpha": r.alpha,
                     "n_warn": len(r.warnings), "n_fringes": r.n_fringes,
                     "vis": r.metrics["vis_median"]})
        if progress:
            progress(i + 1, n)

    good = [i for i, r in enumerate(results) if r is not None]
    if not good:
        raise ValueError("Kein Bild konnte ausgewertet werden.")

    # 1) Qualitaetsfilter
    valid_f = np.array([info[i]["valid"] for i in good])
    resid = np.array([info[i]["resid"] for i in good])
    vis_m = np.array([info[i]["vis"] for i in good])
    q_ok = (valid_f >= 0.6) & (resid <= 3 * np.median(resid) + 0.01) & (vis_m >= 0.1)
    if not q_ok.any():
        q_ok[:] = True
        warns.append("Alle Bilder mit schwacher Qualität – Ergebnis unsicher.")
    cand = [i for i, ok in zip(good, q_ok) if ok]

    # 2) Vorzeichen angleichen (ohne Sollform ist es pro Bild beliebig)
    ref_i = max(cand, key=lambda i: info[i]["valid"])
    ref_c = results[ref_i].coeffs
    C = {}
    for i in cand:
        c = results[i].coeffs.copy()
        if not (p.mirror_d_mm and p.mirror_f_mm):
            if np.dot(c[4:], ref_c[4:]) < 0:
                c = -c
        C[i] = c

    # 3) Ausreisser ueber den Abstand zum Median (j >= 5; orthonormal -> euklidisch sinnvoll)
    M = np.array([C[i] for i in cand])
    med = np.median(M, axis=0)
    dist = np.sqrt(np.sum((M[:, 4:] - med[4:]) ** 2, axis=1))
    mad = 1.4826 * np.median(np.abs(dist - np.median(dist)))
    keep = dist <= np.median(dist) + 3 * mad + 0.05
    used_idx = [i for i, k in zip(cand, keep) if k]
    M = np.array([C[i] for i in used_idx])
    c_med = np.median(M, axis=0)
    c_std = np.std(M, axis=0) if len(used_idx) > 1 else np.zeros_like(c_med)
    rep = float(np.sqrt(np.mean(np.sum((M[:, 4:] - c_med[4:]) ** 2, axis=1)))) if len(used_idx) > 1 else float("nan")

    used = np.zeros(n, bool)
    used[used_idx] = True
    for i in range(n):
        info[i]["used"] = bool(used[i])
    best_i = min(used_idx, key=lambda i: info[i]["resid"])
    alphas = [results[i].alpha for i in used_idx if results[i].alpha is not None]
    if len(used_idx) < max(3, n // 2) and n > 3:
        warns.append(f"Nur {len(used_idx)} von {n} Bildern verwendet – Aufnahmebedingungen prüfen.")
    if rep > 0.1:
        warns.append(f"Große Streuung zwischen den Bildern ({rep:.2f} λ): Messung nicht reproduzierbar "
                     "(Vibration, Luftunruhe, Referenz/Loch?).")

    return SeriesResult(
        coeffs=c_med, coeff_std=c_std, n_total=n, n_used=len(used_idx), used=used, frame_info=info,
        best=results[best_i], pupil=pupil, metrics=summarize(c_med, p),
        alpha=float(np.median(alphas)) if alphas else None,
        alpha_std=float(np.std(alphas)) if len(alphas) > 1 else None,
        repeatability_rms=rep, warnings=warns)
