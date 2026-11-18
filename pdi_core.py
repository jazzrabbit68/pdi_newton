"""pdi_core.py - Auswertekern fuer PDI-Interferogramme (Einzelbild, Fourier-Methode).

Ablauf:
  1. Pupille (Spiegelrand) als Ellipse finden
  2. Traegerfrequenz (Kipp-Streifen) im Spektrum bestimmen
  3. Hintergrund abziehen, demodulieren, tiefpassfiltern  -> komplexes Feld z
  4. Streifenkontrast (Modulation) pro Pixel -> Gueltigkeitsmaske
  5. Phase nur ueber gueltige Pixel entfalten (Unwrapping)
  6. Zernike-Fit (orthonormal, Noll-Indizes) ueber gueltige Pixel, robust gegen Ausreisser

Konvention: Gemessen wird die Wellenfront W (doppelter Weg am Kruemmungsmittelpunkt,
also W = 2 * Oberflaechenfehler). Einheiten intern: Wellenlaengen.
Das Vorzeichen ist aus einem Einzelbild nicht eindeutig (siehe Params.flip_sign).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from math import factorial

import cv2
import numpy as np
from scipy import ndimage as ndi
from skimage.filters import threshold_otsu
from skimage.restoration import unwrap_phase

# ----------------------------------------------------------------------------
# Zernike
# ----------------------------------------------------------------------------

def noll_to_nm(j: int) -> tuple[int, int]:
    """Noll-Index j (ab 1) -> (n, m). m>0: cos-Term, m<0: sin-Term."""
    n = 0
    j1 = j - 1
    while j1 > n:
        n += 1
        j1 -= n
    m = (-1) ** j * ((n % 2) + 2 * int((j1 + ((n + 1) % 2)) / 2))
    return n, m


def _radial(n: int, m: int, rho: np.ndarray) -> np.ndarray:
    m = abs(m)
    out = np.zeros_like(rho, dtype=float)
    for k in range((n - m) // 2 + 1):
        c = ((-1) ** k * factorial(n - k)
             / (factorial(k) * factorial((n + m) // 2 - k) * factorial((n - m) // 2 - k)))
        out += c * rho ** (n - 2 * k)
    return out


def zernike_noll(j: int, rho: np.ndarray, theta: np.ndarray) -> np.ndarray:
    """Orthonormales Zernike-Polynom: RMS ueber die Einheitskreisscheibe = 1."""
    n, m = noll_to_nm(j)
    r = _radial(n, m, rho)
    if m == 0:
        return np.sqrt(n + 1.0) * r
    norm = np.sqrt(2.0 * (n + 1.0))
    return norm * r * (np.cos(m * theta) if m > 0 else np.sin(-m * theta))


ZERNIKE_NAMES = {
    1: "Piston", 2: "Kipp X", 3: "Kipp Y", 4: "Defokus",
    5: "Astigmatismus 45°", 6: "Astigmatismus 0°",
    7: "Koma Y", 8: "Koma X", 9: "Trefoil Y", 10: "Trefoil X",
    11: "Sphaerische Aberration",
    12: "Sek. Astigmatismus 0°", 13: "Sek. Astigmatismus 45°",
}


def zernike_name(j: int) -> str:
    n, m = noll_to_nm(j)
    return ZERNIKE_NAMES.get(j, f"Z{j} (n={n}, m={m})")


def zernike_basis(j_max: int, rho: np.ndarray, theta: np.ndarray) -> np.ndarray:
    return np.stack([zernike_noll(j, rho, theta) for j in range(1, j_max + 1)], axis=-1)


# ----------------------------------------------------------------------------
# Pupille
# ----------------------------------------------------------------------------

@dataclass
class Pupil:
    cx: float
    cy: float
    d1: float      # Durchmesser entlang Achse 1 (px)
    d2: float      # Durchmesser entlang Achse 2 (px)
    angle: float   # Grad (cv2.fitEllipse-Konvention)

    @property
    def diameter(self) -> float:
        return 0.5 * (self.d1 + self.d2)


def to_gray(img: np.ndarray) -> np.ndarray:
    """Beliebiges Bild -> float-Graubild 0..1. Bei Farbbildern der hellste Kanal (roter Laser -> R)."""
    arr = np.asarray(img)
    if arr.dtype == np.uint8:
        scale = 255.0
    elif arr.dtype == np.uint16:
        scale = 65535.0
    else:
        scale = 1.0
    arr = arr.astype(np.float64)
    if arr.ndim == 3:
        ch = int(np.argmax(arr.reshape(-1, arr.shape[2])[:, :3].mean(axis=0)))
        arr = arr[..., ch]
    if scale == 1.0 and arr.max() > 1.5:
        scale = float(arr.max())
    return arr / scale


def find_pupil(gray: np.ndarray) -> Pupil:
    """Spiegelrand als Ellipse: grobe Schwelle, dann Randsuche entlang von Strahlen
    (steilster Helligkeitsabfall) und robuster Ellipsen-Fit. Unempfindlich gegen dunklen Rand und Halo."""
    h, w = gray.shape
    s = max(h, w)
    # Closing fuellt dunkle Streifen, ohne den Rand aufzublasen
    win = max(5, int(0.05 * s)) | 1
    cl = ndi.grey_closing(ndi.gaussian_filter(gray, 1.5), size=(win, win))
    bw = cl > 0.5 * threshold_otsu(cl)
    bw = ndi.binary_fill_holes(bw)
    bw = ndi.binary_opening(bw, iterations=max(2, int(0.01 * s)))
    lab, n = ndi.label(bw)
    if n == 0:
        raise ValueError("Keine helle Fläche gefunden – Bild zu dunkel?")
    sizes = ndi.sum(bw, lab, range(1, n + 1))
    bw = lab == (int(np.argmax(sizes)) + 1)
    cnts, _ = cv2.findContours(bw.astype(np.uint8), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    c = max(cnts, key=cv2.contourArea)
    (cx, cy), (d1, d2), ang = cv2.fitEllipse(c)
    rough = _canonical_pupil(float(cx), float(cy), float(d1), float(d2), float(ang))

    # Randsuche entlang von Strahlen
    cl_s = ndi.gaussian_filter(cl, 1.0)
    R0 = 0.25 * (d1 + d2)
    rs = np.arange(0.6 * R0, 1.3 * R0, 0.5)
    k5 = 10
    pts = []
    for a_ in np.linspace(0, 2 * np.pi, 360, endpoint=False):
        xs = cx + rs * np.cos(a_)
        ys = cy + rs * np.sin(a_)
        prof = ndi.map_coordinates(cl_s, [ys, xs], order=1, mode="nearest")
        dp = np.gradient(prof, rs)
        k = int(np.argmin(dp))
        i0, i1 = max(0, k - k5), min(len(rs) - 1, k + k5)
        if prof[i0] > 1.5 * prof[i1] + 0.02:
            pts.append((xs[k], ys[k]))
    if len(pts) < 40:
        return rough
    pts = np.array(pts, np.float32)
    best = rough
    for _ in range(4):
        (ex, ey), (e1, e2), ea = cv2.fitEllipse(pts)
        cand = Pupil(float(ex), float(ey), float(e1), float(e2), float(ea))
        _, _, rho = _points_rho(pts, cand)
        res = (rho - 1.0) * 0.5 * cand.diameter
        mad = 1.4826 * np.median(np.abs(res - np.median(res)))
        keep = np.abs(res) < max(2.0, 3.0 * mad)
        best = cand
        if keep.all() or keep.sum() < 40:
            break
        pts = pts[keep]
    return _canonical_pupil(best.cx, best.cy, best.d1, best.d2, best.angle)


def pupil_from_points(pts) -> Pupil:
    """Pupille aus per Mausklick gesetzten Randpunkten: 3-5 Punkte -> Kreis, ab 6 Punkten -> Ellipse."""
    pts = np.asarray(pts, float).reshape(-1, 2)
    n = len(pts)
    if n < 3:
        raise ValueError("Mindestens 3 Punkte am Spiegelrand setzen.")
    if n < 6:
        A = np.column_stack([pts[:, 0], pts[:, 1], np.ones(n)])
        b = -(pts[:, 0] ** 2 + pts[:, 1] ** 2)
        sol, *_ = np.linalg.lstsq(A, b, rcond=None)
        cx, cy = -sol[0] / 2.0, -sol[1] / 2.0
        r2 = cx ** 2 + cy ** 2 - sol[2]
        if not np.isfinite(r2) or r2 <= 0:
            raise ValueError("Die Punkte liegen auf einer Geraden – bitte rundum verteilen.")
        r = float(np.sqrt(r2))
        return Pupil(float(cx), float(cy), 2 * r, 2 * r, 0.0)
    (cx, cy), (d1, d2), ang = cv2.fitEllipse(pts.astype(np.float32))
    return _canonical_pupil(float(cx), float(cy), float(d1), float(d2), float(ang))


def _points_rho(pts, p: Pupil):
    dx, dy = pts[:, 0] - p.cx, pts[:, 1] - p.cy
    t = np.deg2rad(p.angle)
    u = dx * np.cos(t) + dy * np.sin(t)
    v = -dx * np.sin(t) + dy * np.cos(t)
    xn = u / (p.d1 / 2.0)
    yn = -v / (p.d2 / 2.0)
    return xn, yn, np.hypot(xn, yn)


def _canonical_pupil(cx, cy, d1, d2, ang) -> Pupil:
    """Achsen eindeutig legen, damit die Zernike-Terme nicht mit zufaelligen Ellipsenwinkeln drehen."""
    if abs(d1 - d2) / max(d1, d2) < 0.04:        # praktisch rund -> Bildachsen
        d = 0.5 * (d1 + d2)
        return Pupil(cx, cy, d, d, 0.0)
    a = ang % 180.0
    if a > 90.0:
        a -= 180.0
    if a > 45.0:
        a -= 90.0
        d1, d2 = d2, d1
    elif a <= -45.0:
        a += 90.0
        d1, d2 = d2, d1
    return Pupil(cx, cy, d1, d2, a)


def pupil_grid(shape: tuple[int, int], p: Pupil):
    """Normierte Koordinaten x, y (Einheitskreis = Pupillenrand) und rho."""
    h, w = shape
    yy, xx = np.mgrid[0:h, 0:w].astype(float)
    dx, dy = xx - p.cx, yy - p.cy
    t = np.deg2rad(p.angle)
    u = dx * np.cos(t) + dy * np.sin(t)
    v = -dx * np.sin(t) + dy * np.cos(t)
    x = u / (p.d1 / 2.0)
    y = -v / (p.d2 / 2.0)
    return x, y, np.hypot(x, y)


# ----------------------------------------------------------------------------
# Auswertung
# ----------------------------------------------------------------------------

@dataclass
class Params:
    min_mod_frac: float = 0.25     # Streifenkontrast-Schwelle relativ zum Median
    bw: float = 0.45               # Filterbreite (1. Durchgang) relativ zur Traegerfrequenz
    iterations: int = 3            # Demodulations-Durchgaenge
    nterms: int = 36               # Zernike-Terme (Noll j = 1..nterms)
    shrink: float = 0.97           # Maske etwas kleiner als der erkannte Rand
    wavelength_nm: float = 650.0   # Laserdiode
    robust: bool = True            # Ausreisser im Fit verwerfen
    flip_sign: bool = False        # Vorzeichen der Wellenfront umkehren
    remove_defocus: bool = True    # RMS/PV ohne Defokus angeben
    mirror_d_mm: float | None = None   # Spiegeldurchmesser -> Sollform (Paraboloid) am Kruemmungsmittelpunkt
    mirror_f_mm: float | None = None   # Brennweite (R = 2 f)
    k_target: float = -1.0             # Konstante der Sollflaeche (-1 = Paraboloid)
    auto_sign: bool = True             # Vorzeichen so waehlen, dass der Korrekturgrad positiv wird


@dataclass
class Result:
    gray: np.ndarray
    pupil: Pupil
    x: np.ndarray
    y: np.ndarray
    rho: np.ndarray
    mask: np.ndarray            # Pupille (shrink)
    valid: np.ndarray           # genug Streifenkontrast, zusammenhaengend
    fit_used: np.ndarray        # tatsaechlich im Fit verwendet
    vis: np.ndarray             # Streifenkontrast (Visibility)
    wrapped: np.ndarray         # gewickelte Phase (rad), ausserhalb nan
    unwrapped: np.ndarray       # entfaltete Phase (Wellen), ausserhalb nan
    fit_map: np.ndarray         # Zernike-Fit (Wellen) ueber ganze Pupille, ohne Piston/Kipp
    residual: np.ndarray        # Daten - Fit (Wellen) auf gueltigen Pixeln, sonst nan
    coeffs: np.ndarray          # orthonormal (RMS-Wellen), Index 0 = j=1
    carrier: tuple[float, float]
    n_fringes: float
    c_raw: np.ndarray | None = None      # gemessene Koeffizienten (ohne Sollform-Abzug)
    target_map: np.ndarray | None = None  # Sollwellenfront (Wellen), falls D und f bekannt
    alpha: float | None = None            # Korrekturgrad: 0 = Kugel, 1 = Sollform
    metrics: dict = field(default_factory=dict)
    warnings: list = field(default_factory=list)


def _nconv(img, mask, sigma):
    m = mask.astype(float)
    num = ndi.gaussian_filter(img * m, sigma, mode="constant")
    den = ndi.gaussian_filter(m, sigma, mode="constant")
    return num / np.maximum(den, 1e-6)


def _estimate_carrier(h, mask, diameter):
    """Traegerfrequenz (cycles/px) als staerkster Peak im Spektrum (Halbebene)."""
    hw = h * ndi.gaussian_filter(mask.astype(float), 2.0)
    N = 1
    while N < 2 * max(h.shape):
        N *= 2
    F = np.fft.fftshift(np.fft.fft2(hw, s=(N, N)))
    P = ndi.gaussian_filter(np.abs(F), 1.0)
    f = np.fft.fftshift(np.fft.fftfreq(N))
    FX, FY = np.meshgrid(f, f)
    fr = np.hypot(FX, FY)
    fmin = 2.5 / diameter
    half = (FX > 0) | ((FX == 0) & (FY > 0))
    sel = (fr > fmin) & (fr < 0.4) & half
    idx = int(np.argmax(np.where(sel, P, 0.0)))
    iy, ix = np.unravel_index(idx, P.shape)
    fx0, fy0 = f[ix], f[iy]
    if 1 <= ix < N - 1 and 1 <= iy < N - 1:
        def para(a, b, c):
            d = a - 2 * b + c
            return 0.0 if abs(d) < 1e-12 else 0.5 * (a - c) / d
        fx0 += para(P[iy, ix - 1], P[iy, ix], P[iy, ix + 1]) * (f[1] - f[0])
        fy0 += para(P[iy - 1, ix], P[iy, ix], P[iy + 1, ix]) * (f[1] - f[0])
    return float(fx0), float(fy0), fmin


def _largest_component(m):
    lab, n = ndi.label(m)
    if n <= 1:
        return m
    sizes = ndi.sum(m, lab, range(1, n + 1))
    return lab == (int(np.argmax(sizes)) + 1)


def fit_zernike(x, y, phi, valid, nterms, robust=True):
    rho = np.hypot(x, y)
    th = np.arctan2(y, x)
    A = zernike_basis(nterms, rho[valid], th[valid])
    b = phi[valid]
    keep = np.ones(b.shape, bool)
    c = np.zeros(nterms)
    for _ in range(5):
        c, *_ = np.linalg.lstsq(A[keep], b[keep], rcond=None)
        if not robust:
            break
        res = b - A @ c
        r0 = res[keep]
        s = 1.4826 * np.median(np.abs(r0 - np.median(r0)))
        new = np.abs(res) < max(0.3, 4.0 * s)
        if new.sum() == keep.sum() or new.sum() < 0.5 * b.size:
            break
        keep = new
    used = np.zeros(valid.shape, bool)
    used[np.where(valid)] = keep
    cond = float(np.sqrt(np.linalg.cond(A[keep].T @ A[keep])))
    return c, used, cond


def analyze(gray: np.ndarray, params: Params | None = None, pupil: Pupil | None = None) -> Result:
    p = params or Params()
    warns: list[str] = []
    gray = np.asarray(gray, float)
    if pupil is None:
        pupil = find_pupil(gray)
    x, y, rho = pupil_grid(gray.shape, pupil)
    mask = rho <= p.shrink
    D = pupil.diameter

    # 1) grober Hintergrund -> Traegerfrequenz
    bg0 = _nconv(gray, mask, 0.1 * D)
    fx0, fy0, fmin = _estimate_carrier((gray - bg0), mask, D)
    f0 = float(np.hypot(fx0, fy0))
    n_fringes = f0 * D
    if f0 < 1.3 * fmin:
        warns.append(f"Kaum Streifenträger ({n_fringes:.1f} Streifen über den Durchmesser): "
                     "Fourier-Methode unsicher – mehr Kippung einstellen.")
    if f0 > 0.25:
        warns.append(f"Streifen sehr eng ({1 / f0:.1f} px Periode): Abtastgrenze, höhere Auflösung/weniger Kippung wählen.")
    period = 1.0 / f0

    # 2) feiner Hintergrund
    bg = _nconv(gray, mask, 1.2 * period)
    h = (gray - bg) * mask
    hh, ww = gray.shape
    yy, xx = np.mgrid[0:hh, 0:ww].astype(float)
    theta = 2 * np.pi * (fx0 * (xx - pupil.cx) + fy0 * (yy - pupil.cy))
    N = 1
    while N < 2 * max(hh, ww):
        N *= 2
    fax = np.fft.fftfreq(N)
    FX, FY = np.meshgrid(fax, fax)
    FR2 = FX ** 2 + FY ** 2
    ang = np.arctan2(y, x)
    ev = rho <= 1.05
    basis_ev = zernike_basis(p.nterms, rho[ev], ang[ev])
    bg_med = np.median(bg[mask])
    core = ndi.binary_erosion(mask, iterations=max(2, int(0.08 * D))) & (bg > 0.2 * bg_med)
    if core.sum() < 100:
        core = mask

    # 3) iterative Demodulation: Traeger + bisher geschaetzte Wellenfront abziehen, schmal filtern
    model = np.zeros(gray.shape)          # Wellen
    c = np.zeros(p.nterms)
    wrapped0 = None
    for it in range(max(1, p.iterations)):
        sig = max((p.bw if it == 0 else 0.5 * p.bw) * f0, 1.5 / D)
        hd = h * np.exp(-1j * (theta + 2 * np.pi * model))
        Hd = np.fft.fft2(hd, s=(N, N))
        z = np.fft.ifft2(Hd * np.exp(-FR2 / (2 * sig ** 2)))[:hh, :ww]
        vis = 2.0 * np.abs(z) / np.maximum(bg, 1e-3)
        vmed = float(np.median(vis[core]))
        valid = mask & (vis > p.min_mod_frac * vmed) & (bg > 0.1 * bg_med)
        valid = ndi.binary_opening(valid, structure=np.ones((3, 3)))
        valid = _largest_component(valid)
        wrapped = np.angle(z)
        if wrapped0 is None:
            wrapped0 = wrapped
        um = unwrap_phase(np.ma.array(wrapped, mask=~valid))
        d = np.ma.filled(um, np.nan) / (2 * np.pi)
        valid = valid & np.isfinite(d)
        ck, used, cond = fit_zernike(x, y, np.nan_to_num(d), valid, p.nterms, p.robust)
        phi = np.where(valid, model + d, np.nan)       # Daten (Wellen, relativ zum Traeger)
        model = model.copy()
        model[ev] += basis_ev @ ck
        c = c + ck
    frac_valid = valid.sum() / max(mask.sum(), 1)
    if vmed < 0.2:
        warns.append(f"Streifenkontrast insgesamt gering (Visibility {vmed:.2f}): Belichtung, Referenz-/Messstrahl-Intensität "
                     "und Streulicht prüfen.")
    if frac_valid < 0.6:
        warns.append(f"Nur {100 * frac_valid:.0f}% der Pupille hat brauchbaren Streifenkontrast.")
    if p.flip_sign:
        c = -c
        phi = -phi
    if cond > 1e3:
        warns.append(f"Zernike-Fit schlecht konditioniert ({cond:.0f}): zu große Lücken in der gültigen Fläche, "
                     "Ergebnis (v.a. Randterme) unsicher.")
    n_trim = int(valid.sum() - used.sum())
    if n_trim > 0.05 * valid.sum():
        warns.append(f"{100 * n_trim / valid.sum():.0f}% der Pixel als Ausreißer verworfen (Unwrapping-Fehler?).")
    wrapped = wrapped0

    th = np.arctan2(y, x)
    full = rho <= 1.0
    basis = zernike_basis(p.nterms, rho[full], th[full])
    c_raw = c.copy()
    alpha = None
    target_map = None
    if p.mirror_d_mm and p.mirror_f_mm:
        # Sollwellenfront gegen die Messkugel (R = 2 f), Doppelweg: T = -K * (rho*r)^4 / (4 R^3 lambda)
        r_mm = p.mirror_d_mm / 2.0
        R_mm = 2.0 * p.mirror_f_mm
        lam_mm = p.wavelength_nm * 1e-6
        t_full = (-p.k_target) * (rho[full] * r_mm) ** 4 / (4.0 * R_mm ** 3 * lam_mm)
        cT, *_ = np.linalg.lstsq(basis, t_full, rcond=None)
        den = float(np.dot(cT[4:], cT[4:]))
        alpha = float(np.dot(c[4:], cT[4:]) / den) if den > 0 else None
        if alpha is not None and alpha < 0 and p.auto_sign:
            c = -c
            phi = -phi
            alpha = -alpha
        c_raw = c.copy()
        c = c - cT                       # Abweichung von der Sollform
        target_map = np.full(gray.shape, np.nan)
        target_map[full] = basis @ cT
        if alpha is not None and not (0.5 < alpha < 1.5):
            warns.append(f"Korrekturgrad {alpha:.2f} ist ungewöhnlich (erwartet 0 bis ~1): "
                         "Durchmesser/Brennweite prüfen oder Auswertung unsicher.")
    c_eval = c.copy()
    c_eval[:3] = 0.0                         # Piston + Kippung nie bewerten
    fit_map = np.full(gray.shape, np.nan)
    fit_map[full] = basis @ c_eval
    eval_c = c_eval.copy()
    if p.remove_defocus:
        eval_c[3] = 0.0
    eval_map = np.full(gray.shape, np.nan)
    eval_map[full] = basis @ eval_c

    resid = np.full(gray.shape, np.nan)
    fitfull = np.full(gray.shape, np.nan)
    fitfull[full] = basis @ c_raw
    resid[valid] = phi[valid] - fitfull[valid]

    rms_wf = float(np.sqrt(np.sum(eval_c ** 2)))
    pv_wf = float(np.nanmax(eval_map) - np.nanmin(eval_map))
    lam = p.wavelength_nm
    metrics = {
        "rms_wf_waves": rms_wf,
        "pv_wf_waves": pv_wf,
        "rms_surface_nm": rms_wf * lam / 2.0,
        "pv_surface_nm": pv_wf * lam / 2.0,
        "strehl": float(np.exp(-(2 * np.pi * rms_wf) ** 2)),
        "fit_residual_rms_waves": float(np.nanstd(resid[used])) if used.any() else float("nan"),
        "valid_fraction": float(frac_valid),
        "n_fringes": float(n_fringes),
        "cond": cond,
        "vis_median": vmed,
        "defocus_removed": p.remove_defocus,
        "alpha": alpha,
        "vs_target": target_map is not None,
    }
    return Result(gray=gray, pupil=pupil, x=x, y=y, rho=rho, mask=mask, valid=valid, fit_used=used,
                  vis=vis, wrapped=np.where(valid, wrapped, np.nan),
                  unwrapped=np.where(valid, phi, np.nan), fit_map=fit_map, residual=resid,
                  coeffs=c, carrier=(fx0, fy0), n_fringes=n_fringes, metrics=metrics, warnings=warns,
                  c_raw=c_raw, target_map=target_map, alpha=alpha)


def load_image(path: str) -> np.ndarray:
    """Bilddatei -> Graubild (float 0..1)."""
    im = cv2.imread(path, cv2.IMREAD_UNCHANGED)
    if im is None:
        raise ValueError(f"Kann Bild nicht lesen: {path}")
    if im.ndim == 3:
        im = cv2.cvtColor(im[..., :3], cv2.COLOR_BGR2RGB)
    return to_gray(im)
