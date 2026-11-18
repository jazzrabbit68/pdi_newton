"""pdi_synth.py - Kuenstliche PDI-Interferogramme mit bekannter Wellenfront (zum Testen der Auswertung)."""
from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi

from pdi_core import zernike_basis


def make_interferogram(
    size: int = 480,
    coeffs: dict[int, float] | None = None,   # Noll-j -> RMS-Wellen
    fringes: tuple[float, float] = (22.0, 6.0),  # Kipp-Streifen ueber den Durchmesser (x, y)
    contrast: float = 0.8,
    falloff: float = 0.35,          # Helligkeitsabfall zum Rand (0 = flach)
    noise: float = 0.01,
    weak_zones: list[tuple[float, float, float, float]] | None = None,  # (x, y, r, faktor) normiert
    gain: float = 1.0,              # >1: Uebersteuerung
    ellipticity: float = 0.0,
    blur_px: float = 0.8,
    halo: float = 0.0,
    seed: int = 0,
    radial_waves=None,              # optional: Funktion rho -> Wellen (rotationssymmetrischer Zusatz, z.B. aus Foucault)
):
    rng = np.random.default_rng(seed)
    coeffs = coeffs or {}
    cx, cy = size / 2 + 6, size / 2 - 4
    a = 0.42 * size
    b = a * (1 - ellipticity)
    yy, xx = np.mgrid[0:size, 0:size].astype(float)
    x = (xx - cx) / a
    y = -(yy - cy) / b
    rho = np.hypot(x, y)
    th = np.arctan2(y, x)

    jmax = max(list(coeffs) + [4])
    basis = zernike_basis(jmax, rho, th)
    W = np.zeros_like(rho)
    for j, c in coeffs.items():
        W += c * basis[..., j - 1]
    if radial_waves is not None:
        W = W + radial_waves(np.minimum(rho, 1.0))
    carrier = 0.5 * (fringes[0] * x + fringes[1] * y)       # Wellen
    phase = 2 * np.pi * (W + carrier)

    V = np.full_like(rho, contrast)
    for (zx, zy, zr, f) in (weak_zones or []):
        d = np.hypot(x - zx, y - zy)
        bump = np.exp(-0.5 * (d / zr) ** 2)
        V = V * (1 - (1 - f) * bump)

    env = np.exp(-falloff * rho ** 2 * 3.0)
    pup = np.clip((1.0 - rho) * a / 1.5 + 0.5, 0, 1)         # weicher Rand ~1.5 px
    I = env * (1 + V * np.cos(phase)) / (1 + contrast) * pup
    if halo > 0:
        I = I + halo * np.exp(-0.5 * ((np.hypot(x - 0.9, y - 0.9)) / 0.5) ** 2) * (rho > 0.8)
    if blur_px > 0:
        I = ndi.gaussian_filter(I, blur_px)
    I = gain * 0.9 * I + noise * rng.standard_normal(I.shape)
    I = np.clip(I, 0, 1)
    img = np.round(I * 255).astype(np.uint8)
    truth = {"W": np.where(rho <= 1, W, np.nan), "coeffs": coeffs, "rho": rho, "cx": cx, "cy": cy, "a": a, "b": b}
    return img, truth
