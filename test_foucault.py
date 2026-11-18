"""Rundlauf: Foucault-Zonenwerte -> Dummy-PDI-Bild -> Auswertung -> Vergleich mit Foucault."""
import sys
import numpy as np
from pdi_core import Params, to_gray, analyze, find_pupil
from pdi_batch import analyze_series
from pdi_synth import make_interferogram
import pdi_foucault as pf

D, F, LAM = 355.6, 1412.0, 650.0
sys.path.insert(0, "/home/claude/spiegelrechner_v86/newton_spiegel_v8.6.0")
kf = pf.load_kern("/home/claude/spiegelrechner_v86/newton_spiegel_v8.6.0")

cases = {
    "Kugel (df=0)": np.zeros(6),
    "Parabel-nah (df = Soll)": None,
    "Parabel + Zonenfehler": None,
}
r = kf.foucault_zonen(D, 6, False)
g = kf.foucault_zonengrenzen(D, 6, False)
def Wz(df):
    return kf.foucault_wellenfront(r, df, D, F, 550.0, best_fit_roc=False, r_grenz=g)
# Messwerte einer perfekten Parabel: W(df) ist linear in df -> Zonen 2..n so waehlen, dass W = 0 (Zone 1 = 0)
A = np.array([Wz(np.eye(6)[k]) - Wz(np.zeros(6)) for k in range(1, 6)]).T
b0 = -Wz(np.zeros(6))
sol, *_ = np.linalg.lstsq(A[1:], b0[1:], rcond=None)
df_par = np.concatenate([[0.0], sol])
print("Parabel-Sollwerte df (mm):", np.round(df_par, 3), " Rest-W:", np.round(Wz(df_par), 4))
for name, df in cases.items():
    if df is None:
        df = df_par.copy()
        if "Zonenfehler" in name:
            df[3] += 0.25
            df[4] -= 0.2
    fp = pf.foucault_profile(D, F, df, 6)
    fn = pf.dummy_wavefront_fn(fp, D, F, LAM)
    img, truth = make_interferogram(size=800, fringes=(40, 8), noise=0.01, radial_waves=fn, blur_px=0.8, seed=3)
    p = Params(mirror_d_mm=D, mirror_f_mm=F, wavelength_nm=LAM)
    res = analyze_series([to_gray(img)], p)
    cen, prof = pf.zone_profile(res.coeffs, p)
    c = pf.compare_profiles(cen, prof, fp, LAM)
    print(f"{name:26s} alpha={res.alpha:.3f}  Foucault RMS={fp['rms']*275:.1f}nm(Oberfl.)  PDI RMS={res.metrics['rms_surface_nm']:.1f}nm  "
          f"Profil-Diff RMS={c['diff_rms_nm']:.1f}nm PV={c['diff_pv_nm']:.1f}nm  valid={res.best.metrics['valid_fraction']:.2f}")
