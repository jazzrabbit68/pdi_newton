"""foucault_dummy.py - aus Foucault-Zonenmessung ein Dummy-PDI-Bild rechnen (Rotationssymmetrie).

Beispiel (Parabel 14"/f 3.97, Zone 1 = 0, 6 Zonen, df in mm):
  python foucault_dummy.py --D 355.6 --f 1412 --df 0 0.93 1.87 2.80 3.73 4.66 --out dummy.png
Leere --df = Kugel. Danach dummy.png im PDI-Tool oeffnen, D/f eintragen, 'Sollform' an.
Optional --foucault-dir <Ordner mit kern_foucault.py>, --csv <Wellenfront-Profil des Foucault-Tests>.
"""
import argparse

import numpy as np
from PIL import Image

import pdi_foucault as pf
from pdi_synth import make_interferogram

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--D", type=float, required=True, help="Spiegeldurchmesser mm")
ap.add_argument("--f", type=float, required=True, help="Brennweite mm")
ap.add_argument("--df", type=float, nargs="*", default=None, help="Foucault-Messwerte je Zone (mm)")
ap.add_argument("--zonen", type=int, default=6)
ap.add_argument("--rand-geteilt", action="store_true")
ap.add_argument("--lam", type=float, default=650.0, help="PDI-Wellenlaenge nm")
ap.add_argument("--size", type=int, default=800)
ap.add_argument("--traeger", type=float, nargs=2, default=(40.0, 8.0), help="Kipp-Streifen x y (bei Parabel >= 40)")
ap.add_argument("--rauschen", type=float, default=0.01)
ap.add_argument("--kein-refokus", action="store_true")
ap.add_argument("--foucault-dir", default=None)
ap.add_argument("--out", default="dummy_pdi.png")
a = ap.parse_args()

df = np.zeros(a.zonen) if not a.df else np.array(a.df)
fp = pf.foucault_profile(a.D, a.f, df, a.zonen, a.rand_geteilt, a.foucault_dir)
fn = pf.dummy_wavefront_fn(fp, a.D, a.f, a.lam, refocus=not a.kein_refokus)
img, _ = make_interferogram(size=a.size, fringes=tuple(a.traeger), noise=a.rauschen, radial_waves=fn, seed=1)
Image.fromarray(img).save(a.out)
print(f"{a.out}: Foucault-Abweichung von der Parabel RMS {fp['rms']:.3f} λ(550) Wellenfront, PV {fp['ptv']:.2f} λ")
