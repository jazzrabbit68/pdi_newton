"""Soll/Ist-Streifenbilder: Modell muss zur Aufnahme passen, Ueberhoehung nur die Abweichung aendern."""
import numpy as np
from pdi_core import Params, to_gray, analyze
from pdi_synth import make_interferogram
import pdi_compare as pc

img, _ = make_interferogram(size=600, coeffs={4: 0.3, 8: 0.05, 11: 0.1}, fringes=(24, 6), noise=0.02, seed=2)
r = analyze(to_gray(img), Params())
o1, o5 = pc.render(r, 1), pc.render(r, 5)
assert o1["corr"] > 0.9, o1["corr"]
assert np.allclose(np.nan_to_num(o1["soll"]), np.nan_to_num(o5["soll"]))        # Soll unabhaengig von der Ueberhoehung
assert np.nanmean(np.abs(o5["ist"] - o5["soll"])) > np.nanmean(np.abs(o1["ist"] - o1["soll"]))
print("test_compare ok, Korrelation %.3f" % o1["corr"])
