"""pdi_report.py - Diagnose-Grafik fuer eine Auswertung (matplotlib, auch fuer die GUI nutzbar)."""
import matplotlib
import numpy as np

from pdi_core import zernike_name


def report_figure(r, truth=None, title=""):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 4, figsize=(17, 8.5))
    ax = ax.ravel()
    ax[0].imshow(r.gray, cmap="gray"); ax[0].set_title("Bild + Pupille")
    t = np.linspace(0, 2 * np.pi, 200)
    p = r.pupil
    ca, sa = np.cos(np.deg2rad(p.angle)), np.sin(np.deg2rad(p.angle))
    ex = p.d1 / 2 * np.cos(t); ey = p.d2 / 2 * np.sin(t)
    ax[0].plot(p.cx + ex * ca - ey * sa, p.cy + ex * sa + ey * ca, "c-", lw=1)
    ax[1].imshow(np.where(r.mask, r.vis, np.nan), cmap="viridis", vmin=0, vmax=1)
    ax[1].set_title("Streifenkontrast")
    ax[2].imshow(np.where(r.valid, 1, 0) + np.where(r.fit_used, 1, 0), cmap="viridis")
    ax[2].set_title(f"gültig ({100 * r.metrics['valid_fraction']:.0f}%) / im Fit")
    ax[3].imshow(r.wrapped, cmap="twilight"); ax[3].set_title("gewickelte Phase")
    im = ax[4].imshow(r.unwrapped, cmap="RdBu_r"); ax[4].set_title("entfaltet (Wellen)")
    plt.colorbar(im, ax=ax[4], fraction=0.046)
    m = r.metrics
    im = ax[5].imshow(r.fit_map, cmap="RdBu_r")
    ax[5].set_title(f"Fit ohne Kipp{'/Defokus' if m['defocus_removed'] else ''}\nRMS {m['rms_wf_waves']:.3f} λ  PV {m['pv_wf_waves']:.2f} λ")
    plt.colorbar(im, ax=ax[5], fraction=0.046)
    im = ax[6].imshow(r.residual, cmap="RdBu_r", vmin=-0.3, vmax=0.3)
    ax[6].set_title(f"Residuum Daten-Fit (RMS {m['fit_residual_rms_waves']:.3f} λ)")
    plt.colorbar(im, ax=ax[6], fraction=0.046)
    js = np.arange(1, 16)
    ax[7].bar(js - 0.2, r.coeffs[:15], width=0.4, label="gemessen")
    if truth is not None:
        ct = np.zeros(15)
        for k, v in truth["coeffs"].items():
            if k <= 15:
                ct[k - 1] = v
        s = np.sign(np.dot(r.coeffs[3:15], ct[3:])) or 1.0
        ax[7].bar(js - 0.2, r.coeffs[:15] * (s - 1) / 2 * 0, width=0.4)
        ax[7].cla(); ax[7].bar(js - 0.2, s * r.coeffs[:15], width=0.4, label="gemessen")
        ax[7].bar(js + 0.2, ct, width=0.4, label="wahr")
    ax[7].set_xticks(js); ax[7].set_ylim(-0.6, 0.6); ax[7].legend(); ax[7].set_title("Zernike (RMS-Wellen)")
    for a in ax[:7]:
        a.axis("off")
    fig.suptitle(title)
    fig.tight_layout()
    return fig
