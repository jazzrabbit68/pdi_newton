#!/usr/bin/env python3
"""pdi_gui.py - Tkinter-Oberflaeche fuer die PDI-Auswertung.

Start:  python pdi_gui.py [bild1.png bild2.png ...]
Benoetigt: numpy, scipy, scikit-image, opencv-python, matplotlib (Tkinter gehoert zu Python).
"""
from __future__ import annotations

import csv
import os
import sys
import time
import tkinter as tk
from tkinter import filedialog, messagebox, simpledialog, ttk

import cv2
import numpy as np

import matplotlib

matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402

from pdi_batch import analyze_series, map_from_coeffs  # noqa: E402
from pdi_core import Params, Pupil, find_pupil, load_image, pupil_from_points, to_gray, zernike_name  # noqa: E402
from pdi_synth import make_interferogram  # noqa: E402

TEST_COEFFS = {4: 0.30, 5: 0.15, 8: 0.20, 11: 0.35, 13: 0.05}


class PupilPicker:
    """Fenster zum Setzen des Spiegelrands per Mausklick (Linksklick: Punkt, Rechtsklick: letzten Punkt loeschen)."""

    def __init__(self, parent, gray, suggestion, on_accept):
        self.gray = gray
        self.suggestion = suggestion
        self.on_accept = on_accept
        self.pts: list[tuple[float, float]] = []
        self.pupil: Pupil | None = None
        self.win = tk.Toplevel(parent)
        self.win.title("Spiegelrand setzen")
        self.win.geometry("760x820")
        ttk.Label(self.win, wraplength=720, text=(
            "Linksklick: Punkt auf den Spiegelrand setzen – mindestens 3, besser 6–8 rundum verteilt "
            "(ab 6 Punkten wird eine Ellipse angepasst). Rechtsklick: letzten Punkt löschen. "
            "Die gestrichelte Linie ist der automatische Vorschlag.")).pack(fill=tk.X, padx=6, pady=4)
        self.fig = Figure(figsize=(7, 7), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.win)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self.canvas.mpl_connect("button_press_event", self.on_click)
        bar = ttk.Frame(self.win)
        bar.pack(fill=tk.X, pady=4)
        ttk.Button(bar, text="Übernehmen", command=self.accept).pack(side=tk.LEFT, padx=6)
        ttk.Button(bar, text="Automatischen Vorschlag übernehmen", command=self.accept_suggestion).pack(side=tk.LEFT, padx=6)
        ttk.Button(bar, text="Punkte löschen", command=self.reset).pack(side=tk.LEFT, padx=6)
        ttk.Button(bar, text="Abbrechen", command=self.win.destroy).pack(side=tk.RIGHT, padx=6)
        self.redraw()

    @staticmethod
    def _outline(p: Pupil):
        t = np.linspace(0, 2 * np.pi, 300)
        ca, sa = np.cos(np.deg2rad(p.angle)), np.sin(np.deg2rad(p.angle))
        ex, ey = p.d1 / 2 * np.cos(t), p.d2 / 2 * np.sin(t)
        return p.cx + ex * ca - ey * sa, p.cy + ex * sa + ey * ca

    def redraw(self):
        self.ax.clear()
        self.ax.imshow(self.gray, cmap="gray")
        if self.suggestion is not None:
            x, y = self._outline(self.suggestion)
            self.ax.plot(x, y, "y--", lw=1)
        if self.pts:
            a = np.array(self.pts)
            self.ax.plot(a[:, 0], a[:, 1], "r+", ms=12, mew=2)
        self.pupil = None
        if len(self.pts) >= 3:
            try:
                self.pupil = pupil_from_points(self.pts)
                x, y = self._outline(self.pupil)
                self.ax.plot(x, y, "c-", lw=1.2)
            except ValueError:
                pass
        self.ax.set_xlim(0, self.gray.shape[1])
        self.ax.set_ylim(self.gray.shape[0], 0)
        self.ax.set_title(f"{len(self.pts)} Punkt(e)" + ("" if self.pupil is None else
                          f" – Durchmesser {self.pupil.diameter:.0f} px"), fontsize=10)
        self.canvas.draw()

    def on_click(self, event):
        if event.inaxes is not self.ax or event.xdata is None:
            return
        if event.button == 1:
            self.pts.append((float(event.xdata), float(event.ydata)))
        elif event.button == 3 and self.pts:
            self.pts.pop()
        self.redraw()

    def reset(self):
        self.pts = []
        self.redraw()

    def accept(self):
        if self.pupil is None:
            messagebox.showinfo("Spiegelrand", "Erst mindestens 3 Punkte setzen.")
            return
        self.on_accept(self.pupil)
        self.win.destroy()

    def accept_suggestion(self):
        if self.suggestion is None:
            messagebox.showinfo("Spiegelrand", "Kein automatischer Vorschlag verfügbar.")
            return
        self.on_accept(self.suggestion)
        self.win.destroy()


class App:
    def __init__(self, root: tk.Tk):
        self.root = root
        root.title("PDI-Auswertung")
        root.geometry("1450x900")
        self.frames: list[np.ndarray] = []
        self.names: list[str] = []
        self.series = None
        self.params: Params | None = None
        self.manual_pupil: Pupil | None = None
        self.manual_sig = None
        self.v = {
            "D": tk.StringVar(value="355.6"), "F": tk.StringVar(value="1412"),
            "lam": tk.StringVar(value="650"), "K": tk.StringVar(value="-1"),
            "use_target": tk.BooleanVar(value=True),
            "mod": tk.StringVar(value="0.25"), "bw": tk.StringVar(value="0.45"),
            "nterms": tk.StringVar(value="36"), "maxdim": tk.StringVar(value="800"),
            "rm_defocus": tk.BooleanVar(value=True), "auto_sign": tk.BooleanVar(value=True),
            "flip": tk.BooleanVar(value=False), "robust": tk.BooleanVar(value=True),
            "pupil_mode": tk.StringVar(value="auto"),
            "pcx": tk.StringVar(value=""), "pcy": tk.StringVar(value=""), "pd": tk.StringVar(value=""),
        }
        self._build_ui()
        self._set_status("Bild öffnen, Video laden oder Webcam-Serie aufnehmen.")

    # ------------------------------------------------------------------ UI
    def _build_ui(self):
        left = ttk.Frame(self.root, padding=6)
        left.pack(side=tk.LEFT, fill=tk.Y)
        right = ttk.Frame(self.root)
        right.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

        box = ttk.LabelFrame(left, text="Bilder", padding=6)
        box.pack(fill=tk.X, pady=3)
        ttk.Button(box, text="Bild(er) öffnen …", command=self.open_images).pack(fill=tk.X)
        ttk.Button(box, text="Video öffnen …", command=self.open_video).pack(fill=tk.X)
        ttk.Button(box, text="Webcam-Serie aufnehmen …", command=self.capture_webcam).pack(fill=tk.X)
        ttk.Button(box, text="Synthetisches Testbild", command=self.make_test).pack(fill=tk.X)
        ttk.Button(box, text="Foucault → Testbild …", command=self.foucault_dialog).pack(fill=tk.X)
        ttk.Button(box, text="Alle entfernen", command=self.clear_frames).pack(fill=tk.X)
        self.lbl_frames = ttk.Label(box, text="0 Bilder geladen")
        self.lbl_frames.pack(anchor=tk.W, pady=(4, 0))
        self._entry(box, "Max. Bildgröße (px)", "maxdim")

        box = ttk.LabelFrame(left, text="Spiegel / Messung", padding=6)
        box.pack(fill=tk.X, pady=3)
        ttk.Checkbutton(box, text="Sollform (Paraboloid) abziehen", variable=self.v["use_target"]).pack(anchor=tk.W)
        self._entry(box, "Durchmesser (mm)", "D")
        self._entry(box, "Brennweite (mm)", "F")
        self._entry(box, "Konische Konstante K", "K")
        self._entry(box, "Wellenlänge (nm)", "lam")

        box = ttk.LabelFrame(left, text="Auswertung", padding=6)
        box.pack(fill=tk.X, pady=3)
        self._entry(box, "Kontrastschwelle (0.1–0.5)", "mod")
        self._entry(box, "Filterbreite (0.3–0.7)", "bw")
        self._entry(box, "Zernike-Terme", "nterms")
        ttk.Checkbutton(box, text="Defokus aus RMS/PV herausrechnen", variable=self.v["rm_defocus"]).pack(anchor=tk.W)
        ttk.Checkbutton(box, text="Vorzeichen automatisch (Sollform)", variable=self.v["auto_sign"]).pack(anchor=tk.W)
        ttk.Checkbutton(box, text="Vorzeichen umkehren", variable=self.v["flip"]).pack(anchor=tk.W)
        ttk.Checkbutton(box, text="Ausreißer im Fit verwerfen", variable=self.v["robust"]).pack(anchor=tk.W)

        box = ttk.LabelFrame(left, text="Pupille (Spiegelrand)", padding=6)
        box.pack(fill=tk.X, pady=3)
        ttk.Button(box, text="Rand per Mausklick setzen …", command=self.pick_pupil).pack(fill=tk.X, pady=(0, 3))
        ttk.Radiobutton(box, text="automatisch", value="auto", variable=self.v["pupil_mode"]).pack(anchor=tk.W)
        ttk.Radiobutton(box, text="manuell (Pixel)", value="manual", variable=self.v["pupil_mode"]).pack(anchor=tk.W)
        self._entry(box, "Mitte x", "pcx")
        self._entry(box, "Mitte y", "pcy")
        self._entry(box, "Durchmesser", "pd")
        ttk.Button(box, text="Erkannte Pupille übernehmen", command=self.adopt_pupil).pack(fill=tk.X, pady=(3, 0))

        ttk.Button(left, text="Auswerten", command=self.evaluate).pack(fill=tk.X, pady=(8, 2), ipady=6)
        ttk.Button(left, text="Soll/Ist-Ansicht …", command=self.open_compare).pack(fill=tk.X, pady=(2, 0))
        ttk.Button(left, text="Report speichern (PNG)", command=self.save_report).pack(fill=tk.X)
        ttk.Button(left, text="Koeffizienten speichern (CSV)", command=self.save_csv).pack(fill=tk.X)

        self.fig = Figure(figsize=(10, 6.4), dpi=100)
        self.canvas = FigureCanvasTkAgg(self.fig, master=right)
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        NavigationToolbar2Tk(self.canvas, right).update()
        self.text = tk.Text(right, height=14, font=("Courier", 10), wrap="none")
        self.text.pack(side=tk.BOTTOM, fill=tk.X)
        self.status = ttk.Label(right, text="", anchor=tk.W, relief=tk.SUNKEN)
        self.status.pack(side=tk.BOTTOM, fill=tk.X)

    def _entry(self, parent, label, key):
        row = ttk.Frame(parent)
        row.pack(fill=tk.X, pady=1)
        ttk.Label(row, text=label, width=24).pack(side=tk.LEFT)
        ttk.Entry(row, textvariable=self.v[key], width=9).pack(side=tk.RIGHT)

    def _set_status(self, msg):
        self.status.config(text=msg)
        self.root.update_idletasks()

    # ------------------------------------------------------------ Eingabe
    def _num(self, key, label):
        try:
            return float(self.v[key].get().replace(",", "."))
        except ValueError:
            raise ValueError(f"Ungültiger Wert bei „{label}“: {self.v[key].get()!r}")

    def get_params(self) -> Params:
        use_t = self.v["use_target"].get()
        return Params(
            min_mod_frac=self._num("mod", "Kontrastschwelle"), bw=self._num("bw", "Filterbreite"),
            nterms=int(self._num("nterms", "Zernike-Terme")), wavelength_nm=self._num("lam", "Wellenlänge"),
            robust=self.v["robust"].get(), flip_sign=self.v["flip"].get(),
            remove_defocus=self.v["rm_defocus"].get(), auto_sign=self.v["auto_sign"].get(),
            mirror_d_mm=self._num("D", "Durchmesser") if use_t else None,
            mirror_f_mm=self._num("F", "Brennweite") if use_t else None,
            k_target=self._num("K", "K"))

    def get_pupil(self):
        if self.v["pupil_mode"].get() != "manual":
            return None
        sig = (self.v["pcx"].get(), self.v["pcy"].get(), self.v["pd"].get())
        if self.manual_pupil is not None and sig == self.manual_sig:
            return self.manual_pupil            # exakt das per Maus/Uebernahme gesetzte (ggf. elliptische) Objekt
        d = self._num("pd", "Pupille Durchmesser")
        return Pupil(self._num("pcx", "Pupille x"), self._num("pcy", "Pupille y"), d, d, 0.0)

    def _set_manual(self, p: Pupil):
        self.v["pcx"].set(f"{p.cx:.1f}")
        self.v["pcy"].set(f"{p.cy:.1f}")
        self.v["pd"].set(f"{p.diameter:.1f}")
        self.v["pupil_mode"].set("manual")
        self.manual_pupil = p
        self.manual_sig = (self.v["pcx"].get(), self.v["pcy"].get(), self.v["pd"].get())
        self._set_status(f"Spiegelrand gesetzt: Durchmesser {p.diameter:.0f} px. „Auswerten“ drücken.")

    def pick_pupil(self):
        if not self.frames:
            messagebox.showinfo("Spiegelrand", "Erst ein Bild laden.")
            return
        gray = np.mean(self.frames, axis=0) if len(self.frames) > 1 else self.frames[0]
        try:
            sug = self.series.pupil if self.series is not None else find_pupil(np.asarray(gray, float))
        except Exception:
            sug = None
        PupilPicker(self.root, gray, sug, self._set_manual)

    def _add_frame(self, gray, name):
        try:
            maxdim = int(self._num("maxdim", "Max. Bildgröße"))
        except ValueError:
            maxdim = 0
        g = np.asarray(gray, np.float32)
        if maxdim > 0 and max(g.shape) > maxdim:
            s = maxdim / max(g.shape)
            g = cv2.resize(g, (int(round(g.shape[1] * s)), int(round(g.shape[0] * s))), interpolation=cv2.INTER_AREA)
        if self.frames and g.shape != self.frames[0].shape:
            messagebox.showwarning("Bildgröße", f"{name}: Größe {g.shape} passt nicht zu den anderen Bildern – übersprungen.")
            return
        self.frames.append(g)
        self.names.append(name)

    def _refresh_count(self):
        self.lbl_frames.config(text=f"{len(self.frames)} Bild(er) geladen")

    def open_images(self):
        paths = filedialog.askopenfilenames(title="Interferogramme wählen",
                                            filetypes=[("Bilder", "*.png *.jpg *.jpeg *.bmp *.tif *.tiff"), ("Alle", "*.*")])
        for pth in paths:
            try:
                self._add_frame(load_image(pth), os.path.basename(pth))
            except Exception as e:
                messagebox.showerror("Fehler", f"{pth}\n{e}")
        self._refresh_count()
        if paths:
            self._set_status(f"{len(paths)} Datei(en) geladen. „Auswerten“ drücken.")

    def open_video(self):
        pth = filedialog.askopenfilename(title="Video wählen", filetypes=[("Video", "*.mp4 *.avi *.mov *.mkv"), ("Alle", "*.*")])
        if not pth:
            return
        step = simpledialog.askinteger("Video", "Nur jedes n-te Bild verwenden:", initialvalue=5, minvalue=1, parent=self.root)
        if not step:
            return
        cap = cv2.VideoCapture(pth)
        i = n = 0
        while n < 40:
            ok, frame = cap.read()
            if not ok:
                break
            if i % step == 0:
                self._add_frame(to_gray(cv2.cvtColor(frame, cv2.COLOR_BGR2RGB)), f"{os.path.basename(pth)}#{i}")
                n += 1
            i += 1
        cap.release()
        self._refresh_count()
        self._set_status(f"{n} Bilder aus dem Video übernommen (max. 40).")

    def capture_webcam(self):
        n = simpledialog.askinteger("Webcam", "Anzahl Bilder in der Serie:", initialvalue=12, minvalue=1, maxvalue=60, parent=self.root)
        if not n:
            return
        idx = simpledialog.askinteger("Webcam", "Kamera-Nummer (0 = erste):", initialvalue=0, minvalue=0, parent=self.root)
        if idx is None:
            return
        cap = cv2.VideoCapture(idx)
        if not cap.isOpened():
            messagebox.showerror("Webcam", "Kamera lässt sich nicht öffnen.")
            return
        cap.set(cv2.CAP_PROP_FRAME_WIDTH, 1280)
        cap.set(cv2.CAP_PROP_FRAME_HEIGHT, 720)
        try:
            # Auto-Belichtung ausschalten (wo der Treiber es erlaubt) - sonst "pumpen" die Streifen
            cap.set(cv2.CAP_PROP_AUTO_EXPOSURE, 0.25)
        except Exception:
            pass
        for _ in range(15):          # einlaufen lassen
            cap.read()
        got = 0
        for k in range(n):
            ok, frame = cap.read()
            if not ok:
                break
            self._add_frame(to_gray(cv2.cvtColor(frame, cv2.COLOR_BGR2RGB)), f"webcam_{time.strftime('%H%M%S')}_{k}")
            got += 1
            self._set_status(f"Webcam: Bild {got}/{n}")
            self.root.update()
            time.sleep(0.15)
        cap.release()
        self._refresh_count()
        self._set_status(f"{got} Webcam-Bilder aufgenommen.")

    def make_test(self):
        img, _ = make_interferogram(size=480, coeffs=TEST_COEFFS, noise=0.02, seed=len(self.frames))
        self._add_frame(to_gray(img), f"test_{len(self.frames)}")
        self._refresh_count()
        self.v["use_target"].set(False)
        messagebox.showinfo("Testbild", "Bekannte Wellenfront (RMS-Wellen):\nj4 Defokus 0.30, j5 Astigm. 0.15,\n"
                            "j8 Koma 0.20, j11 Sphär. Aberr. 0.35, j13 0.05\n\n"
                            "Die Sollform wurde für den Test abgeschaltet. Vorzeichen ist bei einem Einzelbild beliebig.")

    def clear_frames(self):
        self.frames.clear()
        self.names.clear()
        self.series = None
        self._refresh_count()
        self.fig.clear()
        self.canvas.draw()
        self.text.delete("1.0", tk.END)

    def adopt_pupil(self):
        if self.series is None or self.series.pupil is None:
            messagebox.showinfo("Pupille", "Erst auswerten – dann kann die erkannte Pupille übernommen werden.")
            return
        self._set_manual(self.series.pupil)

    # ------------------------------------------------------------ Auswertung
    def evaluate(self):
        if not self.frames:
            messagebox.showinfo("Auswerten", "Keine Bilder geladen.")
            return
        try:
            self.params = self.get_params()
            pupil = self.get_pupil()
        except ValueError as e:
            messagebox.showerror("Eingabe", str(e))
            return

        def progress(i, n):
            self._set_status(f"Auswertung: Bild {i}/{n} …")
            self.root.update()

        self._set_status("Auswertung läuft …")
        self.root.config(cursor="watch")
        self.root.update()
        try:
            self.series = analyze_series(self.frames, self.params, pupil, progress)
        except Exception as e:
            self.root.config(cursor="")
            messagebox.showerror("Auswertung fehlgeschlagen", str(e))
            self._set_status("Fehler.")
            return
        self.root.config(cursor="")
        self.redraw()
        self._set_status(f"Fertig: {self.series.n_used} von {self.series.n_total} Bildern verwendet.")

    # ------------------------------------------------------------ Darstellung
    def _zone_profile(self, m, nbins=24):
        n = m.shape[0]
        g = np.linspace(-1, 1, n)
        X, Y = np.meshgrid(g, g)
        rho = np.hypot(X, Y)
        edges = np.linspace(0, 1, nbins + 1)
        prof, cen = [], []
        for a, b in zip(edges[:-1], edges[1:]):
            sel = (rho >= a) & (rho < b) & np.isfinite(m)
            if sel.any():
                prof.append(np.mean(m[sel]))
                cen.append(0.5 * (a + b))
        return np.array(cen), np.array(prof)

    def redraw(self):
        s, p = self.series, self.params
        r = s.best
        self.fig.clear()
        gs = self.fig.add_gridspec(2, 3, left=0.04, right=0.98, top=0.94, bottom=0.07, wspace=0.28, hspace=0.32)
        ax1, ax2, ax3 = self.fig.add_subplot(gs[0, 0]), self.fig.add_subplot(gs[0, 1]), self.fig.add_subplot(gs[0, 2])
        ax4, ax5, ax6 = self.fig.add_subplot(gs[1, 0]), self.fig.add_subplot(gs[1, 1]), self.fig.add_subplot(gs[1, 2])

        ax1.imshow(r.gray, cmap="gray")
        pu = s.pupil
        t = np.linspace(0, 2 * np.pi, 200)
        ca, sa = np.cos(np.deg2rad(pu.angle)), np.sin(np.deg2rad(pu.angle))
        ex, ey = pu.d1 / 2 * np.cos(t), pu.d2 / 2 * np.sin(t)
        ax1.plot(pu.cx + ex * ca - ey * sa, pu.cy + ex * sa + ey * ca, "c-", lw=1)
        ax1.set_title("Bild und Pupille", fontsize=9)
        ax1.axis("off")

        m = map_from_coeffs(s.coeffs, p)
        lim = max(0.05, float(np.nanmax(np.abs(m))))
        im = ax2.imshow(m, cmap="RdBu_r", vmin=-lim, vmax=lim)
        ax2.set_title(("Abweichung von der Sollform" if p.mirror_d_mm else "Wellenfront")
                      + f" (λ)\nRMS {s.metrics['rms_wf_waves']:.3f}  PV {s.metrics['pv_wf_waves']:.2f}", fontsize=9)
        ax2.axis("off")
        self.fig.colorbar(im, ax=ax2, fraction=0.046)

        v = np.where(r.mask, r.vis, np.nan)
        im = ax3.imshow(v, cmap="viridis", vmin=0, vmax=1)
        ax3.set_title(f"Streifenkontrast (gültig {100 * r.metrics['valid_fraction']:.0f} %)", fontsize=9)
        ax3.axis("off")
        self.fig.colorbar(im, ax=ax3, fraction=0.046)

        im = ax4.imshow(r.residual, cmap="RdBu_r", vmin=-0.3, vmax=0.3)
        ax4.set_title(f"Residuum Daten−Fit (RMS {r.metrics['fit_residual_rms_waves']:.3f} λ)", fontsize=9)
        ax4.axis("off")
        self.fig.colorbar(im, ax=ax4, fraction=0.046)

        js = np.arange(5, min(21, len(s.coeffs) + 1))
        ax5.bar(js, s.coeffs[js - 1], yerr=s.coeff_std[js - 1] if s.n_used > 1 else None, capsize=2)
        ax5.axhline(0, color="k", lw=0.5)
        ax5.set_xticks(js)
        ax5.tick_params(labelsize=7)
        ax5.set_title("Zernike j ≥ 5 (RMS-Wellen)", fontsize=9)

        cen, prof = self._zone_profile(m)
        ax6.plot(cen, prof * p.wavelength_nm / 2.0, "o-")
        ax6.axhline(0, color="k", lw=0.5)
        ax6.set_xlabel("Zonenradius (0 = Mitte, 1 = Rand)", fontsize=8)
        ax6.set_ylabel("Oberfläche (nm)", fontsize=8)
        ax6.tick_params(labelsize=7)
        ax6.set_title("Zonenprofil (Mittel über den Umfang)", fontsize=9)
        ax6.grid(alpha=0.3)
        self.canvas.draw()
        self._fill_text()

    # ------------------------------------------- Foucault-Eingabemaske -> Testbild
    def foucault_dialog(self):
        import pdi_foucault as pf
        try:
            D = self._num("D", "Durchmesser")
            F = self._num("F", "Brennweite")
            lam = self._num("lam", "Wellenlänge")
        except ValueError as e:
            messagebox.showerror("Foucault", str(e))
            return
        win = tk.Toplevel(self.root)
        win.title("Foucault-Messung → Testbild")
        st = {"dir": None, "entries": [], "radii": []}
        nz = tk.StringVar(value="6")
        edge = tk.BooleanVar(value=False)
        car_x, car_y = tk.StringVar(value="40"), tk.StringVar(value="8")
        size = tk.StringVar(value="800")
        top = ttk.Frame(win, padding=6)
        top.pack(fill=tk.X)
        ttk.Label(top, text=f"Spiegel D = {D:g} mm, f = {F:g} mm, λ = {lam:g} nm (aus dem Hauptfenster)").pack(anchor=tk.W)
        row = ttk.Frame(top)
        row.pack(fill=tk.X, pady=3)
        ttk.Label(row, text="Zonen").pack(side=tk.LEFT)
        ttk.Combobox(row, textvariable=nz, values=["4", "6", "8"], width=3).pack(side=tk.LEFT, padx=4)
        ttk.Checkbutton(row, text="Randzone geteilt", variable=edge).pack(side=tk.LEFT, padx=6)
        tab = ttk.Frame(win, padding=6)
        tab.pack(fill=tk.X)

        def kern_dir():
            try:
                pf.load_kern(st["dir"])
            except ImportError:
                d = filedialog.askdirectory(title="Ordner des Spiegelrechners (kern_foucault.py)")
                st["dir"] = d or None
                pf.load_kern(st["dir"])
            return st["dir"]

        def rebuild(*_):
            for w in tab.winfo_children():
                w.destroy()
            try:
                n = int(nz.get())
                st["radii"] = pf.zone_radii(D, n, edge.get(), kern_dir())
            except Exception as e:
                messagebox.showerror("Foucault", str(e))
                return
            st["entries"] = []
            ttk.Label(tab, text="Zone").grid(row=0, column=0)
            ttk.Label(tab, text="Radius (mm)").grid(row=0, column=1)
            ttk.Label(tab, text="Δl (mm)").grid(row=0, column=2)
            for i, r in enumerate(st["radii"]):
                ttk.Label(tab, text=str(i + 1)).grid(row=i + 1, column=0)
                ttk.Label(tab, text=f"{r:.1f}").grid(row=i + 1, column=1)
                v = tk.StringVar(value="0")
                ttk.Entry(tab, textvariable=v, width=9).grid(row=i + 1, column=2, padx=4, pady=1)
                st["entries"].append(v)

        def parabola():
            try:
                df = pf.parabola_readings(D, F, int(nz.get()), edge.get(), kern_dir())
            except Exception as e:
                messagebox.showerror("Foucault", str(e))
                return
            for v, x in zip(st["entries"], df):
                v.set(f"{x:.3f}")

        def readings():
            return [float(v.get().replace(",", ".")) for v in st["entries"]]

        def make(save):
            try:
                df = readings()
                fp = pf.foucault_profile(D, F, df, int(nz.get()), edge.get(), kern_dir())
                fn = pf.dummy_wavefront_fn(fp, D, F, lam)
                img, _ = make_interferogram(size=int(size.get()), fringes=(float(car_x.get().replace(",", ".")),
                                            float(car_y.get().replace(",", "."))), noise=0.01, radial_waves=fn, seed=1)
            except Exception as e:
                messagebox.showerror("Foucault", f"Eingabe oder Berechnung fehlgeschlagen:\n{e}")
                return
            name = "foucault_dummy"
            if save:
                path = filedialog.asksaveasfilename(defaultextension=".png", filetypes=[("PNG", "*.png")],
                                                    initialfile="foucault_dummy.png")
                if not path:
                    return
                cv2.imwrite(path, img)
                name = os.path.basename(path)
            self._add_frame(to_gray(img), name)
            self._refresh_count()
            self.v["use_target"].set(True)
            self._set_status(f"{name}: Foucault-Abweichung von der Parabel RMS {fp['rms']:.3f} λ (Wellenfront bei 550 nm)")

        opt = ttk.Frame(win, padding=6)
        opt.pack(fill=tk.X)
        ttk.Label(opt, text="Träger (Streifen x, y)").grid(row=0, column=0, sticky=tk.W)
        ttk.Entry(opt, textvariable=car_x, width=5).grid(row=0, column=1)
        ttk.Entry(opt, textvariable=car_y, width=5).grid(row=0, column=2, padx=3)
        ttk.Label(opt, text="Bildgröße (px)").grid(row=1, column=0, sticky=tk.W)
        ttk.Entry(opt, textvariable=size, width=7).grid(row=1, column=1, columnspan=2, sticky=tk.W)
        bt = ttk.Frame(win, padding=6)
        bt.pack(fill=tk.X)
        ttk.Button(bt, text="Zonen übernehmen", command=rebuild).pack(fill=tk.X)
        ttk.Button(bt, text="Parabel-Sollwerte einsetzen", command=parabola).pack(fill=tk.X)
        ttk.Button(bt, text="Bild erzeugen (ins Tool laden)", command=lambda: make(False)).pack(fill=tk.X, pady=(6, 0))
        ttk.Button(bt, text="Bild erzeugen und als PNG speichern …", command=lambda: make(True)).pack(fill=tk.X)
        ttk.Label(win, text="Δl: Foucault-Messwert je Zone in mm wie im Spiegelrechner (Zone 1 = 0).\n"
                            "Parabel braucht ≥ 40 Träger-Streifen, sonst wird das Bild nicht auswertbar.",
                  padding=6).pack(anchor=tk.W)
        rebuild()

    # ------------------------------------------------------- Soll/Ist-Fenster
    def open_compare(self):
        if self.series is None or self.series.best is None:
            messagebox.showinfo("Soll/Ist", "Erst auswerten.")
            return
        import pdi_compare
        win = tk.Toplevel(self.root)
        win.title("Soll/Ist-Ansicht")
        win.geometry("1300x820")
        top = ttk.Frame(win, padding=4)
        top.pack(side=tk.TOP, fill=tk.X)
        gain = tk.StringVar(value="5")
        scale = tk.StringVar(value="0.25")
        fdf = tk.StringVar(value="")
        ttk.Label(top, text="Überhöhung der Abweichung ×").pack(side=tk.LEFT)
        ttk.Combobox(top, textvariable=gain, values=["1", "2", "3", "5", "10", "20"], width=4).pack(side=tk.LEFT, padx=4)
        ttk.Label(top, text="Farbskala ± λ").pack(side=tk.LEFT, padx=(12, 0))
        ttk.Combobox(top, textvariable=scale, values=["0.05", "0.1", "0.25", "0.5", "1"], width=5).pack(side=tk.LEFT, padx=4)
        ttk.Label(top, text="Foucault df (mm, Zone 1 = 0, mit Leerzeichen):").pack(side=tk.LEFT, padx=(12, 0))
        ttk.Entry(top, textvariable=fdf, width=34).pack(side=tk.LEFT, padx=4)
        fig = Figure(figsize=(12, 7.4), dpi=100)
        canvas = FigureCanvasTkAgg(fig, master=win)
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        state = {"foucault_dir": None}

        def draw(*_):
            r, p, s = self.series.best, self.params, self.series
            try:
                g = float(gain.get().replace(",", "."))
                sc = float(scale.get().replace(",", "."))
            except ValueError:
                return
            o = pdi_compare.render(r, g)
            cr = o["crop"]
            fig.clear()
            gs = fig.add_gridspec(2, 3, left=0.03, right=0.98, top=0.93, bottom=0.07, wspace=0.30, hspace=0.22)
            for k, (t, im) in enumerate([("Aufnahme", o["real"]), ("Soll (ideale Form)", o["soll"]),
                                         (f"Ist, Abweichung ×{g:g}", o["ist"])]):
                a = fig.add_subplot(gs[0, k])
                a.imshow(im[cr], cmap="gray")
                a.set_title(t, fontsize=10)
                a.axis("off")
            a = fig.add_subplot(gs[1, 0])
            m = pdi_compare.deviation_map_waves(r)
            imh = a.imshow(m, cmap="RdBu_r", vmin=-sc, vmax=sc)
            a.set_title(f"Abweichung von der Sollform (λ), Skala ±{sc:g}", fontsize=10)
            a.axis("off")
            fig.colorbar(imh, ax=a, fraction=0.046)
            a = fig.add_subplot(gs[1, 1:])
            cen, prof = self._zone_profile(map_from_coeffs(s.coeffs, p, remove_defocus=True))
            a.plot(cen, prof * p.wavelength_nm / 2.0, "o-", label="PDI")
            msg = ""
            txt = fdf.get().replace(",", ".").split()
            if txt:
                try:
                    import pdi_foucault as pf
                    df = [float(x) for x in txt]
                    try:
                        fp = pf.foucault_profile_auto(p.mirror_d_mm, p.mirror_f_mm, df, state["foucault_dir"])
                    except ImportError:
                        d = filedialog.askdirectory(title="Ordner des Spiegelrechners (kern_foucault.py)")
                        state["foucault_dir"] = d or None
                        fp = pf.foucault_profile_auto(p.mirror_d_mm, p.mirror_f_mm, df, state["foucault_dir"])
                    c = pf.compare_profiles(cen, prof, fp, p.wavelength_nm)
                    a.plot(cen, c["foucault_nm"], "s--", label="Foucault")
                    msg = f"   Differenz PDI−Foucault: RMS {c['diff_rms_nm']:.1f} nm, PV {c['diff_pv_nm']:.1f} nm"
                except Exception as e:     # falsche Eingabe darf das Fenster nicht zerstoeren
                    msg = f"   Foucault: {e}"
            a.axhline(0, color="k", lw=0.5)
            a.set_xlabel("Zonenradius (0 = Mitte, 1 = Rand)", fontsize=9)
            a.set_ylabel("Oberfläche (nm), ohne Defokus", fontsize=9)
            a.set_title("Zonenprofil" + msg, fontsize=10)
            a.grid(alpha=0.3)
            a.legend(fontsize=8)
            fig.suptitle(f"Modell passt zur Aufnahme (Korrelation {o['corr']:.2f})", fontsize=9)
            canvas.draw()

        ttk.Button(top, text="Neu zeichnen", command=draw).pack(side=tk.LEFT, padx=8)
        draw()

    def _fill_text(self):
        s, p = self.series, self.params
        m = s.metrics
        lam = p.wavelength_nm
        lines = []
        lines.append(f"Bilder: {s.n_used} von {s.n_total} verwendet"
                     + (f"   Wiederholbarkeit: {s.repeatability_rms:.3f} λ" if np.isfinite(s.repeatability_rms) else ""))
        if s.alpha is not None:
            sd = f" ± {s.alpha_std:.3f}" if s.alpha_std is not None else ""
            lines.append(f"Korrekturgrad α = {s.alpha:.3f}{sd}   (0 = Kugel, 1 = Sollform)")
        what = "Abweichung von der Sollform" if p.mirror_d_mm else "Wellenfront"
        dnote = "ohne Kipp/Defokus" if m["defocus_removed"] else "ohne Kipp"
        lines.append(f"{what} ({dnote}):")
        pv_surf_waves = m["pv_wf_waves"] / 2.0
        lines.append(f"  Wellenfront  RMS {m['rms_wf_waves']:.3f} λ = {m['rms_wf_waves'] * lam:.0f} nm    PV {m['pv_wf_waves']:.2f} λ")
        lines.append(f"  Oberfläche   RMS {m['rms_surface_nm']:.0f} nm (λ/{(1 / (m['rms_wf_waves'] / 2.0)) if m['rms_wf_waves'] > 0 else float('inf'):.1f})"
                     f"    PV {m['pv_surface_nm']:.0f} nm (λ/{(1 / pv_surf_waves) if pv_surf_waves > 0 else float('inf'):.1f})")
        lines.append(f"  Strehl (aus Wellenfront-RMS) ≈ {m['strehl']:.2f}")
        for w in s.warnings:
            lines.append("  ! " + w)
        if s.best is not None:
            for w in s.best.warnings:
                lines.append("  ! " + w)
        lines.append("")
        if p.mirror_d_mm:
            lines.append("(Defokus j4 = Fokusversatz relativ zum paraxialen Krümmungsmittelpunkt)")
        lines.append(" j  Name                          RMS(λ)   ±      Oberfl.(nm)")
        for j in range(4, min(16, len(s.coeffs)) + 1):
            lines.append(f"{j:2d}  {zernike_name(j):28s} {s.coeffs[j - 1]:7.3f} {s.coeff_std[j - 1]:6.3f}  {s.coeffs[j - 1] * lam / 2:8.1f}")
        self.text.delete("1.0", tk.END)
        self.text.insert(tk.END, "\n".join(lines))

    # ------------------------------------------------------------ Export
    def save_report(self):
        if self.series is None:
            messagebox.showinfo("Report", "Erst auswerten.")
            return
        pth = filedialog.asksaveasfilename(defaultextension=".png", filetypes=[("PNG", "*.png")])
        if pth:
            self.fig.savefig(pth, dpi=150)
            self._set_status(f"Gespeichert: {pth}")

    def save_csv(self):
        if self.series is None:
            messagebox.showinfo("CSV", "Erst auswerten.")
            return
        pth = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV", "*.csv")])
        if not pth:
            return
        s, lam = self.series, self.params.wavelength_nm
        with open(pth, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh, delimiter=";")
            w.writerow(["j", "Name", "RMS_Wellen", "Streuung_Wellen", "Oberflaeche_nm_RMS"])
            for j in range(1, len(s.coeffs) + 1):
                w.writerow([j, zernike_name(j), f"{s.coeffs[j - 1]:.5f}", f"{s.coeff_std[j - 1]:.5f}", f"{s.coeffs[j - 1] * lam / 2:.2f}"])
        self._set_status(f"Gespeichert: {pth}")


def main():
    root = tk.Tk()
    app = App(root)
    for pth in sys.argv[1:]:
        try:
            app._add_frame(load_image(pth), os.path.basename(pth))
        except Exception as e:
            print(f"{pth}: {e}", file=sys.stderr)
    app._refresh_count()
    root.mainloop()


if __name__ == "__main__":
    main()
