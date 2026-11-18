# PDI-Auswertung (Einzelbild-Interferogramme, Fourier-Methode)

Ausführliche Dokumentation (Prinzip, Zernike, Aufnahme, Ergebnisse lesen): **DOKU.html** bzw. DOKU.md

## Installation
    pip install numpy scipy scikit-image opencv-python matplotlib

## Start
    python pdi_gui.py                 # leer starten
    python pdi_gui.py bild1.png ...   # Bilder direkt laden

Ablauf in der Oberfläche: Bilder laden (Datei, Video oder Webcam-Serie) → Spiegeldaten eintragen
(Durchmesser, Brennweite, Wellenlänge) → "Auswerten". Mehrere Bilder werden einzeln ausgewertet,
schlechte Bilder und Ausreißer fliegen raus, der Rest wird per Median gemittelt.

## Dateien
- pdi_core.py   Auswertekern (Pupille, Fourier-Methode, Unwrapping, Zernike-Fit, Sollform)
- pdi_batch.py  Serienauswertung mit Qualitätsbewertung
- pdi_gui.py    Tkinter-Oberfläche
- pdi_synth.py  künstliche Interferogramme mit bekannter Wellenfront (zum Testen)
- test_core.py, test_batch.py, test_foucault.py   Selbsttests gegen künstliche Bilder
- pdi_compare.py   Soll/Ist-Streifenbilder mit Überhöhung (Fenster 'Soll/Ist-Ansicht')
- pdi_foucault.py, foucault_dummy.py   Brücke zum Foucaulttest (Dummy-PDI-Bild aus Zonenmessung, Profilvergleich)
- pdi_report.py Diagnose-Grafik (Zwischenschritte), z. B. für Fehlersuche

## Konventionen
- Gemessen wird die Wellenfront W = 2 × Oberflächenfehler (Doppelweg am Krümmungsmittelpunkt).
- Zernike-Koeffizienten orthonormal (Noll-Indizes), Einheit: RMS-Wellen. Oberfläche in nm = Koeffizient × λ / 2.
- Sollform: Paraboloid gegen die Messkugel (R = 2 f). Korrekturgrad α: 0 = Kugel, 1 = Paraboloid.
- Das Vorzeichen ist aus einem Einzelbild nicht bestimmbar. Mit Sollform wird es automatisch so gewählt,
  dass α positiv ist.

## Bekannte Grenzen
- Fourier-Methode braucht einen Streifenträger (Kippung, mindestens ca. 8–10 Streifen über den Durchmesser).
  Bilder mit geschlossenen Streifen im Zentrum und kaum Kippung werden unzuverlässig (Warnung erscheint).
- Pupillenerkennung kann bei hellem Halo daneben liegen: in der GUI auf "Rand per Mausklick setzen …" gehen,
  3 bis 8 Punkte rundum auf den Spiegelrand klicken (Rechtsklick löscht den letzten Punkt), "Übernehmen".
  Ab 6 Punkten wird eine Ellipse angepasst (für schräge Aufnahmen).
- Am Spiegelrand ist der Fit weniger zuverlässig (Streifen eng, Kontrast fällt ab).
