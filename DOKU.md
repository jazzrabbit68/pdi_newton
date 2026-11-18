# PDI-Auswertung – Dokumentation

Eigenbau-Auswertung für Einzelbild-Interferogramme von Teleskopspiegeln (Point-Diffraction-Interferometer, PDI).
Gedacht für einen 14"-Spiegel um f/4, Messung am Krümmungsmittelpunkt, Laserdiode, Webcam oder Kamera mit Objektiv.

---

## 1. Das Prinzip des PDI

![Strahlengang](doku/strahlengang.png)

Der Spiegel erzeugt vom Laserpunkt ein Sternbildchen. Genau dort sitzt das **PDI-Plättchen**. Es ist **teildurchlässig**
(abschwächende Schicht) und hat ein **winziges Loch** (wenige µm):

- Das Licht, das durch das Plättchen geht, trägt noch alle Fehler des Spiegels. Das ist der **Messstrahl**.
- Das Licht, das durch das Loch geht, wird gebeugt und kommt als **saubere Kugelwelle** heraus. Das ist die **Referenz**,
  also ein „idealer Spiegel“, den man sich selbst baut.

Beide Strahlen überlagern sich auf der Kamera. Wo Berg auf Berg trifft, ist es hell, wo Berg auf Tal trifft, dunkel.
Das ergibt die Streifen. Zwischen zwei benachbarten Streifen ändert sich die Wellenfront um eine Wellenlänge,
die **Spiegeloberfläche** also um eine halbe (das Licht läuft hin und zurück).

Ein perfekter Spiegel gegen die Idealkugel gibt gerade Linien oder Ringe. Ein Spiegel mit Fehlern gibt verbogene Streifen,
und die Verbiegung ist die Information, die das Tool auswertet.

**Warum das Plättchen teildurchlässig sein muss:** Das Loch lässt nur einen winzigen Bruchteil des Lichts durch. Damit
Referenz und Messstrahl ähnlich hell sind (nur dann haben die Streifen guten Kontrast), wird der Messstrahl durch die
Schicht abgeschwächt. Der Kontrast hängt deshalb an der Durchlässigkeit des Plättchens, der Lochgröße und der Größe des
Sternbildchens, und das ändert sich mit dem Spiegelzustand.

**Warum die Messung zickig sein kann:** Ist das Loch zu groß, trägt die „Referenz“ Spiegelfehler mit und die Messung zeigt
weniger Fehler als vorhanden. Ist es zu klein, kommt zu wenig Licht durch und der Kontrast schwankt. Ein einzelnes schönes
Interferogramm beweist deshalb nicht, dass die Referenz sauber war. Gegenprobe: Spiegel um 90° und 180° drehen. Fehler,
die mitwandern, gehören zum Spiegel; Fehler, die im Bild stehen bleiben, kommen vom Aufbau.

## 2. Was Zernike damit zu tun hat

![Zernike-Tafel](doku/zernike_tafel.png)

Nach der Auswertung liegt eine Höhenkarte der Wellenfront vor, also viele tausend Zahlen. Die sagen noch nicht, was am
Spiegel falsch ist. **Zernike-Polynome** sind Standardformen, aus denen sich jede Höhenkarte zusammensetzen lässt, so wie sich
jede Farbe aus Rot, Grün und Blau mischen lässt. Die Formen haben bekannte Namen:

| Form | Bedeutung am Spiegel |
|---|---|
| Defokus | Fokusabstand/Messabstand falsch (Aufbau, kein Spiegelfehler) |
| Astigmatismus | in einer Richtung stärker gekrümmt: Auflage, Spannung, Politur |
| Koma | einseitig verzerrt: Spiegel schief im Aufbau oder Kantenfehler |
| Sphärische Aberration | Mitte und Rand haben verschiedene Brennweite: Parabel nicht richtig |
| Dreiblatt (Trefoil) | dreieckig: Spannung durch die Auflage |

Das Tool passt diese Formen an die Messung an und gibt pro Form eine Zahl (RMS-Wellen). Aus tausenden Messwerten werden
so etwa zehn Zahlen, die sagen, wo nachgearbeitet werden muss. Kippung und Defokus werden nicht bewertet, sie sind
Einstellung des Aufbaus. Gesamt-RMS und Strehl-Wert lassen sich direkt aus den Zahlen berechnen.

### 2.1 Woher der Name kommt

Die Polynome stammen von dem niederländischen Physiker **Frits Zernike** (1888–1966). Er entwickelte sie 1934 bei der Theorie
des **Foucaulttests**: Er wollte beschreiben, wie sich Wellenfrontfehler im Bild der Messerschneide zeigen, und brauchte
dafür Formen, die auf einer Kreisfläche (der Spiegelöffnung) unabhängig voneinander sind. Bekannter wurde er durch den
Phasenkontrast in der Mikroskopie (Nobelpreis Physik 1953).

Die Durchnummerierung j = 1, 2, 3 … im Tool stammt von Robert Noll (1976), deshalb „Noll-Indizes“. Andere Programme zählen die
Formen teils anders durch und normieren sie anders. Beim Vergleich mit Werten aus anderer Software vorher Zählweise und
Normierung prüfen. Unsere Koeffizienten sind **orthonormal**, ihre Einheit ist RMS-Wellen.

### 2.2 Der Fit ist global

Die Zernike-Formen werden **nicht** auf kleine Spiegelbereiche losgelassen. Jede Form liegt über der **ganzen** Spiegelfläche.
Das Tool sucht die eine Mischung aus (standardmäßig) 36 Formen, die zur gemessenen Wellenfront am besten passt
(Methode der kleinsten Fehlerquadrate). Jedes gültige Pixel liefert eine Gleichung, bei einem Bild sind das einige zehntausend
Gleichungen für 36 Unbekannte.

- **Vorteil:** Rauschen und einzelne schlechte Pixel mitteln sich heraus, und jede Zahl gilt für den ganzen Spiegel.
- **Nachteil:** Lokale Fehler wie Zonen, eine abgesunkene Kante oder ein Kratzer sind für 36 glatte Formen zu kleinteilig.
  Sie landen im **Residuum** (Daten minus Fit), nicht in den Koeffizienten.

Faustregel: Für großräumige Fehler (Astigmatismus, Koma, Parabel) sind die Zernike-Zahlen das richtige Werkzeug. Für lokale
Fehler (Zonen, Kante) schaut man auf **Residuum und Zonenprofil**, ähnlich wie beim Foucaulttest.

### 2.3 Formen von Hand zu erkennen ist schwer, deshalb der Fit

![Zerlegung eines Fehlergemischs](doku/zerlegung.png)

Treten mehrere Fehler gleichzeitig auf, was bei echten Spiegeln der Normalfall ist, überlagern sich ihre Streifenbilder und die
Einzelformen sind von Hand kaum noch zu erkennen (zweites Bild oben, ohne Kippung). Der Fit trennt sie wieder. Im Beispiel
(künstliches Bild mit vier bekannten Fehlern) findet das Tool Astigmatismus 0,22 (wahr 0,22), Koma 0,17 (0,18),
sphärische Aberration 0,27 (0,28) und Dreiblatt 0,10 (0,10), und das Residuum ist praktisch null. Bei echten Bildern ist
das Residuum nie null, und genau dort sieht man, was die Formen nicht erklären.

## 3. Aufbau und Aufnahme

### 3.1 Kamera, Objektiv und Abstände

Das Objektiv wird **auf den Spiegel scharfgestellt**, nicht auf das Plättchen. Das Bild auf dem Sensor ist ein Bild des
Spiegels, und die Streifenlage wird dem Ort auf dem Spiegel zugeordnet. Die Entfernung bis zum Spiegel ist etwa der
Krümmungsradius R = 2·f (bei f = 1412 mm also rund 2,8 m).

**Wichtigster Abstand: Plättchen → Objektiv.** Hinter dem Plättchen läuft das Licht als Kegel auseinander. Der volle
Öffnungswinkel ist D/R, also etwa 0,126 rad (7,2°) bei f/3,97. Der Strahldurchmesser am Objektiv ist

    Strahldurchmesser = 0,126 × Abstand Plättchen–Objektiv      (für einen f/4-Spiegel)

Die Eintrittspupille des Objektivs (Brennweite / Blendenzahl) muss mindestens so groß sein. Sonst werden die Ränder
des Spiegelbildes abgeschnitten (Randabdunkelung, schwacher Kontrast am Rand, unsicherer Spiegelrand).

| Abstand Plättchen–Objektiv | Strahldurchmesser |
|---|---|
| 50 mm | 6,3 mm |
| 100 mm | 12,6 mm |
| 150 mm | 18,9 mm |
| 300 mm | 37,8 mm |

Beispiele (Eintrittspupille = Brennweite / Blende → größter zulässiger Abstand): 135 mm f/8 → 17 mm → bis etwa 134 mm;
135 mm f/11 → 12 mm → bis etwa 97 mm; Webcam-Objektiv 4 mm f/2,8 → 1,4 mm → nur etwa 11 mm.
Also: Objektiv möglichst nah ans Plättchen und lieber etwas weiter aufblenden als zu stark abblenden.

**Der Abstand Plättchen → Sensor ist keine eigene Stellgröße.** Er ergibt sich aus Objektiv und Scharfstellung.
Auf die Streifen selbst hat der Abstand Plättchen–Objektiv keinen direkten Einfluss, denn Mess- und Referenzstrahl
kommen aus demselben Punkt. Die Streifen werden verändert durch:

- **Plättchen seitlich verschieben:** ändert die Kippung (Streifenzahl und Richtung).
- **Plättchen in Achsrichtung verschieben:** ändert Defokus (Ringe).

**Bildmaßstab:** Soll der Spiegel den Sensor füllen, gilt grob  Brennweite ≈ m·u/(1+m)  mit m = Sensorhöhe / Spiegeldurchmesser
und u ≈ Abstand zum Spiegel. Beispiel: APS-C (15,7 mm Höhe) und 356-mm-Spiegel in 2,9 m Entfernung → ca. 125 mm.

**Verzeichnung:** Starke Verzeichnung (typisch bei Weitwinkel-Objektiven, kleinen Webcams) verzerrt die Zuordnung Bild →
Spiegelort und fälscht damit v. a. die rotationssymmetrischen Terme. Ein Teleobjektiv um 100–135 mm ist unkritisch.

### 3.2 Aufnahmeregeln

- **Format:** PNG oder RAW, kein JPEG. Bei Farbkamera und rotem Laser nutzt das Tool automatisch den hellsten (roten) Kanal.
- **Belichtung:** Nichts übersteuern. Auto-Belichtung der Webcam ausschalten, sonst „pumpen“ die Streifen.
- **Streifenträger:** Etwa 8 bis 10 Streifen über den Durchmesser einstellen (Kippung). Bilder mit geschlossenen Ringen
  und kaum Kippung kann die Fourier-Methode nicht auswerten.
- **Nicht zu eng:** Streifenperiode am Rand mindestens etwa 4 Pixel (Warnung des Tools bei Überschreitung).
- **Serie:** 10 bis 20 Bilder vom selben Zustand. Das Tool wertet einzeln aus, sortiert schlechte Bilder und Ausreißer aus
  und mittelt den Rest (Median).
- **Nicht mischen:** Verschiedene Löcher oder Polierstände einzeln auswerten.
- **Scharf auf den Spiegelrand:** Ein unscharfer Rand macht die Pupille ungenau, siehe 5.3.

## 4. Bedienung

Installation:

    pip install numpy scipy scikit-image opencv-python matplotlib

Start:

    python pdi_gui.py                 # leer starten
    python pdi_gui.py bild1.png ...   # Bilder direkt laden

Ablauf:

1. Bilder laden: Datei(en), Video oder Webcam-Serie.
2. Spiegeldaten eintragen: Durchmesser, Brennweite, Wellenlänge der Laserdiode.
3. **Spiegelrand:** „Rand per Mausklick setzen …“, 3 bis 8 Punkte rundum auf den Rand klicken (Rechtsklick löscht den
   letzten), „Übernehmen“. Ab 6 Punkten wird eine Ellipse angepasst (schräge Aufnahme). Der gelbe gestrichelte Kreis ist der
   automatische Vorschlag.
4. „Auswerten“.

Die Oberfläche zeigt Bild mit Pupille, Wellenfront (oder Abweichung von der Sollform), Streifenkontrast, Residuum, Zernike-Balken
und das Zonenprofil. Darunter stehen die Kennzahlen mit Warnungen. Export: Report als PNG, Koeffizienten als CSV.

![Beispielauswertung](doku/beispiel_auswertung.png)

### Soll/Ist-Ansicht

Der Knopf **Soll/Ist-Ansicht …** (nach dem Auswerten) öffnet ein eigenes Fenster:

- **Aufnahme, Soll, Ist:** Das Tool rechnet aus der Auswertung Streifenbilder: Soll = ideale Form (Parabel) mit denselben Kipp- und
  Defokusanteilen wie die Aufnahme, Ist = gemessene Form. Die **Überhöhung ×N** verstärkt nur die Abweichung von der Sollform,
  sodass auch kleine Zonenfehler als verbogene Streifen sichtbar werden. Bei ×1 ist das Ist-Bild das, was der Fit der Aufnahme
  entspricht (Kontrolle: oben im Fenster steht die Korrelation Modell/Aufnahme, typisch über 0,9).
- **Abweichungskarte** mit wählbarer Farbskala (±0,05 bis ±1 λ).
- **Zonenprofil** mit optionaler Foucault-Überlagerung: df-Werte der Zonen (mm, Zone 1 = 0) eintragen, D und f stehen schon im Tool.
  Der Ordner mit `kern_foucault.py` wird gesucht oder per Dialog erfragt.

Beispiel: Ein um 1 mm zu langer Foucault-Randzonenwert ist im Rohbild kaum zu sehen, bei ×5 aber deutlich, und im Zonenprofil
stimmen PDI und Foucault auf etwa 5 nm überein.

## 5. Ergebnisse richtig lesen

### 5.1 Kennzahlen

- **RMS** und **PV** der Wellenfront (Kipp und Defokus herausgerechnet). Oberfläche in nm = Wellenfront / 2.
- **Strehl:** aus der Wellenfront-RMS nach exp(−(2π·RMS)²). Das ist nur für kleine Fehler (unter etwa 0,1 λ RMS) brauchbar;
  darüber nur ein grober Hinweis. Das Tool rechnet mit der Wellenlänge der Laserdiode. Fürs Visuelle (550 nm) wird der Wert
  schlechter.
- **Korrekturgrad α:** 0 = Kugel, 1 = Sollform (Paraboloid). Die Abweichung von der Sollform wird bei der Parabel über
  die konische Konstante K = −1 gerechnet.
- **Residuum:** Daten minus Zernike-Fit. Hohe Werte weisen auf Unwrapping-Fehler oder lokale Strukturen hin, die der Fit nicht
  abbildet.

### 5.2 Vorzeichen

Aus einem Einzelbild ist nicht erkennbar, ob eine Stelle Berg oder Tal ist. Mit Sollform wird das Vorzeichen so gewählt,
dass α positiv wird. Ob „Hügel“ der Richtung entspricht, die man vom Foucaulttest kennt, sollte man einmal an einem bekannten
Spiegel prüfen.

### 5.3 Empfindlichkeit gegenüber dem Spiegelrand (wichtig)

Die Sollform einer Parabel weicht am Krümmungsmittelpunkt bei 14"/f4 um rund 17 Wellen von der Kugel ab und wächst mit der
vierten Potenz des Radius. Deshalb geht ein kleiner Fehler im Pupillenradius stark in die sphärische Aberration ein:

| Pupillenradius gegenüber Erkennung | Korrekturgrad α | Sphär. Aberration (j11) | RMS (λ) | Strehl bei 650 nm |
|---|---|---|---|---|
| −2 % | 0,94 | −0,07 | 0,23 | 0,13 |
| −1 % | 0,99 | −0,02 | 0,22 | 0,14 |
| ± 0 % | 1,02 | +0,03 | 0,24 | 0,11 |
| +1 % | 1,05 | +0,07 | 0,29 | 0,03 |
| +2 % | 1,09 | +0,11 | 0,33 | 0,01 |

(Beispiel: ein Bild des 14"-Spiegels, Durchmesser 355,6 mm, f = 1412 mm.)

Ein Fehler von 1 % im Radius verschiebt α um etwa 0,04 und die sphärische Aberration um etwa 0,05 λ RMS. Das ist so groß wie
der Befund selbst. Konsequenzen:

- Sphärische Aberration und α sind **nur so genau wie der Spiegelrand im Bild**. Auf besser als ±0,05 λ sollte man sie bei
  Klick- oder Auto-Pupille nicht deuten.
- Der **eingetragene Durchmesser** muss zum sichtbaren Rand passen: nicht den Rohlingsdurchmesser, sondern die wirklich
  genutzte Öffnung (Fase, Klammern, Blende).
- **Praktischer Tipp:** Eine Blende oder Maske mit genau bekanntem Durchmesser (z. B. Pappring, 2 bis 3 cm kleiner als der
  Spiegel) direkt vor den Spiegel setzen. Dann ist der Rand im Bild scharf, und der Durchmesser ist bekannt. Diesen Wert in
  das Tool eintragen.
- Die übrigen Terme (Astigmatismus, Koma, Dreiblatt) hängen deutlich weniger am Radius.

## 6. Was das Tool intern macht

1. **Pupille** als Ellipse bestimmen (Automatik oder per Klick).
2. **Streifenträger** (Kippung) im Frequenzspektrum finden.
3. **Hintergrund** abziehen, **demodulieren** und tiefpassfiltern (Fourier-Methode nach Takeda), in mehreren Durchgängen,
   wobei die schon geschätzte Wellenfront mit abgezogen wird.
4. **Streifenkontrast** pro Pixel bestimmen. Wo er zu schwach ist, wird maskiert. Diese Stellen gehen weder ins Entfalten noch
   in den Fit ein. Damit kann ein einzelner schwacher Bereich nicht mehr das ganze Ergebnis verderben.
5. **Entfalten** der Phase nur über die gültige Fläche.
6. **Zernike-Fit** (orthonormal, Noll-Indizes) über die gültigen Pixel, mit Ausreißerabwurf.
7. **Sollform** (Paraboloid gegen die Messkugel R = 2·f) abziehen.
8. **Serie:** je Bild auswerten, schlechte Bilder und Ausreißer aussortieren, Median bilden, Streuung als Wiederholbarkeit ausgeben.

## 7. Selbsttest

`test_core.py` und `test_batch.py` erzeugen künstliche Interferogramme mit bekannter Wellenfront (Rauschen, Übersteuerung,
schwache Zonen, Halo, Ellipse, Serien mit schlechten Bildern) und prüfen, ob die Koeffizienten wiedergefunden werden.
Im Test liegt der Fehler bei etwa 0,01 bis 0,03 Wellen, bei großen Zonenfehlern bei etwa 0,09. `pdi_synth.py` lässt sich auch
für eigene Tests nutzen. Die Tests prüfen nur die Auswertung, nicht den PDI-Aufbau.

## 7.1 Foucaulttest ↔ PDI vergleichen (Brücke zum Spiegelrechner)

Weil ein Spiegel (nahezu) rotationssymmetrisch ist, lässt sich aus einer Foucault-Zonenmessung ein **Dummy-PDI-Bild**
berechnen und mit einer echten Aufnahme vergleichen („Soll und Ist“).

- `pdi_foucault.py` benutzt `kern_foucault.py` direkt aus dem Ordner des Spiegelrechners (keine Kopie, eine Quelle der
  Wahrheit). Ordner: Argument `--foucault-dir`, Umgebungsvariable `PDI_FOUCAULT_DIR` oder ein Ordner `newton_spiegel*` neben dem Tool.
- Konvention (numerisch geprüft): Foucault-W ist die Abweichung von der Parabel in Wellen bei 550 nm (Kugel: W = −T·ρ⁴ mit
  T = D⁴/(512·f³·λ), bei 14″ f/3,97 etwa 17 λ). Das PDI-Tool mit Sollform liefert dieselbe Größe. Die Dummy-Wellenfront gegen die
  Messkugel ist T·ρ⁴ + W_Foucault, umgerechnet auf die PDI-Wellenlänge.
- **Im Tool:** Knopf **Foucault → Testbild …** öffnet eine Eingabemaske: Zonenzahl (4, 6 oder 8, optional Randzone geteilt) wählen,
  die Δl-Werte je Zone in mm eintragen (wie im Spiegelrechner, Zone 1 = 0), **Bild erzeugen** lädt es direkt ins Tool,
  oder als PNG speichern (z. B. für OpenFringe). **Parabel-Sollwerte einsetzen** füllt die Maske mit den Werten einer perfekten Parabel.
- Dasselbe per Kommandozeile: `python foucault_dummy.py --D 355.6 --f 1412 --df 0 0.93 1.87 2.80 3.73 4.66 --out dummy.png`
  (df = Foucault-Messwerte je Zone in mm, Zone 1 = 0; ohne `--df` ergibt sich eine Kugel). Das Bild lässt sich im Tool öffnen.
- Vergleich: `zone_profile()` liefert das Zonenprofil aus den PDI-Koeffizienten, `compare_profiles()` die Differenz zum
  Foucault-Profil (beide ohne Piston und mit bestem Defokus), `save_profile_csv()` schreibt es als CSV.
- Test `test_foucault.py` (Foucault → Dummy → Auswertung → Vergleich): Kugel α = 0,00, Profil-Differenz 3 nm RMS (Oberfläche);
  Parabel α = 0,99, 6 nm; Parabel mit Zonenfehlern α = 0,98, 5 nm.
- **Wichtig:** Eine Parabel hat gegen die Messkugel rund 17 λ Sphärik. Mit wenig Streifenträger (hier 18 Streifen) versagt die
  Auswertung (α 0,07). Mit etwa 40 oder mehr Kippstreifen über den Durchmesser läuft sie sauber. Das gilt auch für echte
  Aufnahmen: bei einer fast fertigen Parabel deutlich kippen, oder das Plättchen leicht defokussieren.

## 8. Bekannte Grenzen

- Ohne Streifenträger (geschlossene Ringe, kaum Kippung) versagt die Fourier-Methode. Dafür wäre Streifenverfolgung nötig
  (noch nicht umgesetzt).
- Systematische Fehler im Aufbau (zu großes Loch, Laserdiode ohne Linse mit Astigmatismus, Verzeichnung) kann kein
  Auswerteprogramm erkennen. Gegenprobe: Spiegel drehen, Foucaulttest.
- Genauigkeit am Spiegelrand ist geringer (enge Streifen, fallender Kontrast), siehe auch 5.3.
- Der Strehl-Wert ist eine grobe Abschätzung, siehe 5.1.

## 9. Dateien

| Datei | Inhalt |
|---|---|
| `pdi_gui.py` | Tkinter-Oberfläche |
| `pdi_core.py` | Auswertekern: Pupille, Fourier-Methode, Unwrapping, Zernike, Sollform |
| `pdi_batch.py` | Serienauswertung mit Qualitätsbewertung |
| `pdi_synth.py` | künstliche Interferogramme zum Testen |
| `pdi_report.py` | Diagnose-Grafik der Zwischenschritte |
| `pdi_compare.py` | Soll/Ist-Streifenbilder (Überhöhung) |
| `pdi_foucault.py`, `foucault_dummy.py` | Brücke zum Foucaulttest (Dummy-Bild, Profilvergleich) |
| `test_core.py`, `test_batch.py`, `test_foucault.py` | Selbsttests |
| `doku/` | Bilder zu dieser Dokumentation |
