# Držadlo s trnem pro zatlačení magnetického držáku hlásičů

Držadlo do pěsti (zploštělá příčná tyč), kterým se dlaní zatlačuje do země díl magnetického držáku rybářských hlásičů: ten se závitem 9,29 mm dole a se zoubky nahoře. Trn zapadne do otvoru v tom dílu a plochá dosedací plocha kolem trnu tlačí na zoubky. Dlaň tlačí shora, prsty tyč obejmou zepředu a zespodu a trn vyčnívá dolů mezi prostředníčkem a prsteníčkem.

**Stav: projekt dočasně uzavřený (6. 10. 2026).** STL jsou hotové a rozměry trnu jsou podle skutečného protikusu, ale nic zatím nebylo vytištěno a vyzkoušeno. Pokračovat se bude po zkušebním tisku (viz „Co zbývá").

## Rozměry (velikost M)
| | |
|---|---|
| Trn | Ø **14,0 mm**, délka **9,4 mm** od dna kapsy, plochá hlava Ø 12 se sražením 1 × 45°, **kolmá pata bez rádiusu** |
| Kapsa kolem trnu | Ø **20,5 mm**, hloubka 1,5 mm, sražení ústí 0,5 × 45°; dno kapsy je dosedací plocha, prstenec široký 3,25 mm tlačí na zoubky |
| Tyč | 105 × 44 mm v půdoryse (půlkruhové konce), výška 30 mm |
| Příčný řez | rovná dlaňová plocha 24 mm, 45° boky, svislé stěny, rovný spodek 30 mm, zaoblení R 4 / R 8 / R 7 |
| Celková výška s trnem | 37,9 mm; trn přesahuje spodek tyče o 7,9 mm |
| Poloha trnu | uprostřed tyče (52,5 mm od konců) |

Velikosti: **S** (95 × 40 × 28), **M** (105 × 44 × 30), **L** (115 × 48 × 32). Trn, kapsa a poloměry detailů jsou u všech stejné, mění se jen tyč.

## Soubory
| Soubor | Obsah |
|---|---|
| `drzadlo_tyc_{S,M,L}.stl` | plná držadla (M: objem 112 ml) |
| `drzadlo_tyc_{S,M,L}_dute.stl` | duté verze pro tisk se **100% výplní** (M: 81,5 ml, asi 101 g PLA) |
| `zkusebni_trny.stl` | destička s pěti trny Ø 13,6 / 13,8 / 14,0 / 14,2 / 14,4 × 9,4 mm k ověření lícování |
| `docs/vykres_{S,M,L}.png` a `.pdf` | výkres s kótami: řez trnem A–A, příčný řez tyčí, půdorys |
| `docs/tyc_dute.png` | řez dutou verzí |
| `drzadlo_tyc.py` | generátor (parametry v `Parametry`) |
| `vykres.py` | výkres s kótami z parametrů modelu |
| `zkusebni_trny.py` | generátor zkušební destičky |
| `drzadlo.py` | pomocná funkce `_zaobli` (zaoblení polyline); původní kulatá verze držadla je zrušená |
| `tests/` | testy (pytest), 61 testů, asi 11 s |

STL jsou v **tiskové orientaci**: dlaňová plocha na stole (z = 0), trn míří nahoru.

## Duté verze
STL nemůže slicer nutit k nějaké hustotě výplně, proto je pevnost zapracovaná do geometrie: dvě uzavřené komory se stěnami 4 mm a plný blok ±17 mm kolem trnu přes celou tloušťku (dráha síly od dlaně k dosedací ploše). Střecha komor je 45° s plochou částí 14 mm (přemostění). Tisk: výplň 100 % globálně, 4 stěny. Žádný modifikátor není potřeba.

## Doporučení k tisku (neověřené)
PETG (houževnatější) nebo PLA, orientace beze změny (dlaň na stole). U plných verzí zvýšit výplň a stěny v oblasti trnu. Zoubky dílu jsou ostré a při plné síle mohou zatlačit do plastu dosedací plochy.

## Použití
```
pip install numpy scipy trimesh manifold3d matplotlib pytest
python drzadlo_tyc.py              # všech šest STL (nebo: python drzadlo_tyc.py M_dute)
python vykres.py M                 # docs/vykres_M.png a .pdf
python zkusebni_trny.py            # zkusebni_trny.stl
python -m pytest tests
```
Skripty se spouštějí ze složky `drzadlo-trn`.

## Historie rozhodnutí
1. **Kulaté držadlo** (puk Ø 80 mm s trnem Ø 9 × 12): nahrazeno. Uživatel navrhl tvar zploštělého válce naležato, který se dá obejmout prsty a trn vyčnívá mezi prsteníčkem a prostředníčkem.
2. **Trn podle fotek:** z fotek protikusu (měřítko závit 9,29 mm) jsem odhadl Ø 15,6 × 8 mm a kapsu Ø 21,2 mm. **Odhad byl špatně** (o 1,6 mm) a trn se nevešel.
3. **Trn podle skutečného protikusu:** Ø 14,0 × 9,4 mm, kapsa Ø 20,5 mm.
4. **Pata trnu** bez zaoblení (původně R 1,2), čistě kolmá.
5. **Duté varianty:** uživatel neuměl udělat modifikátor výplně ve slicerovém programu, proto je pevnost řešená geometrií.

## Poznatky
- Rozměry z fotek jsou nespolehlivé (chyba asi 11 %). Rozměry dílů zadávej z měření skutečného protikusu.
- Plná pevnost se nedá vynutit v STL, jen geometrií (komory plus plný blok) a globálním nastavením výplně.
- Orientace a tvar jsou volené pro tisk bez podpěr: vnější tvar nemá převisy nad 45°, střechy komor mají plochou část nejvýš 14 mm.

## Co zbývá
1. Vytisknout `zkusebni_trny.stl` a zjistit, který průměr (13,6 až 14,4) do otvoru zapadne s lehkým odporem; případně upravit `d_trn`.
2. Změřit skutečný vnější průměr dílu se zoubky a ověřit, že se do kapsy Ø 20,5 vejde; případně upravit `d_kapsa`.
3. Vytisknout držadlo (doporučeno `M_dute` se 100% výplní), vyzkoušet v ruce (poloha trnu mezi prsty, délka tyče) a při zatlačení sledovat dosedací plochu a trn.
4. Podle výsledků případně zaoblit patu trnu (`r_paty`), upravit polohu trnu (`x_trn`) nebo velikost (S/M/L).

Pozn.: vše je na větvi `claude/epic-cannon-x2yxk9` (ne na `main`).
