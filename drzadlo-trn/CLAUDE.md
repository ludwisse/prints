# Projekt drzadlo-trn

## Komunikace
- Piš česky, tykej (uživatel je Petr Ludvík, nejsem programátor, ale rozumí IT).
- Odpovědi stručné a s čísly. Odhady označ jako odhady, měření jako měření. Neověřené věci řekni jako neověřené.

## Cíl
Držadlo do pěsti s trnem pro zatlačení magnetického držáku hlásičů (díl se závitem 9,29 mm a zoubky) do země. Stav, rozměry a historie: `README.md`.

## Pravidla
- Rozměry dílů **vždy ze změřeného skutečného protikusu**, ne z fotek (odhad z fotek byl o 1,6 mm vedle). Když rozměr chybí, požádej uživatele o změření (posuvným měřítkem).
- Na veškerý kód piš testy (pytest, složka `tests/`) a před commitem je spusť. Celá sada trvá asi 11 s.
- Commity jdou na větev `claude/epic-cannon-x2yxk9` v `ludwisse/prints`, pokud uživatel neřekne jinak. Pull requesty nezakládej, pokud o ně uživatel nepožádá.
- Soubory pro uživatele posílej také přímo do chatu (SendUserFile), ať je nemusí hledat v repozitáři ani na větvi.
- Neukládej cache (`__pycache__`, `.pytest_cache`), jsou v `.gitignore`.

## Konvence modelu
- Jednotky: mm. STL je v tiskové orientaci: dlaň na stole (z = 0), trn míří nahoru, osa tyče je x, trn leží uprostřed (`x_trn`).
- Trn: Ø `d_trn` = 14,0 mm, délka `v_trn` = 9,4 mm od dna kapsy; kapsa `d_kapsa` = 20,5 mm × `h_kapsa` = 1,5 mm; pata trnu kolmá (`r_paty` = 0).
- Dosedací plocha = dno kapsy (prstenec mezi trnem a stěnou kapsy).
- Duté varianty (`_dute`): stěna komor `stena_dutiny` = 4 mm (nejméně 2,4), plný blok `r_sloup` = 17 mm kolem trnu, plochá střecha komory nejvýš 2 × `pol_strecha` = 14 mm (přemostění), tisk se 100% výplní.

## Technické pasti
- STL nemůže slicer nutit k hustotě výplně. Pevnost se dá zapracovat jen geometrií (duté komory plus plný blok).
- Při změně parametrů vždy spusť `python drzadlo_tyc.py` (všech šest STL), `python vykres.py <S|M|L>` a testy; kóty ve výkresu se počítají z parametrů.
- V `drzadlo_tyc.py` jsou funkce mezi `trn` a `vyrob` (`plny_rez`, `rez_dutiny`, `komory`, ...). Při hromadných úpravách textu s hledáním podle názvu funkce je nesmaž.
- `prevesy` měří převisy jen na vnější ploše, vnitřní střechy komor se měří zvlášť (`prevesy_vnitrni`).
