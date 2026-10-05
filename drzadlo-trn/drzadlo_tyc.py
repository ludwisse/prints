"""Držadlo do pěsti: zploštělá příčná tyč s trnem pro zatlačení magnetického držáku hlásičů do země.

Trn (Ø 15,6 × 8 mm) zapadá do otvoru v dílu se závitem 9,29 mm a se zoubky nahoře; plochá dosedací plocha
kolem trnu (mělká kapsa Ø 21,2 mm, hloubka 1,5 mm, na dně prstenec mezi trnem a stěnou dílu) tlačí na zoubky.
Dlaň tlačí shora, prsty tyč obejmou zepředu a zespodu, trn vyčnívá dolů z rovného spodku mezi prostředníčkem
a prsteníčkem.

Rozměry trnu a kapsy jsou ODHAD z fotografií (měřítko: závit 9,29 mm): díl se závitem Ø 20,6 mm, krček horní
části Ø 15,8 mm, délka krčku 8,8 mm. Před tiskem ověř na skutečném dílu a případně uprav parametry.
Tyč je v půdoryse „obdélník s půlkruhovými konci"; příčný řez je zploštělý (rovná dlaňová plocha,
45° boky bez převisů, svislé stěny, zaoblená spodní hrana a rovný spodek).

STL je v TISKOVÉ orientaci: dlaň na stole (z = 0), trn míří nahoru, osa tyče je x.
Použití: python drzadlo_tyc.py [S|M|L]
"""
import os
import sys
from dataclasses import dataclass

import manifold3d as mf
import numpy as np
import trimesh

from drzadlo import _zaobli


@dataclass(frozen=True)
class Parametry:
    delka: float = 105.0            # mm, celková délka tyče (šířka ruky + přesah)
    hloubka: float = 44.0           # mm, rozměr tyče ve směru prstů (dlaň → konečky prstů)
    vyska: float = 30.0             # mm, od dlaně ke spodku (rovná plocha s trnem)
    pol_plocha: float = 12.0        # mm, poloviční šířka ploché dlaňové plochy (stůl při tisku)
    r_dlan_zaobleni: float = 4.0    # mm, zaoblení mezi plochou a 45° bokem
    r_bok: float = 8.0              # mm, zaoblení mezi bokem a svislou stěnou
    r_spodek: float = 7.0           # mm, zaoblení spodní hrany (u prstů)
    d_trn: float = 15.6             # mm, průměr trnu (krček dílu ≈ 15,8 mm, minus vůle pro tisk)
    v_trn: float = 8.0              # mm, délka trnu od dna kapsy (krček ≈ 8,8 mm, trn nesmí dosednout na dno otvoru)
    srazeni: float = 1.0            # mm, sražení hrany hlavy trnu (45°)
    r_paty: float = 1.2             # mm, zaoblení trnu u paty (nad dnem kapsy)
    d_kapsa: float = 21.2           # mm, kapsa kolem trnu: díl se závitem Ø ≈ 20,6 mm + vůle
    h_kapsa: float = 1.5            # mm, hloubka kapsy; dno kapsy je dosedací plocha pro zoubky
    sraz_kapsa: float = 0.5         # mm, sražení ústí kapsy (zasouvání)
    x_trn: float = 0.0              # mm, poloha trnu podél tyče (0 = střed)
    segmentu: int = 192


VELIKOSTI = {
    "S": Parametry(delka=95.0, hloubka=40.0, vyska=28.0, pol_plocha=10.0, r_bok=7.0, r_spodek=6.0),
    "M": Parametry(),
    "L": Parametry(delka=115.0, hloubka=48.0, vyska=32.0, pol_plocha=13.0, r_bok=9.0, r_spodek=8.0),
}


def pulprofil(p):
    """Polovina příčného řezu (r, z) od osy dlaně (0, 0) po osu spodku (0, vyska), r = vzdálenost od svislé osy."""
    r_max = p.hloubka / 2
    z_c = r_max - p.pol_plocha                              # 45° bok: dz = dr
    t_c = p.r_bok * np.tan(np.radians(45.0) / 2)
    t_d = p.r_spodek * np.tan(np.radians(90.0) / 2)
    if p.vyska - z_c < t_c + t_d:
        raise ValueError(f"svislá stěna je příliš krátká: {p.vyska - z_c:.1f} mm, potřeba aspoň {t_c + t_d:.1f} mm")
    if p.hloubka <= 2 * p.pol_plocha + 2:
        raise ValueError("plochá dlaňová plocha je téměř stejně široká jako tyč")
    if p.d_kapsa + 2 * p.sraz_kapsa + 3.0 > p.hloubka - 2 * p.r_spodek:
        raise ValueError("rovný spodek je příliš úzký pro kapsu kolem trnu (zvětši hloubku nebo zmenši r_spodek)")
    ostre = [(0, 0), (p.pol_plocha, 0), (r_max, z_c), (r_max, p.vyska), (0, p.vyska)]
    polomery = [0, p.r_dlan_zaobleni, p.r_bok, p.r_spodek, 0]
    return _zaobli(ostre, polomery)


def kapsa(p):
    """Mělká kapsa Ø d_kapsa (osa z) s sraženým ústím; dno je v z = vyska - h_kapsa. Jde přes spodek ven."""
    r = p.d_kapsa / 2
    z0, z1 = p.vyska - p.h_kapsa, p.vyska
    body = [(0, z0), (r, z0), (r, z1 - p.sraz_kapsa), (r + p.sraz_kapsa, z1), (r + p.sraz_kapsa, z1 + 1.0), (0, z1 + 1.0)]
    m = mf.Manifold.revolve(mf.CrossSection([np.array(body, float)]), p.segmentu)
    return m.translate((p.x_trn, 0.0, 0.0))


def trn(p):
    """Trn s plochou sraženou hlavou a zaoblenou patou; stojí na dně kapsy a zasahuje 1 mm pod dno (do těla)."""
    r_t = p.d_trn / 2
    z0 = p.vyska - p.h_kapsa                                   # dno kapsy
    z_hlava = z0 + p.v_trn
    f = np.radians(np.linspace(-90, -180, 24))
    pata = np.column_stack([r_t + p.r_paty + p.r_paty * np.cos(f), z0 + p.r_paty + p.r_paty * np.sin(f)])
    body = np.vstack([[(0, z0 - 1.0)], [(r_t + p.r_paty, z0 - 1.0)], pata,
                      [(r_t, z_hlava - p.srazeni), (r_t - p.srazeni, z_hlava), (0, z_hlava)]])
    m = mf.Manifold.revolve(mf.CrossSection([body]), p.segmentu)
    return m.translate((p.x_trn, 0.0, 0.0))


def vyrob(p=Parametry()):
    """Uzavřená síť držadla s trnem v tiskové orientaci."""
    pul = pulprofil(p)
    # puk: otočený půlprofil kolem svislé osy = konec tyče (půlkruh v půdoryse)
    puk = mf.Manifold.revolve(mf.CrossSection([np.vstack([pul, [(0.0, pul[0][1])]])]), p.segmentu)
    lc = p.delka - p.hloubka                                  # vzdálenost středů koncových půlkruhů
    if lc <= 0:
        raise ValueError("tyč je kratší než její hloubka")
    plny = np.vstack([pul, np.column_stack([-pul[::-1, 0], pul[::-1, 1]])[1:-1]])   # celý příčný řez, souměrný
    rez = mf.Manifold.extrude(mf.CrossSection([plny]), lc)       # osa podél local z, řez v (x, y) = (příčně, výška)
    # (u, v, w) -> (x = w, y = u, z = v): cyklická permutace = vlastní rotace
    rez = rez.transform(np.array([[0, 0, 1, -lc / 2], [1, 0, 0, 0], [0, 1, 0, 0]], float))
    telo = (rez + puk.translate((lc / 2, 0, 0)) + puk.translate((-lc / 2, 0, 0))) - kapsa(p)
    telo = telo + trn(p)
    g = telo.to_mesh()
    return trimesh.Trimesh(np.asarray(g.vert_properties)[:, :3], np.asarray(g.tri_verts), process=True)


def prevesy(m, mez_stupne=45.0, od_stolu_mm=1.5):
    """(podíl, plocha mm²) povrchu nad `od_stolu_mm`, který je převislý víc než `mez_stupne` od svislice."""
    uhel = np.degrees(np.arcsin(np.clip(-m.face_normals[:, 2], -1, 1)))
    zlobi = (uhel > mez_stupne + 0.5) & (m.triangles_center[:, 2] > od_stolu_mm)   # 0,5° tolerance: 45° boky jsou ještě v pořádku
    plocha = float(m.area_faces[zlobi].sum())
    return plocha / m.area, plocha


def main(argv):
    velikosti = argv[1:] or ["S", "M", "L"]
    for v in velikosti:
        p = VELIKOSTI[v.upper()]
        m = vyrob(p)
        out = os.path.join(os.path.dirname(os.path.abspath(__file__)), f"drzadlo_tyc_{v.upper()}.stl")
        m.export(out)
        print(f"{out}: watertight={m.is_watertight}, rozměry {m.extents.round(1)} mm, objem {m.volume / 1000:.1f} ml")


if __name__ == "__main__":
    main(sys.argv)
