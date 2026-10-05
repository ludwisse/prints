"""Držadlo do dlaně s trnem Ø 9 × 12 mm (na zatlačení kolíku do země).

Tvar je rotační těleso, STL je v TISKOVÉ orientaci: dlaňová plocha je na stole (z = 0), trn míří nahoru.
V použití je držadlo otočené: trn dole, dlaň tlačí shora.

Profil (r, z) od osy: plochá dlaňová plocha, zaoblení a 45° bok (tisknutelné bez podpěr), svislý okraj,
kuželová spodní strana pro prsty (30°), plochá opěrná plocha Ø 34 mm pro hlavu kolíku, trn
s přechodovým zaoblením u paty a plochou hlavou se sraženou hranou.
Použití: python drzadlo.py [průměr_dlaně_mm]
"""
import os
import sys
from dataclasses import dataclass

import manifold3d as mf
import numpy as np
import trimesh


@dataclass(frozen=True)
class Parametry:
    d_trn: float = 9.0             # mm, průměr trnu
    v_trn: float = 12.0            # mm, výška trnu nad opěrnou plochou
    srazeni: float = 1.0           # mm, sražení hrany hlavy trnu (45°)
    r_paty: float = 2.0            # mm, zaoblení trnu u paty (proti vrubu)
    d_dlane: float = 80.0          # mm, největší průměr držadla
    r_plocha: float = 28.0         # mm, poloměr ploché dlaňové plochy (stůl při tisku)
    r_zaobleni_dlane: float = 4.0  # mm, zaoblení mezi plochou a 45° bokem
    r_okraj: float = 8.0           # mm, zaoblení mezi bokem a svislým okrajem
    z_okraj_konec: float = 21.0    # mm, výška horního konce svislého okraje
    r_okraj_horni: float = 6.0     # mm, zaoblení mezi svislým okrajem a spodním kuželem
    uhel_kuzele: float = 30.0      # °, sklon spodní strany (pro prsty) vůči vodorovné rovině
    d_opora: float = 34.0          # mm, průměr ploché opěrné plochy kolem trnu
    r_opora: float = 3.0           # mm, zaoblení hrany opěrné plochy
    segmentu: int = 256            # počet segmentů otáčení


def _zaobli(body, polomery, n=24):
    """Polyline `body` (N×2) se zaoblenými rohy: oblouk o poloměru `polomery[i]` v i-tém bodě (0 = ostrý roh).

    Krajní body zůstávají beze změny. Oblouk se vzorkuje n body a je tečný k oběma sousedním úsečkám."""
    body = np.asarray(body, float)
    vysl = [body[0]]
    for i in range(1, len(body) - 1):
        r = polomery[i]
        p0, p1, p2 = body[i - 1], body[i], body[i + 1]
        if r <= 0:
            vysl.append(p1)
            continue
        u = (p0 - p1) / np.linalg.norm(p0 - p1)
        v = (p2 - p1) / np.linalg.norm(p2 - p1)
        uhel = np.arccos(np.clip(u @ v, -1, 1))          # vnitřní úhel v rohu
        stred = p1 + (u + v) / np.linalg.norm(u + v) * (r / np.sin(uhel / 2))
        a = p1 + u * (r / np.tan(uhel / 2))              # tečné body
        b = p1 + v * (r / np.tan(uhel / 2))
        fa = np.arctan2(a[1] - stred[1], a[0] - stred[0])
        fb = np.arctan2(b[1] - stred[1], b[0] - stred[0])
        d = (fb - fa + np.pi) % (2 * np.pi) - np.pi      # kratší oblouk
        vysl.extend(stred + r * np.array([np.cos(fa + d * k), np.sin(fa + d * k)]) for k in np.linspace(0, 1, n))
    vysl.append(body[-1])
    return np.array(vysl)


def profil(p=Parametry()):
    """Uzavřený profil (r, z) pro otočení kolem osy z: vrací (body, z_opora, z_hlava)."""
    r_max = p.d_dlane / 2
    z_c = r_max - p.r_plocha                                  # 45° bok: dz = dr
    z_d = p.z_okraj_konec
    r_n = p.d_opora / 2
    z_e = z_d + (r_max - r_n) * np.tan(np.radians(p.uhel_kuzele))
    z_hlava = z_e + p.v_trn
    r_t = p.d_trn / 2
    # oblouky na svislém okraji se nesmí překrývat (tečné délky: R * tan(úhel_obratu / 2))
    t_c = p.r_okraj * np.tan(np.radians(45.0) / 2)
    t_d = p.r_okraj_horni * np.tan(np.radians(90.0 - p.uhel_kuzele) / 2)
    if z_d - z_c < t_c + t_d:
        raise ValueError(f"svislý okraj je příliš krátký: {z_d - z_c:.1f} mm, potřeba aspoň {t_c + t_d:.1f} mm "
                         "(zvětši z_okraj_konec nebo zmenši poloměry zaoblení)")
    if p.srazeni >= r_t:
        raise ValueError("sražení hrany trnu je větší než jeho poloměr")
    # ostrý obrys: osa, plocha, 45° bok, svislý okraj, kužel pro prsty, opěrná plocha, pata trnu
    ostre = [(0, 0), (p.r_plocha, 0), (r_max, z_c), (r_max, z_d), (r_n, z_e), (r_t, z_e), (r_t, z_hlava - p.srazeni)]
    polomery = [0, p.r_zaobleni_dlane, p.r_okraj, p.r_okraj_horni, p.r_opora, p.r_paty, 0]
    zaobl = _zaobli(ostre, polomery)
    # hlava trnu: sražená hrana a plochá čelní plocha, pak zpět na osu
    body = np.vstack([zaobl, [(r_t - p.srazeni, z_hlava)], [(0.0, z_hlava)]])
    return body, z_e, z_hlava


def vyrob(p=Parametry()):
    """Uzavřená síť držadla s trnem v tiskové orientaci (z = 0 je dlaňová plocha)."""
    body, _, _ = profil(p)
    # na ose nesmí být zdvojené body, jinak revolve vytvoří degenerované trojúhelníky
    prurez = mf.CrossSection([np.asarray(body, float)])
    m = mf.Manifold.revolve(prurez, p.segmentu)
    g = m.to_mesh()
    return trimesh.Trimesh(np.asarray(g.vert_properties)[:, :3], np.asarray(g.tri_verts), process=True)


def preruseni_prevesu(m, mez_stupne=45.0, od_stolu_mm=1.5):
    """Podíl plochy povrchu, který je nad `od_stolu_mm` převislý víc než `mez_stupne` od svislice."""
    n = m.face_normals
    dolu = -n[:, 2]                                  # kosinus úhlu normály od směru dolů
    uhel_od_svisle = np.degrees(np.arcsin(np.clip(dolu, -1, 1)))   # >0 = plocha míří dolů
    vys = m.triangles_center[:, 2]
    zlobi = (uhel_od_svisle > mez_stupne) & (vys > od_stolu_mm)
    return float(m.area_faces[zlobi].sum() / m.area), float(m.area_faces[zlobi].sum())


def main(argv):
    d = float(argv[1]) if len(argv) > 1 else 80.0
    k = d / 80.0
    p = Parametry(d_dlane=d, r_plocha=28.0 * k, z_okraj_konec=21.0 * k)
    m = vyrob(p)
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), f"drzadlo_trn_d{d:g}.stl")
    m.export(out)
    print(f"{out}: watertight={m.is_watertight}, rozměry {m.extents.round(2)} mm, objem {m.volume / 1000:.1f} ml")


if __name__ == "__main__":
    main(sys.argv)
