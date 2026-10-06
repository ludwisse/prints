"""Výkres držadla s kótami (PNG + PDF): řez trnem, příčný řez tyčí a půdorys.

Použití: python vykres.py [S|M|L]   -> docs/vykres_<velikost>.png a .pdf
Kóty jsou počítané z parametrů modelu (drzadlo_tyc.Parametry), ne vypsané ručně, takže se s modelem nerozejdou.
"""
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy.spatial import ConvexHull  # noqa: E402

import drzadlo_tyc as t  # noqa: E402


def koty(p):
    """Všechny kótované rozměry (mm) odvozené z parametrů."""
    z_dno = p.vyska - p.h_kapsa
    return {
        "d_trn": p.d_trn,
        "delka_trnu": p.v_trn,
        "srazeni_hlavy": p.srazeni,
        "d_hlavy_plocha": p.d_trn - 2 * p.srazeni,
        "d_kapsa": p.d_kapsa,
        "h_kapsa": p.h_kapsa,
        "sraz_kapsa": p.sraz_kapsa,
        "presah_pod_spodek": p.v_trn - p.h_kapsa,
        "prstenec_sirka": (p.d_kapsa - p.d_trn) / 2,
        "delka": p.delka,
        "hloubka": p.hloubka,
        "vyska": p.vyska,
        "celkova_vyska": z_dno + p.v_trn,
        "plocha_dlane": 2 * (p.pol_plocha),
        "rovny_spodek": p.hloubka - 2 * p.r_spodek,
        "x_trn": p.x_trn,
        "od_konce_k_trnu": p.delka / 2 + p.x_trn,
        "r_bok": p.r_bok,
        "r_spodek": p.r_spodek,
        "r_dlan_zaobleni": p.r_dlan_zaobleni,
    }


def _kota(ax, a, b, text, odsazeni=0.0, fs=9, barva="tab:blue", pred=0.0):
    """Kóta mezi body a, b; `odsazeni` = posun kótové čáry kolmo ke spojnici, text uprostřed."""
    a, b = np.array(a, float), np.array(b, float)
    d = b - a
    n = np.array([-d[1], d[0]]) / np.linalg.norm(d)
    a2, b2 = a + n * odsazeni, b + n * odsazeni
    if odsazeni:
        for p0, p1 in ((a, a2), (b, b2)):
            ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color=barva, lw=0.6)
    ax.annotate("", xy=a2, xytext=b2, arrowprops=dict(arrowstyle="<->", color=barva, lw=0.9, shrinkA=0, shrinkB=0))
    st = (a2 + b2) / 2
    ang = np.degrees(np.arctan2(d[1], d[0]))
    if ang > 90 or ang < -90:
        ang += 180
    ax.text(st[0] + n[0] * (pred or 1.4), st[1] + n[1] * (pred or 1.4), text, color=barva, fontsize=fs,
            ha="center", va="center", rotation=ang if abs(ang) in (90,) else 0,
            bbox=dict(fc="white", ec="none", pad=0.6))


def kresli(p, soubor_png, titulek="M"):
    k = koty(p)
    m = t.vyrob(p)
    z_hlava = p.vyska - p.h_kapsa + p.v_trn
    fig = plt.figure(figsize=(17, 11))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 1.0], width_ratios=[1.0, 1.0])

    # --- A: řez trnem (pracovní orientace: trn dole) -------------------------------------
    a = fig.add_subplot(gs[0, 0])
    sek = m.section(plane_origin=[p.x_trn, 0, 0], plane_normal=[1, 0, 0])
    for e in sek.entities:
        v = sek.vertices[e.points]
        a.plot(v[:, 1], z_hlava - v[:, 2], "k-", lw=1.5)
    r_t, r_k = p.d_trn / 2, p.d_kapsa / 2
    zdno = p.v_trn                                     # dno kapsy (od hlavy trnu)
    zspod = p.v_trn - p.h_kapsa                        # rovný spodek tyče (od hlavy trnu)
    zvrch = zspod + p.vyska                            # dlaň
    a.axhline(0, color="gray", lw=0.4, ls=":")
    _kota(a, (-r_t, 0), (r_t, 0), f"Ø {k['d_trn']:g}", odsazeni=-3.0, pred=-1.2)
    _kota(a, (-r_t, 0), (-r_t, zdno), f"{k['delka_trnu']:g}", odsazeni=-3.0)
    _kota(a, (r_k + 4, zspod), (r_k + 4, zdno), f"{k['h_kapsa']:g}", odsazeni=0.0)
    _kota(a, (-r_k, zdno), (r_k, zdno), f"Ø {k['d_kapsa']:g}", odsazeni=3.5)
    _kota(a, (r_t + 4, 0), (r_t + 4, zspod), f"{k['presah_pod_spodek']:g}", odsazeni=0.0)
    a.annotate(f"sražení {k['srazeni_hlavy']:g} × 45°", xy=(r_t - 0.4, 0.5), xytext=(r_t + 8, -4.5),
               arrowprops=dict(arrowstyle="->"), fontsize=9)
    a.annotate(f"plochá hlava Ø {k['d_hlavy_plocha']:g}", xy=(-r_t + 1.0, 0), xytext=(-r_t - 4, -5.5), ha="right", fontsize=9)
    a.annotate(f"sražení ústí kapsy {k['sraz_kapsa']:g} × 45°", xy=(r_k + 0.3, zspod + 0.2), xytext=(r_k + 9, zspod + 9),
               arrowprops=dict(arrowstyle="->"), fontsize=9)
    a.annotate("dosedací prstenec (dno kapsy):\n" f"šířka {k['prstenec_sirka']:g}, kolmá pata trnu bez rádiusu",
               xy=(-(r_t + r_k) / 2, zdno), xytext=(-27, zdno + 6), ha="right", fontsize=9, color="crimson",
               arrowprops=dict(arrowstyle="->", color="crimson"))
    a.text(r_t + 8, zspod / 2 - 0.5, f"trn přesahuje spodek\ntyče o {k['presah_pod_spodek']:g}", fontsize=9, va="center")
    a.annotate("dlaň tlačí shora ↓", xy=(0, zvrch), xytext=(0, zvrch + 3), ha="center", fontsize=10)
    a.set_aspect("equal"); a.axis("off")
    a.set_xlim(-46, 46); a.set_ylim(-9, zvrch + 7)
    a.set_title("Řez A–A osou trnu (pracovní orientace: trn dole)", fontsize=12)

    # --- B: příčný řez tyčí --------------------------------------------------------------
    b = fig.add_subplot(gs[0, 1])
    pul = t.pulprofil(p)
    r, z = pul[:, 0], p.vyska - pul[:, 1]           # dlaň nahoře
    xs = np.r_[r, -r[::-1]]; zs = np.r_[z, z[::-1]]
    b.plot(xs, zs, "k-", lw=1.5)
    h2, pp = p.hloubka / 2, p.pol_plocha
    _kota(b, (-pp, p.vyska), (pp, p.vyska), f"plochá dlaň {k['plocha_dlane']:g}", odsazeni=5.0)
    _kota(b, (h2, 0), (h2, p.vyska), f"{k['vyska']:g}", odsazeni=-6.0)
    sp = h2 - p.r_spodek
    _kota(b, (-sp, 0), (sp, 0), f"rovný spodek {k['rovny_spodek']:g}", odsazeni=-5.0)
    _kota(b, (-h2, 0), (h2, 0), f"{k['hloubka']:g}", odsazeni=-13.0)
    b.annotate(f"R {p.r_dlan_zaobleni:g}", xy=(pp + 1, p.vyska - 0.3), xytext=(pp + 8, p.vyska + 4), arrowprops=dict(arrowstyle="->"), fontsize=9)
    b.annotate("45° bok", xy=((pp + h2) / 2 + 1, p.vyska - (h2 - pp) / 2 + 0.5), xytext=(h2 + 7, p.vyska - 4),
               arrowprops=dict(arrowstyle="->"), fontsize=9)
    b.annotate(f"R {p.r_bok:g}", xy=(h2 - 0.4, p.vyska - (h2 - pp) - 2.0), xytext=(h2 + 8, p.vyska - 16), arrowprops=dict(arrowstyle="->"), fontsize=9)
    b.annotate(f"R {p.r_spodek:g}", xy=(h2 - 1.6, 1.8), xytext=(h2 + 8, 5), arrowprops=dict(arrowstyle="->"), fontsize=9)
    b.text(0, p.vyska / 2, "plný řez\n(duté verze viz poznámka)", ha="center", fontsize=9, color="gray")
    b.set_aspect("equal"); b.axis("off"); b.set_xlim(-h2 - 8, h2 + 22); b.set_ylim(-19, p.vyska + 11)
    b.set_title("Příčný řez tyčí mimo trn (dlaň nahoře)", fontsize=12)

    # --- C: půdorys -----------------------------------------------------------------------
    c = fig.add_subplot(gs[1, :])
    v = m.vertices
    hull = ConvexHull(v[:, :2])
    c.fill(v[hull.vertices, 0], v[hull.vertices, 1], fc="#f2f2f2", ec="k", lw=1.5)
    tt = np.linspace(0, 2 * np.pi, 200)
    c.plot(p.x_trn + r_k * np.cos(tt), r_k * np.sin(tt), "--", color="crimson", lw=1.0)
    c.plot(p.x_trn + r_t * np.cos(tt), r_t * np.sin(tt), "-", color="k", lw=1.5)
    c.plot([p.x_trn, p.x_trn], [-p.hloubka / 2 - 3, p.hloubka / 2 + 3], "-.", color="gray", lw=0.6)
    _kota(c, (-p.delka / 2, -p.hloubka / 2), (p.delka / 2, -p.hloubka / 2), f"{k['delka']:g}", odsazeni=-5.0)
    _kota(c, (-p.delka / 2, -p.hloubka / 2), (-p.delka / 2, p.hloubka / 2), f"{k['hloubka']:g}", odsazeni=-5.0)
    _kota(c, (-p.delka / 2, p.hloubka / 2 + 4), (p.x_trn, p.hloubka / 2 + 4), f"{k['od_konce_k_trnu']:g} (trn uprostřed)", odsazeni=0.0)
    c.text(p.x_trn + 12, 12, f"Ø {k['d_trn']:g} trn\nØ {k['d_kapsa']:g} kapsa", fontsize=9, color="crimson")
    c.text(0, p.hloubka / 2 + 11, "Pohled shora na tyč (trn leží mezi prostředníčkem a prsteníčkem)", ha="center", fontsize=10)
    c.text(p.x_trn + 1.0, p.hloubka / 2 + 5.5, "A–A", ha="left", fontsize=10, color="gray")
    c.set_aspect("equal"); c.axis("off"); c.set_xlim(-p.delka / 2 - 15, p.delka / 2 + 15); c.set_ylim(-p.hloubka / 2 - 12, p.hloubka / 2 + 16)
    c.set_title("Půdorys", fontsize=12)

    fig.suptitle(f"Držadlo – příčná tyč {titulek}, trn Ø {k['d_trn']:g} × {k['delka_trnu']:g} mm pro magnetický držák hlásičů   (jednotky: mm)",
                 fontsize=14, y=0.995)
    fig.text(0.01, 0.005,
             f"Celková výška s trnem {k['celkova_vyska']:g} mm. STL je v tiskové orientaci (dlaň na stole, trn nahoru). "
             "Duté varianty (_dute): dvě uzavřené komory, stěny 4 mm, plný blok ±17 mm kolem trnu, tisk se 100% výplní.",
             fontsize=8.5, color="dimgray")
    plt.tight_layout(rect=(0, 0.02, 1, 0.98))
    plt.savefig(soubor_png, dpi=110)
    plt.savefig(os.path.splitext(soubor_png)[0] + ".pdf")
    plt.close(fig)
    return k


if __name__ == "__main__":
    v = (sys.argv[1] if len(sys.argv) > 1 else "M").upper()
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "docs", f"vykres_{v}.png")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    kresli(t.VELIKOSTI_PLNE[v], out, v)
    print(out)
