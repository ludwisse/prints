import os
import sys

import numpy as np
import pytest
import trimesh

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import drzadlo_tyc as t  # noqa: E402


@pytest.fixture(scope="module")
def tvar():
    p = t.Parametry()
    return p, t.vyrob(p)


# --- rozměry zadání ------------------------------------------------------------------------

def test_trn_prumer_15_6_a_delka_8_od_dna_kapsy(tvar):
    p, m = tvar
    z_dno = p.vyska - p.h_kapsa
    assert m.bounds[1][2] == pytest.approx(z_dno + p.v_trn)
    sek = m.section(plane_origin=[0, 0, z_dno + 4.0], plane_normal=[0, 0, 1])
    v = sek.vertices
    assert np.hypot(v[:, 0], v[:, 1]).max() == pytest.approx(p.d_trn / 2, abs=0.001)
    assert v[:, 1].max() - v[:, 1].min() == pytest.approx(p.d_trn, abs=0.01)


def test_hlava_trnu_je_plochá_se_srazenou_hranou(tvar):
    p, m = tvar
    z_hlava = p.vyska - p.h_kapsa + p.v_trn
    horni = m.vertices[np.isclose(m.vertices[:, 2], z_hlava, atol=1e-6)]
    assert len(horni) > 10
    assert np.hypot(horni[:, 0], horni[:, 1]).max() == pytest.approx(p.d_trn / 2 - p.srazeni, abs=0.001)
    pod = m.vertices[np.isclose(m.vertices[:, 2], z_hlava - p.srazeni, atol=1e-6)]
    assert np.hypot(pod[:, 0], pod[:, 1]).max() == pytest.approx(p.d_trn / 2, abs=0.001)


def test_pata_trnu_je_zaoblena(tvar):
    p, m = tvar
    sek = m.section(plane_origin=[0, 0, p.vyska - p.h_kapsa + 0.1], plane_normal=[0, 0, 1])
    assert np.hypot(sek.vertices[:, 0], sek.vertices[:, 1]).max() > p.d_trn / 2 + 0.8     # fillet R1,2 rozšiřuje trn u dna


def test_kapsa_kolem_trnu_ma_prumer_a_hloubku(tvar):
    p, m = tvar
    z_dno = p.vyska - p.h_kapsa
    # těsně nad dnem je kapsa: průřez vnějším obrysem kapsy má poloměr d_kapsa / 2
    sek = m.section(plane_origin=[0, 0, z_dno + 0.8], plane_normal=[0, 0, 1])
    r = np.hypot(sek.vertices[:, 0], sek.vertices[:, 1])
    kolem_kapsy = r[(r > p.d_trn / 2 + 2.0) & (r < 15.0)]                  # jen obrys kapsy, ne vnější obrys tyče
    assert len(kolem_kapsy) > 10
    assert kolem_kapsy.max() < p.d_kapsa / 2 + p.sraz_kapsa + 0.01
    assert kolem_kapsy.min() == pytest.approx(p.d_kapsa / 2, abs=0.01)
    # dno kapsy je vodorovné: v z = dno jsou body mezi trnem a stěnou kapsy
    dno = m.vertices[np.isclose(m.vertices[:, 2], z_dno, atol=1e-6)]
    rd = np.hypot(dno[:, 0], dno[:, 1])
    assert rd.max() == pytest.approx(p.d_kapsa / 2, abs=0.01)


def test_pro_dilu_se_zavitem_je_misto_v_kapse(tvar):
    p, _ = tvar
    assert p.d_kapsa > 20.6 + 0.4                       # díl Ø ≈ 20,6 mm + vůle
    assert p.d_trn < 15.8                               # krček horní části ≈ 15,8 mm: trn je o vůli menší
    assert p.v_trn < 8.8                                # trn nesmí dosednout na dno otvoru


def test_bearing_prstenec_mezi_trnem_a_stenou_dilu(tvar):
    p, _ = tvar
    sirka = (20.6 - p.d_trn) / 2
    assert sirka > 2.0                                  # zoubky dosedají na prstenec široký asi 2,5 mm


def test_rozmery_tyce(tvar):
    p, m = tvar
    assert m.extents[0] == pytest.approx(p.delka, abs=0.05)
    assert m.extents[1] == pytest.approx(p.hloubka, abs=0.05)
    assert m.bounds[0][2] == pytest.approx(0.0, abs=1e-9)                       # dlaň na stole
    assert m.bounds[1][2] == pytest.approx(p.vyska - p.h_kapsa + p.v_trn)


def test_trn_je_uprostred_tyce_a_posun_funguje(tvar):
    p, m = tvar
    hlava = m.vertices[np.isclose(m.vertices[:, 2], m.bounds[1][2], atol=1e-6)]
    assert hlava[:, 0].mean() == pytest.approx(0.0, abs=0.01)
    m2 = t.vyrob(t.Parametry(x_trn=10.0))
    hlava2 = m2.vertices[np.isclose(m2.vertices[:, 2], m2.bounds[1][2], atol=1e-6)]
    assert hlava2[:, 0].mean() == pytest.approx(10.0, abs=0.01)


def test_spodek_je_rovny_kolem_trnu(tvar):
    p, m = tvar
    # v z = vyska je plochý pruh široký 2 * (hloubka/2 - r_spodek); kapsa z něj vyřízla kruh
    na_spodku = m.vertices[np.isclose(m.vertices[:, 2], p.vyska, atol=1e-6)]
    assert len(na_spodku) > 10
    assert np.abs(na_spodku[:, 1]).max() == pytest.approx(p.hloubka / 2 - p.r_spodek, abs=0.01)
    assert (p.hloubka - 2 * p.r_spodek) >= p.d_kapsa + 2 * p.sraz_kapsa + 3.0           # kolem kapsy zbývá rovná plocha


# --- kvalita a tisk ----------------------------------------------------------------------

def test_sit_je_uzavrena_a_kladna(tvar):
    _, m = tvar
    assert m.is_watertight and m.is_winding_consistent and m.volume > 0


def test_objem_odpovida_analyticky(tvar):
    p, m = tvar
    pul = t.pulprofil(p)
    r, z = pul[:, 0], pul[:, 1]
    plocha = 0.5 * abs(np.dot(r, np.roll(z, -1)) - np.dot(z, np.roll(r, -1)))
    cr = (1 / (6 * plocha)) * np.sum((r + np.roll(r, -1)) * (r * np.roll(z, -1) - np.roll(r, -1) * z))
    puk = 2 * np.pi * abs(cr) * plocha
    lc = p.delka - p.hloubka
    ocekavano = 2 * plocha * lc + puk
    trn = np.pi * 4.5 ** 2 * p.v_trn
    assert m.volume == pytest.approx(ocekavano + trn, rel=0.02)


def test_zadne_previsy_nad_45_stupnu(tvar):
    _, m = tvar
    assert t.prevesy(m)[1] == 0.0


def test_dlan_plocha_na_stole_je_dost_siroka(tvar):
    p, m = tvar
    dno = m.vertices[np.isclose(m.vertices[:, 2], 0.0, atol=1e-9)]
    assert np.abs(dno[:, 1]).max() == pytest.approx(p.pol_plocha - p.r_dlan_zaobleni * np.tan(np.radians(22.5)), abs=0.05)


def test_pulprofil_zacina_a_konci_na_ose():
    pul = t.pulprofil(t.Parametry())
    assert pul[0] == pytest.approx([0, 0]) and pul[-1] == pytest.approx([0, 30.0])
    assert (pul[:, 0] >= -1e-9).all() and pul[:, 0].max() == pytest.approx(22.0, abs=1e-6)


def test_prilis_kratka_stena_je_odmitnuta():
    with pytest.raises(ValueError):
        t.pulprofil(t.Parametry(vyska=14.0))


def test_prilis_uzky_spodek_pro_kapsu_je_odmitnut():
    with pytest.raises(ValueError):
        t.pulprofil(t.Parametry(hloubka=34.0, r_spodek=8.0))


def test_tyc_kratsi_nez_hloubka_je_odmitnuta():
    with pytest.raises(ValueError):
        t.vyrob(t.Parametry(delka=30.0))


@pytest.mark.parametrize("velikost", ["S", "L"])
def test_dalsi_velikosti(velikost):
    p = t.VELIKOSTI[velikost]
    m = t.vyrob(p)
    assert m.is_watertight
    assert m.extents[0] == pytest.approx(p.delka, abs=0.05)
    assert t.prevesy(m)[1] == 0.0
    assert m.bounds[1][2] == pytest.approx(p.vyska - p.h_kapsa + p.v_trn)


# --- duté varianty (pro tisk se 100% výplní) -------------------------------------------------

@pytest.fixture(scope="module")
def dute():
    p = t.VELIKOSTI["M_dute"]
    return p, t.vyrob(p)


def test_dute_ma_dve_uzavrene_komory(dute):
    p, m = dute
    casti = m.split(only_watertight=False)
    kladne = [c for c in casti if c.volume > 0]
    zaporne = [c for c in casti if c.volume < 0]
    assert m.is_watertight and len(kladne) == 1 and len(zaporne) == 2
    assert all(-c.volume / 1000 > 10 for c in zaporne)               # každá komora > 10 ml


def test_dute_je_lehci_nez_plne(dute):
    _, m = dute
    plne = t.vyrob(t.VELIKOSTI["M"])
    assert m.volume < 0.8 * plne.volume


def test_dute_stena_je_aspon_3_9_mm(dute):
    p, m = dute
    vnejsi = t.vnejsi_sit(m)
    komory_ = [c for c in m.split(only_watertight=False) if c.volume < 0]
    for k in komory_:
        hloubka = trimesh.proximity.signed_distance(vnejsi, k.vertices)     # kladná = uvnitř tyče
        assert hloubka.min() > p.stena_dutiny - 0.1


def test_dute_plny_blok_kolem_trnu_zustava_plny(dute):
    p, m = dute
    body = np.array([[0, 0, 15.0], [0, 8, 10.0], [0, -8, 20.0], [p.r_sloup - 1.0, 0, 15.0], [-(p.r_sloup - 1.0), 0, 15.0]])
    assert m.contains(body).all()
    # a komory jsou mimo: uprostřed komory vlevo není materiál
    assert not m.contains([[-(p.r_sloup + 8.0), 0, 15.0]])[0]


def test_dute_stropy_komor_jsou_primostitelne(dute):
    _, m = dute
    sirka, plocha_sikma = t.prevesy_vnitrni(m)
    assert 0 < sirka <= 16.0                 # plochá střecha, přemostění do 16 mm
    assert plocha_sikma == 0.0               # ostatní vnitřní převisy jsou ≤ 45°


def test_dute_vnejsi_tvar_nema_previsy(dute):
    _, m = dute
    assert t.prevesy(m)[1] == 0.0


def test_dute_trn_a_kapsa_jsou_stejne_jako_u_plne(dute):
    p, m = dute
    plne = t.vyrob(t.VELIKOSTI["M"])
    assert m.bounds[1][2] == pytest.approx(plne.bounds[1][2])
    assert m.extents == pytest.approx(plne.extents)
    z_hlava = p.vyska - p.h_kapsa + p.v_trn
    horni = m.vertices[np.isclose(m.vertices[:, 2], z_hlava, atol=1e-6)]
    assert np.hypot(horni[:, 0], horni[:, 1]).max() == pytest.approx(p.d_trn / 2 - p.srazeni, abs=0.001)


def test_dute_posun_trnu_posunuje_i_sloup():
    p = t.replace(t.VELIKOSTI["M_dute"], x_trn=10.0)
    m = t.vyrob(p)
    assert m.is_watertight
    assert m.contains([[10.0 + p.r_sloup - 1.0, 0, 15.0], [10.0 - (p.r_sloup - 1.0), 0, 15.0]]).all()


def test_dute_prilis_tenka_stena_je_odmitnuta():
    with pytest.raises(ValueError):
        t.vyrob(t.replace(t.VELIKOSTI["M_dute"], stena_dutiny=2.0))


def test_dute_prilis_siroka_strecha_je_odmitnuta():
    with pytest.raises(ValueError):
        t.vyrob(t.replace(t.VELIKOSTI["M_dute"], pol_strecha=9.0))


def test_dute_prilis_kratky_sloup_je_odmitnut():
    with pytest.raises(ValueError):
        t.vyrob(t.replace(t.VELIKOSTI["M_dute"], r_sloup=10.0))


@pytest.mark.parametrize("klic", ["S_dute", "L_dute"])
def test_dute_dalsi_velikosti(klic):
    m = t.vyrob(t.VELIKOSTI[klic])
    assert m.is_watertight
    assert len([c for c in m.split(only_watertight=False) if c.volume < 0]) == 2
    assert t.prevesy(m)[1] == 0.0
