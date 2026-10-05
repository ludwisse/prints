import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import drzadlo_tyc as t  # noqa: E402


@pytest.fixture(scope="module")
def tvar():
    p = t.Parametry()
    return p, t.vyrob(p)


# --- rozměry zadání ------------------------------------------------------------------------

def test_trn_prumer_9_a_vyska_12_nad_spodkem(tvar):
    p, m = tvar
    assert m.bounds[1][2] == pytest.approx(p.vyska + 12.0)
    sek = m.section(plane_origin=[0, 0, p.vyska + 6.0], plane_normal=[0, 0, 1])
    v = sek.vertices
    assert np.hypot(v[:, 0], v[:, 1]).max() == pytest.approx(4.5, abs=0.001)
    assert v[:, 1].max() - v[:, 1].min() == pytest.approx(9.0, abs=0.01)


def test_hlava_trnu_je_plochá_se_srazenou_hranou(tvar):
    p, m = tvar
    z_hlava = p.vyska + p.v_trn
    horni = m.vertices[np.isclose(m.vertices[:, 2], z_hlava, atol=1e-6)]
    assert len(horni) > 10
    assert np.hypot(horni[:, 0], horni[:, 1]).max() == pytest.approx(4.5 - p.srazeni, abs=0.001)
    pod = m.vertices[np.isclose(m.vertices[:, 2], z_hlava - p.srazeni, atol=1e-6)]
    assert np.hypot(pod[:, 0], pod[:, 1]).max() == pytest.approx(4.5, abs=0.001)


def test_pata_trnu_je_zaoblena(tvar):
    p, m = tvar
    sek = m.section(plane_origin=[0, 0, p.vyska + 0.2], plane_normal=[0, 0, 1])
    assert np.hypot(sek.vertices[:, 0], sek.vertices[:, 1]).max() > 4.9


def test_rozmery_tyce(tvar):
    p, m = tvar
    assert m.extents[0] == pytest.approx(p.delka, abs=0.05)
    assert m.extents[1] == pytest.approx(p.hloubka, abs=0.05)
    assert m.bounds[0][2] == pytest.approx(0.0, abs=1e-9)                       # dlaň na stole
    assert m.bounds[1][2] == pytest.approx(p.vyska + p.v_trn)


def test_trn_je_uprostred_tyce_a_posun_funguje(tvar):
    p, m = tvar
    hlava = m.vertices[np.isclose(m.vertices[:, 2], m.bounds[1][2], atol=1e-6)]
    assert hlava[:, 0].mean() == pytest.approx(0.0, abs=0.01)
    m2 = t.vyrob(t.Parametry(x_trn=10.0))
    hlava2 = m2.vertices[np.isclose(m2.vertices[:, 2], m2.bounds[1][2], atol=1e-6)]
    assert hlava2[:, 0].mean() == pytest.approx(10.0, abs=0.01)


def test_spodek_je_rovny_kolem_trnu(tvar):
    p, m = tvar
    # v z = vyska je plochý pruh široký 2 * (hloubka/2 - r_spodek)
    na_spodku = m.vertices[np.isclose(m.vertices[:, 2], p.vyska, atol=1e-6)]
    assert len(na_spodku) > 10
    assert np.abs(na_spodku[:, 1]).max() == pytest.approx(p.hloubka / 2 - p.r_spodek, abs=0.01)


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
    assert pul[0] == pytest.approx([0, 0]) and pul[-1] == pytest.approx([0, 28.0])
    assert (pul[:, 0] >= -1e-9).all() and pul[:, 0].max() == pytest.approx(20.0, abs=1e-6)


def test_prilis_kratka_stena_je_odmitnuta():
    with pytest.raises(ValueError):
        t.pulprofil(t.Parametry(vyska=14.0))


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
    assert m.bounds[1][2] == pytest.approx(p.vyska + 12.0)
