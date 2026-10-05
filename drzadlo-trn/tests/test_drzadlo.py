import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import drzadlo as d  # noqa: E402


@pytest.fixture(scope="module")
def tvar():
    p = d.Parametry()
    return p, d.vyrob(p), d.profil(p)


# --- pomocná funkce zaoblení ---------------------------------------------------------------

def test_zaobli_pravy_uhel_dava_ctvrtkruh():
    body = d._zaobli([(0, 0), (10, 0), (10, 10)], [0, 2.0, 0], n=50)
    oblouk = body[1:-1]
    stred = np.array([8.0, 2.0])
    assert np.linalg.norm(oblouk - stred, axis=1) == pytest.approx(2.0, abs=1e-9)
    assert oblouk[0] == pytest.approx([8.0, 0.0])        # tečný bod na první úsečce
    assert oblouk[-1] == pytest.approx([10.0, 2.0])      # tečný bod na druhé úsečce


def test_zaobli_nulovy_polomer_necha_roh():
    body = d._zaobli([(0, 0), (5, 0), (5, 5)], [0, 0, 0])
    assert body.tolist() == [[0, 0], [5, 0], [5, 5]]


def test_zaobli_oblouk_je_tecny_u_obou_stran():
    body = d._zaobli([(0, 0), (10, 0), (20, 10)], [0, 3.0, 0], n=60)   # 45° zlom
    oblouk = body[1:-1]
    assert np.allclose(oblouk[0][1], 0, atol=1e-9)                     # začíná na vodorovné úsečce
    sm = oblouk[-1] - oblouk[-2]
    assert np.degrees(np.arctan2(sm[1], sm[0])) == pytest.approx(45, abs=2)   # končí ve směru 45° úsečky


# --- rozměry zadání ------------------------------------------------------------------------

def test_trn_ma_prumer_9_a_vysku_12(tvar):
    p, m, (body, z_opora, z_hlava) = tvar
    assert z_hlava - z_opora == pytest.approx(12.0)
    sek = m.section(plane_origin=[0, 0, z_opora + 6.0], plane_normal=[0, 0, 1])
    x = sek.vertices[:, 0]
    assert x.max() - x.min() == pytest.approx(9.0, abs=0.01)
    assert np.hypot(sek.vertices[:, 0], sek.vertices[:, 1]).max() == pytest.approx(4.5, abs=0.001)


def test_hlava_trnu_je_plochá_se_srazenou_hranou(tvar):
    p, m, (_, z_opora, z_hlava) = tvar
    horni = m.vertices[np.isclose(m.vertices[:, 2], z_hlava, atol=1e-6)]
    assert len(horni) > 10
    assert np.hypot(horni[:, 0], horni[:, 1]).max() == pytest.approx(4.5 - p.srazeni, abs=0.001)   # plocha Ø 7
    pod = m.vertices[np.isclose(m.vertices[:, 2], z_hlava - p.srazeni, atol=1e-6)]
    assert np.hypot(pod[:, 0], pod[:, 1]).max() == pytest.approx(4.5, abs=0.001)    # sražení 1 × 45°


def test_hlava_trnu_je_nejvyssi_bod(tvar):
    _, m, (_, _, z_hlava) = tvar
    assert m.bounds[1][2] == pytest.approx(z_hlava)


def test_pata_trnu_je_zaoblena(tvar):
    p, m, (_, z_opora, _) = tvar
    # těsně nad opěrnou plochou je poloměr větší než 4,5 (zaoblení R2 rozšiřuje trn u paty)
    sek = m.section(plane_origin=[0, 0, z_opora + 0.2], plane_normal=[0, 0, 1])
    assert np.hypot(sek.vertices[:, 0], sek.vertices[:, 1]).max() > 4.9


def test_drzadlo_je_ploche_zespodu_a_mensi_nez_prumer(tvar):
    p, m, _ = tvar
    assert m.bounds[0][2] == pytest.approx(0.0, abs=1e-9)                  # dlaň na stole
    assert m.extents[0] == pytest.approx(p.d_dlane, abs=0.05)
    dno = m.vertices[np.isclose(m.vertices[:, 2], 0.0, atol=1e-9)]
    assert np.hypot(dno[:, 0], dno[:, 1]).max() == pytest.approx(p.r_plocha - p.r_zaobleni_dlane * np.tan(np.radians(22.5)), abs=0.05)


def test_opěrná_plocha_je_vodorovna_o_prumeru_34(tvar):
    p, m, (_, z_opora, _) = tvar
    r = np.hypot(m.vertices[:, 0], m.vertices[:, 1])
    na_opore = m.vertices[np.isclose(m.vertices[:, 2], z_opora, atol=1e-6)]
    assert len(na_opore) > 10
    assert np.hypot(na_opore[:, 0], na_opore[:, 1]).max() > p.d_opora / 2 - 1.0


# --- kvalita sítě a tisk ---------------------------------------------------------------------

def test_sit_je_uzavrena_a_ma_kladny_objem(tvar):
    _, m, _ = tvar
    assert m.is_watertight and m.is_winding_consistent and m.volume > 0


def test_objem_odpovida_rotacnimu_telesu(tvar):
    p, m, (body, _, _) = tvar
    r, z = body[:, 0], body[:, 1]
    plocha_pul = 0.5 * np.abs(np.dot(r, np.roll(z, -1)) - np.dot(z, np.roll(r, -1)))   # shoelace
    teziste_r = (1 / (6 * plocha_pul)) * np.sum((r + np.roll(r, -1)) * (r * np.roll(z, -1) - np.roll(r, -1) * z))
    pappus = 2 * np.pi * abs(teziste_r) * plocha_pul
    assert m.volume == pytest.approx(pappus, rel=0.005)


def test_zadne_previsy_nad_45_stupnu_krome_zaobleni_u_stolu(tvar):
    _, m, _ = tvar
    podil, plocha = d.preruseni_prevesu(m, 45.0, 1.5)
    assert plocha == 0.0


def test_vyska_a_hmotnost_rozumne(tvar):
    _, m, _ = tvar
    assert 40 < m.extents[2] < 55
    assert 100 < m.volume / 1000 < 150            # objem plného tělesa v ml


def test_prilis_kratky_okraj_je_odmitnut():
    with pytest.raises(ValueError):
        d.profil(d.Parametry(d_dlane=100.0, r_plocha=35.0))      # z_c = 15, okraj do 21: nevejdou se oblouky


def test_srazeni_vetsi_nez_polomer_je_odmitnuto():
    with pytest.raises(ValueError):
        d.profil(d.Parametry(srazeni=5.0))


@pytest.mark.parametrize("prumer", [70, 90])
def test_dalsi_velikosti_dlane_jsou_validni(prumer):
    k = prumer / 80
    m = d.vyrob(d.Parametry(d_dlane=prumer, r_plocha=28 * k, z_okraj_konec=21 * k))
    assert m.is_watertight
    assert m.extents[0] == pytest.approx(prumer, abs=0.05)
    assert d.preruseni_prevesu(m)[1] == 0.0
