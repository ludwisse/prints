import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import drzadlo_tyc as t  # noqa: E402
import vykres  # noqa: E402


def test_koty_odpovidaji_zadani():
    k = vykres.koty(t.Parametry())
    assert k["d_trn"] == 14.0 and k["delka_trnu"] == 9.4
    assert k["d_kapsa"] == 20.5
    assert k["prstenec_sirka"] == pytest.approx(3.25)
    assert k["presah_pod_spodek"] == pytest.approx(9.4 - 1.5)
    assert k["celkova_vyska"] == pytest.approx(30 - 1.5 + 9.4)
    assert k["delka"] == 105 and k["hloubka"] == 44 and k["vyska"] == 30
    assert k["plocha_dlane"] == 24 and k["rovny_spodek"] == 30
    assert k["d_hlavy_plocha"] == pytest.approx(12.0)


@pytest.mark.parametrize("velikost", ["S", "L"])
def test_koty_dalsich_velikosti_jsou_konzistentni(velikost):
    p = t.VELIKOSTI_PLNE[velikost]
    k = vykres.koty(p)
    assert k["celkova_vyska"] == pytest.approx(p.vyska - p.h_kapsa + p.v_trn)
    assert k["rovny_spodek"] == pytest.approx(p.hloubka - 2 * p.r_spodek)


def test_vykres_vytvori_png_a_pdf(tmp_path):
    png = tmp_path / "v.png"
    vykres.kresli(t.Parametry(), str(png), "M")
    assert png.exists() and png.stat().st_size > 20000
    assert (tmp_path / "v.pdf").exists() and (tmp_path / "v.pdf").stat().st_size > 5000
