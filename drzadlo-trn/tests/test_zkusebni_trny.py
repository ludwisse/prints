import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import zkusebni_trny as z  # noqa: E402


def test_destička_je_uzavrena_a_ma_pet_trnu():
    m = z.vyrob()
    assert m.is_watertight and m.volume > 0
    sek = m.section(plane_origin=[0, 0, z.DESKA[2] + 4.0], plane_normal=[0, 0, 1])
    assert len(sek.entities) == len(z.PRUMERY)


@pytest.mark.parametrize("i", range(5))
def test_prumery_trnu_odpovidaji(i):
    m = z.vyrob()
    x0 = -z.ROZTEC * (len(z.PRUMERY) - 1) / 2
    stred = x0 + i * z.ROZTEC
    sek = m.section(plane_origin=[0, 0, z.DESKA[2] + 4.0], plane_normal=[0, 0, 1])
    v = sek.vertices
    v = v[np.abs(v[:, 0] - stred) < z.ROZTEC / 2]
    assert np.hypot(v[:, 0] - stred, v[:, 1]).max() == pytest.approx(z.PRUMERY[i] / 2, abs=0.001)


def test_delka_trnu_a_sraz():
    m = z.vyrob()
    assert m.bounds[1][2] == pytest.approx(z.DESKA[2] + z.DELKA)
