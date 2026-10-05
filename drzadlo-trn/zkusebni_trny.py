"""Zkušební destička s trny různých průměrů pro ověření lícování v otvoru dílu se závitem 9,29 mm.

Pět trnů o délce 8 mm, průměry 15,2 / 15,4 / 15,6 / 15,8 / 16,0 mm (zleva doprava), stejně sražené hlavy.
Vyzkoušej, který se zasune do otvoru s lehkým odporem a bez vůle, a ten průměr dej do `d_trn` v drzadlo_tyc.py.
Použití: python zkusebni_trny.py
"""
import os

import manifold3d as mf
import numpy as np
import trimesh

PRUMERY = (15.2, 15.4, 15.6, 15.8, 16.0)
DELKA = 8.0
SRAZENI = 1.0
ROZTEC = 19.0
DESKA = (ROZTEC * len(PRUMERY) + 4.0, 22.0, 3.0)


def trn(d, x):
    r = d / 2
    body = np.array([(0, -1.0), (r, -1.0), (r, DELKA - SRAZENI), (r - SRAZENI, DELKA), (0, DELKA)], float)
    return mf.Manifold.revolve(mf.CrossSection([body]), 192).translate((x, 0, DESKA[2]))


def vyrob():
    deska = mf.Manifold.cube(DESKA, center=True).translate((0, 0, DESKA[2] / 2))
    x0 = -ROZTEC * (len(PRUMERY) - 1) / 2
    m = deska
    for i, d in enumerate(PRUMERY):
        m = m + trn(d, x0 + i * ROZTEC)
    g = m.to_mesh()
    return trimesh.Trimesh(np.asarray(g.vert_properties)[:, :3], np.asarray(g.tri_verts), process=True)


if __name__ == "__main__":
    m = vyrob()
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "zkusebni_trny.stl")
    m.export(out)
    print(out, m.is_watertight, m.extents.round(1))
