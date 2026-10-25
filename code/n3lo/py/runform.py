"""runform.py -- FORM projections for Graph objects (parallel), cached in out/."""
import os, sys, concurrent.futures as cf
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import formtools as ft

def project(g, gi=False):
    body, ks = g.form_body(gi=gi)
    d = ft.run_form(g.name + ('_GI' if gi else ''), body)
    return g.name, d

def project_all(graphs, nproc=8, gi=False):
    res = {}
    with cf.ThreadPoolExecutor(nproc) as ex:
        for name, d in ex.map(lambda g: project(g, gi), graphs):
            res[name] = d
    return res
