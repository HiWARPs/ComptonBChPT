#!/usr/bin/env python3
"""Build db_new.csv: new (2014-2026) proton real-Compton data from A2@MAMI and HIGS.

Run as:  python3 -I build_db_new.py <DATA_DIR>
where DATA_DIR = <release>/code/data

All numerical values are transcribed verbatim (as strings) from the source files;
only theta_cm (resp. theta_lab) is computed, with exact 2-body Compton kinematics.
Sources are read as plain text (no code from downloads is executed).
"""
import csv
import math
import os
import re
import sys

MP = 938.272  # MeV, proton mass used for kinematics

HEADER = ["dataset", "texkey", "obs", "E_lab_MeV", "dE_MeV", "theta_lab_deg",
          "theta_cm_deg", "frame_of_cs", "value", "stat_err", "sys_err_pt",
          "norm_err_frac", "units", "source", "notes"]


def lab_to_cm(E, th_lab_deg):
    """Photon lab scattering angle -> c.m. angle for gamma p -> gamma p at lab photon energy E."""
    beta = E / (E + MP)  # velocity of the c.m. frame
    c = math.cos(math.radians(th_lab_deg))
    ccm = (c - beta) / (1.0 - beta * c)
    return math.degrees(math.acos(max(-1.0, min(1.0, ccm))))


def cm_to_lab(E, th_cm_deg):
    beta = E / (E + MP)
    c = math.cos(math.radians(th_cm_deg))
    cl = (c + beta) / (1.0 + beta * c)
    return math.degrees(math.acos(max(-1.0, min(1.0, cl))))


def fmt(x, nd=3):
    return f"{x:.{nd}f}"


def read(path):
    with open(path, encoding="utf-8", errors="replace") as f:
        return f.read()


rows = []


def add(**kw):
    r = {k: "" for k in HEADER}
    r.update(kw)
    rows.append(r)


def main(D):
    raw = os.path.join(D, "raw", "new")
    ext = os.path.join(D, "extracted")

    # ------------------------------------------------------------------
    # 1) Martel et al. (A2) PRL 114 (2015) 112501 -- Sigma_2x, 273-303 MeV
    #    Values: Martel PhD thesis (UMass 2013) Table 7.3, "Summed" column (phi-summing
    #    method).  The PRL Fig.3 EPS is titled "AsySum.eps: Summed Asymmetry" and its
    #    vector coordinates reproduce these numbers to <=0.001 (checked in verify script).
    # ------------------------------------------------------------------
    txt = read(os.path.join(ext, "Martel2013_thesis.txt"))
    i = txt.index("Angle        Fitted           Summed")
    j = txt.index("Table 7.3: Asymmetry results", i)
    block = txt[i:j]
    pat = re.compile(r"^\s*(\d+)\s+(-?\d\.\d+)\s*±\s*(\d\.\d+)\s+(-?\d\.\d+)\s*±\s*(\d\.\d+)", re.M)
    mrt = pat.findall(block)
    assert len(mrt) == 4, mrt
    E, dE = 288.0, 30.0
    for th, vf, ef, vs, es in mrt:
        thl = float(th)
        add(dataset="A2_Martel2015_S2x", texkey="A2:2014iky", obs="S2X",
            E_lab_MeV=fmt(E, 1), dE_MeV=fmt(dE, 1), theta_lab_deg=fmt(thl, 1),
            theta_cm_deg=fmt(lab_to_cm(E, thl)), frame_of_cs="na",
            value=vs, stat_err=es, sys_err_pt="0", norm_err_frac="0.021", units="1",
            source="Martel PhD thesis UMass 2013 (pmartel_dis.pdf) Table 7.3 'Summed' col.; = points of PRL114,112501 Fig.3 (AsySum.eps)",
            notes=(f"orig frame=lab; Egamma bin 273-303 MeV; theta_lab bin {thl-10:.0f}-{thl+10:.0f} deg; "
                   "stat_err = PRL error bar = stat (+) beam- and target-polarisation syst; "
                   "NORM_ALREADY_IN_STAT_ERR (P_T=81.6+-1.7%, P_e=81.9+-0.1%); "
                   "carbon-subtraction syst (<=10% of total) not included; "
                   f"thesis phi-fit alternative: {vf}+-{ef}"))

    # ------------------------------------------------------------------
    # 3) Paudyal et al. (A2) PRC 102 (2020) 035205 -- Sigma_2z, Table I (TeX source)
    # ------------------------------------------------------------------
    tex = read(os.path.join(raw, "Paudyal2020", "src", "PRC_S2z_Paudyal.tex"))
    i = tex.index(r"\label{table:Asy}")
    tab = tex[tex.rindex(r"\begin{tabular}", 0, i):i]
    lines = [l for l in tab.splitlines() if "^{\\circ}" in l]
    assert len(lines) == 10, len(lines)
    ebins = [(265, 285)] * 5 + [(285, 305)] * 5
    pat = re.compile(r"\$\s*(\d+\.\d)\^\{\\circ\}\$\s*&\s*(-?\d\.\d+)\s*&\s*\$\\pm\s*(\d\.\d+)\$\s*&\s*\$\\pm\s*(\d\.\d+)\$")
    for (e1, e2), l in zip(ebins, lines):
        m = pat.search(l)
        assert m, l
        th, v, rnd, sy = m.groups()
        E = 0.5 * (e1 + e2)
        thl = float(th)
        add(dataset="A2_Paudyal2020_S2z", texkey="A2:2019bqm", obs="S2Z",
            E_lab_MeV=fmt(E, 1), dE_MeV=fmt(e2 - e1, 1), theta_lab_deg=th,
            theta_cm_deg=fmt(lab_to_cm(E, thl)), frame_of_cs="na",
            value=v, stat_err=rnd, sys_err_pt=sy, norm_err_frac="0.104", units="1",
            source="PRC102,035205 Table I (arXiv:1909.02032v3 TeX source PRC_S2z_Paudyal.tex)",
            notes=(f"orig frame=lab; Egamma bin {e1}-{e2} MeV; theta_lab bin {thl-7.5:.0f}-{thl+7.5:.0f} deg; "
                   "stat_err = 'Rand.' column = stat (+) random syst; "
                   "sys_err_pt = 'Syst.' column = CORRELATED syst = target pol 10% (+) beam pol 2.7% (+) carbon subtr. 3-6%; "
                   "NORM_ALREADY_IN_SYS_ERR_PT (norm_err_frac=sqrt(0.10^2+0.027^2))"))

    # ------------------------------------------------------------------
    # 4) Mornacchi et al. (A2) PRL 128 (2022) 132503 -- Supplemental Material
    #    Tables I (dsigma/dOmega, 60 pts) and II (Sigma_3, 36 pts), lab frame
    # ------------------------------------------------------------------
    sup = read(os.path.join(ext, "Mornacchi2022_suppl_265979Suppl.txt"))
    # Table I
    i = sup.index("TABLE I: Measured unpolarized")
    j = sup.index("FIG. 1:", i)
    t1 = sup[i:j]
    hdr = re.search(r"86\.3 − 98\.2\s+98\.1 − 108\.4\s+108\.5 − 118\.7\s+118\.7 − 130\.2\s+130\.3 − 140\.4", t1)
    assert hdr
    cs_bins = [(86.3, 98.2, "92.2"), (98.1, 108.4, "103.2"), (108.5, 118.7, "113.6"),
               (118.7, 130.2, "124.1"), (130.3, 140.4, "135.3")]
    trip = r"(-?\d+\.\d+)\s*±\s*(\d+\.\d+)\s*±\s*(\d+\.\d+)"
    rowpat = re.compile(r"^\D*?(\d{2,3})\s+" + r"\s+".join([trip] * 5) + r"\s*$", re.M)
    found = rowpat.findall(t1)
    assert len(found) == 12, len(found)
    for g in found:
        th = g[0]
        for k, (e1, e2, ec) in enumerate(cs_bins):
            v, st, sy = g[1 + 3 * k: 4 + 3 * k]
            E = float(ec)
            thl = float(th)
            add(dataset="A2_Mornacchi2022_CS", texkey="A2:2021wmm", obs="CS",
                E_lab_MeV=ec, dE_MeV=fmt(e2 - e1, 1), theta_lab_deg=fmt(thl, 1),
                theta_cm_deg=fmt(lab_to_cm(E, thl)), frame_of_cs="lab",
                value=v, stat_err=st, sys_err_pt=sy, norm_err_frac="0.03", units="nb/sr",
                source="PRL128,132503 Suppl. Material Table I (eprints.gla.ac.uk/265979/2/265979Suppl.pdf)",
                notes=(f"orig frame=lab (angles and dsigma/dOmega); Egamma bin {e1}-{e2} MeV "
                       f"(E_lab = bin central value quoted in Mornacchi PhD thesis App.A; bin midpoint {0.5*(e1+e2):.2f}); "
                       f"theta_lab bin {thl-5:.0f}-{thl+5:.0f} deg; combined 2 beamtimes (2018); "
                       "sys_err_pt = point-to-point uncorrelated syst; norm 3% = target density 1% (+) flux 2% (+) cuts/MC 2%"))
    # Table II
    i = sup.index("TABLE II: Measured proton Compton scattering beam asymmetry")
    t2 = sup[i:]
    s3_bins = [(86.3, 98.2, "92.2"), (98.1, 118.7, "108.5"), (118.7, 140.4, "129.5")]
    trip = r"([−-]?\d+\.\d+)\s*±\s*(\d+\.\d+)\s*±\s*(\d+\.\d+)"
    rowpat = re.compile(r"^\D*?(\d{2,3})\s+" + r"\s+".join([trip] * 3) + r"\s*$", re.M)
    found = rowpat.findall(t2)
    assert len(found) == 12, len(found)
    for g in found:
        th = g[0]
        for k, (e1, e2, ec) in enumerate(s3_bins):
            v, st, sy = g[1 + 3 * k: 4 + 3 * k]
            v = v.replace("−", "-")
            E = float(ec)
            thl = float(th)
            add(dataset="A2_Mornacchi2022_S3", texkey="A2:2021wmm", obs="S3",
                E_lab_MeV=ec, dE_MeV=fmt(e2 - e1, 1), theta_lab_deg=fmt(thl, 1),
                theta_cm_deg=fmt(lab_to_cm(E, thl)), frame_of_cs="na",
                value=v, stat_err=st, sys_err_pt=sy, norm_err_frac="0.05", units="1",
                source="PRL128,132503 Suppl. Material Table II (eprints.gla.ac.uk/265979/2/265979Suppl.pdf)",
                notes=(f"orig frame=lab; Egamma bin {e1}-{e2} MeV (E_lab = thesis central value; midpoint {0.5*(e1+e2):.2f}); "
                       f"theta_lab bin {thl-5:.0f}-{thl+5:.0f} deg; combined 2 beamtimes; "
                       "Sigma3=(dsig_par-dsig_perp)/(dsig_par+dsig_perp) [par = photon E-vector in scattering plane] -> NEGATIVE below threshold; "
                       "norm_err_frac 0.05 is a UNIFORMLY distributed correlated scale (half-width; rms = 0.05/sqrt3 = 0.0289) from lin.-pol. determination; "
                       "sys_err_pt = point-to-point uncorrelated (background)"))

    # ------------------------------------------------------------------
    # 5) X. Li et al. (HIGS) PRL 128 (2022) 132502 -- Table I (TeX source)
    # ------------------------------------------------------------------
    tex = read(os.path.join(raw, "Li2022", "src", "prl_proof_Apr2022.tex"))
    i = tex.index(r"\label{tab:DXS}")
    tab = tex[tex.rindex(r"\begin{tabular*}", 0, i):i]
    i_lin = tab.index("Linearly polarized gamma-ray beam")
    pat = re.compile(r"\$(\d+)\^\\circ\$\s*&\s*\$(\d+)\^\\circ\$\s*&(\d+\.\d+)\s*&\$\\pm(\d+\.\d+)\$\s*&\$\\pm(\d+\.\d+)\$\s*&\$\\pm(\d+\.\d+)\$")
    n = 0
    for m in pat.finditer(tab):
        th, phi, v, st, pp, nrm = m.groups()
        lin = m.start() > i_lin
        thl = float(th)
        if lin:
            E, dE = 83.4, 0.027 * 83.4
            obs = "CS_PAR" if phi in ("0", "180") else "CS_PERP"
            beam = "linearly polarised beam (pol. axis horizontal, phi=0/180 in-plane), nearly 100% polarised, no pol. correction quoted"
        else:
            E, dE = 81.3, 0.065 * 81.3
            obs = "CS"
            beam = "circularly polarised beam on unpolarised target = unpolarised dsigma/dOmega"
        add(dataset="HIGS_Li2022_" + ("LIN" if lin else "CIRC"), texkey="Li:2022vnz", obs=obs,
            E_lab_MeV=fmt(E, 1), dE_MeV=fmt(dE, 2), theta_lab_deg=fmt(thl, 1),
            theta_cm_deg=fmt(lab_to_cm(E, thl)), frame_of_cs="lab",
            value=v, stat_err=st, sys_err_pt=pp, norm_err_frac="0.0224", units="nb/sr",
            source="PRL128,132502 Table I (arXiv:2205.10533 TeX source prl_proof_Apr2022.tex)",
            notes=(f"orig frame=lab; detector at theta_lab={th} deg, phi={phi} deg; {beam}; "
                   f"E_lab = intensity-weighted mean beam energy, dE = FWHM spread ({'2.7' if lin else '6.5'}%); "
                   f"tabulated normalisation syst = {nrm} nb/sr (= flux 2% (+) target thickness 1%); "
                   "values bin-centre corrected to detector centre with chiEFT model; "
                   "authors floated circ. and lin. normalisations independently"))
        n += 1
    assert n == 16, n

    out = os.path.join(D, "db_new.csv")
    with open(out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=HEADER, quoting=csv.QUOTE_MINIMAL)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    # summary
    from collections import Counter
    c = Counter((r["dataset"], r["obs"]) for r in rows)
    for k, v in sorted(c.items()):
        print(k, v)
    print("total rows", len(rows), "->", out)


if __name__ == "__main__":
    main(sys.argv[1])
