#!/usr/bin/env python3
"""Independent cross-checks of db_new.csv against second sources.

Run as:  python3 -I verify_db_new.py <DATA_DIR>
Checks (all read-only, plain-text parsing):
  K) kinematics: theta_cm via Lorentz boost vs. via invariant t (lab) = t (cm)
  M1) Martel S2x vs. vector coordinates of PRL Fig.3 (AsySum.eps) -- verification only
  P ) Paudyal S2z vs. pdftotext of arXiv PDF Table I
  Mo) Mornacchi CS/S3 (published, combined) vs. stat-weighted mean of the two
      per-beamtime tables in Mornacchi PhD thesis App. A/B
  L ) HIGS cross sections vs. Li PhD thesis Tables 6.5/6.7 and arXiv PDF text
"""
import csv
import math
import os
import re
import subprocess
import sys

MP = 938.272


def rows_of(D):
    with open(os.path.join(D, "db_new.csv"), encoding="utf-8") as f:
        return list(csv.DictReader(f))


def check_kin(rows):
    worst = 0.0
    for r in rows:
        E = float(r["E_lab_MeV"]); tl = float(r["theta_lab_deg"]); tc = float(r["theta_cm_deg"])
        wp = E / (1 + E / MP * (1 - math.cos(math.radians(tl))))  # scattered photon energy, lab
        t = -2 * E * wp * (1 - math.cos(math.radians(tl)))
        s = MP ** 2 + 2 * MP * E
        wcm = (s - MP ** 2) / (2 * math.sqrt(s))
        ccm = 1 + t / (2 * wcm ** 2)
        tc2 = math.degrees(math.acos(ccm))
        worst = max(worst, abs(tc2 - tc))
    print(f"[K] max |theta_cm(boost) - theta_cm(from t)| = {worst:.4f} deg over {len(rows)} rows")


def check_martel(D, rows):
    eps = open(os.path.join(D, "raw/new/Martel2015/src/Sig2x_gE1E1.eps"), encoding="latin-1").read()
    eps = re.sub(r"\s+", " ", eps)
    assert "AsySum.eps: Summed Asymmetry" in eps
    # y-axis major ticks (PS units): -0.6 .. 0.6 at these y values (from '373 y m -56 X s')
    ticks = [int(y) for y in re.findall(r"373 (\d+) m -56 X s", eps)]
    ticks = sorted(set(ticks))
    vals = [-0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6]
    assert len(ticks) == 7, ticks
    a = (vals[-1] - vals[0]) / (ticks[-1] - ticks[0])
    y0 = ticks[0]
    # x-axis: 0 deg at 317, 180 deg at 2154 (major ticks)
    xs = lambda x: (x - 317) * 180.0 / (2154 - 317)
    pts = re.findall(r"(\d+) (\d+) m20", eps)
    seen = {}
    for x, y in pts:
        seen[(int(x), int(y))] = 1
    up = dict(((int(x), int(y)), int(h)) for x, y, h in re.findall(r"(\d+) (\d+) m (\d+) Y s", eps))
    out = []
    for (x, y) in sorted(seen):
        th = xs(x); v = vals[0] + a * (y - y0)
        e = a * up.get((x, y), float("nan"))
        out.append((round(th, 1), v, e))
    mr = [r for r in rows if r["dataset"] == "A2_Martel2015_S2x"]
    worst = 0
    for (th, v, e), r in zip(out, mr):
        dv = abs(v - float(r["value"])); de = abs(e - float(r["stat_err"]))
        worst = max(worst, dv, de)
        print(f"[M1] th={th:6.1f}  EPS v={v:+.4f} e={e:.4f} | csv th={r['theta_lab_deg']} v={r['value']} e={r['stat_err']}")
    print(f"[M1] max abs diff (value/error) = {worst:.4f}  (EPS resolution ~{a:.4f})")


def pdftext(path):
    return subprocess.run(["/opt/homebrew/bin/pdftotext", "-layout", path, "-"],
                          capture_output=True, text=True).stdout


def check_paudyal(D, rows):
    t = pdftext(os.path.join(D, "raw/new/Paudyal2020/arXiv_1909.02032.pdf"))
    found = re.findall(r"(\d+\.\d)◦ (\d\.\d+) ±(\d\.\d+) ±(\d\.\d+)", t)
    pr = [r for r in rows if r["dataset"] == "A2_Paudyal2020_S2z"]
    ok = len(found) == len(pr) and all(
        (f[0], f[1], f[2], f[3]) == (r["theta_lab_deg"], r["value"], r["stat_err"], r["sys_err_pt"])
        for f, r in zip(found, pr))
    print(f"[P ] Paudyal: {len(found)} rows in arXiv PDF text; exact string match with CSV: {ok}")


def check_mornacchi(D, rows):
    t = open(os.path.join(D, "extracted/Mornacchi2021_thesis.txt"), encoding="utf-8").read()
    def tables(prefix, n):
        res = {}
        for k in range(1, n + 1):
            i = t.index(f"Table {prefix}.{k}: ")
            j = t.find("Table ", i + 20)
            j2 = t.find("Appendix", i + 20)
            seg = t[i: min(x for x in (j, j2) if x > 0)]
            bt = "March" if "March 2018" in seg[:200] else "July"
            for m in re.finditer(r"^\s*(\d+\.\d)\s+(\d+\.\d)\s+(-?\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(-?\d+\.\d+)\s*$", seg, re.M):
                E, th, v, st, sy, co = m.groups()
                res.setdefault((E, th), {})[bt] = (float(v), float(st), float(sy))
        return res
    A = tables("A", 10)
    Bt = tables("B", 6)
    for name, T, ds in (("CS", A, "A2_Mornacchi2022_CS"), ("S3", Bt, "A2_Mornacchi2022_S3")):
        worst_v = worst_e = 0; n = 0; big = []
        for r in rows:
            if r["dataset"] != ds:
                continue
            key = (r["E_lab_MeV"], r["theta_lab_deg"])
            d = T.get(key)
            if not d or len(d) != 2:
                big.append(("missing", key)); continue
            (v1, s1, _), (v2, s2, _) = d["March"], d["July"]
            w1, w2 = 1 / s1 ** 2, 1 / s2 ** 2
            vm = (v1 * w1 + v2 * w2) / (w1 + w2); sm = (w1 + w2) ** -0.5
            dv = abs(vm - float(r["value"])); de = abs(sm - float(r["stat_err"]))
            worst_v = max(worst_v, dv); worst_e = max(worst_e, de); n += 1
            if dv > 0.01 * max(1, abs(float(r["value"]))) or de > 0.01:
                big.append((key, r["value"], round(vm, 4), r["stat_err"], round(sm, 4)))
        print(f"[Mo] {name}: {n} points compared with stat-weighted mean of thesis March/July tables; "
              f"max|dv|={worst_v:.4f} max|dstat|={worst_e:.4f}; flagged={big[:8]}")


def check_li(D, rows):
    t = open(os.path.join(D, "extracted/Li2020_thesis.txt"), encoding="utf-8").read()
    def tab(label):
        i = t.index(label); seg = t[i:i + 2500]
        return re.findall(r"(\d+)◦\s+(\d+)◦\s+\w+\s+(\d+\.\d+)\s+±(\d+\.\d+)\s+±(\d+\.\d+)\s+±(\d+\.\d+)", seg)
    th = tab("Table 6.5: The Compton scattering cross sections of the proton") + \
         tab("Table 6.7: The Compton scattering cross section of the proton using linearly")
    lr = [r for r in rows if r["texkey"] == "Li:2022vnz"]
    ok = len(th) == 16 and all((f[2], f[3], f[4]) == (r["value"], r["stat_err"], r["sys_err_pt"]) for f, r in zip(th, lr))
    print(f"[L ] HIGS: {len(th)} rows in Li thesis Tables 6.5+6.7; exact match with CSV (value, stat, p2p): {ok}")
    pt = pdftext(os.path.join(D, "raw/new/Li2022/arXiv_2205.10533.pdf"))
    vals = re.findall(r"(\d+\.\d\d)\s+±(\d\.\d\d)\s+±(\d\.\d\d)\s+±(\d\.\d\d)", pt)
    ok2 = len(vals) == 16 and all((v[0], v[1], v[2]) == (r["value"], r["stat_err"], r["sys_err_pt"]) for v, r in zip(vals, lr))
    print(f"[L ] HIGS: {len(vals)} rows in arXiv PDF text; exact match with CSV: {ok2}")


if __name__ == "__main__":
    D = sys.argv[1]
    rows = rows_of(D)
    check_kin(rows)
    check_martel(D, rows)
    check_paudyal(D, rows)
    check_mornacchi(D, rows)
    check_li(D, rows)
