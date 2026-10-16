"""Fetch BibTeX entries (INSPIRE texkeys) for the paper's references from the INSPIRE API.
Usage: python3 -I fetch_bib.py out.bib"""
import sys, json, urllib.request, urllib.parse, time
Q = [
 'arxiv:0803.4115','arxiv:0907.0451','arxiv:1510.02794','arxiv:1401.3320','arxiv:1712.05349',
 'arxiv:1508.07952','arxiv:1604.00789','arxiv:1512.03765','arxiv:1203.6834','arxiv:1210.4104',
 'arxiv:hep-ph/0212124','arxiv:hep-ph/0501167','arxiv:1401.0140','arxiv:1903.07952','arxiv:1711.07401',
 'arxiv:2110.15691','arxiv:2205.10533','arxiv:1408.1576','arxiv:1611.03769','arxiv:1909.02032',
 'arxiv:nucl-ex/9507001','arxiv:hep-ph/0212024','arxiv:hep-ph/0508060','arxiv:hep-ph/0512244',
 'arxiv:hep-ph/0609004','arxiv:1506.01343','arxiv:1904.10581','arxiv:1511.01952','arxiv:1711.11546',
 'arxiv:2007.08438','arxiv:2005.09518','arxiv:2006.08626','arxiv:hep-ph/9803347','arxiv:nucl-th/0307070',
 'arxiv:1212.5498','arxiv:2609.15154','arxiv:nucl-ex/0210020','arxiv:nucl-ex/0003002',
 'arxiv:1912.06470',  # Pedroni Sconfietti JPG 47 (2020)?
 't "Bootstrap-based" and a Mornacchi and date 2022',
 'j "Eur.Phys.J.A,10,207"','j "Phys.Rev.Lett.,67,1511"','j "Phys.Lett.B,278,34"','j "Phys.Rev.C,48,1497"',
 'j "Phys.Lett.B,52,122"','j "Sov.J.Nucl.Phys.,21,355"','j "Phys.Rev.,110,733"','j "Phys.Rev.,105,982"','j "Phys.Rev.Lett.,3,93"',
 'j "Nucl.Phys.,18,473"','j "Nuovo Cim.,18,1203"','j "Phys.Rev.C,64,025203"','j "Nucl.Phys.A,620,385"',
 'j "Eur.Phys.J.A,5,269"','j "Phys.Rev.C,65,032202"','j "Eur.Phys.J.A,12,231"','j "Phys.Lett.B,503,245"',
 'j "Nucl.Phys.,18,310"','j "Phys.Rev.D,110,030001"','j "Phys.Rev.D,73,034003"','j "Phys.Rev.Lett.,95,232001"',
 'j "Phys.Rev.Lett.,67,1515"','j "Nucl.Phys.B,373,346"','j "Int.J.Mod.Phys.E,4,193"','j "Phys.Rev.,96,1433"',
 'j "Phys.Rev.Lett.,133,141901"','j "Phys.Rev.D,101,094502"','j "Phys.Rev.C,58,1013"','j "Phys.Rev.D,55,5598"',
 'j "Nucl.Phys.B,80,480"','j "Phys.Rev.C,76,015203"','j "J.Phys.G,47,054001"','j "Eur.Phys.J.A,57,81"',
 'j "Phys.Rev.Lett.,113,262506"','j "Phys.Rev.C,92,024005"',
]
out = []
for q in Q:
    url = 'https://inspirehep.net/api/literature?' + urllib.parse.urlencode({'q': q, 'size': 3, 'format': 'bibtex'})
    try:
        txt = urllib.request.urlopen(url, timeout=30).read().decode()
    except Exception as e:
        txt = ''
        print('FAIL', q, e, file=sys.stderr)
    ents = [e for e in txt.split('\n@') if e.strip()]
    if not ents:
        print('NONE', q, file=sys.stderr)
    else:
        e = ents[0] if ents[0].startswith('@') else '@' + ents[0]
        key = e.split('{', 1)[1].split(',', 1)[0]
        title = [l for l in e.splitlines() if l.strip().startswith('title')]
        print(f'{q:45s} -> {key:28s} {title[0].strip()[:90] if title else ""} {"(+%d more)" % (len(ents)-1) if len(ents) > 1 else ""}', file=sys.stderr)
        out.append(e.strip())
    time.sleep(0.4)
open(sys.argv[1], 'w').write('\n\n'.join(out) + '\n')
