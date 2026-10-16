"""Crude Mathematica-notebook -> linear InputForm text converter (for reading Krupina's PWA notebooks).
Concatenates the string tokens of each Input / Output cell's box structure.  Usage: nb2txt.py in.nb out.txt"""
import re
import sys

TOKEN = re.compile(r'"((?:\\.|[^"\\])*)"')


def unesc(s):
    s = s.replace('\\"', '"').replace('\\<', '').replace('\\>', '').replace('\\[Equal]', '==')
    s = s.replace('\\[IndentingNewLine]', '\n').replace('\\n', '\n').replace('\\[Rule]', '->')
    s = s.replace('\\[Alpha]', 'alpha').replace('\\[Beta]', 'beta').replace('\\[Gamma]', 'gamma')
    s = s.replace('\\[Theta]', 'theta').replace('\\[Pi]', 'Pi').replace('\\[Chi]', 'chi').replace('\\[Sigma]', 'sigma')
    s = s.replace('\\[Omega]', 'omega').replace('\\[Nu]', 'nu').replace('\\[Kappa]', 'kappa').replace('\\[Mu]', 'mu')
    s = s.replace('\\[Infinity]', 'Infinity').replace('\\[LessEqual]', '<=').replace('\\[GreaterEqual]', '>=')
    s = s.replace('\\[Element]', ' in ').replace('\\[Times]', '*').replace('\\[CapitalGamma]', 'Gamma')
    s = s.replace('\\[Lambda]', 'lambda').replace('\\[Epsilon]', 'eps').replace('\\[Delta]', 'delta')
    return s


def cells(text):
    # split at top-level "Cell[" occurrences; keep type tag
    out = []
    for m in re.finditer(r'Cell\[BoxData\[(.*?)\],\s*"(Input|Output|Print)"', text, re.S):
        body, kind = m.group(1), m.group(2)
        # drop CellChangeTimes etc. (not inside BoxData)
        toks = TOKEN.findall(body)
        s = ''.join(unesc(t) for t in toks)
        if kind == 'Output' and len(s) > 3000:
            s = s[:3000] + ' ...[truncated]'
        out.append((kind, s))
    return out


if __name__ == '__main__':
    txt = open(sys.argv[1], encoding='latin-1').read()
    with open(sys.argv[2], 'w') as fh:
        for i, (k, s) in enumerate(cells(txt)):
            fh.write(f'\n### [{i}] {k}\n{s}\n')
