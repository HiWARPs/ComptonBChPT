"""Minimal PDF content-stream interpreter to extract vector graphics (paths with colours) and text
positions from the Mathematica-exported figures of Krupina, Lensky, Pascalutsa PLB 782 (2018) 34.
Used to recover the data points actually used in that paper (digitised from vector coordinates)."""
import re
import zlib
import numpy as np


def streams(path):
    data = open(path, 'rb').read()
    out = []
    for s in re.findall(rb'stream\r?\n(.*?)\r?\nendstream', data, re.S):
        try:
            out.append(zlib.decompress(s).decode('latin-1'))
        except Exception:
            out.append(s.decode('latin-1'))
    return out


TOK = re.compile(r'\((?:\\.|[^\\)])*\)|\[[^\]]*\]|/[^\s/\[\]()<>]+|[^\s\[\]()<>/]+')


def interpret(content):
    toks = TOK.findall(content)
    stack = []
    ctm = np.eye(3)
    gs_stack = []
    state = dict(rg=(0, 0, 0), RG=(0, 0, 0), w=1.0)
    path = []      # list of subpaths, each list of points
    cur = None
    items = []
    texts = []
    tm = np.eye(3)
    in_text = False

    def T(x, y):
        v = np.array([x, y, 1.0]) @ ctm
        return (v[0], v[1])

    for t in toks:
        if re.fullmatch(r'[-+]?(\d+\.?\d*|\.\d+)', t):
            stack.append(float(t))
            continue
        if t.startswith('/') or t.startswith('(') or t.startswith('['):
            stack.append(t)
            continue
        op = t
        if op == 'q':
            gs_stack.append((ctm.copy(), dict(state)))
        elif op == 'Q':
            ctm, state = gs_stack.pop()
        elif op == 'cm':
            a, b, c, d, e, f = stack[-6:]
            M = np.array([[a, b, 0], [c, d, 0], [e, f, 1]])
            ctm = M @ ctm
        elif op == 'rg':
            state['rg'] = tuple(stack[-3:])
        elif op == 'RG':
            state['RG'] = tuple(stack[-3:])
        elif op == 'g':
            state['rg'] = (stack[-1],) * 3
        elif op == 'G':
            state['RG'] = (stack[-1],) * 3
        elif op == 'w':
            state['w'] = stack[-1]
        elif op == 'm':
            cur = [T(*stack[-2:])]
            path.append(cur)
        elif op == 'l':
            if cur is None:
                cur = []
                path.append(cur)
            cur.append(T(*stack[-2:]))
        elif op == 'c':
            if cur is None:
                cur = []
                path.append(cur)
            cur.append(T(*stack[-2:]))
        elif op in ('v', 'y'):
            cur.append(T(*stack[-2:]))
        elif op == 'h':
            pass
        elif op == 're':
            x, y, w_, h_ = stack[-4:]
            cur = [T(x, y), T(x + w_, y), T(x + w_, y + h_), T(x, y + h_)]
            path.append(cur)
        elif op in ('S', 's', 'f', 'f*', 'F', 'B', 'B*', 'b', 'b*', 'n'):
            if op != 'n' and path:
                items.append(dict(op=op, fill=state['rg'], stroke=state['RG'], w=state['w'],
                                  sub=[np.array(sp) for sp in path if len(sp)]))
            path = []
            cur = None
        elif op in ('W', 'W*'):
            pass
        elif op == 'BT':
            in_text = True
            tm = np.eye(3)
        elif op == 'ET':
            in_text = False
        elif op == 'Tm':
            a, b, c, d, e, f = stack[-6:]
            tm = np.array([[a, b, 0], [c, d, 0], [e, f, 1]])
        elif op == 'Td':
            tx, ty = stack[-2:]
            tm = np.array([[1, 0, 0], [0, 1, 0], [tx, ty, 1]]) @ tm
        elif op in ('Tj', 'TJ'):
            s = stack[-1] if stack else ''
            pos = np.array([0, 0, 1.0]) @ tm @ ctm
            texts.append((s, pos[0], pos[1]))
        stack = []
    return items, texts
