"""ltpy.py -- Python access to LoopTools-dq tensor coefficients (via libltwrap.dylib)."""
import ctypes, os, numpy as np
_here = os.path.dirname(os.path.abspath(__file__))
_lib = ctypes.CDLL(os.path.join(_here, 'libltwrap.dylib'))
_names = {}
for line in open(os.path.join(_here, 'coeffs.txt')):
    k, *nm = line.split()
    _names[k] = nm
_d = ctypes.c_double
def _p(x): return ctypes.byref(_d(x))
_lib.lt_init()
_lib.lt_setversion(ctypes.byref(ctypes.c_int(64)))   # Denner D0 (FF fails for D(0,M2,0,M2;u,s;m,m,M,M) below threshold)
def setmudim(m2): _lib.lt_setmudim(_p(m2))
def setdelta(d): _lib.lt_setdelta(_p(d))
def clear(): _lib.lt_clear()
def _get(kind, fn, args):
    n = len(_names[kind])
    fin = np.zeros(n, complex); uv = np.zeros(n, complex)
    fn(*[_p(a) for a in args], fin.ctypes.data_as(ctypes.c_void_p), uv.ctypes.data_as(ctypes.c_void_p))
    return dict(zip(_names[kind], fin)), dict(zip(_names[kind], uv))
def B(p2, m1, m2): return _get('B', _lib.lt_bget, (p2, m1, m2))
def C(p1, p2, p3, m1, m2, m3): return _get('C', _lib.lt_cget, (p1, p2, p3, m1, m2, m3))
def D(p1, p2, p3, p4, p12, p23, m1, m2, m3, m4):
    return _get('D', _lib.lt_dget, (p1, p2, p3, p4, p12, p23, m1, m2, m3, m4))
def setversion(k): _lib.lt_setversion(ctypes.byref(ctypes.c_int(k)))
