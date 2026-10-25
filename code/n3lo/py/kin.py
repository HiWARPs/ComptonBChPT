"""kin.py -- Compton kinematics in nucleon-mass units (M=1), as used by ampl.f90."""
import numpy as np
def from_lab(Elab, th_lab_deg, M=938.272):
    """lab photon energy (MeV), lab angle -> x=(s-M^2)/2, t (units M^2)"""
    w = Elab / M
    c = np.cos(np.radians(th_lab_deg))
    wp = w / (1 + w * (1 - c))
    return w, -2 * w * wp * (1 - c)
def env(x, t, mpi2):
    y = -x - t / 2
    return dict(x=x, t=t, y=y, mpi2=mpi2, ixs=1 / (2 * x), iys=1 / (2 * y), iPq=1 / (2 * x + t / 2), PP=4 - t)
def from_cm(W2, th_cm_deg):
    """s (units M^2), cm angle -> x, t"""
    s = W2; k = (s - 1) / (2 * np.sqrt(s))
    return (s - 1) / 2, -2 * k * k * (1 - np.cos(np.radians(th_cm_deg)))
