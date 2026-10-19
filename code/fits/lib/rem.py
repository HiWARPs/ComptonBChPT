"""E2/M1 ratio of the gamma N Delta transition from the couplings g_M, g_E of the chiral Lagrangian (exact at Q^2 = 0):
R_EM = -G_E*/G_M*, with the Jones-Scadron form factors G_E* = -K g_E, G_M* = g_M + G_E*, K = (M_D - M_N)/(2 (M_D + M_N))
(Pascalutsa and Vanderhaeghen, PRD 73, 034003 (2006), Eqs. for G*, R_EM; author decision 2026-10-20: exact relation everywhere,
replacing the leading-order form K g_E/g_M).  R in %."""
MD, MN, GM0 = 1232.0, 938.272, 2.97
K = (MD - MN) / (2 * (MD + MN))


def rem(gE, gM=GM0):
    """R_EM [%]"""
    return 100 * K * gE / (gM - K * gE)


def ge_of_rem(R, gM=GM0):
    """g_E for R_EM = R [%] (inverse of rem)"""
    r = R / 100
    return r * gM / (K * (1 + r))


def drem_dge(gE, gM=GM0):
    """dR_EM/dg_E [% per unit g_E]"""
    return 100 * K * gM / (gM - K * gE) ** 2


# checks: PV2006 (g_M = 2.9, g_E = -1.0 -> -2.3 %); Phys. Rep. 437 (2007) (2.97, -1.0 -> G_M* = 3.04, G_E* = 0.07, -2.2 %)
assert round(rem(-1.0, 2.9), 1) == -2.3 and round(rem(-1.0), 1) == -2.2
assert round(GM0 + K, 2) == 3.04 and round(K, 2) == 0.07
assert abs(ge_of_rem(rem(-0.77)) + 0.77) < 1e-12 and abs((rem(-1.0 + 1e-6) - rem(-1.0 - 1e-6)) / 2e-6 - drem_dge(-1.0)) < 1e-6
