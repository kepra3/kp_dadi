# -*- coding: utf-8 -*-
"""moments port of custom_model_het_sym.py (8 params). Mirrors the dadi version exactly,
including the two consecutive epochs that both last t1."""
import numpy as np
import moments

def _mig(m12, m21):
    # moments m[i, j] = migration into i from j (same convention as dadi m12/m21)
    return np.array([[0, m12], [m21, 0]])

def _split(ns):
    fs = moments.Demographics1D.snm([ns[0] + ns[1]])
    return fs.split(0, ns[0], ns[1])

def model_func(params, ns):
    nu_1, nu_2, t1, nu11, nu12, me1_12, me1_21, P1 = params
    fsN = _split(ns)
    fsN.integrate([nu_1, nu_2], t1, m=_mig(0, 0))
    fsN.integrate([nu11, nu12], t1, m=_mig(0, 0))
    fsI = _split(ns)
    fsI.integrate([nu_1, nu_2], t1, m=_mig(0, 0))
    fsI.integrate([nu11, nu12], t1, m=_mig(me1_12, me1_21))
    return P1 * fsN + (1 - P1) * fsI
