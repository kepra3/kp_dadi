# -*- coding: utf-8 -*-
"""moments port of custom_model_het_sc.py (11 params)."""
import numpy as np
import moments

def model_func(params, ns):
    nu_1, nu_2, t1, nu11, nu12, t2, nu21, nu22, me2_12, me2_21, P2 = params
    fs = moments.Demographics1D.snm([ns[0] + ns[1]]).split(0, ns[0], ns[1])
    fs.integrate([nu_1, nu_2], t1, m=np.zeros((2, 2)))
    fs.integrate([nu11, nu12], t1, m=np.zeros((2, 2)))
    fsN = fs.copy(); fsN.integrate([nu21, nu22], t2, m=np.zeros((2, 2)))
    fsI = fs.copy(); fsI.integrate([nu21, nu22], t2, m=np.array([[0, me2_12], [me2_21, 0]]))
    return P2 * fsN + (1 - P2) * fsI
