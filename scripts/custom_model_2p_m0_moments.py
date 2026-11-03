# -*- coding: utf-8 -*-
"""moments port of custom_model_2p_m0.py (14 params; the 'het'/'1het' model)."""
import numpy as np
import moments

def model_func(params, ns):
    (nu_1, nu_2, t1, nu11, nu12, me1_12, me1_21,
     t2, nu21, nu22, me2_12, me2_21, P1, P2) = params
    z = np.zeros((2, 2))
    def base():
        fs = moments.Demographics1D.snm([ns[0] + ns[1]]).split(0, ns[0], ns[1])
        fs.integrate([nu_1, nu_2], t1, m=z)
        return fs
    fsN1 = base(); fsN1.integrate([nu11, nu12], t1, m=z)
    fsI1 = base(); fsI1.integrate([nu11, nu12], t1, m=np.array([[0, me1_12], [me1_21, 0]]))
    fs1 = P1 * fsN1 + (1 - P1) * fsI1
    fsN2 = fs1.copy(); fsN2.integrate([nu21, nu22], t2, m=z)
    fsI2 = fs1.copy(); fsI2.integrate([nu21, nu22], t2, m=np.array([[0, me2_12], [me2_21, 0]]))
    return P2 * fsN2 + (1 - P2) * fsI2
