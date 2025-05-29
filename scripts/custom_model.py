# -*- coding: utf-8 -*-

import dadi
from dadi import Numerics, PhiManip, Integration, Spectrum

def model_func(params, ns, pts):
    """
    Heterogeneous model with two populations and migration.
    17 parameters
    """
    nu_1, nu_2, t1, nu11, nu12, m1_12, me1_12, m1_21, me1_21, t2, nu21, nu22, m2_12, m2_21, me2_12, me2_21, P   = params

    _Nanc_size = 1.0  # This value can be used in splits with fractions

    xx = Numerics.default_grid(pts)

    phiN = PhiManip.phi_1D(xx)
    phiN = PhiManip.phi_1D_to_2D(xx, phiN)

    phiN = Integration.two_pops(phiN, xx, T=t1, nu1=nu_1, nu2=nu_2, m12=0, m21=0)
    phiN = Integration.two_pops(phiN, xx, T=t1, nu1=nu11, nu2=nu12, m12=m1_12, m21=m1_21)
    phiN = Integration.two_pops(phiN, xx, T=t2, nu1=nu21, nu2=nu22, m12=m2_12, m21=m2_21)

    fsN = Spectrum.from_phi(phiN, ns, [xx]*len(ns))
    
    phiI = PhiManip.phi_1D(xx)
    phiI = PhiManip.phi_1D_to_2D(xx, phiI)

    phiI = Integration.two_pops(phiI, xx, T=t1, nu1=nu_1, nu2=nu_2, m12=0, m21=0)

    phiI = Integration.two_pops(phiI, xx, T=t1, nu1=nu11, nu2=nu12, m12=me1_12, m21=me1_21)

    phiI = Integration.two_pops(phiI, xx, T=t2, nu1=nu21, nu2=nu22, m12=me2_12, m21=me2_21)
    
    fsI = Spectrum.from_phi(phiI, ns, [xx]*len(ns))


    fs = P * fsN + (1 - P) * fsI
    return fs