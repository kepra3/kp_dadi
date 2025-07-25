# -*- coding: utf-8 -*-

import dadi
from dadi import Numerics, PhiManip, Integration, Spectrum

def model_func(params, ns, pts):
    """
    Heterogeneous model with two populations and secondary contact.
    11 parameters: ... as before ..., P2
    P2: proportion with no migration in second time period
    """
    nu_1, nu_2, t1, nu11, nu12, t2, nu21, nu22, me2_12, me2_21, P2 = params

    _Nanc_size = 1.0  # This value can be used in splits with fractions
     
    xx = Numerics.default_grid(pts)

    # First time period
    phi1 = PhiManip.phi_1D(xx)
    phi1 = PhiManip.phi_1D_to_2D(xx, phiI1)
    phi1 = Integration.two_pops(phiI1, xx, T=t1, nu1=nu_1, nu2=nu_2, m12=0, m21=0)
    phi1 = Integration.two_pops(phiI1, xx, T=t1, nu1=nu11, nu2=nu12, m12=0, m21=0)


    # Second time period: mix no migration at P2 with 1-P2 gene flow
    phiN2 = Integration.two_pops(phi1, xx, T=t2, nu1=nu21, nu2=nu22, m12=0, m21=0)
    phiI2 = Integration.two_pops(phi1, xx, T=t2, nu1=nu21, nu2=nu22, m12=me2_12, m21=me2_21)

    phi = P2 * phiN2 + (1 - P2) * phiI2

    fs = Spectrum.from_phi(phi, ns, [xx]*len(ns))
    return fs
