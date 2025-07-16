# -*- coding: utf-8 -*-

import dadi
from dadi import Numerics, PhiManip, Integration, Spectrum

def model_func(params, ns, pts):
    """
    Heterogeneous model with two populations and migration.
    18 parameters
    P1: proportion with m in first time period (1 - P1) for me
    P2: proportion with m migration in second time period (1 - P2) for me
    """
    nu_1, nu_2, t1, nu11, nu12, m1_12, me1_12, m1_21, me1_21, t2, nu21, nu22, m2_12, me2_12, m2_21, me2_21, P1, P2 = params

    _Nanc_size = 1.0  # This value can be used in splits with fractions
     
    xx = Numerics.default_grid(pts)

    # First time period: mix no migration at P1 with 1-P1 gene flow
    phiN1 = PhiManip.phi_1D(xx)
    phiN1 = PhiManip.phi_1D_to_2D(xx, phiN1)
    phiN1 = Integration.two_pops(phiN1, xx, T=t1, nu1=nu_1, nu2=nu_2, m12=0, m21=0)
    phiN1 = Integration.two_pops(phiN1, xx, T=t1, nu1=nu11, nu2=nu12, m12=m1_12, m21=m1_21)

    phiI1 = PhiManip.phi_1D(xx)
    phiI1 = PhiManip.phi_1D_to_2D(xx, phiI1)
    phiI1 = Integration.two_pops(phiI1, xx, T=t1, nu1=nu_1, nu2=nu_2, m12=0, m21=0)
    phiI1 = Integration.two_pops(phiI1, xx, T=t1, nu1=nu11, nu2=nu12, m12=me1_12, m21=me1_21)

    phi1 = P1 * phiN1 + (1 - P1) * phiI1

    # Second time period: mix no migration at P2 with 1-P2 gene flow
    phiN2 = Integration.two_pops(phi1, xx, T=t2, nu1=nu21, nu2=nu22, m12=m2_12, m21=m2_21)
    phiI2 = Integration.two_pops(phi1, xx, T=t2, nu1=nu21, nu2=nu22, m12=me2_12, m21=me2_21)

    phi = P2 * phiN2 + (1 - P2) * phiI2

    fs = Spectrum.from_phi(phi, ns, [xx]*len(ns))
    return fs
