# -*- coding: utf-8 -*-

import dadi
from dadi import Numerics, PhiManip, Integration, Spectrum

# 1het model without t2 (isolation with migration)
def model_func(params, ns, pts):
    """
    Heterogeneous migration model with two populations P represents propotion of genome with no mig.
    8 parameters
    P1: proportion with no migration in first time period
    """
    nu_1, nu_2, t1, nu11, nu12, me1_12, me1_21, P1 = params

    _Nanc_size = 1.0  # This value can be used in splits with fractions
     
    xx = Numerics.default_grid(pts)

    # First time period: mix no migration at P1 with 1-P1 gene flow
    phiN1 = PhiManip.phi_1D(xx)
    phiN1 = PhiManip.phi_1D_to_2D(xx, phiN1)
    phiN1 = Integration.two_pops(phiN1, xx, T=t1, nu1=nu_1, nu2=nu_2, m12=0, m21=0)
    phiN1 = Integration.two_pops(phiN1, xx, T=t1, nu1=nu11, nu2=nu12, m12=0, m21=0)

    phiI1 = PhiManip.phi_1D(xx)
    phiI1 = PhiManip.phi_1D_to_2D(xx, phiI1)
    phiI1 = Integration.two_pops(phiI1, xx, T=t1, nu1=nu_1, nu2=nu_2, m12=0, m21=0)
    phiI1 = Integration.two_pops(phiI1, xx, T=t1, nu1=nu11, nu2=nu12, m12=me1_12, m21=me1_21)

    phi = P1 * phiN1 + (1 - P1) * phiI1

    fs = Spectrum.from_phi(phi, ns, [xx]*len(ns))
    return fs
