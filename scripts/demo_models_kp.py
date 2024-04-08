#!/usr/anaconda3/env/dadi211/bin/python
# -*- coding: utf-8 -*-

"""
@author: kprata
@date created: 5/5/21
@description: Custom models for demographic analysis using the site frequency spectrum.

Compatible with python 3.6.11 and dadi 2.1.1
"""
from dadi import Numerics, PhiManip, Integration
from dadi.Spectrum_mod import Spectrum
import numpy as np


# Models for testing one population scenarios.

def no_divergence_1d(notused, ns, pts):
    """
    Standard neutral model.

    ns = (n1,)

    n1: Number of samples in resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)

    fs = Spectrum.from_phi(phi, ns, (xx,))
    return fs


def instant_change(params, ns, pts):
    """
    Instantaneous size change some time ago.

    params = (nu,T)
    ns = (n1,)

    nu: Ratio of contemporary to ancient population size
    T: Time in the past at which size change happened (in units of 2*Na
       generations)
    n1: Number of samples in resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nu, T = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)

    phi = Integration.one_pop(phi, xx, T, nu)

    fs = Spectrum.from_phi(phi, ns, (xx,))
    return fs


def bottlegrowth(params, ns, pts):
    """
    Instantanous size change followed by exponential growth.

    params = (nuB,nuF,T)
    ns = (n1,)

    nuB: Ratio of population size after instantanous change to ancient
         population size
    nuF: Ratio of contemporary to ancient population size
    T: Time in the past at which instantaneous change happened and growth began
       (in units of 2*Na generations)
    n1: Number of samples in resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nuB, nuF, T = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)

    nu_func = lambda t: nuB * np.exp(np.log(nuF / nuB) * t / T)
    phi = Integration.one_pop(phi, xx, T, nu_func)

    fs = Spectrum.from_phi(phi, ns, (xx,))
    return fs


def bottleneck(params, ns, pts):

    nuB, nuF, TB, TF = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)

    phi = Integration.one_pop(phi, xx, TB, nuB)
    phi = Integration.one_pop(phi, xx, TF, nuF)

    fs = Spectrum.from_phi(phi, ns, (xx,))
    return fs

# Models for testing two population scenarios.

def no_divergence(notused, ns, pts):
    """
    Standard neutral model, populations never diverge.
    """

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    fs = Spectrum.from_phi(phi, ns, (xx, xx))
    return fs


def no_migration(params, ns, pts):
    """
    Split into two populations, no migration.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T: Time in the past of split (in units of 2*Na generations)
    """
    nu1, nu2, T = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T, nu1, nu2, m12=0, m21=0)

    fs = Spectrum.from_phi(phi, ns, (xx, xx))
    return fs


def sym_migration(params, ns, pts):
    """
    Split into two populations, with symmetric migration.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T: Time in the past of split (in units of 2*Na generations)
    m: Migration rate between populations (2*Na*m)
    """
    nu1, nu2, m, T = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T, nu1, nu2, m12=m, m21=m)

    fs = Spectrum.from_phi(phi, ns, (xx, xx))
    return fs


def asym_migration(params, ns, pts):
    """
    Split into two populations, with different migration rates.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T: Time in the past of split (in units of 2*Na generations)
    m12: Migration from pop 2 to pop 1 (2*Na*m12)
    m21: Migration from pop 1 to pop 2
    """
    nu1, nu2, m12, m21, T = params
    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T, nu1, nu2, m12=m12, m21=m21)
    fs = Spectrum.from_phi(phi, ns, (xx, xx))

    return fs


def iso_inbreeding(params, ns, pts):
    """
    Spilt into two populations, with inbreeding in each population.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    F1: Inbreeding in population 1.
    F2: Inbreeding in population 2.
    T: Time in past of split (in units of 2*Na generations)
    """
    nu1, nu2, F1, F2, T = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T, nu1, nu2)

    fs = Spectrum.from_phi_inbreeding(phi, ns, (xx, xx), (F1, F2), (2, 2))
    return fs


def mig_inbreeding(params, ns, pts):
    """
    Spilt into two populations, with inbreeding in each population.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    F1: Inbreeding in population 1.
    F2: Inbreeding in population 2.
    m: Migration rate between populations (2*Na*m)
    T: Time in past of split (in units of 2*Na generations)
    """
    nu1, nu2, F1, F2, m, T = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T, nu1, nu2, m12=m, m21=m)

    fs = Spectrum.from_phi_inbreeding(phi, ns, (xx, xx), (F1, F2), (2, 2))
    return fs


# Two epoch models

def mig_be_inbred(params, ns, pts):
    """
    Population size change after split with two periods of migration after size change.

    nu1: Size of population 1 after split.
    nu1a: Change of population size 1.
    nu2: Size of population 2 after split.
    nu2a: Change of poulation size 2.
    F1: Inbreeding in population 1.
    F2: Inbreeding in population 2.
    m1: Migration rate between populations (2*Na*m) before size change.
    m2: Migration rate between populations (2*Na*m) after size change.
    T: Time in past of split (in units of 2*Na generations)
    """
    nu1, nu2, nu1a, nu2a, F1, F2, m1, m2, T1, T2 = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, T1, nu1=nu1, nu2=nu2, m12=m1, m21=m1)

    phi = Integration.two_pops(phi, xx, T2, nu1=nu1a, nu2=nu2a, m12=m2, m21=m2)
    fs = Spectrum.from_phi_inbreeding(phi, ns, (xx, xx), (F1, F2), (2, 2))
    return fs


def anc_sym_migration(params, ns, pts):
    """
    Split with symmetric migration followed by isolation.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m: Migration between pop 2 and pop 1.
    T1: The scaled time between the split and the ancient migration (in units of 2*Na generations).
    T2: The scaled time between the ancient migration and present.
    """
    nu1, nu2, m, T1, T2 = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1, nu1, nu2, m12=m, m21=m)

    phi = Integration.two_pops(phi, xx, T2, nu1, nu2, m12=0, m21=0)

    fs = Spectrum.from_phi(phi, ns, (xx, xx))
    return fs


def anc_sym_mig_inbred(params, ns, pts):
    """
    Split with symmetric migration followed by isolation.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m: Migration between pop 2 and pop 1.
    T1: The scaled time between the split and the ancient migration (in units of 2*Na generations).
    T2: The scaled time between the ancient migration and present.
    """
    nu1, nu2, F1, F2, m, T1, T2 = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1, nu1, nu2, m12=m, m21=m)

    phi = Integration.two_pops(phi, xx, T2, nu1, nu2, m12=0, m21=0)

    fs = Spectrum.from_phi_inbreeding(phi, ns, (xx, xx), (F1, F2), (2, 2))
    return fs


def anc_asym_migration(params, ns, pts):
    """
    Split with asymmetric migration followed by isolation.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    T1: The scaled time between the split and the ancient migration (in units of 2*Na generations).
    T2: The scaled time between the ancient migration and present.
    """
    nu1, nu2, m12, m21, T1, T2 = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1, nu1, nu2, m12=m12, m21=m21)

    phi = Integration.two_pops(phi, xx, T2, nu1, nu2, m12=0, m21=0)

    fs = Spectrum.from_phi(phi, ns, (xx, xx))
    return fs


def sec_contact_sym_migration(params, ns, pts):
    """
    Split with no gene flow, followed by period of symmetrical gene flow.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m: Migration between pop 2 and pop 1.
    T1: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    T2: The scaled time between the secondary contact and present.
    """
    nu1, nu2, m, T1, T2 = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1, nu1, nu2, m12=0, m21=0)

    phi = Integration.two_pops(phi, xx, T2, nu1, nu2, m12=m, m21=m)

    fs = Spectrum.from_phi(phi, ns, (xx, xx))
    return fs


def sec_contact_sym_mig_inbred(params, ns, pts):
    """
    Split with no gene flow, followed by period of symmetrical gene flow.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m: Migration between pop 2 and pop 1.
    T1: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    T2: The scaled time between the secondary contact and present.
    """
    nu1, nu2, F1, F2, m, T1, T2 = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1, nu1, nu2, m12=0, m21=0)

    phi = Integration.two_pops(phi, xx, T2, nu1, nu2, m12=m, m21=m)

    fs = Spectrum.from_phi_inbreeding(phi, ns, (xx, xx), (F1, F2), (2, 2))
    return fs


def sec_contact_asym_migration(params, ns, pts):
    """
    Split with no gene flow, followed by period of asymmetrical gene flow.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    T1: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    T2: The scaled time between the secondary contact and present.
    """
    nu1, nu2, m12, m21, T1, T2 = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1, nu1, nu2, m12=0, m21=0)

    phi = Integration.two_pops(phi, xx, T2, nu1, nu2, m12=m12, m21=m21)

    fs = Spectrum.from_phi(phi, ns, (xx, xx))
    return fs


# Three populations models.
# Note: all three population models include inbreeding


def split_nomig(params, ns, pts):
    """
    Model with split between pop 1 and (2,3), then split between 2 and 3.
    Migration does not occur between any population pair.
    nu1: Size of population 1 after split.
    nuA: Size of population (2,3) after split from 1.
    nu2: Size of population 2 after split.
    nu3: Size of population 3 after split.
    T1: The scaled time between the split of pops 1 vs 2 and 3 (in units of 2*Na generations).
    T2: The scaled time between the split of pops 2 and 3 (in units of 2*Na generations).
    """
    # 9 parameters
    nu1, nuA, nu2, nu3, F1, F2, F3, T1, T2 = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1, nu1=nu1, nu2=nuA, m12=0, m21=0)

    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)

    phi = Integration.three_pops(phi, xx, T2, nu1=nu1, nu2=nu2, nu3=nu3, m12=0, m21=0, m23=0, m32=0, m13=0, m31=0)

    fs = Spectrum.from_phi_inbreeding(phi, ns, (xx, xx, xx), (F1, F2, F3), (2, 2, 2))
    return fs


def split_symmig_all(params, ns, pts):
    """
    Model with split between pop 1 and (2,3), then split between 2 and 3.
    Migration is symmetrical between all population pairs (ie 1<->2, 2<->3, and 1<->3).
    nu1: Size of population 1 after split.
    nuA: Size of population (2,3) after split from 1.
    nu2: Size of population 2 after split.
    nu3: Size of population 3 after split.
    mA: Migration rate between population 1 and population (2,3)
    m1: Migration rate between populations 1 and 2 (2*Na*m)
    m2: Migration rate between populations 2 and 3
    m3: Migration rate between populations 1 and 3
    T1: The scaled time between the split of pops 1 vs 2 and 3 (in units of 2*Na generations).
    T2: The scaled time between the split of pops 2 and 3 (in units of 2*Na generations).
    """
    # 13 parameters
    nu1, nuA, nu2, nu3, F1, F2, F3, mA, m1, m2, m3, T1, T2 = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1, nu1=nu1, nu2=nuA, m12=mA, m21=mA)

    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)

    phi = Integration.three_pops(phi, xx, T2, nu1=nu1, nu2=nu2, nu3=nu3, m12=m1, m21=m1, m23=m2, m32=m2, m13=m3, m31=m3)

    fs = Spectrum.from_phi_inbreeding(phi, ns, (xx, xx, xx), (F1, F2, F3), (2, 2, 2))
    return fs


def split_symmig_adjacent(params, ns, pts):
    """
    Model with split between pop 1 and (2,3), then split between 2 and 3. Assume 2 occurs
    in between populations 1 and 3, which do not come in to contact with one another.
    Migration is symmetrical between 'adjacent' population pairs (ie 1<->2, 2<->3, but not 1<->3).
    nu1: Size of population 1 after split.
    nuA: Size of population (2,3) after split from 1.
    nu2: Size of population 2 after split.
    nu3: Size of population 3 after split.
    mA: Migration rate between population 1 and population (2,3)
    m1: Migration rate between populations 1 and 2 (2*Na*m)
    m2: Migration rate between populations 2 and 3
    m3: Migration rate between populations 1 and 3
    T1: The scaled time between the split of pops 1 vs 2 and 3 (in units of 2*Na generations).
    T2: The scaled time between the split of pops 2 and 3 (in units of 2*Na generations).
   """
    # 12 parameters
    nu1, nuA, nu2, nu3, mA, F1, F2, F3, m1, m2, T1, T2 = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1, nu1=nu1, nu2=nuA, m12=mA, m21=mA)

    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)

    phi = Integration.three_pops(phi, xx, T2, nu1=nu1, nu2=nu2, nu3=nu3, m12=m1, m21=m1, m23=m2, m32=m2, m13=0, m31=0)

    fs = Spectrum.from_phi_inbreeding(phi, ns, (xx, xx, xx), (F1, F2, F3), (2, 2, 2))
    return fs


def sec_cont_mig_1(params, ns, pts):
    """
    Model with split between pop 1 and (2,3), gene flow does not occur. Split between pops
    2 and 3, with gene flow. After appearance of 2 and 3, gene flow also occurs between 1
    and 2, and 1 and 3.
    'shorter isolation'
    nu1: Size of population 1 after split.
    nuA: Size of population (2,3) after split from 1.
    nu2: Size of population 2 after split.
    nu3: Size of population 3 after split.
    m1: Migration rate between populations 1 and 2 (2*Na*m)
    m2: Migration rate between populations 2 and 3
    m3: Migration rate between populations 1 and 3
    T1: The scaled time between the split of pops 1 vs 2 and 3 (in units of 2*Na generations).
    T2: The scaled time between the first split and the split of pops 2 and 3 (in units of 2*Na generations).
    """
    # 12 parameters
    nu1, nuA, nu2, nu3, F1, F2, F3, m1, m2, m3, T1, T2 = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1, nu1=nu1, nu2=nuA, m12=0, m21=0)

    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)

    phi = Integration.three_pops(phi, xx, T2, nu1=nu1, nu2=nu2, nu3=nu3, m12=m1, m21=m1, m23=m2, m32=m2, m13=m3, m31=m3)

    fs = Spectrum.from_phi_inbreeding(phi, ns, (xx, xx, xx), (F1, F2, F3), (2, 2, 2))
    return fs


def sec_cont_mig_2(params, ns, pts):
    """
    Model with split between pop 1 and (2,3), gene flow does not occur. Split between pops
    2 and 3, without gene flow between them. After a period of isolation between 1, 2 and 3,
    gene flow resumes between all populations.
    nu1: Size of population 1 after split.
    nuA: Size of population (2,3) after split from 1.
    nu2: Size of population 2 after split.
    nu3: Size of population 3 after split.
    m1: Migration rate between populations 1 and 2 (2*Na*m)
    m2: Migration rate between populations 2 and 3
    m3: Migration rate between populations 1 and 3
    T1: The scaled time between the split of pops 1 vs 2 and 3 (in units of 2*Na generations).
    T2: The scaled time between the first split and the split of pops 2 and 3 (in units of 2*Na generations).
    T3: The scaled time between the T2 and present.
    """
    # 13 parameters
    nu1, nuA, nu2, nu3, F1, F2, F3, m1, m2, m3, T1, T2, T3 = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1, nu1=nu1, nu2=nuA, m12=0, m21=0)

    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)

    phi = Integration.three_pops(phi, xx, T2, nu1=nu1, nu2=nu2, nu3=nu3, m12=0, m21=0, m23=0, m32=0, m13=0, m31=0)

    phi = Integration.three_pops(phi, xx, T3, nu1=nu1, nu2=nu2, nu3=nu3, m12=m1, m21=m1, m23=m2, m32=m2, m13=m2, m31=m2)

    fs = Spectrum.from_phi_inbreeding(phi, ns, (xx, xx, xx), (F1, F2, F3), (2, 2, 2))
    return fs


def sec_cont_mig_3(params, ns, pts):
    """
    Model with split between pop 1 and (2,3), gene flow does not occur. Then gene flow occurs.
    Split between pops 2 and 3, with gene flow, while gene flow also occurs between 1, 2, and 3.

    nu1: Size of population 1 after split.
    nuA: Size of population (2,3) after split from 1.
    nu2: Size of population 2 after split.
    nu3: Size of population 3 after split.
    m1: Migration rate between populations 1 and 2 (2*Na*m)
    m2: Migration rate between populations 2 and 3
    m3: Migration rate between populations 1 and 3
    T1a: The scaled time between the split of pops 1 vs 2 and 3 (in units of 2*Na generations).
    T1b: The scaled time between the first spilt and the spilt of pops 2 and 3.
    T2: The scaled time between T1b (spilt) and present.
    """
    # 14 parameters
    nu1, nuA, nu2, nu3, F1, F2, F3, mA, m1, m2, m3, T1a, T1b, T2 = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1a, nu1=nu1, nu2=nuA, m12=0, m21=0)

    phi = Integration.two_pops(phi, xx, T1a, nu1=nu1, nu2=nuA, m12=mA, m21=mA)

    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)

    phi = Integration.three_pops(phi, xx, T2, nu1=nu1, nu2=nu2, nu3=nu3, m12=m1, m21=m1, m23=m2, m32=m2, m13=m3, m31=m3)

    fs = Spectrum.from_phi_inbreeding(phi, ns, (xx, xx, xx), (F1, F2, F3), (2, 2, 2))
    return fs


def mig_sec_cont23(params, ns, pts):
    """
     Model with split between pop 1 and (2,3), with gene flow. Split
     between pops 2 and 3 without gene flow (while gene flow continues between 1 and 2, and 1 and 3.
     Then all gene flow resumes.
     nu1: Size of population 1 after split.
     nuA: Size of population (2,3) after split from 1.
     nu2: Size of population 2 after split.
     nu3: Size of population 3 after split.
     mA: Migration rate between population 1 and population (2,3)
     m1: Migration rate between populations 1 and 2 (2*Na*m)
     m2: Migration rate between populations 2 and 3
     m3: Migration rate between populations 1 and 3
     T1: The scaled time between the split of pops 1 vs 2 and 3 (in units of 2*Na generations).
     T2: The scaled time between the first split and the split of pops 2 and 3 (in units of 2*Na generations).
     T3: The scaled time between T2 and present (in units of 2*Na generations).
     """
    # 13 parameters
    nu1, nuA, nu2, nu3, F1, F2, F3, mA, m1, m2, T1, T2, T3 = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1, nu1=nu1, nu2=nuA, m12=mA, m21=mA)

    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)

    phi = Integration.three_pops(phi, xx, T2, nu1=nu1, nu2=nu2, nu3=nu3, m12=0, m21=0, m23=m2, m32=m2, m13=m2, m31=m2)

    phi = Integration.three_pops(phi, xx, T3, nu1=nu1, nu2=nu2, nu3=nu3, m12=m1, m21=m1, m23=m2, m32=m2, m13=m2, m31=m2)

    fs = Spectrum.from_phi_inbreeding(phi, ns, (xx, xx, xx), (F1, F2, F3), (2, 2, 2))
    return fs
