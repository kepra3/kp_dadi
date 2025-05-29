#!/usr/anaconda3/env/dadi211/bin/python
# -*- coding: utf-8 -*-

"""
@author: kprata
@date created: 5/5/21
@description: Custom models for demographic analysis using the site frequency spectrum.

Compatible with python 3.6.11 and 2.1.1
"""
from dadi import Numerics, PhiManip, Integration, Spectrum
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


def split_bottlegrowth(params, ns, pts):
    """
    Split into two populations

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T1: Time in the past of split (in units of 2*Na generations) and before growth (T2)

    Instantanous size change followed by exponential growth for both pops.

    nu1, nu2 : Ratio of population size after instantanous change to ancient
         population size
    nu1F, nu2F : Ratio of contemporary to ancient population size
    T2: Time in the past at which growth began
       (in units of 2*Na generations)
    ns: Number of samples in resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nu1, nu2, nu1F, nu2F, T1, T2 = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    nu_func1 = lambda T1: nu1 * np.exp(np.log(nu1F / nu1) * T1 / T2)
    nu_func2 = lambda T1: nu2 * np.exp(np.log(nu2F / nu2) * T1 / T2)

    phi = Integration.two_pops(phi, xx, T1, nu1=nu1, nu2=nu2, m12=0, m21=0)
    phi = Integration.two_pops(phi, xx, T2, nu1=nu_func1, nu2=nu_func2, m12=0, m21=0)

    fs = Spectrum.from_phi(phi, ns, (xx, xx))
    return fs


def split_bottlegrowth_asym_mig(params, ns, pts):
    """
    Split into two populations

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T1: Time in the past of split (in units of 2*Na generations) and before growth (T2)

    Instantanous size change followed by exponential growth for both pops.

    nu1, nu2 : Ratio of population size after instantanous change to ancient
         population size
    nu1F, nu2F : Ratio of contemporary to ancient population size
    T2: Time in the past at which growth began
       (in units of 2*Na generations)
    m12T1: migration rate from pop1 into pop2 over T1
    m21T1: migration rate from pop2 into pop1 over T1
    m12T2: migration rate from pop1 into pop2 over T2
    m21T2: migration rate from pop2 into pop1 over T2
    ns: Number of samples in resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nu1, nu2, nu1F, nu2F, T1, T2, m12T1, m21T1, m12T2, m21T2 = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    nu_func1 = lambda T1: nu1 * np.exp(np.log(nu1F / nu1) * T1 / T2)
    nu_func2 = lambda T1: nu2 * np.exp(np.log(nu2F / nu2) * T1 / T2)

    phi = Integration.two_pops(phi, xx, T1, nu1=nu1, nu2=nu2, m12=m12T1, m21=m21T1)
    phi = Integration.two_pops(phi, xx, T2, nu1=nu_func1, nu2=nu_func2, m12=m12T2, m21=m21T2)

    fs = Spectrum.from_phi(phi, ns, (xx, xx))
    return fs


def split_bottlegrowth_ancient_asym_mig(params, ns, pts):
    """
    Split into two populations

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T1: Time in the past of split (in units of 2*Na generations) and before growth (T2)

    Instantanous size change followed by exponential growth for both pops.

    nu1, nu2 : Ratio of population size after instantanous change to ancient
         population size
    nu1F, nu2F : Ratio of contemporary to ancient population size
    T2: Time in the past at which growth began
       (in units of 2*Na generations)
    m12: migration rate from pop1 into pop2 over T2
    m21: migration rate from pop2 into pop1 over T2
    ns: Number of samples in resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nu1, nu2, nu1F, nu2F, T1, T2, m12, m21 = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    nu_func1 = lambda T1: nu1 * np.exp(np.log(nu1F / nu1) * T1 / T2)
    nu_func2 = lambda T1: nu2 * np.exp(np.log(nu2F / nu2) * T1 / T2)

    phi = Integration.two_pops(phi, xx, T1, nu1=nu1, nu2=nu2, m12=m12, m21=m21)
    phi = Integration.two_pops(phi, xx, T2, nu1=nu_func1, nu2=nu_func2, m12=0, m21=0)

    fs = Spectrum.from_phi(phi, ns, (xx, xx))
    return fs


def split_bottlegrowth_second_asym_mig(params, ns, pts):
    """
    Split into two populations

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T1: Time in the past of split (in units of 2*Na generations) and before growth (T2)

    Instantanous size change followed by exponential growth for both pops.

    nu1, nu2 : Ratio of population size after instantanous change to ancient
         population size
    nu1F, nu2F : Ratio of contemporary to ancient population size
    T2: Time in the past at which growth began
       (in units of 2*Na generations)
    m12: migration rate from pop1 into pop2 over T2
    m21: migration rate from pop2 into pop1 over T2
    ns: Number of samples in resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nu1, nu2, nu1F, nu2F, T1, T2, m12, m21 = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    nu_func1 = lambda T1: nu1 * np.exp(np.log(nu1F / nu1) * T1 / T2)
    nu_func2 = lambda T1: nu2 * np.exp(np.log(nu2F / nu2) * T1 / T2)

    phi = Integration.two_pops(phi, xx, T1, nu1=nu1, nu2=nu2, m12=0, m21=0)
    phi = Integration.two_pops(phi, xx, T2, nu1=nu_func1, nu2=nu_func2, m12=m12, m21=m21)

    fs = Spectrum.from_phi(phi, ns, (xx, xx))
    return fs


def split_sizechange(params, ns, pts):
    """
    Split into two populations

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T1: Time in the past of split (in units of 2*Na generations) and before growth (T2)

    Instantanous size change followed by size change both pops.

    nu1T1, nu2T1 : Ratio of population size after instantanous change to ancient
         population size
    nu1T2, nu2T2 : Ratio of contemporary to ancient population size
    T2: Time in the past at which growth began
       (in units of 2*Na generations)
    ns: Number of samples in resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nu1T1, nu2T1, nu1T2, nu2T2, T1, T2 = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1, nu1=nu1T1, nu2=nu2T1, m12=0, m21=0)
    phi = Integration.two_pops(phi, xx, T2, nu1=nu1T2, nu2=nu2T2, m12=0, m21=0)

    fs = Spectrum.from_phi(phi, ns, (xx, xx))
    return fs


def split_sizechange_asym_mig(params, ns, pts):
    """
    Split into two populations

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T1: Time in the past of split (in units of 2*Na generations) and before growth (T2)

    Instantanous size change followed by size change both pops.

    nu1T1, nu2T1 : Ratio of population size after instantanous change to ancient
         population size
    nu1T2, nu2T2 : Ratio of contemporary to ancient population size
    T2: Time in the past at which growth began
       (in units of 2*Na generations)
    m12T1: migration rate from pop1 into pop2 over T1
    m21T1: migration rate from pop2 into pop1 over T1
    m12T2: migration rate from pop1 into pop2 over T2
    m21T2: migration rate from pop2 into pop1 over T2
    ns: Number of samples in resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nu1T1, nu2T1, nu1T2, nu2T2, T1, T2, m12T1, m21T1, m12T2, m21T2 = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1, nu1=nu1T1, nu2=nu2T1, m12=m12T1, m21=m21T1)
    phi = Integration.two_pops(phi, xx, T2, nu1=nu1T2, nu2=nu2T2, m12=m12T2, m21=m21T2)

    fs = Spectrum.from_phi(phi, ns, (xx, xx))
    return fs


def split_sizechange_ancient_asym_mig(params, ns, pts):
    """
    Split into two populations

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T1: Time in the past of split (in units of 2*Na generations) and before growth (T2)

    Instantanous size change followed by size change both pops.

    nu1T1, nu2T1 : Ratio of population size after instantanous change to ancient
         population size
    nu1T2, nu2T2 : Ratio of contemporary to ancient population size
    T2: Time in the past at which growth began
       (in units of 2*Na generations)
    m12: migration rate from pop1 into pop2 over T1
    m21: migration rate from pop2 into pop1 over T1
    ns: Number of samples in resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nu1T1, nu2T1, nu1T2, nu2T2, T1, T2, m12, m21 = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1, nu1=nu1T1, nu2=nu2T1, m12=m12, m21=m21)
    phi = Integration.two_pops(phi, xx, T2, nu1=nu1T2, nu2=nu2T2, m12=0, m21=0)

    fs = Spectrum.from_phi(phi, ns, (xx, xx))
    return fs


def split_sizechange_second_asym_mig(params, ns, pts):
    """
    Split into two populations

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T1: Time in the past of split (in units of 2*Na generations) and before growth (T2)

    Instantanous size change followed by size change both pops.

    nu1T1, nu2T1 : Ratio of population size after instantanous change to ancient
         population size
    nu1T2, nu2T2 : Ratio of contemporary to ancient population size
    T2: Time in the past at which growth began
       (in units of 2*Na generations)
    m12: migration rate from pop1 into pop2 over T2
    m21: migration rate from pop2 into pop1 over T2
    ns: Number of samples in resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nu1T1, nu2T1, nu1T2, nu2T2, T1, T2, m12, m21 = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1, nu1=nu1T1, nu2=nu2T1, m12=0, m21=0)
    phi = Integration.two_pops(phi, xx, T2, nu1=nu1T2, nu2=nu2T2, m12=m12, m21=m21)

    fs = Spectrum.from_phi(phi, ns, (xx, xx))
    return fs


def hetero_asym_migration(params, ns, pts):
    """
    Split into two populations, with different migration rates.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T: Time in the past of split (in units of 2*Na generations)
    m12: Migration from pop 2 to pop 1 (2*Na*m12)
    m21: Migration from pop 1 to pop 2
    me12: Reduced effective migration rate from pop 2 to pop 1
    me21: Reduced effective migration rate from pop 2 to pop 2
    P: The proportion of the genome evolving neutrally
    """
    nu1, nu2, m12, m21, me12, me21, T, P = params
    xx = Numerics.default_grid(pts)

    phiN = PhiManip.phi_1D(xx)
    phiN = PhiManip.phi_1D_to_2D(xx, phiN)

    phiN = Integration.two_pops(phiN, xx, T, nu1, nu2, m12=m12, m21=m21)
    fsN = Spectrum.from_phi(phiN, ns, (xx, xx))

    phiI = PhiManip.phi_1D(xx)
    phiI = PhiManip.phi_1D_to_2D(xx, phiI)

    phiI = Integration.two_pops(phiI, xx, T, nu1, nu2, m12=me12, m21=me21)
    fsI = Spectrum.from_phi(phiI, ns, (xx, xx))

    fs = P * fsN + (1 - P) * fsI
    return fs


def anc_hetero_asym_migration(params, ns, pts):
    """
    Split with asymmetric migration followed by isolation.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    me12: Effective migration from pop 2 to pop 1 (2*Na*m12).
    me21: Effective migration from pop 1 to pop 2.
    T1: The scaled time between the split and the ancient migration (in units of 2*Na generations).
    T2: The scaled time between the ancient migration and present.
    P: The proportion of the genome evolving neutrally
    """
    nu1, nu2, m12, m21, me12, me21, T1, T2, P = params

    xx = Numerics.default_grid(pts)

    phiN = PhiManip.phi_1D(xx)
    phiN = PhiManip.phi_1D_to_2D(xx, phiN)

    phiN = Integration.two_pops(phiN, xx, T1, nu1, nu2, m12=m12, m21=m21)

    phiN = Integration.two_pops(phiN, xx, T2, nu1, nu2, m12=0, m21=0)

    fsN = Spectrum.from_phi(phiN, ns, (xx, xx))

    phiI = PhiManip.phi_1D(xx)
    phiI = PhiManip.phi_1D_to_2D(xx, phiI)

    phiI = Integration.two_pops(phiI, xx, T1, nu1, nu2, m12=me12, m21=me21)

    phiI = Integration.two_pops(phiI, xx, T2, nu1, nu2, m12=0, m21=0)

    fsI = Spectrum.from_phi(phiI, ns, (xx, xx))

    fs = P * fsN + (1 - P) * fsI
    return fs


def sec_contact_hetero_asym_migration(params, ns, pts):
    """
    Split with no gene flow, followed by period of asymmetrical gene flow.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    me12: Effective migration from pop 2 to pop 1 (2*Na*m12).
    me21: Effective migration from pop 1 to pop 2.
    T1: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    T2: The scaled time between the secondary contact and present.
    P: The proportion of the genome evolving neutrally
    """
    nu1, nu2, m12, m21, me12, me21, T1, T2, P = params

    xx = Numerics.default_grid(pts)

    phiN = PhiManip.phi_1D(xx)
    phiN = PhiManip.phi_1D_to_2D(xx, phiN)

    phiN = Integration.two_pops(phiN, xx, T1, nu1, nu2, m12=0, m21=0)

    phiN = Integration.two_pops(phiN, xx, T2, nu1, nu2, m12=m12, m21=m21)

    fsN = Spectrum.from_phi(phiN, ns, (xx, xx))

    phiI = PhiManip.phi_1D(xx)
    phiI = PhiManip.phi_1D_to_2D(xx, phiI)

    phiI = Integration.two_pops(phiI, xx, T1, nu1, nu2, m12=0, m21=0)

    phiI = Integration.two_pops(phiI, xx, T2, nu1, nu2, m12=me12, m21=me21)

    fsI = Spectrum.from_phi(phiI, ns, (xx, xx))

    fs = P * fsN + (1 - P) * fsI
    return fs


def split_bottlegrowth_hetero_asym_mig(params, ns, pts):
    """
    Split into two populations

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T1: Time in the past of split (in units of 2*Na generations) and before growth (T2)

    Instantanous size change followed by exponential growth for both pops.

    nu1, nu2 : Ratio of population size after instantanous change to ancient
         population size
    nu1F, nu2F : Ratio of contemporary to ancient population size
    T2: Time in the past at which growth began
       (in units of 2*Na generations)
    m12T1: migration rate from pop1 into pop2 over T1
    m21T1: migration rate from pop2 into pop1 over T1
    me12T1: migration rate from pop1 into pop2 over T1
    me21T1: migration rate from pop2 into pop1 over T1
    m12T2: migration rate from pop1 into pop2 over T2
    m21T2: migration rate from pop2 into pop1 over T2
    me12T2: migration rate from pop1 into pop2 over T2
    me21T2: migration rate from pop2 into pop1 over T2
    P: proportion of the genome evolving neutrally
    ns: Number of samples in resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nu1, nu2, nu1F, nu2F, T1, T2, m12T1, m21T1, me12T1, me21T1, m12T2, m21T2, me12T2, me21T2, P = params

    xx = Numerics.default_grid(pts)

    phiN = PhiManip.phi_1D(xx)
    phiN = PhiManip.phi_1D_to_2D(xx, phiN)

    nu_func1 = lambda T1: nu1 * np.exp(np.log(nu1F / nu1) * T1 / T2)
    nu_func2 = lambda T1: nu2 * np.exp(np.log(nu2F / nu2) * T1 / T2)

    phiN = Integration.two_pops(phiN, xx, T1, nu1=nu1, nu2=nu2, m12=m12T1, m21=m21T1)
    phiN = Integration.two_pops(phiN, xx, T2, nu1=nu_func1, nu2=nu_func2, m12=m12T2, m21=m21T2)

    fsN = Spectrum.from_phi(phiN, ns, (xx, xx))

    phiI = PhiManip.phi_1D(xx)
    phiI = PhiManip.phi_1D_to_2D(xx, phiI)

    phiI = Integration.two_pops(phiI, xx, T1, nu1=nu1, nu2=nu2, m12=me12T1, m21=me21T1)
    phiI = Integration.two_pops(phiI, xx, T2, nu1=nu_func1, nu2=nu_func2, m12=me12T2, m21=me21T2)

    fsI = Spectrum.from_phi(phiI, ns, (xx, xx))

    fs = P * fsN + (1 - P) * fsI
    return fs


def split_bottlegrowth_ancient_hetero_asym_mig(params, ns, pts):
    """
    Split into two populations

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T1: Time in the past of split (in units of 2*Na generations) and before growth (T2)

    Instantanous size change followed by exponential growth for both pops.

    nu1, nu2 : Ratio of population size after instantanous change to ancient
         population size
    nu1F, nu2F : Ratio of contemporary to ancient population size
    T2: Time in the past at which growth began
       (in units of 2*Na generations)
    m12: migration rate from pop1 into pop2 over T1
    m21: migration rate from pop2 into pop1 over T1
    me12: effective migration rate from pop1 into pop2 over T1
    me21: effective migration rate from pop2 into pop1 over T1
    P: The proportion of the genome evolving neutrally
    ns: Number of samples in resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nu1, nu2, nu1F, nu2F, T1, T2, m12, m21, me12, me21, P = params

    xx = Numerics.default_grid(pts)

    phiN = PhiManip.phi_1D(xx)
    phiN = PhiManip.phi_1D_to_2D(xx, phiN)

    nu_func1 = lambda T1: nu1 * np.exp(np.log(nu1F / nu1) * T1 / T2)
    nu_func2 = lambda T1: nu2 * np.exp(np.log(nu2F / nu2) * T1 / T2)

    phiN = Integration.two_pops(phiN, xx, T1, nu1=nu1, nu2=nu2, m12=m12, m21=m21)
    phiN = Integration.two_pops(phiN, xx, T2, nu1=nu_func1, nu2=nu_func2, m12=0, m21=0)

    fsN = Spectrum.from_phi(phiN, ns, (xx, xx))

    phiI = PhiManip.phi_1D(xx)
    phiI = PhiManip.phi_1D_to_2D(xx, phiI)

    phiI = Integration.two_pops(phiI, xx, T1, nu1=nu1, nu2=nu2, m12=me12, m21=me21)
    phiI = Integration.two_pops(phiI, xx, T2, nu1=nu_func1, nu2=nu_func2, m12=0, m21=0)

    fsI = Spectrum.from_phi(phiI, ns, (xx, xx))

    fs = P * fsN + (1 - P) * fsI
    return fs


def split_bottlegrowth_second_hetero_asym_mig(params, ns, pts):
    """
    Split into two populations

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T1: Time in the past of split (in units of 2*Na generations) and before growth (T2)

    Instantanous size change followed by exponential growth for both pops.

    nu1, nu2 : Ratio of population size after instantanous change to ancient
         population size
    nu1F, nu2F : Ratio of contemporary to ancient population size
    T2: Time in the past at which growth began
       (in units of 2*Na generations)
    m12: migration rate from pop1 into pop2 over T2
    m21: migration rate from pop2 into pop1 over T2
    me12: effective migration rate from pop1 into pop2 over T2
    me21: effective migration rate from pop2 into pop1 over T2
    P: The proportion of the genome evolving neutrally
    ns: Number of samples in resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nu1, nu2, nu1F, nu2F, T1, T2, m12, m21, me12, me21, P = params

    xx = Numerics.default_grid(pts)

    phiN = PhiManip.phi_1D(xx)
    phiN = PhiManip.phi_1D_to_2D(xx, phiN)

    nu_func1 = lambda T1: nu1 * np.exp(np.log(nu1F / nu1) * T1 / T2)
    nu_func2 = lambda T1: nu2 * np.exp(np.log(nu2F / nu2) * T1 / T2)

    phiN = Integration.two_pops(phiN, xx, T1, nu1=nu1, nu2=nu2, m12=0, m21=0)
    phiN = Integration.two_pops(phiN, xx, T2, nu1=nu_func1, nu2=nu_func2, m12=m12, m21=m21)

    fsN = Spectrum.from_phi(phiN, ns, (xx, xx))

    phiI = PhiManip.phi_1D(xx)
    phiI = PhiManip.phi_1D_to_2D(xx, phiI)

    phiI = Integration.two_pops(phiI, xx, T1, nu1=nu1, nu2=nu2, m12=0, m21=0)
    phiI = Integration.two_pops(phiI, xx, T2, nu1=nu_func1, nu2=nu_func2, m12=me12, m21=me21)

    fsI = Spectrum.from_phi(phiI, ns, (xx, xx))

    fs = P * fsN + (1 - P) * fsI
    return fs


def split_sizechange_hetero_asym_mig(params, ns, pts):
    """
    Split into two populations

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T1: Time in the past of split (in units of 2*Na generations) and before growth (T2)

    Instantanous size change followed by size change both pops.

    nu1T1, nu2T1 : Ratio of population size after instantanous change to ancient
         population size
    nu1T2, nu2T2 : Ratio of contemporary to ancient population size
    T2: Time in the past at which growth began
       (in units of 2*Na generations)
    m12T1: migration rate from pop1 into pop2 over T1
    m21T1: migration rate from pop2 into pop1 over T1
    me12T1: effective migration rate from pop1 into pop2 over T1
    me21T1: effective migration rate from pop2 into pop1 over T1
    m12T2: migration rate from pop1 into pop2 over T2
    m21T2: migration rate from pop2 into pop1 over T2
    me12T2: effective migration rate from pop1 into pop2 over T2
    me21T2: effective migration rate from pop2 into pop1 over T2
    P: The proportion of the genome evolving neutrally
    ns: Number of samples in resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nu1T1, nu2T1, nu1T2, nu2T2, T1, T2, m12T1, m21T1, me12T1, me21T1, m12T2, m21T2, me12T2, me21T2, P = params

    xx = Numerics.default_grid(pts)

    phiN = PhiManip.phi_1D(xx)
    phiN = PhiManip.phi_1D_to_2D(xx, phiN)

    phiN = Integration.two_pops(phiN, xx, T1, nu1=nu1T1, nu2=nu2T1, m12=m12T1, m21=m21T1)
    phiN = Integration.two_pops(phiN, xx, T2, nu1=nu1T2, nu2=nu2T2, m12=m12T2, m21=m21T2)

    fsN = Spectrum.from_phi(phiN, ns, (xx, xx))

    phiI = PhiManip.phi_1D(xx)
    phiI = PhiManip.phi_1D_to_2D(xx, phiI)

    phiI = Integration.two_pops(phiI, xx, T1, nu1=nu1T1, nu2=nu2T1, m12=me12T1, m21=me21T1)
    phiI = Integration.two_pops(phiI, xx, T2, nu1=nu1T2, nu2=nu2T2, m12=me12T2, m21=me21T2)

    fsI = Spectrum.from_phi(phiI, ns, (xx, xx))

    fs = P * fsN + (1 - P) * fsI
    return fs


def split_sizechange_ancient_hetero_asym_mig(params, ns, pts):
    """
    Split into two populations

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T1: Time in the past of split (in units of 2*Na generations) and before growth (T2)

    Instantanous size change followed by size change both pops.

    nu1T1, nu2T1 : Ratio of population size after instantanous change to ancient
         population size
    nu1T2, nu2T2 : Ratio of contemporary to ancient population size
    T2: Time in the past at which growth began
       (in units of 2*Na generations)
    m12: migration rate from pop1 into pop2 over T1
    m21: migration rate from pop2 into pop1 over T1
    me12: effective migration rate from pop1 into pop2 over T1
    me21: effective migration rate from pop2 into pop1 over T1
    P: The proportion of the genome evolving neutrally
    ns: Number of samples in resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nu1T1, nu2T1, nu1T2, nu2T2, T1, T2, m12, m21, me12, me21, P = params

    xx = Numerics.default_grid(pts)

    phiN = PhiManip.phi_1D(xx)
    phiN = PhiManip.phi_1D_to_2D(xx, phiN)

    phiN = Integration.two_pops(phiN, xx, T1, nu1=nu1T1, nu2=nu2T1, m12=m12, m21=m21)
    phiN = Integration.two_pops(phiN, xx, T2, nu1=nu1T2, nu2=nu2T2, m12=0, m21=0)

    fsN = Spectrum.from_phi(phiN, ns, (xx, xx))

    phiI = PhiManip.phi_1D(xx)
    phiI = PhiManip.phi_1D_to_2D(xx, phiI)

    phiI = Integration.two_pops(phiI, xx, T1, nu1=nu1T1, nu2=nu2T1, m12=me12, m21=me21)
    phiI = Integration.two_pops(phiI, xx, T2, nu1=nu1T2, nu2=nu2T2, m12=0, m21=0)

    fsI = Spectrum.from_phi(phiI, ns, (xx, xx))

    fs = P * fsN + (1 - P) * fsI
    return fs

def twoepoch_het_mig(params, ns, pts):
    """
    Heterogeneous model with two populations and migration.
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


def split_sizechange_second_hetero_asym_mig(params, ns, pts):
    """
    Split into two populations

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T1: Time in the past of split (in units of 2*Na generations) and before growth (T2)

    Instantanous size change followed by size change both pops.

    nu1T1, nu2T1 : Ratio of population size after instantanous change to ancient
         population size
    nu1T2, nu2T2 : Ratio of contemporary to ancient population size
    T2: Time in the past at which growth began
       (in units of 2*Na generations)
    m12: migration rate from pop1 into pop2 over T2
    m21: migration rate from pop2 into pop1 over T2
    me12: effective migration rate from pop1 into pop2 over T2
    me21: effective migration rate from pop2 into pop1 over T2
    P: The proportion of the genome evolving neutrally
    ns: Number of samples in resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nu1T1, nu2T1, nu1T2, nu2T2, T1, T2, m12, m21, me12, me21, P = params

    xx = Numerics.default_grid(pts)

    phiN = PhiManip.phi_1D(xx)
    phiN = PhiManip.phi_1D_to_2D(xx, phiN)

    phiN = Integration.two_pops(phiN, xx, T1, nu1=nu1T1, nu2=nu2T1, m12=0, m21=0)
    phiN = Integration.two_pops(phiN, xx, T2, nu1=nu1T2, nu2=nu2T2, m12=m12, m21=m21)

    fsN = Spectrum.from_phi(phiN, ns, (xx, xx))

    phiI = PhiManip.phi_1D(xx)
    phiI = PhiManip.phi_1D_to_2D(xx, phiI)

    phiI = Integration.two_pops(phiI, xx, T1, nu1=nu1T1, nu2=nu2T1, m12=0, m21=0)
    phiI = Integration.two_pops(phiI, xx, T2, nu1=nu1T2, nu2=nu2T2, m12=me12, m21=me21)

    fsI = Spectrum.from_phi(phiI, ns, (xx, xx))

    fs = P * fsN + (1 - P) * fsI
    return fs


# Previous
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

def model_split1_then_23(params, ns, pts):
    """
    A 3-population split: first pop1 vs (pop2, pop3), then pop2 vs pop3.

    params = [nu1, nu2, nu3, T1_23, T23]
        nu1, nu2, nu3: population sizes of pop1, pop2, pop3
        T1_23: time since pop1 split from ancestor of (pop2, pop3)
        T23: time since pop2 and pop3 split from each other
    ns: sample sizes (n1, n2, n3)
    pts: grid sizes for numerical integration
    """
    nu1, nu2, nu3, T1_23, T23 = params

    xx = Numerics.default_grid(pts)

    # Start from 1D ancestral population
    phi = PhiManip.phi_1D(xx)

    # Split into pop1 and ancestral pop23
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    # Integrate until pop2 and pop3 split
    phi = Integration.two_pops(phi, xx, T1_23, nu1=nu1, nu2=1.0)

    # Split ancestral pop23 into pop2 and pop3 (pop1 stays untouched)
    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)

    # Final integration with 3 populations
    phi = Integration.three_pops(phi, xx, T23, nu1=nu1, nu2=nu2, nu3=nu3)

    # Convert phi to frequency spectrum
    fs = Spectrum.from_phi(phi, ns, (xx, xx, xx))
    return fs

def model_split2_then_13(params, ns, pts):
    """
    Three-population model:
        - First, pop2 splits from ancestor of (pop1, pop3)
        - Then, pop1 and pop3 split

    params = [nu1, nu2, nu3, T2_13, T13]
        nu1, nu2, nu3: final population sizes
        T2_13: time since pop2 split from ancestor of (1,3)
        T13: time since pop1 and pop3 split from each other

    ns: sample sizes [n1, n2, n3]
    pts: grid sizes for numerical integration
    """
    nu1, nu2, nu3, T2_13, T13 = params

    xx = Numerics.default_grid(pts)

    # Start from 1D ancestral population
    phi = PhiManip.phi_1D(xx)

    # First split: pop2 and ancestor of (pop1, pop3)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, T2_13, nu1=1.0, nu2=1.0)  # ancestral sizes

    # Second split: ancestor splits into pop1 and pop3
    phi = PhiManip.phi_2D_to_3D_split_1(xx, phi)

    # Final integration of 3 pops
    phi = Integration.three_pops(phi, xx, T13, nu1=nu1, nu2=nu2, nu3=nu3)

    fs = Spectrum.from_phi(phi, ns, (xx, xx, xx))
    return fs

def model_split3_then_12(params, ns, pts):
    """
    Three-population model:
        - First, pop3 splits from ancestor of (pop1, pop2)
        - Then, pop1 and pop2 split

    params = [nu1, nu2, nu3, T3_12, T12]
        nu1, nu2, nu3: final population sizes
        T3_12: time since pop3 split from ancestor of (1,2)
        T12: time since pop1 and pop2 split

    ns: sample sizes [n1, n2, n3]
    pts: grid sizes
    """
    nu1, nu2, nu3, T3_12, T12 = params

    xx = Numerics.default_grid(pts)

    # Start from 1D ancestral population
    phi = PhiManip.phi_1D(xx)

    # First split: pop3 and ancestor of (pop1, pop2)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, T3_12, nu1=1.0, nu2=1.0)  # ancestral sizes

    # Second split: (1,2) split
    phi = PhiManip.phi_2D_to_3D_split_1(xx, phi)

    # Final integration of all three pops
    phi = Integration.three_pops(phi, xx, T12, nu1=nu1, nu2=nu2, nu3=nu3)

    fs = Spectrum.from_phi(phi, ns, (xx, xx, xx))
    return fs

def model_three_way_split(params, ns, pts):
    """
    Simultaneous 3-way split from a common ancestor.

    params = [nu1, nu2, nu3, T]
        nu1, nu2, nu3: final population sizes
        T: time since all 3 populations split simultaneously

    ns: sample sizes [n1, n2, n3]
    pts: grid sizes
    """
    # TODO: work out the correct phi manipulations for this model
    return None
