import dadi
import numpy as np

def model_func(params, ns, pts):
	nu_1, nu_2, t1, nu11, nu12, m1_12, m1_21 = params
	_Nanc_size = 1.0  # This value can be used in splits with fractions
	xx = dadi.Numerics.default_grid(pts)
	phi = dadi.PhiManip.phi_1D(xx)
	phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
	nu1_func = lambda t: nu_1 * (nu11 / nu_1) ** (t / t1)
	nu2_func = lambda t: nu_2 + (nu12 - nu_2) * (t / t1)
	phi = dadi.Integration.two_pops(phi, xx, T=t1, nu1=nu1_func, nu2=nu2_func, m12=m1_12, m21=m1_21)
	sfs = dadi.Spectrum.from_phi(phi, ns, [xx]*len(ns))
	return sfs

data = dadi.Spectrum.from_file('/Users/uqkprat2/git/kp_dadi/data/fs/group1-group2_projected0.8.fs')
data = data.project([16, 16])
pts = [40, 60, 80]
ns = data.sample_sizes

p0 = [11.596867481769653, 0.705547721973308, 1.4573203698391572, 9.831795168485781, 20.73565780796799, 0.0, 0.6484759043326622]
lower_bound = [0.01, 0.01, 1e-15, 0.01, 0.01, 0.0, 0.0]
upper_bound = [100.0, 100.0, 5.0, 100.0, 100.0, 10.0, 10.0]
func_ex = dadi.Numerics.make_extrap_log_func(model_func)
model = func_ex(p0, ns, pts)
ll_model = dadi.Inference.ll_multinom(model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))

theta = dadi.Inference.optimal_sfs_scaling(model, data)
print('Optimal value of theta: {0}'.format(theta))

# As no theta0 or mut. rate + seq. length are not set
theta0 = 1.0
Nanc = int(theta / theta0)
print('Size of ancestral population: {0}'.format(Nanc))
