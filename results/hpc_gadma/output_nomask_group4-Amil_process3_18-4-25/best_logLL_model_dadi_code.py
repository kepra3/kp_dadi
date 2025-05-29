import dadi
import numpy as np

def model_func(params, ns, pts):
	nu_1, nu_2, t1, nu11, nu12, m1_12, m1_21, t2, nu21, nu22, m2_12, m2_21 = params
	_Nanc_size = 1.0  # This value can be used in splits with fractions
	xx = dadi.Numerics.default_grid(pts)
	phi = dadi.PhiManip.phi_1D(xx)
	phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
	phi = dadi.Integration.two_pops(phi, xx, T=t1, nu1=nu11, nu2=nu12, m12=m1_12, m21=m1_21)
	phi = dadi.Integration.two_pops(phi, xx, T=t2, nu1=nu21, nu2=nu22, m12=m2_12, m21=m2_21)
	sfs = dadi.Spectrum.from_phi(phi, ns, [xx]*len(ns))
	return sfs

data = dadi.Spectrum.from_file('/scratch/user/uqkprat2/analysis/kp_dadi/data/fs/group4-Amil_projected0.8.fs')
pts = [50, 70, 90]
ns = data.sample_sizes

p0 = [3.175971533938209, 0.2688140671180461, 4.192680986281217, 1.5185892894556343, 12.667765329402084, 0.0, 0.0, 3.3666774693619104, 30.64860189251811, 11.61541174493469, 0.021510729167047507, 0.20756097798090223]
lower_bound = [0.01, 0.01, 1e-15, 0.01, 0.01, 0.0, 0.0, 1e-15, 0.01, 0.01, 0.0, 0.0]
upper_bound = [100.0, 100.0, 5.0, 100.0, 100.0, 10.0, 10.0, 5.0, 100.0, 100.0, 10.0, 10.0]
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
