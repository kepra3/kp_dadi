import moments
import numpy as np

def model_func(params, ns):
	nu_1, nu_2, t1, nu11, nu12, m1_12, m1_21 = params
	_Nanc_size = 1.0  # This value can be used in splits with fractions
	sts = moments.LinearSystem_1D.steady_state_1D(np.sum(ns))
	fs = moments.Spectrum(sts)
	fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
	migs = np.array([[0, m1_12], [m1_21, 0]])
	fs.integrate(tf=t1, Npop=[nu11, nu12], m=migs, dt_fac=0.01)
	return fs

data = moments.Spectrum.from_file('/scratch/user/uqkprat2/analysis/kp_dadi/data/fs/group1-Amil_projected0.8.fs')
ns = data.sample_sizes

p0 = [1.159455095227221, 1.784372594090734, 4.972450136602187, 14.85079922024904, 11.349182559020194, 0.030355445688900046, 0.09190565296702148]
lower_bound = [0.01, 0.01, 1e-15, 0.01, 0.01, 0.0, 0.0]
upper_bound = [100.0, 100.0, 5.0, 100.0, 100.0, 10.0, 10.0]
model = model_func(p0, ns)
ll_model = moments.Inference.ll_multinom(model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))

theta = moments.Inference.optimal_sfs_scaling(model, data)
print('Optimal value of theta: {0}'.format(theta))

# As no theta0 or mut. rate + seq. length are not set
theta0 = 1.0
Nanc = int(theta / theta0)
print('Size of ancestral population: {0}'.format(Nanc))


plot_ns = [4 for _ in ns]  # small sizes for fast drawing
gen_mod = moments.ModelPlot.generate_model(model_func,
                                           p0, plot_ns)
moments.ModelPlot.plot_model(gen_mod,
                             save_file='model_from_GADMA.png',
                             fig_title='Demographic model from GADMA',
                             draw_scale=False,
                             pop_labels=['group1', 'Amil'],
                             nref=None,
                             gen_time=1.0,
                             gen_time_units='generations',
                             reverse_timeline=True)