import moments
import numpy as np

def model_func(params, ns):
	nu_1, nu_2, t1, nu11, nu12, m1_12, m1_21, t2, nu21, nu22, m2_12, m2_21 = params
	_Nanc_size = 1.0  # This value can be used in splits with fractions
	sts = moments.LinearSystem_1D.steady_state_1D(np.sum(ns))
	fs = moments.Spectrum(sts)
	fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
	migs = np.array([[0, m1_12], [m1_21, 0]])
	fs.integrate(tf=t1, Npop=[nu11, nu12], m=migs, dt_fac=0.01)
	migs = np.array([[0, m2_12], [m2_21, 0]])
	fs.integrate(tf=t2, Npop=[nu21, nu22], m=migs, dt_fac=0.01)
	return fs

data = moments.Spectrum.from_file('../../../data/fs/group1-group2_projected0.8.fs')
ns = data.sample_sizes

p0 = [0.9640108635389305, 0.9913303338879081, 1.5714565152758169, 0.7134859587448723, 14.62973054014538, 0.0, 0.0, 1.9787284083494563, 21.58552343165088, 19.509500885610986, 0.07472220153552123, 0.5127463579074869]
lower_bound = [0.01, 0.01, 1e-15, 0.01, 0.01, 0.0, 0.0, 1e-15, 0.01, 0.01, 0.0, 0.0]
upper_bound = [100.0, 100.0, 5.0, 100.0, 100.0, 10.0, 10.0, 5.0, 100.0, 100.0, 10.0, 10.0]
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
                             pop_labels=['group1', 'group2'],
                             nref=None,
                             gen_time=1.0,
                             gen_time_units='generations',
                             reverse_timeline=True)