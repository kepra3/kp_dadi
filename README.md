# kp_dadi

*This branch is in development - not ready for public use! Please email kepra3@gmail.com for enquiries.*

*An electronic notebook by Katharine Prata (originally posted June 2021)*

A workflow for manually optimising the Diffusion Approximations for Demographic Inference (dadi) 
[(Gutenkunst et el. 2009)](https://dx.plos.org/10.1371/journal.pgen.1000695).

The scripts presented here were altered from the original dadi example [scripts](https://bitbucket.org/gutenkunstlab/dadi/src/master/) and are written for ease of interpretation.
This workflow has followed the advice offered by Ryan Gutenkunst on the [dadi-user google group](https://groups.google.com/g/dadi-user). Code was initially written for the manuscript, [Prata et al., 2022](https://onlinelibrary.wiley.com/doi/10.1111/mec.16391) published in *Molecular Ecology* (checkout: dc97a7f on *3 June 2021*).

This workflow will run through the example from [Prata et al (2022)](https://onlinelibrary.wiley.com/doi/10.1111/mec.16391).

## 1 - Preparation

A conda environment was set up with the versions below to ensure the analysis is compatible with these custom scripts.

```bash
$ conda create -n kp_dadi python=3.10.9 -y
$ conda activate kp_dadi
$ conda install -c conda-forge dadi=2.3.3
```

### 1a - What filters have been used on the data file?

For [Prata et al., (2022)](https://onlinelibrary.wiley.com/doi/10.1111/mec.16391):
A .vcf (v.4.0) file was obtained from an [ipyrad assembly]( https://ipyrad.readthedocs.io/en/latest/index.html ) (v.0.7.23) and was filtered to remove low frequency alleles in [VCFtools]( https://vcftools.github.io/index.html ) (v.0.1.16) as well as multiple other steps (see Methods in ms)

```bash
$ vcftools --vcf agrahno_2c.vcf --min-alleles 2 --max-alleles 2 --mac 3 --recode --stdout > 'AG1-AG2.vcf'
```

This filtering was conducted in order remove low confidence sites from the vcf file. Thus, because this step was done 
the low frequency alleles within the site frequency spectrum were masked in the following analysis.

Individuals within the populations (AG1, AG2, AL1 and AL2) were chosen based upon maximising the number of shared sites 
used in the analysis and which had low to zero admixture.

The all vcf files were named with the population names separated by a hyphen, (e.g., `AG1-AG2.vcf`).

### 1b - Create a frequency spectrum (subsample or projection)

To create the fs for input into dadi the format needs to be changed. Projection or subsampling can be used in order to 
maximise the number of SNP loci (& individuals) in the analysis when there is missing data.

Because inbreeding was used as a parameter in the models the subsample option was used [(Blischak et al., 2020)]( https://doi-org.ezproxy.library.uq.edu.au/10.1093/molbev/msaa042 ).

**To run all scripts in the repository:**

- Change to the script directory

**To run this script:**

-  Make sure the vcf file is named correctly (i.e., name of pop1 and pop2 separated by a hyphen) and is within 'data/vcf/' directory and the corresponding population file is within the 'data/popfile/' directory.

E.g.,
```bash
# Run script
$ python make_fs.py AG1-AG2 folded low subsample 20 9
# Find help
$ python make_fs.py -h
```
Arguments: 

(1) vcf file with the population IDs separated by a hyphen ('-'),

(2) whether spectra is folded or unfolded,

(3) whether spectra is masked or not ('low' - for masking on singletons and doubletons (e.g., if filtered as in step 1a), 'mid' - for masking at 0.5 freq bin if paralogs are an issue, or 'no' - for no masking)

(4) whether the method to use is 'subsample', 'projection' or no, 

(5) the number of genotypes (diploid) to project/subsample down to or use the max if using all individuals within the analysis (i.e., not subsampling/projecting).

This script will provide info about the sfs as well as create and plot the spectrum.

Output files: (1) a fs file (2) spectrum plot(s), 2 x if masked and (3) statistics of the sfs.

Due to a random fs being produced each time the official analysis results for [Prata et al (2022)]() can be found in the following directories: 
(1) `data/official_analysis_fs/` (2) `plots/official_ms_plots/` and (3) `results/official_analysis_results/`.

### 1bi - Plotting, projecting or folding a spectrum from a .fs file

If you already have an sfs output from another program (e.g., ANGSD, easySFS) simply use `plot_fs.py`

```bash
# Run script
$ python plot_fs.py AG1-AG2_subsampled no no low
# NOTE: only use folded if the spectrum is unfolded and you want to fold. Projection will project to 0.8 of sample size.
# Find help
$ python plot_fs.py -h
```

## 2 - Construct models

Models are defined through python functions provided within `demo_models_kp.py`.

Make sure this script is within the same directory you are performing the analysis in (scripts), in order to import the 
script as a module (or by however you like importing modules).

When creating a new model or wanting edit PTS or your upper and lower bounds edit these in `SETTINGS.py`.

## 3 - Optimise models

Using the `optimise_manual.py` script:

E.g.,
```bash
# Run script
$ python optimise_manual.py AG1-AG2 iso_inbred low subsample 3 '../results/test-' -p 1 1 1 1 1 
# Find help
$ python optimise_manual.py -h
```
Arguments: (1) name of your sfs, (2) model, (3) masked low, mid or no, (4) method: subsample or projection or no, (5) number of 
folds, (6) path for output results file and (optional-7) initial parameter values (where you start the optimisation)

If using these scripts please make sure your `PTS` for model extrapolation are defined appropriately in `SETTINGS.py`. You can use 
three grid sizes with the first being much larger (double) than your largest population, and the second and third, each larger than the previous (e.g., if largest population size is 20 haplotypes, the PTS = [40, 50, 60]). The nicknames for the models are also defined within `SETTINGS.py` - please check what they should be and check whether your upper and lower parameter limits are ok.  You may need to change PTS and upper/lower limits if there are issues with optimisation.

The aim of optimisation is to find parameters with the highest likelihood from many runs that have used different starting parameters. It is important to search the parameter space adequately. 

In [Prata et al., 2022](https://onlinelibrary.wiley.com/doi/10.1111/mec.16391) this was achieved through starting with random parameter values and running the optimisation with 3-fold perturbations many times (>100 per model) and then with 2-folds around multiple optimas within each model many times (>100), running from the lower and higher parameter limits, changing the parameter limits to space where the search hasn't been, and finally starting with the parameter values with the lowest log-likelihood scores using 1-fold perturbations (>50). During these optimisations how the likelihood changed with different parameter values was often watched to look for increased likelihood directions and the optimiser getting stuck in parameter space. The optimised parameter values were chosen due to appearing a few times (3x) when starting with different initial values. 

See Methods in [Prata et al. (2022)](https://onlinelibrary.wiley.com/doi/10.1111/mec.16391) for more details and see `dadi_optimisations_official.txt` for final optimisation results.

Official analysis results can be found in `results/official_analysis_results/`.

### 3a -Plot optimisation results - TBD

## 4 - Assess the model fit

### 4a - Qualitatively assess the model

To qualitatively assess the model a simulated fs was created from the parameter values and was compared to the data. 
This is helpful because you can see which snp bins your model is doing badly at estimating. A bad fit is when there is 
bias in certain snp bins - this tells you that another demographic model may be more appropriate.

E.g.,
```bash
# Run script
$ python compare_model.py AG1-AG2 iso_inbred folded 0.05 -o 2.122 25.95 0.0012 0.0455 0.3989
# Find help
$ python compare_model.py -h
```
Arguments: (1) fs, (2) model, (3) whether spectra is 'folded' or 'unfolded', (4) vmin (the scale bar), (5) -o optimised parameter values.

Output: (1) fs plot of the data, (2) fs plot of the simulated model, (3) a residual plot of data - model and (4) all 
previous 3 outputs as well a distribution of the residual scores in a 2x2 figure plot.

Official analysis results can be found in `plots/official_ms_plots/`.

### 4b - Bootstrapping

[*CURRENTLY IN CONSTRUCTION FOR NEW UPDATE*]

Either non-parametric or parametric bootstrapping. Non-parametric is taking random samples from your data 
whereas, parametric bootstrapping is using the model to simulate bootstraps of the data.

Here, we use non-parametric bootstrapping. See manual and dadi-user group for more info on bootstrapping.

```bash
$ python nonparametric_bootstrap_subsample.py AG1-AG2 iso_inbred 100 -g 20 9 -o 2.122 25.95 0.0012 0.0455 0.3989
```
Arguments: (1) fs, (2) model, (3) number of bootstraps, (4) -g number of genotypes in each pop, (5) -o optimised 
parameters.

Official analysis results can be found in `results/official_analaysis_results/bootstrap_vcf_official/`.

### 4c - Goodness-of-fit and parameter confidence intervals (GIM/FIM)

[*CURRENTLY IN CONSTRUCTION FOR NEW UPDATE*]

Using the information gathered from the bootstrapping statistics, bootstrap likelihoods given the optimised parameters
can be plotted next to each other using this Rscript below.

E.g.,
```bash
$ Rscript GOF_plots.R AG1-AG2 iso_inbred -551.02 1040.08
```
Arguments: (1) fs, (2) model, (3) optimised $$\chi^2$$ and (4) optimised log-likelihood.

Official analysis plots can be found in `plots/official_ms_plots/GOF_plots/`.

To calculate the confidence intervals the parameters the Godambe Information Matrix 
[(Coffman et al. 2016)]( https://doi.org/10.1093/molbev/msv255 ) was used. Different values of eps were used to find stable 
parameter confidence intervals.

```bash
$ python confindence_intervals.py AG1-AG2 iso_inbred 100 0.01 -o 2.122 25.95 0.0012 0.0455 0.3989
```

Official analysis results can be found in `results/official_analaysis_results/confidence_intervals_official/`.

