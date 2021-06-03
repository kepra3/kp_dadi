# kp_dadi
*An electronic notebook by Katharine Prata (June 2021)*

A workflow for manually optimising the Diffusion Approximations for Demographic Inference (dadi) 
[(Gutenkunst et el. 2009)](https://dx.plos.org/10.1371/journal.pgen.1000695).

The scripts presented here were altered from the original dadi 
[scripts](https://bitbucket.org/gutenkunstlab/dadi/src/master/). This workflow has followed the advice offered by 
Ryan Gutenkunst on the [dadi-user google group](https://groups.google.com/g/dadi-user).

**Note:** this workflow is for datasets created through reduced representation methods (GBS, ddRAD, RAD, etc.) and which 
do not have a reference genome. The scripts can be used on this type of data or customised to fit other analysis 
requirements.

## 1 - Preparation

A conda environment was setup with the versions below to ensure the analysis works the same every time it's run.

```
$ conda create -n dadi211 python=3.6 -y
$ conda activate dadi211
$ conda install -c conda-forge dadi=2.1.1
```

### 1a - Filter vcf file

A .vcf (v.4.0) file was obtained from an [ipyrad assembly]( https://ipyrad.readthedocs.io/en/latest/index.html ) 
(v.0.7.23) and was filtered to remove low frequency alleles in [VCFtools]( https://vcftools.github.io/index.html ) 
(v.0.1.16) as well as multiple other steps see Prata et al., (xxxx).

```bash
$ vcftools --vcf agrahno_2c.vcf --min-alleles 2 --max-alleles 2 --mac 3 --recode --stdout > 'agrahno_2d-3.vcf'
```

This filtering was conducted in order remove low confidence sites from the vcf file. Thus, because this step was done 
the low frequency alleles within the site frequency spectrum were masked.

Individuals within the populations (AG1, AG2, AL1 and AL2) were chosen based upon maximising the number of shared sites 
used in the analysis and which had low to zero admixture.

The final vcf files was named with the population names separated by a hyphen, (e.g., `AG1-AG2.vcf`).

### 1b - Create a frequency spectrum (subsample or projection)

To create the fs for input into dadi the format needs to be changed. Projection or subsampling can be used in order to 
maximise the number of SNP loci (& individuals) in the analysis when there is missing data.

Because inbreeding was used as a parameter in the models the subsample option was used 
[(Blischak et al., 2020)]( https://doi-org.ezproxy.library.uq.edu.au/10.1093/molbev/msaa042 ).

**To run this script:**

Change to the script directory and make sure the vcf file is named correctly (i.e., name of pop1 and pop2 separated 
by a hyphen) and within the data directory and the population file within the popfile directory.

E.g.,
```bash
$ python make_fs.py AG1-AG2 yes subsample 20 9
```
Arguments: (1) vcf file with the population IDs separated by a hyphen ('-'), (2) yes or no - whether to use 
masking on singletons and doubletons (e.g., if filtered as in step 1a), (3) whether the method to use is subsample, 
projection or neither, (4) the number of genotypes (diploid) to project/subsample down to or use the max if using all 
individuals within the analysis (i.e., not subsampling/projecting).

This script will provide info about the sfs as well as create and plot the spectrum.

Output files: (1) an fs file (2) spectrum plot(s), 2 x if masked and (3) statistics of the sfs.

Due to a random fs being produced each time the official results can be found in the following directories: 
(1) `data/official_analysis_fs/` (2) `plots/official_ms_plots/` and (3) `results/official_analysis_results/`.

## 2 - Construct models

Models are defined through python functions provided within `demo_models_kp.py`.

Make sure this script is within the same directory you are performing the analysis in (scripts), in order to import the 
script as a module (or by however you like importing modules).

## 3 - Optimise models

Using the `optimise_manual.py` script:

E.g.,
```bash
$ python optimise_manual.py AG1-AG2 iso_inbred yes subsample 3
```
Arguments: (1) name of your sfs, (2) model, (3) masked or not, (4) method: subsample or projection, and (5) number of 
folds.

The aim is to find parameters with the highest likelihood from many runs that have used different starting parameters.
It is important to search the parameter space adequately. 

This was achieved through starting with random parameter values and running the optimisation with 3-fold perturbations 
a few times and then with 2-folds a few times, running from the lower and higher parameter limits, changing the 
parameter limits, and starting with the parameter values which the lowest log-likelihood scores using 1-fold 
perturbations. During these optimisations how the likelihood changed with different parameter values was watched. The 
optimised parameter values were chosen due to appearing a few times when starting with different initial values. 
See Prata et al., xxxx for more details and see `dadi_optimisations_official.txt` for optimisation results.

Official results can be found in `results/official_analysis_results/`

## 4 - Assess the model fit

### 4a - Qualitatively assess the model

To qualitatively assess the model a simulated fs was created from the parameter values and was compared to the data. 
This is helpful because you can see which snp bins your model is doing badly at estimating. A bad fit is when there is 
bias in certain snp bins - this tells you that another demographic model may be more appropriate.

E.g.,
```bash
$ python compare_model.py AG1-AG2 iso_inbred 0.05 -o 2.122 25.95 0.0012 0.0455 0.3989
```
Arguments: (1) fs, (2) model, (3) vmin (the scale bar), (4) -o optimised parameter values.

Output: (1) fs plot of the data, (2) fs plot of the simulated model, (3) a residual plot of data - model and (4) all 
previous 3 outputs as well a distribution of the residual scores in a 2x2 figure plot.

Official results can be found in `plots/official_ms_plots/`

### 4b - Bootstrapping

Either non-parametric or parametric bootstrapping. Non-parametric is taking random samples from your data 
whereas, parametric bootstrapping is using the model to simulate bootstraps of the data.

Here, we use non-parametric bootstrapping. See manual and dadi-user group for more info on bootstrapping.

```bash
$ python nonparametric_bootstrap_subsample.py AG1-AG2 iso_inbred 100 -g 20 9 -o 2.122 25.95 0.0012 0.0455 0.3989
```
Arguments: (1) fs, (2) model, (3) number of bootstraps, (4) -g number of genotypes in each pop, (5) -o optimised 
parameters.

Official results can be found in `results/official_analaysis_results/bootstrap_vcf_official/`.

### 4c - Goodness-of-fit and parameter confidence intervals (GIM/FIM)

Using the information gathered from the bootstrapping statistics, bootstrap likelihoods given the optimised parameters
can be plotted next to each other using this Rscript below.

E.g.,
```bash
$ Rscript GOF_plots.R AG1-AG2 iso_inbred -551.02 1040.08
```
Arguments: (1) fs, (2) model, (3) optimised $$\chi^2$$ and (4) optimised log-likelihood.

Official results can be found in `results/official_analaysis_results/bootstrap_vcf_official/GOF_plots/`.

To calculate the confidence intervals the parameters the Godambe Information Matrix 
[(Coffman et al. 2016)]( https://doi.org/10.1093/molbev/msv255 ) was used. Different values of eps were used to find stable 
parameter confidence intervals.

```bash
$ python confindence_intervals.py AG1-AG2 iso_inbred 100 0.01 -o 2.122 25.95 0.0012 0.0455 0.3989
```

Official results can be found in `results/official_analaysis_results/confidence_intervals_official/`.

