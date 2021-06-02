# kp_dadi
*Katharine Prata*

A pipeline for manually optimising the Diffusion Approximations for Demographic Inference (dadi) by Gutenkunst et el. 
(2009): a population genomics analyses.

Using altered [scripts](https://bitbucket.org/gutenkunstlab/dadi/src/master/) from 
[Gutenkunst et al. (2009)](https://dx.plos.org/10.1371/journal.pgen.1000695). This pipeline has followed the help 
offered by Ryan Gutenkunst on the [dadi-user google group](https://groups.google.com/g/dadi-user).

**Note:** the pipeline is for datasets created through reduced representation methods (GBS, ddRAD, RAD, etc.) and which 
do not have a reference genome.

## 1 - Preparation

Best practise to set up a conda environment with the version of the software you are using so that the analysis will 
work the same every time.

```
$ conda create -n dadi211 python=3.6 -y
$ conda activate dadi211
$ conda install -c conda-forge dadi=2.1.1
```

### 1a - Filter your vcf file

.vcf file was obtained from an [ipyrad assembly]( https://ipyrad.readthedocs.io/en/latest/index.html ) and was filtered 
to remove low frequency alleles in vcftools as well as other mutliple other steps see Prata et al., (preprint).

```bash
$ vcftools --vcf agrahno_2c.vcf --min-alleles 2 --max-alleles 2 --mac 3 --recode --stdout > 'agrahno_2d-3.vcf'
```

This filtering was conducted in order remove low confidence sites from the vcf file. However, if this step is done 
then you must mask low frequency alleles within your site frequency spectrum.

Choose individuals within your populations that have the highest number of shared sites in order to maximise the number 
of SNP loci used in your analyses and individuals with lowest (or no) admixture values.

My final vcf file was named `AG1-AG2.vcf` as AG1 is my first population and AG2 is my second population. 

### 1b - Create your frequency spectrum (subsample or projection?)

To create your fs for input into dadi the format will need to be changed and you need to decided where you will 
project or subsample your fs in order to maximise the number of SNP loci (& individuals) in your analysis.

Here change to your script directory and make sure your vcf file is named correctly (i.e., name of pop1 and pop2 separated 
by a hyphen) and within the data directory.

If you have to incorporate inbreeding as a parameter into your model then you can use the subsample option 
within dadi, if not then you can use projection.

```bash
$ python make_fs.py AG1-AG2 yes subsample 10 10
```

First argument is your vcf file with the population IDs separated by a dash ('-').
Second is whether you want to use masking on your singletons and doubletons (e.g., if you have filtered them as in 
step 1a)
Third is whether the method you want to use is subsample, projection or neither.
Fourth is the genotypes (diploid) you want to project/subsample down to, just use the max if you are using all indivs,
(i.e., not subsampling/projecting).

This script will give you info about the number of SNPs in your sfs as well as plot the spectrum.

## 2 - Construct your models

Models are defined through python functions. You can use any of the models provided here (demo_models_kp.py).
Make sure this script is within the same directory you are performing the analysis in (scripts), in order to import the 
script as a module (or by however you like importing modules). The dadi-user group as well as the manual provide plenty 
of instruction in how to create the models.

## 3 - Optimise your models

Using the optimise_manual.py script:

```bash
$ python optimise_manual.py AG1-AG2 iso_inbred yes subsample 3
```
Arguments: name of your sfs, model, masked or not, method: subsample or projection, folds

How many times you run the optimisation and what starting parameters you use is up to you.
The aim is to find parameters with the highest likelihood from many runs that have used different starting parameters.
It is important to search the parameter space adequately. I do this by running from the lower and higher parameter 
limits, random parameters and the optimised parameter values found in previous runs. Increasing the fold to 2 or 3 
allows you to jump around/explore the parameter space easier and using 1-fold allows you to hone in on the optimal 
values. Watching how the likelihood changes when different parameter values are used is important. Good-luck :).

## 4 - Assess your model fit

### 4a - Qualitatively assess your model

To qualitatively assess your model you create a simulated fs and compare it to your data. This is helpful because you 
can see which snp bins your model is doing badly at estimating. A bad fit is when there is bias in certain snp bins - 
this tells you that another demographic model may be more appropriate.

```bash
$ python compare_model.py
```

### 4b - Bootstrapping

You can either use non-parametric or parametric bootstrapping. Non-parametric is taking random samples from your data 
whereas, parametric bootstrapping is using your model to simulate bootstraps of your data.

Here, we use non-parametric bootstrapping. See manual and dadi-user group for more info on bootstrapping.

```bash
$ python nonparametric_bootstrap_subsample.py AG1-AG2 iso_inbred 100 -g 20 9 -o 2.122 25.95 0.0012 0.0455 0.3989
```
Arguments: pop(s), model, number of bootstraps, -g number of genotypes in each pop, -o optimised parameters.

### 4c - Goodness-of-fit and parameter confidence intervals (GIM/FIM)

Using the information gathered from the bootstrapping statistics. Bootstrap likelihoods given the optimised parameters
can be plotted next to each other.

```bash
$ Rscript # TODO
```

To calculate the confidence intervals of the parameters we use the Godambe Information Matrix (Coffman et al.). This 
step is finicky and requires the adjustment of different values of eps (see dadi manual, YRU_CEU.py script and dadi 
user google group).

```bash
$ python confindence_intervals.py AG1-AG2 iso_inbred 100 0.01 -o 2.122 25.95 0.0012 0.0455 0.3989
```