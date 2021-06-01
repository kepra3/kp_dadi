# kp_dadi
*Katharine Prata*

A pipeline for manually optimising the Diffusion Approximations for Demographic Inference (dadi) by Gutenkunst et el. (2009): a population genomics analyses.

Using altered [scripts](https://bitbucket.org/gutenkunstlab/dadi/src/master/) from [Gutenkunst et al. (2009)](https://dx.plos.org/10.1371/journal.pgen.1000695). This pipeline has followed the help offered by Ryan Gutenkunst on the [dadi-user google group](https://groups.google.com/g/dadi-user).

**Note:** the pipeline is for datasets created through reduced representation methods (GBS, ddRAD, RAD, etc.) and which do not have a reference genome.

## 1 - Preparation

.vcf file was obtained from an [ipyrad assembly]( https://ipyrad.readthedocs.io/en/latest/index.html ) and was filtered to remove low frequency alleles in vcftools as well as other mutliple other steps see Prata et al., (preprint).

```bash
$ vcftools --vcf agrahno_2c.vcf --min-alleles 2 --max-alleles 2 --mac 3 --recode --stdout > 'agrahno_2d-3.vcf'
```

This was conducted in order remove low confidence sites from the vcf file. However, if this step is done then you must mask low frequency alleles within your site frequency spectrum.

## 1 - Subsample or projection

Choose individuals within your populations that have the highest number of shared sites, 
lowest (or no) admixture values in order to maxmise the number of SNP loci used in your analyses.
If you do not have to incorporate inbreeding as a parameter into your model then you can use the subsample option within dadi.

You can use this custom script in order to achieve this:
```bash
$ python make_fs.py Pop1-Pop2 yes subsample 10 10
```

First argument is your vcf file with the population IDs separated by a dash ('-').
Second is whether you want to use masking on your singletons and doubletons (e.g., if you have filtered them)
Third is whether the method you want to use is subsample, projection or neither.
Fourth is the genotypes (diploid) you want to project/subsample down to, just use the max if you are using all indivs.

This script will give you info about the number of SNPs in your sfs as well as plot the spectrum.

## 2 - Construct your models

Models are defined through python functions. You can use any of the models provided here (demo_models_kp.py).
Make sure this script is within the same directory you are performing the analysis in, in order to import the script as a module (or by however you like importing modules).
The dadi-user group as well as the manual provide plenty of instruction in how to create the models.

## 3 - Optimise your models

Using the optimise_manual.py script:

```bash
$ python optimise_manual.py 
```

How many times you run the optimisation and what starting parameters you use is up to you.
The aim is to find parameters with the highest likelihood from many runs that have used different starting parameters.
It is important to search the parameter space adequately. I do this by running from the lower and higher parameter limits,
random parameters and the optimised parameter values found in previous runs. Increasing the fold to 2 or 3 allows you to
jump around the parameter space easier.

## 4 - Assess your models

## 5 - Parameter uncertainties