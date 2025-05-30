#!/bin/bash

pop1=$1
pop2=$2
option=$3

if [ $option == "project" ]; then
  echo "Projecting SFS"
  # plot singular and rename
  python plot_fs.py "../data/fs/outQ0.unfold.${pop1}_dadi_${pop2}_dadi_snp_Q5_ind834_sites_qc_dep_ngspara_n80_chrOnly_USE.fs" folded ${option} no
  mv "../data/fs/outQ0.unfold.${pop1}_dadi_${pop2}_dadi_snp_Q5_ind834_sites_qc_dep_ngspara_n80_chrOnly_USE_foldproject.fs" "../data/fs/${pop1}-${pop2}_projected0.8.fs"
  # loop to create bootstraps
  for i in {1..100};
  do python plot_fs.py "../data/fs/outQ0.unfold.${pop1}_dadi_${pop2}_dadi_snp_Q5_ind834_sites_qc_dep_ngspara_n80_chrOnly_USE.saf_boots/outQ0.unfold.${pop1}_dadi_${pop2}_dadi_snp_Q5_ind834_sites_qc_dep_ngspara_n80_chrOnly_USE.saf_${i}_boot.fs" folded projection no -o "projected_${pop1}-${pop2}";
    done
elif [ $option == "no" ]; then
  echo "Not projecting SFS"
  # plot singular and rename
  python plot_fs.py "../data/fs/outQ0.unfold.${pop1}_dadi_${pop2}_dadi_snp_Q5_ind834_sites_qc_dep_ngspara_n80_chrOnly_USE.fs" folded no no
  mv "../data/fs/outQ0.unfold.${pop1}_dadi_${pop2}_dadi_snp_Q5_ind834_sites_qc_dep_ngspara_n80_chrOnly_USE_fold.fs" "../data/fs/${pop1}-${pop2}.fs"

  # loop to create bootstraps
  for i in {1..100};
    do python plot_fs.py "../data/fs/outQ0.unfold.${pop1}_dadi_${pop2}_dadi_snp_Q5_ind834_sites_qc_dep_ngspara_n80_chrOnly_USE.saf_boots/outQ0.unfold.${pop1}_dadi_${pop2}_dadi_snp_Q5_ind834_sites_qc_dep_ngspara_n80_chrOnly_USE.saf_${i}_boot.fs" folded no no -o "${pop1}-${pop2}";
  done
else
  echo "Invalid option. Use 'project' or 'no'."
  exit 1
fi


# Make parameter file and hpc script
for i in {1..8};
  do ./customise_gadma_params.sh ${pop1} ${pop2} ${i} 8 4 "params-full-pop1-pop2-hpc_process" "hpc_files/params-full-${pop1}-${pop2}_process${i}";
done

for i in {1..8};
  do ./customise_gadma_slurm.sh ${pop1} ${pop2} ${i};
done

echo "Files made for ${pop1} and ${pop2}"