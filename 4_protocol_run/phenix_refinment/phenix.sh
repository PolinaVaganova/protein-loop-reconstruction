#!/bin/bash

source /home/polina/tools/phenix-installer-1.19.2-4158-intel-linux-2.6-x86_64-centos6/phenix-1.19.2-4158/phenix_env.sh

# phenix.refine add water
if [ ! -f final_water_1.pdb ]; then
    phenix.refine $1 $2 ../refine_water.eff --unused_ok --overwrite nproc=18
fi


# phenix.refine B-factors
if [ ! -f final_b_factors_1.pdb ]; then
    phenix.refine final_water_1.pdb $2 ../refine_b_factors.eff --unused_ok --overwrite nproc=18
fi



# phenix.molprobity

if [ ! -f validation.out ]; then
    phenix.molprobity $1 $2 output.probe_dots=False output.maps=True output.coot=True output.prefix=validation nproc=18
fi

if [ ! -f validation_final_water.out ]; then
    phenix.molprobity final_water_1.pdb $2 output.probe_dots=False output.maps=True output.coot=True output.prefix=validation_final_water nproc=18
fi

if [ ! -f validation_final_b_factors_1.out ]; then
    phenix.molprobity final_b_factors_1.pdb $2 output.probe_dots=False output.maps=True output.coot=True output.prefix=validation_final_b_factors nproc=18
fi
