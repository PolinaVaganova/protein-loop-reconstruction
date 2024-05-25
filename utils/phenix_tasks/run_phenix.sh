#!/bin/bash

source <path_to_Phenix>/phenix_env.sh

# phenix.refine picks crystallographic water molecules
# Note: to increase the speed the direct summation is turned off in the current .eff-files edition
# Uncomment those lines if needed

if [ ! -f final_water_1.pdb ]; then
    phenix.refine $1 $2 <path_to>/refine_water.eff --unused_ok --overwrite
fi

# phenix.refine B-factors

if [ ! -f final_b_factors_1.pdb ]; then
    phenix.refine final_water_1.pdb $2 <path_to>/refine_b_factors.eff --unused_ok --overwrite
fi


# phenix.molprobity calls
# Note: MolProbity uses n_gaussian scattering table, thus, R-values might deviate from the above

if [ ! -f validation.out ]; then
    phenix.molprobity $1 $2 output.prefix=validation
fi

if [ ! -f validation_final_water.out ]; then
    phenix.molprobity final_water_1.pdb $2 output.prefix=validation_final_water
fi

if [ ! -f validation_final_b_factors_1.out ]; then
    phenix.molprobity final_b_factors_1.pdb $2 output.prefix=validation_final_b_factors
fi
