#!/bin/bash

# loop over all parameter settings
# print the filename (so we have some visual progress indicator)
# then submit the tophat jobs to SLURM

## working directory
BASEDIR=/scratch/hpc-prf-pose/peter
APPL=dafd_project/DAFD
WORKDIR=${BASEDIR}/${APPL}

## go to application directory
cd $WORKDIR

## directories for output
RES=${WORKDIR}/results
## make dir if necessary
if [[ ! -d ${RES} ]]; then
	mkdir -p ${RES}
fi

#============================================================

# all parameter settings
arr_inst=( "zip" "grid" ) # 2 3 )
arr_maxwait=( 1200 ) # "1500" )
arr_ameth=( "w" ) # "n" "r" "s" )
arr_rmeth=( "l" ) # "s" "n" )
arr_bundl=( "b" ) # "" )
arr_alpha=( 0.6 ) # 0.8 )
arr_beta=( 0.9 ) # 0.8 )

# loop over all combinations
for inst in "${arr_inst[@]}"; do
	for wait in "${arr_maxwait[@]}"; do
		for ameth in "${arr_ameth[@]}"; do
			for rmeth in "${arr_rmeth[@]}"; do
				for bundle in "${arr_bundl[@]}"; do
					for alpha in "${arr_alpha[@]}"; do
						for beta in "${arr_beta[@]}"; do
							## specify output files
							OUT_FILE=${RES}/dafd_${inst}_${wait}_${ameth}_${rmeth}_${bundle}_${alpha}_${beta}.res
							ERR_FILE=${RES}/dafd_${inst}_${wait}_${ameth}_${rmeth}_${bundle}_${alpha}_${beta}.err
							
							sbatch --output=${OUT_FILE} --error=${ERR_FILE} "$WORKDIR/dafd.sbatch" ${inst} ${wait} ${ameth} ${rmeth} ${bundle} ${alpha} ${beta}
							echo "setting (script): ${inst} ${wait} ${ameth} ${rmeth} ${bundle} ${alpha} ${beta}"
			
							sleep 1 # pause to be kind to the scheduler
						done
					done
				done
			done
		done
	done
done
