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
arr_maxwait=( 1200 1500 )
arr_ameth=( "w" "n" "r" )
arr_rmeth=( "l" "s" "n" )
arr_bundl=( "b" "" )
arr_alpha=( 0.60 0.70 )
arr_beta=( 0.80 0.90 ) # 0.8 )

# loop over all combinations
for inst in "${arr_inst[@]}"; do
	grid_inst="false"
	if [ $inst = "grid" ]; then
		grid_inst="true"
	fi
	for wait in "${arr_maxwait[@]}"; do
		for ameth in "${arr_ameth[@]}"; do
			ameth_w=0
			if [ ${ameth} = "w" ]; then
				ameth_w=1
			fi
			for rmeth in "${arr_rmeth[@]}"; do
				rmeth_l=0
				if [ ${rmeth} = "l" ]; then
					rmeth_l=1
				fi
				for bundle in "${arr_bundl[@]}"; do
					use_bundling="false"
					if [ ! -z ${bundle} ]; then
						use_bundling="true"
					fi
					for alpha in "${arr_alpha[@]}"; do
						for beta in "${arr_beta[@]}"; do
							stat_file="${grid_inst}_${wait}_${ameth}_${rmeth}_${use_bundling}_${alpha}_${beta}.txt"
							if [[ -e "${RES}/${stat_file}" ]]; then
								echo "skipped, ${stat_file} exists"
								continue
							fi
							echo $stat_file
							## specify output files
							OUT_FILE=${RES}/dafd_${inst}_${wait}_${ameth}_${rmeth}_${bundle}_${alpha}_${beta}.res
							ERR_FILE=${RES}/dafd_${inst}_${wait}_${ameth}_${rmeth}_${bundle}_${alpha}_${beta}.err
							
							#sbatch --output=${OUT_FILE} --error=${ERR_FILE} "$WORKDIR/dafd.sbatch" ${inst} ${wait} ${ameth} ${rmeth} ${bundle} ${alpha} ${beta}
							#echo "setting (script): ${inst} ${wait} ${ameth} ${rmeth} ${bundle} ${alpha} ${beta}"
			
							sleep 0.0 # pause to be kind to the scheduler
						done
					done
				done
			done
		done
	done
done
