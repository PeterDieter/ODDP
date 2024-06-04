#!/bin/bash

# loop over all parameter settings
# print the filename (so we have some visual progress indicator)
# then submit the jobs to SLURM

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
## directories for logs
LOG=${WORKDIR}/log
## make dir if necessary
if [[ ! -d ${LOG} ]]; then
	mkdir -p ${LOG}
fi

#============================================================

# all parameter settings
arr_inst=( "zip" "grid" )
arr_maxwait=( 1200 1500 1800 )
arr_ameth=( "n" "w" "r" )
arr_rmeth=( "l" "s" )
arr_bundl=( "" "b" )
arr_alpha=( 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 1.00 )
arr_beta=( 0.50 0.60 0.70 0.80 0.90 )

#============================================================

# loop over all combinations
for inst in "${arr_inst[@]}"; do
	grid_inst="false"
	if [ $inst = "grid" ]; then
		grid_inst="true"
	fi
	for wait in "${arr_maxwait[@]}"; do
		for ameth in "${arr_ameth[@]}"; do
			for rmeth in "${arr_rmeth[@]}"; do
				for bundle in "${arr_bundl[@]}"; do
					use_bundling="false"
					bundl_arg="-"
					if [ ! -z ${bundle} ]; then
						use_bundling="true"
						bundl_arg="${bundle}"
					fi
					for alpha in "${arr_alpha[@]}"; do
						use_alpha="${alpha}"
						if [ ! ${ameth} = "w" ]; then
							use_alpha=""
						fi
						for beta in "${arr_beta[@]}"; do
							# --------------------------------------------------------
							use_beta="${beta}"
							if [ ! ${rmeth} = "l" ]; then
								use_beta=""
							fi
							# skip if result file already exits
							stat_file="${grid_inst}_${wait}_${ameth}_${rmeth}_${use_bundling}_${use_alpha}_${use_beta}.txt"
							if [[ -e "${RES}/${stat_file}" ]]; then
								echo "SKIPPED: \"${stat_file}\" already exists"
								continue
							fi
							# skip unnecessary method/alpha combinations
							if [[  -z ${use_alpha} && "${alpha}" != "${arr_alpha[0]}" ]]; then
								echo "schedule: ${inst} ${wait} ${ameth} ${rmeth} ${bundl_arg} ${alpha} ${beta} --- SKIPPED: Non-w AMeth with an alpha already started"
								continue
							fi
							# skip unnecessary method/beta combinations
							if [[  -z ${use_beta} && "${beta}" != "${arr_beta[0]}" ]]; then
								echo "schedule: ${inst} ${wait} ${ameth} ${rmeth} ${bundl_arg} ${alpha} ${beta} --- SKIPPED: Non-l RMeth with an beta already started"
								continue
							fi
							# --------------------------------------------------------
							## specify log files
							OUT_FILE=${LOG}/dafd_${inst}_${wait}_${ameth}_${rmeth}_${use_bundling}_${alpha}_${beta}.res
							ERR_FILE=${LOG}/dafd_${inst}_${wait}_${ameth}_${rmeth}_${use_bundling}_${alpha}_${beta}.err
							
							sbatch --output=${OUT_FILE} --error=${ERR_FILE} "$WORKDIR/dafd.sbatch" ${inst} ${wait} ${ameth} ${rmeth} ${bundl_arg} ${alpha} ${beta}
							echo "schedule: ${inst} ${wait} ${ameth} ${rmeth} ${bundl_arg} ${alpha} ${beta}"
			
							sleep 0.5 # pause to be kind to the scheduler
							# --------------------------------------------------------
						done
					done
				done
			done
		done
	done
done
