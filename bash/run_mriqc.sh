#!/bin/bash -x
# Written by Anders Lillevik Thorsen, August 2021
# This scripts splits the use of mriqc for OBIC into sets of four subjects

headdir=/data/OBIC/BIDS
mriqcdir=/data/OBIC/BIDS/mriqc
listfile=${mriqcdir}/list.txt

	cd ${mriqcdir}

# Sets up performace parameters for mriqc NOT CURRENTLY USED

	nthreads=18 # Usable vCPUs
	mem=25 #Usable memory in Gb

# Sets block size of participants to be run at once

	blocksize=20

# Getting number of participants

	nrows=$(wc -l ${listfile} | awk '{print $1}')
	echo 'nrows is '${nrows}

# Iterates over blocks of participants

	for (( i=1; i <= ${nrows}; i=i+${blocksize} ))

			do
			
			echo 'line number/iteration number is '${i} # Prints line number/iteration number

			subjects=$(sed -n ''${i},$(( ${i}+${blocksize} ))p'' ${listfile}) # Gets participant names
			
			echo 'Now running '${subjects} # Prints subjects for troubleshooting

		# Runs mriqc
			docker run -it --rm -v ${headdir}:/data:ro -v ${mriqcdir}:/out poldracklab/mriqc /data /out participant --participant_label ${subjects} --no-sub #--n_proc ${nthreads} --mem_gb $mem

		#break # Breaks loop after one iteration
	done




