#!/bin/bash
# Written by Anders Lillevik Thorsen, December 2021
# A port from tcsh to bash of UiO's script for calculating gray-white matter contrast
# Takes a single subject's processed Freesurfer directory as input

# Setup relevant variables
export SUBJECTS_DIR=/data/OBIC/Freesurfer
dcm=nu
s=$1
target=fsaverage
fwhm=10

echo "Subject is $s"
echo "SUBJECTS_DIR is $SUBJECTS_DIR"

cd $SUBJECTS_DIR
cd $s

# Check if file already exists, if so don't rerun GWC estimation

if [ ! -f "surf/lh.nu.avg.gm.mgh" ]; then

# extract WM/GM intensity from dcm.mgz
for hemi in lh rh
	do

	# compute points for gray matter
        mri_vol2surf \
                --mov $SUBJECTS_DIR/${s}/mri/$dcm.mgz \
                --hemi $hemi --noreshape --interp trilinear \
                --projfrac-avg 0.1 0.6 0.1 \
                --o $SUBJECTS_DIR/$s/surf/${hemi}.${dcm}.avg.gm.mgh --regheader $s

	# compute points for white matter
        mri_vol2surf \
                --mov $SUBJECTS_DIR/$s/mri/$dcm.mgz \
                --hemi $hemi --noreshape --interp trilinear \
                --projdist-avg -1.5 -0.15 0.15 \
                --o $SUBJECTS_DIR/$s/surf/${hemi}.${dcm}.avg.wm.mgh --regheader $s


        # computing the contrast
        mri_concat \
                --i $SUBJECTS_DIR/${s}/surf/${hemi}.${dcm}.avg.wm.mgh \
                --i $SUBJECTS_DIR/${s}/surf/${hemi}.${dcm}.avg.gm.mgh \
                --paired-diff-norm --mul 100 \
               	--o $SUBJECTS_DIR/${s}/surf/${hemi}.${dcm}.w-g.avg.mgh

	done

# surf -> fsaverage 
# Register surface from native space to fsaverage space
for meas in w-g.avg avg.wm avg.gm
do
	        for hemi in lh rh
		do
			# Tranform to fsaverage
	                mri_surf2surf \
	                        --srcsubject $s \
	                        --srchemi $hemi \
	                        --srcsurfreg sphere.reg \
	                        --trgsubject fsaverage \
	                        --trghemi $hemi \
	                        --trgsurfreg sphere.reg \
	                        --tval $SUBJECTS_DIR/${s}/surf/${hemi}.${dcm}.${meas}_fsaverage.mgh \
	                        --sval $SUBJECTS_DIR/${s}/surf/${hemi}.${dcm}.${meas}.mgh \
	                        --noreshape \
	                        --no-cortex

        # get surface data from aparc
                #mris_anatomical_stats \
                        #-t $SUBJECTS_DIR/${s}/surf/${hemi}.${dcm}.${meas}.mgh \
                        #-f $SUBJECTS_DIR/${s}/stats/${hemi}.${dcm}.${meas}.aparc.stats \
                        #-a aparc $s $hemi

       	 	done
done


# smooth surface
for meas in w-g.avg
do
        for hemi in lh rh
	do

                mri_surf2surf \
                        --hemi $hemi \
                        --s $target \
                        --sval $SUBJECTS_DIR/${s}/surf/${hemi}.${dcm}.${meas}_fsaverage.mgh \
                        --fwhm $fwhm \
                        --cortex \
                        --tval $SUBJECTS_DIR/${s}/surf/${hemi}.${dcm}.${meas}_fsaverage.${fwhm}.mgh \


        done
done
else
echo "The subjec's GWC file already exists, exiting script"
fi
