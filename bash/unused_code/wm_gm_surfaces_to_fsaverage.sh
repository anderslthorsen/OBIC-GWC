#!/bin/bash
# Written by Anders Lillevik Thorsen, August 2022
# This script registers each subject's mean gray and white matter signal intensity surfaces to fsaverage and smoothes them with a 10mm sphere.

export SUBJECTS_DIR=/data/OBIC/Freesurfer
dcm=nu
target=fsaverage
fwhm=10

cd /data/OBIC/Freesurfer

for subject in `cat GWC_include20Dec2021.txt`
	do
		bash /data/OBIC/scripts/extract_intensity_contrast_bash.sh ${subject}
	done


echo "Subject is $s"
echo "SUBJECTS_DIR is $SUBJECTS_DIR"

cd $SUBJECTS_DIR
cd $s

# Check if file already exists, if so don't rerun GWC estimation

if [ ! -f "surf/lh.nu.avg.gm.mgh" ]; then

#smooth
for meas in avg.wm avg.gm

do

        for hemi in lh rh
	do

                mri_surf2surf \
                        --hemi $hemi \
                        --s $target \
                        --sval $SUBJECTS_DIR/${subject}/surf/${hemi}.${dcm}.${meas}_fsaverage.mgh \
                        --fwhm $fwhm \
                        --cortex \
                        --tval $SUBJECTS_DIR/${subject}/surf/${hemi}.${dcm}.${meas}_fsaverage.${fwhm}.mgh \


        done
done
else
echo "The subjec's GWC file already exists, exiting script"
fi
