#!/bin/bash
# Written by Anders Lillevik Thorsen, August 2022
# Concatenate GM and WM signal intensity maps and then calculate the mean value over all subjects

cd /data/OBIC/Freesurfer/estimates

# Concatenate left GM surfaces
mri_concat \
--f /data/OBIC/Freesurfer/estimates/lh.GM_nu.avg.gm_fsaverage.txt \
--o /data/OBIC/Freesurfer/estimates/GWC_lh_mean.mgz \
--mean

# Estimate mean of left GM surfaces


/data/OBIC/Freesurfer/ICA/lh.nu.w-g.avg_fsaverage.mgh # concatenated file for lh

# Calculate mean of GWC for lh
mri_concat \
--f /data/OBIC/Freesurfer/GWC_w-g.avg_fsaverage.txt \
--o /data/OBIC/Freesurfer/estimates/GWC_lh_mean.mgz \
--mean
