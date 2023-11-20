#!/bin/bash
# Written by Anders Lillevik Thorsen, December 2021
# Estimates mean and standard deviation of GWC maps from all subjects

cd /data/OBIC/Freesurfer/estimates

# Calculate mean of GWC for lh
mri_concat \
--f /data/OBIC/Freesurfer/GWC_w-g.avg_fsaverage.txt \
--o /data/OBIC/Freesurfer/estimates/GWC_lh_mean.mgz \
--mean

/data/OBIC/Freesurfer/ICA/lh.nu.w-g.avg_fsaverage.mgh # concatenated file for lh

# Calculate mean of GWC for lh
mri_concat \
--f /data/OBIC/Freesurfer/GWC_w-g.avg_fsaverage.txt \
--o /data/OBIC/Freesurfer/estimates/GWC_lh_mean.mgz \
--mean


# Calculate standard deviation of GWC
mri_concat \
--f /data/OBIC/Freesurfer/GWC_w-g.avg_fsaverage.txt \
--o /data/OBIC/Freesurfer/estimates/GWC_lh_SD.mgz \
--std

# Reverse values of GWC to make map appear as in Norbom et al., 2020, Biol Psychiatry
# How to reverse view in Freeview?
mri_concat \
--i /data/OBIC/Freesurfer/estimates/GWC_lh_mean.mgz \
--o /data/OBIC/Freesurfer/estimates/GWC_lh_mean_reversed.mgz \
--mul -1

# Compute difference between mean maps of smoothed GWC at gm/wm stage or after calculcating the difference (GWC)
mri_concat \
--i /data/OBIC/Freesurfer/estimates/GWC_lh_mean.mgz \
--i /data/OBIC/Freesurfer/estimates/GWC_lh_smoothed_mean.mgz \
--o /data/OBIC/Freesurfer/estimates/GWC_lh_mean_diff_norm1.mgz \
--paired-diff-norm1

# Multiply difference by 100 to get percentage difference
mri_concat \
--i /data/OBIC/Freesurfer/estimates/GWC_lh_mean_diff_norm1.mgz \
--o /data/OBIC/Freesurfer/estimates/GWC_lh_mean_diff_norm1_mul100.mgz \
--mul 100

# Calculate mean GWC for lh and rh based on multi-subject volume

mri_concat \
--i /data/OBIC/Freesurfer/ICA/lh.nu.w-g.avg_fsaverage.mgh \
--o /data/OBIC/Freesurfer/estimates/GWC_lh_mean.mgz \
--mean

mri_concat \
--i /data/OBIC/Freesurfer/ICA/rh.nu.w-g.avg_fsaverage.mgh \
--o /data/OBIC/Freesurfer/estimates/GWC_rh_mean.mgz \
--mean

# Invert values so that more myelin/lower GWC has higher values
mri_concat \
--i /data/OBIC/Freesurfer/estimates/GWC_lh_mean.mgz \
--o /data/OBIC/Freesurfer/estimates/GWC_lh_mean_inverted.mgz \
--mul -1

mri_concat \
--i /data/OBIC/Freesurfer/estimates/GWC_rh_mean.mgz \
--o /data/OBIC/Freesurfer/estimates/GWC_rh_mean_inverted.mgz \
--mul -1

# Load mean GWC overlaid on pial surface from fsaverage
freeview -f /data/OBIC/Freesurfer/fsaverage/surf/lh.pial:overlay=/data/OBIC/Freesurfer/estimates/GWC_lh_mean.mgz -f ../fsaverage/surf/rh.pial:overlay=/data/OBIC/Freesurfer/estimates/GWC_rh_mean.mgz

# Load mean inversed GWC overlaid on pial surface from fsaverage
freeview -f /data/OBIC/Freesurfer/fsaverage/surf/lh.pial:overlay=/data/OBIC/Freesurfer/estimates/GWC_lh_mean_inverted.mgz -f /data/OBIC/Freesurfer/fsaverage/surf/rh.pial:overlay=/data/OBIC/Freesurfer/estimates/GWC_rh_mean_inverted.mgz


freeview -f /data/OBIC/Freesurfer/fsaverage/surf/lh.pial:overlay=/data/OBIC/Freesurfer/estimates/GWC_lh_mean.mgz -f /data/OBIC/Freesurfer/fsaverage/surf/lh.pial:overlay=/data/OBIC/Freesurfer/estimates/GWC_lh_mean_inverted.mgz

freeview -f /data/OBIC/Freesurfer/fsaverage/surf/lh.pial:overlay=/data/OBIC/Freesurfer/ICA/icasso_results_rerun_d7/lh.IC_03.mgh -f /data/OBIC/Freesurfer/fsaverage/surf/rh.pial:overlay=/data/OBIC/Freesurfer/ICA/icasso_results_rerun_d7/rh.IC_03.mgh

# Load mean GM and WM signal intensity on pial surface from fsaverage
# GM
freeview -f /data/OBIC/Freesurfer/fsaverage/surf/lh.pial:overlay=/data/OBIC/Freesurfer/estimates/GM_lh_mean.mgz -f ../fsaverage/surf/rh.pial:overlay=/data/OBIC/Freesurfer/estimates/GM_rh_mean.mgz

# WM
freeview -f /data/OBIC/Freesurfer/fsaverage/surf/lh.pial:overlay=/data/OBIC/Freesurfer/estimates/WM_lh_mean.mgz -f ../fsaverage/surf/rh.pial:overlay=/data/OBIC/Freesurfer/estimates/WM_rh_mean.mgz

# GM and WM for left hemisphere
freeview -f /data/OBIC/Freesurfer/fsaverage/surf/lh.pial:overlay=/data/OBIC/Freesurfer/estimates/GM_lh_mean.mgz -f /data/OBIC/Freesurfer/fsaverage/surf/lh.pial:overlay=/data/OBIC/Freesurfer/estimates/WM_lh_mean.mgz
