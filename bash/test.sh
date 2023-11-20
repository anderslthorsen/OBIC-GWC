#!/bin/bash
# Written by Vilde
# This script takes the input of a subject ID and loads the relevant T1.mgz, skullstripp output (brainmask.mgz), white matter volume parcellation (wm.mgz), volume-based anatomical parcellation (aparc + aseg), and the intensity-normalized image that is input to surfaces generation (brain.finalsurfs.mgz) and white (blue outlines) and pial (red outlines) furface renderings for both hemispheres, based on Savilia et al 2015


outputdir=/data/OBIC/Freesurfer # Where you want the Freesurfer data to end up
cd ${outputdir}

freeview -v $1/mri/T1.mgz \
$1/mri/wm.mgz \
$1/mri/brainmask.mgz \
$1/mri/finalsurfs-mgz \
$1/mri/aparc+aseg.mgz:colormap=lut:opacity=0.75 \
-f $1/surf/lh.white:edgecolor=blue \
$1/surf/lh.pial:edgecolor=red \
$1/surf/rh.white:edgecolor=blue \
$1/surf/rh.pial:edgecolor=red
