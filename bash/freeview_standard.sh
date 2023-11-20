#!/bin/bash
# Written by Anders Lillevik Thorsen, September 2021
# This script takes the input of a subject ID and loads the relevant T1, brainmask, segmentation, pial and white matter surfaces from recon-all in Freesurfer


#outputdir=/data/OBIC/Freesurfer # Where you want the Freesurfer data to end up
#cd ${outputdir}

freeview -v $1/mri/T1.mgz \
$1/mri/wm.mgz \
$1/mri/brainmask.mgz \
$1/mri/aseg.mgz:colormap=lut:opacity=0.2 \
-f $1/surf/lh.white:edgecolor=blue \
$1/surf/lh.pial:edgecolor=red \
$1/surf/rh.white:edgecolor=blue \
$1/surf/rh.pial:edgecolor=red
