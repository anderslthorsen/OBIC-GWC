#!/bin/bash
# Written by Vilde Brecke, Dec 2021
# This script takes the input of a subject ID and loads the relevant Hippocampus and amygdala segmentation (HBT)with brainmask an


outputdir=/data/ENIGMA/amyg_hippo2021/BIDS/Test1/ # Where you want the Freesurfer data to end up
cd ${outputdir}

freeview -v $1/mri/T1.mgz \
$1/mri/brainmask.mgz \
$1/mri/rh.hippoAmygLabels-T1.v21.HBT.mgz:colormap=lut:opacity=0.2 \
$1/mri/lh.hippoAmygLabels-T1.v21.HBT.mgz:colormap=lut:opacity=0.2 \


