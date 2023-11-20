#!/bin/tcsh
#fs
setenv SUBJECTS_DIR /data/OBIC/Freesurfer
set dcm="nu"
set s=$1
set target="fsaverage"
set fwhm="10"

echo $s
echo $SUBJECTS_DIR

cd $SUBJECTS_DIR
cd $s


# extract WM/GM intensity from dcm.mgz
foreach hemi (lh rh)

	# compute points for gray matter - UNSMOOTHED
        mri_vol2surf \
                --mov $SUBJECTS_DIR/${s}/mri/$dcm.mgz \
                --hemi $hemi --noreshape --interp trilinear \
                --projfrac-avg 0.1 0.6 0.1 \
                --o $SUBJECTS_DIR/$s/surf/${hemi}.${dcm}.avg.gm.mgh --regheader $s

	# compute points for white matter - UNSMOOTHED
        mri_vol2surf \
                --mov $SUBJECTS_DIR/$s/mri/$dcm.mgz \
                --hemi $hemi --noreshape --interp trilinear \
                --projdist-avg -1.5 -0.15 0.15 \
                --o $SUBJECTS_DIR/$s/surf/${hemi}.${dcm}.avg.wm.mgh --regheader $s


        # computing the contrast - UNSMOOTHED

        mri_concat \
                --i $SUBJECTS_DIR/${s}/surf/${hemi}.${dcm}.avg.wm.mgh \
                --i $SUBJECTS_DIR/${s}/surf/${hemi}.${dcm}.avg.gm.mgh \
                --paired-diff-norm --mul 100 \
               --o $SUBJECTS_DIR/${s}/surf/${hemi}.${dcm}.w-g.avg.mgh


	# compute points for gray matter - SMOOTHED

        mri_vol2surf \

                --mov $SUBJECTS_DIR/${s}/mri/$dcm.mgz \
                --hemi $hemi --noreshape --interp trilinear \
                --projfrac-avg 0.1 0.6 0.1 \
		--fwhm $fwhm \ 
                --o $SUBJECTS_DIR/$s/surf/${hemi}.${dcm}.avg.gm.smoothed.mgh --regheader $s

	# compute points for white matter - SMOOTHED

        mri_vol2surf \
                --mov $SUBJECTS_DIR/$s/mri/$dcm.mgz \
                --hemi $hemi --noreshape --interp trilinear \
                --projdist-avg -1.5 -0.15 0.15 \
		--fwhm $fwhm \ 
                --o $SUBJECTS_DIR/$s/surf/${hemi}.${dcm}.avg.wm.smoothed.mgh --regheader $s


        # computing the contrast - SMOOTHED
        mri_concat \
                --i $SUBJECTS_DIR/${s}/surf/${hemi}.${dcm}.avg.wm.smoothed.mgh \
                --i $SUBJECTS_DIR/${s}/surf/${hemi}.${dcm}.avg.gm.smoothed.mgh \
                --paired-diff-norm --mul 100 \
               --o $SUBJECTS_DIR/${s}/surf/${hemi}.${dcm}.w-g.avg.smoothed.mgh




end


# surf -> fsaverage 
foreach meas (w-g.avg avg.wm avg.gm avg.wm.smoothed avg.gm.smoothed w-g.avg.smoothed) # Now also includes seperately smoothed wm, gm and wm-gm contrast


        foreach hemi (lh rh)

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




        end


end


#smooth
foreach meas (w-g.avg) 


        foreach hemi (lh rh)


                mri_surf2surf \
                        --hemi $hemi \
                        --s $target \
                        --sval $SUBJECTS_DIR/${s}/surf/${hemi}.${dcm}.${meas}_fsaverage.mgh \
                        --fwhm $fwhm \
                        --cortex \
                        --tval $SUBJECTS_DIR/${s}/surf/${hemi}.${dcm}.${meas}_fsaverage.${fwhm}.mgh \


        end
end
