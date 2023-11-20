%% By Anders L.Thorsen and Vilde Brecke, August 2022
% This script runs vertex-wise group comparisons using FSL Palm.
% It is currently set up for running a combined analysis for both hemispheres together.

cd('/data/OBIC/Freesurfer/PALM/symptom_dims/')

%% Both hemispheres together for GWC
% Runs Palm on both hemispheres at the same time. Consider using palm_hemisplit and palm_hemimerge instead to perform one analysis which might be faster.
palm -i symptom_dims_rh.nu.w-g.avg_fsaverage10.mgh ... % 4D surface file for all subjects
-i symptom_dims_lh.nu.w-g.avg_fsaverage10.mgh ... % 4D surface file for all subjects
-s /opt/freesurfer/subjects/fsaverage/surf/rh.white /opt/freesurfer/subjects/fsaverage/surf/rh.white.avg.area.mgh ... % Specifies surface and surface area, needed for TFCE
-s /opt/freesurfer/subjects/fsaverage/surf/lh.white /opt/freesurfer/subjects/fsaverage/surf/lh.white.avg.area.mgh ... % Specifies surface and surface area, needed for TFCE
-d symptom_dim_design.mat ...  % Design matrix
-t symptom_dim_contrasts.con ... % Contrast matrix of directional t-tests. Only specify contrasts that you really need, as it takes a long time to perform permutations per contrast
-eb symptom_dim_exchangeability_block.csv ... % Ensures that permutations are only done within each site
-T ... % Enables TFCE correction of p-values
-tfce2D ... % Specifices that TFCE is done on surface/2D data
-n 1000 ... % Number of permutations. Should be at least 1000 for the final analysis
-corrmod ... % Bonferroni corrects for the two hemispheres
-corrcon ... % Bonferroni correct for the two one-sided t-tests (positive and negative association with sex/rel symptoms)
-accel tail ... % Enables tail approximation for faster/fewer permutations
-save1-p ... % P-values are saved as 1-p (so that 0.05 is 0.95). Useful for plotting
-o /data/OBIC/Freesurfer/PALM/symptom_dims/GWC_symptomdims_2hemispheres_2contrasts

%% Threshold results using Freesurfer's mri_surfcluster
% Thresholds the right hemisphere (m1) and OCD>HC cotrast (c1) which is FWE-corrected for two
% hemispheres (mfwep) at the vertex-wise level (dpv)
%!mri_surfcluster --in GWC_2hemispheres_2contrasts_dpv_tstat_mfwep_m1_c1.mgz --subject fsaverage --hemi rh --surf white --annot aparc --thmin 0.95 --thmax 1 --sign pos --no-adjust --sd /data/OBIC/Freesurfer --sum result_mfwep_m1_c1.txt 