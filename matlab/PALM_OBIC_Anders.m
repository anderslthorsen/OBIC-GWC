%% By Anders L.Thorsen and Vilde Brecke, August 2022
% This script runs vertex-wise group comparisons using FSL Palm.
% It is currently set up for running a combined analysis for both hemispheres together.

cd('/data/OBIC/Freesurfer/PALM')

%% Both hemispheres together for GWC
% Runs Palm on both hemispheres at the same time. Consider using palm_hemisplit and palm_hemimerge instead to perform one analysis which might be faster.
palm -i /data/OBIC/Freesurfer/ICA/rh.nu.w-g.avg_fsaverage.mgh ... % 4D surface file for all subjects
-i /data/OBIC/Freesurfer/ICA/lh.nu.w-g.avg_fsaverage.mgh ... % 4D surface file for all subjects
-s /opt/freesurfer/subjects/fsaverage/surf/rh.white /opt/freesurfer/subjects/fsaverage/surf/rh.white.avg.area.mgh ... % Specifies surface and surface area, needed for TFCE
-s /opt/freesurfer/subjects/fsaverage/surf/lh.white /opt/freesurfer/subjects/fsaverage/surf/lh.white.avg.area.mgh ... % Specifies surface and surface area, needed for TFCE
-d /data/OBIC/Freesurfer/PALM/design_OCD_vs_HC_Vest.mat ...  % Design matrix
-t /data/OBIC/Freesurfer/PALM/design_OCD_vs_HC.con ... % Contrast matrix of directional t-tests. Only specify contrasts that you really need, as it takes a long time to perform permutations per contrast
-eb /data/OBIC/Freesurfer/PALM/ExchangeabilityBlocks.csv ... % Ensures that permutations are only done within each site
-T ... % Enables TFCE correction of p-values
-tfce2D ... % Specifices that TFCE is done on surface/2D data
-n 1000 ... % Number of permutations. Should be at least 1000 for the final analysis
-corrmod ... % Bonferroni corrects for the two hemispheres
-corrcon ... % Bonferroni correct for the two one-sided t-tests (OCD>HC and HC>OCD)
-accel tail ... % Enables tail approximation for faster/fewer permutations
-save1-p ... % P-values are saved as 1-p (so that 0.05 is 0.95). Useful for plotting
-o /data/OBIC/Freesurfer/PALM/GWC_2hemispheres_2contrasts

%% Both hemispheres together for gray matter signal intensity
% Runs Palm on both hemispheres at the same time. Consider using palm_hemisplit and palm_hemimerge instead to perform one analysis which might be faster.
palm -i /data/OBIC/Freesurfer/ICA/rh.nu.gm.avg_fsaverage.mgh ... % 4D surface file for all subjects
-i /data/OBIC/Freesurfer/ICA/lh.nu.gm.avg_fsaverage.mgh ... % 4D surface file for all subjects
-s /opt/freesurfer/subjects/fsaverage/surf/rh.white /opt/freesurfer/subjects/fsaverage/surf/rh.white.avg.area.mgh ... % Specifies surface and surface area, needed for TFCE
-s /opt/freesurfer/subjects/fsaverage/surf/lh.white /opt/freesurfer/subjects/fsaverage/surf/lh.white.avg.area.mgh ... % Specifies surface and surface area, needed for TFCE
-d /data/OBIC/Freesurfer/PALM/design_OCD_vs_HC_Vest.mat ...  % Design matrix
-t /data/OBIC/Freesurfer/PALM/design_OCD_vs_HC.con ... % Contrast matrix of directional t-tests. Only specify contrasts that you really need, as it takes a long time to perform permutations per contrast
-eb /data/OBIC/Freesurfer/PALM/ExchangeabilityBlocks.csv ... % Ensures that permutations are only done within each site
-T ... % Enables TFCE correction of p-values
-tfce2D ... % Specifices that TFCE is done on surface/2D data
-n 1000 ... % Number of permutations. Should be at least 1000 for the final analysis
-corrmod ... % Bonferroni corrects for the two hemispheres
-corrcon ... % Bonferroni correct for the two one-sided t-tests (OCD>HC and HC>OCD)
-accel tail ... % Enables tail approximation for faster/fewer permutations
-save1-p ... % P-values are saved as 1-p (so that 0.05 is 0.95). Useful for plotting
-o /data/OBIC/Freesurfer/PALM/GM_2hemispheres_2contrasts

%% Both hemispheres together for white matter signal intensity
% Runs Palm on both hemispheres at the same time. Consider using palm_hemisplit and palm_hemimerge instead to perform one analysis which might be faster.
palm -i /data/OBIC/Freesurfer/ICA/rh.nu.wm.avg_fsaverage.mgh ... % 4D surface file for all subjects
-i /data/OBIC/Freesurfer/ICA/lh.nu.wm.avg_fsaverage.mgh ... % 4D surface file for all subjects
-s /opt/freesurfer/subjects/fsaverage/surf/rh.white /opt/freesurfer/subjects/fsaverage/surf/rh.white.avg.area.mgh ... % Specifies surface and surface area, needed for TFCE
-s /opt/freesurfer/subjects/fsaverage/surf/lh.white /opt/freesurfer/subjects/fsaverage/surf/lh.white.avg.area.mgh ... % Specifies surface and surface area, needed for TFCE
-d /data/OBIC/Freesurfer/PALM/design_OCD_vs_HC_Vest.mat ...  % Design matrix
-t /data/OBIC/Freesurfer/PALM/design_OCD_vs_HC.con ... % Contrast matrix of directional t-tests. Only specify contrasts that you really need, as it takes a long time to perform permutations per contrast
-eb /data/OBIC/Freesurfer/PALM/ExchangeabilityBlocks.csv ... % Ensures that permutations are only done within each site
-T ... % Enables TFCE correction of p-values
-tfce2D ... % Specifices that TFCE is done on surface/2D data
-n 1000 ... % Number of permutations. Should be at least 1000 for the final analysis
-corrmod ... % Bonferroni corrects for the two hemispheres
-corrcon ... % Bonferroni correct for the two one-sided t-tests (OCD>HC and HC>OCD)
-accel tail ... % Enables tail approximation for faster/fewer permutations
-save1-p ... % P-values are saved as 1-p (so that 0.05 is 0.95). Useful for plotting
-o /data/OBIC/Freesurfer/PALM/WM_2hemispheres_2contrasts

%% Both hemispheres together for Regional Vulnearablity Index
% Runs Palm on both hemispheres at the same time. Consider using palm_hemisplit and palm_hemimerge instead to perform one analysis which might be faster.
palm -i /data/OBIC/Freesurfer/ICA/rh.nu.w-g.avg_fsaverage.mgh ... % 4D surface file for all subjects
-i /data/OBIC/Freesurfer/ICA/lh.nu.w-g.avg_fsaverage.mgh ... % 4D surface file for all subjects
-s /opt/freesurfer/subjects/fsaverage/surf/rh.white /opt/freesurfer/subjects/fsaverage/surf/rh.white.avg.area.mgh ... % Specifies surface and surface area, needed for TFCE
-s /opt/freesurfer/subjects/fsaverage/surf/lh.white /opt/freesurfer/subjects/fsaverage/surf/lh.white.avg.area.mgh ... % Specifies surface and surface area, needed for TFCE
-d /data/OBIC/Freesurfer/PALM/design_RVI_Anders_Vest_reanalysis.mat ...  % Design matrix
-t /data/OBIC/Freesurfer/PALM/contrasts_RVI_Anders_Vest.con ... % Contrast matrix of directional t-tests. Only specify contrasts that you really need, as it takes a long time to perform permutations per contrast
-eb /data/OBIC/Freesurfer/PALM/ExchangeabilityBlocks.csv ... % Ensures that permutations are only done within each site
-T ... % Enables TFCE correction of p-values
-tfce2D ... % Specifices that TFCE is done on surface/2D data
-n 1000 ... % Number of permutations. Should be at least 1000 for the final analysis
-corrmod ... % Bonferroni corrects for the two hemispheres
-corrcon ... % Bonferroni correct for the two one-sided t-tests (OCD>HC and HC>OCD)
-accel tail ... % Enables tail approximation for faster/fewer permutations
-save1-p ... % P-values are saved as 1-p (so that 0.05 is 0.95). Useful for plotting
-o /data/OBIC/Freesurfer/PALM/GWC_RVI_2hemispheres_2contrasts_reanalysis

%% Threshold results using Freesurfer's mri_surfcluster
% Thresholds the right hemisphere (m1) and OCD>HC cotrast (c1) which is FWE-corrected for two
% hemispheres (mfwep) at the vertex-wise level (dpv)
!mri_surfcluster --in GWC_2hemispheres_2contrasts_dpv_tstat_mfwep_m1_c1.mgz --subject fsaverage --hemi rh --surf white --annot aparc --thmin 0.95 --thmax 1 --sign pos --no-adjust --sd /data/OBIC/Freesurfer --sum result_mfwep_m1_c1.txt 

% Reanalysis of RVI on 8 June 2023
    % rh
    !mri_surfcluster --in GWC_RVI_2hemispheres_2contrasts_reanalysis_tfce_tstat_mfwep_m1_c1.mgz --subject fsaverage --hemi rh --surf white --annot aparc --thmin 0.95 --thmax 1 --sign pos --no-adjust --sd /data/OBIC/Freesurfer --sum RVI_reanalysis_result_mfwep_m1_c1.txt 
    % lh
    !mri_surfcluster --in GWC_RVI_2hemispheres_2contrasts_reanalysis_tfce_tstat_mfwep_m2_c1.mgz --subject fsaverage --hemi lh --surf white --annot aparc --thmin 0.95 --thmax 1 --sign pos --no-adjust --sd /data/OBIC/Freesurfer --sum RVI_reanalysis_result_mfwep_m2_c1.txt 

  %% Visualize resuts in Freeview
  !freeview -f ../fsaverage/surf/lh.pial:overlay=GWC_RVI_2hemispheres_2contrasts_reanalysis_tfce_tstat_mfwep_m2_c1.mgz -f ../fsaverage/surf/rh.pial:overlay=GWC_RVI_2hemispheres_2contrasts_reanalysis_tfce_tstat_mfwep_m1_c1.mgz