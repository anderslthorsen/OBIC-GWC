%% DA2021
% Modified by Anders Lillevik Thorsen, December 2021
% Changes: Edited fsdata_rh to "rh"
% Changed path to data to path on our server
% Changed line 77, lh to rh
%% Clears variables and console
clc
clear

%% Toolboxes
addpath(genpath('/opt/freesurfer/matlab')); %located under freesurfer installation path / $FREESURFER_HOME
addpath(genpath('toolbox/FastICA_25/')) %https://research.ics.aalto.fi/ica/fastica/code/dlcode.shtml
addpath(genpath('toolbox/icasso122/')) %https://research.ics.aalto.fi/ica/icasso/about+download.shtml
addpath(genpath('toolbox/save2pdf')); %https://www.mathworks.com/matlabcentral/fileexchange/16179-save2pdf

%% Paths/data
cd('/data/OBIC/Freesurfer/ICA') % Path to folder containing this script and freesurfer data
fsdata_lh='lh.nu.w-g.avg_fsaverage.mgh'; % Participant x Vertex input data, left hemisphere (e.g. G/W-ratios)
fsdata_rh='rh.nu.w-g.avg_fsaverage.mgh'; % Participant x Vertex input data, right hemisphere (e.g. G/W-ratios)

%% Read surface data
% Read left hemisphere data
[sL, mL, pL, vL] =load_mgh(fsdata_lh); %LH Vertex data asssigned to sL
size(sL) % = Vertices x 1 x 1 x Participants
surfLL=MRIread(fsdata_lh); %FS-matlab object with header, for later writing of component maps

% Read right hemisphere data
[sR, mR, pR, vR] =load_mgh(fsdata_rh); %RH Vertex data asssigned to sR
surfRR=MRIread(fsdata_rh);

% Transpose (participant x vertices), combine LH and RH, do columnwise z-normalization (across participants, within vertices)
%% Turn Z-score normalization on or off

%mergeLR = [squeeze(sL)', squeeze(sR)']; % When used, this option does not use z-score normalisation
% Til Vilde 24.02.22: Det kan være at klammene [] må taes vekk for at programmet skal kjøre, kan du
% sjekke?
% Fra Vilde: Ingen behov for å fjerne klammene

mergeLR = zscore([squeeze(sL)', squeeze(sR)']); % ORIGINAL CODE: When used, this option uses z-score normalisation 

size(mergeLR) %participants x vertices (lh + rh)

%% Loop through extracting 16-20 components  
for i = 7
    
    %% Run ICA
    Nics=i; % number of components to estimate
    Mtimes=100; %Iterations (for stability/reliability measures)
    visopt='basic'; %Plotting: basic/off
    [ica.Iq,ica.A,ica.W,ica.S,ica.sR]=icasso(mergeLR,Mtimes,'approach','symm','g','tanh','epsilon',1e-10,'maxNumIterations', 1000,'lastEig',Nics,'vis',visopt);
    % OUTPUT:
    %   Iq    (vector) quality index of each estimate
    %   A     (matrix) estimated columns of the mixing matrix (A) = pinv(W)
    %   W     (matrix) estimated rows of the demixing matrix (W)
    %   S     (matrix) estimated independent components
    %   sR    (struct) Icasso result data structure
    
    % Create results folder for chosen model order
    outpath=[outpath_pref,'_d', num2str(Nics)];
    if ~exist(outpath, 'dir')
        mkdir(outpath)
    end
    
    % Save ICA results
    save([outpath filesep 'ICA_Results_d',num2str(Nics),'.mat'],'ica', '-v7.3'); % Save the ICA results
    
    % Save some icasso plots if requested
    if strcmp(visopt,'basic')
        movefile('tmp1.pdf',[outpath filesep 'Clusterindex.pdf']);
        movefile('tmp2.pdf',[outpath filesep 'Stability.pdf']);
        movefile('tmp3.pdf',[outpath filesep 'Dendrogram.pdf']);
        movefile('tmp4.pdf',[outpath filesep 'SimilarityCentrotype.pdf']);
        movefile('tmp5.pdf',[outpath filesep 'Signalplot.pdf']);
    end
    
    %% Save ICA spatial component maps
    size(ica.S) % Components x Vertices, these are the spatial components
    
    % Split hemispheres
    icaSl=ica.S(:,(1:size(ica.S,2)/2));
    icaSr=ica.S(:,(1+size(ica.S,2)/2):end);
    
    % Save ICA spatial maps
    for I=1:size(icaSl,1)
        save_mgh(icaSl(I,:)', [outpath filesep 'lh.IC_' sprintf('%02d',I) '.mgh'], mL, pL);
        save_mgh(icaSr(I,:)', [outpath filesep 'rh.IC_' sprintf('%02d',I) '.mgh'], mR, pR); % 20. December - shouldnt this be rh.IC?
    end
    
end
% %% Do stats on ICA subject weights
% % Prepare data table
% X=readtable('data.csv'); % Read predictor matrix (behavioral/clincal data, participants x variables)
% size(ica.A) % Participants x Components, these are the subject weights
% Y=array2table(ica.A);   % Table with ICA subject weights (participants x components), can be submited to linear models to investigate effects of scanner, age, sex, diagnosis etc.
% grot=strsplit(sprintf('IC%02d\n',1:Nics)); % Make ICA variable labels
% Y.Properties.VariableNames=grot(1:Nics); % Set ICA variable labels
% datatable=[X,Y]; % Merge predictor and ICA-weight tables
%
% % Linear model
% modelspec='IC01 ~ Scanner + Euler + QC + Age + Sex + Diagnosis';
% themodel=fitlm(datatable,modelspec);
% results_f=anova(themodel);
% results_f % Main effects of predictors
% themodel  % Model estimates and t-statistics !! correction for multiple testing when looking at several ICs!!

