%% PAPER - "Association between working memory abilities and brain activity underlying long-term recognition of auditory sequences" (2022)

%Gemma Fernández Rubio (gemmafr@clin.au.dk)

%If you find this script useful, please cite the following papers:
%Fernández-Rubio, G., Carlomagno, F., Vuust, P., Kringelbach, M. L., & Bonetti, L. (2022). Associations between abstract working memory abilities and brain activity underlying long-term recognition of auditory sequences. bioRxiv.
%https://www.biorxiv.org/content/10.1101/2022.05.19.492607v1.abstract
%Bonetti, L., Brattico, E., Carlomagno, F., Cabral, J., Stevner, A., Deco, G., Whybrow, P.C., Pearce, M., Pantazis, D., Vuust, P., & Kringelbach, M.L. (2022). Spatiotemporal whole-brain dynamics of auditory patterns recognition. bioRxiv.
%https://www.biorxiv.org/content/10.1101/2020.06.23.165191v3

%To use this script, you will need to download the following:
%(1) LBPD functions (https://github.com/leonardob92/LBPD-1.0.git) provided by Dr. Leonardo Bonetti (leonardo.bonetti@clin.au.dk)
%(2) FieldTrip (http://www.fieldtriptoolbox.org/)
%(3) SPM (https://www.fil.ion.ucl.ac.uk/spm/)
%(4) OSL toolbox (https://ohba-analysis.github.io/osl-docs/)

%% LBPD_startup_D

%adding some required paths and starting up OSL

pathl = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD'; %path to stored functions
addpath(pathl);
LBPD_startup_D(pathl);

%% PREPROCESSING

%here we follow these steps:
%(1) copying and pasting epoched and continuous files
%(2) band-pass filtering
%(3) (re)epoching
%(4) merging the files

%% Copy epoched and continuous files

%we copy and paste continuous files for all blocks
%we copy and paste epoched files only for block minor

for ii = 10:71 %over subjects
    copyfile(['/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/after_maxfilter_mc/dffspmeeg_SUBJ00' num2str(ii) 'recog*_tsssdsm*'],'/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma'); %continuous files (major and atonal)
    copyfile(['/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/after_maxfilter_mc/e80dffspmeeg_SUBJ00' num2str(ii) 'recogminor_tsssdsm*'],'/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma'); %epoched files (minor)
    disp(ii)
end

%% Band-pass filtering

%we filter the continuous data into delta (slow) frequency band (0.1 - 1 Hz)

addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing'); %path to cluster functions
block = 3; %1 = atonal; 2 = major, 3 = minor

for ii = 10:71 %over subjects
    S3 = []; %empty structure
    if block == 1 %atonal
        S3.D = ['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/dffspmeeg_SUBJ00' num2str(ii) 'recogatonal_tsssdsm.mat']; %atonal continuous data
    elseif block == 2 %major
        S3.D = ['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/dffspmeeg_SUBJ00' num2str(ii) 'recogmajor_tsssdsm.mat']; %major continuous data    
    else %minor
        S3.D = ['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/dffspmeeg_SUBJ00' num2str(ii) 'recogminor_tsssdsm.mat']; %minor continuous data    
    end
    S3.freq = [0.1 1]; %frequency band
    S3.band = 'bandpass';
    S.prefix = 'f_LBPD';
%     D = spm_eeg_filter(S3); %actual function
    jobid = job2cluster(@filtering,S3); %running with parallel computing
    disp(ii)
end

filter_list = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/fdffspmeeg*.mat'); %checking that all the files are there

%% (Re)Epoching

%we redo the epoching using the filtered continuous data
%we use the information from the previously (unfiltered) epoched data and perform the epoching again

addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing'); %path to cluster functions
block = 3; %1 = atonal; 2 = major, 3 = minor

for ii = 1:71 %over subjects
    if block == 1 %previously epoched data (atonal)
        if ii <= 9
            D_e = spm_eeg_load(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/e80dffspmeeg_SUBJ000' num2str(ii) 'recogatonal_tsssdsm.mat']); %loading previously computed epoched data to get epochinfo
        else
            D_e = spm_eeg_load(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/e80dffspmeeg_SUBJ00' num2str(ii) 'recogatonal_tsssdsm.mat']); %loading previously computed epoched data to get epochinfo  
        end
    elseif block == 2 %previously epoched data (major)
        if ii <= 9
            D_e = spm_eeg_load(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/e80dffspmeeg_SUBJ000' num2str(ii) 'recogmajor_tsssdsm.mat']); %loading previously computed epoched data to get epochinfo
        else
            D_e = spm_eeg_load(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/e80dffspmeeg_SUBJ00' num2str(ii) 'recogmajor_tsssdsm.mat']); %loading previously computed epoched data to get epochinfo
        end
    else %previously epoched data (minor)
        if ii <= 9 
            D_e = spm_eeg_load(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/e80dffspmeeg_SUBJ000' num2str(ii) 'recogminor_tsssdsm.mat']); %loading previously computed epoched data to get epochinfo
        else
            D_e = spm_eeg_load(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/e80dffspmeeg_SUBJ00' num2str(ii) 'recogminor_tsssdsm.mat']); %loading previously computed epoched data to get epochinfo
        end
    end
    epochinfo = D_e.epochinfo; %extracting the epochinfo (information about the previously computed epoching that was stored inside the SPM object)
    S3 = []; %empty structure
    if block == 1 %continuous data (atonal)
        if ii <= 9 
            S3.D = ['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/fdffspmeeg_SUBJ000' num2str(ii) 'recogatonal_tsssdsm.mat']; %atonal filtered continuous data
        else
            S3.D = ['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/fdffspmeeg_SUBJ00' num2str(ii) 'recogatonal_tsssdsm.mat']; %atonal filtered continuous data
        end
    elseif block == 2 %continuous data (major)
        if ii <= 9
            S3.D = ['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/fdffspmeeg_SUBJ000' num2str(ii) 'recogmajor_tsssdsm.mat']; %major filtered continuous data
        else
            S3.D = ['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/fdffspmeeg_SUBJ00' num2str(ii) 'recogmajor_tsssdsm.mat']; %major filtered continuous data
        end
    else %continuous data (minor)
        if ii <= 9 
            S3.D = ['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/fdffspmeeg_SUBJ000' num2str(ii) 'recogminor_tsssdsm.mat']; %minor filtered continuous data
        else
            S3.D = ['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/fdffspmeeg_SUBJ00' num2str(ii) 'recogminor_tsssdsm.mat']; %minor filtered continuous data
        end
    end
    if block == 3
        epochinfo.trl(:,2) = epochinfo.trl(:,1) + 525; %the time dimension is different for epoched minor data, so we change it manually
    end
    S3.trl = epochinfo.trl; %trial information from epochinfo
    S3.conditionlabels = epochinfo.conditionlabels; %condition labels
    S3.prefix = 'er'; %epoched redone
%     D = spm_eeg_epochs(S3); %actual function
    jobid = job2cluster(@cluster_epoch,S3); %running with parallel computing
    disp(ii)
end

epoch_filter_list = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/erfdffspmeeg*.mat'); %checking that all the files are there

%% Merging the files

%we combine the atonal, major and minor epoched files into one file per subject

epoch_list = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/erfdffspmeeg*mat'); %epoched files
subj_list = cell(71,3); %store 3 epoched files per subject (atonal, major and minor)

for ii = 3:3:length(epoch_list) %over epoched files
    subj_list(ii/3,:) = {epoch_list(ii-2),epoch_list(ii-1),epoch_list(ii)}; %every 3 files to store atonal, major and minor per subject
end

addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing'); %path to cluster functions

for ii = 65:71%length(subj_list) %over subjects
    S = []; %empty structure
    S.recode = 'same';
    S.prefix = 'm'; %merged
    for ss = 1:3 %over blocks for subject ii
        S.D{ss} = spm_eeg_load([subj_list{ii,ss}.folder '/' subj_list{ii,ss}.name]);
    end
    Dout = spm_eeg_merge(S); %actual function
    disp(ii)
end

merged_list = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/merfdffspmeeg*.mat'); %checking that all the files are there

%% COREGISTRATION

%here we follow these steps:
%(1) creating MRI folders
%(2) copying MRI files
%(3) rhino coregistration
%(4) checking rhino coregistration

%% New MRI folders

%we create one folder for each subject to save the nifti files

for ii = 1:71 %over subjects
    if ii <= 9
        mkdir(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/MRI_nifti_images_all_blocks/000' num2str(ii)])
    else
        mkdir(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/MRI_nifti_images_all_blocks/00' num2str(ii)])
    end
end

%% Moving the MRI nifti files

%we copy the MRI nifti files from the original subjects' folders to the newly created folders

MRI_path = '/projects/MINDLAB2017_MEG-LearningBach/scratch/Leonardo/LearningBach/after_maxfilter_mc/NIFTI_INV2'; %path to subjects' MRI folders
old_folders = dir([MRI_path '/0*']); %subjects' MRI folders
output_path = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/MRI_nifti_images_all_blocks'; %path to subjects' new MRI folders
new_folders = dir([output_path '/0*']); %new subjects' MRI folders

%OBS! make sure that MRI_path and output_path have the same number of folders!

%OBS! if the number of folders is not the same: delete the extra folders from output_path (for subjects that use a template) and manually add them after running the loop

for ii = 1:length(old_folders) %over subjects' MRI folders
    if ii ~= 36 %no nifti file for this subject
    old_name = [old_folders(ii).folder '/' old_folders(ii).name]; %old folder
    file_to_copy = dir([old_name '/*.nii']); %list nifti files
    for jj = 1:length(file_to_copy)
        if ~strcmp(file_to_copy(jj).name(1), 'y') %not the nifti files that start with 'y'
            nifti_file = [file_to_copy(jj).folder '/' file_to_copy(jj).name]; %nifti file we need
        end
    end
    new_name = [new_folders(ii).folder '/' new_folders(ii).name]; %new folder
    if old_name(94:95) == new_name(90:91) %compare the name of the folders to make sure they are the same
        copyfile(nifti_file,new_name) %copy the nifti files into the new folders
    end
    disp(ii)
    end
end

%after this, we manually created new folders for subjects 0007, 0033, 0038 and 0041, and copied the MRI template into them

%% Settings for cluster (parallel computing)

clusterconfig('scheduler', 'cluster'); %if you do not want to submit to the cluster, but simply want to test the script on the hyades computer, you can instead of 'cluster', write 'none'
clusterconfig('long_running', 1); %this is the cue we want to use for the clsuter. There are 3 different cues
clusterconfig('slot', 1); %slot is memory, and 1 memory slot is about 8 GB
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing')

%% RHINO coregistration

%we coregister the data

list = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/merfdffspmeeg*.mat'); %dir to epoched/merged files
a = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/MRI_nifti_images_all_blocks'; %set path to MRI subjects' folders
%OBS! check that all MEG data are in the same order and number as MRI nifti files!

for ii = 7:length(list) %over epoched/merged files
    S = []; %empty structure
    S.ii = ii; %subjects
    S.D = [list(ii).folder '/' list(ii).name]; %path to epoched/merged files
    D = spm_eeg_load(S.D);
    dummyname = D.fname;
    if 7 == exist([a '/' dummyname(19:22)],'dir') %if you have the MRI folder
        dummymri = dir([a '/' dummyname(19:22) '/*.nii']); %path to nifti files (ending with .nii)
        if ~isempty(dummymri)
            S.mri = [dummymri(1).folder '/' dummymri(1).name];
            %standard parameters
            S.useheadshape = 1;
            S.use_rhino = 1; %set 1 for rhino, 0 for no rhino
            S.forward_meg = 'Single Shell';
            S.fid.label.nasion = 'Nasion';
            S.fid.label.lpa = 'LPA';
            S.fid.label.rpa = 'RPA';
            jobid = job2cluster(@coregfunc,S); %running with parallel computing
        else
            warning(['subject ' dummyname(19:22) ' does not have the MRI'])
        end
    end
    disp(ii)
end

%% Checking RHINO coregistration

copy_label = 0; %1 = pasting inv RHINO from epoched data (where it was computed) to continuous data; 0 = showing RHINO coregistration
list = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/merfdffspmeeg*.mat'); %dir to epoched/merged files

for ii = 1:length(list)
    D = spm_eeg_load([list(ii).folder '/' list(ii).name]);
    if isfield(D,'inv')
        if copy_label == 0 %display RHINO coregistration
            if isfield(D,'inv') %checking if the coregistration was already run
                rhino_display(D)
            end
        else %pasting inv RHINO from epoched data (where it was computed) to continuous data
            inv_rhino = D.inv;
            D2 = spm_eeg_load([list(ii).folder '/' list(ii).name(2:end)]);
            D2.inv = inv_rhino;
            D2.save();
        end
    end
    disp(['Subject ' num2str(ii)])
end

%% BEAMFORMING

%here we follow these steps:
%(1) settings for cluster (parallel computing)
%(2) source reconstruction
%(3) statistics over participants

%% Settings for cluster (parallel computing)

clusterconfig('scheduler', 'cluster'); %if you do not want to submit to the cluster, but simply want to test the script on the hyades computer, you can instead of 'cluster', write 'none'
clusterconfig('long_running', 1); %this is the cue we want to use for the clsuter. There are 3 different cues
clusterconfig('slot', 1); %slot is memory, and 1 memory slot is about 8 GB
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing')

%% Source reconstruction

%OBS! check that you don't run this for subjects with no D.inv (RHINO)

%user settings
clust_l = 1; %1 = parallel computing; 0 = running locally
timek = 1:526; %time points
freqq = []; %frequency range (empty [] for broad band)
sensl = 1; %1 = magnetometers only; 2 = gradiometers only; 3 = both magnetometers and gradiometers (OBS! suggested 1!)
workingdir = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/Beam_AllBlocks_Prefiltered_0.1_1'; %high-order working directory (a subfolder for each analysis with information about frequency, time and absolute value will be created)

if isempty(freqq)
    absl = 1; % 1 = absolute value of sources; 0 = not
else
    absl = 0;
end

%actual computation
list = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/merfdffspmeeg*.mat'); %dir to epoched/merged files with RHINO coregistration
condss = {'Old_Correct','New_Correct'}; %conditions
load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/time_normal.mat');
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing');

if ~exist(workingdir,'dir') %creating working folder if it does not exist
    mkdir(workingdir)
end

for ii = 2:length(list) %over subjects
    S = []; %empty structure
    S.Aarhus_cluster = clust_l; %1 = parallel computing; 0 = running locally
    S.norm_megsensors.zscorel_cov = 1; % 1 for zscore normalization; 0 otherwise
    S.norm_megsensors.workdir = workingdir;
    S.norm_megsensors.MEGdata_e = [list(ii).folder '/' list(ii).name];
    S.norm_megsensors.freq = freqq; %frequency range
    S.norm_megsensors.forward = 'Single Shell'; %forward solution
    S.beamfilters.sensl = sensl; %1 = magnetometers; 2 = gradiometers; 3 = both MEG sensors (mag and grad) (SUGGESTED 3!)
    S.beamfilters.maskfname = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_T1_8mm_brain.nii.gz'; %path to brain mask
    S.inversion.znorml = 0; %1 = invert MEG data using the zscored normalized one
    %0 = normalize the original data with respect to maximum and minimum of the experimental conditions if you have both mag and grad
    %0 = use original data in the inversion if you have only mag or grad
    %OBS! suggested 0 in both cases
    S.inversion.timef = timek; %data points to be extracted; leave empty [] for the full length of the epoch
    S.inversion.conditions = condss; %cell with characters for the labels of the experimental conditions
    S.inversion.bc = [1 15]; %extreme time-samples for baseline correction (leave empty [] if you do not want to apply it)
    S.inversion.abs = absl; %1 = absolute values of sources time-series (recommended 1!)
    S.inversion.effects = 1;
    S.smoothing.spatsmootl = 0; %1 = spatial smoothing; 0 = otherwise
    S.smoothing.spat_fwhm = 100; %spatial smoothing fwhm (suggested = 100)
    S.smoothing.tempsmootl = 0; %1 = temporal smoothing; 0 = otherwise
    S.smoothing.temp_param = 0.01; %temporal smoothing parameter (suggested = 0.01)
    S.smoothing.tempplot = [1 2030 3269]; %vector with source indices to be plotted (original vs temporally smoothed timeseries); leave empty [] for no plots
    S.nifti = 1; %1 = plotting nifti images of the reconstructed sources of the experimental conditions
    S.out_name = ['SUBJ_' list(ii).name(19:22)]; %name (character) for output nifti images (conditions name is automatically detected and added)
    if clust_l ~= 1 %useful  mainly for debugging purposes
        MEG_SR_Beam_LBPD(S);
    else
        jobid = job2cluster(@MEG_SR_Beam_LBPD,S); %running with parallel computing
    end
end

%% Statistics over participants

%OBS! parallel computing, specially useful if you have several contrasts

clust = 1; %1 = parallel computing; 0 = running locally
freqq = 1; %1 = 2-8Hz; 2 = broadband (currently only for Block 3)

%building structure
asd = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/Beam*');
S = []; %empty structure
S.workingdir = [asd(freqq).folder '/' asd(freqq).name]; %path where the data from MEG_SR_Beam_LBPD.m is stored
S.sensl = 1; %1 = magnetometers only; 2 = gradiometers only; 3 = both magnetometers and gradiometers
S.plot_nifti = 1; %1 to plot nifti images; 0 otherwise
S.plot_nifti_name = []; %character with name for nifti files (it may be useful if you run separate analysis);lLeave empty [] to not  specify any name
S.contrast = [1 -1; -1 1]; %specify contrasts
S.effects = 1; %mean over subjects for now (but t-tests for contrasts)
if clust == 1
    S.Aarhus_clust = 1; %1 to use paralle computing
    jobid = job2cluster(@MEG_SR_Stats1_Fast_LBPD,S); %running with parallel computing
else
    S.Aarhus_clust = 0; %1 to use paralle computing
    MEG_SR_Stats1_Fast_LBPD(S)
end

%% CONTRAST MEMORIZED VS NOVEL SEQUENCES

%here we follow these steps:
%(1) t-test memorized vs novel sequences
%(2) defining clusters on 3D brain voxel statistics
%(3) combining images with significant clusters
%(4) extracting information about the clusters at source level and reporting it in xlsx files

%% T-test memorized vs novel sequences

%we contrast memorized vs novel sequences for all blocks for each of the time windows corresponding to the tones

subj_old = zeros(23,27,23,526,71); %5D matrix with sources and timepoints for all subjects
subj_new = zeros(23,27,23,526,71); %5D matrix with sources and timepoints for all subjects

for ii = 1:71 %over subjects
    if ii <= 9
        dum_old = load_nii(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/Beam_AllBlocks_Prefiltered_0.1_1/SUBJ_000' num2str(ii) '_cond_Old_Correct_norm0_spatsmoot_0_tempsmoot_0_abs_1.nii.gz']); %load file
        subj_old(:,:,:,:,ii) = dum_old.img; %store image
        dum_new = load_nii(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/Beam_AllBlocks_Prefiltered_0.1_1/SUBJ_000' num2str(ii) '_cond_New_Correct_norm0_spatsmoot_0_tempsmoot_0_abs_1.nii.gz']); %load file
        subj_new(:,:,:,:,ii) = dum_new.img; %store image
    else
        dum_old = load_nii(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/Beam_AllBlocks_Prefiltered_0.1_1/SUBJ_00' num2str(ii) '_cond_Old_Correct_norm0_spatsmoot_0_tempsmoot_0_abs_1.nii.gz']); %load file
        subj_old(:,:,:,:,ii) = dum_old.img; %store image
        dum_new = load_nii(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/Beam_AllBlocks_Prefiltered_0.1_1/SUBJ_00' num2str(ii) '_cond_New_Correct_norm0_spatsmoot_0_tempsmoot_0_abs_1.nii.gz']); %load file
        subj_new(:,:,:,:,ii) = dum_new.img; %store image
    end
    disp(ii)
end

onset = [16 54 91 129 166 203 240 277 314 351 388 425]; %onset of tones (and extra time points) in time samples

for jj = 2:length(onset) %over onset of tones 
    if onset(jj) == 16 || 129 %add 36 or 37 to the onset to obtain the full time window for each tone
        x = 37;
    else
        x = 36;
    end
    
    %time window for memorized sequences
    time_window = subj_old(:,:,:,onset(jj):onset(jj) + x,:); %time window of tone x
    dpm1 = squeeze(mean(time_window,4)); %average of tone x %the fourth dimension (time) is eliminated using 'squeeze' since it corresponds to only one point in time
    
    %time window for memorized sequences
    time_window = subj_new(:,:,:,onset(jj):onset(jj) + x,:); %time window of tone x
    dpm2 = squeeze(mean(time_window,4)); %average of tone x %the fourth dimension (time) is eliminated using 'squeeze' since it corresponds to only one point in time
    
    %t-test
    sz = size(dpm1);
    T = zeros(sz(1),sz(2),sz(3)); %t-value
    P = zeros(sz(1),sz(2),sz(3));%p-value
    for pp = 1:sz(1) %first dimension of voxel coordinates
        for yy = 1:sz(2) %second dimension of voxel coordinates
            for zz = 1:sz(3) %third dimension of voxel coordinates
                if dpm1(pp,yy,zz,1) ~= 0 && dpm2(pp,yy,zz,1) ~= 0 %only for significant brain values
                    a = squeeze(dpm1(pp,yy,zz,:)); %memorized
                    b = squeeze(dpm2(pp,yy,zz,:)); %novel
                    [h,p,ci,stats] = ttest(a,b); %actual t-test for tone x
                    T(pp,yy,zz) = stats.tstat; %t-value
                    P(pp,yy,zz) = 1-p; %p-value
                end
            end
        end
    end
    
    %storing the images
    pathcontrast = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/ContrastAllBlocksPrefiltered';
    dum2 = load_nii('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/FirLev_Block_1_freq_0.1_1.oat/wholebrain_first_level_BC_coape_subj_group_level_everybody_randomise_c1_dir_tone_5/tstat1_gc1_8mm.nii.gz');
    %use dum2 to index a file with only four dimensions (3 voxel dimensions and 1 subject dimension)
    dum2.img = T; %write the image
    save_nii(dum2,[pathcontrast '/Tval_tone_' num2str(jj) '.nii.gz']); %save the nifti image for the t-values
    dum2.img = P; %write the image
    save_nii(dum2,[pathcontrast '/Pval_tone_' num2str(jj) '.nii.gz']); %save the nifti image for the p-values
    disp(jj)
end

clear

%% Defining clusters on 3D brain voxel statistics

%loading p-values and t-values
clear DATA
tone = 1; %set to 1 2 3 4 or 5
tvalbin = 3.5; %value to binarize the data with
T = load_nii(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/ContrastAllBlocksPrefiltered/Tval_tone_' num2str(tone) '.nii.gz']);

%extracting matrix with statistics
T2 = T.img;

%mask for non-0 voxels in brain imges (basically building a layout for actual brain voxels)
mask = zeros(size(T2,1),size(T2,2),size(T2,3));
mask(T2~=0) = 1; %assigning 1 when you have real brain voxels

%memorized vs novel (removing values below t-value = 2)
data = T2;
data(data==0) = NaN; %removing non-brain voxels
data(T2>tvalbin) = 1; %significant brain voxels (positive)
data(T2<-tvalbin) = 1; %significant brain voxels (negative)
data(data~=1) = NaN; %assigning 1 to significant voxels
DATA{1} = data; %storing data

%novel vs memorized (removing values below t-value = 2)
data = T2;
data(data==0) = NaN; %removing non-brain voxels
data(T2>-2) = NaN; %removing non-significant voxels
data(~isnan(data)) = 1; %assigning 1 to significant voxels
DATA{2} = data; %storing data

%getting MNI coordinates
%OBS! you may get a warning since the skeletonized image is not exactly in MNI space, but close enough
[ mni_coords, xform ] = osl_mnimask2mnicoords('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_T1_8mm_brain.nii.gz');

%preparation of information and data for the actual function
for ii = 1%:2 %over directions of the contrast (cond1>cond2 and cond2>cond1)
    S = [];
    S.T = ['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/ContrastAllBlocksPrefiltered/Tval_tone_' num2str(tone) '.nii.gz'];
    S.outdir = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/ContrastAllBlocksPrefiltered/MCS'; %output path
    S.parcelfile = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/aal_8mm_try5.nii.gz';
    S.labels = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/AAL_labels.mat';
    S.MNIcoords = mni_coords; %MNI coordinates of 8mm MNI152T1 brain
    S.data = DATA{ii};
    S.mask = mask; %mask of the brain layout you have your results in
    S.permut = 1000; %number of permutations for Monte Carlo simulation (MCS)
    S.clustmax = 1; %1 = only max cluster size of each permutation MCS (more strict); 0 = every size of each cluster detected for each permutation MCS (less strict)
    S.permthresh = 0.05; %threhsold for MCS
    
    %final names
    if ii == 1
        S.anal_name = ['Tval_tone_' num2str(tone) '_Cond1vsCond2']; %name for the analysis (used to identify and save image and results)
    else
        S.anal_name = ['Tval_tone_' num2str(tone) '_Cond2vsCond1']; %name for the analysis (used to identify and save image and results)
    end
    
    %actual function
    PP = BrainSources_MonteCarlosim_3D_LBPD_D(S);
end

%% Combining images with significant clusters

%we combine images with more than one cluster

tone = 1; %set to 1 2 3 4 or 5
path = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/ContrastAllBlocksPrefiltered/MCS';

%tone x cluster 1
name1 = ['Tval_tone_' num2str(tone) '_Cond1vsCond2_SignClust_1_Tvals.nii.gz'];

%tone x cluster 2
name2 = ['Tval_tone_' num2str(tone) '_Cond1vsCond2_SignClust_2_Tvals.nii.gz'];

%tone x cluster 3
name3 = ['Tval_tone_' num2str(tone) '_Cond1vsCond2_SignClust_3_Tvals.nii.gz'];

output = ['Tval_tone_' num2str(tone) '_All_Clusters.nii.gz'];

cmd = ['fslmaths ' path '/' name1 ' -add ' path '/' name2 ' -add ' path '/' name3 ' ' path  '/' output]; %OBS! be careful with the spacing
system(cmd)

%after this, we created the final images using Workbench

%% Extracting information about the clusters at source level and reporting it in xlsx files

%we obtain information about the brain regions forming the clusters
%the tables can be found in SUPPLEMENTARY MATERIALS

path = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/ContrastAllBlocksPrefiltered/MCS/FinalImages';

for tonebytone = 1:5 %over tones
    fname = [path '/Tval_tone_' num2str(tonebytone) '_All_Clusters.nii.gz']; %tone x cluster 1
    [ mni_coords, xform ] = osl_mnimask2mnicoords(fname); %getting MNI coordinates of significant voxels within the provided image
    V = nii.load(fname); %loading the image
    VV = V(V~=0); %extracting statistics
    VI = find(V~=0); %indices of non-zero values of nifti image
    parcelfile = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/aal_8mm_try5.nii.gz'; %path to AAL template
    load('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/AAL_labels.mat'); %loading AAL labels
    K = nii.load(parcelfile); %extracting AAL coordinates information
    [VV2, II] = sort(VV,'descend'); %sorting results in order to have strongest voxels at the top (positive t-values) or at the bottom (negative t-values)
    VI = VI(II);
    mni_coords = mni_coords(II,:);
    PD = cell(length(VV2),4);  %final cell
    ROI = zeros(length(VI),1); %getting AAL indices
    cnt = 0;
    for ii = 1:length(VI)
        ROI(ii) = K(VI(ii));
        if ROI(ii) > 0 && ROI(ii) < 91
            cnt = cnt + 1;
            PD(cnt,1) = {lab(ROI(ii),3:end)}; %storing ROI
            PD(cnt,4) = {mni_coords(ii,:)}; %storing MNI coordinates
            if mni_coords(ii,1) > 0 %storing hemisphere
                PD(cnt,2) = {'R'};
            else
                PD(cnt,2) = {'L'};
            end
            PD(cnt,3) = {round(VV2(ii),2)}; %storing t-statistics
        end
    end
    PDn = cell2table(PD(~any(cellfun('isempty',PD),2),:)); %remove the possible empty cell
    writetable(PDn,'/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/ContrastAllBlocksPrefiltered/SupplementaryTablesSourceClusters/AllBlocks_Prefiltered_0.1_1.xlsx','Sheet',tonebytone); %printing excel file
end

%% CORRELATION WITH WM SCORES

%here we follow these steps:
%(1) correlation with WM scores
%(2) defining clusters on 3D brain voxel statistics
%(3) combining images with significant clusters
%(4) extracting information about the clusters at source level and reporting it in xlsx files

%% Correlation with WM scores

contrast = 3; %set contrast to 1 (memorized), 2 (novel) or 3 (memorized vs. novel)
list_WM = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/behavioral_data/WM.xlsx'); %dir to working memory excel file
[~,~,WM] = xlsread([list_WM(1).folder '/' list_WM(1).name]); %loads working memory data
SUBJ = zeros(23,27,23,526,71); %5D matrix with sources and timepoints for all subjects

for ii = 1:71 %over subjects
    if ii <= 9
        if contrast == 1 %memorized
            dum2 = load_nii(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/Beam_AllBlocks_Prefiltered_0.1_1/SUBJ_000' num2str(ii) '_cond_Old_Correct_norm0_spatsmoot_0_tempsmoot_0_abs_1.nii.gz']);
            SUBJ(:,:,:,:,ii) = dum2.img; %store image
        elseif contrast == 2 %novel
            dum2 = load_nii(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/Beam_AllBlocks_Prefiltered_0.1_1/SUBJ_000' num2str(ii) '_cond_New_Correct_norm0_spatsmoot_0_tempsmoot_0_abs_1.nii.gz']);
            SUBJ(:,:,:,:,ii) = dum2.img; %store image
        elseif contrast == 3 %memorized vs. novel
            dum_old = load_nii(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/Beam_AllBlocks_Prefiltered_0.1_1/SUBJ_000' num2str(ii) '_cond_Old_Correct_norm0_spatsmoot_0_tempsmoot_0_abs_1.nii.gz']);
            dum_new = load_nii(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/Beam_AllBlocks_Prefiltered_0.1_1/SUBJ_000' num2str(ii) '_cond_New_Correct_norm0_spatsmoot_0_tempsmoot_0_abs_1.nii.gz']);
            dum2.img = dum_old.img - dum_new.img; %contrast memorized and novel images
            SUBJ(:,:,:,:,ii) = dum2.img; %store image
        end
    else
        if contrast == 1 %memorized
            dum2 = load_nii(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/Beam_AllBlocks_Prefiltered_0.1_1/SUBJ_00' num2str(ii) '_cond_Old_Correct_norm0_spatsmoot_0_tempsmoot_0_abs_1.nii.gz']);
            SUBJ(:,:,:,:,ii) = dum2.img; %store image
        elseif contrast == 2 %novel
            dum2 = load_nii(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/Beam_AllBlocks_Prefiltered_0.1_1/SUBJ_00' num2str(ii) '_cond_New_Correct_norm0_spatsmoot_0_tempsmoot_0_abs_1.nii.gz']);
            SUBJ(:,:,:,:,ii) = dum2.img; %store image
        elseif contrast == 3 %memorized vs. novel
            dum_old = load_nii(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/Beam_AllBlocks_Prefiltered_0.1_1/SUBJ_00' num2str(ii) '_cond_Old_Correct_norm0_spatsmoot_0_tempsmoot_0_abs_1.nii.gz']);
            dum_new = load_nii(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/Beam_AllBlocks_Prefiltered_0.1_1/SUBJ_00' num2str(ii) '_cond_New_Correct_norm0_spatsmoot_0_tempsmoot_0_abs_1.nii.gz']);
            dum2.img = dum_old.img - dum_new.img; %contrast memorized and novel images
            SUBJ(:,:,:,:,ii) = dum2.img; %store image
        end
    end
    disp(ii)
end

onset = [16 54 91 129 166 203 240 277 314 351 388 425]; %onset of tones (and extra time points) in time samples

for jj = 1:length(onset) %over onset of tones 
    if onset(jj) == 16 || 129 %add 36 or 37 to the onset to obtain the full time window for each tone
        x = 37;
    else
        x = 36;
    end
    time_window = SUBJ(:,:,:,onset(jj):onset(jj) + x,:); %time window of tone x
    dpm1 = squeeze(mean(time_window,4)); %average of tone x %the fourth dimension (time) is eliminated using 'squeeze' since it corresponds to only one point in time
    
    %correlation
    sz = size(dpm1);
    T = zeros(sz(1),sz(2),sz(3)); %rho-value (correlation coefficient)
    P = zeros(sz(1),sz(2),sz(3)); %p-value
    wm = cell2mat(WM(2:end,2)); %WM scores
    idxnonans = ~isnan(wm); %1 for numbers, 0 for non-values
    for pp = 1:sz(1) %first dimension of voxel coordinates
        for yy = 1:sz(2) %second dimension of voxel coordinates
            for zz = 1:sz(3) %third dimension of voxel coordinates
                if dpm1(pp,yy,zz,1) ~= 0 %only for significant brain values
                    a = squeeze(dpm1(pp,yy,zz,idxnonans));
                    [rho,p] = corr(a,wm(idxnonans)); %correlation for tone x
                    T(pp,yy,zz) = rho; %rho-value
                    P(pp,yy,zz) = 1-p; %p-value
                end
            end
        end
    end
    
    %storing the images
    pathcontrast = ['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/ContrastWM/AllBlocks_Prefiltered_Contrast' num2str(contrast)];
    if ~exist(pathcontrast,'dir')
        mkdir(pathcontrast) %create a new folder for each contrast
    end
    dum2 = load_nii(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/FirLev_Block_1_freq_0.1_1.oat/wholebrain_first_level_BC_coape_subj_group_level_everybody_randomise_c1_dir_tone_5/tstat1_gc1_8mm.nii.gz']);
    %use dum2 to index a file with only four dimensions (3 voxel dimensions and 1 subject dimension)
    dum2.img = T; %write the image
    save_nii(dum2,[pathcontrast '/Rho_WM_tone_' num2str(jj) '.nii.gz']); %save the nifti image for the t-values
    dum2.img = P; %write the image
    save_nii(dum2,[pathcontrast '/Pval_WM_tone_' num2str(jj) '.nii.gz']); %save the nifti image for the p-values
    disp(jj)
end

clear

%% Defining clusters on 3D brain voxel statistics

%loading p-values and t-values
clear DATA
tone = 1; %set to 1 2 3 4 or 5
tvalbin = 0.25; %value to binarize the data with
T = load_nii(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/ContrastWM/AllBlocks_Prefiltered_Contrast3/Rho_WM_tone_' num2str(tone) '.nii.gz']);

%extracting matrix with statistics
T2 = T.img;

%mask for non-0 voxels in brain imges (basically building a layout for actual brain voxels)
mask = zeros(size(T2,1),size(T2,2),size(T2,3));
mask(T2~=0) = 1; %assigning 1 when you have real brain voxels

%preparing data
data = T2;
data(data==0) = NaN; %removing non-brain voxels
data(T2>tvalbin) = 1; %significant brain voxels (positive)
data(T2<-tvalbin) = 1; %significant brain voxels (negative)
data(data~=1) = NaN; %assigning 1 to significant voxels
DATA{1} = data; %storing data

%getting MNI coordinates
%OBS! you may get a warning since the skeletonized image is not exactly in MNI space, but close enough
[ mni_coords, xform ] = osl_mnimask2mnicoords('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_T1_8mm_brain.nii.gz');

%preparation of information and data for the actual function
for ii = 1 %over directions of the contrast (cond1>cond2 and cond2>cond1)
    S = []; %empty structure
    S.T = ['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/ContrastWM/AllBlocks_Prefiltered_Contrast3/Rho_WM_tone_' num2str(tone) '.nii.gz'];
    S.outdir = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/ContrastWM/AllBlocks_Prefiltered_Contrast3/MCS'; %output path
    S.parcelfile = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/aal_8mm_try5.nii.gz';
    S.labels = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/AAL_labels.mat';
    S.MNIcoords = mni_coords; %MNI coordinates of 8mm MNI152T1 brain
    S.data = DATA{ii};
    S.mask = mask; %mask of the brain layout you have your results in
    S.permut = 1000; %number of permutations for Monte Carlo simulation (MCS)
    S.clustmax = 1; %1 = only max cluster size of each permutation MCS (more strict); 0 = every size of each cluster detected for each permutation MCS (less strict)
    S.permthresh = 0.05; %threhsold for MCS
    
    %final names
    if ii == 1
        S.anal_name = ['Rval_tone_' num2str(tone) '_Cond1vsCond2']; %name for the analysis (used to identify and save image and results)
    else
        S.anal_name = ['Rval_tone_' num2str(tone) '_Cond2vsCond1']; %name for the analysis (used to identify and save image and results)
    end
    
    %actual function
    PP = BrainSources_MonteCarlosim_3D_LBPD_D(S);
end

%% Combining images with significant clusters

%we combine images with more than one cluster

tone = 5; %set to 1 2 3 4 or 5

path = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/ContrastWM/AllBlocks_Prefiltered_Contrast3/MCS';

%tone x cluster 1
name1 = ['Rval_tone_' num2str(tone) '_Cond1vsCond2_SignClust_1_Tvals.nii.gz'];

%tone x cluster 2
name2 = ['Rval_tone_' num2str(tone) '_Cond1vsCond2_SignClust_2_Tvals.nii.gz'];

output = ['Rval_tone_' num2str(tone) '_All_Clusters.nii.gz'];

cmd = ['fslmaths ' path '/' name1 ' -add ' path '/' name2 ' ' path  '/' output]; %OBS! be careful with the spacing
system(cmd)

%after this, we created the final images using Workbench

%% Extracting information about the clusters at source level and reporting it in xlsx files

%we obtain information about the brain regions forming the clusters
%the tables can be found in SUPPLEMENTARY MATERIALS

path = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/ContrastWM/AllBlocks_Prefiltered_Contrast3/MCS/FinalImages';

for tonebytone = 1:5
    fname = [path '/Rval_tone_' num2str(tonebytone) '_All_Clusters.nii.gz']; %tone x cluster 1
    [ mni_coords, xform ] = osl_mnimask2mnicoords(fname); %getting MNI coordinates of significant voxels within the provided image
    V = nii.load(fname); %loading the image
    VV = V(V~=0); %extracting statistics
    VI = find(V~=0); %indices of non-zero values of nifti image
    parcelfile = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/aal_8mm_try5.nii.gz'; %path to AAL template
    load('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/AAL_labels.mat'); %loading AAL labels %load this from the provided codes folder
    K = nii.load(parcelfile); %extracting AAL coordinates information
    [VV2, II] = sort(VV,'descend'); %sorting results in order to have strongest voxels at the top (positive t-values) or at the bottom (negative t-values)
    VI = VI(II);
    mni_coords = mni_coords(II,:);
    PD = cell(length(VV2),4);  %final cell
    ROI = zeros(length(VI),1); %getting AAL indices
    cnt = 0;
    for ii = 1:length(VI)
        ROI(ii) = K(VI(ii));
        if ROI(ii) > 0 && ROI(ii) < 91
            cnt = cnt + 1;
            PD(cnt,1) = {lab(ROI(ii),3:end)}; %storing ROI
            PD(cnt,4) = {mni_coords(ii,:)}; %storing MNI coordinates
            if mni_coords(ii,1) > 0 %storing hemisphere
                PD(cnt,2) = {'R'};
            else
                PD(cnt,2) = {'L'};
            end
            PD(cnt,3) = {round(VV2(ii),2)}; %storing t-statistics
        end
    end
    PDn = cell2table(PD(~any(cellfun('isempty',PD),2),:)); %remove the possible empty cell
    writetable(PDn,'/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/ContrastWM/AllBlocks_Prefiltered_Contrast3/SupplementaryTablesSourceClusters/Correlation_WM_AllBlocks_Prefiltered_0.1_1.xlsx','Sheet',tonebytone); %printing excel file
end