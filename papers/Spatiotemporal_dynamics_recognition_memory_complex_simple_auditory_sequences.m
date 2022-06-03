%% PAPER - "The spatiotemporal dynamics of recognition memory for complex versus simple auditory sequences" (2022)

%Gemma Fernández Rubio (gemmafr@clin.au.dk)

%If you find this script useful, please cite the following papers:
%Fernández Rubio, G., Brattico, E., Kotz, S. A., Kringelbach, M. L., Vuust, P., & Bonetti,  L. (2022). The spatiotemporal dynamics of recognition memory for complex versus simple auditory sequences. bioRxiv.
%https://www.biorxiv.org/content/10.1101/2022.05.15.492038v1
%Bonetti, L., Brattico, E., Carlomagno, F., Cabral, J., Stevner, A., Deco, G., Whybrow, P.C., Pearce, M., Pantazis, D., Vuust, P., & Kringelbach, M.L. (2022). Spatiotemporal whole-brain dynamics of auditory patterns recognition. bioRxiv.
%https://www.biorxiv.org/content/10.1101/2020.06.23.165191v3

%To use this script, you will need to download the following:
%(1) LBPD functions (https://github.com/leonardob92/LBPD-1.0.git) provided by Dr. Leonardo Bonetti (leonardo.bonetti@clin.au.dk)
%(2) FieldTrip (http://www.fieldtriptoolbox.org/)
%(3) SPM (https://www.fil.ion.ucl.ac.uk/spm/)
%(4) OSL toolbox (https://ohba-analysis.github.io/osl-docs/)

%% Settings for cluster (parallel computing)

addpath('/projects/MINDLAB2020_MEG-AuditoryPatternRecognition/scripts/leonardo')
clusterconfig('slot', 2); %set manually the job cluster slots
clusterconfig('scheduler', 'cluster');
clusterconfig('long_running', 1); %0 = short queue, 1 = all queue, 2 = long queue

%% LBPD_startup_D

pathl = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD'; %path to stored functions
addpath(pathl);
LBPD_startup_D(pathl);
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/scripts_osl_learningbach/Cluster') %path to the function that submits jobs to the server

%% PREPROCESSING

%here we follow these steps:
%(1) maxfilter
%(2) converting .fif files into SPM objects
%(3) creating a list of recog files
%(4) filtering and downsampling
%(5) epoching
%(6) defining the conditions
%(7) checking answers
%(8) averaging

%% Maxfilter

%OBS! before running maxfilter: (1) close matlab, (2) open terminal and write 'use anaconda', (3) open matlab, (4) run maxfilter

addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/osl/osl-core'); %path to OSL functions
osl_startup %starting up OSL
basefifdatadir = '/raw/sorted/MINDLAB2017_MEG-LearningBach';
basefifdatadirraw = '/raw/sorted/MINDLAB2017_MEG-LearningBach';
datadir = '/projects/MINDLAB2017_MEG-LearningBach/scratch/Leonardo/LearningBach/maxfilter_preproc/';
workingdir = datadir;
subnum = '1'; %set the sub ID for starting from a specific chosen subject in the temp1fif list
temp1fif = dir([basefifdatadir '/0*']);

for ii = str2double(subnum):length(temp1fif) %this iterates over subjects (remember to set only the subjects that you want)
    temp2fif = [basefifdatadirraw '/' temp1fif(ii).name];
    dummy3fif = dir([temp2fif '/20*']);
    dummy2_3fif = dir([temp2fif '/' dummy3fif(1).name '/MEG/0*']);
    if isempty(dummy2_3fif)
        dummy2_3fif = dir([temp2fif '/' dummy3fif(2).name '/MEG/0*']);
    end
    if isempty(dummy2_3fif)
        display('problem in getting the raw MEG data!');
    end
    fif_files{ii,1} = dummy2_3fif;
    for k = 1:length(dummy2_3fif) %this iterates over files for each subject
        temp4fif = [dummy2_3fif(k).folder '/' dummy2_3fif(k).name '/files/' dummy2_3fif(k).name(5:end) '.fif'];
        a{ii,k} = temp4fif; %cell for storing all the paths
        S = [];
        S.dataset = temp4fif;
        S.outfile = ['spmeeg_SUBJ' temp1fif(ii).name dummy2_3fif(k).name(5:end)];
    end
end

%new proper lines for maxfilter
maxfilter_path = '/neuro/bin/util/maxfilter';
maxDir = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach_Fsample1000Hz'; %with movement compensation
project = 'MINDLAB2017_MEG-LearningBach';
for ii = str2double(subnum):length(temp1fif) %iterates over subjects
    fif_files_dummy = fif_files{ii,1};
    spm_files_dummy = fif_files{ii,1};
    for j = 1:length(spm_files_dummy) %iterates over experimental blocks
        if ~isempty(strfind(spm_files_dummy(j).name,'rest10')) || ~isempty(strfind(spm_files_dummy(j).name,'learminor')) || ~isempty(strfind(spm_files_dummy(j).name,'recogminor'))  %here for now I want only learminor, recogminor and rest10
            rawName{j} = [fif_files_dummy(j).folder '/' fif_files_dummy(j).name '/files' '/' fif_files_dummy(j).name(5:end) '.fif'];
            [~,n,~] = fileparts(rawName{j});
            maxfName = ['SUBJ00' num2str(ii) rawName{j}(regexp(rawName{j},'files/')+6:regexp(rawName{j},'.fif')-1)];
            badchans = [];
            %commands for submitting the job to the cluster
            cmd = ['submit_to_cluster -q maxfilter.q -n 4 -p ',num2str(project),' "',maxfilter_path,' -f ',[fif_files_dummy(j).folder '/' fif_files_dummy(j).name '/files' '/' fif_files_dummy(j).name(5:end) '.fif'],' -o ' [maxDir '/' maxfName '_tsssdsm.fif'] ' -st 4 -corr 0.98 -movecomp -ds 3 ',' -format float -v | tee ' [maxDir '/log_files/' maxfName '_tsssdsm.log"']];
            system(cmd);
        end
    end
end


%% Converting .fif files into SPM objects and creating some file names for later steps

%setting directories
clear
datadir = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/after_maxfilter_mc'; %with movement compensation in maxfilter
workingdir = datadir;
basefifdatadir = datadir;

%converting .fif files into SPM objects
temp1fif = dir([basefifdatadir '/*.fif']);
for ii = 1:length(temp1fif) %this iterates over blocks through all subjects (remember to set only the subjects that you want)
    temp2fif = [basefifdatadir '/' temp1fif(ii).name];
    fif_files{ii,1} = temp2fif;
    S = [];
    S.dataset = fif_files{ii,1};
    D = spm_eeg_convert(S); 
end

%create the spm object files path and name
for ii = 1:length(temp1fif)
    spm_files_basenames{ii,1} = ['spmeeg_' temp1fif(ii).name(1:end-4) '.mat'];
    spm_files{ii,1} = [datadir '/' 'spmeeg_' temp1fif(ii).name(1:end-4) '.mat'];
end

%% Creating a list of recog files only (SPM objects)

k = 0;
kk = 0;
jj = 0;
clear spm_files_lear_basen spm_files_recog_basen spm_files_lear spm_files_recog xlsx_basenames

for ii = 1:length(temp1fif)
    if strcmp(temp1fif(ii).name(9:12),'lear')
        k = k + 1;
        %spm_files_recog{k,1} = spm_files{ii,1};
        spm_files_lear_basen{k,1} = spm_files_basenames{ii,1};
        spm_files_lear{k,1}=[workingdir '/dff' spm_files_lear_basen{k}];
%         xlsx_basenames{k,1} = [spm_files_lear_basen{k,1}(12:15) spm_files_lear_basen{k,1}(21:23) '.xlsx'];
    end
    if strcmp(temp1fif(ii).name(9:12),'reco')
        kk = kk + 1;
        %spm_files_recog{k,1} = spm_files{ii,1};
        spm_files_recog_basen{kk,1} = spm_files_basenames{ii,1};
        spm_files_recog{kk,1}=[workingdir '/dff' spm_files_recog_basen{kk}];
        xlsx_basenames{kk,1} = [spm_files_recog_basen{kk,1}(12:15) spm_files_recog_basen{kk,1}(21:23) '.xlsx'];
    end
    if strcmp(temp1fif(ii).name(9:12),'rest')
        jj = jj + 1;
        spm_files_rest10_basen{jj,1} = spm_files_basenames{ii,1};
        spm_files_rest10{jj,1}=[workingdir '/dff' spm_files_rest10_basen{jj}];
    end
end

%% Filtering and downsampling

%high-pass filter
%preffix for these files: 'f'
for ii = 1:length(spm_files) % iterates over blocks !OBS!
    spm_files{ii} = [workingdir '/' spm_files_basenames{ii}];
    S2 = [];
    S2.D = spm_files{ii};
    S2.band = 'high';
    S2.freq = 0.1;
    D = spm_eeg_filter(S2);
end

%notch filter
%preffix for these files: 'ff'
for ii = 1:length(spm_files)
    spm_files{ii} = [workingdir '/f' spm_files_basenames{ii}];
    S2 = [];
    S2.D = [spm_files{ii}];
    S2.band = 'stop';
    S2.freq = [48 52]; %this defines the notch filter frequency range (i.e., around 50Hz)
    D = spm_eeg_filter(S2);
end

%donwsampling
for ii = 1:length(spm_files) % iterates over experimental blocks
    spm_files{ii} = [workingdir '/ff' spm_files_basenames{ii}];
    S = [];
    S.D = spm_files{ii};
    S.fsample_new = 150; % in Hz
    D = spm_eeg_downsample(S);    
end

%removing bad trials (if any) from the data by performing a visual inspection
for ii = 1:length(spm_fils)
    spm_files{ii}=[workingdir '/dff' spm_files_basenames{ii}];
    D = spm_eeg_load(spm_files{ii});
    D = oslview(D);
    D.save(); %this is needed to save the marked bad events (and/or channels) in OSLview
end

%removing eyeblink and heartbeat artifacts by computing ICA
for ii = 1:length(spm_files)
    spm_files{ii}=[workingdir '/dff' spm_files_basenames{ii}];
    D = spm_eeg_load(spm_files{ii});
    D = osl_africa(D,'do_ident',false,'do_remove',false,'used_maxfilter',true); %independent components
    D = osl_africa(D,'do_remove',false); 
    D = osl_africa(D,'do_ident',false,'do_remove',true);
    D.save();
    disp(['just done subject number ' num2str(ii)])
end

%% Epoching

%one epoch per old/new excerpt
%baseline = (-)100ms

prefix_tobeadded = 'e80';

%creates a list of only the recog files (subgroup of the whole set of data files)
k = 0;
for ii = 1:length(temp1fif)
    if strcmp(temp1fif(ii).name(9:12),'reco')
        k = k + 1;
        %spm_files_recog{k,1} = spm_files{ii,1};
        spm_files_recog_basen{k,1} = spm_files_basenames{ii,1};
        spm_files_recog{k,1}=[workingdir '/dff' spm_files_recog_basen{k}];
        xlsx_basenames{k,1} = [spm_files_recog_basen{k,1}(12:15) spm_files_recog_basen{k,1}(21:23) '.xlsx'];
    end
end

%iterates for the recognition files only
for ii = 3:3:length(spm_files_recog_basen) %indexing only the recognition files
    D = spm_eeg_load(spm_files_recog{ii,1}); %load the spm object
    events = D.events;
    count_evval = 0;
    for ieve = 1:length(events)
        if strcmp(events(ieve).type,'STI101_up')
            if events(ieve).value == 5
                count_evval = count_evval + 1;
                trigcor(count_evval,1) = events(ieve).time + 0.010; %this takes the correct triggers and adds 10ms of delay of the sound travelling into the tubes
            end
        end
    end
    if count_evval ~= 4 && count_evval ~= 80
        disp('warning.. there is something wrong with the triggers');
        disp(spm_files_reco_basen{ii}(8:end-12));
    end
    trl_sam = zeros(length(trigcor),3);
    trl_sec = zeros(length(trigcor),3);
    deftrig = zeros(length(trigcor),1);
    for k = 1:length(trigcor)
        deftrig(k,1) = 0.012 + trigcor(k,1); %adding a 0.012 seconds delay to the triggers sent during the experiment (this delay was due to technical reasons related to the stimuli)
        trl_sec(k,1) = deftrig(k,1) - 0.1000; %beginning time-window epoch in s (please note that we computed this operation two times, obtaining two slightly different pre-stimulus times
        %this was done because for some computations was convenient to have a slightly longer pre-stimulus time
        trl_sec(k,2) = deftrig(k,1) + 6.000; %end time-window epoch in s
        trl_sec(k,3) = trl_sec(k,2) - trl_sec(k,1); %range time-windows in s
        trl_sam(k,1) = round(trl_sec(k,1) * 150); %beginning time-window epoch in samples
        trl_sam(k,2) = round(trl_sec(k,2) * 150); %end time-window epoch in samples
        trl_sam(k,3) = -15; %sample before the onset of the stimulus (corresponds to 0.100ms)
    end
    %create the epochinfo structure that is required for the source reconstruction later
    epochinfo.trl = trl_sam;
    epochinfo.time_continuous = D.time;
    %switch the montage to 0 (OSL does epoching with the not denoised data)
    D = D.montage('switch',0);
    %build structure for spm_eeg_epochs
    S = [];
    S.D = D;
    S.trl = trl_sam;
    S.prefix = prefix_tobeadded;
    D = spm_eeg_epochs(S);
    %store the epochinfo structure inside the D object
    D.epochinfo = epochinfo;
    D.save();
    %take bad segments registered in OSLview and check if they overlap with the trials. If so, it gives the number of overlapped trials that will be removed later   
    count = 0;
    Bad_trials = zeros(length(trigcor),1);
    for kkk = 1:length(events) %over events
        if strcmp(events(kkk).type,'artefact_OSL')
            for k = 1:length(trl_sec) %over trials
                if events(kkk).time - trl_sec(k,2) < 0 %if end of trial is > than beginning of artifact
                    if trl_sec(k,1) < (events(kkk).time + events(kkk).duration) %if beginning of trial is < than end of artifact
                        Bad_trials(k,1) = 1; %it is a bad trial (stored here)
                        count = count + 1;
                    end
                end                  
            end
        end
    end
    %if bad trials were detected, their indices are stored within D.badtrials field
    if count == 0
        disp('there are no bad trials marked in oslview for');
        disp(spm_files_recog_basen(ii));
    else
        D = badtrials(D,find(Bad_trials),1); %get the indices of the badtrials marked as '1' (that means bad)
        D.save(); %saving on disk
        disp('bad trials are ')
        D.badtrials
    end
end


%% Defining the conditions

%here we label the conditions for the atonal and tonal trials

atonal_label = 1; % 1 = atonal; 0 = tonal

%define conditions - only 1 epoch for each memorized/novel sequence (baseline = (-)100ms)
basnam_con_oc = {'Old_Correct'}; %memorized correct
basnam_con_nc = {'New_Correct'}; %novel correct
basnam_con_ouc = {'Old_Incorrect'}; %memorized incorrect
basnam_con_nuc = {'New_Incorrect'}; %novel incorrect

%dir to MEG behavioral results and epoched files
listXLSX_atonal = dir('/scratch7/MINDLAB2017_MEG-LearningBach/BehavioralResponses/*ato.xlsx'); %dir to MEG behavioral results (atonal)
listXLSX_major = dir('/scratch7/MINDLAB2017_MEG-LearningBach/BehavioralResponses/*maj.xlsx'); %dir to MEG behavioral results (tonal)
listMEG_atonal = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/e80*atonal_tsssdsm.mat'); %dir to epoched files (atonal)
listMEG_major = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/e80*major_tsssdsm.mat'); %dir to epoched files (tonal)

for ii = 1:length(listMEG_atonal) %over epoched files
    if atonal_label == 0 %OBS! remember to set this to atonal (1) or tonal (0)
        D = spm_eeg_load([listMEG_atonal(ii).folder '/' listMEG_atonal(ii).name]); %loads epoched files
        [~,~,raw_recog] = xlsread([listXLSX_atonal(ii).folder '/' listXLSX_atonal(ii).name]); %loads behavioral results
    else
        D = spm_eeg_load([listMEG_major(ii).folder '/' listMEG_major(ii).name]); %loads epoched files
        [~,~,raw_recog] = xlsread([listXLSX_major(ii).folder '/' listXLSX_major(ii).name]); %loads behavioral results
    end
    for k = 1:length(D.trialonset) %over onset of each trial
        if strcmp(raw_recog{(k + 1),3},'No_response') %if there was no response
            D = D.conditions(k,'No_response'); %assign 'no response'
        elseif raw_recog{(k + 1),2}(1) == raw_recog{(k + 1),3}(1) %if the response was correct
            if raw_recog{(k + 1),2}(1) == 'O' %if the trial was memorized
                D = D.conditions(k,basnam_con_oc); %assign 'memorized correct'
            else
                D = D.conditions(k,basnam_con_nc); %otherwise assign 'novel correct'
            end
        else %otherwise the response was incorrect
            if raw_recog{(k + 1),2}(1) == 'O' %if the trial was memorized
                D = D.conditions(k,basnam_con_ouc); %assign 'memorized incorrect'
            else
                D = D.conditions(k,basnam_con_nuc); %otherwise assign 'novel incorrect'
            end
        end
    end
    if ~isempty(D.badtrials) %overwriting badtrials (if any) on condition labels
        BadTrials = D.badtrials;
        for badcount = 1:length(BadTrials)
            D = D.conditions(BadTrials(badcount),'Bad_trial');
        end
    end
    D = D.montage('switch',1); %AFRICA denoised montage
    D.epochinfo.conditionlabels = D.conditions; %added for later use in the source reconstruction
    D.save(); %saving data on disk
    disp(num2str(ii))
end

%%  Checking answers

%this counts how many memorized correct and novel correct answers we have for each condition (atonal and tonal) per subject

listMEG_atonal = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/e80*atonal_tsssdsm.mat'); %dir to epoched files (atonal)
listMEG_major = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/e80*major_tsssdsm.mat'); %dir to epoched files (tonal)

%atonal
res_ato = zeros(length(listMEG_atonal),2);
for ii = 1:length(listMEG_atonal) %over epoched files
    D = spm_eeg_load([listMEG_atonal(ii).folder '/' listMEG_atonal(ii).name]); %loads epoched files
    cntOC = 0; %memorized correct %starts the counter from 0
    cntNC = 0; %novel correct %starts the counter from 0
    for jj = 1:length(D.conditions) %over conditions (memorized correct, novel correct, memorized incorrect, novel incorrect)
       if strcmp('Old_Correct',D.conditions{jj})
           cntOC = cntOC + 1;
       end
       if strcmp('New_Correct',D.conditions{jj})
           cntNC = cntNC + 1;
       end
    end
    res_ato(ii,1) = cntOC; %total number of memorized correct per subject
    res_ato(ii,2) = cntNC; %total number of novel correct per subject
    disp(ii)
end

%tonal
res_maj = zeros(length(listMEG_major),2);
for ii = 1:length(listMEG_major) %over epoched files
    D = spm_eeg_load([listMEG_major(ii).folder '/' listMEG_major(ii).name]);
    cntOC = 0; %memorized correct %starts the counter from 0
    cntNC = 0; %novel correct %starts the counter from 0
    for jj = 1:length(D.conditions) %over conditions (memorized correct, novel correct, memorized incorrect, novel incorrect)
       if strcmp('Old_Correct',D.conditions{jj})
           cntOC = cntOC + 1;
       end
       if strcmp('New_Correct',D.conditions{jj})
           cntNC = cntNC + 1;
       end
    end
    res_maj(ii,1) = cntOC; %total number of memorized correct per subject
    res_maj(ii,2) = cntNC; %total number of novel correct per subject
    disp(ii)
end

%% Averaging

%settings for cluster (parallel computing)
clusterconfig('scheduler', 'cluster'); %set automatically the long run queue
clusterconfig('long_running', 1); %set automatically the long run queue
clusterconfig('slot', 1); %set manually the job cluster slots
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing')

%averaging
output_prefix_to_be_set = 'm'; %after doing the averaging
epoch_list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/e80*.mat'); %dir to epoched files (atonal and tonal)
for ii = 1:length(epoch_list) %over epoched files
    input = [];
    input.D = [epoch_list(ii).folder '/' epoch_list(ii).name];
    input.prefix = output_prefix_to_be_set;
    jobid = job2cluster(@sensor_average, input);
end
average_list = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/me80*.mat');

%combining planar gradiometers
for ii = 1:length(average_list) %over averaged files
    input = [];
    input.D = [average_list(ii).folder '/' average_list(ii).name];
    D = spm_eeg_load(input.D);
    D = D.montage('switch',1); %AFRICA denoised data
    D.save();
    jobid = job2cluster(@combining_planar_cluster, input);
end

%% MEG SENSOR DATA ANALYSES

%here we follow these steps:
%(1) extracting MEG sensor data
%(2) loading t-test results and reshaping them
%(3) getting information about the timepoints and plotting

%% Extracting MEG sensor data

block = 1; %1 = atonal; 2 = tonal
load_data = 1; %set to 1 if you want to load the data instead of extracting it from SPM objects
subjects = 1:71;
outdir = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/MEG_sensors'; %output path

if block == 1 %atonal
    list = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/Pme80*recogatonal_tsssdsm.mat');
else %tonal
    list = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/Pme80*recogmajor_tsssdsm.mat');
end

v = 1:71; %subjects
S = [];
S.data = [];
if load_data == 1 %if you already computed and saved the t-tests you can load them here
    if block == 1 %atonal
        load([outdir '/recogatonal_all_conditions.mat']);
    else %tonal
        load([outdir '/recogmajor_all_conditions.mat']);
    end
    S.data = data_mat(:,:,subjects,:);
    S.chanlabels = chanlabels;
    S.time_real = time_sel;
else %otherwise you can extract the data from SPM MEG objects (one for each subject)
    S.spm_list = cell(1,length(v));
    for ii = 1:length(v)
        S.spm_list(ii) = {[list(v(ii)).folder '/' list(v(ii)).name]};
    end
end
S.conditions = {'Old_Correct','New_Correct'}; %set the conditions
S.timeextract = []; %time-points to be extracted
S.centerdata0 = 0; %1 to make data starting at 0
S.save_data = 1; %only meaningful if you read data from SPM objects saved on disk
if block == 1
    S.save_name_data = 'recogatonal_all_conditions';
else
    S.save_name_data = 'recogmajor_all_conditions';
end

%individual waveform plotting
S.waveform_singlechannels_label = 0; %1 = plot single channel waveforms
S.wave_plot_conditions_together = 0; %1 = plot the average of all
S.mag_lab = 1; %1 = magnetometers; 2 = gradiometers
S.x_lim_temp_wave = []; %limits for time (in secs) (e.g. [-0.1 3.4])
S.y_lim_ampl_wave = []; %limit for amplitude (e.g. [0 120] magnetometes, [0 6] gradiometers)

%averaged waveform plotting
S.waveform_average_label = 0; %average of some channels
S.legc = 0; %1 = legend
S.left_mag = 13; %13 %37 (visual) %43 (visual) %199 (visual) %203 (visual) %channels for averaging
S.signtp = {[]};
% S.sr = 150; %sampling rate (Hz)
S.avewave_contrast = 0; %1 = plot the contrast between conditions (averaged waveform)
S.save_label_waveaverage = 0;
S.label_plot = 'c';
S.color_line = [1 0 0; 0 0 1];

%t-tests
S.t_test_for_permutations = 1;
S.cond_ttests_tobeplotted_topoplot = [2 1]; %this is for both topoplot and t-tests; here [1 2] means cond1 vs cond2

%topoplotting
S.topoplot_label = 0; %1 = plot topoplots
S.fieldtrip_mask = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External';
S.topocontr = 0;
S.topocondsing = [2]; %condition for topoplot
%S.xlim = [1.05 1.15]; %time topolot
S.xlim = [0.913 1.160]; %time topolot (cluster I)
S.zlimmag = [-3.38 3.38]; %magnetometers amplitude topoplot limits
S.zlimgrad = [-2.51 2.51]; %gradiometers amplitude topoplot limits
S.colormap_spec = 0;
% x = []; x.bottom = [0 0 1]; x.botmiddle = [0 0.5 1]; x.middle = [1 1 1]; x.topmiddle = [1 1 0.5]; x.top = [1 0.95 0]; %yellow - blue
x = []; x.bottom = [0 0 0.5]; x.botmiddle = [0 0.5 1]; x.middle = [1 1 1]; x.topmiddle = [1 0 0]; x.top = [0.6 0 0]; %red - blue
S.colormap_spec_x = x;
S.topoplot_save_label = 0;
S.outdir = outdir; %output path

[out] = MEG_sensors_plotting_ttest_LBPD_D2(S); %open this for more information on the function and parameters


%% Loading t-test results and reshaping them

%the analyses have been conducted for planar gradiometers only

%input information
contrn = 0; %1 = load memorized vs novel t-tests; 0 = load novel vs memorized t-tests
block = 2; %1 = atonal; 2 = tonal

grad = 1; %1 = gradiometers; 0 = magnetometers
p_thresh = 0.01; %for binarising p-values matrices

%actual computation
if grad == 1 %gradiometers
    min_time_point = 16; %16 = 0 seconds (first 15 points are pre-stimulus time) %%time-points to be selected (in time samples)
    max_time_point = 399;
else %magnetometers on the basis of gradiometers results
    min_time_point = 16;
    max_time_point = 399;
end
clear DATAP2 TSTAT2
outdir = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/MEG_sensors'; %path where t-test results are stored in
if contrn == 1
    if block == 1 %atonal
        load([outdir '/recogatonal_all_conditions_OUT_Old_Correct_vs_New_Correct.mat']);
    else %tonal
        load([outdir '/recogmajor_all_conditions_OUT_Old_Correct_vs_New_Correct.mat']);
    end
else
    if block == 1 %atonal
        load([outdir '/recogatonal_all_conditions_OUT_New_Correct_vs_Old_Correct.mat']);
    else %tonal
        load([outdir '/recogmajor_all_conditions_OUT_New_Correct_vs_Old_Correct.mat']);
    end
end
chanlabels = OUT.chanlabels;
time_sel = OUT.time_sel;
%here gradiometers and magnetometers are always extracted, but in the following steps only the requested channels (either magnetometers or gradiometers) are used
DATAP2(:,:,1) = OUT.TSTATP_mag(:,min_time_point:max_time_point);
DATAP2(:,:,2) = OUT.TSTATP_grad(:,min_time_point:max_time_point);
TSTAT2(:,:,1) = abs(OUT.TSTAT_mag(:,min_time_point:max_time_point));
TSTAT2(:,:,2) = OUT.TSTAT_grad(:,min_time_point:max_time_point);
chanlab = chanlabels(1:2:204)';
label = zeros(length(chanlab),1); %channels label
for ii = 1:length(chanlab)
    label(ii,1) = str2double(chanlab{ii,1}(4:end));
end

%individuating positive vs negative t-values
P = DATAP2; %trick for looking only into the positive t-values
P(P < p_thresh) = 1; %binarizing p-values according to threshold
P(P < 1) = 0;
%old (pos) > new (neg)
P(TSTAT2 < 0) = 0; %deleting p-values for negative contrasts
TSTAT_mag_pos = P(:,:,1);
TSTAT_grad = P(:,:,2);
%negative vs positive p-values
P = DATAP2;
P(P < p_thresh) = 1; %binarizing p-values according to threshold
P(P < 1) = 0;
P(TSTAT2 > 0) = 0; %deleting p-values for positive contrasts
TSTAT_mag_neg = P(:,:,1);
%load a 2D approximation of the MEG channels location
[~,~,raw_channels] = xlsread('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MatrixMEGChannelLayout_With0_2.xlsx');
% reshaping data for computational purposes
S = [];
S.label = label;
S.TSTAT_mag_pos = TSTAT_mag_pos;
S.TSTAT_mag_neg = TSTAT_mag_neg;
S.TSTAT_grad = TSTAT_grad;
S.raw_channels = raw_channels;

[MAG_data_pos, MAG_data_neg, GRAD_data] = MEG_sensors_MCS_reshapingdata_LBPD_D(S);

%actual Monte Carlo simulations
S = [];
if grad == 1 %actual gradiometers data
    S.data(:,:,:,1) = GRAD_data;
    S.data(:,:,:,2) = zeros(size(MAG_data_pos,1),size(MAG_data_pos,2),size(MAG_data_pos,3));
    S.data(:,:,:,3) = zeros(size(MAG_data_pos,1),size(MAG_data_pos,2),size(MAG_data_pos,3));
else %actual magnetometers data
    S.data(:,:,:,1) = zeros(size(MAG_data_pos,1),size(MAG_data_pos,2),size(MAG_data_pos,3));
    S.data(:,:,:,2) = MAG_data_pos;
    S.data(:,:,:,3) = MAG_data_neg;
end
S.sensortype = [];
S.MEGlayout = cell2mat(raw_channels);
S.permut = 1000;
S.clustmax = 1;
S.permthresh = 0.001;

[MAG_clust_pos, MAG_clust_neg, GRAD_clust] = MEG_sensors_MonteCarlosim_LBPD_D(S);

%% Getting information about the time-points and plotting

%here we create the images for the significant clusters identified in the previous section
%the figures can be found in SUPPLEMENTARY MATERIALS

%OBS! check cluster number, block (tonal or atonal), MEG sensor type, condition (memorized vs. novel or novel vs. memorized), and time interval
clustnum = 4;
block = 1; %1 = atonal; 2 = tonal
typec = 1; %1 = gradiometers; 2 = magnetometers positive; 3 = magnetometers negative
cond = 2; %1 = memorized vs. novel; 2 = novel vs. memorized
outdir = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/MEG_sensors';

if block == 1 %atonal
    sensl = 'recogatonal_all_conditions';
    load([outdir '/' sensl '.mat']);
    if cond == 1
        load([outdir '/Final_sensor_clusters/atonal_oldvsnew_p01_p001_perm1000.mat']);
    else
        load([outdir '/Final_sensor_clusters/atonal_newvsold_p01_p001_perm1000.mat']);
    end
else %tonal
    sensl = 'recogmajor_all_conditions';
    load([outdir '/' sensl '.mat']);
    if cond == 1
        load([outdir '/Final_sensor_clusters/major_oldvsnew_p01_p001_perm1000.mat']);
    else
        load([outdir '/Final_sensor_clusters/major_newvsold_p01_p001_perm1000.mat']);
    end
end
S = [];

%computing data
S.outdir = outdir;
S.data = [];
S.data = data_mat;
S.chanlabels = chanlabels;
S.time_real = time_sel;
S.conditions = {'Old_Correct','New_Correct'};
S.timeextract = []; %time-points to be extracted
S.centerdata0 = 0; %1 = data starting at 0
S.save_data = 0; %only meaningfull if you read data from SPM objects saved on disk
S.save_name_data = sensl;

%individual waveform plotting (settings for individual channel plots are not used in the paper..)
S.waveform_singlechannels_label = 0; %1 = plot single channel waveforms
S.wave_plot_conditions_together = 0; %1 = plot the average of all
S.mag_lab = 2; %1 = magnetometers; 2 = gradiometers
S.x_lim_temp_wave = [0 2.5]; %limits for time (in secs) (e.g. [-0.1 3.4])
S.y_lim_ampl_wave = []; %limit for amplitude (e.g. [0 120] magnetometes, [0 6] gradiometers)

%averaged waveform plotting
S.waveform_average_label = 1;
S.legc = 0; %set 1 for legend
if typec == 1
    clustplot = GRAD_clust{clustnum,3};
else
    clustplot = MAG_clust_pos{clustnum,3};
end
S.left_mag = [];
for ii = 1:size(clustplot,1)
    S.left_mag(ii) = find(cellfun(@isempty,strfind(chanlabels,clustplot{ii,1})) == 0);
end
S.signtp(1) = {[0.653 0.74]}; %in seconds
S.sr = 150; %sampling rate (Hz)
S.avewave_contrast = 1; %1 to plot the contrast between conditions (averaged waveform)
S.save_label_waveaverage = 0;
S.label_plot = 'block_minor';

%t-tests
S.t_test_for_permutations = 0;
S.cond_ttests_tobeplotted_topoplot = [1 2]; %this is for both topoplot and t-tests; here [1 2] means cond1 vs cond2

%topoplotting
S.topoplot_label = 0; %1 = plot topoplot
S.fieldtrip_mask = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External';
S.topocontr = 1;
S.topocondsing = [1];
S.xlim = [0.1333 0.1867]; %time topolot
S.zlimmag = [-5.3 5.3]; %magnetometers amplitude topoplot limits
S.zlimgrad = [-4.6 4.6]; %gradiometers amplitude topoplot limits
S.colormap_spec = 0;
% x = []; x.bottom = [0 0 1]; x.botmiddle = [0 0.5 1]; x.middle = [1 1 1]; x.topmiddle = [1 1 0.5]; x.top = [1 0.95 0]; %yellow - blue
x = []; x.bottom = [0 0 0.5]; x.botmiddle = [0 0.5 1]; x.middle = [1 1 1]; x.topmiddle = [1 0 0]; x.top = [0.6 0 0]; %red - blue
S.colormap_spec_x = x;
S.topoplot_save_label = 0;

[out] = MEG_sensors_plotting_ttest_LBPD_D2(S);



%% Extracting first and last significant time-point for each channel, converting them into seconds (from time-samples) and printing it as xlsx file 

%here you can print the statistics in tables/excel files that can be convenient to read
%you can simply select the significant cluster number that you want by specifying it in the line below
%the tables can be found in SUPPLEMENTARY MATERIALS

clustnum = 1; %number of cluster
block = 4; % 1 = atonal novel vs memorized; 2 = atonal memorized vs novel; 3 = tonal novel vs memorized; 4 = tonal memorized vs novel

list = dir(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/MEG_sensors/Final_sensor_clusters/*mat']);
load([list(block).folder '/' list(block).name])

for clustnum = 1:length(GRAD_clust)
    hh = GRAD_clust{clustnum,3};
    load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/MEG_sensors/time.mat')
    min_time_point = 16; %16 = 0 seconds (first 15 points are pre-stimulus time)
    max_time_point = 399;
    time_sel2 = time_sel(1,min_time_point:max_time_point);
    clear ff2
    ff2(:,1) = hh(:,1); %extracting channel names
    for ii = 1:size(ff2,1)
        ff2(ii,2) = {time_sel2(hh{ii,2}(1))}; %first significant time-point
        ff2(ii,3) = {time_sel2(hh{ii,2}(end))}; %last significant time-point
    end
    PDn = cell2table(ff2); %converting cell to table
    writetable(PDn,['Block_' num2str(block) '.xlsx'],'Sheet',clustnum) %saving xlsx file
end

%% COREGISTRATION

%(1) creating MRI folders
%(2) moving the MRI nifti files
%(2) actual rhino coregistration
%(3) checking rhino coregistration

%% Creating MRI subjects' folders

%here we create one folder for each subject to save the nifti files

v_ID = 1:71; %vector with subjects' ID

for ii = 1:length(v_ID)
    mkdir(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/MRI_nifti_images_check/00' num2str(v_ID(ii))])
end

%after this, I manually added one '0' to the first nine folders (in order to have four digits on all the folders' names)

%% Moving the MRI nifti files

%we copy the MRI nifti files from the original subjects' folders to the newly created folders

MRI_path = '/projects/MINDLAB2017_MEG-LearningBach/scratch/Leonardo/LearningBach/after_maxfilter_mc/NIFTI_INV2'; %path to subjects' MRI folders
original_folders = dir([MRI_path '/0*']); %subjects' MRI folders
output_path = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/MRI_nifti_images'; %path to subjects' new MRI folders
new_folders = dir([output_path '/0*']); %new subjects' MRI folders

%OBS! make sure that MRI_path and output_path have the same number of folders!

%OBS! if the number of folders is not the same: delete the extra folders from output_path (for subjects that use a template) and manually add them after running the loop

for ii = 1:length(original_folders) %over folders
    name_original = [original_folders(ii).folder '/' original_folders(ii).name];
    file_to_copy = dir([name_original '/*.nii']); %list nifti files
    for jj = 1:length(file_to_copy)
        if ~strcmp(file_to_copy(jj).name(1), 'y') %not the nifti files that start with 'y'
            nifti_file = [file_to_copy(jj).folder '/' file_to_copy(jj).name]; %nifti file we need
        end
    end
    name_new = [new_folders(ii).folder '/' new_folders(ii).name];
    if name_original(94:95) == name_new(79:80) %compare the name of the folders
        copyfile(nifti_file,name_new) %copy the nifti files into the new folders
    end
    disp(ii)
end

%after this, we manually created new folders for subjects 0007, 0033, 0036 and 0041 and copied the MRI template into them

%% Actual RHINO coregistration

addpath('/projects/MINDLAB2020_MEG-AuditoryPatternRecognition/scripts/gemma')

%lists of new epoched files in which the coregistration needs to be pasted in:
list_atonal = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/e80*recogatonal_tsssdsm.mat'); %atonal files
list_major = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/e80*recogmajor_tsssdsm.mat'); %tonal files

%running rhino
%OBS! check that all MEG data are in the same order and number as the MRI nifti files!
a = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/MRI_nifti_images/00*'); %set path to subjects' MRI folders
for ii = 1:length(list_major) %OBS! set this to atonal or tonal
    S = [];
    S.ii = ii;
    %S.D = [list_atonal(ii).folder '/' list_atonal(ii).name]; %path to atonal files
    S.D = [list_major(ii).folder '/' list_major(ii).name]; %path to tonal files
    dummymri = dir([a(ii).folder '/' a(ii).name '/*.nii']); %path to nifti files
    S.mri = [dummymri(1).folder '/' dummymri(1).name];
    %standard parameters
    S.useheadshape = 1;
    S.use_rhino = 1; %set 1 for rhino, 0 for no rhino
    S.forward_meg = 'Single Shell';
    S.fid.label.nasion = 'Nasion';
    S.fid.label.lpa = 'LPA';
    S.fid.label.rpa = 'RPA';
    jobid = job2cluster(@coregfunc_g,S); %running with parallel computing
end

%% Checking RHINO coregistration

%check whether the brain is centered
%green dots represent the MEG sensors
%pink rhombi represent the fiducials (nasion and ears)
%colored dots represent the pen digitization

%check this for all subjects and both atonal and major

%OBS! consider removing subjects 0053, 0054 and 0057 in the future

for ii = 1:length(list_major)
    D = spm_eeg_load([list_major(ii).folder '/' list_major(ii).name]);
    rhino_display(D)
end

%% SOURCE RECONSTRUCTION

%(1) beamforming
%(2) moving files into a common folder
%(3) first level, subject level, group level, and creating nifti images

%% Overview

%actual source reconstruction
%independent sections for different blocks (atonal and tonal) and frequency bands (0.1-1Hz and 2-8Hz)

%four steps:
%1) beamforming
%2) first level (independent analyses for each subject)
%3) subject level (not relevant, but necessary for OSL to properly function)
%4) group level (all subjects or groups of subjects together)

%settings for cluster (parallel computing)
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/scripts_osl_learningbach/Cluster')
clusterconfig('scheduler', 'cluster'); %set automatically the long run queue
clusterconfig('long_running', 1); %set automatically the long run queue
clusterconfig('slot', 2); %set manually the job cluster slots


%% Beamforming

block = 2; %1 = atonal; 2 = tonal
freqq = 2; %1 = 0.1-1Hz; 2 = 2-8Hz; 3 = broadband (0.1 - 40 Hz)

addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing') %adds path to function that submits jobs to the server
if block == 1 %atonal
    %list of the epoched files:
    list_block = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/e80*recogatonal_tsssdsm.mat');
    s = 1:71; %subjects in list_block
else %tonal
    list_block = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/e80*recogmajor_tsssdsm.mat');
    s = 1:71;
end

for ii = 1:length(s)
    processed_file = [list_block(s(ii)).folder '/' list_block(s(ii)).name]; %epoched file
    D_epoched = spm_eeg_load(processed_file); %load epoched file
    D_epoched = D_epoched.montage('switch',1); %montage 1 (AFRICA denoised data)
    D_epoched.save();
    %Beamform
    oat = [];
    oat.source_recon.D_epoched(str2double(list_block(s(ii)).name(18:21))) = {processed_file};
    pca_dim_1 = 50; %number of principal components
    oat.source_recon.pca_dim = pca_dim_1(ones(length(oat.source_recon.D_epoched),1)); %this seems to be necessary %it is for having the same number of values than the number of input D objects
    oat.source_recon.modalities = {'MEGMAG'; 'MEGPLANAR'};
    oat.source_recon.conditions = {'Old_Correct','New_Correct'};
    oat.source_recon.gridstep = 8; %in mm
    oat.source_recon.time_range = [-0.1 3.4]; %time range in secs
    if freqq == 1
        oat.source_recon.freq_range = [0.1 1]; %frequency range in Hz
    elseif freqq == 2
        oat.source_recon.freq_range = [2 8]; %frequency range in Hz
    else
        oat.source_recon.freq_range = [0.1 40]; %frequency range in Hz
    end
    %S.source_recon.pca_order         = 250;
    oat.source_recon.type = 'Scalar';
    oat.source_recon.method = 'beamform';
    oat.source_recon.normalise_method = 'mean_eig';
    oat.source_recon.forward_meg = 'MEG Local Spheres';
    %S.source_recon.prefix            = '';
    oat.source_recon.report.do_source_variance_maps = 1;
    oat.source_recon.sessions_to_do = [];
    oat.source_recon.sessions_to_do = str2double(list_block(s(ii)).name(18:21)); %sessions to do among the file_list
    oat.source_recon.dirname = [list_block(1).folder '/source/Block_' num2str(block) '_freq_' num2str(oat.source_recon.freq_range(1)) '_' num2str(oat.source_recon.freq_range(2)) '_SUBJ_' list_block(s(ii)).name(18:21)];
    
    jobid = job2cluster(@cluster_beamforming,oat); %running with parallel computing
end

%% Moving files into a common folder

%the files that were computed independently for each subject in Beamforming are moved into a common new folder
%this way, we go back to the "usual" way of running oat analysis pipeline in OSL

%OBS! this is run only the first time
%each block is repeated three times: one for 0.1-1Hz, one for 2-8Hz, and one for 0.1-40Hz
%block 1 = atonal; block 2 = tonal

%Block 1 - freq 0.1 - 2Hz
list_block = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/Block_1_freq_0.1_1*');
mkdir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source', 'FirLev_Block_1_freq_0.1_1.oat');
bs = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/FirLev_Block_1_freq_0.1_1.oat';
for ii = 1:length(list_block)
    %concat file .dat
    as = [list_block(ii).folder '/' list_block(ii).name '/concatMfsession' num2str(str2double(list_block(ii).name(25:28))) '_spm_meeg.dat'];
    status = movefile(as,bs)
    %concat file .mat
    as = [list_block(ii).folder '/' list_block(ii).name '/concatMfsession' num2str(str2double(list_block(ii).name(25:28))) '_spm_meeg.mat'];
    status = movefile(as,bs)
    %session file
    as = [list_block(ii).folder '/' list_block(ii).name '/session' num2str(str2double(list_block(ii).name(25:28))) '_recon.mat'];
    status = movefile(as,bs)
end

%Block 2 - freq 0.1 - 1Hz
list_block = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/Block_2_freq_0.1_1*');
mkdir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source', 'FirLev_Block_2_freq_0.1_1.oat');
bs = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/FirLev_Block_2_freq_0.1_1.oat';
for ii = 1:length(list_block)
    %concat file .dat
    as = [list_block(ii).folder '/' list_block(ii).name '/concatMfsession' num2str(str2double(list_block(ii).name(25:28))) '_spm_meeg.dat'];
    status = movefile(as,bs)
    %concat file .mat
    as = [list_block(ii).folder '/' list_block(ii).name '/concatMfsession' num2str(str2double(list_block(ii).name(25:28))) '_spm_meeg.mat'];
    status = movefile(as,bs)
    %session file
    as = [list_block(ii).folder '/' list_block(ii).name '/session' num2str(str2double(list_block(ii).name(25:28))) '_recon.mat'];
    status = movefile(as,bs)
end

%Block 1 - freq 2 - 8Hz
list_block = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/Block_1_freq_2*');
mkdir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source','FirLev_Block_1_freq_2_8.oat');
bs = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/FirLev_Block_1_freq_2_8.oat';
for ii = 1:length(list_block)
    %concat file .dat
    as = [list_block(ii).folder '/' list_block(ii).name '/concatMfsession' num2str(str2double(list_block(ii).name(23:26))) '_spm_meeg.dat'];
    status = movefile(as,bs)
    %concat file .mat
    as = [list_block(ii).folder '/' list_block(ii).name '/concatMfsession' num2str(str2double(list_block(ii).name(23:26))) '_spm_meeg.mat'];
    status = movefile(as,bs)
    %session file
    as = [list_block(ii).folder '/' list_block(ii).name '/session' num2str(str2double(list_block(ii).name(23:26))) '_recon.mat'];
    status = movefile(as,bs)
end

%Block 2 - freq 2 - 8Hz
list_block = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/Block_2_freq_2*');
mkdir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source','FirLev_Block_2_freq_2_8.oat');
bs = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/FirLev_Block_2_freq_2_8.oat';
for ii = 1:length(list_block)
    %concat file .dat
    as = [list_block(ii).folder '/' list_block(ii).name '/concatMfsession' num2str(str2double(list_block(ii).name(23:26))) '_spm_meeg.dat'];
    status = movefile(as,bs)
    %concat file .mat
    as = [list_block(ii).folder '/' list_block(ii).name '/concatMfsession' num2str(str2double(list_block(ii).name(23:26))) '_spm_meeg.mat'];
    status = movefile(as,bs)
    %session file
    as = [list_block(ii).folder '/' list_block(ii).name '/session' num2str(str2double(list_block(ii).name(23:26))) '_recon.mat'];
    status = movefile(as,bs)
end

%Block 1 - freq 0.1 - 40Hz
list_block = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/Block_1_freq_0.1_40*');
mkdir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source', 'FirLev_Block_1_freq_0.1_40.oat');
bs = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/FirLev_Block_1_freq_0.1_40.oat';
for ii = 1:length(list_block)
    %concat file .dat
    as = [list_block(ii).folder '/' list_block(ii).name '/concatMfsession' num2str(str2double(list_block(ii).name(26:29))) '_spm_meeg.dat'];
    status = movefile(as,bs)
    %concat file .mat
    as = [list_block(ii).folder '/' list_block(ii).name '/concatMfsession' num2str(str2double(list_block(ii).name(26:29))) '_spm_meeg.mat'];
    status = movefile(as,bs)
    %session file
    as = [list_block(ii).folder '/' list_block(ii).name '/session' num2str(str2double(list_block(ii).name(26:29))) '_recon.mat'];
    status = movefile(as,bs)
end

%Block 2 - freq 0.1 - 40Hz
list_block = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/Block_2_freq_0.1_40*');
mkdir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source', 'FirLev_Block_2_freq_0.1_40.oat');
bs = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/FirLev_Block_2_freq_0.1_40.oat';
for ii = 1:length(list_block)
    %concat file .dat
    as = [list_block(ii).folder '/' list_block(ii).name '/concatMfsession' num2str(str2double(list_block(ii).name(26:29))) '_spm_meeg.dat'];
    status = movefile(as,bs)
    %concat file .mat
    as = [list_block(ii).folder '/' list_block(ii).name '/concatMfsession' num2str(str2double(list_block(ii).name(26:29))) '_spm_meeg.mat'];
    status = movefile(as,bs)
    %session file
    as = [list_block(ii).folder '/' list_block(ii).name '/session' num2str(str2double(list_block(ii).name(26:29))) '_recon.mat'];
    status = movefile(as,bs)
end

%checking whether all the files were moved
bs_1_1 = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/FirLev_Block_1_freq_0.1_1.oat'); %block 1 freq 1
bs_2_1 = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/FirLev_Block_2_freq_0.1_1.oat'); %block 2 freq 1
bs_1_2 = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/FirLev_Block_1_freq_2_8.oat'); %block 1 freq 2
bs_2_2 = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/FirLev_Block_2_freq_2_8.oat'); %block 2 freq 2
bs_1_3 = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/FirLev_Block_1_freq_0.1_40.oat'); %block 1 freq 3
bs_2_3 = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/FirLev_Block_2_freq_0.1_40.oat'); %block 2 freq 3

%% First level, subject level, group level, and creating nifti images

block = 2; %1 = atonal; 2 = tonal
freqq = 3; %1 = 0.1-1Hz; 2 = 2-8Hz; 3 = 0.1-40Hz
first_lev = 0; %1 to run first-level analysis
N100_coape = 1; %1 to run first-level analysis in 2-8Hz with coape
%absolute value after the reconstruction and before contrasts); set to 1
%also if you want to run the subsequent steps on this particular data
subj_lev = 0; %1 to run subject-level analysis
group_lev = 1; %1 to run group-level analysis
subj_nifti = 0; %1 to run subject-level nifti images creation
group_nifti = 0; %1 to run group-level nifti images creation

%OBS! do not run the next two steps together: (1) cluster stats, (2) fsl_maths
cluster_based_fsl_statistics = 0; % 1 to run cluster-based statistics by using randomise (FSL) through OSL
fsl_maths_cluster = 0; % 1 to run fsl_maths
%fasl_maths uses the significant clusters as masks for the voxel statistics (always run locally)
%this is crucial to have proper nifti images to be used in Workbench

%codes to get a "normal" oat for source reconstruction (all subjects in the same folder)
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing') %path to the function that submits the jobs to the server
if block == 1
    if freqq == 1
        outdir = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/FirLev_Block_1_freq_0.1_1.oat';
    elseif freqq == 2
        outdir = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/FirLev_Block_1_freq_2_8.oat';
    else
        outdir = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/FirLev_Block_1_freq_0.1_40.oat';
    end
    list_block = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/e80*recogatonal_tsssdsm.mat');
    s = 1:71;
elseif block == 2
    if freqq == 1
        outdir = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/FirLev_Block_2_freq_0.1_1.oat';
    elseif freqq == 2
        outdir = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/FirLev_Block_2_freq_2_8.oat';
    else
        outdir = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/FirLev_Block_2_freq_0.1_40.oat';
    end
    list_block = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/e80*recogmajor_tsssdsm.mat');
    s = 1:71;
else
    warning('block number is wrong.. check it better..')
end

%actual creation of oat structure
oat = [];
%cnt = 0;
for ii = 1:length(s)
    oat.source_recon.D_epoched(ii) = {[list_block(s(ii)).folder '/' list_block(s(ii)).name]};
    %oat.source_recon.sessions_to_do(ii) = {ii}; 
    oat.source_recon.sessions_to_do(ii) = {[str2double(list_block(s(ii)).name(18:21))]};
    %oat.source_recon.results_fnames(ii) = {['session' num2str(str2double(list_cont(s(ii)).name(12:15))) '_recon']};
end
for ii = 1:(str2double(list_block(end).name(18:21)))
    oat.source_recon.results_fnames(ii) = {['session' num2str(ii) '_recon']};
end
%Beamform
pca_dim_1 = 50;
oat.source_recon.pca_dim = pca_dim_1(ones(length(oat.source_recon.D_epoched),1)); %this seems to be necessary. It is for having the same number of values than the number of input D objects
oat.source_recon.modalities = {'MEGMAG'; 'MEGPLANAR'};
oat.source_recon.conditions = {'Old_Correct','New_Correct'};
oat.source_recon.gridstep = 8; %in mm
oat.source_recon.time_range = [-0.1 3.4]; %time range in secs
if freqq == 1
    oat.source_recon.freq_range = [0.1 1]; %frequency range in Hz
elseif freqq == 2
    oat.source_recon.freq_range = [2 8]; %frequency range in Hz
else
    oat.source_recon.freq_range = [0.1 40]; %frequency range in Hz
end
%S.source_recon.pca_order = 250;
oat.source_recon.type = 'Scalar';
oat.source_recon.method = 'beamform';
oat.source_recon.normalise_method = 'mean_eig';
oat.source_recon.forward_meg = 'MEG Local Spheres';
%S.source_recon.prefix            = '';
oat.source_recon.report.do_source_variance_maps = 1;
oat.source_recon.dirname = [];
oat.source_recon.dirname = outdir(1:end-4);

%first level
%each experimental block for each subject, independently
%same for both blocks
design_matrix_summary = {};
design_matrix_summary{1} = [1 0];design_matrix_summary{2} = [0 1];
%contrasts to be calculated:
oat.first_level.contrast = {};
%contrast design matrix
oat.first_level.contrast{1} = [1 0]'; %old
oat.first_level.contrast{2} = [0 1]'; %new
oat.first_level.contrast{3} = [1 -1]'; %old - new
%contrast names
oat.first_level.contrast_name = {};
oat.first_level.contrast_name{1} = 'old';
oat.first_level.contrast_name{2} = 'new';
oat.first_level.contrast_name{3} = 'old-new';
%baseline correction
oat.first_level.bc = [1 1 0];
%design matrix summary
oat.first_level.design_matrix_summary = design_matrix_summary;
%original contrasts
oat.first_level.report.first_level_cons_to_do = [1:length(oat.first_level.contrast)]; %[1 2 3]; %better to do 3 2 1 in order to get the information for the peak value of the contrast old-new
oat.first_level.time_range = [-0.1 3.4]; %time range in secs
oat.first_level.post_tf_downsample_factor = 1;
if freqq == 1
    %slow negativity
    oat.first_level.cope_type = 'coape'; %this should calculate absolute values of glm copes before comparing different conditions (otherwise, by default the absolute value is done on the copes after comparing different conditions); the absolute value here seems to be necessary because of the sign ambiguity issue introduced by source reconstruction
    oat.first_level.name = ['wholebrain_first_level_BC_coape']; %OBS! remember to check this name
else
    %N100
    if N100_coape == 1
        oat.first_level.cope_type = 'coape'; %this should calculate absolute values of glm copes before comparing different conditions (otherwise, by default the absolute value is done on the copes after comparing different conditions); the absolute value here seems to be necessary because of the sign ambiguity issue introduced by source reconstruction
        oat.first_level.name = ['wholebrain_first_level_BC_N100_coape_real']; %OBS! remember to check this name
    else
        oat.first_level.cope_type = 'none';
        oat.first_level.name = ['wholebrain_first_level_BC_coape'];
    end
end

if first_lev == 1
    %running first level on parallel computing
    for ii = 1:length(s)
        oat.first_level.sessions_to_do = [];
        oat.first_level.sessions_to_do = [str2double(list_block(s(ii)).name(18:21))]; %here it seems that the session indexes oat.source_recon.results_fnames{ii} is directly related to the sequential
        %oat.first_level.sessions_to_do = [str2double(list(ii).name(1:4))];
        jobid = job2cluster(@cluster_beamfirstlevel,oat);
    end
end

%subject level
%this does not do anything since we have only one experimental block for each participant (however, for computational reasons, we have to run it)
%this is needed to read the proper subjects
for ii = 1:(str2double(list_block(end).name(18:21))) %everybody
    oat.first_level.results_fnames(ii) = {['session' num2str(ii) '_' oat.first_level.name '.mat']};
end
oat.subject_level.name = 'subj';
oat.subject_level.session_index_list = cell(1,(str2double(list_block(end).name(18:21))));
for ii = 1:(str2double(list_block(end).name(18:21)))
    oat.subject_level.session_index_list{ii} = ii;
    oat.subject_level.results_fnames(ii) = {['subject' num2str(ii) '_' oat.first_level.name '_' oat.subject_level.name '.mat']};
end
if subj_lev == 1
    for ii = 1:length(s)
        oat.subject_level.subjects_to_do = [];
        oat.subject_level.subjects_to_do = str2double(list_block(s(ii)).name(18:21));
        jobid = job2cluster(@cluster_beamsubjlevel,oat);
    end
end

%group level
for ii = 1:(str2double(list_block(end).name(18:21))) %everybody
    oat.subject_level.subjects_to_do(ii) = ii;
end
oat.group_level = [];
oat.group_level.name = 'group_level_everybody';
oat.group_level.subjects_to_do = []; %original IDs
for ii = 1:length(s)
    oat.group_level.subjects_to_do(ii) = str2double(list_block(s(ii)).name(18:21));
end
%results name
oat.group_level.results_fnames = [oat.first_level.name '_' oat.subject_level.name '_' oat.group_level.name '.mat'];
if ~exist([oat.source_recon.dirname '/' oat.group_level.results_fnames],'var')
    copyfile(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter/source/mask.nii.gz'],[oat.source_recon.dirname '.oat/' oat.first_level.name '_' oat.subject_level.name '_' oat.group_level.name '_mask.nii.gz'])
end
%Spatial and temporal averaging options
oat.group_level.time_range = [-0.1 3.4]; %time range in secs
oat.group_level.space_average = 0;
oat.group_level.time_average = 0;
oat.group_level.time_smooth_std = 0; %secs
oat.group_level.use_tstat = 0;
%path to AAL template
%oat.group_level.mask_fname = '/projects/MINDLAB2017_MEG-LearningBach/scripts/osl/AAL_Nifti/aal_8mm_try5.nii.gz';
%Spatial and temporal smoothing options
oat.group_level.spatial_smooth_fwhm = 0; %mm
oat.group_level.group_varcope_time_smooth_std = 100;
oat.group_level.group_varcope_spatial_smooth_fwhm = 100; %smooths the variance of the group copes. It is recommended to do this
%store copes (useful for doing the permutation test later, otherwise it needs to compute again the group level analysis)
oat.group_level.store_lower_level_copes = 1;
%Set up design matrix and contrasts
%this is if you have only the general mean across the all participants
oat.group_level.group_design_matrix = ones(1,length(oat.group_level.subjects_to_do)); %if you want all of the participants
oat.group_level.group_contrast = [];
oat.group_level.group_contrast{1} = [1];
oat.group_level.group_contrast_name = {};
oat.group_level.group_contrast_name{1} = 'mean';
%oat.group_level.glm_method='ols'; %ols or fixed-effects
%Define which contrasts to perform for the report
oat.group_level.first_level_contrasts_to_do = [1:length(oat.first_level.contrast)]; %list of first level contrasts to run the group analysis on
oat.group_level.report.first_level_cons_to_do = [1]; %purely used to determine which contrasts are shown in the report, the first one in list determines the contrast used to find the vox, time, freq with the max stat
oat.group_level.report.group_level_cons_to_do = [1]; %purely used to determine which contrasts are shown in the report, the first one in list determines the contrast used to find the vox, time, freq with the max stat
oat.group_level.report.show_lower_level_copes = 0;
oat.group_level.report.show_lower_level_cope_maps = 0;
if group_lev == 1
    jobid = job2cluster(@cluster_beamgrouplevel,oat);
end

%subject-level - nifti
if subj_nifti == 1
    oat.source_recon.dirname = [oat.source_recon.dirname '.oat']; %inserting ".oat" to allow the function to read the file
    if ~exist([oat.source_recon.dirname '/' oat.group_level.results_fnames],'var')
        copyfile(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter/source/mask.nii.gz'],[oat.source_recon.dirname '/' oat.first_level.name '_mask.nii.gz'])
    end
    for ii = 1:length(oat.group_level.subjects_to_do)
        addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/scripts_osl_learningbach/Cluster') %path to function that submits the jobs to the server
        S2 = [];
        S2.oat = oat;
        S2.stats_fname = oat.subject_level.results_fnames{oat.group_level.subjects_to_do(ii)};
        %S2.stats_fname = oat.subject_level.results_fnames{find(oat.group_level.subjects_to_do==7)};
        S2.first_level_contrasts = 1:length(oat.first_level.contrast); %remember that in this way you define the order of the output (this numbers refers to the order (numbers) defined in the contrasts; for example here tstat3 refers to the contrast memorized-novel)
        S2.resamp_gridstep = oat.source_recon.gridstep;
        jobid = job2cluster(@cluster_oat_save_nii_stats,S2);
    end
end

%group-level - nifti
if group_nifti == 1
    S2 = [];
    oat.source_recon.dirname = [oat.source_recon.dirname '.oat']; %inserting ".oat" to allow the function to read the file
    S2.oat = oat;
    S2.stats_fname = oat.group_level.results_fnames;
    S2.first_level_contrasts = 1:length(oat.first_level.contrast); %remember that in this way you define the order of the output (this numbers refers to the order (numbers) defined in the contrasts; for example here tstat3 refers to the contrast memorized-novel)
    S2.group_level_contrasts = [1];
    S2.resamp_gridstep = oat.source_recon.gridstep;
    jobid = job2cluster(@cluster_oat_save_nii_stats,S2);
end

%clusters (OSL - FSL randomise)
if cluster_based_fsl_statistics == 1
    t_val = 2; %t-value to binarize statistics
    oat.source_recon.dirname = [oat.source_recon.dirname '.oat']; %inserting ".oat" to allow the function to read the file
    %same for both blocks
    %positive contrasts
    time_tones_356 = [0 0.250 0.500 0.750 1.000]; %onset of tones (blocks 1 2)
    for ii = 1:5 %over musical tones
        S = [];
        oat.clustname = ['tone_' num2str(ii)];
        S.oat = oat; %oat computed in this section
        S.cluster_stats_thresh = t_val; %t-value to binarize statistics
        S.cluster_stats_nperms = 5000; %permutations
        S.first_level_copes_to_do = oat.group_level.first_level_contrasts_to_do; %all the contrasts
        S.group_level_copes_to_do = 1;
        S.group_varcope_spatial_smooth_fwhm = oat.group_level.group_varcope_spatial_smooth_fwhm;
        S.write_cluster_script = 0;
        S.time_range = [time_tones_356(ii) (time_tones_356(ii) + 0.250)];
        S.time_average=1;
        %Run the permutations (on cluster.. parallel computing)
        jobid = job2cluster(@clusterbasedpermutation_osl,S);
    end
    %doing clusters for reversed contrasts (e.g. cond2>cond1)
    %reversing statistics by multiplying it by -1
    if N100_coape ~= 1
        warning('you are computing cluster statistics on contrasted data with polarity.. thus the interpretation of the results/figures that I am outputting now may be highly controversial..')
    end
    load([oat.source_recon.dirname '/' oat.group_level.results_fnames]); %loading data
    oat.group_level.results_fnames = [oat.group_level.results_fnames(1:end-4) '_reversed.mat']; %name for the new data that is being reversed
    oat_stage_results.fname = oat.group_level.results_fnames(1:end-4); %this is important for later operations in the cluster function
    for cc = 1:3 %all contrasts
        ghi = oat_stage_results.lower_level_copes{cc}; %extracting data of interest
        ghi2 = ghi * (-1); %reversing the sign of contrast
        oat_stage_results.lower_level_copes{cc} = ghi2; %storing data
    end
    save([oat.source_recon.dirname '/' oat.group_level.results_fnames],'oat_stage_results','-v7.3'); %saving new reversed data
    %computing clusters as done before; only for contrasts between experimental conditions
    for ii = 1:5 %over musical tones
        S = [];
        oat.clustname = ['tone_' num2str(ii)];
        S.oat = oat; %oat computed in this section
        S.cluster_stats_thresh = t_val; %t-value to binarize statistics
        S.cluster_stats_nperms = 5000; %permutations
        S.first_level_copes_to_do = oat.group_level.first_level_contrasts_to_do;
        S.group_level_copes_to_do = 1;
        S.group_varcope_spatial_smooth_fwhm = oat.group_level.group_varcope_spatial_smooth_fwhm;
        S.write_cluster_script = 0;
        S.time_range = [time_tones_356(ii) (time_tones_356(ii) + 0.250)];
        S.time_average=1;
        jobid = job2cluster(@clusterbasedpermutation_osl,S);
    end
end

%using significant clusters as mask for the voxel statistics
if fsl_maths_cluster == 1
    if cluster_based_fsl_statistics == 1
        hj = oat.source_recon.dirname;
        hjname = oat.group_level.results_fnames(1:end-13);
    else
        hj = [oat.source_recon.dirname '.oat'];
        hjname = oat.group_level.results_fnames(1:end-4);
    end
    listmath = dir([hj '/' hjname '_randomise*']);
    blk = 1:3; %contrasts 1:3
    %adding positive and negative contrasts masks
    for new_folders = blk %over relevant contrasts (two conditions against each other)
        for pp = 1:5 %over musical tones
            %getting filename for reversed
            corrp_rev = dir([hj '/' hjname '_reversed_randomise_c' num2str(new_folders) '_dir_tone_' num2str(pp) '/clustere_corrp*8mm.nii.gz']);
            %getting filename for original
            corrp_ori = dir([hj '/' hjname '_randomise_c' num2str(new_folders) '_dir_tone_' num2str(pp) '/clustere_corrp*8mm.nii.gz']);
            %putting together the two masks (significant clusters for both original and reversed contrasts)
            cmd = ['fslmaths ' corrp_ori(1).folder '/' corrp_ori(1).name ' -add ' corrp_rev(1).folder '/' corrp_rev(1).name ' ' corrp_ori(1).folder '/' corrp_ori(1).name];
            system(cmd)
        end
    end
    %doing the mask operation (fslmath)
    for ii = 1:length(listmath) %over contrasts and tones 
        if strcmp(listmath(ii).name(end-5:end-2),'tone') %only if the contrast if for one of the musical tones
            path_tstat = dir([listmath(ii).folder '/' listmath(ii).name '/tst*8mm.nii.gz']); %looking for proper t-stats
            path_corrp = dir([listmath(ii).folder '/' listmath(ii).name '/clustere_corrp*8mm.nii.gz']); %looking for proper significant clusters (corrected p-vals)
            if isnan(str2double(listmath(ii).name(end-12:end-11))) %this a trick to get the contrast label "c" also when the contrast is 10 (so two characters); the idea is that if we get also "c" in the conversion to number that produces a NaN..
                output_name = [hj '/' listmath(ii).name(end-5:end) '_' listmath(ii).name(end-12:end-11)];
            else
                output_name = [hj '/' listmath(ii).name(end-5:end) '_' listmath(ii).name(end-13:end-11)];
            end
            cmd = ['fslmaths ' path_tstat(1).folder '/' path_tstat(1).name ' -mas ' path_corrp(1).folder '/' path_corrp(1).name ' ' output_name]; %apply fslmaths to create a mask from the corrected b0
            system(cmd)
        end
        disp([num2str(ii) ' / ' num2str(length(listmath))])
    end
end
%after that, we use Workbench for producing the images


%% CLUSTER-BASED MONTE-CARLO SIMULATIONS WITH IN-HOUSE-BUILT FUNCTIONS

%here we follow these steps:
%(1) defining clusters on 3D brain voxel statistics
%(2) combining images with significant clusters

%% Defining clusters on 3D brain voxel statistics

%loading p-values and t-values
clear DATA
freq = 1; %1 = 0.1-1Hz; 2 = 2-8Hz
block = 2; %1 = atonal; 2 = tonal
tone = 5; %set to 1 2 3 4 or 5
tvalbin = 2.3; %value to binarize the data with
contrast = 3; %set to 1 2 or 3
if freq == 1
    T = load_nii(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/FirLev_Block_' num2str(block) '_freq_0.1_1.oat/wholebrain_first_level_BC_coape_subj_group_level_everybody_randomise_c' num2str(contrast) '_dir_tone_' num2str(tone) '/tstat3_gc1_8mm.nii.gz']);
else
    T = load_nii(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/FirLev_Block_' num2str(block) '_freq_2_8.oat/wholebrain_first_level_BC_N100_coape_real_subj_group_level_everybody_randomise_c' num2str(contrast) '_dir_tone_' num2str(tone) '/tstat3_gc1_8mm.nii.gz']);
end

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

%novel vs memorized (removing values below t-value = 2
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
    if freq == 1
        S.T = ['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/FirLev_Block_' num2str(block) '_freq_0.1_1.oat/wholebrain_first_level_BC_coape_subj_group_level_everybody_randomise_c' num2str(contrast) '_dir_tone_' num2str(tone) '/tstat3_gc1_8mm.nii.gz'];
    else
        S.T = ['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/FirLev_Block_' num2str(block) '_freq_2_8.oat/wholebrain_first_level_BC_N100_coape_real_subj_group_level_everybody_randomise_c' num2str(contrast) '_dir_tone_' num2str(tone) '/tstat3_gc1_8mm.nii.gz'];
    end
    S.outdir = ['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/OldvsNew_MCS/' FREQ{freq} '/Block_' num2str(block)]; %output path
    S.parcelfile = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/aal_8mm_try5.nii.gz';
    S.labels = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/AAL_labels.mat';
    S.MNIcoords = mni_coords; %MNI coordinates of 8mm MNI152T1 brain
    S.data = DATA{ii};
    S.mask = mask; %mask of the brain layout with the results
    S.permut = 1000; %number of permutations for Monte Carlo simulation (MCS)
    S.clustmax = 1; %set 1 for only max cluster size of each permutation MCS (more strict); set 0 for every size of each cluster detected for each permutation MCS (less strict)
    S.permthresh = 0.001; %threhsold for MCS
    
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

freq = 1; %1 = 0.1-1Hz; 2 = 2-8Hz
block = 2; %1 = atonal; 2 = tonal
tone = 3; %set to 1 2 3 4 or 5

if freq == 1
    if block == 1
        path = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/OldvsNew_MCS/0.1_1/Block_1';
    else
        path = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/OldvsNew_MCS/0.1_1/Block_2';
    end
else
    if block == 1
        path = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/OldvsNew_MCS/2_8/Block_1';
    else
        path = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/OldvsNew_MCS/2_8/Block_2';
    end
end

%tone x cluster 1
name1 = ['Tval_tone_' num2str(tone) '_Cond1vsCond2_SignClust_1_Tvals.nii.gz'];

%tone x cluster 2
name2 = ['Tval_tone_' num2str(tone) '_Cond1vsCond2_SignClust_2_Tvals.nii.gz'];

%tone x cluster 3
name3 = ['Tval_tone_' num2str(tone) '_Cond1vsCond2_SignClust_3_Tvals.nii.gz'];

output = ['Tval_tone_' num2str(tone) '_All_Clusters.nii.gz'];

% cmd = ['fslmaths ' path '/' name1 ' -add ' path '/' name2 ' ' path '/' output]; %only two clusters
cmd = ['fslmaths ' path '/' name1 ' -add ' path '/' name2 ' -add ' path '/' name3 ' ' path  '/' output]; %OBS! be careful with the spacing
system(cmd)

%after this, we created the final images using Workbench


%% Extracting information about the clusters at source level and reporting it in xlsx files

%here we obtain information about the brain regions forming the clusters
%the tables can be found in SUPPLEMENTARY MATERIALS

freq = 2; %1 = 0.1-1Hz; 2 = 2-8Hz
block = 3; %1 = atonal; 2 = tonal; 3 = tonal vs atonal
tonebytone = 1; %set to 1 2 3 4 or 5

if freq == 1
    if block == 1
        path = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/OldvsNew_MCS/0.1_1/Block_1/FinalImages';
    elseif block == 2
        path = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/OldvsNew_MCS/0.1_1/Block_2/FinalImages';
    else
        path = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/ContrastMajorAtonal/Freq_01_1/Contrast1/Cluster_Permutation_test/FinalImages';
    end
else
    if block == 1
        path = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/OldvsNew_MCS/2_8/Block_1/FinalImages';
    elseif block == 2
        path = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/OldvsNew_MCS/2_8/Block_2/FinalImages';
    else
        path = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/ContrastMajorAtonal/Freq_2_8/Contrast1/Cluster_Permutation_test/FinalImages';
    end
end


for tonebytone = 1:5
    %tone x cluster 1
    if block < 3
        fname = [path '/Tval_tone_' num2str(tonebytone) '_All_Clusters.nii.gz'];
    else
        fname = [path '/Tval_MajvsAto_tone_' num2str(tonebytone) '_All_Clusters.nii.gz'];
    end
    %getting MNI coordinates of significant voxels within the provided image
    [ mni_coords, xform ] = osl_mnimask2mnicoords(fname);
    %loading the image
    V = nii.load(fname);
    %extracting statistics
    VV = V(V~=0);
    %indices of non-zero values of nifti image
    VI = find(V~=0);
    %path to AAL template
    parcelfile = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/aal_8mm_try5.nii.gz'; %load this from the provided codes folder
    %loading AAL labels
    load('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/AAL_labels.mat'); %load this from the provided codes folder
    %extracting AAL coordinates information
    K = nii.load(parcelfile);
    %sorting results in order to have strongest voxels at the top (positive t-values) or at the bottom (negative t-values)
    [VV2, II] = sort(VV,'descend');
    VI = VI(II);
    mni_coords = mni_coords(II,:);
    %final cell
    PD = cell(length(VV2),4);
    %getting AAL indices
    ROI = zeros(length(VI),1);
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
    writetable(PDn,['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/OldvsNew_MCS/SupplementaryTablesSourceClusters/Block_' num2str(block) '_freq_' num2str(freq) '.xlsx'],'Sheet',tonebytone) %printing excel file
end


%% COMPARING TONAL VS. ATONAL

%here we follow these steps:
%(1) t-tests
%(2) cluster permutations test
%(3) combining images with significant clusters

%% T-tests

subj_ato = zeros(23,27,23,526,71); %5D matrix with voxel coordinates and timepoints for all subjects  
subj_maj = zeros(23,27,23,526,71); %5D matrix with voxel coordinates and timepoints for all subjects 

freq = 1; %1 = 0.1-1Hz; 2 = 2-8Hz
contrast = 2; %set contrast to 1 2 or 3

for ii = 1:71 %over subjects
    if freq == 1
        dum2 = load_nii(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/FirLev_Block_1_freq_0.1_1.oat/subject' num2str(ii) '_wholebrain_first_level_BC_coape_subj_dir/tstat' num2str(contrast) '_8mm.nii.gz']);
        subj_ato(:,:,:,:,ii) = dum2.img; %4D matrix for each subject
        dum2 = load_nii(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/FirLev_Block_2_freq_0.1_1.oat/subject' num2str(ii) '_wholebrain_first_level_BC_coape_subj_dir/tstat' num2str(contrast) '_8mm.nii.gz']);
        subj_maj(:,:,:,:,ii) = dum2.img; %4D matrix for each subject
    else
        dum2 = load_nii(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/FirLev_Block_1_freq_2_8.oat/subject' num2str(ii) '_wholebrain_first_level_BC_N100_coape_real_subj_dir/tstat' num2str(contrast) '_8mm.nii.gz']);
        subj_ato(:,:,:,:,ii) = dum2.img;
        dum2 = load_nii(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/FirLev_Block_2_freq_2_8.oat/subject' num2str(ii) '_wholebrain_first_level_BC_N100_coape_real_subj_dir/tstat' num2str(contrast) '_8mm.nii.gz']);
        subj_maj(:,:,:,:,ii) = dum2.img;
    end
    disp(ii)
end

onset = [16 54 91 129 166]; %onset of tones in time samples

for jj = 1:length(onset) %over onset of tones
    if onset(jj) == 16 || 129 %add 36 or 37 to the onset to obtain the full time window for each tone
        x = 37;
    else
        x = 36;
    end
    time_window_maj = subj_maj(:,:,:,onset(jj):onset(jj) + x,:); %time window of tone x for major
    dpm1 = squeeze(mean(time_window_maj,4)); %average of tone x for major
    %the fourth dimension (time) is eliminated using 'squeeze' since it corresponds to only one point in time
    time_window_ato = subj_ato(:,:,:,onset(jj):onset(jj) + x,:); %time window of tone x for atonal
    dpm2 = squeeze(mean(time_window_ato,4)); %average of tone x for atonal
    
    %t-test
    sz = size(dpm1);
    T = zeros(sz(1),sz(2),sz(3));
    P = zeros(sz(1),sz(2),sz(3));
    for pp = 1:sz(1) %first dimension of voxel coordinates
        for yy = 1:sz(2) %second dimension of voxel coordinates
            for zz = 1:sz(3) %third dimension of voxel coordinates
                if dpm1(pp,yy,zz,1) ~= 0 %only for values that correspond to the brain (not the 'box' that MATLAB creates)
                    a = squeeze(dpm1(pp,yy,zz,:)); %major
                    b = squeeze(dpm2(pp,yy,zz,:)); %atonal                   
                    [h,p,ci,stats] = ttest(a,b); %actual t-test for tone x
                    T(pp,yy,zz) = stats.tstat; %t-value
                    P(pp,yy,zz) = 1-p; %p-value
                end
            end
        end
    end
    if freq == 1
        pathfreq = ['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/ContrastMajorAtonal/Freq_01_1/Contrast' num2str(contrast)];
        if ~exist(pathfreq,'dir')
            mkdir(pathfreq) %create a new folder for each contrast
        end
    else
        pathfreq = ['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/ContrastMajorAtonal/Freq_2_8/Contrast' num2str(contrast)];
        if ~exist(pathfreq,'dir')
            mkdir(pathfreq) %create a new folder for each contrast
        end
    end
    
    dum2 = load_nii(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/FirLev_Block_1_freq_0.1_1.oat/wholebrain_first_level_BC_coape_subj_group_level_everybody_randomise_c1_dir_tone_5/tstat1_gc1_8mm.nii.gz']);
    %use dum2 to index a file with only four dimensions (3 voxel dimensions and 1 subject dimension)
    dum2.img = T; %write the image
    save_nii(dum2,[pathfreq '/Tval_MajvsAto_tone_' num2str(jj) '.nii.gz']); %save the nifti image for the t-values
    dum2.img = P; %write the image
    save_nii(dum2,[pathfreq '/Pval_MajvsAto_tone_' num2str(jj) '.nii.gz']); %save the nifti image for the p-values
    disp(jj)
end

%% Cluster permutations test

%defining clusters on 3D brain voxel statistics

%loading p-values and t-values
clear DATA

freq = 1; %1 = 0.1-1Hz; 2 = 2-8Hz
contrast = 1; %set to 1 2 or 3
tone = 5; %set to 1 2 3 4 or 5

if freq == 1
    P = load_nii(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/ContrastMajorAtonal/Freq_01_1/Contrast' num2str(contrast) '/Pval_MajvsAto_tone_' num2str(tone) '.nii.gz']);
    T = load_nii(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/ContrastMajorAtonal/Freq_01_1/Contrast' num2str(contrast) '/Tval_MajvsAto_tone_' num2str(tone) '.nii.gz']);
else
    P = load_nii(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/ContrastMajorAtonal/Freq_2_8/Contrast' num2str(contrast) '/Pval_MajvsAto_tone_' num2str(tone) '.nii.gz']);
    T = load_nii(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/ContrastMajorAtonal/Freq_2_8/Contrast' num2str(contrast) '/Tval_MajvsAto_tone_' num2str(tone) '.nii.gz']);
end

%extracting matrix with statistics
P2 = P.img;
T2 = T.img;

%mask for non-0 voxels in brain imges (basically building a layout for actual brain voxels)
mask = zeros(size(T2,1),size(T2,2),size(T2,3));
mask(T2~=0) = 1; %assigning 1 when you have real brain voxels

%preparing data
data = T2;
data(data==0) = NaN; %removing non-brain voxels
P2(P2==0) = NaN; %removing non-brain voxels
data(P2<0.95) = NaN; %removing non-significant voxels

%tonal vs. atonal
%cond1>cond2 (so t-values are positive)
data1 = data;
data1(T2<0) = NaN; %removing negative t-values (so cond1>cond2)
data1(~isnan(data1)) = 1; %assigning 1 to significant voxels
DATA{1} = data1; %storing data

%atonal vs. tonal
%cond2>cond1 (so t-values are negative)
data2 = data;
% data2(P2<0.95) = NaN; %removing non-significant voxels (MORE STRICT THRESHOLD HERE FOR COND2>COND1)
data2(T2>0) = NaN; %removing positive t-values (so cond2>cond1)
data2(~isnan(data2)) = 1; %assigning 1 to significant voxels
DATA{2} = data2; %storing data

%getting MNI coordinates
%OBS! you may get a warning since the skeletonized image is not exactly in MNI space, but close enough
[ mni_coords, xform ] = osl_mnimask2mnicoords('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_T1_8mm_brain.nii.gz');

%preparation of information and data for the actual function
for ii = 1:2 %over directions of the contrast (cond1>cond2 and cond2>cond1)
    S = [];
    if freq == 1
        S.T = ['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/ContrastMajorAtonal/Freq_01_1/Contrast' num2str(contrast) '/Tval_MajvsAto_tone_' num2str(tone) '.nii.gz']; %path to image with t-values
        S.outdir = ['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/ContrastMajorAtonal/Freq_01_1/Contrast' num2str(contrast) '/Cluster_Permutation_test']; %output path
    else
        S.T = ['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/ContrastMajorAtonal/Freq_2_8/Contrast' num2str(contrast) '/Tval_MajvsAto_tone_' num2str(tone) '.nii.gz']; %path to image with t-values
        S.outdir = ['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/ContrastMajorAtonal/Freq_2_8/Contrast' num2str(contrast) '/Cluster_Permutation_test']; %output path
    end
    S.parcelfile = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/aal_8mm_try5.nii.gz';
    S.labels = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/AAL_labels.mat';
    S.MNIcoords = mni_coords; %MNI coordinates of 8mm MNI152T1 brain
    S.data = DATA{ii};
    S.mask = mask; %mask of the brain layout you have your results in
    S.permut = 10; %number of permutations for Monte Carlo simulation (MCS)
    S.clustmax = 1; %set 1 for only max cluster size of each permutation MCS (more strict); set 0 for every size of each cluster detected for each permutation MCS (less strict)
    S.permthresh = 0.001; %threhsold for MCS
    
    if ii == 1
        S.anal_name = ['Tval_MajvsAto_tone_' num2str(tone) '_Cond1vsCond2']; %name for the analysis (used to identify and save image and results)
    else
        S.anal_name = ['Tval_MajvsAto_tone_' num2str(tone) '_Cond2vsCond1']; %name for the analysis (used to identify and save image and results)
    end
    
    %actual function
    PP = BrainSources_MonteCarlosim_3D_LBPD_D(S);
end

%% Combining images (Cond1>Cond2 and Cond2>Cond1) with significant clusters

%here we obtain information about the brain regions forming the clusters
%the tables can be found in SUPPLEMENTARY MATERIALS

freq = 1; %1 = 0.1-1Hz; 2 = 2-8Hz
contrast = 1; %set to 1 2 or 3
tone = 4; %set to 1 2 3 4 or 5

if freq == 1
    %path to cluster permutation test images
    path = ['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/ContrastMajorAtonal/Freq_01_1/Contrast' num2str(contrast) '/Cluster_Permutation_test'];
else
    path = ['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/source/ContrastMajorAtonal/Freq_2_8/Contrast' num2str(contrast) '/Cluster_Permutation_test'];
end

%major vs. atonal (cond1>con2)
name1 = ['Tval_MajvsAto_tone_' num2str(tone) '_Cond1vsCond2_SignClust_1_Tvals.nii.gz'];

%atonal vs. tonal (cond2>cond1)
name2 = ['Tval_MajvsAto_tone_' num2str(tone) '_Cond2vsCond1_SignClust_1_Tvals.nii.gz'];

output = ['Tval_MajvsAto_tone_' num2str(tone) '_All_Clusters.nii.gz'];

cmd = ['fslmaths ' path '/' name1 ' -add ' path '/' name2 ' ' path '/' output]; %OBS! be careful with the spacing
system(cmd)

%after this, we created the final images using Workbench

%%


