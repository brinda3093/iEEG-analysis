%% Assessing Sleep, Selecting Cycles to be analyzed, Calculating FFT for Analysis
% Last Edited: 1-8-2020; Brinda Sevak

% The following script can be used to Assess the sleep which helps classify
% the given night of data based on the output of the classifier after the
% bad channels and bad stretches of the data are removed. Further we can
% select the sections of the night: 
% a) Sleep 
% b) Wake Before Sleep & 
% c) Wake After Sleep 
% using EEGLAB which are to be used for further analysis and calculations 

% The script is divided into the following steps
% Step 1.1
% Step 1.2
% Step 1.3
% Step 2.1
% Step 2.2
% Step 2.3

 
%% Step 1.1
% Import the data to remove the bad sections of the data

eeglab;

Patient = 'Patient_2';
Night = 'Night_3';
pID = 'P2';
nID = 'N3';

filename = [pID,nID];
cd(['/Users/bsevak/Documents/Merged Data_BF/Merged_Data/']);
load([Patient,'/',Night,'/',Patient,'_',Night,'_Transpose.mat']);
EEG = pop_importdata([Patient,'/',Night,'/',Patient,'_',Night,'_Transpose.mat']);
fs = 200;
EEG.data = EEGdata;
EEG = pop_editset(EEG, 'srate',fs); 

EEG.filename = filename;

%% Step 1.2
% Filtering the data 
cd(['/Users/bsevak/Documents/Merged Data_BF/Merged_Data/',Patient,'/',Night,'/']);

LOW_CUTOFF = 0.5;
HIGH_CUTOFF = [];

EEG1 = pop_eegfiltnew(EEG, LOW_CUTOFF, HIGH_CUTOFF, [], 0, [], 0);

EEG1 = csc_eeg_plotter(EEG1);

EEG1.badchannels = EEG1.csc_hidden_channels;
EEG1 = pop_select(EEG1, 'nochannel', EEG1.badchannels); % remove

%
eegplot(EEG1.data,...
  'srate', EEG1.srate,...
  'eloc_file', EEG1.chanlocs,...
  'spacing', 5000,...
  'winlength', 300,...
  'dispchans', EEG1.nbchan,...
  'command', 'disp(''No data rejected. Use pop_select for this.'')',...
  'butlabel', 'MARK', ...
  'events', EEG1.event);
%% Step 1.3 - 

EEG = epi_log(@pop_select, EEG1, 'nopoint', TMPREJ(:,1:2));

EEG.removed_stretches = TMPREJ(:,1:2);

EEG = pop_saveset(EEG, 'filename', strcat([EEG.filename '_HP_nobadstr.set']));
save([pID,nID,'_data_assess_sleep.mat'],'EEG');

load(['/Users/bsevak/Documents/Merged Data_BF/Merged_Data/',Patient,'/',Night,'/',Patient,'_',Night,'_200Hz_resampled.mat']);

Data = double(EEG.data);
El_number = double(1:length(Data(:,1)));
El_name = {header_ns3(:).Label};

fs = 200;
cd(['/Users/bsevak/Documents/Merged Data_BF/Merged_Data/',Patient,'/',Night,'/',pID,nID]);
save([pID,nID,'_data.mat'],'Data','El_number','El_name','fs','-v7.3');


%% Step 2.1 - Classifying into sleep-wake based on Delta-Beta ratio

cd(['/Users/bsevak/Documents/Merged Data_BF/Merged_Data/',Patient,'/',Night,'/'])
global subject_id z;
subject_id  = {[pID,nID]};
z = 1;

Assess_the_sleep_NST();
global Db;
db_sleep = [];
db_sleep_x = [];
db_wake = [];
db_wake_x = [];
for i = 1:length(Db)
    if (Db(i) > prctile(Db,90))
        db_sleep = [db_sleep Db(i)];
        db_sleep_x = [db_sleep_x i];
    elseif (Db(i) < prctile(Db, 5))
        db_wake = [db_wake Db(i)];
        db_wake_x = [db_wake_x i];
    end
end

save(['Db_',pID,nID,'.mat'],'Db');
disp('Move to step 2.2 to select the section of wake or sleep you want to analyze based on the CLassifier output');

%% Step 2.2 Select the valid night for sleep from the Delta Beta ratio on the EEGLAB

cd(['/Users/bsevak/Documents/Merged Data_BF/Merged_Data/',Patient,'/',Night,'/'])

type = 'sleep';
eeglab; 
    EEG = pop_importdata(['Db_',pID,nID,'.mat']);
    EEG.srate = 1;
    EEG.data = Db;
    EEG.xmax = length(EEG.data);
    disp('Mark the sleep cycle to be analyzed');
    eegplot(EEG.data,...
        'srate', EEG.srate,...
        'eloc_file', EEG.chanlocs,...
        'spacing', 110,...
        'winlength', length(EEG.data),...
        'dispchans', EEG.nbchan,...
        'command', 'disp(''No data rejected. Use pop_select for this.'')',...
        'butlabel', 'MARK', ...
        'events', EEG.event);
    waitfor(gcf);
    
sleep_epoch = data_epoch(EEG, pID, nID, type, Db, TMPREJ);

%% Step 2.3 Select the valid wake before sleep from the Delta Beta ratio on the EEGLAB
% Skip if no data is available

eeglab;
type = 'WakeBS';
EEG = pop_importdata(['Db_',pID,nID,'.mat']);
    EEG.srate = 1;
    EEG.data = Db;
    EEG.xmax = length(EEG.data);
    disp('Mark the Wake before sleep to be analyzed');
    eegplot(EEG.data,...
        'srate', EEG.srate,...
        'eloc_file', EEG.chanlocs,...
        'spacing', 110,...
        'winlength', length(EEG.data),...
        'dispchans', EEG.nbchan,...
        'command', 'disp(''No data rejected. Use pop_select for this.'')',...
        'butlabel', 'MARK', ...
        'events', EEG.event);
    waitfor(gcf);

wakeBS_epoch = data_epoch(EEG, pID, nID, type, Db, TMPREJ);


%% Step 2.4 Select the valid wake after sleep from the Delta Beta ratio on the EEGLAB
% Skip if no data is available

eeglab;
EEG = pop_importdata(['Db_',pID,nID,'.mat']);
    EEG.srate = 1;
    EEG.data = Db;
    EEG.xmax = length(EEG.data);
    disp('Mark the Wake after sleep to be analyzed');
    eegplot(EEG.data,...
        'srate', EEG.srate,...
        'eloc_file', EEG.chanlocs,...
        'spacing', 110,...
        'winlength', length(EEG.data),...
        'dispchans', EEG.nbchan,...
        'command', 'disp(''No data rejected. Use pop_select for this.'')',...
        'butlabel', 'MARK', ...
        'events', EEG.event);
    waitfor(gcf);

wakeAS_epoch = data_epoch(EEG, pID, nID, type, Db, TMPREJ);
%% Step 2.6 Saving the sleep, wakeBS and wakeAS data for FFT Calculations
%  Sleep Data for FFT calculation

load([pID,nID,'_data_assess_sleep.mat']);

Data = EEG.data;
data_fft = Data(:,sleep_epoch(1)*30*200:sleep_epoch(end)*30*200);

save([pID,nID,'sleep_data_fft.mat'],'data_fft');

eeglab;
EEG = pop_importdata([pID,nID,'sleep_data_fft.mat']);
EEG.data = data_fft;
EEG.srate = 200; %Sampling rate

EEG = pop_saveset(EEG, 'filename', [pID,nID,'_sleepData']);

%% WakeBS Data for FFT Calculation

load([pID,nID,'_data_assess_sleep.mat']);

Data = EEG.data;
data_fft = Data(:,wakeBS_epoch(1)*30*200:wakeBS_epoch(end)*30*200);

save([pID,nID,'wakeBS_data_fft.mat'],'data_fft');

eeglab;
EEG = pop_importdata([pID,nID,'wakeBS_data_fft.mat']);
EEG.data = data_fft;
EEG.srate = 200; %Sampling rate

EEG = pop_saveset(EEG, 'filename', [pID,nID,'_wakeBSData']);

%% Wake AS Data for FFT Calculation

load([pID,nID,'_data_assess_sleep.mat']);

Data = EEG.data;
data_fft = Data(:,wakeAS_epoch(1)*30*200:wakeAS_epoch(end)*30*200);

save([pID,nID,'wakeAS_data_fft.mat'],'data_fft');

eeglab;
EEG = pop_importdata([pID,nID,'wakeAS_data_fft.mat']);
EEG.data = data_fft;
EEG.srate = 200; %Sampling rate

EEG = pop_saveset(EEG, 'filename', [pID,nID,'_wakeASData']);

%%
% FFT time course - overview
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EEG_fft = pop_loadset();
filtered_file = EEG_fft.filename;

% TODO: Should we give this file a more meaningful name?
[path, name, ~] = fileparts(filtered_file);
mat_path = path;
mat_name = [name '.mat'];
options = struct(...
    'save_file',        0            ,... % TODO: What does this do?
    'save_name',        mat_name     ,...
    'save_path',        mat_path     , ...
    'epoch_length',     6            ,...  % window length in seconds
    'freq_limit',       240          ,...  % number of bands - to go till 40Hz
    'ylimitmax',        1            , ... % TODO: What does this do?
    'fft_bands',        1)           ;     % TODO: What does this do?

% Calculate the FFT
% Stores output variables and input options in specified file.

[fft_all, freq_range] = csc_FFT_last(EEG_fft, options);
% freq_range = 240 equally-spaced frequencies from 0 to 40Hz.
% fft_all = power in each channel, in each frequency given by 'freq_range',
% for each NON-OVERLAPPING time window.

% "Calculate bands' ffts".
% Stores output variables and input options in specified file.
[fft_bands, freq_bands, freq] = csc_calculate_freq_bands_last_modified(fft_all, freq_range, options);
% freq_bands = Internally defined bands, presumably not data-dependent.
% freq = strings with common names (e.g. "theta") for each of these bands.
% fft_bands = power in each of these bands, in each channel, at each time
% window.

fft_bands = fft_bands(1:8,:,:); % remove channel 10 which is zero;

% Visually identify the bad epochs using thresholding of the band(s) of interest (2nd parameters)
% This uses power in the delta and gamma bands.
% The automatic threshold is set so that 1% of epochs are rejected.
% Stores output variables and input options in specified file.
[channel_thresholds, fft_SWA_threshold, fft_HF_threshold] = ...
    csc_artifact_detection_fft_last(fft_bands       , ...
    [1 6]           , ...
    'semi_automatic', ...
    options);

% fft_SWA_threshold = the power of delta in each channel at each time
% will now also display beta/delta ratio as 2nd band.

fft_bands(:,:,7) = fft_bands(:,:,1)./fft_bands(:,:,5); % delta/beta ratio;

[channel_thresholds, fft_SWA_threshold, fft_HF_threshold] = ...
    csc_artifact_detection_fft_last(fft_bands       , ...
    [1 7]           , ...
    'semi_automatic', ...
    options);

% "calculate the percentage good epochs and good fft with sliding mean"
% Stores output variables and input options in specified file.
good_epochs = calculate_data_saved(fft_bands       , ...
    fft_SWA_slided  , ...
    fft_HF_slided   , ...
    options);

% Stores output variables and input options in specified file.
good_epochs_thresholds = calculate_data_saved_thresholds(fft_bands          , ...
    fft_SWA_threshold  , ...
    fft_HF_threshold   , ...
    options);

% "plot the spectra"
Spectra_Plotter(fft_all,good_epochs, options);

save([filtered_file(1:end-4) '_fft_wholedata.mat'], 'fft_all', 'fft_bands',... %'fft_HF_slided', 'fft_SWA_slided', 
    'fft_HF_threshold', 'fft_SWA_threshold', 'options', 'freq_range',... %'good_epochs' ,
    'good_epochs_thresholds', 'filtered_file' , '-mat', '-v7.3');


