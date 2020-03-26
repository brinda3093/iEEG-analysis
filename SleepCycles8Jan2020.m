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
% Step 1.1 - Importing the data into EEGLAB 
% Step 1.2 - Filtering the data followed by bad channel rejection in csc_eeg_plotter to remove
% and using EEGLab to remove bad stretches of the data
% Step 1.3 - Saving the filtered data in .mat and .set files for future
% use
% Step 2.1 - Classifying the nights data into sleep-wake states using the
% classifier provided by Vaclav Kremen - Mayo Clinic and modified to fit
% the data (Assess_the_sleep_NST). Saving the Delta-Beta ratio as obtained
% from the classifier
% Step 2.2 - Using the Delta-Beta Ratio obtained from the classifier to
% further select the sections of the sleep from the night to be analyzed.
% Step 2.3  - Using the Delta-Beta Ratio obtained from the classifier to
% further select the sections of the wake before sleep from the night to be analyzed.
% Step 2.4  - Using the Delta-Beta Ratio obtained from the classifier to
% further select the sections of the wake after sleep from the night to be analyzed.
% Step 2.5 -  Saving the sleep, wakeBS and wakeAS data for FFT Calculations
% Step 3 - Calculate the FFT power for the different bands for the different sections of the night.

% Repeat the step 3 for each section of the night by selecting the
% appropriate set file

% Functions required with the script are data_epochs.m, Assess_the_sleep_NST.m eeglab
% The following script saves the .mat files of the FFT data in the end
% for each electrode which can be used for the analysis. 
%% Step 1.1
% Import the data to remove the bad sections of the data
clear;

eeglab;

Patient = 'Patient_6';
Night = 'Night_5';
pID = 'P6';
nID = 'N5';

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
HIGH_CUTOFF = 40;

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

cd(['/Users/bsevak/Documents/Merged Data_BF/Merged_Data/',Patient,'/',Night,'/']);
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
    EEG = pop_importdata();%['Db_',pID,nID,'.mat']);
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
    
%sleep_epoch = data_epoch(EEG, pID, nID, type, Db, TMPREJ);
    
sleep_epoch = data_epoch(EEG, x, type, Db, TMPREJ);
%% Step 2.3 Select the valid wake before sleep from the Delta Beta ratio on the EEGLAB
% Skip if no data is available

eeglab;
type = 'WakeBS';
EEG = pop_importdata(['Db_',x,'.mat']);
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

%wakeBS_epoch = data_epoch(EEG, pID, nID, type, Db, TMPREJ);
wakeBS_epoch = data_epoch(EEG, x, type, Db, TMPREJ);

%% Step 2.4 Select the valid wake after sleep from the Delta Beta ratio on the EEGLAB
% Skip if no data is available

eeglab;
type = 'WakeAS';
EEG = pop_importdata(['Db_',x,'.mat']);
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

%wakeAS_epoch = data_epoch(EEG, pID, nID, type, Db, TMPREJ);
wakeAS_epoch = data_epoch(EEG,x, type, Db, TMPREJ);
%% Step 2.5 Saving the sleep, wakeBS and wakeAS data for FFT Calculations
%  Sleep Data for FFT calculation

load([pID,nID,'_data_assess_sleep.mat']);

%Data = EEG1.data;
Data = EEG_fft.data;
last = sleep_epoch(end);
data_fft = Data(:,sleep_epoch(1)*30*200:last*30*200);

%save([pID,nID,'sleep_data_fft.mat'],'data_fft');
save([x,'sleep_data_fft.mat'],'data_fft');

eeglab;
%EEG = pop_importdata([pID,nID,'sleep_data_fft.mat']);
EEG = pop_importdata([x,'sleep_data_fft.mat']);
EEG.data = data_fft;
EEG.srate = 400; %Sampling rate

%EEG = pop_saveset(EEG, 'filename', [pID,nID,'_sleepData']);
EEG = pop_saveset(EEG, 'filename', [x,'_sleepData']);
%% WakeBS Data for FFT Calculation

load([pID,nID,'_data_assess_sleep.mat']);

%Data = EEG1.data;
Data = EEG_fft.data;
last = wakeBS_epoch(end);
data_fft = Data(:,wakeBS_epoch(1)*30*200:last*30*200);

%save([pID,nID,'wakeBS_data_fft.mat'],'data_fft');
save([x,'wakeBS_data_fft.mat'],'data_fft');

eeglab;
EEG = pop_importdata([x,'wakeBS_data_fft.mat']);
EEG.data = data_fft;
EEG.srate = 400; %Sampling rate

%EEG = pop_saveset(EEG, 'filename', [pID,nID,'_wakeBSData']);
EEG = pop_saveset(EEG, 'filename', [x,'_wakeBSData']);
%% Wake AS Data for FFT Calculation

load([pID,nID,'_data_assess_sleep.mat']);

%Data = EEG1.data;
Data = EEG_fft.data;
last = wakeAS_epoch(end);
data_fft = Data(:,wakeAS_epoch(1)*30*200:last*30*200);

%save([pID,nID,'wakeAS_data_fft.mat'],'data_fft');
save([x,'wakeAS_data_fft.mat'],'data_fft');

eeglab;
EEG = pop_importdata([x,'wakeAS_data_fft.mat']);
EEG.data = data_fft;
EEG.srate = 400; %Sampling rate

%EEG = pop_saveset(EEG, 'filename', [pID,nID,'_wakeASData']);
EEG = pop_saveset(EEG, 'filename', [x,'_wakeASData']);
%% Step 3 - Calculate the FFT power for the different bands for the different sections of the night.
% Repeat the step 3 for each section of the night by selecting the
% appropriate set file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

fft_SWA_2_4 = fft_bands(:,:,1);
fft_SWA_4_8 = fft_bands(:,:,2);
fft_Spindles = fft_bands(:,:,5);
fft_SWA_DB = fft_bands(:,:,3)./fft_bands(:,:,6);

save([filtered_file(1:end-4) '_fft_wholedata.mat'], 'fft_all', 'fft_bands','fft_SWA_2_4','fft_SWA_4_8','fft_Spindles','fft_SWA_DB',... %'fft_HF_slided', 'fft_SWA_slided', 
    'options', 'freq_range',... %'good_epochs' ,
    'filtered_file' , '-mat', '-v7.3');

disp('FFT data saved');





