%% Sleep Scoring - iEEG - Columbia
% For CUBFXX patients or any .mat file with data

eeglab;
pID = 'CUBF10';
nID = '0208';
cd(['/Volumes/Data Backup Epilepsy/',pID,'/',pID,'_',nID,'/']);
fname = ['/Volumes/Data Backup Epilepsy/',pID,'/',pID,'_',nID,'/',pID,'_',nID,'.mat']; 
load(fname);
load(['/Volumes/Data Backup Epilepsy/',pID,'/',pID,'_',nID,'/',pID,'_',nID,'_channel_labels.mat']);
fs = 200;
EEG = pop_importdata(fname);
EEG.data = data;
EEG = pop_editset(EEG,'srate',fs);
%EEG.chanlocs =  channel_labels;

% Bandpass filtering the signal
%l_freq = 1;
%h_freq = 40;
%EEG = pop_eegfiltnew(EEG,'locutoff',l_freq,'hicutoff',h_freq);

EEGcopy = EEG;

% Remove bad channels before average referencing the data
EEG = csc_eeg_plotter(EEG);
% h_channels = EEG.csc_hidden_channels;
% EEG.csc_hidden_channels = [];
% for i = 1:length(h_channels)
%     EEG.csc_hidden_channels(i) = EEG.chanlocs(i).labels;
% end
% 
% 
% EEG.csc_hidden_channels = (EEG.chanlocs(h_channels).labels)';
hidden_chans = EEG.csc_hidden_channels;
save([pID,'_',nID,'_hiddenchans.mat'],'hidden_chans');

EEG = pop_select(EEG,'nochannel',EEG.csc_hidden_channels);

%% Average reference the data to remove noise
EEG_avgref = EEG.data - mean(EEG.data);
EEG.data = EEG_avgref;

%% Evaluating the signal for sleep scoring

EEG = csc_eeg_plotter(EEG);

cd(['/Volumes/Data Backup Epilepsy/',pID,'/',pID,'_',nID,'/']);
pop_saveset(EEG,[pID,'_',nID,'_filtered']);

%%
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
    'epoch_length',     6           ,...  % window length in seconds
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

figure;
plot(mean(fft_SWA_DB));
EEG_fft = csc_eeg_plotter(EEG_fft);

%% After scoring 

load('sleepscore_checkedMB.mat');

% Sleep Epoch
for i = 1:length(event_data)
    if (event_data{i,3} == 1)||(event_data{i,3} == 2)||(event_data{i,3} == 3)
        event_val(i,1) = 0;
    else
        event_val(i,1) = 1;
    end  
end

% Finding the corresponding sleep stages
sleep = find(event_val == 0);
wake = find(event_val == 1);
D = diff([0,diff(sleep')==1,0]);
D1 = diff([0,diff(wake')==1,0]);

% Start and end points of a consecutive block of sleep
event_points_sleep(:,1) = sleep(D>0);
event_points_sleep(:,2) = sleep(D<0);

event_points_wake(:,1) = wake(D1>0);
event_points_wake(:,2) = wake(D1<0);


% If the wake between sleep is greater than 5 mins discard it else it
% can be considered as sleep vice versa for sleep 
sleep_dur = (event_points_sleep(:,2) - event_points_sleep(:,1))*0.5;
wake_dur = (event_points_wake(:,2) - event_points_wake(:,1))*0.5;


for i = 1:length(wake_dur)
    if wake_dur(i) < 5
        event_val(event_points_wake(i,1):event_points_wake(i,2)) = 0;
    end
end
     
ref_sleep = find(event_val == 0);
ref_wake = find(event_val == 1);

D = diff([0,diff(ref_sleep')==1,0]);
D1 = diff([0,diff(ref_wake')==1,0]);

event_points_sleep2(:,1) = ref_sleep(D>0);
event_points_sleep2(:,2) = ref_sleep(D<0);
event_points_wake2(:,1) = ref_wake(D1>0);
event_points_wake2(:,2) = ref_wake(D1<0);

% Sleep epoch in sample / time points
sleep_epoch = event_points_sleep2(1,1)*30:event_points_sleep2(end,2)*30;
wakeBS_epoch = event_points_wake2(1,1)*30:event_points_wake2(1,2)*30;
wakeAS_epoch = event_points_wake2(end,1)*30:event_points_wake2(end,2)*30;

EEG_sleep = EEG_fft.data(:,sleep_epoch);
EEG_wakeBS = EEG_fft.data(:,wakeBS_epoch);
EEG_wakeAS = EEG_fft.data(:,wakeAS_epoch);

save('CUBF10_0208_sleepData.mat','EEG_sleep');
save('CUBF10_0208_wakeBSData.mat','EEG_wakeBS');
save('CUBF10_0208_wakeASData.mat','EEG_wakeAS');

%% Step 3 - Calculate the FFT power for the different bands for the different sections of the night.
% Repeat the step 3 for each section of the night by selecting the
% appropriate set file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FFT time course - overview
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fname = 'CUBF10_0208_wakeBSData.mat';

EEG_fft = pop_importdata(fname);
EEG_fft.data = EEG_wakeBS;
fs = 200;
EEG_fft = pop_editset(EEG_fft,'srate',fs);

filtered_file = 'wakeBSData';

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

save([filtered_file '_fft_wholedata.mat'], 'fft_all', 'fft_bands','fft_SWA_2_4','fft_SWA_4_8','fft_Spindles','fft_SWA_DB',... %'fft_HF_slided', 'fft_SWA_slided', 
    'options', 'freq_range',... %'good_epochs' ,
    'filtered_file' , '-mat', '-v7.3');

disp('FFT data saved');
