%% Sleep Scoring - iEEG - Columbia
% For CUBFXX patients or any .mat file with data

eeglab;
pID = 'CUBF27';
nID = '0430';
cd(['/Volumes/Data Backup Epilepsy/',pID,'/',pID,'_',nID,'/']);
fname = ['/Volumes/Data Backup Epilepsy/',pID,'/',pID,'_',nID,'/',pID,'_',nID,'.mat']; 
load(fname);
load(['/Volumes/Data Backup Epilepsy/',pID,'/',pID,'_',nID,'/',pID,'_',nID,'_channel_labels.mat']);
fs = 200;
EEG = pop_importdata(fname);
EEG.data = data;%merged_matfile;
EEG = pop_editset(EEG,'srate',fs);
%EEG.chanlocs =  channel_labels;

% Bandpass filtering the signal
l_freq = 1;
h_freq = 40;
EEG = pop_eegfiltnew(EEG,'locutoff',l_freq,'hicutoff',h_freq);
%EEG = pop_cleanline(EEG);
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

% Removing bad sections of data using EEGLAB before sleep scoring
% Section added 6/3/2020
eegplot(EEG.data,...
  'srate', EEG.srate,...
  'eloc_file', EEG.chanlocs,...
  'spacing', 500,...
  'winlength', 30,...
  'dispchans', EEG.nbchan,...
  'command', 'disp(''No data rejected. Use pop_select for this.'')',...
  'butlabel', 'MARK', ...
  'events', EEG.event);

EEG = epi_log(@pop_select, EEG, 'nopoint', TMPREJ(:,1:2));
EEG.removed_stretches = TMPREJ(:,1:2);
EEG.filename = [pID,'_',nID];
EEG = pop_saveset(EEG, 'filename', strcat([EEG.filename '_HP_nobadstr.set']));
save([pID,nID,'_data_assess_sleep.mat'],'EEG','TMPREJ');



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

%Load the sleep scoring file
load('sleepscore_checkedMB.mat');
%%
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
%%
% Sleep epoch in sample / time points
%sleep_epoch_first = [event_data{event_points_sleep2(1,1),2}*200:event_data{event_points_sleep2(end,2),2}*200];
sleep_epoch = [event_data{event_points_sleep2(2,1),2}*200:event_data{event_points_sleep2(5,2),2}*200];
%wakeAS_epoch = [event_data{event_points_wake2(2,1),2}*200:event_data{event_points_wake2(2,2),2}*200];
%wakeBS_epoch = [event_data{event_points_wake2(1,1),2}*200:event_data{event_points_wake2(1,2),2}*200];
wakeBS_epoch = [event_data{event_points_wake2(1,1),2}*200:event_data{event_points_wake2(1,2),2}*200];
wakeAS_epoch = [event_data{event_points_wake2(2,1),2}*200:event_data{event_points_wake2(2,2),2}*200];
%sleep_epoch_last = [event_data{event_points_sleep2(5,1),2}*200:event_data{event_points_sleep2(5,2),2}*200];

wakeBS_epoch(1) = [];
%EEG_sleep_first = EEG_fft.data(:,sleep_epoch_first);
%EEG_sleep_last = EEG_fft.data(:,sleep_epoch_last);
%%
EEG_sleep = EEG_fft.data(:,sleep_epoch);

EEG_wakeBS = EEG_fft.data(:,wakeBS_epoch);
EEG_wakeAS = EEG_fft.data(:,wakeAS_epoch);
%%


%save('CUBF09_0123_sleepData_first.mat','EEG_sleep_first');
%save('CUBF09_0123_sleepData_last.mat','EEG_sleep_last');
save('CUBF09_0123_sleepData_full.mat','EEG_sleep');
save('CUBF10_0208_wakeBSData.mat','EEG_wakeBS');
save('CUBF10_0208_wakeASData.mat','EEG_wakeAS');

%% Step 3 - Calculate the FFT power for the different bands for the different sections of the night.
% Repeat the step 3 for each section of the night by selecting the
% appropriate set file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FFT time course - overview
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fname = 'CUBF09_0123_sleepData_full.mat';
load(fname);

EEG_fft = pop_importdata(fname);
EEG_fft.data = EEG_sleep;
fs = 200;
EEG_fft = pop_editset(EEG_fft,'srate',fs);

filtered_file = 'sleepData_full';

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

%% Spike Removal

eeglab
EEG = csc_eeg_plotter(EEG);

load('/Volumes/Data Backup Epilepsy/CUBF10/CUBF10_0208/spikes_cubf208.mat'); % Respective spike file

EEG.csc_event_data = event_data;

EEG = pop_saveset(EEG);

EEG2 = EEG;

%% Check if the spikes lie in the bad stretches.

EEG_nobdstr = pop_loadset();

removed_stretches = EEG_nobdstr.event;
for i = 1:length(event_data)
    spike_locs(i,1) = event_data{i,2}*200;
end

spike_new = spike_locs;
sub_latency = 0;

for j = 1:length(removed_stretches)
    i = 1;
    while i <= length(spike_locs)
        sub_latency = removed_stretches(j).duration;
        if removed_stretches(j).latency <= spike_locs(i,1)
            spike_new(i,1) = spike_new(i,1) - sub_latency;
            i = i+1;
        else
            i = i +1;
        end
    end
end
event_data_new = event_data;

for i = 1:length(event_data_new)
     event_data_new{i,2}= spike_new(i,1)/200;
end
%% Frequency of Spikes in Sleep vs Wake
% for others

% for i = 1:length(event_points_wake2)
%     x = [event_data{event_points_wake2(i,1),2},event_data{event_points_wake2(i,2),2}];
%     fr = 0;
%   
%         for k = 1:length(event_data_new)
%             
%             if x(1,1) <= event_data_new{k,2} && event_data_new{k,2} < (x(1,2))
%                 fr = fr +1;
%                 freq_wake(i,1) = fr*30/(x(1,2) - x(1,1) + 1);
%             end
%         end
% end

% for CUBF07_111
% 
% x = [19920,26550;
%     event_data{1945,2}, event_data{2039,2} ];
% for i = 1:length(x)
% fr = 0;
%   
%         for k = 1:length(event_data_new)
%             
%             if x(i,1) <= event_data_new{k,2} && event_data_new{k,2} < (x(i,2))
%                 fr = fr +1;
%                 freq_wake(i,1) = fr*30/(x(i,2) - x(i,1) + 1);
%             end
%         end
% end

for i = 1:length(event_points_sleep2)
    x = [event_data{event_points_sleep2(i,1),2},event_data{event_points_sleep2(i,2),2}];
    fr = 0;
  
        for k = 1:length(event_data_new)
            if x(1,1) <= event_data_new{k,2} && event_data_new{k,2} < (x(1,2))
                fr = fr +1;
                freq_sleep(i,1) = fr*30/(x(1,2) - x(1,1) + 1);
            end
        end
end
       
        
tot_sleep = mean(freq_sleep);
tot_wake = mean(freq_wake);

per_diff = (tot_wake - tot_sleep)/(tot_sleep+tot_wake)*100;


%% Spike frequency individual cycles
%Wake 


    x1 = [event_data{event_points_wake(1,2),2} - 3600, event_data{event_points_wake2(1,2),2}]; % Cycle 1 
    x2 = [event_data{event_points_wake2(end,1),2} ,event_data{event_points_wake2(end,1),2}+3600]; % Cycle 2
    fr = 0;
  
        for k = 1:length(event_data_new)
            if x1(1,1) <= event_data_new{k,2} && event_data_new{k,2} < (x1(1,2))
                fr = fr +1;
                freq_wake_1(1,1) = fr*30/(x1(1,2) - x1(1,1) + 1);
            end
        end
       
        fr= 0;
         for k = 1:length(event_data_new)
            if x2(1,1) <= event_data_new{k,2} && event_data_new{k,2} < (x2(1,2))
                fr = fr +1;
                freq_wake_2(1,1) = fr*30/(x2(1,2) - x2(1,1) + 1);
            end
        end
        


SpikeFreq_W = 100*(freq_wake(6,1) - freq_wake(2,1))/tot_wake;

%%
% Sleep Frequency

    x1 = [event_data{event_points_sleep2(1,1),2}, event_data{event_points_sleep2(1,1),2}+ 3600]; % Cycle 1 
    x2 = [event_data{event_points_sleep2(end,2),2} - 3600 ,event_data{event_points_sleep2(end,2),2}]; % Cycle 2
    fr = 0;
  
        for k = 1:length(event_data_new)
            if x1(1,1) <= event_data_new{k,2} && event_data_new{k,2} < (x1(1,2))
                fr = fr +1;
                freq_sleep_1(1,1) = fr*30/(x1(1,2) - x1(1,1) + 1);
            end
        end
%%        
        fr= 0;
         for k = 1:length(event_data_new)
            if x2(1,1) <= event_data_new{k,2} && event_data_new{k,2} < (x2(1,2))
                fr = fr +1;
                freq_sleep_2(1,1) = fr*30/(x2(1,2) - x2(1,1) + 1);
            end
        end
        
        
SpikeFreq_S = 100*(freq_sleep_2(1,1) - freq_sleep_1(1,1))/tot_sleep;
%SpikeFreq_S = 100*(freq_sleep(5,1)-freq_sleep(2,1))/tot_sleep;
%% Remove the 1s before and after the marked spike 

spiketime  = [];

for ispike = 1:length(EEG.csc_event_data)
    spiketime(ispike,1) = round(EEG.csc_event_data{ispike,2}*EEG.srate);
end


% Spikes to delete - 1s before and after the marked event on EEGplotter
del_spikes = [];
del_spikes(:,1) = spiketime - EEG.srate;
del_spikes(:,2) = spiketime + EEG.srate;

% 
EEG_spikes = pop_select(EEG,'point',del_spikes); % Spikes
EEG_nospikes = pop_select(EEG,'nopoint',del_spikes); % Data without the event markers


EEG_spikes = csc_eeg_plotter(EEG_spikes); % Helps visualize deleted spikes
EEG_nospikes = csc_eeg_plotter(EEG_nospikes); % Visualizing data without spikes

EEG_nospikes = pop_saveset(EEG_nospikes);


%% Track sleep 
% This section is only to be used after removing the spikes and then
% tracking the sleep and wake epoch changes accordingly

% To track sleep after removing the spikes, check the EEG.event for the filtered file after spike removal 
% You will need the event_points_sleep2 and event_points_wake2 from the
% sleep marking

removed_spikes = EEG.event;

% Latency and Durations are in samples
for i = 1:length(removed_spikes)
    spikeEvents(i,1:2) = [removed_spikes(i).latency , [removed_spikes(i).latency + removed_spikes(i).duration]];
end

% Load the sleep epochs as marked


%load('sleepscore_BS_withSpikes_6142020.mat'); % CUBF09_123
sleep_epoch1_start = event_data{event_points_sleep2(2,1),2}*200;
%sleep_epoch2_start = event_data{event_points_sleep2(5,1),2}*200;

sleep_epoch1_end = event_data{event_points_sleep2(5,2),2}*200;
%sleep_epoch2_end = event_data{event_points_sleep2(5,2),2}*200;



%% New Sleep Epoch Values
% Epoch 1
duration_delete.before_sleepep1 = 0;
for i = 1:length(spikeEvents)
    duration_delete.before_sleepep1 = duration_delete.before_sleepep1 + removed_spikes(i).duration;
    if spikeEvents(i,1) >= sleep_epoch1_start
        spike_sleep1_start = i;
        break;
    end
end

new_sleep_epoch1_start = sleep_epoch1_start - duration_delete.before_sleepep1;

duration_delete.during_sleepep1 = 0;
for i = spike_sleep1_start:length(spikeEvents)
    duration_delete.during_sleepep1 = duration_delete.during_sleepep1 + removed_spikes(i).duration;
    if spikeEvents(i,1) >= sleep_epoch1_end
        spike_sleep1_end = i;
        break;
    end
end
del_dur =  duration_delete.during_sleepep1 + duration_delete.before_sleepep1;
new_sleep_epoch1_end = sleep_epoch1_end - del_dur;

% Epoch 2
duration_delete.before_sleepep2 = 0;
for i = spike_sleep1_end:length(spikeEvents)
    duration_delete.before_sleepep2 = duration_delete.before_sleepep2 + removed_spikes(i).duration;
    if spikeEvents(i,1) >= sleep_epoch2_start
        spike_sleep2_start = i;
        break;
    end
end

del_dur =  duration_delete.during_sleepep1 + duration_delete.before_sleepep1 + duration_delete.before_sleepep2;
new_sleep_epoch2_start = sleep_epoch2_start - del_dur;

duration_delete.during_sleepep2 = 0;
for i = spike_sleep2_start:length(spikeEvents)
    duration_delete.during_sleepep2 = duration_delete.during_sleepep2 + removed_spikes(i).duration;
    if spikeEvents(i,1) >= sleep_epoch1_end
        spike_sleep2_end = i;
        break;
    end
end
del_dur =  duration_delete.during_sleepep1 + duration_delete.before_sleepep1 + duration_delete.before_sleepep2+duration_delete.during_sleepep2;
new_sleep_epoch2_end = sleep_epoch2_end - del_dur;

%% Wake Epochs

wake_epoch1_start = event_data{event_points_wake2(1,1),2}*200;
wake_epoch2_start = event_data{event_points_wake2(2,1),2}*200;

wake_epoch1_end = event_data{event_points_wake2(1,2),2}*200;
wake_epoch2_end = event_data{event_points_wake2(2,2),2}*200;
del_dur = 0;
% Epoch1
% Epoch 1
duration_delete.before_wakeep1 = 0;
for i = 1:length(spikeEvents)
    duration_delete.before_wakeep1 = duration_delete.before_wakeep1 + removed_spikes(i).duration;
    if spikeEvents(i,1) >= wake_epoch1_start
        spike_wake1_start = i;
        break;
    end
end
del_dur = duration_delete.before_wakeep1;
new_wake_epoch1_start = wake_epoch1_start - del_dur;

duration_delete.during_wakeep1 = 0;
for i = spike_wake1_start:length(spikeEvents)
    duration_delete.during_wakeep1 = duration_delete.during_wakeep1 + removed_spikes(i).duration;
    if spikeEvents(i,1) >= wake_epoch1_end
        spike_wake1_end = i;
        break;
    end
end
del_dur = del_dur + duration_delete.during_wakeep1;
new_wake_epoch1_end = wake_epoch1_end - del_dur;

% Epoch 2
duration_delete.before_wakeep2 = 0;
for i = spike_wake1_end:length(spikeEvents)
    duration_delete.before_wakeep2 = duration_delete.before_wakeep2 + removed_spikes(i).duration;
    if spikeEvents(i,1) >= wake_epoch2_start
        spike_wake2_start = i;
        break;
    end
end
del_dur = del_dur + duration_delete.before_wakeep2;
new_wake_epoch2_start = wake_epoch2_start - del_dur;

duration_delete.during_wakeep2 = 0;
for i = spike_wake2_start:length(spikeEvents)
    duration_delete.during_wakeep2 = duration_delete.during_wakeep2 + removed_spikes(i).duration;
    if spikeEvents(i,1) >= wake_epoch1_end
        spike_wake2_end = i;
        break;
    end
end
del_dur = del_dur + duration_delete.during_wakeep2;
new_wake_epoch2_end = wake_epoch2_end - del_dur;

%% Defining the EEG values for new epochs defined
EEG_nospikes_sleep = EEG.data(:,new_sleep_epoch1_start:new_sleep_epoch1_end);
%EEG_nospikes_sleep_first = EEG.data(:,new_sleep_epoch1_start:new_sleep_epoch1_end);
%EEG_nospikes_sleep_last = EEG.data(:,new_sleep_epoch2_start:new_sleep_epoch2_end);

EEG_nospikes_wakeBS = EEG.data(:,1:new_wake_epoch1_end);%new_wake_epoch1_start:new_wake_epoch1_end);
EEG_nospikes_wakeAS = EEG.data(:,new_wake_epoch2_start:length(EEG.data));%new_wake_epoch2_end);

% save('CUBF09_0123_sleepFirst_WOSpikes.mat','EEG_nospikes_sleep_first');
% save('CUBF09_0123_sleepLast_WOSpikes.mat','EEG_nospikes_sleep_last');
save('CUBF09_0123_sleep_full_WOSpikes.mat','EEG_nospikes_sleep');
save('CUBF10_0208_wakeBS_WOSpikes.mat','EEG_nospikes_wakeBS');
save('CUBF10_0208_wakeAS_WOSpikes.mat','EEG_nospikes_wakeAS');


%% Step 3 - Calculate the FFT power for the different bands for the different sections of the night.
% Repeat the step 3 for each section of the night by selecting the
% appropriate set file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FFT time course - overview
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fname = 'CUBF09_0123_sleep_full_WOSpikes.mat';
load(fname);

EEG_fft = pop_importdata(fname);
EEG_fft.data = EEG_nospikes_sleep;
fs = 200;
EEG_fft = pop_editset(EEG_fft,'srate',fs);

filtered_file = 'sleepData_full';

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


