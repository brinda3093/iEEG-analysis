
eeglab
EEG = csc_eeg_plotter(EEG);

load('/Volumes/Data Backup Epilepsy/CUBF10/CUBF10_0208/spikes_cubf208.mat'); % Respective spike file

EEG.csc_event_data = event_data;

EEG = pop_saveset(EEG);

EEG2 = EEG;

%% 

p_sleep = [790,  1534]*30*200;
p_wakeBS = [274 , 789]*30*200;
p_wakeAS = [1535 , 1603]*30*200;

for i = 1:length(event_data)
    if event_data{i,1} == 'event 1'
        spike_val(i,1) = floor(event_data{i,2}*200);
    end
end

sleep_ep = [];

for i = 1:length(event_data)
    if spike_val(i,1) > p_sleep(1,1) && spike_val(i,1) < p_sleep(1,2)
        sleep_ep = [sleep_ep , spike_val(i,1)];
    end
end


%% Spike Frequency



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

