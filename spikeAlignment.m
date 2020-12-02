
eeglab
EEG = csc_eeg_plotter(EEG);

load('/Volumes/Data Backup Epilepsy/CUBF10/CUBF10_0208/spikes_cubf208.mat'); % Respective spike file

EEG.csc_event_data = event_data;

EEG = pop_saveset(EEG);

EEG2 = EEG;

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




