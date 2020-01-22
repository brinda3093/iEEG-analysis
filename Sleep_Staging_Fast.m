% Script for analyzing the night faster

clear;

eeglab;

EEG = pop_importdata();

LOW_CUTOFF = 0.5;
HIGH_CUTOFF = [];

EEG1 = pop_eegfiltnew(EEG, LOW_CUTOFF, HIGH_CUTOFF, [], 0, [], 0);
EEG1 = csc_eeg_plotter(EEG1);
EEG1.badchannels = EEG1.csc_hidden_channels;
EEG1 = pop_select(EEG1, 'nochannel', EEG1.badchannels); % remove

eegplot(EEG1.data,...
  'srate', EEG1.srate,...
  'eloc_file', EEG1.chanlocs,...
  'spacing', 5000,...
  'winlength', 300,...
  'dispchans', EEG1.nbchan,...
  'command', 'disp(''No data rejected. Use pop_select for this.'')',...
  'butlabel', 'MARK', ...
  'events', EEG1.event);
%% 

% Remove bad stretches of data
EEG = epi_log(@pop_select, EEG1, 'nopoint', TMPREJ(:,1:2));

EEG.removed_stretches = TMPREJ(:,1:2);

Data1 = double(EEG.data);
Data = mean(Data1,1);
El_number = 1;
El_name{1,1} = '1'; 

fs = 200;
x = 'trialFile';

% Create a directory in any of the working folders
mkdir('/Users/bsevak/Documents/',x);

% Change the directory according to the folder created
cd(['/Users/bsevak/Documents/',x]);
save([x,'_data.mat'],'Data','El_number','El_name','fs','-v7.3');

% Directory minus the final folder (here /trialFile)
cd('/Users/bsevak/Documents/');
global subject_id z;
subject_id  = {'trialFile'};
z = 1;

Assess_the_sleep_NST();