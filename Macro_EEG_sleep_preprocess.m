% Once night is selected, let's do these steps:

%%
% Import in EEGlab

% there you do command line or gui, as you like - need to have only 1
% variable in .mat file for GUI to work

eeglab


%%
% take out bad channels - plotter will automatically save bad channels you mask out in EEG
% structure in workspace but not on computer

EEG = csc_eeg_plotter(EEG);

EEG = pop_saveset(EEG, 'filename', EEG.filename);

EEG = pop_select(EEG, 'nochannel', EEG.csc_hidden_channels );

EEG = pop_saveset(EEG, 'filename', strcat([EEG.filename(1:end-4) '_nobadchan.set']));


%%
% take out bad epochs

LOW_CUTOFF = 0.5;
HIGH_CUTOFF = [];

EEG1 = pop_eegfiltnew(EEG, LOW_CUTOFF, HIGH_CUTOFF, [], 0, [], 0)

eegplot(EEG1.data,...
  'srate', EEG1.srate,...
  'eloc_file', EEG1.chanlocs,...
  'spacing', 110,...
  'winlength', 120,...
  'dispchans', EEG1.nbchan,...
  'command', 'disp(''No data rejected. Use pop_select for this.'')',...
  'butlabel', 'MARK', ...
  'events', EEG1.event);

% eegplot will now save the stretches to be rejected in a workspace
% variable called TMPREJ. 

% Remove the marked stretches
% EEG = epi_log(@pop_select, EEG, 'nopoint', [120*EEG.srate, EEG.pnts]);
EEG = epi_log(@pop_select, EEG1, 'nopoint', TMPREJ(:,1:2));

%Let's keep in mind what we removed
EEG.removed_stretches = TMPREJ(:,1:2);


EEG = pop_saveset(EEG, 'filename', strcat([EEG.filename(1:end-4) '_HP_nobadstr.set']));

%%
% run classifier again
% using popp_select, select 1st and last cycles as the epoch from 1st to
% last black dot marked by classifier
Assess_the_sleep_modified()
global db_extracted;
db_sleep = [];
db_sleep_x = [];
db_wake = [];
db_wake_x = [];
for i = 1:length(db_extracted)
    if (db_extracted(i) > prctile(db_extracted,90))
        db_sleep = [db_sleep db_extracted(i)];
        db_sleep_x = [db_sleep_x i];
    elseif (db_extracted(i) < prctile(db_extracted, 5))
        db_wake = [db_wake db_extracted(i)];
        db_wake_x = [db_wake_x i];
    end
end




EEG_cycle = pop_select(EEG, 'time', [db_sleep_x(1)*30 db_sleep_x(end)*30]);
EEG_fft = pop_saveset(EEG_cycle, 'filename', strcat('C1_0710_fft.set'));


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
[fft_bands, freq_bands, freq] = csc_calculate_freq_bands_last(fft_all, freq_range, options);
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
Spectra_Plotter(fft_all, good_epochs, options);

save([filtered_file(1:end-4) '_fft_wholedata.mat'], 'fft_all', 'fft_bands', 'fft_HF_slided', 'fft_SWA_slided', 'fft_HF_threshold', 'fft_SWA_threshold', 'options', 'freq_range', 'good_epochs' , 'good_epochs_thresholds', 'filtered_file' , '-mat', '-v7.3');


