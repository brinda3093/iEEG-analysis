%% .mat to .set conversion for using csc_eeg_plotter()

% Selecting the directory where the files are supposed to be stored
cd('/Users/bsevak/Documents/Merged Data_BF/Merged_Data/');

% Path of the excel file containing the list of the files/folders in the
% Patient and Night format which are to be converted to .set files from
% .mat files
xfile = readtable('/Users/bsevak/Documents/convertToSet.xlsx');

% For loop across all the files in the excel file
for i = 37:length(xfile.Patient) 
    
    patient = char(xfile.Patient(i)); % Patient ID
    night = char(xfile.Night(i)); % Night 
    
    disp([patient,' ',night]);
    
    %Loading the mat file to be converted 
    data = load([patient,'/',night,'/',patient,'_',night,'_200Hz_resampled.mat']);
    
    % EEG data file was saved in either variable name merged_matfile or file_new, hence using the if loop for the same 
    if isfield(data,'merged_matfile')
        EEGdata = data.merged_matfile';
    elseif isfield(data,'file_new')
        EEGdata = data.file_new';
    end
    
    % Transposing the data to be used for EEGLAB such that the data is in
    % the format --> Elec X Data in the .mat file
    save([patient,'/',night,'/',patient,'_',night,'_Transpose.mat'],'EEGdata');

%% Using EEGLAB functions 
    eeglab;
    EEG = pop_importdata([patient,'/',night,'/',patient,'_',night,'_Transpose.mat']);
    fs = 200;
    EEG.data = EEGdata;
    EEG = pop_editset(EEG, 'srate',fs); 
    % Filter parameters
    % Bandpass FIR filtering 
    
    locutoff = 0.1;  % Low CutOff Frequency
    hicutoff = 40;   % High CutOff Frequency
    filtorder = 4;  % Filter Order

    % Filtering the dataset
    EEG = pop_eegfilt( EEG, locutoff, hicutoff, filtorder);
    
    % Saving the data in .set file in the corresponding folder similar to
    % the dataset 
    EEG = pop_saveset(EEG, 'filename', [patient,'_',night], ...
    'filepath',['/Users/bsevak/Documents/Merged Data_BF/Merged_Data/',patient,'/',night,'/']);

    %csc_eeg_plotter(EEG);
    
    clearvars -except xfile; 
end
    
    
    