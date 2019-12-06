%% .mat to .set conversion for using csc_eeg_plotter()

cd('/Users/bsevak/Documents/Merged Data_BF/Merged_Data/');
xfile = readtable('/Users/bsevak/Documents/convertToSet.xlsx');
for i = 37:length(xfile.Patient) 
    patient = char(xfile.Patient(i));
    night = char(xfile.Night(i));
    disp([patient,' ',night]);
    data = load([patient,'/',night,'/',patient,'_',night,'_200Hz_resampled.mat']);
    if isfield(data,'merged_matfile')
        EEGdata = data.merged_matfile';
    elseif isfield(data,'file_new')
        EEGdata = data.file_new';
    end
    save([patient,'/',night,'/',patient,'_',night,'_Transpose.mat'],'EEGdata');

%%
    eeglab;
    EEG = pop_importdata([patient,'/',night,'/',patient,'_',night,'_Transpose.mat']);
% Filter parameters
    locutoff = 0.1;
    hicutoff = 40;
    filtorder = 4;

    EEGOUT = pop_eegfilt( EEG, locutoff, hicutoff, filtorder);
    EEG = pop_saveset(EEGOUT, 'filename', [patient,'_',night], ...
    'filepath',['/Users/bsevak/Documents/Merged Data_BF/Merged_Data/',patient,'/',night,'/']);
    clearvars -except xfile; 
end