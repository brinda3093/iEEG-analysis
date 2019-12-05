
x = readtable('/Users/bsevak/Documents/TrialAutomate.xlsx');

for filedet = 25:length(x.File)
    pathP = (['/Volumes/EPILEPSY/BF/Data folders/',char(x.File(filedet)),'/']);
cd(pathP);
merged_matfile = [];
fileList_resampled = swa_getFiles(pwd,'_resampled_200Hz.mat');
data_resampled_merged = [];
for i_resampledfile = 1:size(fileList_resampled,1)
   filename_resampled = fileList_resampled{i_resampledfile};
   load(fileList_resampled{i_resampledfile});
   data_resampled_merged = [data_resampled_merged; data_resampled]; 
end

merged_matfile = [merged_matfile;   data_resampled_merged];
%merged_matfile = [merged_matfile; zeros(1,size(merged_matfile,2))];

%save([filename_resampled(:,1:end-25) '_200Hz_resampled_merged.mat'],...
%    'merged_matfile','header_ns5','headerinfo', '-mat', '-v7.5');

patient = [char(x.Patient(filedet))]; %'Patient_6';
night =  [char(x.Night(filedet))];%'Night_2_4';
                                                                                                  
cd(['/Users/bsevak/Documents/Merged Data_BF/Merged_Data/',patient,'/',night]);

save([patient,'_',night,'_200Hz_resampled.mat'],...
    'merged_matfile','header_ns3','headerinfo', '-mat', '-v7.3');
disp('Files Merged');
end