%% Merging the .mat files
% Nov 12, 2019

%% Merge multiple files with the same names
% clearvars;
% clc;
% 
% 
% %fileName = num2str(filel(c,5))
%     
% cd('/Volumes/EPILEPSY/BF/Data folders/20150715-152219');
% merged_matfile = [];
% fileList_resampled = swa_getFiles(pwd,'_resampled_200Hz.mat');
% data_resampled_merged = [];
% for i_resampledfile = 1:size(fileList_resampled,1)
%    filename_resampled = fileList_resampled{i_resampledfile};
%    load(fileList_resampled{i_resampledfile});
%    data_resampled_merged = [data_resampled_merged; data_resampled]; 
% end
% 
% merged_matfile = [merged_matfile;   data_resampled_merged];
% %merged_matfile = [merged_matfile; zeros(1,size(merged_matfile,2))];
% 
% %save([filename_resampled(:,1:end-25) '_200Hz_resampled_merged.mat'],...
% %    'merged_matfile','header_ns5','headerinfo', '-mat', '-v7.5');
% 
% patient = 'Patient_3';
% night = 'Night_1_2';
% 
% cd(['/Users/bsevak/Documents/Merged Data_BF/Merged_Data/',patient,'/',night]);
% 
% save([patient,'_',night,'_200Hz_resampled.mat'],...
%     'merged_matfile','header_ns3','headerinfo', '-mat', '-v7.3');
% disp('Files Merged');
% % % 


%% Merge files with different file names and same night
% 
clear;
clc;
patient = 'Patient_3';
night = 'Night_1';

cd(['/Users/bsevak/Documents/Merged Data_BF/Merged_Data/',patient,'/',night,'/']);

% Loading the .mat files for the same day different names
file1 = load(['/Users/bsevak/Documents/Merged Data_BF/Merged_Data/',patient,'/',night,'_1/',patient,'_',night,'_1_200Hz_resampled.mat']);
file2 = load(['/Users/bsevak/Documents/Merged Data_BF/Merged_Data/',patient,'/',night,'_2/',patient,'_',night,'_2_200Hz_resampled.mat']);
%file3 = load(['/Users/bsevak/Documents/Merged Data_BF/Merged_Data/',patient,'/',night,'_3/',patient,'_',night,'_3_200Hz_resampled.mat']);
%file4 = load(['/Users/bsevak/Documents/Merged Data_BF/Merged_Data/',patient,'/',night,'_4/',patient,'_',night,'_4_200Hz_resampled.mat']);
%file5 = load(['/Users/bsevak/Documents/Merged Data_BF/Merged_Data/',patient,'/',night,'_5/',patient,'_',night,'_5_200Hz_resampled.mat']);
% %file6 = load(['/Users/bsevak/Documents/Merged Data_BF/Merged_Data/',patient,'/',night,'_6/',patient,'_',night,'_6_200Hz_resampled.mat']);
% file7 = load(['/Users/bsevak/Documents/Merged Data_BF/Merged_Data/',patient,'/',night,'_7/',patient,'_',night,'_7_200Hz_resampled.mat']);
% file8 = load(['/Users/bsevak/Documents/Merged Data_BF/Merged_Data/',patient,'/',night,'_8/',patient,'_',night,'_8_200Hz_resampled.mat']);

merged_matfile = [file1.merged_matfile ; file2.merged_matfile];%file3.merged_matfile ];%file4.merged_matfile];%...
    %file5.merged_matfile ];% file7.merged_matfile ;file8.merged_matfile];
header_ns3 = file1.header_ns3;


save([patient,'_',night,'_200Hz_resampled.mat'],'merged_matfile','header_ns3','-mat','-v7.3');
disp(['Files Merged for ',night]);