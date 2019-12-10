%% Convert .ns3 to .mat
% Oct 30, 2019

% .csv file with a list of all the folders to be checked for .ns3 to .mat
% file conversion
x = readtable('/Users/bsevak/Documents/FilesToConvert.xlsx');

len = 33;
while len <= length(x.Folder)
    disp(char(x.Folder(len)));
    cd(['/Volumes/EPILEPSY/BF/Data folders/',char(x.Folder(len))])
    %cd(['/Volumes/EPILEPSY/BF/Data folders/20150722-120659'])

    % Obtaining all the files with .ns3 extension
    fileList_ns3 = swa_getFiles(pwd, '.*.ns3');

    % Converting all the files to .mat files
    for i_ns3file = 1:size(fileList_ns3,1)
        disp(['File :', num2str(i_ns3file)]);
        data_resampled = [];
        filename_ns3 = fileList_ns3{i_ns3file};
        [data_info_ns3] = openNSx(filename_ns3, 'read');
        original_Fs = data_info_ns3.MetaTags.SamplingFreq;
        resampled_Fs = 200;
        headerinfo = data_info_ns3.MetaTags;
    
        data_ns3 = data_info_ns3.Data;
        header_ns3 = data_info_ns3.ElectrodesInfo;
        
        % Resampling the data to 200 Hz from 2000 Hz
        for ichannel = 1:size(data_info_ns3.Data,1)
            [data_ns3_chan] = data_ns3(ichannel,:)';
            data_channel_resampled = resample(double(data_ns3_chan),resampled_Fs,original_Fs);
            data_resampled = [data_resampled data_channel_resampled];
        end
        save([filename_ns3(:,51:end-4) '_resampled_200Hz.mat'], 'data_resampled', 'original_Fs', 'resampled_Fs', 'header_ns3','headerinfo', '-mat'),
    end
    
    len = len + 1;
end