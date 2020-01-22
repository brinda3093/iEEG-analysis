%% Analysis of different cycles of sleep and wake

%% Step 1

% Comparison between the first and last stage of the sleep cycle
patient = 'Patient_6';
night = 'Night_4';

pID = 'P6';
nID = 'N4';

cd(['/Users/bsevak/Documents/Merged Data_BF/Merged_Data/',patient,'/',night,'/']);

sleep = load([pID,nID,'_sleepData_fft_wholedata.mat']);
wakeBS = load([pID,nID,'_wakeBSData_fft_wholedata.mat']);
wakeAS = load([pID,nID,'_wakeASData_fft_wholedata.mat']);


%% Calculating the length of the sleep, wakeBS and wakeAS segment 
% 1 FFT window = 6s of data
% 2 hour of data = 7200 s of data = 1200 FFT windows
% 1 hour of data = 3600 s of data  = 600 FFT windows
% 1/2 hour of data = 1800 s of data = 300 FFT windows 

sleep_len = length(sleep.fft_SWA_2_4);
wakeBS_len = length(wakeBS.fft_SWA_2_4);
wakeAS_len = length(wakeAS.fft_SWA_2_4);
n_chan = length(sleep.fft_SWA_2_4(:,1));

%% Sleep FFT for comparison

if sleep_len > 1200
    s_epoch_len = 600;
    sleep_first_2_4 = sleep.fft_SWA_2_4(:,1:s_epoch_len);
    sleep_last_2_4 = sleep.fft_SWA_2_4(:,sleep_len-s_epoch_len+1:sleep_len);
    
    sleep_first_4_8 = sleep.fft_SWA_4_8(:,1:s_epoch_len);
    sleep_last_4_8 = sleep.fft_SWA_4_8(:,sleep_len-s_epoch_len+1:sleep_len);
else
    s_epoch_len = sleep_len/2;
    sleep_first_2_4 = sleep.fft_SWA_2_4(:,1:s_epoch_len);
    sleep_last_2_4 = sleep.fft_SWA_2_4(:,s_epoch_len+1:sleep_len);
    
    sleep_first_4_8 = sleep.fft_SWA_4_8(:,1:s_epoch_len);
    sleep_last_4_8 = sleep.fft_SWA_4_8(:,s_epoch_len+1:sleep_len);
end
    
%% Wake FFT for comparison

if (wakeBS_len >= 300 && wakeAS_len >= 300)
    w_epoch_len = 300;
    wakeBS_2_4 = wakeBS.fft_SWA_2_4(:,wakeBS_len-w_epoch_len+1:wakeBS_len);
    wakeAS_2_4 = wakeAS.fft_SWA_2_4(:,1:w_epoch_len);
    
    wakeBS_4_8 = wakeBS.fft_SWA_4_8(:,wakeBS_len-w_epoch_len+1:wakeBS_len);
    wakeAS_4_8 = wakeAS.fft_SWA_4_8(:,1:w_epoch_len);
else
    wake_len = min(wakeBS_len,wakeAS_len);
    w_epoch_len = wake_len;
    wakeBS_2_4 = wakeBS.fft_SWA_2_4(:,1:wake_len);
    wakeAS_2_4 = wakeAS.fft_SWA_2_4(:,1:wake_len);
    
    wakeBS_4_8 = wakeBS.fft_SWA_4_8(:,1:wake_len);
    wakeAS_4_8 = wakeAS.fft_SWA_4_8(:,1:wake_len);
end


%% Comparing the sleep and wake states in each channel

for i = 1:n_chan
    % Sleep - 2_4
    [h,p] = ttest2(sleep_first_2_4(i,:),sleep_last_2_4(i,:));
    barplot = figure;    
    barplot; 
    type = categorical({'First','Last'});
   
    bar(type, [mean(sleep_first_2_4(i,:)) mean(sleep_last_2_4(i,:))]); hold on
    
    title(['SWA Delta(2-4 Hz) FFT comparison in electrode ',num2str(i)]);
    xlabel('FFT power - Sleep');
    errorbar(type, [mean(sleep_first_2_4(i,:)) mean(sleep_last_2_4(i,:))] , [ std(sleep_first_2_4(i,:))  std(sleep_last_2_4(i,:))  ]/sqrt(s_epoch_len), '.','color','k') 
    
    saveas(barplot,['FirstvsLastSleep_2_4_E',num2str(i),'.png']);
    close(barplot);
    
    % Sleep - 4_8
    [h,p] = ttest2(sleep_first_4_8(i,:),sleep_last_4_8(i,:));
    barplot = figure;
    barplot; 
    type = categorical({'First','Last'});
   
    bar(type, [mean(sleep_first_4_8(i,:)) mean(sleep_last_4_8(i,:))]); hold on
    
    title(['SWA Theta(4-8 Hz) FFT comparison in electrode ',num2str(i)]);
    xlabel('FFT power - Sleep');
    errorbar(type, [mean(sleep_first_4_8(i,:)) mean(sleep_last_4_8(i,:))] , [ std(sleep_first_4_8(i,:))  std(sleep_last_4_8(i,:))  ]/sqrt(s_epoch_len), '.','color','k') 
    
    saveas(barplot,['FirstvsLastSleep_4_8_E',num2str(i),'.png']);
    close(barplot);
    
    
%     % Wake - 2_4
    [h,p] = ttest2(wakeBS_2_4(i,:),wakeAS_2_4(i,:));
    barplot = figure;
    barplot; 
    type = categorical({'Before','After'});
   
    bar(type, [mean(wakeBS_2_4(i,:)) mean(wakeAS_2_4(i,:))]); hold on
    title(['SWA Delta(2-4 Hz) FFT comparison in electrode ',num2str(i)]);
    xlabel('FFT power - Wake');
    errorbar(type, [mean(wakeBS_2_4(i,:)) mean(wakeAS_2_4(i,:))] , [ std(wakeBS_2_4(i,:))  std(wakeAS_2_4(i,:))  ]/sqrt(s_epoch_len), '.','color','k') 
    
    saveas(barplot,['BeforevsAfterWake_2_4_E',num2str(i),'.png']);
    close(barplot);
    
     % Wake - 4_8
    [h,p] = ttest2(wakeBS_4_8(i,:),wakeAS_4_8(i,:));
    barplot = figure;
    barplot; 
    type = categorical({'Before','After'});
   
    bar(type, [mean(wakeBS_4_8(i,:)) mean(wakeAS_4_8(i,:))]); hold on
    title(['SWA Theta(4-8 Hz) FFT comparison in electrode ',num2str(i)]);
    xlabel('FFT power - Wake');
    errorbar(type, [mean(wakeBS_4_8(i,:)) mean(wakeAS_4_8(i,:))] , [ std(wakeBS_4_8(i,:))  std(wakeAS_4_8(i,:))  ]/sqrt(s_epoch_len), '.','color','k') 
   
    saveas(barplot,['BeforevsAfterWake_4_8_E',num2str(i),'.png']);
    close(barplot);
end




