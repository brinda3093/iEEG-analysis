%% Analysis of different cycles of sleep and wake

%% Step 1

% Comparison between the first and last stage of the sleep cycle
%patient = 'Patient_6';
%night = 'Night_4';

%pID = 'P6';
%nID = 'N4';

%cd(['/Users/bsevak/Documents/Merged Data_BF/Merged_Data/',patient,'/',night,'/']);

%sleep =load('_fft_wholedata.mat');
sleep = load(['sleepData_fft_wholedata.mat']);
wakeBS = load(['wakeBSData_fft_wholedata.mat']);
wakeAS = load(['wakeASData_fft_wholedata.mat']);
load('CUBF10_0208_channel_labels.mat');
load('CUBF10_208_hiddenchans.mat');

chan_lab = [];
for i = 1:length(channel_labels)
    if ~ismember(i,hidden_chans)
        chan_lab = [chan_lab, channel_labels(i).labels];
    end
end

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
    dataS2_4.(['E',num2str(i)]) = [sleep_first_2_4(i,:)', sleep_last_2_4(i,:)'];
    bar(type, [mean(sleep_first_2_4(i,:)) mean(sleep_last_2_4(i,:))]); hold on
    
    title(['SWA Delta(2-4 Hz) FFT comparison in electrode ',chan_lab(i)]);
    xlabel('FFT power - Sleep');
    errorbar(type, [mean(sleep_first_2_4(i,:)) mean(sleep_last_2_4(i,:))] , [ std(sleep_first_2_4(i,:))  std(sleep_last_2_4(i,:))  ]/sqrt(s_epoch_len), '.','color','k') 
    
    saveas(barplot,['FirstvsLastSleep_2_4_E',char(chan_lab(i)),'.png']);
    close(barplot);
    
    % Sleep - 4_8
    [h,p] = ttest2(sleep_first_4_8(i,:),sleep_last_4_8(i,:));
    barplot = figure;
    barplot; 
    type = categorical({'First','Last'});
    dataS4_8.(['E',num2str(i)]) = [sleep_first_4_8(i,:)', sleep_last_4_8(i,:)'];
    
    bar(type, [mean(sleep_first_4_8(i,:)) mean(sleep_last_4_8(i,:))]); hold on
    
    title(['SWA Theta(4-8 Hz) FFT comparison in electrode ',chan_lab(i)]);
    xlabel('FFT power - Sleep');
    errorbar(type, [mean(sleep_first_4_8(i,:)) mean(sleep_last_4_8(i,:))] , [ std(sleep_first_4_8(i,:))  std(sleep_last_4_8(i,:))  ]/sqrt(s_epoch_len), '.','color','k') 
    
    saveas(barplot,['FirstvsLastSleep_4_8_E',char(chan_lab(i)),'.png']);
    close(barplot);
    
    
%     % Wake - 2_4
    [h,p] = ttest2(wakeBS_2_4(i,:),wakeAS_2_4(i,:));
    barplot = figure;
    barplot; 
    type = categorical({'Before','After'});
    dataW2_4.(['E',num2str(i)]) = [wakeBS_2_4(i,:)', wakeAS_2_4(i,:)'];
    
    bar(type, [mean(wakeBS_2_4(i,:)) mean(wakeAS_2_4(i,:))]); hold on
    title(['SWA Delta(2-4 Hz) FFT comparison in electrode ',chan_lab(i)]);
    xlabel('FFT power - Wake');
    errorbar(type, [mean(wakeBS_2_4(i,:)) mean(wakeAS_2_4(i,:))] , [ std(wakeBS_2_4(i,:))  std(wakeAS_2_4(i,:))  ]/sqrt(s_epoch_len), '.','color','k') 
    
    saveas(barplot,['BeforevsAfterWake_2_4_E',char(chan_lab(i)),'.png']);
    close(barplot);
    
     % Wake - 4_8
    [h,p] = ttest2(wakeBS_4_8(i,:),wakeAS_4_8(i,:));
    barplot = figure;
    barplot; 
    type = categorical({'Before','After'});
    dataW4_8.(['E',num2str(i)]) = [wakeBS_4_8(i,:)', wakeAS_2_4(i,:)'];
    
    bar(type, [mean(wakeBS_4_8(i,:)) mean(wakeAS_4_8(i,:))]); hold on
    title(['SWA Theta(4-8 Hz) FFT comparison in electrode ',chan_lab(i)]);
    xlabel('FFT power - Wake');
    errorbar(type, [mean(wakeBS_4_8(i,:)) mean(wakeAS_4_8(i,:))] , [ std(wakeBS_4_8(i,:))  std(wakeAS_4_8(i,:))  ]/sqrt(s_epoch_len), '.','color','k') 
   
    saveas(barplot,['BeforevsAfterWake_4_8_E',char(chan_lab(i)),'.png']);
    close(barplot);
end


%% Different type of graph
LAMY_B = [];
for elecs = 1:7
    LAMY_B = [wakeBS_2_4(elecs,:)',LAMY_B];
end


LAMY_A = [];
for elecs = 1:7
    LAMY_A = [wakeAS_2_4(elecs,:)',LAMY_A];
end

%% 
dataSBefore2_4 = [];

for i = 1:n_chan
    dataSBefore2_4 = [dataS2_4.(['E',num2str(i)])(:,1), dataSBefore2_4];
end

dataSAfter2_4 = [];

for i = 1:n_chan
    dataSAfter2_4 = [dataS2_4.(['E',num2str(i)])(:,2), dataSAfter2_4];
end

xylabel = repmat('ba',47,1);
eleclabel = repmat(1:47,1,2);
h = figure;
boxplot([dataSBefore2_4 dataSAfter2_4], {eleclabel' xylabel(:)},'factorgap',1);
title('Before vs After Sleep 2-4 Hz');
%%

dataSBefore4_8 = [];

for i = 1:n_chan
    dataSBefore4_8 = [dataS4_8.(['E',num2str(i)])(:,1), dataSBefore4_8];
end

dataSAfter4_8 = [];

for i = 1:n_chan
    dataSAfter4_8 = [dataS4_8.(['E',num2str(i)])(:,2), dataSAfter4_8];
end

xylabel = repmat('ba',47,1);
eleclabel = repmat(1:47,1,2);
h = figure;
boxplot([dataSBefore4_8 dataSAfter4_8], {eleclabel' xylabel(:)},'factorgap',1);
title('Before vs After Sleep 4-8 Hz');

%%
dataWBefore4_8 = [];

for i = 1:n_chan
    dataWBefore4_8 = [dataW4_8.(['E',num2str(i)])(:,1), dataWBefore4_8];
end

dataWAfter4_8 = [];

for i = 1:n_chan
    dataWAfter4_8 = [dataW4_8.(['E',num2str(i)])(:,2), dataWAfter4_8];
end

xylabel = repmat('ba',47,1);
eleclabel = repmat(1:47,1,2);
h = figure;
boxplot([dataWBefore4_8 dataWAfter4_8], {eleclabel' xylabel(:)},'factorgap',1);
title('Before vs After Wake 4-8 Hz');
%%

dataWBefore2_4 = [];

for i = 1:n_chan
    dataWBefore2_4 = [dataW2_4.(['E',num2str(i)])(:,1), dataWBefore2_4];
end

dataWAfter2_4 = [];

for i = 1:n_chan
    dataWAfter2_4 = [dataW2_4.(['E',num2str(i)])(:,2), dataWAfter2_4];
end

xylabel = repmat('ba',47,1);
eleclabel = repmat(1:47,1,2);
h = figure;
boxplot([dataWBefore2_4 dataWAfter2_4], {eleclabel' xylabel(:)},'factorgap',1);
title('Before vs After Wake 2-4 Hz');