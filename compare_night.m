%% Analysis of different cycles of sleep and wake

%% Step 1

% Comparison between the first and last stage of the sleep cycle
%patient = 'Patient_6';
%night = 'Night_4';

pID = 'CUBF09';
nID = '0123';

%cd(['/Users/bsevak/Documents/Merged Data_BF/Merged_Data/',patient,'/',night,'/']);

%
sleep =load('sleepData_full_fft_wholedata.mat');
spike = 'WOSpikes';
%sleep_first = load(['sleepData_first_fft_wholedata.mat']);
%sleep_last = load(['sleepData_last_fft_wholedata.mat']);

wakeBS = load(['wakeBSData_fft_wholedata.mat']);
wakeAS = load(['wakeASData_fft_wholedata.mat']);
load([pID,'_',nID,'_channel_labels.mat']);
load([pID,'_',nID,'_hiddenchans.mat']);

chan_lab = [];
for i = 1:length(channel_labels)
   if ~ismember(i,hidden_chans)
        chan_lab = [chan_lab, channel_labels(i)];%.labels];
   end
end

%% Calculating the length of the sleep, wakeBS and wakeAS segment 
% 1 FFT window = 6s of data
% 2 hour of data = 7200 s of data = 1200 FFT windows
% 1 hour of data = 3600 s of data  = 600 FFT windows
% 1/2 hour of data = 1800 s of data = 300 FFT windows 

sleep_len = length(sleep.fft_SWA_2_4(1,:));
%sleepFirst_len = length(sleep_first.fft_SWA_2_4(1,:));
%sleepLast_len = length(sleep_last.fft_SWA_2_4(1,:));
wakeBS_len = length(wakeBS.fft_SWA_2_4(1,:));
wakeAS_len = length(wakeAS.fft_SWA_2_4(1,:));
n_chan = length(sleep.fft_SWA_2_4(:,1));

%% Complete Sleep Data

sleep_full_2_4 = sleep.fft_SWA_2_4;
sleep_full_4_8 = sleep.fft_SWA_4_8;
sleep_full_Delta = sleep.fft_bands(:,:,3);
sleep_full_Beta = sleep.fft_bands(:,:,6);
sleep_full_DB = sleep.fft_SWA_DB;

for i = 1:n_chan
    dataSDelta_full.(['E',num2str(i)]) = [sleep_full_Delta(i,:)'];
end

for i = 1:n_chan
    DeltaS_full(i) = mean(dataSDelta_full.(['E',num2str(i)])(:,1)) ;
end



%% Sleep FFT for comparison

if sleep_len > 1200
    s_epoch_len = 600;
    sleep_first_2_4 = sleep.fft_SWA_2_4(:,1:s_epoch_len);
    sleep_last_2_4 = sleep.fft_SWA_2_4(:,sleep_len-s_epoch_len+1:sleep_len);
    
    sleep_first_4_8 = sleep.fft_SWA_4_8(:,1:s_epoch_len);
    sleep_last_4_8 = sleep.fft_SWA_4_8(:,sleep_len-s_epoch_len+1:sleep_len);
    
    sleep_first_Delta = sleep.fft_bands(:,1:s_epoch_len,3);
    sleep_last_Delta = sleep.fft_bands(:,sleep_len-s_epoch_len+1:sleep_len,3);
    
    sleep_first_Beta = sleep.fft_bands(:,1:s_epoch_len,6);
    sleep_last_Beta = sleep.fft_bands(:,sleep_len-s_epoch_len+1:sleep_len,6);
    
    sleep_first_DB = sleep.fft_SWA_DB(:,1:s_epoch_len);
    sleep_last_DB = sleep.fft_SWA_DB(:,sleep_len-s_epoch_len+1:sleep_len);
else
    s_epoch_len = sleep_len/2;
    sleep_first_2_4 = sleep.fft_SWA_2_4(:,1:s_epoch_len);
    sleep_last_2_4 = sleep.fft_SWA_2_4(:,s_epoch_len+1:sleep_len);
    
    sleep_first_4_8 = sleep.fft_SWA_4_8(:,1:s_epoch_len);
    sleep_last_4_8 = sleep.fft_SWA_4_8(:,s_epoch_len+1:sleep_len);
    
    sleep_first_Delta = sleep.fft_bands(:,1:s_epoch_len,3);
    sleep_last_Delta = sleep.fft_bands(:,s_epoch_len+1:sleep_len,3);
    
    sleep_first_Beta = sleep.fft_bands(:,1:s_epoch_len,6);
    sleep_last_Beta = sleep.fft_bands(:,s_epoch_len+1:sleep_len,6);
    
    sleep_first_DB = sleep.fft_SWA_DB(:,1:s_epoch_len);
    sleep_last_DB = sleep.fft_SWA_DB(:,s_epoch_len+1:sleep_len);
end

%% Sleep & Wake First and Last - FFT
if sleepFirst_len > sleepLast_len
    s_epoch_len = sleepLast_len;
else
    s_epoch_len = sleepFirst_len;
end

sleep_first_2_4 = sleep_first.fft_SWA_2_4(:,1:s_epoch_len);
sleep_last_2_4 = sleep_last.fft_SWA_2_4(:,1:s_epoch_len);
    
sleep_first_4_8 = sleep_first.fft_SWA_4_8(:,1:s_epoch_len);
sleep_last_4_8 = sleep_last.fft_SWA_4_8(:,1:s_epoch_len);

sleep_first_Delta = sleep_first.fft_bands(:,1:s_epoch_len,3);
sleep_last_Delta = sleep_last.fft_bands(:,1:s_epoch_len,3);

sleep_first_Beta = sleep_first.fft_bands(:,1:s_epoch_len,6);
sleep_last_Beta = sleep_last.fft_bands(:,1:s_epoch_len,6);
    
sleep_first_DB = sleep_first.fft_SWA_DB(:,1:s_epoch_len);
sleep_last_DB = sleep_last.fft_SWA_DB(:,1:s_epoch_len);

%% Wake FFT for comparison

if (wakeBS_len >= 300 && wakeAS_len >= 300)
    w_epoch_len = 300;
    wakeBS_2_4 = wakeBS.fft_SWA_2_4(:,wakeBS_len-w_epoch_len+1:wakeBS_len);
    wakeAS_2_4 = wakeAS.fft_SWA_2_4(:,1:w_epoch_len);
    
    wakeBS_4_8 = wakeBS.fft_SWA_4_8(:,wakeBS_len-w_epoch_len+1:wakeBS_len);
    wakeAS_4_8 = wakeAS.fft_SWA_4_8(:,1:w_epoch_len);
    
    wakeBS_Delta = wakeBS.fft_bands(:,wakeBS_len-w_epoch_len+1:wakeBS_len,3);
    wakeAS_Delta = wakeAS.fft_bands(:,1:w_epoch_len,3);
    
    wakeBS_Beta = wakeBS.fft_bands(:,wakeBS_len-w_epoch_len+1:wakeBS_len,6);
    wakeAS_Beta = wakeAS.fft_bands(:,1:w_epoch_len,6);
    
    wakeBS_DB = wakeBS.fft_SWA_DB(:,wakeBS_len-w_epoch_len+1:wakeBS_len);
    wakeAS_DB = wakeAS.fft_SWA_DB(:,1:w_epoch_len);
else
    wake_len = min(wakeBS_len,wakeAS_len);
    w_epoch_len = wake_len;
    wakeBS_2_4 = wakeBS.fft_SWA_2_4(:,1:wake_len);
    wakeAS_2_4 = wakeAS.fft_SWA_2_4(:,1:wake_len);
    
    wakeBS_4_8 = wakeBS.fft_SWA_4_8(:,1:wake_len);
    wakeAS_4_8 = wakeAS.fft_SWA_4_8(:,1:wake_len);
    
    wakeBS_Delta = wakeBS.fft_bands(:,1:wake_len,3);
    wakeAS_Delta = wakeAS.fft_bands(:,1:wake_len,3);
    
    wakeBS_Beta = wakeBS.fft_bands(:,1:wake_len,6);
    wakeAS_Beta = wakeAS.fft_bands(:,1:wake_len,6);
    
    wakeBS_DB = wakeBS.fft_SWA_DB(:,1:wake_len);
    wakeAS_DB = wakeAS.fft_SWA_DB(:,1:wake_len);
end


%% Comparing the sleep and wake states in each channel


for i = 1:n_chan
    dataSDelta.(['E',num2str(i)]) = [sleep_first_Delta(i,:)', sleep_last_Delta(i,:)'];
    dataWDelta.(['E',num2str(i)]) = [wakeBS_Delta(i,:)', wakeAS_Delta(i,:)'];
end

for i = 1:n_chan
    % Sleep - 2_4
    %[h,p] = ttest2(sleep_first_2_4(i,:),sleep_last_2_4(i,:));
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
    %[h,p] = ttest2(sleep_first_4_8(i,:),sleep_last_4_8(i,:));
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
    
    % Sleep - DB
    %[h,p] = ttest2(sleep_first_DB(i,:),sleep_last_DB(i,:));
    barplot = figure;
    barplot; 
    type = categorical({'First','Last'});
    dataSBeta.(['E',num2str(i)]) = [sleep_first_Beta(i,:)', sleep_last_Beta(i,:)'];
    
    bar(type, [mean(sleep_first_Beta(i,:)) mean(sleep_last_Beta(i,:))]); hold on
    
    title(['SWA Beta FFT comparison in electrode ',chan_lab(i)]);
    xlabel('FFT power - Sleep');
    errorbar(type, [mean(sleep_first_Beta(i,:)) mean(sleep_last_Beta(i,:))] , [ std(sleep_first_Beta(i,:))  std(sleep_last_Beta(i,:))  ]/sqrt(s_epoch_len), '.','color','k') 
    
    saveas(barplot,['FirstvsLastSleep_Beta_E',char(chan_lab(i)),'.png']);
    close(barplot);
    
    % Sleep - DB
    %[h,p] = ttest2(sleep_first_DB(i,:),sleep_last_DB(i,:));
    barplot = figure;
    barplot; 
    type = categorical({'First','Last'});
    dataSDB.(['E',num2str(i)]) = [sleep_first_DB(i,:)', sleep_last_DB(i,:)'];
    
    bar(type, [mean(sleep_first_DB(i,:)) mean(sleep_last_DB(i,:))]); hold on
    
    title(['SWA Delta Beta FFT comparison in electrode ',chan_lab(i)]);
    xlabel('FFT power - Sleep');
    errorbar(type, [mean(sleep_first_DB(i,:)) mean(sleep_last_DB(i,:))] , [ std(sleep_first_DB(i,:))  std(sleep_last_DB(i,:))  ]/sqrt(s_epoch_len), '.','color','k') 
    
    saveas(barplot,['FirstvsLastSleep_DB_E',char(chan_lab(i)),'.png']);
    close(barplot);
    
    
% %     % Wake - 2_4
%    [h,p] = ttest2(wakeBS_2_4(i,:),wakeAS_2_4(i,:));
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
%    [h,p] = ttest2(wakeBS_4_8(i,:),wakeAS_4_8(i,:));
    barplot = figure;
    barplot; 
    type = categorical({'Before','After'});
    dataW4_8.(['E',num2str(i)]) = [wakeBS_4_8(i,:)', wakeAS_4_8(i,:)'];
    
    bar(type, [mean(wakeBS_4_8(i,:)) mean(wakeAS_4_8(i,:))]); hold on
    title(['SWA Theta(4-8 Hz) FFT comparison in electrode ',chan_lab(i)]);
    xlabel('FFT power - Wake');
    errorbar(type, [mean(wakeBS_4_8(i,:)) mean(wakeAS_4_8(i,:))] , [ std(wakeBS_4_8(i,:))  std(wakeAS_4_8(i,:))  ]/sqrt(s_epoch_len), '.','color','k') 
   
    saveas(barplot,['BeforevsAfterWake_4_8_E',char(chan_lab(i)),'.png']);
    close(barplot);
    
    % Wake - Beta ratio
%    [h,p] = ttest2(wakeBS_DB(i,:),wakeAS_DB(i,:));
    barplot = figure;
    barplot; 
    type = categorical({'Before','After'});
    dataWBeta.(['E',num2str(i)]) = [wakeBS_Beta(i,:)', wakeAS_Beta(i,:)'];
    
    bar(type, [mean(wakeBS_Beta(i,:)) mean(wakeAS_Beta(i,:))]); hold on
    title(['SWA Beta Ratio FFT comparison in electrode ',chan_lab(i)]);
    xlabel('FFT power - Wake');
    errorbar(type, [mean(wakeBS_Beta(i,:)) mean(wakeAS_Beta(i,:))] , [ std(wakeBS_Beta(i,:))  std(wakeAS_Beta(i,:))  ]/sqrt(s_epoch_len), '.','color','k') 
   
    saveas(barplot,['BeforevsAfterWake_Beta_E',char(chan_lab(i)),'.png']);
    close(barplot);  
    
    % Wake - DB ratio
%    [h,p] = ttest2(wakeBS_DB(i,:),wakeAS_DB(i,:));
    barplot = figure;
    barplot; 
    type = categorical({'Before','After'});
    dataWDB.(['E',num2str(i)]) = [wakeBS_DB(i,:)', wakeAS_DB(i,:)'];
    
    bar(type, [mean(wakeBS_DB(i,:)) mean(wakeAS_DB(i,:))]); hold on
    title(['SWA Delta Beta Ratio FFT comparison in electrode ',chan_lab(i)]);
    xlabel('FFT power - Wake');
    errorbar(type, [mean(wakeBS_DB(i,:)) mean(wakeAS_DB(i,:))] , [ std(wakeBS_DB(i,:))  std(wakeAS_DB(i,:))  ]/sqrt(s_epoch_len), '.','color','k') 
   
    saveas(barplot,['BeforevsAfterWake_DB_E',char(chan_lab(i)),'.png']);
    close(barplot);  
    
    
end

%% Comparing the different electrodes for sleep and wake

for i = 1:n_chan
   
   compareS2_4(i) = (mean(dataS2_4.(['E',num2str(i)])(:,2)) - mean(dataS2_4.(['E',num2str(i)])(:,1)))/...
       (mean(dataS2_4.(['E',num2str(i)])(:,1))); %+ mean(dataS2_4.(['E',num2str(i)])(:,2)));
   
   compareS4_8(i) = (mean(dataS4_8.(['E',num2str(i)])(:,2)) - mean(dataS4_8.(['E',num2str(i)])(:,1)))/...
       (mean(dataS4_8.(['E',num2str(i)])(:,1))); %+ mean(dataS4_8.(['E',num2str(i)])(:,2)));
   
   compareSDelta(i) = (mean(dataSDelta.(['E',num2str(i)])(:,2)) - mean(dataSDelta.(['E',num2str(i)])(:,1)))/...
       (mean(dataSDelta.(['E',num2str(i)])(:,1))); %+ mean(dataSDelta.(['E',num2str(i)])(:,2)));
   
   compareSBeta(i) = (mean(dataSBeta.(['E',num2str(i)])(:,2)) - mean(dataSBeta.(['E',num2str(i)])(:,1)))/...
       (mean(dataSBeta.(['E',num2str(i)])(:,1))); %+ mean(dataSBeta.(['E',num2str(i)])(:,2)));
   
   compareSDB(i) = (mean(dataSDB.(['E',num2str(i)])(:,2)) - mean(dataSDB.(['E',num2str(i)])(:,1)))/...
       (mean(dataSDB.(['E',num2str(i)])(:,1))); %+ mean(dataSDB.(['E',num2str(i)])(:,2)));
   
   compareW2_4(i) = (mean(dataW2_4.(['E',num2str(i)])(:,2)) - mean(dataW2_4.(['E',num2str(i)])(:,1)))/...
       (mean(dataW2_4.(['E',num2str(i)])(:,1)));% + mean(dataW2_4.(['E',num2str(i)])(:,2)));
   
   compareW4_8(i) = (mean(dataW4_8.(['E',num2str(i)])(:,2)) - mean(dataW4_8.(['E',num2str(i)])(:,1)))/...
       (mean(dataW4_8.(['E',num2str(i)])(:,1)));% + mean(dataW4_8.(['E',num2str(i)])(:,2)));
   
    compareWDelta(i) = (mean(dataWDelta.(['E',num2str(i)])(:,2)) - mean(dataWDelta.(['E',num2str(i)])(:,1)))/...
       (mean(dataWDelta.(['E',num2str(i)])(:,1)));% + mean(dataWDelta.(['E',num2str(i)])(:,2)));
   
   compareWBeta(i) = (mean(dataWBeta.(['E',num2str(i)])(:,2)) - mean(dataWBeta.(['E',num2str(i)])(:,1)))/...
       (mean(dataWBeta.(['E',num2str(i)])(:,1)));% + mean(dataWBeta.(['E',num2str(i)])(:,2)));
   
   compareWDB(i) = (mean(dataWDB.(['E',num2str(i)])(:,2)) - mean(dataWDB.(['E',num2str(i)])(:,1)))/...
       (mean(dataWDB.(['E',num2str(i)])(:,1)));% + mean(dataWDB.(['E',num2str(i)])(:,2)));
   
   % Absolute Values
   compareS2_4_abs(i) = (mean(dataS2_4.(['E',num2str(i)])(:,2)) - mean(dataS2_4.(['E',num2str(i)])(:,1)));
      
   compareS4_8_abs(i) = (mean(dataS4_8.(['E',num2str(i)])(:,2)) - mean(dataS4_8.(['E',num2str(i)])(:,1)));
   
   compareSDelta_abs(i) = (mean(dataSDelta.(['E',num2str(i)])(:,2)) - mean(dataSDelta.(['E',num2str(i)])(:,1)));
       
   compareSBeta_abs(i) = (mean(dataSBeta.(['E',num2str(i)])(:,2)) - mean(dataSBeta.(['E',num2str(i)])(:,1)));
   
   compareSDB_abs(i) = (mean(dataSDB.(['E',num2str(i)])(:,2)) - mean(dataSDB.(['E',num2str(i)])(:,1)));
      
   compareW2_4_abs(i) = (mean(dataW2_4.(['E',num2str(i)])(:,2)) - mean(dataW2_4.(['E',num2str(i)])(:,1)));
       
   compareW4_8_abs(i) = (mean(dataW4_8.(['E',num2str(i)])(:,2)) - mean(dataW4_8.(['E',num2str(i)])(:,1)));
   
   compareWDelta_abs(i) = (mean(dataWDelta.(['E',num2str(i)])(:,2)) - mean(dataWDelta.(['E',num2str(i)])(:,1)));
   
   compareWBeta_abs(i) = (mean(dataWBeta.(['E',num2str(i)])(:,2)) - mean(dataWBeta.(['E',num2str(i)])(:,1)));
   
   compareWDB_abs(i) = (mean(dataWDB.(['E',num2str(i)])(:,2)) - mean(dataWDB.(['E',num2str(i)])(:,1)));  
   
end

save([pID,'_',nID,'_forComparison_',spike,'.mat'],'compareS2_4','compareS4_8','compareSBeta','compareSDB','compareW2_4','compareW4_8','compareWBeta','compareWDB',...
    'compareS2_4_abs','compareS4_8_abs','compareSBeta_abs','compareSDB_abs','compareW2_4_abs','compareW4_8_abs','compareWBeta_abs','compareWDB_abs',...
    'dataS2_4','dataS4_8','dataSBeta','dataSDB','dataW2_4','dataW4_8','dataWBeta','dataWDB','chan_lab','hidden_chans','dataSDelta','dataWDelta',...
    'compareSDelta','compareWDelta','compareSDelta_abs','compareWDelta_abs');

%% S_2_4

figure;
sgtitle(['Sleep Percentage',spike]);

subplot(1,4,1)
imagesc(compareS2_4');
caxis([-1 1]);
yticks(1:n_chan)
yticklabels(chan_lab);
title(['2-4 Hz']);
colormap(jet);
colorbar;

% figure;
% title('Sleep 2-4 Hz');
% hintonDiagram(compareS2_4');
% yticks([1:n_chan]);
% yticklabels(chan_lab);

%figure;
subplot(1,4,2);
imagesc(compareS4_8');
caxis([-1 1]);
yticks(1:n_chan);
yticklabels(chan_lab);
title(['4-8 Hz']);
colorbar;

subplot(1,4,3);
imagesc(compareSBeta');
caxis([-1 1]);
yticks(1:n_chan);
yticklabels(chan_lab);
title(['Beta']);
colorbar;
   
%figure;
subplot(1,4,4);
imagesc(compareSDB');
caxis([-1 1]);
yticks(1:n_chan);
yticklabels(chan_lab);
title(['Delta-Beta']);
colorbar;


% 


%% Wake

figure;
sgtitle(['Wake Percentage',spike]);

subplot(1,4,1);
imagesc(compareW2_4');
caxis([-0.5 0.5]);
yticks(1:n_chan);
yticklabels(chan_lab);
title('2-4 Hz');
colormap(jet);
colorbar;


subplot(1,4,2);
imagesc(compareW4_8');
caxis([-0.5 0.5]);
yticks(1:n_chan);
yticklabels(chan_lab);
title(['4-8 Hz']);
colorbar;

subplot(1,4,3);
imagesc(compareWBeta');
caxis([-0.5 0.5]);
yticks(1:n_chan);
yticklabels(chan_lab);
title(['Beta']);
colorbar;

subplot(1,4,4);
imagesc(compareWDB');
caxis([-0.5 0.5]);
yticks(1:n_chan);
yticklabels(chan_lab);
title(['Delta-Beta']);
colorbar;


%% Absolute Values

%Sleep
figure;
sgtitle(['Sleep Absolute',spike]);

subplot(1,4,1)
imagesc(compareS2_4_abs');
caxis([-500 500]);
yticks(1:n_chan)
yticklabels(chan_lab);
title(['2-4 Hz']);
colormap(jet);
colorbar;

subplot(1,4,2);
imagesc(compareS4_8_abs');
caxis([-500 500]);
yticks(1:n_chan);
yticklabels(chan_lab);
title(['4-8 Hz']);
colorbar;

subplot(1,4,3);
imagesc(compareSBeta_abs');
caxis([-500 500]);
yticks(1:n_chan);
yticklabels(chan_lab);
title(['Beta']);
colorbar;
   
%figure;
subplot(1,4,4);
imagesc(compareSDB_abs');
caxis([-500 500]);
yticks(1:n_chan);
yticklabels(chan_lab);
title(['Delta-Beta']);
colorbar;


figure;
sgtitle(['Wake Absolute',spike]);

subplot(1,4,1);
imagesc(compareW2_4_abs');
caxis([-500 500]);
yticks(1:n_chan);
yticklabels(chan_lab);
title('2-4 Hz');
colormap(jet);
colorbar;


subplot(1,4,2);
imagesc(compareW4_8_abs');
caxis([-500 500]);
yticks(1:n_chan);
yticklabels(chan_lab);
title(['4-8 Hz']);
colorbar;

subplot(1,4,3);
imagesc(compareWBeta_abs');
caxis([-500 500]);
yticks(1:n_chan);
yticklabels(chan_lab);
title(['Beta']);
colorbar;

subplot(1,4,4);
imagesc(compareWDB_abs');
caxis([-500 500]);
yticks(1:n_chan);
yticklabels(chan_lab);
title(['Delta-Beta']);
colorbar;

%% Plotting Delta band

figure;
sgtitle(['Sleep ',spike]);

subplot(1,2,1)
imagesc(compareSDelta');
caxis([-1 1]);
yticks(1:n_chan)
yticklabels(chan_lab);
title(['Delta']);
colormap(jet);
colorbar;

subplot(1,2,2)
imagesc(compareSDelta_abs');
%caxis([-1000 1000]);
yticks(1:n_chan)
yticklabels(chan_lab);
title(['Delta Abs']);
colormap(jet);
colorbar;
%%
for i = 1:n_chan
    DeltaS_first(i) = mean(dataSDelta.(['E',num2str(i)])(:,1)) ;
    DBS_first(i) = mean(dataSDB.(['E',num2str(i)])(:,1));
    DeltaS_last(i) = mean(dataSDelta.(['E',num2str(i)])(:,2)) ;
    DBS_last(i) = mean(dataSDB.(['E',num2str(i)])(:,2));
end

subplot(1,3,3)
imagesc(DeltaS_first');
%caxis([-1000 1000]);
yticks(1:n_chan)
yticklabels(chan_lab);
title(['Delta First']);
colormap(jet);
colorbar;
%%
subplot(1,4,4)
imagesc(DBS_first');
%caxis([-1000 1000]);
yticks(1:n_chan)
yticklabels(chan_lab);
title(['Delta-Beta First']);
colormap(jet);
colorbar;


%%
figure;

subplot(1,2,1)
imagesc(DeltaS_last');
%caxis([-1000 1000]);
yticks(1:n_chan)
yticklabels(chan_lab);
title(['Delta Last']);
colormap(jet);
colorbar;

subplot(1,2,2)
imagesc(DBS_last');
%caxis([-1000 1000]);
yticks(1:n_chan)
yticklabels(chan_lab);
title(['Delta-Beta Last']);
colormap(jet);
colorbar;


%% 

figure

subplot(1,3,1)
imagesc(DeltaS_full');
caxis([0 2500]);
yticks(1:n_chan)
yticklabels(chan_lab);
title(['Delta Whole night Absolute']);
colormap(jet);
colorbar;

subplot(1,3,2)
imagesc(DeltaS_first');
caxis([0 2500]);
yticks(1:n_chan)
yticklabels(chan_lab);
title(['Delta first night Absolute']);
colormap(jet);
colorbar;

subplot(1,3,3)
imagesc(DeltaS_last');
caxis([0 2500]);
yticks(1:n_chan)
yticklabels(chan_lab);
title(['Delta last night Absolute']);
colormap(jet);
colorbar;