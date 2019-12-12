% Load the file of the patient;
%x = readtable('/Users/bsevak/Documents/AssessSleep.xlsx');

%for p = 7%:length(x.Patient)


Patient = 'Patient_1' %char(x.Patient(p));
Night = 'Night_1'%char(x.Night(p));

cd(['/Users/bsevak/Documents/Merged Data_BF/Merged_Data/',Patient,'/',Night,'/'])
load([Patient,'_',Night,'_200Hz_resampled.mat']);

% Setting Parameters according to the Assess the Sleep algorithm
if exist('merged_matfile')
    Data = merged_matfile';
elseif exist('file_new')
    Data = file_new';
end

El_number = double(1:length(Data(:,1)));
El_name = {header_ns3(:).Label};

fs = 200;

cd(['/Users/bsevak/Documents/Merged Data_BF/Merged_Data/',Patient,'/',Night,'/',Patient,'_',Night]);
save([Patient,'_',Night,'_data.mat'],'Data','El_number','El_name','fs','-v7.3');


%%
cd(['/Users/bsevak/Documents/Merged Data_BF/Merged_Data/',Patient,'/',Night,'/'])
global subject_id z;
subject_id  = {[Patient,'_',Night]};
z = 1;

Assess_the_sleep_NST();

clearvars -except x

%end