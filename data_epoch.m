function [select_epoch, EEG] = data_epoch(EEG, pID, nID, type, Db, TMPREJ)
    EEG.trials = 1;
    EEG.nbchan = length(EEG.data(:,1));
    EEG = epi_log(@pop_select, EEG, 'point', TMPREJ(:,1:2));

    EEG.removed_stretches = TMPREJ(:,1:2);

    DBVals = EEG.data;
    save(['DBVals_',type,pID,nID,'.mat'],'DBVals'); 
    
    %%  Selecting the data from the epochs
    
    Db_sin = single(Db);
    select_epoch = [];

    % Sleep Epoch
    for k = 1:length(DBVals)
        for l = 1:length(Db_sin)
            if Db_sin(l) == DBVals(k)
                select_epoch = [select_epoch, l];
            end
        end
    end

    
end