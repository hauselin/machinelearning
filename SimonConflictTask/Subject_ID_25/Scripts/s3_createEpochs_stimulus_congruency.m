function s3_createEpochs_stimulus_congruency(subject)
%% Epoch data and detect artifact in epochs using EEGLAB artifact detection procedures
% Last modified by Hause Lin 01-08-18 22:08 hauselin@gmail.com
% all runs
% exclude error trials and trials after error trials

%% Set up parameters and read data

% subject = 'AB_06'; % subject id

binlister = 'Stimulus_congruency.txt';
epochRange = [-2500 2500]; % ms
epochBaselineCorrect = [-200 0]; % ms
plotSave = true; % plot and save figures
epochSave = true; % save epoched data?
tfrSave = false; % save time-frequency data?
dataDirectory = 'EEGData'; % directory containing subject/participant data

currentDir = pwd; % current directory
parentDir = currentDir(1:end-8); % remove 'Scripts' from path

% try reading reward csv
try
    readtable(fullfile(parentDir, 'avg_rew_appended', [subject '_averageReward.csv']),'TreatAsEmpty','NA');
catch
    disp(['No reward csv for subject ' subject, '!']);
    dlmwrite(fullfile(parentDir,'Messages',['No reward csv for subject ' subject ' ' datestr(datetime('now')) '.txt']),['No reward csv for subject ' subject ' ' datestr(datetime('now'))],'delimiter','');
    return
end

% add paths
addpath('/psyhome/u4/linhause/matlabtoolboxes/eeglab14_1_2b'); % addpath on UTSC cluster
addpath('/users/hause/dropbox/Apps/MATLAB Toolboxes and Packages/eeglab14_1_2b') % add path on local machine
addpath('/Users/Hause/Dropbox/Apps/MATLAB Toolboxes and Packages/fieldtrip-20180723/')
ft_defaults

%% start analysing data

try

    dataDirectoryAbsPath = fullfile(parentDir, dataDirectory); % full path to data directory
    currentSubjectDirectory = fullfile(dataDirectoryAbsPath, subject); % full path to current subject's data directory

    clc;
    if ~exist(fullfile(currentSubjectDirectory),'dir') % if directory folder doesn't exist, skip this subject
        disp(['Subject ' subject, ' directory is unavailable']);
        return
    end
    disp(['Subject ' subject, '']);

    directoryToReadFrom = fullfile(currentSubjectDirectory,'continuous'); % directory containing correct dataset
    filesInDirectory = dir(directoryToReadFrom); % files in data directory
    fileIdx = find(~cellfun(@isempty, strfind({filesInDirectory.name}, 'ICAPruned_Reref.set'))); % index of file to read

    if length(fileIdx) ~= 1 % if more than one raw data found, error
        error('Check number of files in data directory!');
    end

    % Read data with ICA components removed
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab; % run eeglab
    EEG = pop_loadset('filename',filesInDirectory(fileIdx).name,'filepath',directoryToReadFrom);

    % eeglab redraw
    % check events
    % pop_squeezevents(EEG); % summarize events (ERPLAB)
    % eeg_eventtypes(EEG) % summarize events (EEGLAB)
    close all

    %% Clear evenlist and remove events with codelabel '5sec' (created by adjust plugin)

    EEG.EVENTLIST = [];
    try
        toIncludeEventsIdx = find(cellfun(@isempty,strfind({EEG.event.type},'5sec')));
        EEG.event = EEG.event(toIncludeEventsIdx);
        EEG = eeg_checkset(EEG);
        [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG,EEG,CURRENTSET);
        disp('Removed events with "5sec" label');
    end
    
    %% Reorder channels
    
    newChanOrder = {'Fp1','Fpz','Fp2','AF8','AF4','AFz','AF3','AF7','F7','F5','F3','F1','Fz','F2','F4','F6','F8','FC6','FC4','FC2','FCz','FC1','FC3','FC5','C5','C3','C1','Cz','C2','C4','C6','FT8','T8','TP8','FT7','T7','TP7','CP5','CP3','CP1','CPz','CP2','CP4','CP6','P10','P8','P6','P4','P2','Pz','P1','P3','P5','P7','P9','PO7','PO3','POz','PO4','PO8','O2','Oz','O1','Iz','M1','M2','SO1','LO1','IO1','LO2'};
    EEG = reorderchans(EEG,newChanOrder,true);
    close all

    %% Create eventlist, assign bins, and extract epochs using ERPLAB

    % create eventlist with ERPLAB and save it in results folder
    outPath = fullfile(parentDir,'Results','Eventlist_erplab'); mkdir(outPath);
    EEG = pop_creabasiceventlist(EEG,'AlphanumericCleaning','on','BoundaryNumeric',{-99},'BoundaryString',{'boundary'},'Eventlist',fullfile(outPath,[EEG.subject '_eventlist_raw.txt']));

    % assign bins with ERPLAB binlister
    EEG = pop_binlister(EEG,'BDF',fullfile(parentDir,'Binlister',binlister),'ExportEL',fullfile(outPath,[EEG.subject '_' binlister]),'IndexEL',1,'SendEL2','EEG','UpdateEEG','on','Voutput','EEG');
    
    % get nobaseline correction data to get CNV amplitude
    EEG_nobaselinecorrection = pop_epochbin(EEG,epochRange,'none'); % no baseline correction

    % extract epochs with ERPLAB
    if ~isempty(epochBaselineCorrect)
        EEG = pop_epochbin(EEG,epochRange,epochBaselineCorrect); % extract epochs with baseline correction
    else
        EEG = pop_epochbin(EEG,epochRange,'none'); % no baseline correction
    end
    
    EEG.bindescr = {EEG.EVENTLIST.bdf.description}; % bin description
    EEG = eeg_checkset(EEG);
    [ALLEEG,EEG,CURRENTSET] = eeg_store(ALLEEG,EEG,CURRENTSET);
    disp('Extracted epochs!');

    % Check if time-locking events in epochs match: EEG.EVENTLIST.eventinfo.bini & EEG.EVENTLIST.eventinfo.bepoch
    timelockeventsbepoch = [[EEG.EVENTLIST.eventinfo.bepoch] > 0];
    if sum([EEG.EVENTLIST.eventinfo(timelockeventsbepoch).bini] == -1)
        error('Error! Check EEG.EVENTLIST.eventinfo.bini and EEG.EVENTLIST.eventinfo.bepoch!');
    else
        disp('Epoched data and events match!')
    end
    
    %% Get single-trial CNV amplitude (-200 to 0ms)
    
    cfg = [];
    cfg.chan = 'Cz';
    cfg.chanind = find(strcmpi({EEG.chanlocs.labels},'Cz'));
    cfg.times = [-200 0];
    cfg.timeind = dsearchn(EEG.times',cfg.times');
    cnvAmplitude = squeeze(mean(EEG_nobaselinecorrection.data(cfg.chanind,cfg.timeind(1):cfg.timeind(2),:),2));
    clear cfg
    
    %% Artifact detection on epoched data with EEGLAB

    if exist(fullfile(currentSubjectDirectory,'parameters',[EEG.subject '_bin_' strrep(binlister,'.txt','') '_EEGRejectField.mat']),'file')

        % load artifact detection and rejected epoch information from .mat file
        disp('Not rerunning artifact detection...');
        disp('Overwriting existing EEG.reject with previously saved EEG.reject stored in .mat file...');
        disp([EEG.subject '_bin_' strrep(binlister,'.txt','') '_EEGRejectField.mat']);
        load(fullfile(currentSubjectDirectory,'parameters',[EEG.subject '_bin_' strrep(binlister,'.txt','') '_EEGRejectField.mat']));
        EEG.reject = EEGRejectField;
        EEG = eeg_checkset(EEG);

        % find([EEG.EVENTLIST.eventinfo.flag]);
        EEG = pop_syncroartifacts(EEG,'Direction','eeglab2erplab'); % transfer EEG.reject (eeglab) info to EEG.EVENTLIST (erplab)
        EEG = eeg_checkset(EEG);
        [ALLEEG,EEG,CURRENTSET] = eeg_store(ALLEEG,EEG,CURRENTSET);

    else % run artifact detection

        % Citation: Delorme, A., Sejnowski, T., & Makeig, S. (2007). Enhanced detection of artifacts in EEG data using higher-order statistics and independent component analysis. NeuroImage, 34(4), 1443-1449. doi:10.1016/j.neuroimage.2006.11.004
        disp('Detecting artifacts...');

        % exclude eye and EMG channels from artifact detection
        eyeEmgChans = {'SO1','IO1','SO2','IO2','LO1','LO2','CorsOut','CorsIns','CORRins','CORRout','ZYGup','ZYGlow','COORins','COORout'}; % chan to exclude
        allChannels = {EEG.chanlocs.labels}; % all channel labels
        channelsToDetectArtifact = find(~ismember(allChannels, eyeEmgChans)); % exclude channels above

        % reject by linear trend/variance (max slope, uV/epoch: 100; R-squred limit: 0.5)
        EEG = pop_rejtrend(EEG,1,channelsToDetectArtifact,EEG.pnts,100,0.5,0,0,0);
        disp(['Trials rejected via linear trend/variance: ' num2str(find(EEG.reject.rejconst))]);

        % reject by probability (5 SD)
        EEG = pop_jointprob(EEG,1,channelsToDetectArtifact,5,5,0,0,0,[],0);
        disp(['Trials rejected via probability (5 SD): ' num2str(find(EEG.reject.rejjp))]);

        % reject by spectra (detect muscle and eye movement)
        % muscle: -100 to 100 dB, 20 to 40 Hz
        % eye: -50 to 50 dB, 0 to 2 Hz
        EEG = pop_rejspec(EEG,1,'elecrange',channelsToDetectArtifact,'method','fft','threshold',[-50 50;-100 25],'freqlimits',[0 2;20 40],'eegplotcom','','eegplotplotallrej',0,'eegplotreject',0);
        disp(['Trials rejected via spectra: ' num2str(find(EEG.reject.rejfreq))]);

        % update EEG.reject.rejglobal and EEG.reject.rejglobalE fields with all rejected epochs
        EEG = eeg_rejsuperpose(EEG,1,1,1,1,1,1,1,1);

        % save EEG.reject field as a .mat file so don't need to rerun artifact detection again in the future
        EEGRejectField = EEG.reject;
        save(fullfile(currentSubjectDirectory,'parameters',[EEG.subject '_bin_' strrep(binlister,'.txt','') '_EEGRejectField.mat']),'EEGRejectField'); % save just EEG.reject (maybe this will
        disp('Finished artifact detection and saved EEG.reject field to parameters folder');

        % find([EEG.EVENTLIST.eventinfo.flag]);
        EEG = pop_syncroartifacts(EEG,'Direction','eeglab2erplab'); % transfer EEG.reject (eeglab) info to EEG.EVENTLIST (erplab)
        EEG = eeg_checkset(EEG);
        [ALLEEG,EEG,CURRENTSET] = eeg_store(ALLEEG,EEG,CURRENTSET);
    end

    % save subject's artifact information into subject's parameters directory
    colNames = {'subject','epochs','rejLinear','rejProb','rejSpec','rejTotal','rejPerc','acceptTotal','acceptPerc','binlister'};
    values = {EEG.subject,EEG.trials,length(find(EEG.reject.rejconst)),length(find(EEG.reject.rejjp)),length(find(EEG.reject.rejfreq)),length(find(EEG.reject.rejglobal)),round(length(find(EEG.reject.rejglobal))/EEG.trials*100,2),EEG.trials-length(find(EEG.reject.rejglobal)),round((EEG.trials-length(find(EEG.reject.rejglobal)))/EEG.trials*100,2),binlister};
    T = cell2table(values,'VariableNames',colNames);
    outPath = fullfile(parentDir,'Results','ArtifactInfo'); mkdir(outPath);
    writetable(T, fullfile(outPath,[EEG.subject '_bin_artifactN_' binlister]),'delimiter',',');
    disp('Saved artifact information to parameters directory');

    % gather all subject's artifact detection summary in a table as save as csv file in Binlister directory
    % tempFiles = dir(fullfile(dataDirectoryAbsPath,'**','parameters','',['*_bin_artifactN_' binlister])); % find matching files recursively (only works for MATLAB 2016b onwards)
    % gather all subject's artifact detection summary in a table as save as csv file in Binlister directory
    artifactTableAllSubjects = table();
    tempFiles = dir(fullfile(outPath,['*_bin_artifactN_' binlister])); % find matching files recursively
    filePaths = {tempFiles.name};
    for i=1:length(filePaths)
        tempData = readtable(fullfile(outPath,filePaths{i}));
        artifactTableAllSubjects = [artifactTableAllSubjects; tempData];
    end
    artifactTableAllSubjects.reject = artifactTableAllSubjects.rejPerc > 25; % if > 25% trials rejected, mark subject to reject
    artifactTableAllSubjects = [table([1:height(artifactTableAllSubjects)]') artifactTableAllSubjects]; % add subject number count
    artifactTableAllSubjects.Properties.VariableNames{1} = 'subjectN';

    writetable(artifactTableAllSubjects,fullfile(outPath,['_artifactSummary_' strrep(binlister,'.txt','') '.csv']));
    disp('Saved all subject''s artifact information to directory');

    %% Save single-trial event information

    el = EEG.EVENTLIST.eventinfo; % save ERPLAB eventlist structure
    timeLockedEventsIdx = [el.bepoch] ~= 0; % find indices with time-locked events (bepoch is not 0)
    allEvents = {el.binlabel}; % all event labels (B1, B2 etc.)
    epochEvent = allEvents(timeLockedEventsIdx); % event/code for every epoch (including epochs with artifacts)
    % find artifact and no-artifact epochs
    artifactFlags = [el.flag];
    cleanEpochs = find(~artifactFlags(timeLockedEventsIdx)); % epoch (indices) without artifacts

    if ((length(cleanEpochs) + sum(artifactFlags)) ~= EEG.trials) || (length(epochEvent) ~= EEG.trials) % my manual count of trials and EEG.trials should match
        error('Error! Trial numbers don''t match! Double check!');
    end

    % Save single-trial info (design matrix) as table
    elT = struct2table(el);
    elT.enable = []; elT.dura = []; elT.codelabel = []; % remove variables
    % add subject id to table
    C = cell(size(elT,1),1); % empty cell to store subject id
    C(:) = {EEG.subject}; % fill empty cell with subject id
    elT = [C elT]; % join cell with table
    elT.Properties.VariableNames = {'subject','eventN','eventCode','binlabel','timeS','samplingPoint','artifactFlag','binIndicator','epochN'}; % rename variable names

    % add average reward rate info
    rewardRate = readtable(fullfile(parentDir, 'avg_rew_appended', [subject '_averageReward.csv']),'TreatAsEmpty','NA');
    rewardRate = rewardRate(:, {'subject','eventN','eventCode','samplingPoint','avg_reward','reward','iti','trial_num','run_num','trial_in_run','rt','congruencyDC','accDC','previousCongruencyDC','previousAcc','trialType','keyRep'});
    elT = outerjoin(elT,rewardRate,'Type','left','Mergekeys',true);  % join

    elT.bindescr = elT.binlabel; % assign bindescr to each epoch
    for i=1:size(elT,1) %  for each event...
        bindescrTemp = elT{i,'bindescr'}{1}; % get binlabel
        if strcmpi(bindescrTemp,'""') % if empty, it's not time-locking event
            elT{i,'bindescr'}{1} = '';
        else  % if not empty, fill in bindescr from ERP
            elT{i,'bindescr'}{1} = EEG.bindescr{elT{i,'binIndicator'}};
        end
    end

    % elT.rt = [NaN; diff(elT.timeS)]; % compute RT (change this line accordingly for different designs/experiments)
    % elT.rt(elT.artifactFlag == 1) = NaN; % convert artifact epochs to NaN
    el_timelockevent = elT((elT.epochN ~= 0),:); % eventlist with just time-locked event info
    
    % add CNV amplitude to eventinfo
    if size(el_timelockevent,1) == length(cnvAmplitude)
        el_timelockevent = [el_timelockevent table(cnvAmplitude)];
    else
        error('Check CNV amplitude dimensions!');
    end
    
    % save data
    outPath = fullfile(parentDir,'Results','TrialEventInfo',strrep(binlister,'.txt','')); mkdir(outPath);
    writetable(elT,fullfile(outPath,[EEG.subject '_eventAllEvents.csv'])) % all events
    writetable(el_timelockevent,fullfile(outPath,[EEG.subject '_eventTimeLockedEvents.csv'])) % only time-locked events
    disp('Saved event list as txt and csv files');
    
    %% Create design matrix or array of events

    % variable names to save together with matfile later on
    % only variables of integer or double type!
    epochs_allVars = {'epochN','artifactFlag','binIndicator','avg_reward','reward','iti','trial_num','run_num','trial_in_run','rt','congruencyDC','accDC','previousCongruencyDC','previousAcc','trialType','keyRep','cnvAmplitude'};
    epochs_all = el_timelockevent(:,epochs_allVars);
    % ensure all variables are double
    for i = 1:length(epochs_allVars)
       if class(epochs_all.(i)) ~= 'double'
           epochs_all.(i) = double(epochs_all.(i));
       end
    end

    % create a matrix indicating bin for each trial
    epochs_all = table2array(epochs_all);
    if sum(epochs_all(:,1) < 0) || sum(epochs_all(:,3) < 0)
        error('Error! Check design matrix epoch number and bin number!');
    end

    epochs_all(:,2) = epochs_all(:,2) == 0; % convert artifacts rejected to artifacts accepted
    epochs_clean = epochs_all(epochs_all(:,2) == 1,:); % select only epochs without artifacts
    epochs_clean = epochs_clean((sum(isnan(epochs_clean),2) == 0),:); % select only epochs without NaN

    EEG.epochs_all = epochs_all;
    EEG.epochs_clean = epochs_clean;
    EEG.epochs_allVars = epochs_allVars;

    % epochs_all and epochs_clean
    outPath = fullfile(parentDir,'Results','DesignMatrix',strrep(binlister,'.txt','')); mkdir(outPath);
    dlmwrite(fullfile(outPath,[EEG.subject '_epochsAll_designMat.txt']),epochs_all,'delimiter','\t')
    save(fullfile(outPath,[EEG.subject '_epochsAll_designMat.mat']),'epochs_all','epochs_allVars');
    dlmwrite(fullfile(outPath,[EEG.subject '_epochsClean_designMat.txt']),epochs_clean,'delimiter','\t')
    save(fullfile(outPath,[EEG.subject '_epochsClean_designMat.mat']),'epochs_clean','epochs_allVars');
    disp('Saved design matrix as txt and mat files');

    %% Select only clean and relevant epochs

    % select only relevant trials/epochs
    % blocks1to5Idx = find(ismember(epochs_clean(:,ismember(epochs_allVars,'block')),[1:10]));
    accurateTrialsIdx = find(epochs_clean(:,ismember(epochs_allVars,'accDC')) == 0); % current trial correct
    % preAccTrialsIdx = find(epochs_clean(:,ismember(epochs_allVars,'previousAcc')) == 0); % previous trial correct
    % relevantIdx = intersect(intersect(blocks1to5Idx,accurateTrialsIdx),preAccTrialsIdx);
    % relevantIdx = intersect(preAccTrialsIdx,accurateTrialsIdx);
    relevantIdx = accurateTrialsIdx;
    epochsRelevant = epochs_clean(relevantIdx,:);
    epochsRelevantIdx = epochsRelevant(:,1);
    
    EEG.epochs_allVars = epochs_allVars;
    EEG.epochs_all = epochs_all;
    EEG.epochs_clean = epochsRelevant;
    EEG.history = [];
    
    ALLEEG = [];
    ALLEEG = EEG;
    ALLEEG(2) = pop_select(ALLEEG(1),'trial',epochsRelevantIdx); % select only clean epochs
    
    ALLEEG(1).setname = [EEG.subject '_all_epochs'];
    ALLEEG(1).filename = ALLEEG(1).setname;
    ALLEEG(2).setname = [EEG.subject '_clean_epochs'];
    ALLEEG(2).filename = ALLEEG(2).setname;
    
    if epochSave
        outPath = fullfile(parentDir,'Results','Epoch_stimulus_correctTrials'); mkdir(outPath); %!!! REMEMBER TO CHANGE THIS!
        pop_saveset(ALLEEG(1),'filename',[EEG.subject '_epochsAll.set'],'filepath',outPath);
    end

    %% Create and save ERPset using ERPLAB

    ERP = pop_averager(ALLEEG(2),'Criterion','good','SEM','on');
    outPath = fullfile(parentDir,'Results','ERP_stimulus_correctTrials'); mkdir(outPath); %!!! REMEMBER TO CHANGE THIS!
    pop_savemyerp(ERP,'erpname',[EEG.subject '_bin_' strrep(binlister,'.txt','')],'filename',[EEG.subject '.erp'],'filepath',outPath,'Warning','off');
    % eeglab redraw
    % erplab redraw

    %% Convert all epochs to fieldtrip format for further processing

    EEGft = eeglab2fieldtrip(ALLEEG(2),'preprocessing');
    epochs_clean = epochsRelevant;

    %% Compute and plot ERPs in fieldtrip

    % average of all trials
    cfg = [];
    erp_alltrials = ft_timelockanalysis(cfg,EEGft);

    if true
        % erp for each bin
        binsUnique = sort(unique(epochs_clean(:,3)));
        erp_bins = cell(1,length(binsUnique));
        for binI = 1:length(erp_bins)
            cfg = [];
            cfg.trials = find(epochs_clean(:,3) == binsUnique(binI)); % select trials
            erp_bins{binI} = ft_timelockanalysis(cfg,EEGft); % compute mean
        end

        if plotSave
            cfg = [];
            cfg.xlim = [-0.5 1.0];
            cfg.showlabels = 'yes';
            cfg.comment = 'stimulus-locked congruent vs. incongruent (correct-only trials)';
            cfg.showcomment = 'yes';
            cfg.rotate = 90; % rotate plot
            cfg.fontsize = 16;
            cfg.linewidth = 0.5;
            cfg.layout = 'biosemi64.lay';
            figure(5)
            ft_multiplotER(cfg,erp_bins{:});
            set(gcf,'units','normalized','outerposition',[0 0 1 1],'PaperPositionMode','auto','DefaultTextInterpreter','none','PaperOrientation','portrait'); % maximize figure
            figPath = fullfile(parentDir,'Figures','ERP_stimulus_correctTrials'); mkdir(figPath); %!!! REMEMBER TO CHANGE THIS!
            print(gcf,'-djpeg','-r200', fullfile(figPath,[EEG.subject '_conditions.jpg']));
        end
        close all

    end
    
    %% Single-trial erp regression
    
    % create design matrix
    epochs_allVars;
    regressors = {'avg_reward' 'reward' 'iti' 'trial_num' 'congruencyDC' 'previousCongruencyDC' 'previousAcc' 'keyRep'};
    colIdx = ismember(epochs_allVars,regressors);
    designMat = epochs_clean(:,colIdx);
    % z-score
    designMat(:,1) = zscore(designMat(:,1));
    designMat(:,2) = zscore(designMat(:,2));
    designMat(:,3) = zscore(designMat(:,3));
    designMat(:,4) = zscore(designMat(:,4));
    
    cfg = [];
    cfg.intercept = 'yes';
    cfg.times = [-0.2:0.2:1.0];
    cfg.n_permutes = 500;
    cfg.model = 'y ~ b0 + avg_reward + reward + iti + trial_num + congruencyDC + previousCongruencyDC + previousAcc + keyRep';
    cfg.type = 'zpermute';
    cfg.designmatrix = designMat;
    erpBeta = ft_regressERF_permute(cfg,ALLEEG(2));
    erpBeta
    
    regressorPlot = [2];
    regressorNames = {'avgReward'};
    figure(101), clf
    colormap viridis
    set(gcf,'units','normalized','outerposition',[0 0 1 1],'PaperPositionMode','auto','DefaultTextInterpreter','none','PaperOrientation','portrait'); % maximize figure
    for regressorI = 1:length(regressorPlot)
        subplot(1,2,regressorI)
        imagesc(erpBeta.cfg.times,1:length(erpBeta.cfg.chan),squeeze(erpBeta.realbeta(regressorPlot(regressorI),:,:)))
        colorbar
        set(gca,'clim',[-2 2])
        set(gca,'ytick',1:length(erpBeta.cfg.chan),'yticklabel',erpBeta.cfg.chan)
        title(['beta ' regressorNames{regressorI}])
        
        subplot(1,2,regressorI+1)
        imagesc(erpBeta.cfg.times,1:length(erpBeta.cfg.chan),squeeze(erpBeta.zmap(regressorPlot(regressorI),:,:)))
        colorbar
        set(gca,'clim',[-2 2])
        set(gca,'ytick',1:length(erpBeta.cfg.chan),'yticklabel',erpBeta.cfg.chan)
        title(['permuted z map ' regressorNames{regressorI}])
    end
    
    figPath = fullfile(parentDir,'Figures','ERP_stimulus_correctTrials'); mkdir(figPath);
    print(gcf,'-djpeg','-r200', fullfile(figPath,[EEG.subject '_beta_avgReward.jpg']));
    
    regressorPlot = [3];
    regressorNames = {'reward'};
    figure(101), clf
    colormap viridis
    set(gcf,'units','normalized','outerposition',[0 0 1 1],'PaperPositionMode','auto','DefaultTextInterpreter','none','PaperOrientation','portrait'); % maximize figure
    for regressorI = 1:length(regressorPlot)
        subplot(1,2,regressorI)
        imagesc(erpBeta.cfg.times,1:length(erpBeta.cfg.chan),squeeze(erpBeta.realbeta(regressorPlot(regressorI),:,:)))
        colorbar
        set(gca,'clim',[-2 2])
        set(gca,'ytick',1:length(erpBeta.cfg.chan),'yticklabel',erpBeta.cfg.chan)
        title(['beta ' regressorNames{regressorI}])
        
        subplot(1,2,regressorI+1)
        imagesc(erpBeta.cfg.times,1:length(erpBeta.cfg.chan),squeeze(erpBeta.zmap(regressorPlot(regressorI),:,:)))
        colorbar
        set(gca,'clim',[-2 2])
        set(gca,'ytick',1:length(erpBeta.cfg.chan),'yticklabel',erpBeta.cfg.chan)
        title(['permuted z map ' regressorNames{regressorI}])
    end
    
    figPath = fullfile(parentDir,'Figures','ERP_stimulus_correctTrials'); mkdir(figPath);
    print(gcf,'-djpeg','-r200', fullfile(figPath,[EEG.subject '_beta_reward.jpg']));
    
    
    regressorPlot = [6];
    regressorNames = {'congruency'};
    figure(101), clf
    colormap viridis
    set(gcf,'units','normalized','outerposition',[0 0 1 1],'PaperPositionMode','auto','DefaultTextInterpreter','none','PaperOrientation','portrait'); % maximize figure
    for regressorI = 1:length(regressorPlot)
        subplot(1,2,regressorI)
        imagesc(erpBeta.cfg.times,1:length(erpBeta.cfg.chan),squeeze(erpBeta.realbeta(regressorPlot(regressorI),:,:)))
        colorbar
        set(gca,'clim',[-2 2])
        set(gca,'ytick',1:length(erpBeta.cfg.chan),'yticklabel',erpBeta.cfg.chan)
        title(['beta ' regressorNames{regressorI}])
        
        subplot(1,2,regressorI+1)
        imagesc(erpBeta.cfg.times,1:length(erpBeta.cfg.chan),squeeze(erpBeta.zmap(regressorPlot(regressorI),:,:)))
        colorbar
        set(gca,'clim',[-2 2])
        set(gca,'ytick',1:length(erpBeta.cfg.chan),'yticklabel',erpBeta.cfg.chan)
        title(['permuted z map ' regressorNames{regressorI}])
    end
    
    figPath = fullfile(parentDir,'Figures','ERP_stimulus_correctTrials'); mkdir(figPath);
    print(gcf,'-djpeg','-r200', fullfile(figPath,[EEG.subject '_beta_congruency.jpg']));
    
    % set up strucutre to plot stuff
    erpBetaCell = {};
    tempBeta = erp_alltrials;
    
    tempBeta.avg = squeeze(erpBeta.zmap(2,:,:));
    tempBeta.time = erpBeta.cfg.times
    tempBeta.var = [];
    tempBeta.dof = [];
    erpBetaCell{1} = tempBeta;
    
    tempBeta.avg = squeeze(erpBeta.zmap(3,:,:));
    erpBetaCell{2} = tempBeta;
    
    tempBeta.avg = squeeze(erpBeta.zmap(6,:,:));
    erpBetaCell{3} = tempBeta;
    
    cfg = [];
    cfg.xlim = [-0.2 0.8];
    cfg.showlabels = 'yes';
    cfg.comment = 'stimulusReward locked (correct trials): avgReward, reward, congruency';
    cfg.showcomment = 'yes';
    cfg.rotate = 90; % rotate plot
    cfg.fontsize = 16;
    cfg.linewidth = 0.5;
    cfg.layout = 'biosemi64.lay';
    figure(5)
    ft_multiplotER(cfg,erpBetaCell{:});
    set(gcf,'units','normalized','outerposition',[0 0 1 1],'PaperPositionMode','auto','DefaultTextInterpreter','none','PaperOrientation','portrait'); % maximize figure
    figPath = fullfile(parentDir,'Figures','ERP_stimulus_correctTrials'); mkdir(figPath);
    print(gcf,'-djpeg','-r200', fullfile(figPath,[EEG.subject '_betaERPs.jpg']));
    
    close all
    
    %% Save ERP regresion beta coefficients
    
    erpBetaCell = {};
    erpBetaCell{1} = erpBeta;
    erpBetaCell
    erpBetaCell{1}
    outPath = fullfile(parentDir,'Results','ERP_stimulus_beta'); mkdir(outPath);
    save(fullfile(outPath,[EEG.subject '.mat']),'erpBetaCell');
    
    %% Time-frequency decomposition in fieldtrip (keep single trials)

    cfg = [];
    % cfg.channel = {'all','-M1','-M2','-SO1','-LO1','-IO1','-LO2'}; % all or channels in cell array
    cfg.channel = {'FCz','Cz','Pz','Fpz','Fz','T7','T8','PO3','PO4','Oz','CP3','CP4','F5','F6'}; % all or channels in cell array
    cfg.method = 'wavelet'; % time-freq decomposition method
    cfg.foi = 1:40;	% freq of interest
    cfg.width = logspace(log10(3),log10(12),length(cfg.foi)); % cycles
    cfg.output = 'fourier';	% fourier, pow, or powandcsd
    cfg.ds = 0.02; % downsampling spacing in seconds
    cfg.ds2 = (1/EEGft.fsample)*(round(0.02/(1/EEGft.fsample))); % actual downsampling spacing based on sampling rate
    cfg.toistartidx = dsearchn(EEGft.time{1}',-1.0); % specify output time begin in seconds
    cfg.toiendidx = dsearchn(EEGft.time{1}',1.0); % specify output time end in seconds
    cfg.toi = EEGft.time{1}(cfg.toistartidx):cfg.ds2:EEGft.time{1}(cfg.toiendidx); % downsampled timepoints to return 
    cfg.keeptrials = 'yes'; % return single trials (yes) or average (no)
    cfg.trials = 'all'; % all trials
    cfg.pad = 'nextpow2'; % pad data to up to next power of 2 (faster more efficient fft)
    % tf_fourierSpec = ft_freqanalysis(cfg,EEGft);

    % run or reload previously saved data
    outPath = fullfile(parentDir,'Results','TFR_fourierSpec_stimulus_correctTrials'); % directory to save data to or read data from
    outFile = fullfile(outPath,[EEG.subject '_singleTrials.mat']); % output file name 
    if exist(outFile) % if exist, load it
        disp(['Loading data...']);
        disp(outFile);
        load(outFile);
    else % if don't exist, run analysis
        tf_fourierSpec = ft_freqanalysis(cfg,EEGft); % run analysis
        disp('Saving data...');
        disp(outFile);
        mkdir(outPath);
        save(outFile,'tf_fourierSpec','epochs_allVars','epochs_all','epochs_clean');
    end

    % prepare objects for each single-trial regression model later on
    tfPowerZmap = cell(0); % store regression models later on
    
    %% Convert fourier spectrum to power
    
    tf_pow = tf_fourierSpec;
    tf_pow.powspctrm = abs(tf_fourierSpec.fourierspctrm).^2;
    tf_pow = rmfield(tf_pow,'fourierspctrm');
   
    %% Single-trial regression: y ~ b0 + avg_reward + reward + congruencyDC + iti + trial_num + keyRep + previousCongruencyDC

    % create design matrix
    epochs_allVars;
    regressors = {'avg_reward' 'reward' 'iti' 'trial_num' 'congruencyDC' 'previousCongruencyDC' 'previousAcc' 'keyRep'};
    % regressors = {'avg_reward' 'reward' 'iti' 'trial_num' 'congruencyDC' 'previousCongruencyDC' 'keyRep'};
    colIdx = ismember(epochs_allVars,regressors);
    designMat = epochs_clean(:,colIdx);
    
    % z-score
    designMat(:,1) = zscore(designMat(:,1));
    designMat(:,2) = zscore(designMat(:,2));
    designMat(:,3) = zscore(designMat(:,3));
    designMat(:,4) = zscore(designMat(:,4));

    % fit model
    cfg = [];
    cfg.designmatrix = designMat;
    cfg.model = 'y ~ avg_reward + reward + iti + trial_num + congruencyDC + previousCongruencyDC + previousAcc + keyRep';
    cfg.intercept = 'no';
    cfg.transformY = 'tiedrank';
    cfg.type = 'zpermute';
    cfg.n_permutes = 500;
    regressOutput = ft_regressTFRPower_permute(cfg,tf_pow);
   
    % cluster correction
%     cfg = [];
%     cfg.chan = 'FCz';
%     cfg.regressorlabels = {'avgReward' 'reward' 'iti' 'trialNum' 'congruency' 'previousCongruency' 'previousAcc' 'keyRep'}
%     ft_regressTFRPower_permute_clusterplot(cfg,regressOutput);
%     
%     outPath = fullfile(parentDir,'Figures','TFR_pow_stimulus_correctTrials_beta'); mkdir(outPath);
%     tempFigname = fullfile(outPath,[EEG.subject '_FCz.jpg']);
%     print(gcf,'-djpeg','-r100', tempFigname);

    tfPowerZmap{1} = regressOutput;
    clear regressOutput
    
    %% Save z/beta values

    tfPowerZmap;
    outPath = fullfile(parentDir,'Results','TFR_pow_stimulus_correctTrials_beta_avgReward'); mkdir(outPath);
    save(fullfile(outPath,[EEG.subject '.mat']),'tfPowerZmap');

    %% MultiplotTFR
    
    if plotSave
        fignames = {'avgReward'};
        for j = 1:length(tfPowerZmap) % for each fitted model, plot z/beta values
            cfg = [];
            cfg.zMax = 2;
            cfg.zlim = [-cfg.zMax cfg.zMax];
            cfg.comment = fignames{j};
            cfg.showlabels = 'yes';
            cfg.parameter = 'powspctrm';
            cfg.layout = 'biosemi64.lay';
            cfg.colormap = 'viridis';
            cfg.colorbar = 'yes';
            cfg.fontsize = 14;
            cfg.trials = 1; % coefficient number
            close all; figure(201); clf
            ft_multiplotTFR(cfg,tfPowerZmap{j});
            set(gcf,'units','normalized','outerposition',[0 0 1 1],'PaperPositionMode','auto','DefaultTextInterpreter','none','PaperOrientation','portrait'); % maximize figure
            outPath = fullfile(parentDir,'Figures','TFR_pow_stimulus_correctTrials_beta'); mkdir(outPath);
            tempFigname = fullfile(outPath,[EEG.subject '_multiplotTF_' fignames{j} '.jpg']);
            print(gcf,'-djpeg','-r100', tempFigname);
            close all
        end
    end
    
    if plotSave
        fignames = {'congruency'};
        for j = 1:length(tfPowerZmap) % for each fitted model, plot z/beta values
            cfg = [];
            cfg.zMax = 2;
            cfg.zlim = [-cfg.zMax cfg.zMax];
            cfg.comment = fignames{j};
            cfg.showlabels = 'yes';
            cfg.parameter = 'powspctrm';
            cfg.layout = 'biosemi64.lay';
            cfg.colormap = 'viridis';
            cfg.colorbar = 'yes';
            cfg.fontsize = 14;
            cfg.trials = 5; % coefficient number
            figure(201); clf
            ft_multiplotTFR(cfg,tfPowerZmap{j});
            set(gcf,'units','normalized','outerposition',[0 0 1 1],'PaperPositionMode','auto','DefaultTextInterpreter','none','PaperOrientation','portrait'); % maximize figure
            outPath = fullfile(parentDir,'Figures','TFR_pow_stimulus_correctTrials_beta'); mkdir(outPath);
            tempFigname = fullfile(outPath,[EEG.subject '_multiplotTF_' fignames{j} '.jpg']);
            print(gcf,'-djpeg','-r100', tempFigname);
            close all
        end
    end

    %% plot FCz

    if plotSave
        figure(101); clf
        set(gcf,'units','normalized','outerposition',[0 0 1 1],'PaperPositionMode','auto','DefaultTextInterpreter','none','PaperOrientation','portrait'); % maximize figure
        for j = 1:length(tfPowerZmap)
            cfg = [];
            cfg.chan = 'FCz';
            cfg.trials = 1; % coefficient number
            cfg.subplot = [1,2,1];
            cfg.fontsize = 14;
            cfg.colormap = 'viridis';
            cfg.title = {cfg.chan,tfPowerZmap{j}.cfg.regressparams.model 'avgReward'};
            ft_contourfTFR(cfg,tfPowerZmap{j});
            set(0,'DefaultTextInterpreter','none') % disable TEX
            
            cfg = [];
            cfg.chan = 'FCz';
            cfg.trials = 5; % coefficient number
            cfg.subplot = [1,2,2];
            cfg.fontsize = 14;
            cfg.colormap = 'viridis';
            cfg.title = {cfg.chan,tfPowerZmap{j}.cfg.regressparams.model 'congruency'};
            ft_contourfTFR(cfg,tfPowerZmap{j});
            set(0,'DefaultTextInterpreter','none') % disable TEX
        end
        outPath = fullfile(parentDir,'Figures','TFR_pow_stimulus_correctTrials_beta'); mkdir(outPath);
        tempFigname = fullfile(outPath,[EEG.subject '_FCz_avgReward_congruency.jpg']);
        print(gcf,'-djpeg','-r100', tempFigname);
        close all
    end

    %% Save parameters

    params.workingdirectory = parentDir;
    params.binlister = fullfile(parentDir,'Binlister',binlister);
    params.epochRange = epochRange; % ms
    params.epochBaselineCorrect = epochBaselineCorrect;

    save(fullfile(parentDir,'Binlister',[strrep(binlister,'.txt','') '_epochingParams.mat']),'params');

    %% Done

    disp(['Finished subject ' subject, '!']);
    dlmwrite(fullfile(parentDir,'Messages',['Finished_Subject_' subject ' ' datestr(datetime('now')) '.txt']),['Finished!' ' ' datestr(datetime('now'))],'delimiter','');

catch ME % if errors when runnig script, catch them
    disp(['Error with subject ' subject, '!']);
    save(fullfile(parentDir,'Messages',['Error_MException_' subject ' ' datestr(datetime('now')) '.mat']),'ME'); % save MException object
end

clear
close all

end % end function
