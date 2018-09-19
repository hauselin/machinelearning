function s4_stimulus(subject)

%% set up parameters and read data

% subject = 'AB_06'; % subject id
dataDir = 'TFR_fourierSpec_stimulus_correctTrials';
currentDir = pwd; % current directory
parentDir = currentDir(1:end-8); % remove 'Scripts' from path
dataAbsDir = fullfile(parentDir,'Results',dataDir);
filename = dir(fullfile(dataAbsDir,[subject '*'])); % find subject data file

% load data
load(fullfile(dataAbsDir,filename.name))
tf_fourierSpec
size(tf_fourierSpec.fourierspctrm,1) % trials

%% load design matrix

load(fullfile(parentDir,'Results','DesignMatrix','Stimulus_congruency',[subject '_epochsClean_designMat.mat']))
epochs_allVars
size(epochs_clean)

% select only correct trials
accurateTrialsIdx = find(epochs_clean(:,ismember(epochs_allVars,'accDC')) == 0); % current trial correct
epochs_clean = epochs_clean(accurateTrialsIdx,:);
size(epochs_clean)

%% Convert fourier spectrum to power
    
if false
    tf_pow = tf_fourierSpec;
    tf_pow.powspctrm = abs(tf_fourierSpec.fourierspctrm).^2;
    tf_pow = rmfield(tf_pow,'fourierspctrm');
end
    
%% single-trial regression: power ~ congruencyDC

if false 
    % create design matrix
    epochs_allVars;
    regressors = {'congruencyDC'};
    colIdx = ismember(epochs_allVars,regressors);
    designMat = epochs_clean(:,colIdx);

    % fit model
    cfg = [];
    cfg.designmatrix = designMat;
    cfg.model = 'y ~ congruencyDC';
    cfg.intercept = 'no';
    cfg.transformY = 'tiedrank';
    cfg.type = 'zpermute';
    cfg.n_permutes = 500;
    regressOutput = ft_regressTFRPower_permute(cfg,tf_pow);
    regressOutput

    %% save results

    outPath = fullfile(parentDir,'Results','TFR_pow_stimulus_correctTrials_beta_congruency'); mkdir(outPath);
    save(fullfile(outPath,[subject '.mat']),'regressOutput');

    %% plot

    close all
    cfg = [];
    cfg.zMax = 2;
    cfg.zlim = [-cfg.zMax cfg.zMax];
    cfg.comment = regressOutput.cfg.regressparams.model;
    cfg.showlabels = 'yes';
    cfg.parameter = 'powspctrm';
    cfg.layout = 'biosemi64.lay';
    cfg.colormap = 'viridis';
    cfg.colorbar = 'yes';
    cfg.fontsize = 14;
    cfg.trials = 1; % coefficient number
    close all; figure(201); clf
    ft_multiplotTFR(cfg,regressOutput);

    set(gcf,'units','normalized','outerposition',[0 0 1 1],'PaperPositionMode','auto','DefaultTextInterpreter','none','PaperOrientation','portrait'); % maximize figure
    outPath = fullfile(parentDir,'Figures','TFR_pow_stimulus_correctTrials_beta'); mkdir(outPath);
    tempFigname = fullfile(outPath,[subject '_multiplotTF_congruency.jpg']);
    print(gcf,'-djpeg','-r100', tempFigname);
    close all


    cfg = [];
    cfg.chan = 'FCz';
    cfg.trials = 1; % coefficient number
    cfg.fontsize = 14;
    cfg.colormap = 'viridis';
    cfg.title = {cfg.chan,regressOutput.cfg.regressparams.model 'congruency'};
    ft_contourfTFR(cfg,regressOutput);

    set(0,'DefaultTextInterpreter','none') % disable TEX
    set(gcf,'units','normalized','outerposition',[0 0 1 1],'PaperPositionMode','auto','DefaultTextInterpreter','none','PaperOrientation','portrait'); % maximize figure
    outPath = fullfile(parentDir,'Figures','TFR_pow_stimulus_correctTrials_beta'); mkdir(outPath);
    tempFigname = fullfile(outPath,[subject '_FCz_congruency.jpg']);
    print(gcf,'-djpeg','-r100', tempFigname);
    close all
end

%% run analyses for only first 5 runs

if false 
    % select first half of data
    epochs_allVars
    runs1to5 = find(ismember(epochs_clean(:,ismember(epochs_allVars,'run_num')),[1:5])); % indices
    cfg = [];
    cfg.trials = runs1to5;
    tf_pow50 = ft_selectdata(cfg,tf_pow);

    % create design matrix
    epochs_allVars;
    regressors = {'avg_reward' 'reward' 'iti' 'trial_num' 'congruencyDC' 'previousCongruencyDC' 'previousAcc' 'keyRep'};
    % regressors = {'avg_reward' 'reward' 'iti' 'trial_num' 'congruencyDC' 'previousCongruencyDC' 'keyRep'};
    colIdx = ismember(epochs_allVars,regressors);
    designMat = epochs_clean(runs1to5,colIdx);

    % z-score
    designMat(:,1) = zscore(designMat(:,1));
    designMat(:,2) = zscore(designMat(:,2));
    designMat(:,3) = zscore(designMat(:,3));
    designMat(:,4) = zscore(designMat(:,4));

    % fit model
    cfg = [];
    cfg.designmatrix = designMat;
    cfg.model = 'y ~ avg_reward + reward + iti + trial_num + congruencyDC + previousCongruencyDC + previousAcc + keyRep (50%)';
    cfg.intercept = 'no';
    cfg.transformY = 'tiedrank';
    cfg.type = 'zpermute';
    cfg.n_permutes = 500;
    regressOutput = ft_regressTFRPower_permute(cfg,tf_pow50);
     
    tfPowerZmap{1} = regressOutput;
    clear regressOutput
    
    %% Save z/beta values

    tfPowerZmap;
    outPath = fullfile(parentDir,'Results','TFR_pow_stimulus_correctTrials_beta_avgReward50'); mkdir(outPath);
    save(fullfile(outPath,[subject '.mat']),'tfPowerZmap');
   
end

%% single-trial regression: power ~ rt

if false
    % create design matrix
    epochs_allVars;
    regressors = {'rt'};
    colIdx = ismember(epochs_allVars,regressors);
    designMat = epochs_clean(:,colIdx);
    
    % z-score
    designMat(:,1) = zscore(designMat(:,1));
    
    % fit model
    cfg = [];
    cfg.designmatrix = designMat;
    cfg.model = 'y ~ rt';
    cfg.intercept = 'no';
    cfg.transformY = 'tiedrank';
    cfg.type = 'zpermute';
    cfg.n_permutes = 500;
    regressOutput = ft_regressTFRPower_permute(cfg,tf_pow);
    regressOutput
    
    %% save results

    outPath = fullfile(parentDir,'Results','TFR_pow_stimulus_correctTrials_beta_rt'); mkdir(outPath);
    save(fullfile(outPath,[subject '.mat']),'regressOutput');
    
    %% plot

    close all
    cfg = [];
    cfg.zMax = 2;
    cfg.zlim = [-cfg.zMax cfg.zMax];
    cfg.comment = regressOutput.cfg.regressparams.model;
    cfg.showlabels = 'yes';
    cfg.parameter = 'powspctrm';
    cfg.layout = 'biosemi64.lay';
    cfg.colormap = 'viridis';
    cfg.colorbar = 'yes';
    cfg.fontsize = 14;
    cfg.trials = 1; % coefficient number
    close all; figure(201); clf
    ft_multiplotTFR(cfg,regressOutput);

    set(gcf,'units','normalized','outerposition',[0 0 1 1],'PaperPositionMode','auto','DefaultTextInterpreter','none','PaperOrientation','portrait'); % maximize figure
    outPath = fullfile(parentDir,'Figures','TFR_pow_stimulus_correctTrials_beta'); mkdir(outPath);
    tempFigname = fullfile(outPath,[subject '_multiplotTF_rt.jpg']);
    print(gcf,'-djpeg','-r100', tempFigname);
    close all

    cfg = [];
    cfg.chan = 'FCz';
    cfg.trials = 1; % coefficient number
    cfg.fontsize = 14;
    cfg.colormap = 'viridis';
    cfg.title = {cfg.chan,regressOutput.cfg.regressparams.model 'rt'};
    ft_contourfTFR(cfg,regressOutput);

    set(0,'DefaultTextInterpreter','none') % disable TEX
    set(gcf,'units','normalized','outerposition',[0 0 1 1],'PaperPositionMode','auto','DefaultTextInterpreter','none','PaperOrientation','portrait'); % maximize figure
    outPath = fullfile(parentDir,'Figures','TFR_pow_stimulus_correctTrials_beta'); mkdir(outPath);
    tempFigname = fullfile(outPath,[subject '_FCz_rt.jpg']);
    print(gcf,'-djpeg','-r100', tempFigname);
    close all
    
end

%% single-trial regression: power ~ rt + congruency

if false
    % create design matrix
    epochs_allVars;
    regressors = {'rt' 'congruencyDC'};
    colIdx = ismember(epochs_allVars,regressors);
    designMat = epochs_clean(:,colIdx);
    
    % z-score
    designMat(:,1) = zscore(designMat(:,1));
    
    % fit model
    cfg = [];
    cfg.designmatrix = designMat;
    cfg.model = 'y ~ rt + congruency';
    cfg.intercept = 'no';
    cfg.transformY = 'tiedrank';
    cfg.type = 'zpermute';
    cfg.n_permutes = 500;
    regressOutput = ft_regressTFRPower_permute(cfg,tf_pow);
    regressOutput
    
    %% save results

    outPath = fullfile(parentDir,'Results','TFR_pow_stimulus_correctTrials_beta_rtCongruency'); mkdir(outPath);
    save(fullfile(outPath,[subject '.mat']),'regressOutput');
    
    %% plot

    % rt
    close all
    cfg = [];
    cfg.zMax = 2;
    cfg.zlim = [-cfg.zMax cfg.zMax];
    cfg.comment = [regressOutput.cfg.regressparams.model ' (rt)'];
    cfg.showlabels = 'yes';
    cfg.parameter = 'powspctrm';
    cfg.layout = 'biosemi64.lay';
    cfg.colormap = 'viridis';
    cfg.colorbar = 'yes';
    cfg.fontsize = 14;
    cfg.trials = 1; % coefficient number
    close all; figure(201); clf
    ft_multiplotTFR(cfg,regressOutput);

    set(gcf,'units','normalized','outerposition',[0 0 1 1],'PaperPositionMode','auto','DefaultTextInterpreter','none','PaperOrientation','portrait'); % maximize figure
    outPath = fullfile(parentDir,'Figures','TFR_pow_stimulus_correctTrials_beta'); mkdir(outPath);
    tempFigname = fullfile(outPath,[subject '_multiplotTF_rt_controlCongruency.jpg']);
    print(gcf,'-djpeg','-r100', tempFigname);
    close all

    % congruency
    close all
    cfg = [];
    cfg.zMax = 2;
    cfg.zlim = [-cfg.zMax cfg.zMax];
    cfg.comment = [regressOutput.cfg.regressparams.model ' (congruency)'];
    cfg.showlabels = 'yes';
    cfg.parameter = 'powspctrm';
    cfg.layout = 'biosemi64.lay';
    cfg.colormap = 'viridis';
    cfg.colorbar = 'yes';
    cfg.fontsize = 14;
    cfg.trials = 2; % coefficient number
    close all; figure(201); clf
    ft_multiplotTFR(cfg,regressOutput);

    set(gcf,'units','normalized','outerposition',[0 0 1 1],'PaperPositionMode','auto','DefaultTextInterpreter','none','PaperOrientation','portrait'); % maximize figure
    outPath = fullfile(parentDir,'Figures','TFR_pow_stimulus_correctTrials_beta'); mkdir(outPath);
    tempFigname = fullfile(outPath,[subject '_multiplotTF_congruency_controlRt.jpg']);
    print(gcf,'-djpeg','-r100', tempFigname);
    close all
    
    % rt and congruency
    close all
    cfg = [];
    cfg.subplot = [1 2 1];
    cfg.chan = 'FCz';
    cfg.trials = 1; % coefficient number
    cfg.fontsize = 14;
    cfg.colormap = 'viridis';
    cfg.title = {cfg.chan,regressOutput.cfg.regressparams.model 'rt'};
    ft_contourfTFR(cfg,regressOutput);
    
    cfg.subplot = [1 2 2];
    cfg.chan = 'FCz';
    cfg.trials = 2; % coefficient number
    cfg.fontsize = 14;
    cfg.colormap = 'viridis';
    cfg.title = {cfg.chan,regressOutput.cfg.regressparams.model 'congruency'};
    ft_contourfTFR(cfg,regressOutput);

    set(0,'DefaultTextInterpreter','none') % disable TEX
    set(gcf,'units','normalized','outerposition',[0 0 1 1],'PaperPositionMode','auto','DefaultTextInterpreter','none','PaperOrientation','portrait'); % maximize figure
    outPath = fullfile(parentDir,'Figures','TFR_pow_stimulus_correctTrials_beta'); mkdir(outPath);
    tempFigname = fullfile(outPath,[subject '_FCz_rtCongruency.jpg']);
    print(gcf,'-djpeg','-r100', tempFigname);
    close all
    
end

%% single-trial regression: power ~ CNV

if false
    
    % create design matrix
    epochs_allVars;
    regressors = {'cnvAmplitude'};
    colIdx = ismember(epochs_allVars,regressors);
    designMat = epochs_clean(:,colIdx);
    
    % z-score
    designMat(:,1) = zscore(designMat(:,1));
    
    % fit model
    cfg = [];
    cfg.designmatrix = designMat;
    cfg.model = 'y ~ cnv';
    cfg.intercept = 'no';
    cfg.transformY = 'tiedrank';
    cfg.type = 'zpermute';
    cfg.n_permutes = 500;
    regressOutput = ft_regressTFRPower_permute(cfg,tf_pow);
    regressOutput
    
    %% save results

    outPath = fullfile(parentDir,'Results','TFR_pow_stimulus_correctTrials_beta_cnv'); mkdir(outPath);
    save(fullfile(outPath,[subject '.mat']),'regressOutput');
    
    %% plot

    close all
    cfg = [];
    cfg.zMax = 2;
    cfg.zlim = [-cfg.zMax cfg.zMax];
    cfg.comment = regressOutput.cfg.regressparams.model;
    cfg.showlabels = 'yes';
    cfg.parameter = 'powspctrm';
    cfg.layout = 'biosemi64.lay';
    cfg.colormap = 'viridis';
    cfg.colorbar = 'yes';
    cfg.fontsize = 14;
    cfg.trials = 1; % coefficient number
    close all; figure(201); clf
    ft_multiplotTFR(cfg,regressOutput);

    set(gcf,'units','normalized','outerposition',[0 0 1 1],'PaperPositionMode','auto','DefaultTextInterpreter','none','PaperOrientation','portrait'); % maximize figure
    outPath = fullfile(parentDir,'Figures','TFR_pow_stimulus_correctTrials_beta'); mkdir(outPath);
    tempFigname = fullfile(outPath,[subject '_multiplotTF_cnv.jpg']);
    print(gcf,'-djpeg','-r100', tempFigname);
    close all

    cfg = [];
    cfg.chan = 'FCz';
    cfg.trials = 1; % coefficient number
    cfg.fontsize = 14;
    cfg.colormap = 'viridis';
    cfg.title = {cfg.chan,regressOutput.cfg.regressparams.model 'cnv'};
    ft_contourfTFR(cfg,regressOutput);

    set(0,'DefaultTextInterpreter','none') % disable TEX
    set(gcf,'units','normalized','outerposition',[0 0 1 1],'PaperPositionMode','auto','DefaultTextInterpreter','none','PaperOrientation','portrait'); % maximize figure
    outPath = fullfile(parentDir,'Figures','TFR_pow_stimulus_correctTrials_beta'); mkdir(outPath);
    tempFigname = fullfile(outPath,[subject '_FCz_cnv.jpg']);
    print(gcf,'-djpeg','-r100', tempFigname);
    close all
    
end







%% 

clear
close all

end % end function
