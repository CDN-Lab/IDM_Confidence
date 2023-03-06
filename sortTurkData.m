function [] = sortTurkData()
    %% sortTurkData
    
    % This script sorts the data from their Experiment 1. This sorting script
    % provides a general template for organizing data into the format used for
    % fitting CASANDRE.
    
    % Citation: WT Adler, WJ Ma. (2018). Comparing Bayesian and non-Bayesian
    % accounts of human confidence reports. PLOS Computational Biology. 14(11):
    % e1006572. https://doi.org/10.1371/journal. pcbi.1006572
    
    % Stimulus and task:
    % During each session, each subject completed two orientation
    % categorization tasks, Tasks A and B. The stimuli were drifting Gabors for
    % some subjects and ellipses for other subjects. On each trial, a category
    % C was selected randomly (both categories were equally probable), and a
    % stimulus s was drawn from the corresponding stimulus distribution and
    % displayed. The subject categorized the stimulus and simultaneously
    % reported their confidence on a 4-point scale, with a single button press.
    % The categories were defined by normal distributions on orientation, which
    % differed by task. In Task A, the distributions had different means (±mu_C)
    % and the same standard deviation (sigma_C); leftward-tilting stimuli were more
    % likely to be from category 1. In Task B, the distributions had the same
    % mean (0°) and different standard deviations (sigma_1, sigma_2); stimuli around
    % the horizontal were more likely to be from category 1.
    
    % Please see these links for more information about these data:
    % https://github.com/wtadler/confidence/
    % https://osf.io/s46pr/
    
    % Prelim Subjects:
    % turkPrelim
    % 1-all_data
    
    % Start with clean slate
    clearvars
    
    % Set Paths
    thisPath    = fullfile(fileparts(mfilename('fullpath')));
    addpath(genpath(fullfile(thisPath,'..')));
    % Load data
    cd(thisPath)

    rootPath    ='/Users/pizarror/Library/CloudStorage/OneDrive-SharedLibraries-NationalInstitutesofHealth/NIMH CDN lab - Documents/Experiments/1-Data/IDM/IDM-CloudResearch/';
    dataPath    =fullfile(rootPath,"IDM_Model_Output/");
    outPath     =fullfile(rootPath,"Confidence_analysis/CASENDRE_output/");
    batches     = dir(fullfile(dataPath,'Batch*'));
    tasks       = {'cpdm'};
    
    for ib = 1:numel(batches)
        batch           = batches(ib).name;
        disp(['Working on Batch: ',batch])

        folderList      = dir(fullfile(dataPath,batch));
        folderNames     = {folderList.name};
        folderNames     = folderNames(~ismember(folderNames ,{'.','..'}));
        folderInd       = 1;
        for iF = 1:numel(folderNames)
            tmpFileList     = dir(fullfile(dataPath,batch,folderNames{iF},tasks{1},'*.csv'));
            if ~isempty(tmpFileList)
                subList(folderInd)    = tmpFileList(end);
                folderInd       = folderInd + 1;
            end
        end
        % % Prelim data:
        % relColNames     = {'gabor_orient', 'gabor_contrast', 'run_dimension', 'cpdm_trial_resp_keys','cpdm_trial_resp_rt','cpdm_trials_thisTrialN'};
        condLabels{1}   = {'low_vol_low_risk', 'low_vol_high_risk', 'high_vol_low_risk', 'high_vol_high_risk'};

        for iO = 1:numel(subList)
            disp(['Processing subject: ',num2str(iO),' of ',num2str(numel(subList))])
            % tmpData     = readtable(subList(iO).name);
            tmpData     = readtable(fullfile(subList(iO).folder,subList(iO).name),'PreserveVariableNames',true);
            % tmpData     = readtable(fullfile(subList(iO).folder,subList(iO).name));
            colNames    = tmpData.Properties.VariableNames;
            % disp(colNames(:))
            
            try
                colInd = match_index_col(colNames);
                dataRaw     = table2cell(tmpData(2:end,colInd));
                for iBlock  = 1:4
                    % Step 1: Pull out first block condition of trials
                    blockEnd    = find(strcmp(condLabels{1}{iBlock},dataRaw(:,3)));
                    if numel(blockEnd) == 1 % Early prelim data
                        trialInds   = blockEnd-200:blockEnd-1;
                    elseif numel(blockEnd) == 200 % Cloud Research
                        trialInds   = blockEnd;%(1:end-1);
                    else
                        warning('Block trials not aligned');
                    end
                    
                    % Step 2: Re-code responses numerically
                    tmpResp    = strrep(dataRaw(trialInds,4),'q','-2');
                    tmpResp    = strrep(tmpResp,'p','2');
                    tmpResp    = strrep(tmpResp,'a','-1');
                    tmpResp    = strrep(tmpResp,'l','1');
                    tmpResp    = str2double(tmpResp);
                    
                    % Step 3: Reorganize reliability/contrast levels
                    dataMat             = nan*ones(numel(trialInds),4);
                    dataMat(:,[1 2 4])  = cell2mat(dataRaw(trialInds,[1 2 5]));
                    dataMat(:,3)        = tmpResp;
                    
                    condLabels{2}{iBlock}   = unique(dataMat(:,2));
                    
                    % Step 4: convert to one 2-D matrix per observer and reliability level(trial x [stimulus strength, 2-D choice])
                    relList = unique(dataMat(:,2)); % n levels of stimulus reliability
                    
                    for iR = 1:numel(relList)
                        trialList   = find(dataMat(:,2) == relList(iR));
                        
                        trials{iO}{iBlock}{iR}(:,1) = dataMat(trialList, 1);
                        trials{iO}{iBlock}{iR}(:,2) = dataMat(trialList, 3);
                        trials{iO}{iBlock}{iR}(:,3) = dataMat(trialList, 4);
                    end
                end
            catch e %e is an MException struct
                fprintf(1,'The identifier was: %s \n',e.identifier);
                fprintf(1,'There was an error! The message was: %s \n',e.message);
                ProblemFiles{iO} = subList(iO).name;
                keyboard
            end
        end
        
        
        %% Quality Check
        
        complInds = find(~cellfun(@isempty,trials));
        
        totalCatchCorrect       = nan*zeros(1,numel(trials));
        totalCatchTrials        = nan*zeros(1,numel(trials));
        totalTrials             = nan*zeros(1,numel(trials));
        totalResponded          = nan*zeros(1,numel(trials));
        
        for iO = 1:numel(trials)
            fileName(iO)    = string([subList(iO).name]);
            if ~isempty(trials{iO})
                tmpResponded    = 0;
                tmpTotal        = 0;
                tmpCatchCorrect = 0;
                tmpCatchTotal   = 0;
                
                for iB = 1:4 % number of blocks
                    if numel(trials{iO}{iB}) == 2
                        catchTrialInds      = find(abs(trials{iO}{iB}{end}(:,1))>4);
                        tmpCatchCorrect     = tmpCatchCorrect + sum(sign(trials{iO}{iB}{end}(catchTrialInds,1)) == sign(trials{iO}{iB}{end}(catchTrialInds,2)));
                        tmpCatchTotal       = tmpCatchTotal + numel(catchTrialInds);
                    elseif numel(trials{iO}{iB}) == 6
                        catchTrialInds      = find(abs(trials{iO}{iB}{end}(:,1))>4);
                        tmpCatchCorrect     = tmpCatchCorrect + sum(sign(trials{iO}{iB}{end}(catchTrialInds,1)) == sign(trials{iO}{iB}{end}(catchTrialInds,2)));
                        tmpCatchTotal       = tmpCatchTotal + numel(catchTrialInds);
                        
                        catchTrialInds      = find(abs(trials{iO}{iB}{end-1}(:,1))>4);
                        tmpCatchCorrect     = tmpCatchCorrect + sum(sign(trials{iO}{iB}{end-1}(catchTrialInds,1)) == sign(trials{iO}{iB}{end-1}(catchTrialInds,2)));
                        tmpCatchTotal       = tmpCatchTotal + numel(catchTrialInds);
                    end
                    
                    for iR = 1:numel(trials{iO}{iB})
                        tmpOData        = trials{iO}{iB}{iR};
                        tmpResponded    = tmpResponded + sum(~isnan(tmpOData(:,2)));
                        tmpTotal        = tmpTotal + numel(tmpOData(:,2));
                    end
                end
                
                totalCatchCorrect(iO)   = tmpCatchCorrect;
                totalCatchTrials(iO)    = tmpCatchTotal;
                
                totalTrials(iO)         = tmpTotal;
                totalResponded(iO)      = tmpResponded;
            end
        end
        
        catchCorrectProportion = totalCatchCorrect./totalCatchTrials;
        
        % whether they get the bonus
        bonus                       = zeros(size(fileName));
        bonus(totalResponded>780)   = 1;
        bonus(catchCorrectProportion>.75)     = 2;
        
        % whether I was able to process the data for mTurk data
        % processed   = ~cellfun(@isempty,trials);
        % bonusTable = table(fileName',bonus',processed',totalResponded',catchCorrectProportion');
        % saveNameTable = fullfile(dataPath,folderID,[folderID '.csv']);
        % writetable(bonusTable,saveNameTable);
            
        %% Save sorted data
        saveName = fullfile(outPath,batch,[tasks{1},'_trials_subList_' date '.mat']);
        % saveName = fullfile(dataPath,folderID,[date '.mat']);
        save(saveName, 'trials','condLabels','bonus','subList');
        clear trials condLabels bonus subList;
    end
end


function thisTrialN_col = get_thisTrialN(colNames)
    for iC = 1:numel(colNames)
        if endsWith(colNames{iC},'thisTrialN')
            thisTrialN_col = colNames{iC};
        end
    end
end


function colInd = match_index_col(colNames)
    % disp(colNames(:))
    thisTrialN_col = get_thisTrialN(colNames);
    relColNames     = {'cpdm_gabor_orient', 'cpdm_gabor_contrast', 'cpdm_run_dimension', 'cpdm_trial_resp.keys','cpdm_trial_resp.rt',thisTrialN_col};
    % Get indices of relevant column names
    colInd = zeros(numel(relColNames),1);
    for iCol = 1:numel(relColNames)
        comp_result = strcmp(relColNames{iCol},colNames);
        if ~sum(comp_result)
            fprintf('Unlucky with this subject...')
            disp(colNames(:))
                % exit
        end
        colInd(iCol)    = find(comp_result);
    end


end

