function [] = fitTurkValTasks()

    % Start with clean slate
    clearvars;
    clc;
    
    % Set Paths
    thisPath    = fullfile(fileparts(mfilename('fullpath')));
    addpath(genpath(fullfile(thisPath,'..')));
    % Load data
    cd(thisPath)
    
    rootPath    ='/Users/pizarror/Library/CloudStorage/OneDrive-SharedLibraries-NationalInstitutesofHealth/NIMH CDN lab - Documents/Experiments/1-Data/IDM/IDM-CloudResearch/';
    % dataPath    =fullfile(rootPath,"IDM_Model_Output/");
    outPath     =fullfile(rootPath,"Confidence_analysis/CASENDRE_output/");
    batches     = dir(fullfile(outPath,'Batch*'));
    tasks       = {'crdm','cdd'};
    
    % Set parameters
    nRuns         = 25;        % Number of restarts for model fit, let's make it bigger
    attemptsTotal = 10;        % Number of attempts at improving the fit
    sampleRate    = 100;       % Higher values produce slower, more precise estimates.
    delta         = 5;         % Number of standard deviations below and above mean, used to compute confidence variable distributions
    calcPrecision = [sampleRate, delta];
    asymFlag      = 0;         % if 1, fit CASANDRE with asymmetrical confidence criteria; if 0, fit CASANDRE with symmetrical confidence criteria
    
    % Set options
    options               = optimoptions('fmincon');
    options.MaxIterations = 25;
    options.Display       = 'off';
    
    
    for ib = 1:numel(batches) % batches 1-4
        disp(['Working on Batch: ',batches(ib).name])
        batch           = batches(ib).name;
        for iT = 1:numel(tasks) % crdm, cdd
            fileList    = dir(fullfile(outPath,batch,[tasks{iT},'_trials_subList_*.mat'])); 
            % Load data
            load(fullfile(fileList.folder,fileList.name), 'trials','subList');%,'condLabels','bonus');
    
            % Step 1: unpack data
            nSubjects      = numel(trials);
            nConfCrit = 1;              % Number of confidence criteria if asymFlag = 0;
            tic;
            for iO = 1:nSubjects
                improveAttempts = 0;
                for iL = 1:nRuns
                    % taskID = iT;
                    nRel      = numel(trials{1}{1});
                    % subject = strrep(subList(iO).name,sprintf('_%s_SV_hat.csv',tasks{iT}),'');
                    fprintf('Fitting data for subject %d of %d, task %s, loop %d... \n', iO, nSubjects, tasks{iT}, iL);
                    clear stimValue choiceVec
        
                    nParams = 3 + nRel + nConfCrit; % [Guess rate, meta-noise, stimulus criterion], [stimulus sensitivity], [confidence criteria]
                    
                    if asymFlag == 1
                       nParams = nParams + nConfCrit; %account for additional parameters if using asymmetrical confidence critiera 
                    end
        
                    for iR = 1:nRel                
                        % Set experiment parameters
                        stimValue{iR} = trials{iO}{1}{iR}(:,1);   % The different stimulus conditions in units of stimulus magnitude (e.g., orientation in degrees)
                        choice        = trials{iO}{1}{iR}(:,2);   % Subject choice (sign = perceptual choice; value = confidence rating)
                        
                        % Recode responses
                        respOptions   = [-(nConfCrit+1):-1, 1:(nConfCrit+1)];
                        choiceVec{iR} = zeros(2*(nConfCrit+1), numel(stimValue{iR}));
                        
                        for iC = 1:numel(respOptions)
                            selTrials = (choice == respOptions(iC));
                            choiceVec{iR}(iC, selTrials) = 1;
                        end
                    end
                    
                    % Initialize parameters
                    saveName = fullfile(outPath,batch,[tasks{iT},'_bestParams_subjects_' date '.mat']);
                    try load(saveName, 'bestParamEst');
                        startVec = bestParamEst{iO}{1};
                    catch
                        startVec = [.01 abs(randn(1,nRel)) -0.1 0.5 sort(2*rand(1,nConfCrit))];
                        if asymFlag == 1
                            startVec = [startVec, rand(1, nConfCrit)];
                        end
                    end
        
                    % Search bounds:
                    LB          = zeros(numel(startVec),1);
                    UB          = zeros(numel(startVec),1);
                    
                    LB(1,1)             = 0;            UB(1,1)                 = 0.05;                  % Guess rate
                    LB(2:1+nRel,1)      = 0.01;         UB(2:1+nRel,1)          = 10;                   % Stimulus sensitivity
                    LB(2+nRel,1)        = -3;           UB(2+nRel,1)            = 3;                    % Stimulus criterion
                    LB(3+nRel,1)        = 0.1;          UB(3+nRel,1)            = 5;                    % Meta uncertainty
                    LB(4+nRel:end,1)    = 0;            UB(4+nRel:end,1)        = 5;                    % Confidence criteria
                    
                    % Define objective function
                    obFun = @(paramVec) giveNLL(paramVec, stimValue, choiceVec, calcPrecision, asymFlag);
                    % Introduce some noise
                    startVec = startVec .* .9+.5*rand(size(startVec));
                    % Make sure start is within bounds
                    startVec = max(min(startVec, UB' - 0.001), LB' + 0.001);
        
                    % Fit model
                    [paramEst{iO}{1}, NLL{iO}{1}] = fmincon(obFun, startVec, [], [], [], [], LB, UB, [], options);
                                
                    % Compare current fits with best fits
                    try load(saveName, 'bestParamEst', 'bestNLL');
                        if bestNLL{iO}{1} <= NLL{iO}{1}
                            % check if we are not improving after attemptsTotal
                            if improveAttempts > attemptsTotal
                                fprintf('We tried improving for %d loops without success. \n',improveAttempts);
                                fprintf('**WARNING** We will save and move on to the next subject \n');
                                break
                            end
                            improveAttempts = improveAttempts+1;
                        else
                            bestNLL{iO}{1}      = NLL{iO}{1};
                            bestParamEst{iO}{1} = paramEst{iO}{1};
                            fprintf('New best fit found for full model! \n')
                            % Save fits to new .mat file
                            save(saveName, 'trials', 'bestParamEst', 'bestNLL','subList');
                            % Found improvement so let's reset
                            improveAttempts = 0;
                        end
                    catch
                        bestNLL{iO}{1}      = NLL{iO}{1};
                        bestParamEst{iO}{1} = paramEst{iO}{1};
                        % Save fits to new .mat file
                        save(saveName, 'trials', 'bestParamEst', 'bestNLL','subList');
                    end
                    
                end
            end
            tElapsed = toc/60.0;
            fprintf('Finished, for all %d subjects, it took %0.2f minutes',nSubjects,tElapsed)
        end
    end
end

function [NLL] = giveNLL(paramVec, stimValue, choiceVec, calcPrecision, asymFlag)
    % Unpack data
    nRel = numel(stimValue);  % n levels of stimulus uncertainty
    
    for iR = 1:nRel
        % First select appropriate subset of parameters, then get NLL for each level of stimulus uncertainty
        % Required order for getLlhChoice: [guess rate, stim sens, stim crit, meta-uncertainty, conf criteria]
        params      = paramVec([1, 2+(iR-1), 2+nRel, 3+nRel, 4+nRel:end]);
        choiceLlh   = getLlhChoice(stimValue{iR}, params, calcPrecision, asymFlag);
        
        NLLvec(iR)  = -sum(sum(choiceVec{iR}.*log(choiceLlh)));
    end
    
    % Compute NLL across entire data-set
    NLL = sum(NLLvec);
end

