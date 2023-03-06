function [] = fitTurkDataComb()

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
    tasks       = {'cpdm'};
    
    
    % Set parameters
    nRuns         = 25;         % Number of restarts for model fit (no need to run 25 loops, try 3 next time)
    attemptsTotal = 5;        % Number of attempts at improving the fit
    sampleRate    = 100;       % Higher values produce slower, more precise estimates.
    delta         = 5;         % Number of standard deviations below and above mean, used to compute confidence variable distributions
    calcPrecision = [sampleRate, delta];
    asymFlag      = 0;         % if 1, fit CASANDRE with asymmetrical confidence criteria; if 0, fit CASANDRE with symmetrical confidence criteria
    
    % Set options
    options               = optimoptions('fmincon');
    options.MaxIterations = 25;
    options.Display       = 'off';
    


    for ib = 1:numel(batches)
        batch           = batches(ib).name;
        disp(['Working on Batch: ',batch])

        fileList    = dir(fullfile(outPath,batch,[tasks{1},'_trials_subList_*.mat']));
        % Load data
        load(fullfile(fileList.folder,fileList.name), 'trials','condLabels','bonus','subList');
    
        % Step 1: unpack data
        nSubjects      = numel(trials);
        nTasks    = numel(trials{1});
        nConfCrit = 1;              % Number of confidence criteria if asymFlag = 0;
        for iO = 1:nSubjects
            tic;
            improveAttempts = 0;
            for iL = 1:nRuns
                if bonus(iO)>0
                    % fprintf('Fitting observer %d, loop %d... \n', iO, iL)
                    fprintf('Fitting data for subject %d of %d, task %s, loop %d... \n', iO, nSubjects, tasks{1}, iL);
                    clear stimValue choiceVec
                    try
                        for iT = 1:nTasks
                            totRel      = numel(condLabels{2}{end});
                            totRelLevs  = condLabels{2}{end};
                            nRel        = numel(trials{1}{iT});
                            nRelLevls   = condLabels{2}{iT};
                            numTasks    = numel(condLabels{1});
                            
                            nParams = 2 + numTasks + totRel + numTasks*nConfCrit; % [Guess rate, stimulus criterion], [meta-uncertainty], [stimulus sensitivity], [confidence criteria]
                            
                            if asymFlag == 1
                                warning('This doesnt work yet in this program');
                                keyboard
                                nParams = nParams + numTasks*nConfCrit; %account for additional parameters if using asymmetrical confidence critiera
                            end
                            
                            for j = 1:numel(nRelLevls)
                                relInds(j) = find(nRelLevls(j)==totRelLevs);
                            end
                            
                            stimValue{iT} = cell(1,totRel);
                            choiceVec{iT} = cell(1,totRel);
                            for iR = 1:nRel
                                % Set experiment parameters
                                stimValue{iT}{relInds(iR)}  = trials{iO}{iT}{iR}(:,1);   % The different stimulus conditions in units of stimulus magnitude (e.g., orientation in degrees)
                                choice                      = trials{iO}{iT}{iR}(:,2);   % Subject choice (sign = perceptual choice; value = confidence rating)
                                
                                % Recode responses
                                respOptions                 = [-(nConfCrit+1):-1, 1:(nConfCrit+1)];
                                choiceVec{iT}{relInds(iR)}  = zeros(2*(nConfCrit+1), numel(stimValue{iT}{relInds(iR)}));
                                
                                for iC = 1:numel(respOptions)
                                    selTrials                                   = (choice == respOptions(iC));
                                    choiceVec{iT}{relInds(iR)}(iC, selTrials)   = 1;
                                end
                            end
                        end
                        
                        
                        % Initialize parameters
                        saveName = fullfile(outPath,batch,[tasks{1},'_bestParams_subjects_' date '.mat']);
                        try load(saveName, 'bestParamEstComb');
                            startVec = bestParamEstComb{iO}{iT};
                        catch
                            startVec = [.01 abs(randn(1,totRel)) -0.1 abs(randn(1,numTasks)) sort(2*rand(1,numTasks*nConfCrit))];
                            if asymFlag == 1
                                startVec = [startVec, rand(1, numTasks*nConfCrit)];
                            end
                        end
                        
                        % Search bounds:
                        LB          = zeros(numel(startVec),1);
                        UB          = zeros(numel(startVec),1);
                        
                        LB(1,1)                           = 0;    UB(1,1)                           = 0.075;          % Guess rate
                        LB(2:1+totRel,1)                  = 0.01; UB(2:1+totRel,1)                  = 10;            % Stimulus sensitivity
                        LB(2+totRel,1)                    = -3;   UB(2+totRel,1)                    = 3;             % Stimulus criterion
                        LB(3+totRel:3+totRel+numTasks,1)  = 0.1;  UB(3+totRel:3+totRel+numTasks,1)  = 3;             % Meta uncertainty
                        LB(end-numTasks+1:end,1)          = 0.01;   UB(end-numTasks+1:end,1)        = 5.1;             % Confidence criteria
                        
                        % Define objective function
                        obFun = @(paramVec) giveNLL(paramVec, stimValue, choiceVec, calcPrecision, asymFlag);
                        
                        % Introduce some noise
                        startVec = startVec .* .9+.5*rand(size(startVec));
                        
                        % Make sure start is within bounds
                        startVec = max(min(startVec, UB' - 0.001), LB' + 0.001);
                        
                        % Fit model
                        [paramEst{iO}, NLL{iO}] = fmincon(obFun, startVec, [], [], [], [], LB, UB, [], options);
                        
                        %             % Fit model without meta-uncertainty
                        %             LB(2,1) = 0.001;    UB(2,1) = 0.002;   startVec(2) = 0.0015;
                        %             [paramEst_NM{iO}{iT}, NLL_NM{iO}{iT}] = fmincon(obFun, startVec, [], [], [], [], LB, UB, [], options);
                        
                        % Compare current fits with best fits
                        try load(saveName, 'bestParamEstComb', 'bestNLLComb');
                            if bestNLLComb{iO} <= NLL{iO}
                                if improveAttempts > attemptsTotal
                                    fprintf('We tried improving for %d loops without success. \n',improveAttempts);
                                    fprintf('**WARNING** We will save and move on to the next subject \n');
                                    break
                                end
                                improveAttempts = improveAttempts+1;
                            else
                                bestNLLComb{iO}      = NLL{iO};
                                bestParamEstComb{iO} = paramEst{iO};
                                fprintf('New best fit found for full model! \n')
                                % Save fits to new .mat file
                                save(saveName, 'trials', 'bestParamEstComb', 'bestNLLComb','condLabels','bonus','subList');
                                % Found improvement so let's reset
                                improveAttempts = 0;
                            end
                        catch
                            bestNLLComb{iO}      = NLL{iO};
                            bestParamEstComb{iO} = paramEst{iO};
                        end
                        
                        % Save fits to new .mat file
                        save(saveName, 'trials', 'bestParamEstComb', 'bestNLLComb','condLabels','bonus','subList');
                        fprintf('Success fitting subject %d, loop %d... \n', iO, iL)
                    catch
                        warning('Unsuitable subject');
                    end
                end
            end
            tElapsed = toc/60.0;
            fprintf('Finished, for subject %d of %d, it took %0.2f minutes \n',iO,nSubjects,tElapsed)
        end
        
    end
end


function [NLL] = giveNLL(paramVec, stimValue, choiceVec, calcPrecision, asymFlag)
% Unpack data
numTasks    = numel(stimValue);
nRel        = numel(stimValue{1});  % n levels of total stimulus uncertainty
NLLvec      = nan*zeros(numTasks,nRel);

for iT  = 1:numTasks
    for iR = 1:nRel
        if ~isempty(stimValue{iT}{iR})
            % First select appropriate subset of parameters, then get NLL for each level of stimulus uncertainty
            % Required order for getLlhChoice: [guess rate, stim sens, stim crit, meta-uncertainty, conf criteria]
            params      = paramVec([1, 2+(iR-1), 2+nRel, 3+nRel+(iT-1), 2+nRel+numTasks+iT]);
            choiceLlh   = getLlhChoice(stimValue{iT}{iR}, params, calcPrecision, asymFlag);
            
            NLLvec(iT,iR)  = -sum(sum(choiceVec{iT}{iR}.*log(choiceLlh)));
        end
    end
end
% Compute NLL across entire data-set
NLL = nansum(NLLvec(:));
end
