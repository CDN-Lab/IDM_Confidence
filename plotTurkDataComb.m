function [] = plotTurkDataComb()
    %% plotAdlerMa2018
    
    % Code to  visualize CASANDRE fits to Adler & Ma (2018)
    % Expt 1 data.
    
    % Start with clean slate
    clearvars;
    clc;
    close all;
    
    % Set Paths
    thisPath    = fullfile(fileparts(mfilename('fullpath')));
    addpath(genpath(fullfile(thisPath,'..')));

    rootPath    ='/Users/pizarror/Library/CloudStorage/OneDrive-SharedLibraries-NationalInstitutesofHealth/NIMH CDN lab - Documents/Experiments/1-Data/IDM/IDM-CloudResearch/';
    % dataPath    =fullfile(rootPath,"IDM_Model_Output/");
    outPath     =fullfile(rootPath,"Confidence_analysis/CASENDRE_output/");
    batches     = dir(fullfile(outPath,'Batch*'));
    tasks       = {'cpdm'};

    % Specify some plotting variables
    nSamples    = 100;      % Resolution of predicted psychometric and confidence functions
    nPoints     = [8, 8, 5, 5];   % Number of points in the observed psychometric and confidence functions [task 1, task 2]
    
    % Calculation precision
    sampleRate    = 100;        % Higher values produce slower, more precise estimates.
    delta         = 5;          % Number of standard deviations below and above mean, used to compute confidence variable distributions
    calcPrecision = [sampleRate, delta];
    asymFlag      = 0;
    
    for ib = 1:numel(batches) % batches 1-4
        disp(['Working on Batch: ',batches(ib).name])
        batch           = batches(ib).name;

        fileList    = dir(fullfile(outPath,batch,[tasks{1},'_bestParams_subjects_*.mat'])); 
        % Load data
        load(fullfile(fileList.folder,fileList.name),'trials', 'bestParamEstComb', 'bestNLLComb','condLabels','bonus','subList');
        nSubjects      = numel(trials);

        for iT = 1:4 %plotTask %blocks

            subject_metauncertainty = cell(nSubjects,2);
            tic;

            for iO = 1:nSubjects
                subject = strrep(subList(iO).name,sprintf('_%s.csv',tasks{1}),'');
                fprintf('Plotting data for subject %d of %d, block %d, task %s... \n', iO, nSubjects, iT, tasks{1});
    
                set(figure(iT), 'OuterPosition', [100 100 2000 1000])
                % set(figure(2), 'OuterPosition', [100 100 2000 1000])
                % iT = plotTask;
                % iO = plotObs;
                                
                % Step 1: unpack data
                % nTasks      = numel(trials{1});
                nRel        = numel(trials{1}{iT});
                totRel      = numel(condLabels{2}{end});
                totRelLevs  = condLabels{2}{end};
                nRelLevls   = condLabels{2}{iT};
                numTasks    = numel(condLabels{1});
                for j = 1:numel(nRelLevls)
                    relInds(j) = find(nRelLevls(j)==totRelLevs);
                end
                if asymFlag == 1
                    warning('This doesnt work yet in this program');
                    keyboard
                end
                
                nConfCrit = 3;
                %%
                
                for iR = 1:nRel
                    
                    % Set plot color
                    col = [1-iR/nRel 0 iR/nRel];
                    
                    % Set experiment parameters
                    stimValue{iR} = trials{iO}{iT}{iR}(:,1);   % The different stimulus conditions in units of stimulus magnitude (e.g., orientation in degrees)
                    choice        = trials{iO}{iT}{iR}(:,2);
                    
                    % Compute observed psychometric and confidence functions (this requires grouping of stimuli)
                    edges = prctile(stimValue{iR}, linspace(0, 100, nPoints(iT) + 1));
                    [nTrials, edges, binInd] = histcounts(stimValue{iR}, edges);
                    for iP = 1:nPoints(iT)
                        oriPF(iP) = nanmean(stimValue{iR}(binInd == iP));
                        obsPF(iP) = nanmean((choice(binInd == iP)) > 0);
                        obsCF(iP) = nanmean(abs(choice(binInd == iP)));
                        
                        % Psychometric for low and high confidence trials
                        indHighConf    = abs(choice) >= 2;
                        obsPF_HC(iP)   = nanmean((choice(binInd == iP & indHighConf)) > 0);
                        obsPF_LC(iP)   = nanmean((choice(binInd == iP & ~indHighConf)) > 0);
                        nTrials_HC(iP) = sum(binInd == iP & indHighConf);
                        nTrials_LC(iP) = sum(binInd == iP & ~indHighConf);
                    end
                    
                    % First specify function arguments
                    stimPlot = linspace(min(stimValue{iR}), max(stimValue{iR}), nSamples);
                    
                    % Required order for getLlhChoice: [guess rate, stim sens, stim crit, meta uncertainty, conf criteria]
                    params   = bestParamEstComb{iO}([1, 2+(relInds(iR)-1), 2+totRel, 3+totRel+(iT-1), 2+totRel+numTasks+iT]);
                    paramsNM = [params(1:3), 0.01, params(5:end)];
                    paramsNG = [0, params(2:end)];
                    
                    % Get model predictions
                    fitLlh   = getLlhChoice(stimPlot, params, calcPrecision, asymFlag);
                    fitLlhNG = getLlhChoice(stimPlot, paramsNG, calcPrecision, asymFlag);
                    
                    % The PF and CF predicted on the basis of the likelihood functions
                    predPF = sum(fitLlh(size(fitLlh, 1)/2+1:end, :));
                    predCF = fitLlh' * [2 1 1 2]'; % multiply choice-likelihood by confidence ratings
                    
                    % Predicted CF without meta-uncertainty
                    predCFnoMN  = getLlhChoice(stimPlot, paramsNM, calcPrecision, asymFlag)' * [2 1 1 2]';
                    
                    % Predicted PF without guesses, split out for low and high confidence trials
                    predPF_HC = fitLlhNG(4, :)./sum(fitLlhNG([1, 4], :));
                    predPF_LC = fitLlhNG(3, :)./sum(fitLlhNG([2, 3], :));
                    
                    subplotColumnNumber = 6;
                    
                    % Plot PF
                    figure(iT)
                    subplot(4,subplotColumnNumber,iR)
                    plot(stimPlot, predPF, '-', 'linewidth', 2, 'color', col)
                    hold on, box off, axis square
                    axis([-5 5 0 1])
                    xlabel('Stimulus value')
                    ylabel('Proportion category 1')
                    
                    for iP = 1:nPoints(iT)
                        plot(oriPF(iP), obsPF(iP), 'ko', 'markerfacecolor', col, 'markersize', round(nTrials(iP))+1)
                    end
                    
                    if iR == nRel
                        legend('Full model', 'location', 'best')
                    end
                    
                    % Plot CF
                    figure(iT)
                    subplot(4,subplotColumnNumber, subplotColumnNumber+iR)
                    plot(stimPlot, predCF, '-', 'linewidth', 2, 'color', col)
                    hold on, box off, axis square
                    plot(stimPlot, predCFnoMN, 'k--', 'linewidth', 2)
                    axis([-5 5 1 2])
                    xlabel('Stimulus value')
                    ylabel('Mean confidence level')
                    
                    for iP = 1:nPoints(iT)
                        plot(oriPF(iP), obsCF(iP), 'ko', 'markerfacecolor', col, 'markersize', round(nTrials(iP))+1)
                    end
                    
                    if iR == nRel
                        legend('Full model', 'meta uncertainty = 0', 'location', 'best')
                    end
                    
                    % Plot PF split out for high and low confidence trials
                    figure(iT)
                    subplot(4,subplotColumnNumber,2*subplotColumnNumber+iR)
                    plot(stimPlot, predPF_HC, '-', 'linewidth', 2, 'color', [0 1 0])
                    hold on, box off, axis square
                    plot(stimPlot, predPF_LC, '-', 'linewidth', 2, 'color', [1 0 0])
                    axis([-5 5 0 1])
                    xlabel('Stimulus value')
                    ylabel('Proportion category 1')
                    
                    for iP = 1:nPoints(iT)
                        plot(oriPF(iP), obsPF_HC(iP), 'ko', 'markerfacecolor', [0 1 0], 'markersize', round(nTrials_HC(iP))+1)
                        plot(oriPF(iP), obsPF_LC(iP), 'ko', 'markerfacecolor', [1 0 0], 'markersize', round(nTrials_LC(iP))+1)
                    end
                    
                    if iR == nRel
                        legend('High confidence', 'Low confidence', 'location', 'best')
                    end
                    
                    % Plot confidence vs performance
                    figure(iT)
                    subplot(4,subplotColumnNumber,3*subplotColumnNumber+1)
                    plot(predPF, predCF, '-', 'linewidth', 2, 'color', col)
                    hold on, box off, axis square
                    for iP = 1:nPoints(iT)
                        plot(obsPF(iP), obsCF(iP), 'ko', 'markerfacecolor', col, 'markersize', round(nTrials(iP))+1)
                    end
                    axis([0 1 1 2])
                    xlabel('Proportion category 1')
                    ylabel('Mean confidence level')
                    
                    % Plot param values
                    figure(iT)
                    subplot(4,subplotColumnNumber,3*subplotColumnNumber+2);

                    metauncertainty = bestParamEstComb{iO}(3+totRel+(iT-1)); %bestParamEstComb{iO}(2+totRel);
                    % subject = subList(iO).name(1:30);
                    subject_metauncertainty(iO,:) = {subject,metauncertainty};
                    
                    text(0,.5,sprintf('Lapse rate = %.2g \nDecision criterion = %.2g \nMeta-uncertainty = %.2g \nConfidence criterion = %.2g \n%s\nNLL = %.4g\n ',bestParamEstComb{iO}(1),bestParamEstComb{iO}(2+totRel),bestParamEstComb{iO}(3+totRel+(iT-1)),bestParamEstComb{iO}(2+totRel+numTasks+iT),condLabels{1}{iT},bestNLLComb{iO}));
                    %     text(0,.85,sprintf('Decision criterion = %.2g \nMeta-uncertainty = %.2g \nConfidence criterion = %.2g \nLapse rate = %.2g\n%s\n ',bestParamEstComb{iO}(2+nRel),bestParamEstComb{iO}(3+nRel),bestParamEstComb{iO}(4+nRel),bestParamEstComb{iO}(1),condLabels{1}{plotTask}));
                    % Required order for getLlhChoice: [guess rate, stim sens, stim crit, meta uncertainty, conf criteria]
                    
                    
                    interLine       = .1;
                    linePosition    = .75;
                    for itextR = 1:nRel
                        text(.85,linePosition,sprintf('Sensitivity%d = %.2g',itextR,bestParamEstComb{iO}(2+(itextR-1))));
                        linePosition = linePosition-interLine;
                    end
                    
                    interLine       = .1;
                    linePosition    = .75;
                    for itextT = 1:numTasks
                        text(1.6,linePosition,sprintf('Meta-uncertainty%d = %.2g',itextT,bestParamEstComb{iO}(3+totRel+(itextT-1))));
                        linePosition = linePosition-interLine;
                    end
                    
                    interLine       = .1;
                    linePosition    = .75;
                    for itextT = 1:numTasks
                        text(2.3,linePosition,sprintf('Confidence criterion%d = %.2g',itextT,bestParamEstComb{iO}(2+totRel+numTasks+itextT)));
                        linePosition = linePosition-interLine;
                    end
                    
                    axis([0 1 0 1]);
                    axis off
                    
                end
                
                figFolder = fullfile(outPath,batch,'figures');
                if ~exist(figFolder, 'dir')
                    mkdir(figFolder)
                end
                % have to add .eps to end of it because of the periods in subject names
                saveas(gcf,fullfile(figFolder,sprintf('%s_%s_block_%d.eps',tasks{1},subject,iT)),'epsc');
                close
                % saveas(gcf,fullfile(thisPath,'figures',sprintf('S%s_block%d',subList(plotObs).name(1:10),plotTask)),'epsc');
            end

            subj_meta_fn = fullfile(outPath,batch,sprintf('%s_block_%d_metauncertainty_by_subject.csv',tasks{1},iT));
            sub_meta_table = cell2table(subject_metauncertainty,"VariableNames",["subject" "meta_uncertainty"]);
            writetable(sub_meta_table,subj_meta_fn)
            tElapsed = toc/60.0;
            fprintf('Finished, for all %d subjects, it took %0.2f minutes \n',nSubjects,tElapsed)
        end
    end

