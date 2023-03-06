function [] = sortTurkValData()
    %% sortTurkValData
    
    
    % This script sorts the data from their version 1 of our IDM turk task.
    % specifically the 56 identified subjects who did well in the percpetual
    % decision-making task
    
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
    tasks       = {'crdm','cdd'};
    
    for ib = 1:numel(batches)
        batch           = batches(ib).name;
        disp(['Working on Batch: ',batch])
        for iT = 1:numel(tasks)
            folderList      = dir(fullfile(dataPath,batch));
            folderNames     = {folderList.name};
            folderNames     = folderNames(~ismember(folderNames ,{'.','..'}));
            folderInd       = 1;
            for iF = 1:numel(folderNames)
                tmpFileList     = dir(fullfile(dataPath,batch,folderNames{iF},tasks{iT},'*SV_hat.csv'));
                if ~isempty(tmpFileList)
                    subList(folderInd)    = tmpFileList;
                    folderInd       = folderInd + 1;
                end
            end
                
            for iO = 1:numel(subList)
                disp(['Processing subject: ',num2str(iO),' of ',num2str(numel(subList))])
                clear dataMat
                try
                    tmpData     = readtable(fullfile(subList(iO).folder,subList(iO).name));
                    % do not use ambig_trial column here
                    tmpData     = tmpData(:,{'SV_delta','confidence'});
                    % colNames    = tmpData.Properties.VariableNames;
                    dataRaw     = table2cell(tmpData);
                    for iTrial = 1:size(dataRaw,1)
                        dataMat(iTrial,1) = double(dataRaw{iTrial,1});
                        if ~ischar(dataRaw{iTrial,2})
                            dataMat(iTrial,2) = double(dataRaw{iTrial,2});
                        elseif ~isempty(str2num(dataRaw{iTrial,2}))
                            dataMat(iTrial,2) = str2num(dataRaw{iTrial,2});
                        else
                            dataMat(iTrial,2) = nan;
                        end
                    end
                                
                    % These lines recode the confidence responses from 4 point
                    % scale to 2 point scale
                    for iC = 1:3 %criterion
                        confVariation(iC) = mean(abs(dataMat(:,2))>iC);
                    end
                    [~,confCriterion]   = min(abs(confVariation-.5));
                    
                    highInds            = find(abs(dataMat(:,2))>confCriterion);
                    lowInds             = find(abs(dataMat(:,2))<=confCriterion);
                    dataMat(highInds,2) = 2*sign(dataMat(highInds,2));
                    dataMat(lowInds,2)  = 1*sign(dataMat(lowInds,2));
                                                
                    trials{iO}{1}{1} = dataMat;
                    
                catch e %e is an MException struct
                    fprintf(1,'The identifier was:\n%s',e.identifier);
                    fprintf(1,'There was an error! The message was:\n%s',e.message);
                    ProblemFiles{iO} = subList(iO).name;
                    keyboard
                end
            end
            saveName = fullfile(outPath,batch,[tasks{iT},'_trials_subList_' date '.mat']);
            save(saveName, 'trials','subList');
            clear trials subList;
        end
    
    end
end
