%% Ngo et al., eLife 2020: Sleep spindles mediate hippocampal-neocortical coupling during long-duration ripples
%
% Calculate histogram of NC- and HIPP-spindle onsets time-locked to ripples
% and NREM-control events
%
% Requirements:
% - Fieldtrip added to Matlab search path
% - line 28: "additional_functions" folder added to search path
% - line 30: specification of root path containing "EEGs", "Ripples" and
%   "Control events" folder (see https://osf.io/3hpvr/)
%
% created by H.-V.V. Ngo

clear
close all


%% bookkeeping
%% check if fieldtrip is available
if contains(path,'fieldtrip')
    ft_defaults;
else
    error('FieldTrip not in path');
end


%% directories etc.
addpath('C:\Ngo_et_al_eLife2020\additional_functions')  %% specify path

dirRoot = 'C:\Ngo_et_al_eLife2020';                     %% specify path

svName  = 'spindleOnsetHist';                           %% filename for saved data


%% fundamental variables
numType = 2;
numPat  = 14;
numRep  = 100;
numCh   = 2;

fsample = 1000;


%% results structure
out = [];

out.def.type    = {'Ripple'; 'Ctrl'};
out.def.label   = {'NC'; 'HIPP'};

out.dim.spiOnsetHist  = 'Channel x Time';

out.param.artfctPad     = [-0.25 0.25];     %% padding applied to artifacts
out.param.stageoi       = [2 3 4];          %% sleep stages of interest, i.e. NREM.
out.param.minNREMLen    = 3;                %% minimal required length of NREM to be analyzed
out.param.timePad       = [-0.5 0.5];       %% temporal padding for segmentation
out.param.binSize       = 0.05;             %% bin size of histogram (in s)
out.param.smooth        = 3;                %% smoothing factor for histogram = 3rd moving average

out.dof         = nan(numPat,numRep,numCh);
out.time        = out.param.timePad(1) : out.param.binSize : out.param.timePad(2)-out.param.binSize;
out.histogram   = zeros(numRep,numType,numCh,numel(out.time));


%% timekeeping
scrptSta = tic;


%% loop across all repetitions
for iRep = 1 : numRep
    %% Loop across all patients
    for iPat = 1 : numPat
        fprintf('analyze patient %d (Rep %d)\n', iPat, iRep);

        %% get useful parameters
        tmplt   = load(fullfile(dirRoot,'EEGs',sprintf('pat%02d_NC_supplement.mat',iPat)));
        datalen = tmplt.datalen;
        scoring = ismember(tmplt.scoring,out.param.stageoi);

        clear tmplt
        

        %% gather artifacts
        artfctfltr  = zeros(numCh,datalen);

        for iCh = 1 : numCh
            %% load data
            inSupplmt   = load(fullfile(dirRoot,'EEGs',sprintf('pat%02d_%s_supplement.mat',iPat,out.def.label{iCh})));


            %% prepare artifacts
            %.. add padding, ensure padding is within data range
            artfct                          = inSupplmt.artifacts + round(out.param.artfctPad * fsample);
            artfct(artfct(:,1) < 1,1)       = 1;
            artfct(artfct(:,2) > datalen,2) = datalen;

            %.. create binary artifact filter
            artfctfltr(iCh,:) = hvn_createBnrySignal(artfct,datalen);

            clear artfct

            %.. remove NREM intervals shorter than minSlpLen
            cleanNREM                    = hvn_extrctBnryBouts(~artfctfltr);
            rmvIdx                      = diff(cleanNREM,1,2) < round(out.param.minNREMLen * fsample);
            rmvSamples                  = cell2mat(reshape(arrayfun(@(x,y)  x:y, cleanNREM(rmvIdx,1), cleanNREM(rmvIdx,2),'UniformOutput',0),1,sum(rmvIdx)));
            artfctfltr(iCh,rmvSamples)  = 1;

            clear inSupplmt artfct cleanNREM rmvIdx rmvSamples
        end
    

        %% load ripple and surrogate events
        evtTms      = cell(numType,1);
        inRipple    = load(fullfile(dirRoot,'Ripples',sprintf('pat%02d_HIPP_ripples.mat',iPat)));                           %% ripples
        inCtrl      = load(fullfile(dirRoot,'Control events',num2str(iRep),sprintf('pat%02d_HIPP_controlNREM.mat',iPat)));  %% control events
        
        evtTms{1} = inRipple.EvtInfo(1).maxTime(~inRipple.rejects)';
        evtTms{2} = inCtrl.nonEvtInfo(1).maxTime';

        clear inRipple inCtrl


        for iCh = 1 : numCh
            %% Select artifacts of current channel
            tmpfltr = all([scoring; ~artfctfltr(iCh,:)]);

            %% Prepare temporary config
            tmp_trl = cell(numType,1);
            for iType = 1 : numType
                tmp_trl{iType} = [];
                tmp_trl{iType} = [evtTms{iType} + round(out.param.timePad(1) * fsample),...
                                  evtTms{iType} + round(out.param.timePad(2) * fsample),...
                                  ones(numel(evtTms{iType}),1) * round(out.param.timePad(1) * fsample), evtTms{iType}];

                tmp_trl{iType}(tmp_trl{iType}(:,1) < 1 | tmp_trl{iType}(:,2) > datalen,:)                       = [];   %% remove trials outside data range
                tmp_trl{iType}(arrayfun(@(x,y) ~all(tmpfltr(x:y)),tmp_trl{iType}(:,1),tmp_trl{iType}(:,2)),:)   = [];   %% remove trials overlapping with artifacts
            end

            
            %% Adjust trial numbers between ripples and surrogates
            numType1 = size(tmp_trl{1},1);
            numType2 = size(tmp_trl{2},1);

            if numType1 > numType2
                tmp_trl{1} = tmp_trl{1}(randperm(numType1,numType2),:);
            elseif numType2 > numType1
                tmp_trl{2} = tmp_trl{2}(randperm(numType2,numType1),:);
            end

            out.dof(iPat,iRep,iCh) = size(tmp_trl{1},1);              %% save degrees of freedom


            %% compute histogram
            for iType = 1 : numType
                inSpindle   = load(fullfile(dirRoot,'Spindles',sprintf('pat%02d_%s_spindles.mat',iPat,out.def.label{iCh})));    %% Spindle events
                spiOnset    = {inSpindle.EvtInfo.staTime};

                tmp_histo = hvn_slythm_evtHistogram(spiOnset,out.def.label(iCh),tmp_trl{iType}(:,4),datalen,fsample,out.param.timePad,out.param.binSize,0,out.param.smooth);

                out.hist(iRep,iType,iCh,:)  = squeeze(out.spiOnsetHist.hist(iRep,iType,iCh,:))' + tmp_histo.avg(1,:);
                
                clear inSpindle tmp_histo
            end

            clear tmpfltr tmp_trl

        end
    end         %% iPat
end             %% iRep


%% normalise data on the number of captured events
tmpHisto = 100 * out.histo ./ repmat(sum(out.histo,4),[1 1 1 numel(out.time)]);


%% perform statistics
out.stats.ripple    = squeeze(mean(tmpHisto(:,1,:,:),1,'omitnan'));
out.stats.ctrl.mean = squeeze(mean(tmpHisto(:,2,:,:),1,'omitnan'));
out.stats.ctrl.std  = squeeze(std(tmpHisto(:,2,:,:),1,'omitnan'));
out.stats.ctrl.sem  = squeeze(std(tmpHisto(:,2,:,:),1,'omitnan')) ./ sqrt(numRep);
out.stats.zval      = (out.stats.ripple - out.stats.ctrl.mean) ./ out.stats.ctrl.std;
out.stats.pval      = 2*normcdf(-abs(out.stats.zval));


%% save results
save(fullfile(dirSave,[filSave '.mat']),'-struct','out');


%% timekeeping
fprintf('This analysis took %.2f s\n', toc(scrptSta));
