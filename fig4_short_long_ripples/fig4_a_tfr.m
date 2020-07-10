%% Ngo et al., eLife 2020: Sleep spindles mediate hippocampal-neocortical coupling during long-duration ripples
%
% Calculate time-frequency representations for NC and HIPP time-locked to
% ripples and NREM-control events and based on a tertial split of event
% duration
%
% Requirements:
% - Fieldtrip added to Matlab search path
% - line 29: "additional_functions" folder added to search path
% - line 31: specification of root path containing "EEGs", "Ripples" and
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

svName  = 'tfr_ripple_nrem_durationSplit';              %% filename for saved data


%% fundamental parameters
numType = 2;
numSplt = 3;
numPat  = 14;
numCh   = 2;
numRep  = 100;

fsample = 1000;


%% results structure
out = [];

out.def.type    = {'Ripple'; 'Ctrl'};
out.def.label  = {'NC'; 'HIPP'};

out.dim.indiv   = 'Type X Split: Trial x Channel x Freq x Time';
out.dim.grdavg  = 'Type x Split x Channel x Freq x Time';

out.param.artfctPad     = [-0.25 0.25];     %% padding applied to artifacts
out.param.stageoi       = [2 3 4];          %% sleep stages of interest, i.e. NREM.
out.param.minNREMLen    = 3;                %% minimal required length of NREM to be analyzed
out.param.timePad       = [-3.1 2.1];      %% temporal padding for segmentation

out.time = -2.0 : 0.02 : 1.0;   numTime = numel(out.time);      %% time vector
out.freq = 1 : 0.5 : 20;        numFreq = numel(out.freq);      %% frequency vector

out.param.mwCycl                        = ceil(out.freq * 0.5);
out.param.mwCycl(out.param.mwCycl < 5)  = 5;
out.param.mwCycl(1:8)                   = [2 3 3 3 3 3 4 4];

out.dof = nan(numPat,numSplt);

out.tfr.indiv   = cell(numType,numSplt);
out.tfr.grdavg  = nan(numType,numSplt,numCh,numFreq,numTime);

out.splt.ripple = nan(numPat,numSplt-1);
out.splt.ctrl   = nan(numPat,numRep,numSplt-1);


%% timekeeping
scrptSta = tic;


%% main loop across patients
for iPat = 1 : numPat
    fprintf('... analyze patient %d\n', iPat);

    %% get useful parameters
    tmplt   = load(fullfile(dirRoot,'EEGs',sprintf('pat%02d_NC_supplement.mat',iPat)));
    datalen = tmplt.datalen;
    scoring = ismember(tmplt.scoring,out.param.stageoi);
    
    clear tmplt
    

    %% gather data and artifacts
    inEEG       = cell(numCh,1);
    artfctfltr  = zeros(numCh,datalen);

    for iCh = 1 : numCh
        %% load data
        inEEG{iCh}  = load(fullfile(dirRoot,'EEGs',sprintf('pat%02d_%s.mat',iPat,out.def.label{iCh})));                 %% EEG
        inSupplmt   = load(fullfile(dirRoot,'EEGs',sprintf('pat%02d_%s_supplement.mat',iPat,out.def.label{iCh})));      %% supplement
        

        %% prepare artifacts
        %.. add padding, ensure padding is within data range
        artfct                          = inSupplmt.artifacts + round(out.param.artfctPad * fsample);
        artfct(artfct(:,1) < 1,1)       = 1;
        artfct(artfct(:,2) > datalen,2) = datalen;
        
        %.. create binary artifact filter
        artfctfltr(iCh,:) = hvn_createBnrySignal(artfct,datalen);
        
        clear artfct

        %.. remove NREM intervals shorter than minNREMLen
        cleanNREM                   = hvn_extrctBnryBouts(~artfctfltr);
        rmvIdx                      = diff(cleanNREM,1,2) < round(out.param.minNREMLen * fsample);
        rmvSamples                  = cell2mat(reshape(arrayfun(@(x,y)  x:y, cleanNREM(rmvIdx,1), cleanNREM(rmvIdx,2),'UniformOutput',0),1,sum(rmvIdx)));
        artfctfltr(iCh,rmvSamples)  = 1;
        
        clear inSupplmt artfct cleanNREM rmvIdx rmvSamples
    end

    
    %% append data and create binary vector indexing clean NREM data points
    godData         = ft_appenddata([],inEEG{:});
    godData.fsample = fsample;
    
    godfltr = all([scoring; ~artfctfltr]);

    clear inEEG scoring artfctfltr

    
    %% process ripples
    inRipple = load(fullfile(dirRoot,'Ripples',sprintf('pat%02d_HIPP_ripples.mat',iPat)));      %% ripple events
    
    %.. segmentation
    tfg     = [];
    tfg.trl = [inRipple.EvtInfo.maxTime' + round(out.param.timePad(1) * fsample),...
               inRipple.EvtInfo.maxTime' + round(out.param.timePad(2) * fsample),...
               ones(inRipple.EvtInfo.numEvt,1) * round(out.param.timePad(1) * fsample),...
               inRipple.EvtInfo.duration'];
    
    tfg.trl(tfg.trl(:,1) < 1 | tfg.trl(:,2) > datalen | inRipple.rejects,:)     = [];     %% remove trials outside data range
    tfg.trl(arrayfun(@(x,y) ~all(godfltr(x:y)),tfg.trl(:,1),tfg.trl(:,2)),:)    = [];     %% remove trials overlapping with artifacts

    rippleTrl = ft_redefinetrial(tfg,godData);

    clear inRipple

    %.. calculate tfr
    tfg             = [];
    tfg.output      = 'pow';
    tfg.method      = 'wavelet';
    tfg.taper       = 'hanning';
    tfg.toi         = out.time;
    tfg.foi         = out.freq;
    tfg.pad         = 'nextpow2';
    tfg.width       = out.param.mwCycl;
    tfg.keeptrials  = 'yes';
    tfg.polyremoval = 1;
    tmp_tfr         = ft_freqanalysis(tfg,rippleTrl);

    rippleTFR = squeeze(tmp_tfr.powspctrm);
    rippleDur = tmpTrl.trialinfo;

    clear rippleTrl tmp_tfr
    
    
    %% process control events
    ctrlTFR = cell(1,numRep);
    ctrlDur = cell(1,numRep);
    
    %.. loop across realizations
    for iRep = 1 : numRep
        inCtrl    = load(fullfile(dirRoot,'Control events',num2str(iRep),sprintf('pat%02d_HIPP_controlNREM.mat',iPat)));        
            
        %.. segmentation
        tfg     = [];
        tfg.trl = [inCtrl.nonEvtInfo.maxTime' + round(out.param.timePad(1) * fsample),...
                   inCtrl.nonEvtInfo.maxTime' + round(out.param.timePad(2) * fsample),...
                   ones(inCtrl.nonEvtInfo.numEvt,1) * round(out.param.timePad(1) * fsample),...
                   (inCtrl.nonEvtInfo(1).endTime' - inCtrl.nonEvtInfo(1).staTime') / fsample];

        tfg.trl(tfg.trl(:,1) < 1 | tfg.trl(:,2) > datalen,:)                        = [];   %% remove trials outside data range
        tfg.trl(arrayfun(@(x,y) ~all(godfltr(x:y)),tfg.trl(:,1),tfg.trl(:,2)),:)    = [];   %% remove trials overlapping with artifacts

        ctrlTrl = ft_redefinetrial(tfg,godData);
        
        clear inCtrl

        %.. calculate tfr
        tfg             = [];
        tfg.output      = 'pow';
        tfg.method      = 'wavelet';
        tfg.keeptrials  = 'yes';
        tfg.toi         = out.time;
        tfg.foi         = out.freq;
        tfg.pad         = 'nextpow2';
        tfg.width       = out.param.mwCycl;
        tfg.polyremoval = 1;
        tmp_tfr         = ft_freqanalysis(tfg,ctrlTrl);

        ctrlTFR{iRep}  = squeeze(tmp_tfr.powspctrm);
        ctrlDur{iRep}  = ctrlTrl.trialinfo;

        clear ctrlTrl tmp_tfr
            
    end     %% iRep

    clear godData
    
    
    %% balance trial number and average surrogate repetitions
    for iSplt = 1 : numSplt-1
        out.splt.ripple(iSub,iSplt) = prctile(rippleDur,100*iSplt/numSplt);
        out.splt.ctrl(iSub,:,iSplt) = cellfun(@(x) prctile(x,100*iSplt/numSplt),ctrlDur);
    end
        
    for iSplt = 1 : numSplt
        %... determine trial numbers
        idxRipple   = [];
        idxCtrl     = cell(1,numRep);
        
        switch iSplt
            case 1
                idxRipple = find(rippleDur < out.splt.ripple(iPat,1));
                for iRep = 1 : numRep
                    idxCtrl{iRep} = find(ctrlDur{iRep} < out.splt.ctrl(iPat,iRep,1));
                end
            case 2
                idxRipple = find((rippleDur >= out.splt.ripple(iPat,1)) & ...
                                 (rippleDur <= out.splt.ripple(iPat,2)));
                for iRep = 1 : numRep
                    idxCtrl{iRep} = find((ctrlDur{iRep} >= out.splt.ctrl(iPat,iRep,1)) & ...
                                         (ctrlDur{iRep} <= out.splt.ctrl(iPat,iRep,2)));
                end
            case 3
                idxRipple = find(rippleDur > out.splt.ripple(iPat,2));
                for iRep = 1 : numRep
                    idxCtrl{iRep} = find(ctrlDur{iRep} > out.splt.ctrl(iPat,iRep,2));
                end
        end
            
        numRipple   = numel(idxRipple);
        numCtrl     = cellfun(@(x) size(x,1),idxCtrl);

        %... extract minimal trial number
        minTrl              = min([numRipple numCtrl]);
        out.dof(iPat,iSplt) = minTrl;

        %... adjust trial numbers and average control sets
        if minTrl > 0
            %... ripples
            out.tfr.indiv{1,iSplt} = cat(1,out.tfr.indiv,rippleTFR(idxRipple(randperm(numRipple,minTrl)),:,:,:));
            
            %... surrogates
            tmp_ctrl_mat = nan(numRep,minTrl,numCh,numFreq,numTime);
            for iRep = 1 : numRep
                tmp_ctrl_mat(iRep,:,:,:,:) = ctrlTFR{iRep}(idxCtrl{iRep}(randperm(numCtrl(iRep),minTrl)),:,:,:);
            end

            out.tfr.indiv{2,iSplt} = cat(1,out.tfr.indiv{2,iSplt},squeeze(mean(tmp_ctrl_mat,1)));
            
            clear tmp_ctrl_mat
        end

        clear idxRipple idxCtrl
    end     %% iSplt
end         %% iSub


%% grand averages and baseline correction
for iType = 1 : numType
    for iSplt = 1 : numSplt
        out.tfr.grdavg(iType,iSplt,:,:,:) = squeeze(nanmean(out.tfr.indiv{iType,iSplt},1));
    end     %% iSplt
end         %% iType


%% statistical analysis: ripple vs. NREM-control events
out.stats.time = -1 : 0.02 : 1;     numTimeStats = numel(out.stats.time);
out.stats.freq = out.freq;

out.stats.tval = zeros(numSplt,numCh,numFreq,numTimeStats);
out.stats.pval = ones(numSplt,numCh,numFreq,numTimeStats);
out.stats.mask = zeros(numSplt,numCh,numFreq,numTimeStats,'logical');


for iCh = 1 : numCh
    for iSplt = 1 : numSplt
        %.. statistic parameters
        sfg                     = [];
        sfg.parameter           = 'powspctrm';
        sfg.latency             = [out.stats.time(1) out.stats.time(end)];
        sfg.method              = 'montecarlo';
        sfg.statistic           = 'indepsamplesT';
        sfg.correctm            = 'cluster';
        sfg.numrandomization    = 1000;
        sfg.neighbours          = [];
        sfg.clusteralpha        = 0.05;
        sfg.clusterstatistic    = 'maxsum';
        sfg.tail                = 0;
        sfg.alpha               = 0.025;

        sfg.design(1,1:(2*numTrl)) = [ones(1,numTrl) 2*ones(1,numTrl)];
        sfg.ivar = 1; % the 1st row in cfg.design contains the independent variable
        
        %.. prepare data
        numTrl  = sum(out.dof(:,iSplt),1);
        tmpData = cell(numType,numTrl);
        
        for iType = 1 : numType
            for iTrl = 1 : numTrl
                tmpData{iType,iTrl}.time                = out.time;
                tmpData{iType,iTrl}.freq                = out.freq;
                tmpData{iType,iTrl}.label               = {'dummy'};
                tmpData{iType,iTrl}.dimord              = 'chan_freq_time';
                tmpData{iType,iTrl}.powspctrm           = nan(1,numFreq,numTime);
                tmpData{iType,iTrl}.powspctrm(1,:,:)    = squeeze(out.tfr.indiv{iType,iSplt}(iTrl,iCh,:,:));
            end
        end
        
        %.. perform statistics
        tmpStats = ft_freqstatistics(sfg,tmpData{1,1:numTrl},tmpData{2,1:numTrl});

        out.stats.pval(iSplt,iCh,:,:) = tmpStats.prob;
        out.stats.tval(iSplt,iCh,:,:) = tmpStats.stat;
        out.stats.mask(iSplt,iCh,:,:) = tmpStats.mask;

        clear tmpStats tmpData
    end
end


%% save results
save(fullfile(dirSave,[filSave '.mat']),'-struct','out');


%% timekeeping
fprintf('This analysis took %.2f s\n', toc(statime));
