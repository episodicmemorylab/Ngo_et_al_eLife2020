%% Ngo et al., eLife 2020: Sleep spindles mediate hippocampal-neocortical coupling during long-duration ripples
%
% Calculate NC-HIPP partial directed coherence time-locked to ripples and
% NREM-control events
%
% Requirements:
% - Fieldtrip added to Matlab search path
% - line 30: "additional_functions" folder added to search path
% - line 32: specification of root path containing "EEGs", "Ripples" and
%   "Control events" folder (see https://osf.io/3hpvr/)
% - "mask_spindleCluster.mat" obtained from
%   'fig2_a_tfr_ripple_nremCtrl'-script accessible in root folder
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

svName  = 'partialdirectedcoh';                         %% filename for saved data


%% fundamental parameters
numType = 2;
numPat  = 14;
numCh   = 2;
numRep  = 100;

fsample = 1000;


%% results structure
out = [];

out.def.type   = {'Ripple'; 'Ctrl'};
out.def.label  = {'NC'; 'HIPP'};

out.dim.pdc = '(Rep) x Freq x Time';

out.param.artfctPad     = [-0.25 0.25];     %% padding applied to artifacts
out.param.stageoi       = [2 3 4];          %% sleep stages of interest, i.e. NREM.
out.param.minNREMLen    = 3;                %% minimal required length of NREM to be analyzed
out.param.timePad       = [-2.1 2.1];       %% temporal padding for segmentation
out.param.fftPad        = [-0.256 0.255];
out.param.do_zscore     = 1;

out.time = -1.0 : 0.02 : 1.0;       numTime = numel(out.time);          %% time vector
out.freq = 0 : 1/0.512 : 21;        numFreq = numel(out.freq);          %% frequency vector

out.dof         = [];
out.pdc.ripple  = nan(numCh,numFreq,numTime);
out.pdc.ctrl    = nan(numRep,numCh,numFreq,numTime);

out.stats.map.zval = zeros(numFreq,numTime);
out.stats.map.pval = ones(numFreq,numTime);
out.stats.map.mask = zeros(numFreq,numTime,'logical');

out.stats.spindleCluster = cell(1,numCh);


%% timekeeping
scrptSta = tic;


%% step 1: loop through all patients and repetitions to acquire minimal trial number
tmp_trl = cell(numPat,1);

for iPat = 1 : numPat
    fprintf('...step 1: patient %d\n', iPat); 

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

        %.. remove NREM intervals shorter than minNREMLen
        cleanNREM                    = hvn_extrctBnryBouts(~artfctfltr);
        rmvIdx                      = diff(cleanNREM,1,2) < round(out.param.minNREMLen * fsample);
        rmvSamples                  = cell2mat(reshape(arrayfun(@(x,y)  x:y, cleanNREM(rmvIdx,1), cleanNREM(rmvIdx,2),'UniformOutput',0),1,sum(rmvIdx)));
        artfctfltr(iCh,rmvSamples)  = 1;
        
        clear inSupplmt artfct cleanNREM rmvIdx rmvSamples
    end
    
    %% create binary vector indexing clean sleep data points
    godfltr = all([scoring; ~artfltr]);
    
    clear scoring artfctfltr
    
    
    %% prepare ripple events
    inRipple = load(fullfile(dirRoot,'Ripples',sprintf('pat%02d_HIPP_ripples.mat',iPat)));      %% ripple events
    
    %.. trial definition
    tfg     = [];
    tfg.trl = [inRipple.EvtInfo.maxTime' + round(out.param.timePad(1) * fsample),...
               inRipple.EvtInfo.maxTime' + round(out.param.timePad(2) * fsample),...
               ones(inRipple.EvtInfo.numEvt,1) * round(out.param.timePad(1) * fsample)];
    
    tfg.trl(tfg.trl(:,1) < 1 | tfg.trl(:,2) > datalen | inRipple.rejects,:)     = [];     %% remove trials outside data range
    tfg.trl(arrayfun(@(x,y) ~all(godfltr(x:y)),tfg.trl(:,1),tfg.trl(:,2)),:)    = [];     %% remove trials overlapping with artifacts

    numRipple               = size(tfg.trl,1);
    tmp_trl{iPat}.ripple    = tfg.trl;
    
    clear inRipple
   
    
    %% process control events
    tmp_trl{iPat}.ctrl  = cell(1,numRep);
    numCtrl             = nan(1,numRep);
    
    for iRep = 1 : numRep
        inCtrl = load(fullfile(dirRoot,'Control events',num2str(iRep),sprintf('pat%02d_HIPP_controlNREM.mat',iPat)));
            
        %.. trial definition
        tfg     = [];
        tfg.trl = [inCtrl.nonEvtInfo.maxTime' + round(out.param.timePad(1) * fsample),...
                   inCtrl.nonEvtInfo.maxTime' + round(out.param.timePad(2) * fsample),...
                   ones(inCtrl.nonEvtInfo.numEvt,1) * round(out.param.timePad(1) * fsample)];

        tfg.trl(tfg.trl(:,1) < 1 | tfg.trl(:,2) > datalen,:)                        = [];   %% remove trials outside data range
        tfg.trl(arrayfun(@(x,y) ~all(godfltr(x:y)),tfg.trl(:,1),tfg.trl(:,2)),:)    = [];   %% remove trials overlapping with artifacts

        tmp_trl{iPat}.ctrl{iRep}    = tfg.trl;
        numCtrl(iRep)               = size(tfg.trl,1);
        
        clear inCtrl
    end     %% iRep

    
    %% save minimal trial number
    numCtrl         = cellfun(@(x) size(x,1),ctrlTrl);
    minTrl          = min([numRipple numSurrgt]);
    out.dof{iPat}   = minTrl;
    
    
    %% balance out trials
    tmp_trl{iPat}.ripple = tmp_trl{iPat}.ripple(randperm(numRipple,minTrl),:);
    
    for iRep = 1 : numRep
        tmp_trl{iPat}.ctrl{iRep} = tmp_trl{iPat}.ctrl{iRep}(randperm(numCtrl(iRep),minTrl),:);
    end
                
end     %% iSub


%% step 2: calculate pdc for ripples
rippleFFT = cell(numPat,1);

for iPat = 1 : numPat
    fprintf('... step 2a: patient %d\n', iPat);
    
    %% gather data
    inEEG  = cell(numCh,1);
    
    for iCh = 1 : numTrgt
        inEEG{iCh}  = load(fullfile(dirRoot,'EEGs',sprintf('pat%02d_%s.mat',iPat,out.def.label{iCh})));     %% EEG
    end

    godData         = ft_appenddata([],inEEG{:});
    godData.fsample = fsample;
   
    clear inEEG

    
    %% prepare ripple trials
    tfg         = [];
    tfg.trl     = tmp_trl{iPat}.ripple;
    rippleTrl   = ft_redefinetrial(tfg,godData);
    
    clear godData
    
    
    %% (optional) z-score
    if out.param.do_zscore == 1
        tmpMat  = cell2mat(rippleTrl.trial);
        tmpMean = mean(tmpMat,2);
        tmpStd  = std(tmpMat,1,2);
        rippleTrl.trial = cellfun(@(x) (x-tmpMean) ./ tmpStd, rippleTrl.trial,'UniformOutput',0);
    end
    
    
    %% calculate short-time fft
    rippleFFT{iPat} = nan(out.dof(iPat),numTime,numCh,numFreq);

    for iTime = 1 : numTime
        tfg         = [];
        tfg.toilim  = out.time(iTime) + out.param.fftPad;
        tmp_trl2    = ft_redefinetrial(tfg,rippleTrl);

        tfg                 = [];
        tfg.output          = 'fourier';
        tfg.method          = 'mtmfft';
        tfg.foi             = out.freq;
        tfg.pad             = 'nextpow2';
        tfg.taper           = 'hanning';
        tfg.keeptrials      = 'yes';
        tfg.polyremoval     = 1;
        tmp_fft             = ft_freqanalysis(tfg,tmp_trl2);

        rippleFFT{iPat}(:,iTime,:,:) = tmp_fft.fourierspctrm;

        clear tmp_trl2 tmp_fft
    end

    clear rippleTrl
    
end


%% pdc calculation
fprintf('... step 2b: calculate pdc for ripples\n');

rippleFFT = cell2mat(rippleFFT);
for iTime = 1 : numTime
    tmpData                 = [];
    tmpData.label           = out.def.label;
    tmpData.dimord          = 'rpttap_chan_freq';
    tmpData.freq            = out.freq;
    tmpData.fourierspctrm   = squeeze(rippleFFT(:,iTime,:,:));

    tfg                 = [];
    tfg.method          = 'pdc';
    tfg.pdc.channelcmb  = out.def.label';
    tmp_pdc             = ft_connectivityanalysis(tfg,tmpData);

    out.pdc.ripple(1,iTime,:) = squeeze(tmp_pdc.pdcspctrm(1,2,:));
    out.pdc.ripple(2,iTime,:) = squeeze(tmp_pdc.pdcspctrm(2,1,:)); 

    clear tmpPDC

end     %% iTime

clear rippleFFT


%% Step 3: process control events
for iRep = 1 : numRep
    
    ctrlFFT = cell(numPat,1);
    
    for iPat = 1 : numPat
        fprintf('... step 3a (rep %d): patient %d\n', iRep, iPat);

        %% gather data
        inEEG  = cell(numCh,1);

        for iCh = 1 : numTrgt
            inEEG{iCh}  = load(fullfile(dirRoot,'EEGs',sprintf('pat%02d_%s.mat',iPat,out.def.label{iCh})));
        end

        godData         = ft_appenddata([],inEEG{:});
        godData.fsample = fsample;

        clear inEEG
        
        
        %% prepare control trials
        tfg     = [];
        tfg.trl = tmp_trl{iPat}.ctrl{iRep};
        ctrlTrl = ft_redefinetrial(tfg,godData);
        
        clear godData
        
        
        %% (optional) z-score
        if out.param.do_zscore == 1
            tmpMat  = cell2mat(ctrlTrl.trial);
            tmpMean = mean(tmpMat,2);
            tmpStd  = std(tmpMat,1,2);
            ctrlTrl.trial = cellfun(@(x) (x-tmpMean) ./ tmpStd, ctrlTrl.trial,'UniformOutput',0);
            
            clear tmpMat tmpMean tmpStd
        end
        

        %% calculate short-time fft
        ctrlFFT{iPat} = nan(out.dof(iPat),numTime,numCh,numFreq);

        for iTime = 1 : numTime
            tfg         = [];
            tfg.toilim  = out.time(iTime) + out.param.fftPad;
            tmp_trl2    = ft_redefinetrial(tfg,ctrlTrl);

            %--- Perform FFT
            tfg                = [];
            tfg.output         = 'fourier';
            tfg.method         = 'mtmfft';
            tfg.foi            = out.freq;
            tfg.pad            = 'nextpow2';
            tfg.taper          = 'hanning';
            tfg.keeptrials     = 'yes';
            tfg.polyremoval    = 1;
            tmp_fft            = ft_freqanalysis(tfg,tmp_trl2);

            ctrlFFT{iPat}(:,iTime,:,:) = tmp_fft.fourierspctrm;

            clear tmp_trl2 tmp_fft
        end

        clear ctrlTrl
        
    end     %% iPat
    
    
    %% calculate coherence
    fprintf('... step 3b (rep %d): calculate pdc for control events\n', iRep);

    ctrlFFT = cell2mat(ctrlFFT);
    for iTime = 1 : numTime
        tmpData                 = [];
        tmpData.label           = out.def.label;
        tmpData.dimord          = 'rpttap_chan_freq';
        tmpData.freq            = out.freq;
        tmpData.fourierspctrm   = squeeze(ctrlFFT(:,iTime,:,:));

        tfg                 = [];
        tfg.method          = 'pdc';
        tfg.pdc.channelcmb  = out.def.label';
        tmp_pdc             = ft_connectivityanalysis(tfg,tmpData);

        out.pdc.ctrl(iRep,1,iTime,:) = squeeze(tmp_pdc.pdcspctrm(1,2,:));
        out.pdc.ctrl(iRep,2,iTime,:) = squeeze(tmp_pdc.pdcspctrm(2,1,:)); 

        clear tmpPDC

    end     %% iTime

    clear ctrlFFT
    
end     %% iRep


%% perform statistics
%% time-resolved coherence map
tmpRipple       = squeeze(out.pdc.ripple(1,:,:)-out.pdc.ripple(2,:,:))';
tmpCtrl_mean    = squeeze(mean(out.pdc.ctrl(:,1,:,:)-out.pdc.ctrl(:,2,:,:),1))';
tmpCtrl_std     = squeeze(std(out.pdc.ctrl(:,1,:,:)-out.pdc.ctrl(:,2,:,:),1,1))';

out.stats.map.zval = (tmpRipple - tmpCtrl_mean) ./ tmpCtrl_std;
out.stats.map.pval = 2*normcdf(-abs(out.stats.map.zval));
out.stats.map.mask = outStats.map.pval < 0.05;


%% significant spindle cluster
inMask = load(fullfile(dirRoot,'mask_spindleCluster.mat'));

idxTime = arrayfun(@(x) nearest(out.time,x),[inMask.time(1) inMask.time(end)]);

pdc_freq    = ceil(out.freq);                    %% round up tpdc freq
intrsctFreq = intersect(pdc_freq,inMask.freq);   %% intersection frequencies between mask and tpdc

%... remove uncommon frequency components
pdc_ripple  = out.pdc.ripple(:,idxTime(1):idxTime(2),ismember(pdc_freq,intrsctFreq));
pdc_ctrl    = out.pdc.ctrl(:,:,idxTime(1):idxTime(2),ismember(pdc_freq,intrsctFreq));

tmp_mask = nan(numel(inMask.time),numel(intrsctFreq));
tmp_mask(logical(inMask.spindle(ismember(inMask.freq,intrsctFreq),:))') = 1;
tmp_mask = permute(repmat(tmp_mask,[1 1 numRep]),[3 1 2]);

for iCh = 1 : numCh
    tmpRipple = squeeze(pdc_ripple(iCh,:,:)) .* squeeze(tmp_mask(1,:,:));
    tmpCtrl = squeeze(pdc_ctrl(:,iCh,:,:)) .* tmp_mask;

    out.stats.spindleCluster{iCh}.ripple    = mean(tmpRipple,'all','omitnan');
    out.stats.spindleCluster{iCh}.ctrl.mean = mean(tmpCtrl,'all','omitnan');
    out.stats.spindleCluster{iCh}.ctrl.dist = mean(tmpCtrl,[2 3],'omitnan');
    out.stats.spindleCluster{iCh}.ctrl.std  = std(mean(tmpCtrl,[2 3],'omitnan'),1,1,'omitnan');
    out.stats.spindleCluster{iCh}.zval      = (out.stats.spindleCluster{iCh}.ripple - out.stats.spindleCluster{iCh}.ctrl.mean) ./ out.stats.spindleCluster{iCh}.ctrl.std;
end


%% save results
save(fullfile(dirRoot,[svName '.mat']),'-struct','out');


%% timekeeping
fprintf('This analysis took %.2f s\n', toc(scrptSta));
