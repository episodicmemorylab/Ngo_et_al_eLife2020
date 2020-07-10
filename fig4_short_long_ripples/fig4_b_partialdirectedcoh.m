%% Ngo et al., eLife 2020: Sleep spindles mediate hippocampal-neocortical coupling during long-duration ripples
%
% Calculate NC-HIPP partial directed coherence time-locked to
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

svName  = 'partialdirectedcoh_durationSplit';           %% filename for saved data


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
out.def.label   = {'NC'; 'HIPP'};

out.dim.pdc = '(Rep) x Freq x Time';

out.param.artfctPad     = [-0.25 0.25];     %% padding applied to artifacts
out.param.stageoi       = [2 3 4];          %% sleep stages of interest, i.e. NREM.
out.param.minNREMLen    = 3;                %% minimal required length of NREM to be analyzed
out.param.timePad       = [-2.1 2.1];      %% temporal padding for segmentation
out.param.fftPad        = [-0.256 0.255];
out.param.do_zscore     = 1;

out.time = -1.0 : 0.02 : 1.0;    numTime = numel(out.time);
out.freq = 0 : 1/0.512 : 21;     numFreq = numel(out.freq);

out.dof = nan(numPat,numSplt);

out.pdc.ripple = nan(numSplt,numCh,numTime,numFreq);
out.pdc.surrgt = nan(numRep,numSplt,numCh,numTime,numFreq);

out.splt.ripple     = nan(numPat,numSplt-1);
out.splt.ctrl       = nan(numPat,numRep,numSplt-1);

out.stats.spindleBand = cell(1,numSplt);


%% timekeeping
scrptSta = tic;


%% step 1: loop through all patients and repetitions to acquire minimal trial number
tmp_trl = cell(numPat,numSplt);

for iPat = 1 : numPat
    fprintf('... step 1: patient %d\n', iPat);
  

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
    
    %% create binary vector indexing clean sleep data points
    godfltr = all([scoring; ~artfltr]);
    
    clear scoring artfctfltr
    
    
    %% process ripple events
    inRipple = load(fullfile(dirRoot,'Ripples',sprintf('pat%02d_HIPP_ripples.mat',iPat)));      %% ripple events
        
    tfg     = [];
    tfg.trl = [inRipple.EvtInfo.maxTime' + round(out.param.timePad(1) * fsample),...
               inRipple.EvtInfo.maxTime' + round(out.param.timePad(2) * fsample),...
               ones(inRipple.EvtInfo.numEvt,1) * round(out.param.timePad(1) * fsample),...
               inRipple.EvtInfo.duration'];

    tfg.trl(tfg.trl(:,1) < 1 | tfg.trl(:,2) > datalen | inRipple.rejects,:)     = [];   %% rmv trials outside data range
    tfg.trl(arrayfun(@(x,y) ~all(godfltr(x:y)),tfg.trl(:,1),tfg.trl(:,2)),:)    = [];   %% rmv trials overlapping with artifacts
    
    numRipple  = size(tfg.trl,1);
    rippleTrl  = tfg.trl;
    
    clear inRipple
   
    
    %% process surrogates
    ctrlTrl = cell(1,numRep);
    numCtrl = nan(1,numRep);
    
    for iRep = 1 : numRep
        inCtrl = load(fullfile(dirRoot,'Control events',num2str(iRep),sprintf('pat%02d_HIPP_controlNREM.mat',iPat)));
                
        tfg     = [];
        tfg.trl = [inCtrl.nonEvtInfo.maxTime' + round(out.param.timePad(1) * fsample),...
                   inCtrl.nonEvtInfo.maxTime' + round(out.param.timePad(2) * fsample),...
                   ones(inCtrl.nonEvtInfo.numEvt,1) * round(out.param.timePad(1) * fsample),...
                   (inCtrl.nonEvtInfo(1).endTime' - inCtrl.nonEvtInfo(1).staTime') / fsample];

        tfg.trl(tfg.trl(:,1) < 1 | tfg.trl(:,2) > datalen,:) = [];                              %% Remove trials outside data range
        tfg.trl(arrayfun(@(x,y) ~all(godfltr(x:y)),tfg.trl(:,1),tfg.trl(:,2)),:) = [];    %% rmv trials overlapping with artifacts
        
        numCtrl(iRep) = size(tfg.trl,1);
        ctrlTrl{iRep} = tfg.trl;

        clear inCtrl
    end     %% iRep
    
    
    %% 3-way split based on ripple duration & balancing of trials
    for iSplt = 1 : numSplt-1
        out.splt.ripple(iPat,iSplt) = prctile(rippleTrl(:,4),100*iSplt/numSplt);
        out.splt.ctrl(iPat,:,iSplt) = cellfun(@(x) prctile(x(:,4),100*iSplt/numSplt),ctrlTrl);
    end
        
    for iSplt = 1 : numSplt
        
        idxRipple   = [];
        idxCtrl     = cell(1,numRep);
        
        switch iSplt
            case 1
                idxRipple   = find(rippleTrl(:,4) < out.splt.ripple(iPat,1));
                idxCtrl     = cellfun(@(x) find(x(:,4) < out.splt.ctrl(iPat,iRep,1)),ctrlTrl,'UniformOutput',0);
                
            case 2
                idxRipple   = find((rippleTrl(:,4) >= out.splt.ripple(iPat,1)) & ...
                                   (rippleTrl(:,4) <= out.splt.ripple(iPat,2)));
                idxCtrl     = cellfun(@(x) find((x(:,4) >= out.splt.ctrl(iPat,iRep,1)) & ...
                                                (x(:,4) <=out.splt.ctrl(iPat,iRep,2))),ctrlTrl,'UniformOutput',0);         
                
            case 3
                idxRipple   = find(rippleTrl(:,4) > out.splt.ripple(iPat,2));
                idxCtrl     = cellfun(@(x) find(x(:,4) > out.splt.ctrl(iPat,iRep,2)),ctrlTrl,'UniformOutput',0);
        end

        %... extract minimal trial number
        numRipple           = numel(idxRipple);
        numCtrl             = cellfun(@(x) size(x,1),idxCtrl);
        minTrl              = min([numRipple numCtrl]);
        out.dof(iPat,iSplt) = minTrl;

        
        %% balance out trials
        tmp_trl{iPat,iSplt}.ripple = rippleTrl(idxRipple(randperm(numRipple,minTrl)),1:3);
           
        for iRep = 1 : numRep
            tmp_trl{iPat,iSplt}.surrgt{iRep} = ctrlTrl{iRep}(idxCtrl{iRep}(randperm(numCtrl(iRep),minTrl)),1:3);
        end

        clear idxRipple idxCtrl

    end     %% iSplt       
    
    clear rippleTrl ctrlTrl
    
end         %% iPat


%% step 2: calculate pdc for ripples
for iSplt = 1 : numSplt
    rippleFFT = cell(numPat,1);

    %% tfr calculation
    for iPat = 1 : numPat
        fprintf('... step 2a (split %d): patient %d\n', iSplt, iPat);


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
        tfg.trl     = tmp_trl{iPat,iSplt}.ripple;
        rippleTrl   = ft_redefinetrial(tfg,godData);

        clear godData

        
        %% (optional) z-score
        if out.param.do_zscore == 1
            tmpMat          = cell2mat(rippleTrl.trial);
            tmpMean         = mean(tmpMat,2);
            tmpStd          = std(tmpMat,1,2);
            rippleTrl.trial = cellfun(@(x) (x-tmpMean) ./ tmpStd, rippleTrl.trial,'UniformOutput',0);
        end


        %% calculate short-time fft
        rippleFFT{iPat} = nan(out.dof(iPat,iSplt),numTime,numCh,numFreq);
        
        for iTime = 1 : numTime
            tfg        = [];
            tfg.toilim = out.time(iTime) + out.param.pdcPad;
            tmp_trl2   = ft_redefinetrial(tfg,rippleTrl);

            %--- perform fft
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

            tmp_trl2 = []; tmp_fft = [];
        end

        clear rippleTrl

    end     %% iPat

    %% partial directed coherence calculation
    fprintf('... step 2b (split %d): calculate pdc for ripples\n', iSplt);
    
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

        out.pdc.ripple(iSplt,1,iTime,:) = squeeze(tmp_pdc.pdcspctrm(1,2,:));
        out.pdc.ripple(iSplt,2,iTime,:) = squeeze(tmp_pdc.pdcspctrm(2,1,:)); 

        clear tmpData tmp_pdc

    end     %% iTime

    clear rippleFFT
    
end     %% iSplt


%% Step 3: process surrogates
for iSplt = 1 : numSplt
    for iRep = 1 : numRep

        ctrlFFT = cell(numPat,1);
        
        %% tfr calculation
        for iPat = 1 : numPat
            fprintf('... step 3a (split %d, rep %d): patient %d\n', iSplt, iRep, iPat);


            %% gather data
            inEEG  = cell(numCh,1);

            for iCh = 1 : numTrgt
                inEEG{iCh}  = load(fullfile(dirRoot,'EEGs',sprintf('pat%02d_%s.mat',iPat,out.def.label{iCh})));
            end

            godData         = ft_appenddata([],inEEG{:});
            godData.fsample = fsample;

            clear inEEG

            
            %% prepare control events
            tfg     = [];
            tfg.trl = tmp_trl{iPat,iSplt}.ctrl{iRep};
            ctrlTrl = ft_redefinetrial(tfg,godData);
        
            clear godData
            
            
            %% (optional) z-score
            if out.param.do_zscore == 1
                tmpMat          = cell2mat(ctrlTrl.trial);
                tmpMean         = mean(tmpMat,2);
                tmpStd          = std(tmpMat,1,2);
                ctrlTrl.trial   = cellfun(@(x) (x-tmpMean) ./ tmpStd, ctrlTrl.trial,'UniformOutput',0);
            end


            %% calculate short-time fft
            ctrlFFT{iPat} = nan(out.dof(iPat,iSplt),numTime,numCh,numFreq);
            
            for iTime = 1 : numTime
                tfg        = [];
                tfg.toilim = out.time(iTime) + out.param.pdcPad;
                tmp_trl2    = ft_redefinetrial(tfg,ctrlTrl);

                %--- perform fft
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

                clear tmp_trl2 tmp_fft;
            end

            clear ctrlTrl
            
        end     %% iPat


        %% temporal partial directed coherence calculation
        fprintf('... step 3b (split %d): calculate pdc for control events\n', iSplt);

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
            tmp_pdc              = ft_connectivityanalysis(tfg,tmpData);

            out.pdc.ctrl(iRep,iSplt,1,iTime,:) = squeeze(tmp_pdc.pdcspctrm(1,2,:));
            out.pdc.ctrl(iRep,iSplt,2,iTime,:) = squeeze(tmp_pdc.pdcspctrm(2,1,:)); 

            clear tmpData tmp_pdc

        end     %% iTime

        clear ctrlFFT

    end     %% iRep
end         %% iSplt


%% perform statistics
%% spindle band time-series
idxFreq = arrayfun(@(x) nearest(out.freq,x),[12 16]);

for iSplt = 1 : numSplt
    tmpRipple   = squeeze(nanmean(out.pdc.ripple(iSplt,:,:,idxFreq(1):idxFreq(2)),4));
    tmpCtrl     = squeeze(nanmean(out.pdc.ctrl(:,iSplt,:,:,idxFreq(1):idxFreq(2)),5));
    
    out.stats.spindleBand{iSplt}.ripple     = squeeze(tmpRipple(1,:)-tmpRipple(2,:));
    out.stats.spindleBand{iSplt}.ctrl.dist  = squeeze(tmpCtrl(:,1,:)-tmpCtrl(:,2,:));
    out.stats.spindleBand{iSplt}.ctrl.mean  = squeeze(nanmean(tmpCtrl(:,1,:)-tmpCtrl(:,2,:)))';
    out.stats.spindleBand{iSplt}.ctrl.std   = squeeze(nanstd(tmpCtrl(:,1,:)-tmpCtrl(:,2,:),1,1))';
    out.stats.spindleBand{iSplt}.zval       = (out.stats.spindleBand{iSplt}.ripple - out.stats.spindleBand{iSplt}.ctrl.mean) ./ out.stats.spindleBand{iSplt}.ctrl.std;
end


%% save results
save(fullfile(dirRoot,[svName '.mat']),'-struct','out');


%% timekeeping
fprintf('This analysis took %.2f s\n', toc(scrptSta));
