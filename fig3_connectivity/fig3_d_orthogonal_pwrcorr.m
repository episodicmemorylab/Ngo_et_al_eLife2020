%% Ngo et al., eLife 2020: Sleep spindles mediate hippocampal-neocortical coupling during long-duration ripples
%
% Calculate NC-HIPP orthogonalized power correlation time-locked to ripples
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

svName  = 'powercorrelation';                           %% filename for saved data


%% fundamental parameters
numType = 2;
numPat  = 14;
numCh   = 2;
numRep  = 100;

fsample = 1000;


%% results structure
out = [];

out.def.type    = {'Ripple'; 'Ctrl'};
out.def.label   = {'NC'; 'HIPP'};

out.dim.pwrcorr = '(Rep) x Freq x Time';

out.param.artfctPad     = [-0.25 0.25];     %% padding applied to artifacts
out.param.stageoi       = [2 3 4];          %% sleep stages of interest, i.e. NREM.
out.param.minNREMLen    = 3;                %% minimal required length of NREM to be analyzed
out.param.timePad       = [-2.1 2.1];       %% temporal padding for segmentation

out.time = -1.0 : 0.02 : 1.0;   numTime = numel(out.time);          %% time vector
out.freq = 1 : 0.5 : 20;        numFreq = numel(out.freq);          %% frequency vector

out.param.mwCycl                        = ceil(out.freq * 0.5);     %% number of cycles per frequency for morlet wavelets
out.param.mwCycl(out.param.mwCycl < 5)  = 5;
out.param.mwCycl(1:8)                   = [2 3 3 3 3 3 4 4];

out.dof         = [];
out.pwrcorr.ripple  = nan(numFreq,numTime);
out.pwrcorr.ctrl    = nan(numRep,numFreq,numTime);

out.stats.map.zval   = zeros(numFreq,numTime);
out.stats.map.pval   = ones(numFreq,numTime);
out.stats.map.mask   = zeros(numFreq,numTime,'logical');


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


%% step 2: calculate coherence for ripples
ripple_x   = cell(numPat,1); 
ripple_y   = cell(numPat,1);
ripple_x_y = cell(numPat,1);
ripple_y_x = cell(numPat,1);

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
    
    
    %% compute complex tfr
    tfg             = [];
    tfg.output      = 'fourier';
    tfg.method      = 'wavelet';
    tfg.toi         = out.time;
    tfg.foi         = out.freq;
    tfg.pad         = 'nextpow2';
    tfg.width       = out.param.mwCycl;
    tfg.keeptrials  = 'yes';
    tfg.polyremoval = 1;
    tmp_tfr         = ft_freqanalysis(tfg,rippleTrl);
    
    clear rippleTrl
    
    %% orthogonalise 
    x   = squeeze(tmp_tfr.fourierspctrm(:,1,:,:));
    y   = squeeze(tmp_tfr.fourierspctrm(:,2,:,:));
    x_y = nan(size(x));
    y_x = nan(size(x));
    
    for iTrl = 1 : size(tmpTFR.fourierspctrm,1)
        [x_y(iTrl,:,:), y_x(iTrl,:,:)]  = mt_pwrcorr_ortho('hipp',squeeze(x(iTrl,:,:)),squeeze(y(iTrl,:,:)),0); 
    end
    
    %.. absolute value
    ripple_x{iPat}   = abs(x).^2;
    ripple_y{iPat}   = abs(y).^2;
    ripple_x_y{iPat} = abs(x_y).^2;
    ripple_y_x{iPat} = abs(y_x).^2;
    
    clear rippleTrl tmp_tfr x y x_y y_x
    
end


%% power correlation calculation
fprintf('... step 2b: calculate power correlation for ripples\n');

ripple_x   = cell2mat(ripple_x);
ripple_y   = cell2mat(ripple_y);
ripple_x_y = cell2mat(ripple_x_y);
ripple_y_x = cell2mat(ripple_y_x);

for iFreq = 1 : numFreq
    for iTime = 1 : numTime
        tmpCorr     = corr(ripple_x(:,iFreq,iTime),ripple_y_x(:,iFreq,iTime));
        tmpCorr2    = corr(ripple_x_y(:,iFreq,iTime),ripple_y(:,iFreq,iTime));

        out.pwrcorr.ripple(iFreq,iTime) = mean([tmpCorr,tmpCorr2]);            
    end
end

clear ripple_x ripple_y ripple_x_y ripple_y_x


%% Step 3: process 
for iRep = 1 : numRep
    
    ctrl_x   = cell(numPat,1); 
    ctrl_y   = cell(numPat,1);
    ctrl_x_y = cell(numPat,1);
    ctrl_y_x = cell(numPat,1);
    
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
        

        %% compute complex TFR
        tfg             = [];
        tfg.output      = 'fourier';
        tfg.method      = 'wavelet';
        tfg.keeptrials  = 'yes';
        tfg.toi         = out.time;
        tfg.foi         = out.freq;
        tfg.pad         = 'nextpow2';
        tfg.width       = out.param.mwCycl;
        tfg.polyremoval = 1;
        tmp_tfr         = ft_freqanalysis(tfg,ctrlTrl);

        %% orthogonalise 
        x   = squeeze(tmp_tfr.fourierspctrm(:,1,:,:));
        y   = squeeze(tmp_tfr.fourierspctrm(:,2,:,:));
        x_y = nan(size(x));
        y_x = nan(size(x));

        for iTrl = 1 : size(tmpTFR.fourierspctrm,1)
            [x_y(iTrl,:,:), y_x(iTrl,:,:)]  = mt_pwrcorr_ortho('hipp',squeeze(x(iTrl,:,:)),squeeze(y(iTrl,:,:)),0); 
        end

        %.. absolute value
        ctrl_x{iPat}   = abs(x).^2;
        ctrl_y{iPat}   = abs(y).^2;
        ctrl_x_y{iPat} = abs(x_y).^2;
        ctrl_y_x{iPat} = abs(y_x).^2;

        clear ctrlTrl tmp_tfr x y x_y y_x
        
    end     %% iPat
    
    
    %% calculate coherence
    ctrl_x   = cell2mat(ctrl_x);
    ctrl_y   = cell2mat(ctrl_y);
    ctrl_x_y = cell2mat(ctrl_x_y);
    ctrl_y_x = cell2mat(ctrl_y_x);

    for iFreq = 1 : numFreq
        for iTime = 1 : numTime
            tmpCorr     = corr(ctrl_x(:,iFreq,iTime),ctrl_y_x(:,iFreq,iTime));
            tmpCorr2    = corr(ctrl_x_y(:,iFreq,iTime),ctrl_y(:,iFreq,iTime));

            out.pwrcorr.ctrl(iRep,iFreq,iTime) = mean([tmpCorr,tmpCorr2]);            
        end
    end

    clear ctrl_x ctrl_y ctrl_x_y ctrl_y_x
    
end     %% iRep


%% perform statistics
%% time-resolved coherence map
tmpRipple       = out.pwrcorr.ripple;
tmpCtrl_mean    = squeeze(nanmean(out.pwrcorr.ctrl,1));
tmpCtrl_std     = squeeze(nanstd(out.pwrcorr.ctrl,1,1));

out.stats.map.zval = (tmpRipple - tmpCtrl_mean) ./ tmpCtrl_std;
out.stats.map.pval = 2*normcdf(-abs(out.stats.map.zval));
out.stats.map.mask = outStats.map.pval < 0.05;


%% save results
save(fullfile(dirSave,[filSave '.mat']),'-struct','out');


%% timekeeping
fprintf('This analysis took %.2f s\n', toc(scrptSta));
