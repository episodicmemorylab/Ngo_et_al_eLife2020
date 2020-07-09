%% Ngo et al., eLife 2020: Sleep spindles mediate hippocampal-neocortical coupling during long-duration ripples
%
% Average HIPP activity time-locked to offline identified ripple events
%
% Requirements:
% - Fieldtrip added to Matlab search path
% - line 26: "additional_functions" folder added to search path
% - line 28: specification of root path containing "EEGs" and "Ripples"
%   folder (see https://osf.io/3hpvr/)
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
addpath('E:\GitHub\Ngo_et_al_eLife2020\additional_functions')   %% specify path

dirRoot = 'E:\work_uob\sleesio\upload';                         %% specify path

svName  = 'rippleAverage';                                      %% filename for saved data


%% fundamental parameters
numPat  = 14;   %% number of patients
numCh   = 1;    %% number of channels


%% output structure
out = [];

out.def.label  = {'HIPP'};

out.param.artfctPad     = [-0.25 0.25];     %% padding applied to artifacts
out.param.stageoi       = [2 3 4];          %% sleep stages of interest, i.e. NREM.
out.param.minNREMLen    = 3;                %% minimal required length of NREM to be analyzed
out.param.avgPad        = [-0.5 0.5];       %% temporal padding around minimum of ripple events
out.param.baseln        = [-0.5 -0.4];      %% interval for baseline correction

out.time = out.param.avgPad(1) : 0.001 : out.param.avgPad(2);   numTime = numel(out.time);  %% time vector

out.dof    = nan(numPat,numCh);             %% degrees of freedom, i.e. numbers of analyzed ripples per patient and channel
out.indiv  = nan(numPat,numCh,numTime);     %% individual averages per patient
out.grdavg = nan(2,numCh,numTime);          %% mean and sem across patients


%% timekeeping
scrptSta = tic;


%% main loop across patients
for iPat = 1 : numPat
    fprintf('analyse patient %d\n', iPat);
    
            
    %% get useful parameters
    tmplt   = load(fullfile(dirRoot,'EEGs',sprintf('pat%02d_NC_supplement.mat',iPat)));
    fsample = tmplt.fsample;
    datalen = tmplt.datalen;
    scoring = ismember(tmplt.scoring,out.param.stageoi);
    
    clear tmplt
    
    
    %% loop across channels
    for iCh = 1 : numCh
        %% load data
        inEEG       = load(fullfile(dirRoot,'EEGs',sprintf('pat%02d_%s.mat',iPat,out.def.label{iCh})));  %% EEG
        inRipple  = load(fullfile(dirRoot,'Ripples',sprintf('pat%02d_%s_ripples.mat',iPat,out.def.label{iCh})));                              %% Spindle events
        
        
        %% filter data
        tfg                     = [];
        tfg.hpfilter            = 'yes';
        tfg.hpfreq              = 0.3;
        tfg.hpfiltord           = 3;
        tfg.lpfilter            = 'yes';
        tfg.lpfreq              = 150;
        tfg.lpfiltord           = 5;
        tfg.lpinstabilityfix    = 'reduce';
        tfg.bsfilter            = 'yes';
        tfg.bsfreq              = [48 52];

        inEEG = ft_preprocessing(tfg,inEEG);
        
        
        %% prepare binary vector indexing clean NREM data points
        inSupplmt = load(fullfile(dirRoot,'EEGs',sprintf('pat%02d_%s_supplement.mat',iPat,out.def.label{iCh})));

        %.. add padding, ensure padding is within data range
        artfct                          = inSupplmt.artifacts + round(out.param.artfctPad * fsample);
        artfct(artfct(:,1) < 1,1)       = 1;
        artfct(artfct(:,2) > datalen,2) = datalen;
        
        %.. create binary artifact filter
        artfctfltr = hvn_createBnrySignal(artfct,datalen);
        
        clear artfct

        %.. remove NREM intervals shorter than minNREMLen
        cleanNREM                = hvn_extrctBnryBouts(~artfctfltr);
        rmvIdx                   = diff(cleanNREM,1,2) < round(out.param.minNREMLen * fsample);
        rmvSamples               = cell2mat(reshape(arrayfun(@(x,y)  x:y, cleanNREM(rmvIdx,1), cleanNREM(rmvIdx,2),'UniformOutput',0),1,sum(rmvIdx)));
        artfctfltr(rmvSamples)   = 1;

        godfltr = all([scoring; ~artfctfltr]);

        clear inSupplmt artfctfltr artfct cleanNREM rmvIdx rmvSamples


        %% segmentation
        tfg     = [];
        tfg.trl = [inRipple.EvtInfo.maxTime' + round(out.param.avgPad(1) * fsample),...
                   inRipple.EvtInfo.maxTime' + round(out.param.avgPad(2) * fsample),...
                   ones(inRipple.EvtInfo.numEvt,1) * round(out.param.avgPad(1) * fsample)];

        tfg.trl(tfg.trl(:,1) < 1 | tfg.trl(:,2) > datalen | inRipple.rejects,:) = [];          %% remove trials outside data range
        tfg.trl(arrayfun(@(x,y) ~all(godfltr(x:y)),tfg.trl(:,1),tfg.trl(:,2)),:) = [];          %% remove trials overlapping with artifacts

        rippleTrl = ft_redefinetrial(tfg,inEEG);
        
        clear inEEG
        
        
        %% average across trials
        rippleAvg = ft_timelockanalysis([],rippleTrl);
        
        clear rippleTrl


        %% baseline correction
        tfg             = [];
        tfg.baseline    = out.param.baseln;
        rippleAvg       = ft_timelockbaseline(tfg,rippleAvg);


        %% save results in output structure
        out.dof(iPat,iCh)      = rippleAvg.dof(1);
        out.indiv(iPat,iCh,:)  = squeeze(rippleAvg.avg);
        
        clear rippleAvg
    end    
end


%% grand average across patients
out.grdavg(1,:,:) = nanmean(out.indiv,1);                   %% mean
out.grdavg(2,:,:) = nanstd(out.indiv,1,1) ./ sqrt(numPat);  %% sem


%% save results
save(fullfile(dirRoot,[svName '.mat']),'-struct','out');


%% plot grand average
plt_xlim = [-0.5 0.5];
plt_ylim = [-80 70];

figure('Name','Ripple average');
    
boundedline(out.time,squeeze(out.grdavg(1,iCh,:)),squeeze(out.grdavg(2,iCh,:)))
    
xlabel('Time (s)');
xlim(plt_xlim);

ylabel('Amplitude (uV)');
ylim(plt_ylim)


%% timekeeping
fprintf('This analysis took %.2f s\n', toc(scrptSta));
