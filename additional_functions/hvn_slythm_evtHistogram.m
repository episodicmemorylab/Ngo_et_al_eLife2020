%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Slythm_CrossEvtHistogram
% Created by H.V.-V. Ngo
%
% Calculates histogram of (discrete) events time-locked to a specific event/trigger
% - Requires Fieldtrip toolbox
% - Time-bin i summarises activity from [t_i, t_i+1), i.e. up to time-bin
%   i+1 excluding that time bin itself. Only the last bin comprises
%   [t_end-1 t_end]
%
% Usage: outHist = Slythm_CrossEvtHistogram(inEvent, inLabel, inTrigger, datalen, fsample, window, binsize, doNorm)
%
% Parameters:
% inEvent   = (N x 1) cell containing the discrete events times (as row
%             vectors) to be examined, N = number of channels
% inLabel   = (N x 1) cell with N channel labels
% inTrigger = (M x 1) array with trigger/events the inEvent data will be
%             examined relative to, M is the number of events
% datalen   = length of original data recording
% fsample   = sampling rate of original data recording
% window    = [preTime, postTime], interval to be analysed centered on
%             inTrigger data, preTime < 0 means interval begins before the
%             corresponding trigger
% binsize   = resolution of Cross-Event-Histogram
% doNorm    = 0 | 1 | 2, indicating whether results should be normalized
%             0 = deactivated
%             1 = normalise by number of triggers * 100
%             2 = normalise by number of captured events * 100
% doSmooth  = odd integer, for smoothing of resulting histogram,
%             default = 1 i.e. no smoothing
%
% To-do
% [ ] Complete bookkeeping section
%
% Last update: 18-06-12 by HVN


function outHist = hvn_slythm_evtHistogram(inEvent, inLabel, inTrigger, datalen, fsample, window, binsize, doNorm, doSmooth)

%% timekeeping
statime = tic;


%% bookkeeping
numCh = size(inEvent,1);

switch nargin
    case 7
        doNorm = 0;
        doSmooth = 1;
    case 8
        doSmooth = 1;
end


%% Prepare outHist variable
outHist         = [];
outHist.time    = window(1,1): binsize : window(1,2);
outHist.avg     = nan(numCh,numel(outHist.time));
outHist.dof     = [];
outHist.label   = inLabel;
outHist.binsize = binsize;

if~any(ismember(outHist.time,0))
    warning('Specified histogram bins do not include zero-bin, this might cause distortions');
end


%% create fieldtrip structure representing inEvents data as binary events
dataEvent               = [];
dataEvent.fsample       = fsample;
dataEvent.label         = inLabel;
dataEvent.trial         = {zeros(numCh,datalen)};
dataEvent.time          = {(0:datalen-1) ./ fsample};
dataEvent.sampleinfo    = [1 datalen];

for iCh = 1 : numCh
    dataEvent.trial{1,1}(iCh,inEvent{iCh,1}) = 1;
end


%% Prepare trial config for segmentation
tfg             = [];
tfg.minlength   = window(1,2) - window(1,1);
tfg.trl         = [inTrigger+round(window(1,1) * fsample), inTrigger+round(window(1,2) * fsample), ones(size(inTrigger,1),1) * round(window(1,1) * fsample)];

tfg.trl(tfg.trl(:,1) < 1 | tfg.trl(:,2) > datalen,:) = [];                  % ensure trials do not exceed data range

tmpTrl = ft_redefinetrial(tfg,dataEvent);                                   % create trial data


%% perform averaging time-locked to specified trigger
tmpAvg        = ft_timelockanalysis([],tmpTrl);
tmpAvg.avg    = tmpAvg.avg .* tmpAvg.dof;           %% multiply with degrees of freedom to get absolute number


%% Create histogram
tZeroAvg    = nearest(tmpAvg.time,0);
tZeroHist   = nearest(outHist.time,0);

%--- Left tail, if existing
if tZeroHist > 1
    numBin = tZeroHist-1;
    tmpBin = [arrayfun(@(x) nearest(tmpAvg.time,x),outHist.time(1:numBin)) tZeroAvg];
    
    for iBin = 1 : numBin-1
        outHist.avg(:,iBin) = sum(tmpAvg.avg(:,tmpBin(iBin):tmpBin(iBin+1)-1),2);
    end
    
    outHist.avg(:,numBin) = sum(tmpAvg.avg(:,tmpBin(numBin):tmpBin(numBin+1)),2);   %% last left tail bin including t = 0
end

%-- Right tail, if existing
if tZeroHist <= numel(outHist.time)
        numBin = numel(outHist.time)-tZeroHist;
        tmpBin = [tZeroAvg arrayfun(@(x) nearest(tmpAvg.time,x),outHist.time(tZeroHist+1:end))];
        
    for iBin = 1 : numBin
        outHist.avg(:,iBin+tZeroHist-1) = sum(tmpAvg.avg(:,tmpBin(iBin)+1:tmpBin(iBin+1)),2);
    end
end

%% Remove last bin
outHist.avg(:,end)  = [];
outHist.time(:,end) = [];


%% Smooth data
outHist.avg = smoothdata(outHist.avg,2,'movmean',doSmooth);


%% Normalise, if desired
switch doNorm
    case 1      %% Normalise by number of trigger events
        outHist.dof = repmat(tmpAvg.dof(1,1),size(outHist.avg));
        outHist.avg = 100 * outHist.avg ./ outHist.dof;
    case 2      %% Normalise by number of captured events
        outHist.dof = repmat(sum(outHist.avg,2),1,numel(outHist.time));
        outHist.avg = 100 * outHist.avg ./ outHist.dof;
    otherwise
        outHist.dof = repmat(tmpAvg.dof(1,1),size(outHist.avg));
end


% for iCh = 1 : numCh
%     %--- resample data for histcounts function
%     tmpHist = [];
%     for iTime = 1 : size(tmpAvg.time,2)
%         if tmpAvg.avg(iCh,iTime) ~= 0
%            tmpHist = [tmpHist, ones(1,round(tmpAvg.avg(iCh,iTime))) * tmpAvg.time(iTime)];
%         end
%     end
%     [resHist, ~] = histcounts(tmpHist,window(1) : binsize : window(2));
%     
%     %--- Optional: Normalize data
%     switch doNorm
%         case 1
%             outHist.dof = [outHist.dof; unique(tmpAvg.dof)];
%             resHist     = 100 * resHist ./ unique(tmpAvg.dof);
%         case 2
%             outHist.dof = [outHist.dof; sum(resHist)];
%             resHist     = 100 * resHist ./ sum(resHist);
%     end
%     
%     outHist.avg = [outHist.avg; resHist];
% end


fprintf('This shit got done in only %.2f s\n', toc(statime));
end
