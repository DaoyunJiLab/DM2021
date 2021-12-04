function [pinfo,data] = DataManager_FindWaveProp(pinfo,data, cellind, vv)
%%find firing properties for each spike: channel 1,2,3,4 amplitudes, maximum height, maximum width, all average values
%%all amplitude values in AD units (need gain to switch to mV), time values in data point (need sampling frequency to switch to second)
if (~isempty(cellind))
   nspike = numel(pinfo.general.animalname);  if (~isfield(pinfo, 'firing')) pinfo.firing = []; end
   if (~isfield(pinfo.firing, 'champ')) pinfo.firing.champ = cell(1,nspike); end
   if (~isfield(pinfo.firing, 'maxwd')) pinfo.firing.maxwd = cell(1,nspike); end
   if (~isfield(pinfo.firing, 'maxht'))    pinfo.firing.maxht = cell(1,nspike); end
   if (~isfield(pinfo.firing, 'posamp'))    pinfo.firing.posamp = cell(1,nspike); end
   if (~isfield(pinfo.firing, 'negamp'))    pinfo.firing.negamp = cell(1,nspike); end
   if (~isfield(pinfo.firing, 'npampratio'))    pinfo.firing.npampratio = cell(1,nspike); end
   if (~isfield(pinfo.firing, 'hfwd'))    pinfo.firing.hfwd = cell(1,nspike); end
   if (~isfield(pinfo.firing, 'hlwd'))    pinfo.firing.hlwd = cell(1,nspike); end
   if (~isfield(pinfo.firing, 'hlht'))    pinfo.firing.hlht = cell(1,nspike); end
   if (~isfield(pinfo.firing, 'maxposamp'))    pinfo.firing.maxposamp = cell(1,nspike); end
   if (~isfield(pinfo.firing, 'maxnegamp'))    pinfo.firing.maxnegamp = cell(1,nspike); end
   if (~isfield(pinfo.firing, 'maxnpampratio'))    pinfo.firing.maxnpampratio = cell(1,nspike); end
   if (~isfield(pinfo.firing, 'maxhfwd'))    pinfo.firing.maxhfwd = cell(1,nspike); end
   if (~isfield(pinfo.firing, 'maxhlwd'))    pinfo.firing.maxhlwd = cell(1,nspike); end
   if (~isfield(pinfo.firing, 'maxhlht'))    pinfo.firing.maxhlht = cell(1,nspike); end
end

for (jjjk = 1:numel(cellind))
    i = cellind(jjjk);
    disp(strcat('-----> wave prob ---', pinfo.general.parmfile{i}));
    gain = pinfo.general.gain{i}; fs = pinfo.parm.fs(i); spikebasepoint = pinfo.parm.spikebasepoint(i);
    if (~isempty(pinfo.general.wavefile{i}))
        wavefile = pinfo.general.wavefile{i};
      if exist(wavefile, 'file') 
        [spiketime, wv] = ReadSpikeWaveSPW(wavefile);
        %%%%%spiketime in second, wv = [32, 4, nspike]
        nspike = numel(spiketime);
        if (nspike >0)
           hfwd = zeros(4,1); hlwd = zeros(4,1); hlht = zeros(4,1); posamp = zeros(4,1); negamp = zeros(4,1); ampratio = zeros(4,1); 
           %calculate average firing
           wvavg = sum(wv,3)/nspike;
           if ~isempty(wvavg)
              for (j = 1:4)
                  [hfwd(j), hlwd(j), hlht(j), posamp(j), negamp(j), ampratio(j)] = findparam(wvavg(:,j)*gain(j), fs, spikebasepoint);
              end
           end
           if (vv == 1) plotspikewaveavg(wavefile, wvavg, gain); end
           wvavg = []; wv = []; spiketime = [];
           pinfo.firing.posamp{i} = posamp;
           pinfo.firing.negamp{i} = negamp;
           pinfo.firing.npampratio{i} = ampratio;
           pinfo.firing.hfwd{i} = hfwd;
           pinfo.firing.hlwd{i} = hlwd;
           pinfo.firing.hlht{i} = hlht;
           [pinfo.firing.maxposamp{i}, ikk] = max(posamp);
           pinfo.firing.maxnegamp{i} = negamp(ikk);
           pinfo.firing.maxnpampratio{i} = ampratio(ikk);
           pinfo.firing.maxhfwd{i} = hfwd(ikk);
           pinfo.firing.maxhlwd{i} = hlwd(ikk);
           pinfo.firing.maxhlht{i} = hlht(ikk);
        end
      end
    end
    if (~isempty(pinfo.general.parmfile{i}))
        %get spike parm file (pinfo.filename.parmfile{i})
        parmfile = pinfo.general.parmfile{i};
        S = load(parmfile, '-ascii'); %this loading method suppose not sensitive to different parameter columns
        if (~isempty(S))
           pinfo.firing.champ{i} = [mean(S(:,2))*gain(1) mean(S(:,3))*gain(2) mean(S(:,4))*gain(3) mean(S(:,5))*gain(4)];
           pinfo.firing.maxwd{i} = mean(S(:,6))/fs;
           pinfo.firing.maxht{i} = mean(S(:,7))*max(gain);
           S = [];
        end
    end
end

function [hfwd, hlwd, hlht, posamp, negamp, ampratio] = findparam(wvavg, fs, spikebasepoint)
if (spikebasepoint >= 1)
    baseline = mean(wvavg(1:spikebasepoint)); %take the first xx points as baseline
else
    baseline = 0;
end
wvavg = wvavg - baseline;
[posamp, postime] = max(wvavg); [negamp, negtime] = min(wvavg); 
if (posamp ~= 0) 
    ampratio = abs(negamp/posamp); 
else
    ampratio = NaN;
end
hlht = posamp - negamp;
hlwd = (negtime - postime)/fs; 

%%%pos peak half width in s
hfwd = NaN; hftime1 = []; hftime2 = [];
if (negtime > postime)
    for (i = 1:postime)
        if (wvavg(i) >= posamp/2)
            hftime1 = i; break
        end 
    end
    for (i = postime:negtime)
        if (wvavg(i) <= posamp/2)
            hftime2 = i; break
        end
   end
   if ((~isempty(hftime1)) && (~isempty(hftime2)))
      hfwd = (hftime2 - hftime1)/fs;
   end
end

function plotspikewaveavg(wavefile, wvavg, gain)
[npoint,nch] = size(wvavg);
tlim = [-2 npoint+2]/32.556; maxV = max([max(gain(1)*wvavg(:,1)) max(gain(2)*wvavg(:,2)) max(gain(3)*wvavg(:,3)) max(gain(4)*wvavg(:,4))]);
minV = min([min(gain(1)*wvavg(:,1)) min(gain(2)*wvavg(:,2)) min(gain(3)*wvavg(:,3)) min(gain(4)*wvavg(:,4))]);
vlim = [minV-5 maxV+5];
hf = figure('Name', strcat(wavefile, '_Averagefirings'), 'NumberTitle', 'off');
ha1 = axes('Parent', hf, 'NextPlot', 'add', 'Position', [0.1 0.1 0.15 0.9],...
    'XLim', tlim, 'YLim', vlim); xlabel('ms'); ylabel('uV');
ha2 = axes('Parent', hf, 'NextPlot', 'add', 'Position', [0.35 0.1 0.15 0.9],...
    'XLim', tlim, 'YLim', vlim);
ha3 = axes('Parent', hf, 'NextPlot', 'add', 'Position', [0.6 0.1 0.15 0.9],...
    'XLim', tlim, 'YLim', vlim);
ha4 = axes('Parent', hf, 'NextPlot', 'add', 'Position', [0.85 0.1 0.15 0.9],...
    'XLim', tlim, 'YLim', vlim);
plot ([1:npoint]/32.556, gain(1)*wvavg(:,1), 'Parent', ha1, 'LineWidth', 2, 'Color', [1 0 0]);
xlabel('ms'); ylabel('uV');
plot ([1:npoint]/32.556, gain(2)*wvavg(:,2), 'Parent', ha2, 'LineWidth', 2, 'Color', [1 0 0]);
plot ([1:npoint]/32.556, gain(3)*wvavg(:,3), 'Parent', ha3, 'LineWidth', 2, 'Color', [1 0 0]);
plot ([1:npoint]/32.556, gain(4)*wvavg(:,4), 'Parent', ha4, 'LineWidth', 2, 'Color', [1 0 0]);


