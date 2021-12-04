function DataManager_EEGReCompute_Callback
%%Add or change variable values of the selected group, and save into a new database

hf = gcbf; pinfo = getappdata(hf, 'eeg'); data = getappdata(hf, 'eegdata'); tagmark = get(gcbo, 'Tag');
hgroup = getappdata(hf, 'hgroup'); hfield = getappdata(hf, 'hfield');
ok = 1; plotparm = getappdata(hf, 'plotparm'); vv = plotparm.showdetail; %if vv=1, shoe details of computation
cc = plotparm.compute; %if cc=0, assign parameters for computation only; if 1, do the real computation
ow = plotparm.overwrite; %if ow=1, plot into the current database
[writefilename, okk] = getoutputfile(hf, ow);

%get selected cellind
groupselection = getappdata(hgroup, 'selection'); cellind = []; grpind = find(groupselection == 1); 
for (kk = 1:numel(grpind)) cellind = union(cellind, data.grouplist.groupindex{grpind(kk)}); end
if (~isempty(cellind)) && okk
    if (strcmp(tagmark, 'recomputeinit'))
        [pinfo,data] = DataManager_EEGComputeInit_Callback(pinfo,data, cellind, vv);
    elseif (strcmp(tagmark, 'recomputeall'))
        [pinfo,data] = DataManager_EEGComputeAll_Callback(pinfo,data, cellind, vv);
    elseif (strcmp(tagmark, 'spectral'))
        if (cc==0) pinfo = assignspectralparm(pinfo, cellind);  end
        if (cc == 1) [pinfo, data] = DataManager_FindSpectralProp(pinfo, data, cellind, vv); end
    elseif (strcmp(tagmark, 'ripples'))
        if (cc==0) pinfo = assignrippleparm(pinfo, cellind);  end
        if (cc == 1) [pinfo, data] = DataManager_FindRippleProp(pinfo, data, cellind, vv); end
    elseif (strcmp(tagmark, 'theta'))
        if (cc==0) pinfo = assignthetaparm(pinfo, cellind);  end
        if (cc == 1) [pinfo, data] = DataManager_FindThetaProp(pinfo, data, cellind, vv); end
    elseif (strcmp(tagmark, 'spindle'))
        if (cc==0) pinfo = assignspindleparm(pinfo, cellind); end
        if (cc == 1) [pinfo, data] = DataManager_FindSpindleProp(pinfo, data, cellind, vv); end
    elseif (strcmp(tagmark, 'hvs'))
        if (cc==0) pinfo = assignhvsparm(pinfo, cellind); end
        if (cc == 1) [pinfo, data] = DataManager_FindHVSProp(pinfo, data, cellind, vv); end
    elseif (strcmp(tagmark, 'cluster0'))
        if (cc==0) pinfo = assigncluster0parm(pinfo, cellind); end
        if (cc == 1) [pinfo, data] = DataManager_FindCluster0Prop(pinfo, data, cellind, vv); end
    elseif (strcmp(tagmark, 'slowtheta'))
        if (cc==0) pinfo = assignslowthetaparm(pinfo, cellind); end
        if (cc == 1) [pinfo, data] = DataManager_FindSlowThetaProp(pinfo, data, cellind, vv); end    
    end
    %%%%%%%%save and plot the new database
    if (~isempty(pinfo.general.finaldir))
        if (ok)
            eeg = pinfo; eegdata = data; save(writefilename, 'eeg', 'eegdata');
            if (ow)
               iii = get(hf, 'Children'); delete(iii(~strcmp(get(iii, 'Type'), 'uimenu'))); 
               fname = get(hf, 'Name'); ftt = strfind(fname, '__'); 
               newname = fname(1:ftt-1); set(hf, 'Name', newname); hmain = hf; 
            else
               hmain = DataManager_DataManager_Callback; plotparm.linkspike = 0; plotparm.linkeeg = 1; plotparm.linkbehav = 0;
            end
            setappdata(hmain,'plotparm', plotparm);
            %DataManager_PlotEEGDatabase(hmain, eeg, eegdata);
            DataManager_PlotSpikeDatabase(hmain, eeg, eegdata, 'EEG File', 'eegID', '.eegdb');
            set(hmain, 'Name', strcat(get(hmain, 'Name'), '__', writefilename));
            setappdata(hmain, 'eeg', eeg); setappdata(hmain, 'eegdata', eegdata); pinfo = []; data = [];
        end
    else
        disp('-------------> no items in the database!');
    end
else
    disp('--------------> no groups selected, groups do not contain any items, or action cancelled');
end
disp('**********************');

function pinfo = assignslowthetaparm(pinfo, cellind)
pp = {'Enter trough threshold (0-1):'; 'Enter start threshold (0-1):'; 'Enter min trough duration (s):';...
      'Enter max trough duration (s):'; 'Enter min number of peaks:'; 'Enter min event gap (s):'}; %; 'Peak mode? (negative/positive/absolute):'};
norm = 'self'; normmeth = '0-1';
if isfield(pinfo, 'cluster0')
   norm = pinfo.parm.cl0NormMode{cellind(1)}; normmeth = pinfo.parm.cl0NormMeth{cellind(1)};
end
if (strcmp(norm, 'self')) && (strcmp(normmeth, '0-1'))
    def = {'0.02'; '0.03'; '0.05'; '0.12'; '3'; '0.2'}; %; 'negative'};
end
III=inputdlg(pp, 'Parameters for computing slow theta properties', 6, def, 'on');
if (~isempty(III))
   pthr = str2num(III{1}); sthr = str2num(III{2}); mindur = str2num(III{3}); 
   maxdur = str2num(III{4}); minpeaknum = str2num(III{5}); minevgap = str2num(III{6}); 
   nspike = numel(pinfo.general.eegfile); 
   if (~isfield(pinfo.parm, 'slowthetaTroughThres')) pinfo.parm.slowthetaTroughThres = NaN*ones(1, nspike); end
   if (~isfield(pinfo.parm, 'slowthetaStartThres')) pinfo.parm.slowthetaStartThres = NaN*ones(1, nspike); end
   if (~isfield(pinfo.parm, 'slowthetaMinDur')) pinfo.parm.slowthetaMinDur = NaN*ones(1, nspike); end
   if (~isfield(pinfo.parm, 'slowthetaMaxDur')) pinfo.parm.slowthetaMaxDur = NaN*ones(1, nspike); end
   if (~isfield(pinfo.parm, 'slowthetaMinPeakNum')) pinfo.parm.slowthetaMinPeakNum = NaN*ones(1, nspike); end  
   if (~isfield(pinfo.parm, 'slowthetaMinEventGap')) pinfo.parm.slowthetaMinEventGap = NaN*ones(1, nspike); end
   ncell = numel(cellind);
   %for (i = 1:ncell) pinfo.parm.slowthetaPeakMode{cellind(i)} = peakmode; end
   pinfo.parm.slowthetaTroughThres(cellind) = pthr*ones(1,ncell);
   pinfo.parm.slowthetaStartThres(cellind) = sthr*ones(1,ncell);
   pinfo.parm.slowthetaMinDur(cellind) = mindur*ones(1,ncell); 
   pinfo.parm.slowthetaMaxDur(cellind) = maxdur*ones(1,ncell);
   pinfo.parm.slowthetaMinPeakNum(cellind) = minpeaknum*ones(1,ncell); 
   pinfo.parm.slowthetaMinEventGap(cellind) = minevgap*ones(1,ncell);
end

function pinfo = assigncluster0parm(pinfo, cellind)
p = {'Normalization mode (self/run/none'; 'Normalization method(Z-transform/0-1/mean-subtraction):'; ... 
     'Smoothing (yes/no):'; 'Smoothing sigma:'; 'Smoothing N of sigmas:';'Time bin size (s):';}; 
d = {'none'; '0-1'; 'no'; '2'; '5'; '0.00003'}; ook = 1;
II = inputdlg(p, 'Binning/normalization/smoothing parameters', 6, d, 'on'); %%%resizable window
if ~isempty(II)
  mode = 'none';
  if strncmp(II{1}, 'run', 1)
    mode = 'run';
  elseif strncmpi(II{1}, 'self', 1)
    mode = 'self';
  end
  norm = 'Z-transform';
  if strncmp(II{2}, '0', 1) || strncmp(II{2}, '1', 1)
    norm = '0-1';
  elseif strncmpi(II{1}, 'mean', 1)
    norm = 'mean-subtraction';
  end
  smooth = 'no';
  if strncmpi(II{3}, 'yes', 1)
    smooth = 'yes';
  end
  sigma = str2num(II{4}); Nsigma = str2num(II{5}); binsize = str2num(II{6});
else 
  ook = 0;
end

pp = {'Start threshold:'; 'Peak threshold:'; 'Max gap (s):'; 'Min duration (s)'; 'Max duration (s)'};
if (strcmp(mode, 'self'))
    def = {'0.15'; '0.6'; '0.02'; '0.05'; '0.2'};
elseif (strcmp(mode, 'run'))
    def = {'0.15'; '0.6'; '0.02'; '0.05'; '0.2'};
else
    def = {'0.05'; '0.3'; '0.005'; '0.05'; '0.2'};
end
III=inputdlg(pp, 'Cluster0 event detection parameters', 2, def, 'on');
if (~isempty(III))
   pthr = str2num(III{2}); sthr = str2num(III{1}); mgap = str2num(III{3}); mindur = str2num(III{4}); maxdur = str2num(III{5});
else
    ook = 0;
end
if ook
   nspike = numel(pinfo.general.eegfile);
   if (~isfield(pinfo.parm, 'cl0NormMode')) pinfo.parm.cl0NormMode = cell(1, nspike); end
   if (~isfield(pinfo.parm, 'cl0NormMeth')) pinfo.parm.cl0NormMeth = cell(1, nspike); end
   if (~isfield(pinfo.parm, 'cl0Smooth')) pinfo.parm.cl0Smooth = cell(1, nspike); end
   if (~isfield(pinfo.parm, 'cl0Sigma')) pinfo.parm.cl0Sigma = NaN*ones(1, nspike); end
   if (~isfield(pinfo.parm, 'cl0NSigma')) pinfo.parm.cl0NSigma = NaN*ones(1, nspike); end
   if (~isfield(pinfo.parm, 'cl0Binsize')) pinfo.parm.cl0Binsize = NaN*ones(1, nspike); end
   if (~isfield(pinfo.parm, 'cl0StartThres')) pinfo.parm.cl0StartThres = NaN*ones(1, nspike); end
   if (~isfield(pinfo.parm, 'cl0PeakThres')) pinfo.parm.cl0PeakThres = NaN*ones(1, nspike); end
   if (~isfield(pinfo.parm, 'cl0MaxGap')) pinfo.parm.cl0MaxGap = NaN*ones(1, nspike); end
   if (~isfield(pinfo.parm, 'cl0MinDur')) pinfo.parm.cl0MinDur = NaN*ones(1, nspike); end
   if (~isfield(pinfo.parm, 'cl0MaxDur')) pinfo.parm.cl0MaxDur = NaN*ones(1, nspike); end
   ncell = numel(cellind);
   for (i = 1:ncell) pinfo.parm.cl0NormMode{cellind(i)} = mode; end
   for (i = 1:ncell) pinfo.parm.cl0NormMeth{cellind(i)} = norm; end
   for (i = 1:ncell) pinfo.parm.cl0Smooth{cellind(i)} = smooth; end
   pinfo.parm.cl0Sigma(cellind) = sigma*ones(1,ncell);
   pinfo.parm.cl0NSigma(cellind) = Nsigma*ones(1,ncell);
   pinfo.parm.cl0Binsize(cellind) = binsize*ones(1,ncell);
   pinfo.parm.cl0StartThres(cellind) = sthr*ones(1,ncell);
   pinfo.parm.cl0PeakThres(cellind) = pthr*ones(1,ncell);
   pinfo.parm.cl0MaxGap(cellind) = mgap*ones(1,ncell);
   pinfo.parm.cl0MinDur(cellind) = mindur*ones(1,ncell);
   pinfo.parm.cl0MaxDur(cellind) = maxdur*ones(1,ncell);
end

function pinfo = assignspindleparm(pinfo, cellind)
ok = 1; p{1} = 'Normalization mode (self/run/none)?'; d{1} = 'self'; II = inputdlg(p, 'Select a normalization mode', 1, d);
if strncmpi(II{1}, 'self', 2)
    norm = 'self';
elseif strncmpi(II{1}, 'run', 2)
    norm = 'run';
elseif strncmpi(II{1}, 'none', 2)
    norm = 'none';
elseif (~isempty(II))
    disp('-------> unknown mode, set to default (none)');
    norm = 'none';
else
    ok = 0;
end
if (ok)
pp = {'Enter peak threshold (std/mV):'; 'Enter start threshold (std/mV):'; 'Enter max gap (s):';...
      'Enter min duration (s):'; 'Enter max duration (s):'};
if (strcmp(norm, 'self'))
    def = {'6'; '4'; '0.15'; '0.3'; '10'};
elseif (strcmp(norm, 'run'))
    def = {'6'; '4'; '0.15'; '0.3'; '10'};
else
    def = {'0.3'; '0.2'; '0.15'; '0.3'; '10'};
end
III=inputdlg(pp, 'Parameters for computing spindle properties', 4, def, 'on');
if (~isempty(III))
   pthr = str2num(III{1}); sthr = str2num(III{2}); mgap = str2num(III{3}); mindur = str2num(III{4}); maxdur = str2num(III{5}); 
   nspike = numel(pinfo.general.eegfile);
   if (~isfield(pinfo.parm, 'spinNorm')) pinfo.parm.spinNorm = cell(1, nspike); end
   if (~isfield(pinfo.parm, 'spinPeakThres')) pinfo.parm.spinPeakThres = NaN*ones(1, nspike); end
   if (~isfield(pinfo.parm, 'spinStartThres')) pinfo.parm.spinStartThres = NaN*ones(1, nspike); end
   if (~isfield(pinfo.parm, 'spinMaxGap')) pinfo.parm.spinMaxGap = NaN*ones(1, nspike); end
   if (~isfield(pinfo.parm, 'spinMinDur')) pinfo.parm.spinMinDur = NaN*ones(1, nspike); end
   if (~isfield(pinfo.parm, 'spinMaxDur')) pinfo.parm.spinMaxDur = NaN*ones(1, nspike); end
   ncell = numel(cellind);
   for (i = 1:ncell) pinfo.parm.spinNorm{cellind(i)} = norm; end
   pinfo.parm.spinPeakThres(cellind) = pthr*ones(1,ncell);
   pinfo.parm.spinStartThres(cellind) = sthr*ones(1,ncell);
   pinfo.parm.spinMaxGap(cellind) = mgap*ones(1,ncell);
   pinfo.parm.spinMinDur(cellind) = mindur*ones(1,ncell); pinfo.parm.spinMaxDur(cellind) = maxdur*ones(1,ncell);
end
end

function pinfo = assignhvsparm(pinfo, cellind)
ok = 1; p{1} = 'Normalization mode (self/run/none)?'; d{1} = 'self'; II = inputdlg(p, 'Select a normalization mode', 1, d);
if strncmpi(II{1}, 'self', 2)
    norm = 'self';
elseif strncmpi(II{1}, 'run', 2)
    norm = 'run';
elseif strncmpi(II{1}, 'none', 2)
    norm = 'none';
elseif (~isempty(II))
    disp('-------> unknown mode, set to default (none)');
    norm = 'none';
else
    ok = 0;
end
if (ok)
pp = {'Enter peak threshold (std/mV):'; 'Enter start threshold (std/mV):'; 'Enter max gap (s):';...
      'Enter min number of peaks:'; 'Enter max duration (s):'; 'Peak mode? (negative/positive/absolute):'};
if (strcmp(norm, 'self'))
    def = {'6'; '2.5'; '0.5'; '5'; '100'; 'negative'};
elseif (strcmp(norm, 'run'))
    def = {'6'; '2.5'; '0.5'; '5'; '100'; 'negative'};
else
    def = {'0.3'; '0.12'; '0.5'; '5'; '100', 'negative'};
end
III=inputdlg(pp, 'Parameters for computing HVS properties', 6, def, 'on');
if (~isempty(III))
   pthr = str2num(III{1}); sthr = str2num(III{2}); mgap = str2num(III{3}); mindur = str2num(III{4}); maxdur = str2num(III{5}); 
   nspike = numel(pinfo.general.eegfile); peakmode = III{6};
   if (~isfield(pinfo.parm, 'hvsNorm')) pinfo.parm.hvsNorm = cell(1, nspike); end
   if (~isfield(pinfo.parm, 'hvsPeakThres')) pinfo.parm.hvsPeakThres = NaN*ones(1, nspike); end
   if (~isfield(pinfo.parm, 'hvsStartThres')) pinfo.parm.hvsStartThres = NaN*ones(1, nspike); end
   if (~isfield(pinfo.parm, 'hvsMaxGap')) pinfo.parm.hvsMaxGap = NaN*ones(1, nspike); end
   if (~isfield(pinfo.parm, 'hvsMinPeakNum')) pinfo.parm.hvsMinPeakNum = NaN*ones(1, nspike); end
   if (~isfield(pinfo.parm, 'hvsMaxDur')) pinfo.parm.hvsMaxDur = NaN*ones(1, nspike); end
   if (~isfield(pinfo.parm, 'hvsPeakMode')) pinfo.parm.hvsPeakMode = cell(1, nspike); end
   ncell = numel(cellind);
   for (i = 1:ncell) pinfo.parm.hvsNorm{cellind(i)} = norm; pinfo.parm.hvsPeakMode{cellind(i)} = peakmode; end
   pinfo.parm.hvsPeakThres(cellind) = pthr*ones(1,ncell);
   pinfo.parm.hvsStartThres(cellind) = sthr*ones(1,ncell);
   pinfo.parm.hvsMaxGap(cellind) = mgap*ones(1,ncell); 
   pinfo.parm.hvsMinPeakNum(cellind) = mindur*ones(1,ncell); pinfo.parm.hvsMaxDur(cellind) = maxdur*ones(1,ncell);
end
end

function pinfo = assignthetaparm(pinfo, cellind)
disp('-----------> No parameters are needed for this analysis');

function pinfo = assignspectralparm(pinfo, cellind)
pp = {'[Min Max] frequency (Hz)'; 'Frequency step (Hz):'; 'Window length (s; enter 0 if no spectrogram)';...
    'Window shift time (s)'; 'Normalization frequency range (Hz)'}; 
def = {'0.5 400'; '0.5'; '5'; '1'; '2 400'};
III=inputdlg(pp, 'PSD parameters', 5, def, 'on');
if (~isempty(III))
   ff = str2num(III{1}); minf = min(ff); maxf = max(ff); 
   fstep = str2num(III{2}); ws = str2num(III{3}); wshift = str2num(III{4});
   Nff = str2num(III{5}); normminf = min(Nff); normmaxf = max(Nff);
   pp2 = {'Event selection keyword'; 'Event selection type'};
   def2 = {'Track'; 'run'};
   III2=inputdlg(pp2, 'Spectrogram parameters', 2, def2, 'on');
   if ~isempty(III2)
      evkeyword = III2{1}; evtype = III2{2};
      nspike = numel(pinfo.general.eegfile);
      if (~isfield(pinfo.parm, 'specMinFreq')) pinfo.parm.specMinFreq = NaN*ones(1, nspike); end
      if (~isfield(pinfo.parm, 'specMaxFreq')) pinfo.parm.specMaxFreq = NaN*ones(1, nspike); end
      if (~isfield(pinfo.parm, 'specFreqStep')) pinfo.parm.specFreqStep = NaN*ones(1, nspike); end
      if (~isfield(pinfo.parm, 'specWinSize')) pinfo.parm.specWinSize = NaN*ones(1, nspike); end
      if (~isfield(pinfo.parm, 'specWinShift')) pinfo.parm.specWinShift = NaN*ones(1, nspike); end
      if (~isfield(pinfo.parm, 'specNormMinFreq')) pinfo.parm.specNormMinFreq = NaN*ones(1, nspike); end
      if (~isfield(pinfo.parm, 'specNormMaxFreq')) pinfo.parm.specNormMaxFreq = NaN*ones(1, nspike); end
      if (~isfield(pinfo.parm, 'specEvtKeyword')) pinfo.parm.specEvtKeyword = cell(1, nspike); end
      if (~isfield(pinfo.parm, 'specEvtType')) pinfo.parm.specEvtType = cell(1, nspike); end
      ncell = numel(cellind);
      pinfo.parm.specMinFreq(cellind) = minf*ones(1,ncell);
      pinfo.parm.specMaxFreq(cellind) = maxf*ones(1,ncell);
      pinfo.parm.specFreqStep(cellind) = fstep*ones(1,ncell);
      pinfo.parm.specWinSize(cellind) = ws*ones(1,ncell);
      pinfo.parm.specWinShift(cellind) = wshift*ones(1,ncell);
      pinfo.parm.specNormMinFreq(cellind) = normminf*ones(1,ncell);
      pinfo.parm.specNormMaxFreq(cellind) = normmaxf*ones(1,ncell);
      for (i = 1:numel(cellind))
          pinfo.parm.specEvtKeyword{cellind(i)} = evkeyword; pinfo.parm.specEvtType{cellind(i)} = evtype;
      end
   end
end

function pinfo = assignrippleparm(pinfo, cellind)
ok = 1; p{1} = 'Normalization mode (self/run/none)?'; d{1} = 'run'; II = inputdlg(p, 'Select a normalization mode', 1, d);
if strncmpi(II{1}, 'self', 2)
    norm = 'self';
elseif strncmpi(II{1}, 'run', 2)
    norm = 'run';
elseif strncmpi(II{1}, 'none', 2)
    norm = 'none';
elseif (~isempty(II))
    disp('-------> unknown mode, set to default (none)');
    norm = 'none';
else
    ok = 0;
end
if (ok)
pp = {'Enter [peak start] threshold (std/mV):'; 'Enter max gap (s):';...
      'Enter [min max] duration (s):'; 'Enter max amplitude (std):'; 'Peak mode? (positive/negative/absolute/both):'};
if (strcmp(norm, 'self'))
    def = {'6 2.5'; '0.03'; '0.02 0.4'; '20'; 'negative'};
elseif (strcmp(norm, 'run'))
    def = {'6 2.5'; '0.03'; '0.02 0.4'; '20'; 'negative'};
else
    def = {'0.15 0.07'; '0.03'; '0.02 0.4'; '0.5', 'negative'};
end
III=inputdlg(pp, 'Parameters for computing ripple properties', 5, def, 'on');
if (~isempty(III))
   thr = str2num(III{1}); pthr = thr(1); sthr = thr(2); mgap = str2num(III{2}); 
   dur = str2num(III{3}); mindur = dur(1); maxdur = dur(2); maxamp = str2num(III{4}); peakmode = III{5};
   nspike = numel(pinfo.general.eegfile);
   if (~isfield(pinfo.parm, 'rippNorm')) pinfo.parm.rippNorm = cell(1, nspike); end
   if (~isfield(pinfo.parm, 'rippPeakThres')) pinfo.parm.rippPeakThres = NaN*ones(1, nspike); end
   if (~isfield(pinfo.parm, 'rippStartThres')) pinfo.parm.rippStartThres = NaN*ones(1, nspike); end
   if (~isfield(pinfo.parm, 'rippMaxGap')) pinfo.parm.rippMaxGap = NaN*ones(1, nspike); end
   if (~isfield(pinfo.parm, 'rippMinDur')) pinfo.parm.rippMinDur = NaN*ones(1, nspike); end
   if (~isfield(pinfo.parm, 'rippMaxDur')) pinfo.parm.rippMaxDur = NaN*ones(1, nspike); end
   if (~isfield(pinfo.parm, 'rippMaxAmp')) pinfo.parm.rippMaxAmp = NaN*ones(1, nspike); end
   if (~isfield(pinfo.parm, 'rippPeakMode')) pinfo.parm.rippPeakMode = cell(1, nspike); end
   ncell = numel(cellind);
   for (i = 1:ncell) pinfo.parm.rippNorm{cellind(i)} = norm;  pinfo.parm.rippPeakMode{cellind(i)} = peakmode; end
   pinfo.parm.rippPeakThres(cellind) = pthr*ones(1,ncell);
   pinfo.parm.rippStartThres(cellind) = sthr*ones(1,ncell);
   pinfo.parm.rippMaxGap(cellind) = mgap*ones(1,ncell);
   pinfo.parm.rippMinDur(cellind) = mindur*ones(1,ncell); 
   pinfo.parm.rippMaxDur(cellind) = maxdur*ones(1,ncell);
   pinfo.parm.rippMaxAmp(cellind) = maxamp*ones(1,ncell);
end
end

function [writefilename, okk] = getoutputfile(hf, ow)
okk = 1; writefilename = [];
if (ow == 0)
   [fname, pname] = uiputfile(fullfile(cd, '*.eegdb'), 'Write the new spike database to:');
   if (numel(fname)>1)
      writefilename = fullfile(pname, fname);
   else
      okk = 0;
   end
else
   %input = questdlg('The current database will be altered and overwritten. Are you sure?', 'Overwrite?', 'Yes');
   %if (strcmp(input, 'Yes'))
      fname = get(hf, 'Name'); ftt = strfind(fname, '__'); writefilename = fname(ftt+2:numel(fname));
   %else
   %   okk = 0;
   %end
end

