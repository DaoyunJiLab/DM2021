function DataManager_ReCompute_Callback
%%Add or change variable values of the selected group, and save into a new database

hf = gcbf; pinfo = getappdata(hf, 'pinfo'); data = getappdata(hf, 'data'); tagmark = get(gcbo, 'Tag');
hgroup = getappdata(hf, 'hgroup');hfield = getappdata(hf, 'hfield');
ok = 1; plotparm = getappdata(hf, 'plotparm'); vv = plotparm.showdetail; %if vv=1, shoe details of computation
cc = plotparm.compute; %if cc=0, assign parameters for computation only; if 1, do the real computation
ow = plotparm.overwrite; %if ow=1, plot into the current database
[writefilename, ok] = getoutputfile(hf, ow);

%get selected cellind
groupselection = getappdata(hgroup, 'selection'); cellind = []; grpind = find(groupselection == 1); 
for (kk = 1:numel(grpind)) cellind = union(cellind, data.grouplist.groupindex{grpind(kk)}); end
if (~isempty(cellind)) && ok
    if (strcmp(tagmark, 'recomputeinit'))
        [pinfo,data] = DataManager_ComputeInit_Callback(pinfo,data, cellind, vv);
    elseif (strcmp(tagmark, 'recomputeall'))
        [pinfo,data] = DataManager_ComputeAll_Callback(pinfo,data, cellind, vv);
    elseif (strcmp(tagmark, 'firing'))
        if (cc==0) [pinfo, ok] = assignfiringparm(pinfo, cellind);  end
        if (cc == 1) [pinfo, data] = DataManager_FindSpikeFiring(pinfo, data, cellind, vv); end
    elseif (strcmp(tagmark, 'waveform'))
        if (cc==0) [pinfo, ok] = assignwaveformparm(pinfo, cellind); end
        if (cc ==1) [pinfo, data] = DataManager_FindWaveProp(pinfo, data, cellind, vv); end
    elseif (strcmp(tagmark, 'celltype'))
        if (cc==0) [pinfo, ok] = assigncelltypeparm(pinfo, cellind); end
        if (cc==1) [pinfo, data] = DataManager_FindNeuronType(pinfo, data, cellind, vv); end
    elseif (strcmp(tagmark, 'clusterqua'))
        if (cc==0) [pinfo, ok] = assignclusterquaparm(pinfo, cellind); end
        if (cc==1) [pinfo, data] = DataManager_FindClusterQuality(pinfo, data, cellind, vv); end
    elseif (strcmp(tagmark, 'pfprop'))
        if (cc ==0) [pinfo, ok] = assignfieldparm(pinfo, cellind); end 
        if (cc==1) 
            if plotparm.linkbehav == 1
                behav = getappdata(hf, 'behav'); bhdata = getappdata(hf, 'bhdata');
                [pinfo,data] = DataManager_FindSpatialInfoFieldProp(pinfo,data,behav,bhdata,cellind,vv);
            else
                msgbox('No behaviroal database linked; Fields not computed'); ok =0;
            end
        end
    elseif (strcmp(tagmark, 'pfdynam'))
        if (cc ==0) [pinfo, ok] = assignpfdynamparm(pinfo, cellind); end 
        if (cc==1) 
            if plotparm.linkbehav == 1
                behav = getappdata(hf, 'behav'); bhdata = getappdata(hf, 'bhdata');
                [pinfo,data] = DataManager_FindFieldDynam(pinfo,data,behav,bhdata,cellind,vv);
            else
                msgbox('No behaviroal database linked; Field dynamics cannot computed'); ok = 0;
            end
        end
    elseif (strcmp(tagmark, 'thetaphase'))
        if (cc ==0) [pinfo, ok] = assignthetaphaseparm(pinfo, cellind); end 
        if (cc==1) 
            if plotparm.linkeeg == 1  %%check eeg
                eeg = getappdata(hf, 'eeg'); eegdata = getappdata(hf, 'eegdata');
            else
                msgbox('No EEG database linked; Theta phases cannot computed'); ok = 0;
            end
            if ok %%%check behav
            if plotparm.linkbehav == 1
                behav = getappdata(hf, 'behav'); bhdata = getappdata(hf, 'bhdata');
                [pinfo,data] = DataManager_FindPhaseProperties(pinfo,data,eeg,eegdata,behav,bhdata,cellind);
            else
                bb = questdlg('No behaviroal database linked; No phase precession will be computed. Continue?'); 
                if strcmp(bb, 'Yes')
                    behav = []; bhdata = [];
                    [pinfo,data] = DataManager_FindPhaseProperties(pinfo,data,eeg,eegdata,behav,bhdata,cellind);
                else
                    ok = 0;
                end
            end
            end
        end
    end
 
    %%%%%%%%save and plot the new database
    if (~isempty(pinfo.general.parmfile))
        if (ok)
            s = whos('data'); Gb1 = s.bytes/(1024^3);
            s = whos('pinfo'); Gb2 = s.bytes/(1024^3);
            if Gb1 + Gb2 < 2
               save(writefilename, 'pinfo', 'data'); %disp('Small one');
            else
               save(writefilename, 'pinfo', 'data', '-mat', '-v7.3'); %disp('Big one');
            end
            %save(writefilename, 'pinfo', 'data');
            %save(writefilename, 'pinfo', 'data', '-mat', '-v7.3'); %%% This allows to save files bigger than 2GB on 64-bit systems
            if (ow)
               iii = get(hf, 'Children'); delete(iii(~strcmp(get(iii, 'Type'), 'uimenu'))); 
               fname = get(hf, 'Name'); ftt = strfind(fname, '__'); 
               newname = fname(1:ftt-1); set(hf, 'Name', newname); hmain = hf; 
            else
               hmain = DataManager_DataManager_Callback; plotparm.linkspike = 1; plotparm.linkeeg = 0; plotparm.linkbehav = 0;
            end
            setappdata(hmain,'plotparm', plotparm);
            DataManager_PlotSpikeDatabase(hmain, pinfo, data, 'Cell', 'clname', '.spikedb');
            set(hmain, 'Name', strcat(get(hmain, 'Name'), '__', writefilename));
            setappdata(hmain, 'pinfo', pinfo); setappdata(hmain, 'data', data); pinfo = []; data = [];
            %%%%%now reset all the selection options and update the display
            %   ----------do not do this.
        end
    else
        disp('-------------> no cells in the database!');
    end
else
    disp('--------------> no groups selected, groups do not contain any cells, or action cancelled');
end
disp('**********************');

function [pinfo, ok] = assignpfdynamparm(pinfo, cellind)
nspike = numel(pinfo.general.parmfile); ncell = numel(cellind); ok = 1;
%%%%%assign parameters for computing place field dynamics
pp = {'Session segment duration (s)'; 'Baseline session for 2D map dynamics (session name)?'; 'Baseline segment for 2D map dynamics (first, last, average, number(s))?';...
      'Baseline session for 1D map dynamics (session name)?'; 'Baseline lap for 1D map dynamics (first, last, average, number(s))?'};
def = {'300'; 'Open1'; 'last'; 'Track1'; 'last'};
III=inputdlg(pp, 'Parameters for computing field dynamics', 4, def, 'on');
if (~isempty(III))
   if (~isfield(pinfo.parm, 'fd2DSessSegTime')) pinfo.parm.fd2DSessSegDur = NaN*ones(1, nspike); end 
   if (~isfield(pinfo.parm, 'fd2DBaseSess')) pinfo.parm.fd2DBaseSess = cell(1, nspike); end
   if (~isfield(pinfo.parm, 'fd2DBaseSeg')) pinfo.parm.fd2DBaseSeg = cell(1, nspike); end
   if (~isfield(pinfo.parm, 'fd1DBaseSess')) pinfo.parm.fd1DBaseSess = cell(1, nspike); end
   if (~isfield(pinfo.parm, 'fd1DBaseLap')) pinfo.parm.fd1DBaseLap = cell(1, nspike); end
   for (i = 1:ncell) pinfo.parm.fd2DSessSegDur(cellind(i)) = str2num(III{1}); end
   for (i = 1:ncell) pinfo.parm.fd2DBaseSess{cellind(i)} = III{2}; end
   for (i = 1:ncell) pinfo.parm.fd2DBaseSeg{cellind(i)} = III{3}; end
   for (i = 1:ncell) pinfo.parm.fd1DBaseSess{cellind(i)} = III{4}; end
   for (i = 1:ncell) pinfo.parm.fd1DBaseLap{cellind(i)} = III{5}; end
else
   ok = 0;
end

function [pinfo, ok] = assignfieldparm(pinfo, cellind)
nspike = numel(pinfo.general.parmfile); ncell = numel(cellind); ok = 1;
%%%%%assign parameters for computing place field properties
pp = {'Smooth 2D rate maps?'; 'Enter 2D smoothing sigma (cm);'; 'Enter 2D N sigmas for smoothing:';...
    'Smooth 1D rate maps?'; 'Enter 1D smoothing sigma (cm):'; 'Enter 1D N sigmas for smoothing:'};
def = {'yes'; '5'; '5'; 'yes'; '5'; '5'};
III=inputdlg(pp, 'Parameters for computing rate maps', 6, def, 'on');
if (~isempty(III))
   s2d = III{1}; sig2d = str2num(III{2}); Nsig2d = str2num(III{3}); 
   if (~isfield(pinfo.parm, 'fSmooth2D')) pinfo.parm.fSmooth2D = cell(1, nspike); end
   if (~isfield(pinfo.parm, 'fSmooth2DSigma')) pinfo.parm.fSmooth2DSigma = NaN*ones(1, nspike); end
   if (~isfield(pinfo.parm, 'fSmooth2DNSigma')) pinfo.parm.fSmooth2DNSigma = NaN*ones(1, nspike); end
   for (i = 1:ncell) pinfo.parm.fSmooth2D{cellind(i)} = s2d; end
   pinfo.parm.fSmooth2DSigma(cellind) = sig2d*ones(1,ncell); pinfo.parm.fSmooth2DNSigma(cellind) = Nsig2d*ones(1,ncell); 
   s1d = III{4}; sig1d = str2num(III{5}); Nsig1d = str2num(III{6}); 
   if (~isfield(pinfo.parm, 'fSmooth1D')) pinfo.parm.fSmooth1D = cell(1, nspike); end
   if (~isfield(pinfo.parm, 'fSmooth1DSigma')) pinfo.parm.fSmooth1DSigma = NaN*ones(1, nspike); end
   if (~isfield(pinfo.parm, 'fSmooth1DNSigma')) pinfo.parm.fSmooth1DNSigma = NaN*ones(1, nspike); end
   for (i = 1:ncell) pinfo.parm.fSmooth1D{cellind(i)} = s1d; end
   pinfo.parm.fSmooth1DSigma(cellind) = sig1d*ones(1,ncell); pinfo.parm.fSmooth1DNSigma(cellind) = Nsig1d*ones(1,ncell); 
else
    ok = 0;
end
if (ok)
pp = {'Enter 1D grid N for computing baseline rate (x total grid N):'; 'Enter 1D minimum peak rate (Hz):'; ...
    'Enter 1D minimum in-field rate (x peak rate):'; 'Enter 1D maximum in-field gap (cm):'; 'Enter 1D minimum laps'};
def = {'0.3'; '3'; '0.1'; '5'; '4'};
III=inputdlg(pp, 'Parameters for defining 1D fields', 5, def, 'on');
if (~isempty(III))
   gridN = str2num(III{1}); minpeak = str2num(III{2}); thres = str2num(III{3}); gap = str2num(III{4}); minlap = str2num(III{5});
   if (~isfield(pinfo.parm, 'f1DBaseRateGridN')) pinfo.parm.f1DBaseRateGridN = NaN*ones(1, nspike); end
   if (~isfield(pinfo.parm, 'f1DMinPeakRate')) pinfo.parm.f1DMinPeakRate = NaN*ones(1, nspike); end
   if (~isfield(pinfo.parm, 'f1DThresholdRate')) pinfo.parm.f1DThresholdRate = NaN*ones(1, nspike); end
   if (~isfield(pinfo.parm, 'f1DMaxGap')) pinfo.parm.f1DMaxGap = NaN*ones(1, nspike); end
   if (~isfield(pinfo.parm, 'f1DMinLapNum')) pinfo.parm.f1DMinLapNum = NaN*ones(1, nspike); end
   pinfo.parm.f1DBaseRateGridN(cellind) = gridN*ones(1,ncell); pinfo.parm.f1DMinPeakRate(cellind) = minpeak*ones(1,ncell);
   pinfo.parm.f1DThresholdRate(cellind) = thres*ones(1,ncell); pinfo.parm.f1DMaxGap(cellind) = gap*ones(1,ncell); 
   pinfo.parm.f1DMinLapNum(cellind) = minlap*ones(1,ncell);
else
    ok = 0;
end
end
if (ok)
pp = {'Enter 2D grid N for computing baseline rate (x total grid N):'; 'Enter 2D minimum peak rate (Hz):'; ...
    'Enter 2D minimum in-field rate (x peak rate):'; 'Enter 2D maximum in-field gap (cm):'};
def = {'0.3'; '5'; '0.1'; '5'};
III=inputdlg(pp, 'Parameters for defining 2D fields', 4, def, 'on');
if (~isempty(III))
   gridN = str2num(III{1}); minpeak = str2num(III{2}); thres = str2num(III{3}); gap = str2num(III{4}); 
   if (~isfield(pinfo.parm, 'f2DBaseRateGridN')) pinfo.parm.f2DBaseRateGridN = NaN*ones(1, nspike); end
   if (~isfield(pinfo.parm, 'f2DMinPeakRate')) pinfo.parm.f2DMinPeakRate = NaN*ones(1, nspike); end
   if (~isfield(pinfo.parm, 'f2DThresholdRate')) pinfo.parm.f2DThresholdRate = NaN*ones(1, nspike); end
   if (~isfield(pinfo.parm, 'f2DMaxGap')) pinfo.parm.f2DMaxGap = NaN*ones(1, nspike); end
   pinfo.parm.f2DBaseRateGridN(cellind) = gridN*ones(1,ncell); pinfo.parm.f2DMinPeakRate(cellind) = minpeak*ones(1,ncell);
   pinfo.parm.f2DThresholdRate(cellind) = thres*ones(1,ncell); pinfo.parm.f2DMaxGap(cellind) = gap*ones(1,ncell); 
else
    ok = 0;
end
end
if (ok)
pp = {'Max diretional shifts (cm; -1 for no directionality)?'; 'Correlate rotated rate maps?'; 'Enter rotating center X (cm):';...
    'Entering rotating center Y (cm):'; 'Enter rorating angle (degree):'};
def = {'25'; 'no'; '100'; '45'; '90'};
III=inputdlg(pp, 'Parameters for computing directional/ratational properties', 5, def, 'on');
if (~isempty(III))
   rot = III{2}; dirpix = str2num(III{1}); xx = str2num(III{3}); yy = str2num(III{4}); ang = str2num(III{5}); 
   if (~isfield(pinfo.parm, 'fDirShift')) pinfo.parm.fDirShift = NaN*ones(1, nspike); end
   if (~isfield(pinfo.parm, 'fRotate')) pinfo.parm.fRotate = cell(1, nspike); end
   if (~isfield(pinfo.parm, 'fRotate2DCenX')) pinfo.parm.fRotate2DCenX = NaN*ones(1, nspike); end
   if (~isfield(pinfo.parm, 'fRotate2DCenY')) pinfo.parm.fRotate2DCenY = NaN*ones(1, nspike); end
   if (~isfield(pinfo.parm, 'fRotate2DAngle')) pinfo.parm.fRotate2DAngle = NaN*ones(1, nspike); end
   pinfo.parm.fDirShift(cellind) = dirpix*ones(1,ncell); for (i = 1:ncell) pinfo.parm.fRotate{cellind(i)} = rot; end
   pinfo.parm.fRotate2DCenX(cellind) = xx*ones(1,ncell); pinfo.parm.fRotate2DCenY(cellind) = yy*ones(1,ncell); 
   pinfo.parm.fRotate2DAngle(cellind) = ang*ones(1,ncell); 
else
    ok = 0;
end
end

function [pinfo, ok] = assignclusterquaparm(pinfo, cellind)
ok = 1;
pp = {'Enter minimum number of cluster points'; 'Enter minimum standard deviation of amplitudes (AD units):'}; %; 'Enter number of bins:'};
def = {'100'; '100'}; %; '100'};
III=inputdlg(pp, 'Parameters for assessing cluster quality', 2, def, 'on');
if (~isempty(III))
   minN = str2num(III{1}); mindev = str2num(III{2}); %numbin = str2num(III{3}); 
   nspike = numel(pinfo.general.parmfile);
   if (~isfield(pinfo.parm, 'clustMinNpoints')) pinfo.parm.clustMinNpoints = NaN*ones(1, nspike); end
   if (~isfield(pinfo.parm, 'clustMinDev')) pinfo.parm.clustMinDev = NaN*ones(1, nspike); end
   if (~isfield(pinfo.parm, 'clustNumBin')) pinfo.parm.clustNumBin = NaN*ones(1, nspike); end
   ncell = numel(cellind);
   pinfo.parm.clustMinNpoints(cellind) = minN*ones(1,ncell);
   pinfo.parm.clustMinDev(cellind) = mindev*ones(1,ncell);
   pinfo.parm.clustNumBin(cellind) = 0*ones(1,ncell);
else
    ok = 0;
end

function [pinfo, ok] = assignwaveformparm(pinfo, cellind)
pp = {'Enter number of points for waveform baseline'}; ok = 1;
def = {'0'};
III=inputdlg(pp, 'Parameters for computing waveform properties', 3, def, 'on');
if (~isempty(III))
   baseN = str2num(III{1});  
   nspike = numel(pinfo.general.parmfile);
   if (~isfield(pinfo.parm, 'spikebasepoint')) pinfo.parm.spikebasepoint = NaN*ones(1, nspike); end
   ncell = numel(cellind);
   pinfo.parm.spikebasepoint(cellind) = baseN*ones(1,ncell);
else
   ok = 0;
end

function [pinfo, ok] = assignfiringparm(pinfo, cellind)
pp = {'Enter binsize(s) for computing max rate'; 'Enter maximum spike burst interval(s):'; 'Enter spike refractory time(s):'};
def = {'1'; '0.01'; '0.002'}; ok = 1;
III=inputdlg(pp, 'Parameters for computing firing properties', 3, def, 'on');
if (~isempty(III))
   bin = str2num(III{1}); Bint = str2num(III{2}); Rint = str2num(III{3}); 
   nspike = numel(pinfo.general.parmfile);
   if (~isfield(pinfo.parm, 'maxratetimebin')) pinfo.parm.maxratetimebin = NaN*ones(1, nspike); end
   if (~isfield(pinfo.parm, 'spikeburstint')) pinfo.parm.spikeburstint = NaN*ones(1, nspike); end
   if (~isfield(pinfo.parm, 'spikerefracint')) pinfo.parm.spikerefracint = NaN*ones(1, nspike); end
   ncell = numel(cellind);
   pinfo.parm.maxratetimebin(cellind) = bin*ones(1,ncell);
   pinfo.parm.spikeburstint(cellind) = Bint*ones(1,ncell);
   pinfo.parm.spikerefracint(cellind) = Rint*ones(1,ncell);
else
    ok = 0;
end

function [pinfo, ok] = assignthetaphaseparm(pinfo, cellind)
pp = {'Reference EEG record area'; 'Reference EEG band'; 'Reference EEG keyword'};
def = {'CA1'; 'theta'; 'smooth'}; ok = 1;
III=inputdlg(pp, 'Parameters for reference EEG', 3, def, 'on');
if (~isempty(III))
   nspike = numel(pinfo.general.parmfile);
   if (~isfield(pinfo.parm, 'phaseEEGarea')) pinfo.parm.phaseEEGarea = cell(1, nspike); end
   if (~isfield(pinfo.parm, 'phaseEEGband')) pinfo.parm.phaseEEGband = cell(1, nspike); end
   if (~isfield(pinfo.parm, 'phaseEEGkeyword')) pinfo.parm.phaseEEGkeyword = cell(1, nspike); end
   ncell = numel(cellind);
   for (i = 1:ncell)
       pinfo.parm.phaseEEGarea{cellind(i)} = III{1};
       pinfo.parm.phaseEEGband{cellind(i)} = III{2};
       pinfo.parm.phaseEEGkeyword{cellind(i)} = III{3};
   end
else
   ok = 0;
end
if ok
pp = {'Session keyword'; 'Session keytype'};
def = {'Track'; 'linear'};
III=inputdlg(pp, 'Parameters for session selection', 2, def, 'on');
if (~isempty(III))
   if (~isfield(pinfo.parm, 'phaseSessKeyword')) pinfo.parm.phaseSessKeyword = cell(1, nspike); end
   if (~isfield(pinfo.parm, 'phaseSessKeytype')) pinfo.parm.phaseSessKeytype = cell(1, nspike); end
   for (i = 1:ncell)
       pinfo.parm.phaseSessKeyword{cellind(i)} = III{1};
       pinfo.parm.phaseSessKeytype{cellind(i)} = III{2};
   end
else
   ok = 0;
end
end
if ok
pp = {'Event keyword'; 'Event keytype'};
def = {'Track'; 'run'};
III=inputdlg(pp, 'Parameters for event selection', 2, def, 'on');
if (~isempty(III))
   if (~isfield(pinfo.parm, 'phaseEvtKeyword')) pinfo.parm.phaseEvtKeyword = cell(1, nspike); end
   if (~isfield(pinfo.parm, 'phaseEvtKeytype')) pinfo.parm.phaseEvtKeytype = cell(1, nspike); end
   for (i = 1:ncell)
       pinfo.parm.phaseEvtKeyword{cellind(i)} = III{1};
       pinfo.parm.phaseEvtKeytype{cellind(i)} = III{2};
   end
else
   ok = 0;
end
end
% if ok
% pp = {'Min slope'; 'Max slope'; 'Slope step'};
% def = {'-20'; '10'; '0.1'};
% III=inputdlg(pp, 'Parameters for slope search', 3, def, 'on');
% if (~isempty(III))
%    if (~isfield(pinfo.parm, 'phaseMinSlope')) pinfo.parm.phaseMinSlope = NaN*ones(1, nspike); end
%    if (~isfield(pinfo.parm, 'phaseMaxSlope')) pinfo.parm.phaseMaxSlope = NaN*ones(1, nspike); end
%    if (~isfield(pinfo.parm, 'phaseSlopeStep')) pinfo.parm.phaseSlopeStep = NaN*ones(1, nspike); end
%    pinfo.parm.phaseMinSlope(cellind) = str2num(III{1})*ones(1, ncell);
%    pinfo.parm.phaseMaxSlope(cellind) = str2num(III{2})*ones(1, ncell);
%    pinfo.parm.phaseSlopeStep(cellind) = str2num(III{3})*ones(1, ncell);
% else
%    ok = 0;
% end
% end
% if ok
% pp = {'Min offset'; 'Max offset'; 'Offset step'};
% def = {'-1000'; '1000'; '10'};
% III=inputdlg(pp, 'Parameters for offset search', 3, def, 'on');
% if (~isempty(III))
%    if (~isfield(pinfo.parm, 'phaseMinOffset')) pinfo.parm.phaseMinOffset = NaN*ones(1, nspike); end
%    if (~isfield(pinfo.parm, 'phaseMaxOffset')) pinfo.parm.phaseMaxOffset = NaN*ones(1, nspike); end
%    if (~isfield(pinfo.parm, 'phaseOffsetStep')) pinfo.parm.phaseOffsetStep = NaN*ones(1, nspike); end
%    pinfo.parm.phaseMinOffset(cellind) = str2num(III{1})*ones(1, ncell);
%    pinfo.parm.phaseMaxOffset(cellind) = str2num(III{2})*ones(1, ncell);
%    pinfo.parm.phaseOffsetStep(cellind) = str2num(III{3})*ones(1, ncell);
% else
%    ok = 0;
% end
% end

function [pinfo, ok] = assigncelltypeparm(pinfo, cellind)
%%%%first get recarea and determine there are hp and/or ctx cells
recarea = pinfo.general.recarea; nspike = numel(recarea); ok = 1;
hpcellind = find( (strncmpi(recarea,'CA',2)) | (strncmpi(recarea,'HP',2)) | (strncmpi(recarea,'DG',2)) | (strncmpi(recarea,'hip',3)) );
hpcellind = intersect(hpcellind, cellind);
ctxcellind = find( (strncmpi(recarea,'ctx',3)) | (strncmpi(recarea,'V',1)) | (strncmpi(recarea,'PFC',3)) | (strncmpi(recarea,'PC',1)) );
ctxcellind = intersect(ctxcellind, cellind);
if (~isempty(hpcellind))
    pp = {'Enter HP runmeanrate(Hz) threshold:'; 'Enter HP maxnpampratio threshold:'; 'Enter HP maxhlwd(s) threshold:'};
    def = {'8'; '0.72'; '0.0002'};
    III=inputdlg(pp, 'Parameters for classifying neuron types', 3, def, 'on');
    if (~isempty(III))
        hpTmeanrate = str2num(III{1}); hpTratio = str2num(III{2}); hpTwd = str2num(III{3});
    else
        ok = 0;
    end
end
if (ok) & (~isempty(ctxcellind))
    pp = {'Enter CTX runmeanrate(Hz) threshold:';'Enter CTX maxnpampratio threshold:'; 'Enter CTX maxhlwd(s) threshold:'};
    def = {'10'; '0.72'; '0.0002'};
    III=inputdlg(pp, 'Input for CTX classification parmameters', 3, def, 'on');
    if (~isempty(III))
        ctxTmeanrate = str2num(III{1}); ctxTratio = str2num(III{2}); ctxTwd = str2num(III{3});
    else
        ok = 0;
    end
end
if (ok)
    if (~isfield(pinfo.parm, 'CelltypeTmeanrate')) pinfo.parm.CelltypeTmeanrate = NaN*ones(1,nspike); end 
    if (~isfield(pinfo.parm, 'CelltypeThlwd')) pinfo.parm.CelltypeThlwd = NaN*ones(1,nspike); end
    if (~isfield(pinfo.parm, 'CelltypeTampratio')) pinfo.parm.CelltypeTampratio = NaN*ones(1,nspike); end
end
if (ok) & (~isempty(hpcellind))
    if (~isempty(hpTmeanrate))
        pinfo.parm.CelltypeTmeanrate(hpcellind) = hpTmeanrate*ones(1,numel(hpcellind));
    end
    if (~isempty(hpTwd))
        pinfo.parm.CelltypeThlwd(hpcellind) = hpTwd*ones(1,numel(hpcellind));
    end
    if (~isempty(hpTratio))
        pinfo.parm.CelltypeTampratio(hpcellind) = hpTratio*ones(1,numel(hpcellind));
    end
end
if (ok) & (~isempty(ctxcellind))
    if (~isempty(ctxTmeanrate))
        pinfo.parm.CelltypeTmeanrate(ctxcellind) = ctxTmeanrate*ones(1,numel(ctxcellind));
    end
    if (~isempty(ctxTwd))
        pinfo.parm.CelltypeThlwd(ctxcellind) = ctxTwd*ones(1,numel(ctxcellind));
    end
    if (~isempty(ctxTratio))
        pinfo.parm.CelltypeTampratio(ctxcellind) = ctxTratio*ones(1,numel(ctxcellind));
    end
end

function [writefilename, okk] = getoutputfile(hf, ow)
okk = 1; writefilename = [];
if (ow == 0)
   [fname, pname] = uiputfile(fullfile(cd, '*.spikedb'), 'Write the new spike database to:');
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


