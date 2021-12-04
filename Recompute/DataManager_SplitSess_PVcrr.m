function DataManager_SplitSess_PVcrr
%Compute population vector: - split laps in a session into
%subgroups (early/late) and then compute PV crr among groups

%For (open field) session: original PFDynam can be used to do the same;
%just need to set the parameter behav.parm.sessSegTime

hf = gcbf; pinfo = getappdata(hf, 'pinfo'); data = getappdata(hf, 'data'); %tagmark = get(gcbo, 'Tag');
hgroup = getappdata(hf, 'hgroup'); %hfield = getappdata(hf, 'hfield');
ok = 1; oc = 0; plotparm = getappdata(hf, 'plotparm'); 
disp('Split track sessions to groups of laps and re-compute spatial variables...');
%get selected cellind
groupselection = getappdata(hgroup, 'selection'); cellind = []; grpind = find(groupselection == 1); 
for (kk = 1:numel(grpind)) cellind = union(cellind, data.grouplist.groupindex{grpind(kk)}); end
if isempty(cellind) 
    ok = 0; disp('-----> no cells selected; aborted'); 
end 
if ok
    if plotparm.linkbehav == 1
       behav = getappdata(hf, 'behav'); bhdata = getappdata(hf, 'bhdata');
    else
       disp('-----> no behaviroal database linked; aborted'); ok = 0;
    end
end
if ok
    if ~isfield(pinfo, 'PVcrr')
        disp('-----> PVcrr variables need to be computed first; aborted'); ok = 0;
    end
end
if ok
    if (~isfield(pinfo, 'fielddynam')) % || (~isfield(pinfo.fielddynam, 'lapgrpPairname'))
        disp('-----> field dynamics need to be performed first; aborted'); ok = 0;
    end
end
if ok
    input = inputdlg({'Split by time/lap number'; 'Duration/lap number'; 'Min time/lap number'}, 'Split method', 3, {'lap'; '5'; '5'}); 
    if (~isempty(input))
        parm.splitmethod = input{1}; parm.grouplength = str2num(input{2}); parm.grpminlapnum = str2num(input{3}); 
    else
        ok = 0;
    end
end 
if ok
   ow = plotparm.overwrite; %if ow=1, plot into the current database
   [writefilename, ok] = getoutputfile(hf, ow);
end
if ok
   [pinfo, data] = splitsession_spPVcrr(pinfo, data, behav, bhdata, cellind, parm); oc = 1;
end
if ok
    if oc
        save(writefilename, 'pinfo', 'data', '-mat', '-v7.3');
        if (ow)
            iii = get(hf, 'Children'); delete(iii(~strcmp(get(iii, 'Type'), 'uimenu'))); 
            fname = get(hf, 'Name'); ftt = strfind(fname, '__'); 
            newname = fname(1:ftt-1); set(hf, 'Name', newname); hmain = hf; 
        else
            hmain = DataManager_DataManager_Callback; plotparm.linkspike = 1; plotparm.linkeeg = 0; plotparm.linkbehav = 0;
        end
        setappdata(hmain, 'plotparm', plotparm);
        DataManager_PlotSpikeDatabase(hmain, pinfo, data, 'Cell', 'clname', '.spikedb');
        set(hmain, 'Name', strcat(get(hmain, 'Name'), '__', writefilename));
        setappdata(hmain, 'pinfo', pinfo); setappdata(hmain, 'data', data); pinfo = []; data = [];
    else
        disp('-----> database not altered; no new data to save');
    end 
end
disp('****************');

function [pinfo, data] = splitsession_spPVcrr(pinfo, data, behav, bhdata, cellind, parm)
%%%%%% Spatial field dynamics
if (~isempty(cellind))
   nspike = numel(pinfo.general.parmfile);
   %%%%assign parm variables
   if (~isfield(pinfo.parm, 'SplitMethod')) pinfo.parm.pvSplitMethod = cell(1, nspike); end   
   if (~isfield(pinfo.parm, 'lapgrpLength')) pinfo.parm.pvlapgrpLength = cell(1, nspike); end
   if (~isfield(pinfo.parm, 'lapgrpMinNumber')) pinfo.parm.pvlapgrpMinNumber = cell(1, nspike); end
   for (ttt = 1:numel(cellind))
       pinfo.parm.pvSplitMethod{cellind(ttt)} = parm.splitmethod;
       pinfo.parm.pvlapgrpLength{cellind(ttt)} = parm.grouplength; 
       pinfo.parm.pvlapgrpMinNumber{cellind(ttt)} = parm.grpminlapnum;
   end
   %%% define session spPVcrr variables
   if (~isfield(pinfo, 'spPVcrr')) pinfo.spPVcrr = []; end 
   if (~isfield(pinfo.spPVcrr, 'sessPairName')) pinfo.spPVcrr.sessPairName = cell(1, nspike); end
   if (~isfield(pinfo.spPVcrr, 'sessPairType')) pinfo.spPVcrr.sessPairType = cell(1, nspike); end
   if (~isfield(pinfo.spPVcrr, 'sessDgnlMeanCrr')) pinfo.spPVcrr.sessDgnlMeanCrr = cell(1, nspike); end %%%real values only assigned to the first cell of a day; others set as NaN
   if (~isfield(pinfo.spPVcrr, 'sessDgnlMedCrr')) pinfo.spPVcrr.sessDgnlMedCrr = cell(1, nspike); end
   if (~isfield(pinfo.spPVcrr, 'sessDgnlVarCrr')) pinfo.spPVcrr.sessDgnlVarCrr = cell(1, nspike); end
   %%% define event spPVcrr variables - along real diagonal line 
   if (~isfield(pinfo.spPVcrr, 'evtPairName')) pinfo.spPVcrr.evtPairName = cell(1, nspike); end
    if (~isfield(pinfo.spPVcrr, 'evtPairType')) pinfo.spPVcrr.evtPairType = cell(1, nspike); end
   if (~isfield(pinfo.spPVcrr, 'evtMeanPeakCrr')) pinfo.spPVcrr.evtMeanPeakCrr = cell(1, nspike); end
   if (~isfield(pinfo.spPVcrr, 'evtMedPeakCrr')) pinfo.spPVcrr.evtMedPeakCrr = cell(1, nspike); end
   if (~isfield(pinfo.spPVcrr, 'evtVarPeakCrr')) pinfo.spPVcrr.evtVarPeakCrr = cell(1, nspike); end 
   if (~isfield(pinfo.spPVcrr, 'evtMeanPeakLoc')) pinfo.spPVcrr.evtMeanPeakLoc = cell(1, nspike); end 
   if (~isfield(pinfo.spPVcrr, 'evtMedPeakLoc')) pinfo.spPVcrr.evtMedPeakLoc = cell(1, nspike); end
   if (~isfield(pinfo.spPVcrr, 'evtVarPeakLoc')) pinfo.spPVcrr.evtVarPeakLoc = cell(1, nspike); end 
   if (~isfield(pinfo.spPVcrr, 'evtMeanDgnlCrr')) pinfo.spPVcrr.evtMeanDgnlCrr = cell(1, nspike); end
   if (~isfield(pinfo.spPVcrr, 'evtMedDgnlCrr')) pinfo.spPVcrr.evtMedDgnlCrr = cell(1, nspike); end
   if (~isfield(pinfo.spPVcrr, 'evtVarDgnlCrr')) pinfo.spPVcrr.evtVarDgnlCrr = cell(1, nspike); end  
   if (~isfield(pinfo.spPVcrr, 'evtMeanNgbrCrr')) pinfo.spPVcrr.evtMeanNgbrCrr = cell(1, nspike); end 
   if (~isfield(pinfo.spPVcrr, 'evtMedNgbrCrr')) pinfo.spPVcrr.evtMedNgbrCrr = cell(1, nspike); end
   if (~isfield(pinfo.spPVcrr, 'evtVarNgbrCrr')) pinfo.spPVcrr.evtVarNgbrCrr = cell(1, nspike); end     
   if (~isfield(pinfo.spPVcrr, 'evtMeanDistCrr')) pinfo.spPVcrr.evtMeanDistCrr = cell(1, nspike); end 
   if (~isfield(pinfo.spPVcrr, 'evtMedDistCrr')) pinfo.spPVcrr.evtMedDistCrr = cell(1, nspike); end
   if (~isfield(pinfo.spPVcrr, 'evtVarDistCrr')) pinfo.spPVcrr.evtVarDistCrr = cell(1, nspike); end  
   %%% define event spPVcrr variables - along flipped diagonal line   
   if (~isfield(pinfo.spPVcrr, 'evtMeanFPeakCrr')) pinfo.spPVcrr.evtMeanFPeakCrr = cell(1, nspike); end 
   if (~isfield(pinfo.spPVcrr, 'evtMedFPeakCrr')) pinfo.spPVcrr.evtMedFPeakCrr = cell(1, nspike); end
   if (~isfield(pinfo.spPVcrr, 'evtVarFPeakCrr')) pinfo.spPVcrr.evtVarFPeakCrr = cell(1, nspike); end 
   if (~isfield(pinfo.spPVcrr, 'evtMeanFPeakLoc')) pinfo.spPVcrr.evtMeanFPeakLoc = cell(1, nspike); end 
   if (~isfield(pinfo.spPVcrr, 'evtMedFPeakLoc')) pinfo.spPVcrr.evtMedFPeakLoc = cell(1, nspike); end
   if (~isfield(pinfo.spPVcrr, 'evtVarFPeakLoc')) pinfo.spPVcrr.evtVarFPeakLoc = cell(1, nspike); end 
   if (~isfield(pinfo.spPVcrr, 'evtMeanFDgnlCrr')) pinfo.spPVcrr.evtMeanFDgnlCrr = cell(1, nspike); end 
   if (~isfield(pinfo.spPVcrr, 'evtMedFDgnlCrr')) pinfo.spPVcrr.evtMedFDgnlCrr = cell(1, nspike); end
   if (~isfield(pinfo.spPVcrr, 'evtVarFDgnlCrr')) pinfo.spPVcrr.evtVarFDgnlCrr = cell(1, nspike); end  
   if (~isfield(pinfo.spPVcrr, 'evtMeanFNgbrCrr')) pinfo.spPVcrr.evtMeanFNgbrCrr = cell(1, nspike); end 
   if (~isfield(pinfo.spPVcrr, 'evtMedFNgbrCrr')) pinfo.spPVcrr.evtMedFNgbrCrr = cell(1, nspike); end
   if (~isfield(pinfo.spPVcrr, 'evtVarFNgbrCrr')) pinfo.spPVcrr.evtVarFNgbrCrr = cell(1, nspike); end     
   if (~isfield(pinfo.spPVcrr, 'evtMeanFDistCrr')) pinfo.spPVcrr.evtMeanFDistCrr = cell(1, nspike); end 
   if (~isfield(pinfo.spPVcrr, 'evtMedFDistCrr')) pinfo.spPVcrr.evtMedFDistCrr = cell(1, nspike); end
   if (~isfield(pinfo.spPVcrr, 'evtVarFDistCrr')) pinfo.spPVcrr.evtVarFDistCrr = cell(1, nspike); end  
   data.spPVcrr.sessXYbin = cell(1, nspike);  
   data.spPVcrr.sessCrr = cell(1, nspike); data.spPVcrr.sessCrrP = cell(1, nspike); 
   data.spPVcrr.evtXbin = cell(1, nspike); data.spPVcrr.evtXavgjoint = cell(1, nspike); 
   data.spPVcrr.evtCrr = cell(1, nspike); data.spPVcrr.evtCrrP = cell(1, nspike);
   %%%%read in spPVcrr parameters
   parm.peakrange(1) = pinfo.parm.pvPeakRangeMin{cellind(1)}; parm.peakrange(2) = pinfo.parm.pvPeakRangeMax{cellind(1)};       
   parm.ngbrrange(1) = pinfo.parm.pvNgbrRangeMin{cellind(1)}; parm.ngbrrange(2) = pinfo.parm.pvNgbrRangeMax{cellind(1)};
   parm.distrange(1) = pinfo.parm.pvDistRangeMin{cellind(1)}; parm.distrange(2) = pinfo.parm.pvDistRangeMax{cellind(1)}; 
   parm.mincellnum = pinfo.parm.pvMinCellnum{cellind(1)};
   parm.evkeyword = pinfo.parm.pvevKeyword{cellind(1)}; parm.evkeytype = pinfo.parm.pvevKeytype{cellind(1)};
   parm.evkeynoword = pinfo.parm.pvevKeynoword{cellind(1)}; parm.evkeynotype = pinfo.parm.pvevKeynotype{cellind(1)};  
end
%%%%%determine how many datasets to compute
allfinaldir = pinfo.general.finaldir(cellind);
uniqfinaldir = unique(allfinaldir);
for (jjjk = 1:numel(uniqfinaldir))
    cellindnow = cellind(find(strcmp(allfinaldir, uniqfinaldir{jjjk}))); %%%all cellind in this dataset
if numel(cellindnow) >= parm.mincellnum
    i = cellindnow(1); %%%only assign variables for the first cell of the group
    disp(['-----> compute split-session PV crr in ', pinfo.general.finaldir{i}]);
    [XYbin, sessvector, sessoccutime, sessname, Xbin, evtvector, evtoccutime, evtname, fjoint] = getPVvectors(pinfo, data, behav, bhdata, cellindnow, parm);
       %%%XYbin{1}(x), XYbin{2}(y), sessvector{nsess}(ncell, y, x); sessoccutime{tt}(y,x)
       %%%Xbin(x), evtvector{nev}(ncell, x), evtname{nev}, fjoint(xjoint)
    %%%%%%%Now: pair up and compute - seesions - only compute crr's for diagonal locations
    nsess = numel(sessname); ip = nsess*(nsess-1)/2; ipair =0;
    pinfo.spPVcrr.sessPairName{i} = cell(1, ip); data.spPVcrr.sessXYbin{i} = XYbin; 
    data.spPVcrr.sessCrr{i} = cell(1, ip); data.spPVcrr.sessCrrP{i} = cell(i, ip);
    pinfo.spPVcrr.sessDgnlMeanCrr{i} = NaN*ones(1,ip); pinfo.spPVcrr.sessDgnlMedCrr{i} = NaN*ones(1,ip); pinfo.spPVcrr.sessDgnlVarCrr{i} = NaN*ones(1,ip); 
    for (m = 1:nsess)
        for (n = m+1:nsess)
            ipair = ipair + 1;
            pinfo.spPVcrr.sessPairName{i}{ipair} = strcat(sessname{m}, '+', sessname{n});
            pinfo.spPVcrr.sessPairType{i}{ipair} = findevtpairtype(sessname{m}, sessname{n});
            [sess2DCrr, sess2DP] = find2DspPVcrr(XYbin, sessvector{m}, sessvector{n});
            data.spPVcrr.sessCrr{i}{ipair} = sess2DCrr; data.spPVcrr.sessCrrP{i}{ipair} = sess2DP;
            aa = reshape(sess2DCrr, [1, numel(XYbin{1})*numel(XYbin{2})]);
            aa = aa(~isnan(aa)); pinfo.spPVcrr.sessDgnlMeanCrr{i}(ipair) = mean(aa); 
            pinfo.spPVcrr.sessDgnlMedCrr{i}(ipair) = median(aa); pinfo.spPVcrr.sessDgnlVarCrr{i}(ipair) = var(aa);
        end
    end
    %%%%%%% pair up and compute - events - compute crr's for all spatial locations
    nev = numel(evtname); ip = nev + nev*(nev-1)/2; ipair =0;
    data.spPVcrr.evtXbin{i} = Xbin; data.spPVcrr.evtXavgjoint{i} = fjoint;
    pinfo.spPVcrr.evtPairName{i} = cell(1, ip); data.spPVcrr.evtCrr{i} = cell(1, ip); data.spPVcrr.evtCrrP{i} = cell(i, ip);
    for (m = 1:numel(evtname))
        for (n = m:numel(evtname))
            ipair = ipair + 1;
            pinfo.spPVcrr.evtPairName{i}{ipair} = strcat(evtname{m}, '+', evtname{n});
            pinfo.spPVcrr.evtPairType{i}{ipair} = findevtpairtype(evtname{m}, evtname{n});
            [evt1DCrr, evt1DP] = find1DnowspPVcrr(Xbin, evtvector{m}, evtvector{n});
            data.spPVcrr.evtCrr{i}{ipair} = evt1DCrr; data.spPVcrr.evtCrrP{i}{ipair} = evt1DP;
        end
    end
    %            evt measures: peak, diagonal, neighbor, distal
    pinfo.spPVcrr.evtMeanPeakCrr{i} = NaN*ones(1,ip); pinfo.spPVcrr.evtMedPeakCrr{i} = NaN*ones(1,ip); pinfo.spPVcrr.evtVarPeakCrr{i} = NaN*ones(1,ip);
    pinfo.spPVcrr.evtMeanPeakLoc{i} = NaN*ones(1,ip); pinfo.spPVcrr.evtMedPeakLoc{i} = NaN*ones(1,ip); pinfo.spPVcrr.evtVarPeakLoc{i} = NaN*ones(1,ip);
    pinfo.spPVcrr.evtMeanDgnlCrr{i} = NaN*ones(1,ip); pinfo.spPVcrr.evtMedDgnlCrr{i} = NaN*ones(1,ip); pinfo.spPVcrr.evtVarDgnlCrr{i} = NaN*ones(1,ip);
    pinfo.spPVcrr.evtMeanNgbrCrr{i} = NaN*ones(1,ip); pinfo.spPVcrr.evtMedNgbrCrr{i} = NaN*ones(1,ip); pinfo.spPVcrr.evtVarNgbrCrr{i} = NaN*ones(1,ip);
    pinfo.spPVcrr.evtMeanDistCrr{i} = NaN*ones(1,ip); pinfo.spPVcrr.evtMedDistCrr{i} = NaN*ones(1,ip); pinfo.spPVcrr.evtVarDistCrr{i} = NaN*ones(1,ip);
    for m = 1:ip
       crrnow = data.spPVcrr.evtCrr{i}{m}; [peakV, peakloc, DiaV, NeiV, DisV] = findtargetvaluesincrr(crrnow, Xbin, parm);
       peakV = peakV(~isnan(peakV)); peakloc = peakloc(~isnan(peakloc));
       pinfo.spPVcrr.evtMeanPeakCrr{i}(m) = mean(peakV); pinfo.spPVcrr.evtMedPeakCrr{i}(m) = median(peakV); pinfo.spPVcrr.evtVarPeakCrr{i}(m) = var(peakV);
       pinfo.spPVcrr.evtMeanPeakLoc{i}(m) = mean(peakloc); pinfo.spPVcrr.evtMedPeakLoc{i}(m) = median(peakloc); pinfo.spPVcrr.evtVarPeakLoc{i}(m) = var(peakloc);
       DiaV = DiaV(~isnan(DiaV)); NeiV = NeiV(~isnan(NeiV)); DisV = DisV(~isnan(DisV));
       pinfo.spPVcrr.evtMeanDgnlCrr{i}(m) = mean(DiaV); pinfo.spPVcrr.evtMedDgnlCrr{i}(m) = median(DiaV); pinfo.spPVcrr.evtVarDgnlCrr{i}(m) = var(DiaV);
       pinfo.spPVcrr.evtMeanNgbrCrr{i}(m) = mean(NeiV); pinfo.spPVcrr.evtMedNgbrCrr{i}(m) = median(NeiV); pinfo.spPVcrr.evtVarNgbrCrr{i}(m) = var(NeiV);       
       pinfo.spPVcrr.evtMeanDistCrr{i}(m) = mean(DisV); pinfo.spPVcrr.evtMedDistCrr{i}(m) = median(DisV); pinfo.spPVcrr.evtVarDistCrr{i}(m) = var(DisV);
    end
    %            evt measures: fliped peak, diagonal, neighbor, distal
    pinfo.spPVcrr.evtMeanFPeakCrr{i} = NaN*ones(1,ip); pinfo.spPVcrr.evtMedFPeakCrr{i} = NaN*ones(1,ip); pinfo.spPVcrr.evtVarFPeakCrr{i} = NaN*ones(1,ip);
    pinfo.spPVcrr.evtMeanFPeakLoc{i} = NaN*ones(1,ip); pinfo.spPVcrr.evtMedFPeakLoc{i} = NaN*ones(1,ip); pinfo.spPVcrr.evtVarFPeakLoc{i} = NaN*ones(1,ip);
    pinfo.spPVcrr.evtMeanFDgnlCrr{i} = NaN*ones(1,ip); pinfo.spPVcrr.evtMedFDgnlCrr{i} = NaN*ones(1,ip); pinfo.spPVcrr.evtVarFDgnlCrr{i} = NaN*ones(1,ip);
    pinfo.spPVcrr.evtMeanFNgbrCrr{i} = NaN*ones(1,ip); pinfo.spPVcrr.evtMedFNgbrCrr{i} = NaN*ones(1,ip); pinfo.spPVcrr.evtVarFNgbrCrr{i} = NaN*ones(1,ip);
    pinfo.spPVcrr.evtMeanFDistCrr{i} = NaN*ones(1,ip); pinfo.spPVcrr.evtMedFDistCrr{i} = NaN*ones(1,ip); pinfo.spPVcrr.evtVarFDistCrr{i} = NaN*ones(1,ip);
    for m = 1:ip
       crrnow = data.spPVcrr.evtCrr{i}{m}; crrnow = flip(crrnow); %%%flip the crr matrix
       [peakV, peakloc, DiaV, NeiV, DisV] = findtargetvaluesincrr(crrnow, Xbin, parm);
       peakV = peakV(~isnan(peakV)); peakloc = peakloc(~isnan(peakloc));
       pinfo.spPVcrr.evtMeanFPeakCrr{i}(m) = mean(peakV); pinfo.spPVcrr.evtMedFPeakCrr{i}(m) = median(peakV); pinfo.spPVcrr.evtVarFPeakCrr{i}(m) = var(peakV);
       pinfo.spPVcrr.evtMeanFPeakLoc{i}(m) = mean(peakloc); pinfo.spPVcrr.evtMedFPeakLoc{i}(m) = median(peakloc); pinfo.spPVcrr.evtVarFPeakLoc{i}(m) = var(peakloc);
       DiaV = DiaV(~isnan(DiaV)); NeiV = NeiV(~isnan(NeiV)); DisV = DisV(~isnan(DisV));
       pinfo.spPVcrr.evtMeanFDgnlCrr{i}(m) = mean(DiaV); pinfo.spPVcrr.evtMedFDgnlCrr{i}(m) = median(DiaV); pinfo.spPVcrr.evtVarFDgnlCrr{i}(m) = var(DiaV);
       pinfo.spPVcrr.evtMeanFNgbrCrr{i}(m) = mean(NeiV); pinfo.spPVcrr.evtMedFNgbrCrr{i}(m) = median(NeiV); pinfo.spPVcrr.evtVarFNgbrCrr{i}(m) = var(NeiV);       
       pinfo.spPVcrr.evtMeanFDistCrr{i}(m) = mean(DisV); pinfo.spPVcrr.evtMedFDistCrr{i}(m) = median(DisV); pinfo.spPVcrr.evtVarFDistCrr{i}(m) = var(DisV);
    end
else
    disp(strcat('----------> warning: not enough cells: ', pinfo.general.finaldir{i}));
end
end

function [crr, pp] = find2DspPVcrr(XYbin, V1, V2)
%%%XYbin{1}(x), XYbin{2}(y), sessvector = V1/V2(ncell, y, x)
nx = numel(XYbin{1}); ny = numel(XYbin{2}); crr = NaN*ones(ny,nx); pp = NaN*ones(ny, nx);
for m = 1:ny
    for n = 1:nx
        [R,P] = corrcoef(V1(:,m,n), V2(:,m,n)); crr(m,n) = R(1,2); pp(m,n) = P(1,2);
    end
end
function [crr, pp] = find1DnowspPVcrr(Xbin, V1, V2)
%%%Xbin(x), evtvector = V1/V2(ncell, x)
nx = numel(Xbin); crr = NaN*ones(nx,nx); pp = NaN*ones(nx, nx);
for m = 1:nx
    for n = 1:nx
        [R,P] = corrcoef(V1(:,m), V2(:,n)); crr(m,n) = R(1,2); pp(m,n) = P(1,2);
    end
end

function [XYbin, sessvector, sessoccutime, sessname, Xbin, evtvector, evtoccutime, evtname, fjoint] = getPVvectors(pinfo, data, behav, bhdata, cellind, parm)
i = cellind(1); finaldirnow = pinfo.general.finaldir{i}; %sessnamenow = pinfo.general.sessionname{i};%%%%here already selected just one day to analyze
%%%%%%%%%%%%%%%%%%%%%%%% 2D PV vectors %%%%%%%%%%%%%%%%%%%%%%%%%
%%%% read in split-session 2D spatial bins and ratemaps FROM field dynam
%%%%      data: -- does not work - field dynamics only computed on cells with place fields - not all selected cells (cellind) have fields
%%%% Decide to change in fielddynam - since cells can be selected there and it is not necessay to impose,  unlike the 1D dynamics 
nseg = numel(data.fielddynam.sess2DRateMaps{cellind(1)});
if nseg == 0
    disp(['---------->  warning: session variables not computed; segmentation not found in ', finaldirnow]);
end
XY = cell(1, nseg); sessname = cell(1, nseg); ratemaps = cell(1, nseg); gridoccup = cell(1, nseg); segvalid = zeros(1, nseg);
withinsegnum = pinfo.fielddynam.sessWithinSegNum{i}; 
%%%only select the first and last of each session
ii = find(withinsegnum == 1); segvalid(ii) = ones(1, numel(ii)); %%%this contains all first seg's of all sessions
for (jk = 2:numel(ii)) segvalid(ii(jk)-1) = 1; end
segvalid(numel(segvalid)) = 1; %%%these contaon all last seg's of all sessions 
for (tt = 1:nseg) 
    withinnumnow = withinsegnum(tt);
    if segvalid(tt) %%%only do this for the first and last of a session
       sessname{tt} = strcat(pinfo.fielddynam.sessName{i}{tt}, '_', num2str(withinnumnow));
       sessnum = data.fielddynam.sess2DsegSessNum{i}(tt); 
       XY{tt} = data.fielddynam.sess2DXYbin{i}{sessnum}; %%%this is for each session
       gridoccup{tt} = data.fielddynam.sess2DOccuTime{i}{tt}; %%%%This is for each segment
       ratemaps{tt} = cell(1, numel(cellind));
       for (j = 1:numel(cellind))
           ratemaps{tt}{j} = data.fielddynam.sess2DRateMaps{cellind(j)}{tt};
       end
    end
end
ii = find(segvalid); XY = XY(ii); ratemaps = ratemaps(ii); gridoccup = gridoccup(ii); sessname = sessname(ii);
[XYbin, sessvector, sessoccutime] = alignandget2Dvectors(XY, ratemaps, gridoccup, finaldirnow); %%% assuming origin alignment
%%%%%%%%%%%%%%%%%%%%%% 1D PV vectors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% cannot read 1D spatial bins and ratemaps FROM field dynam or split-session field dynam data (only done on ONE with dominant field):
%%%% so this need to be recomputed
smParm.d1sigma = pinfo.parm.fSmooth1DSigma(i); smParm.d1Nsig = pinfo.parm.fSmooth1DNSigma(i); smParm.d1sm = pinfo.parm.fSmooth1D{i};
timeunit = pinfo.parm.timeunit(i); spiketime = cell(1, numel(cellind));
for (m=1:numel(cellind)) 
    spiketime{m} = timeunit*data.spike.spiketime{cellind(m)}; 
end
evName = pinfo.general.eventname{i}; evTime = data.events.eventtimes{i}; evType = pinfo.parm.eventtype{i}; 
evsel = checkifcorrectevents(evName, evType, parm.evkeyword, parm.evkeytype, parm.evkeynoword, parm.evkeynotype); %%%select events according to user defined criteria
evind = find(evsel==1); evName = evName(evind); evTime = evTime(evind); evType = evType(evind); 
nev = 0; X = []; ratemaps = []; gridoccup = []; xjoint = []; evtname = []; %%%%%% valid event groups will accumulate as each event is parsed to early/late groups (group 1/2)
for j = 1:numel(evName)
    evSess = identifysession(evTime{j}, pinfo.general.sessionname{i}, pinfo.general.sessionstartT{i}, pinfo.general.sessionendT{i});
    posid = []; evid = [];
    if (~isempty(evSess))
        posid = find( strcmp(behav.general.finaldir, finaldirnow) & strcmp(behav.general.sessname, evSess) );
    end
    if numel(posid) == 1
        if (isfield(behav.general, 'eventname'))
            evid = find(strcmp(behav.general.eventname{posid}, evName{j}));
        else
            evid = find(strcmp(behav.behavior.eventname{posid}, evName{j}));
        end
    end
    if (numel(posid)~=1)||(numel(evid)~=1)
        disp(['---------->  not computed: no or more than 1 positon/event files match the session: ', finaldirnow, '___', evName{j}]);
    else       
        xbin = bhdata.event.Xbin{posid}{evid}; xjointnow = bhdata.event.Xjoint{posid}{evid};
        if ~isempty(xjointnow)
           lapoccutime = bhdata.event.LapOccuptime{posid}{evid}; framerate = behav.parm.framerate(posid);
           lappostime = bhdata.event.LapAllPostimestamp{posid}{evid}; lapx = bhdata.event.LapAllX{posid}{evid};
           [lapind, grpname] = findlapgroups(evName{j}, evTime{j}, numel(evTime{j}.start), parm);
           for tt = 1:numel(grpname)
               evok.start = evTime{j}.start(lapind{tt}); evok.ent = evTime{j}.ent(lapind{tt}); occutime = zeros(numel(xbin), 1);
               for (m = 1:numel(lapind{tt}))
                    occutime = occutime + lapoccutime{lapind{tt}(m)};
               end
               nev = nev + 1; X{nev} = xbin; gridoccup{nev} = occutime; xjoint{nev} = xjointnow; 
               evtname{nev} = strcat(evName{j}, '_', num2str(tt)); ratemaps{nev} = cell(1, numel(cellind));
               for m = 1:numel(cellind)
                   [D1ratenow, ~] = getlinearmap(spiketime{m}, evok, lappostime(lapind{tt}), lapx(lapind{tt}), xbin, occutime, smParm, 1/framerate);
                   ratemaps{nev}{m} =  D1ratenow;
               end
           end
        end
    end
end
[Xbin, evtvector, evtoccutime, evtname, fjoint] = alignget1Dvectors(X, ratemaps, gridoccup, xjoint, evtname, finaldirnow);

function [peakV, peakloc, DiaV, NeiV, DisV] = findtargetvaluesincrr(crrnow, Xbin, parm)
%%% parm.pvngbrrange[min max]; parm.pvdistrange[min max]; parm.pvsearchpeakrange[min max];
%%%% search for peak value and peak location first; 
nx = numel(Xbin); peakV = NaN*ones(1, nx); peakloc = NaN*ones(1, nx); binsize = median(diff(Xbin));
P1 = ceil(parm.peakrange(1)/binsize); P2 = ceil(parm.peakrange(2)/binsize);
for i = 1:nx
    ind = intersect([i+P1:i+P2], 1:nx); 
    aa = crrnow(i,ind); ii = find(~isnan(aa));
    if ~isempty(ii) 
        aa = aa(ii); peakV(i) = max(aa); tt=find(aa==peakV(i)); 
        if numel(tt)>1
           inow = randperm(numel(tt)); tt = tt(inow(1)); %%%take one of peak locations if multiple
        end
        peakloc(i) = Xbin(ind(ii(tt)))-Xbin(i);
        %%peakloc(i) = median(Xbin(ind(ii(tt)))-Xbin(i));
    end
end
%%% compute values in different catogeries
N1 = ceil(parm.ngbrrange(1)/binsize); N2 = ceil(parm.ngbrrange(2)/binsize); %%%neighbor bourdaries
D1 = ceil(parm.distrange(1)/binsize); D2 = ceil(parm.distrange(2)/binsize); %%%distal boundaries
isdia = false(nx); %%%be careful here: false(n) = zerso(n,n) a square array
isnei = false(nx); isdist = false(nx); %%% select if 2d indices are true for a category
for m = 1:nx %%%for every row
    for n = 1:nx %%% every column
        aa = abs(n-m);
        if aa == 0 %%% diagonal
            isdia(m,n) = 1; 
        end
        if (aa>=N1)&&(aa<=N2)  %%% neighbor
            isnei(m,n) = 1; 
        end 
        if (aa>=D1)&&(aa<=D2)  %%% distal
            isdist(m,n) = 1; 
        end 
    end
end
DiaV = crrnow(isdia); NeiV = crrnow(isnei); DisV = crrnow(isdist);

function [XYbin, sessvector, sessoccutime] = alignandget2Dvectors(XY, ratemaps, gridoccup, finaldir)
%%%% check 2D grid consistency, arrange into PVs - assuming same 2D space
%%%% across sessions; otherwise aligned at origin [0,0] with same binsize
XYbin{1} = []; XYbin{2} = []; sessvector = []; sessoccutime = [];
nsess = numel(gridoccup); 
if nsess>= 1
   nX = zeros(1, nsess); nY = zeros(1, nsess); ncell = zeros(1, nsess);
   for i = 1:nsess
       nX(i) = numel(XY{i}{1}); nY(i) = numel(XY{i}{2}); ncell(i) = numel(ratemaps{i});
   end
   minx = min(nX); miny = min(nY); minc = min(ncell);
   if (~isempty(find(nX~=minx, 1))) || (~isempty(find(nY~=miny, 1))) 
       disp(['----------> warning: 2D spatial dimensions do not matach across sessions: ' finaldir]); 
   end
   if (~isempty(find(ncell~=minc, 1))) 
       disp(['----------> warning: number of 2D cell rate maps do not matach across sessions: ' finaldir]); 
   end
   XYbin{1} = XY{1}{1}(1:minx); XYbin{2} = XY{1}{2}(1:miny);
   sessvector = cell(1, nsess); sessoccutime = cell(1, nsess);
   for tt = 1:nsess
       sessvector{tt} = NaN*ones(minc, miny, minx); %%%%initial vector assignment: rate of cell minc at [minx miny]
       sessoccutime{tt} = gridoccup{tt}(1:miny, 1:minx);
       for k = 1:minc
           sessvector{tt}(k,:,:) = ratemaps{tt}{k}(1:miny, 1:minx);
       end
   end
end

function [Xbin, evtvector, evtoccutime, evtname, fjoint] = alignget1Dvectors(X, ratemaps, gridoccup, xjoint, evtname, msgflag)
%%%% align at all middle joints, arrange into PVs 
%%%% the issue here is that some trajectories have oppsoite directions (joints) or non-overlapping portions
%%%% Assumption: all trajs have similar lenght and same number of joints
%            This is ensured by event/session selection at the beginning
Xbin = []; evtvector = []; evtoccutime = []; nevt = numel(gridoccup); %disp(evtname);
if nevt == 1
    Xbin = X{1}; fjoint = xjoint{1}; evtoccutime{1} = gridoccup{1};
    evtvector = getratevectors(ratemaps{1}, numel(ratemaps{1}), numel(Xbin)); 
end
if nevt>1
%%%%first need to determine which events are forward or backward (flip)
%%%%%%%% there is no other way around, has to rely on format Track1_leftrgiht, etc..    
   [~, nameword] = strtok(evtname, '_'); flipsel = settrajdir(nameword);
   for i = 1:nevt
       if flipsel(i) %%%%assuming rate/refrate xx/refxx are opposite trajectories
          disp(['----------> flipping trajectory: ' evtname{i}]);
          xjoint{i} = max(xjoint{i}) - flip(xjoint{i}); 
          gridoccup{i} = flip(gridoccup{i}); X{i} = max(xjoint{i})-flip(X{i}); %%%need to be very careful here, some xx not start at 0 or end at max(xjoint)
          for (m = 1:numel(ratemaps{i}))
              ratemaps{i}{m} = flip(ratemaps{i}{m}); 
          end
       end
   end
%%%%% now get the reference xbin: get the binsize and then use the average joints   
   binsize = 0; nj = zeros(1, nevt); evok = true(1, nevt);
   for i = 1:nevt
       binsize = binsize + X{i}(2) - X{i}(1); nj(i) = numel(xjoint{i});
   end
   binsize = binsize/nevt; nj = min(nj); fjoint = zeros(nevt,nj);
   for i = 1:nevt
       for j = 1:nj
           fjoint(i,j) = xjoint{i}(j);
       end
   end
   fjoint = mean(fjoint); %%%now refjoint(1, nj)
   Xbin = [];
   for i=2:nj
       Xbin = [Xbin fjoint(i-1):binsize:fjoint(i)]; 
   end
   %%%%%%%%%%%%%%%%now align each traj to reference: Xbin, fjoint
   %        First: make sure cell numbers are the same
   ncell = zeros(1, nevt);
   for i = 1:nevt ncell(i) = numel(ratemaps{i}); end
   nc = min(ncell); 
   if (~isempty(find(ncell~=nc, 1))) 
       disp(['----------> warning: number of 1D cell rate maps do not matach across sessions: ' msgflag]); 
   end
   for i = 1:nevt ratemaps{i} = ratemaps{i}(1:nc); end
   %        align to reference
   isok = true(1,nevt); 
   for i = 1:nevt
       [ratemaps{i}, gridoccup{i}, isok(i)] = realignjointtoref(xjoint{i}, fjoint, ratemaps{i}, gridoccup{i}, X{i}, Xbin, msgflag);
   end
   ratemaps = ratemaps(isok); evtoccutime = gridoccup(isok); evtname = evtname(isok);
   evtvector = getratevectors(ratemaps, nc, numel(Xbin)); 
end
function evtvector = getratevectors(ratemaps, nc, nbin)
nev = numel(ratemaps); evtvector = cell(1, nev); 
for tt = 1:nev
    evtvector{tt} = NaN*ones(nc, nbin); %%%%initial vector assignment: rate of cell minc at [minx miny]
    for k = 1:nc
        evtvector{tt}(k,:) = ratemaps{tt}{k};
    end
end
function [rr, occut, isok] = realignjointtoref(xjoint, refjoint, rate, occutime, xx, refxx, trajflag)
%%%%realign trajectories assuming the center landmarks are matched
nc = numel(rate); rr = cell(1, nc); occut = []; isok = true(1); 
occutime = occutime'; 
for m = 1:nc 
    rate{m} = rate{m}'; 
    if isempty(rate{m}) isok = 0; break; end %%%% sometimes rate maps are not computed for certain sessions when not enough laps were run
end %%%%there are multiple cells here, but only one occutime 
if numel(xjoint) ~= numel(refjoint)
    disp(['----------> warning: trajectory joint points do not match; split PVcrr not computed for ', trajflag]); isok = 0;
end
if isok
   occut = zeros(numel(refxx), 1); 
   for m = 1:nc 
       rate{m} = rate{m}'; rr{m} = NaN*ones(1, numel(refxx)); 
   end 
   nj = numel(xjoint); 
   if (nj == 2)
       if (numel(xx) == numel(refxx))
           rr = rate; occut = occutime;
       else
           occut = interpn(xx, occutime, refxx); 
           for m = 1:nc rr{m} = interpn(xx, rate{m}, refxx); end 
       end
   else
       for (i = 2:numel(xjoint)) 
           jj1 = find( (xx>=xjoint(i-1)) & (xx<xjoint(i))); jj2 = find( (refxx>=refjoint(i-1)) & (refxx<refjoint(i))); 
           if (numel(jj1) == numel(jj2))
               occut(jj2) = occutime(jj1); for m=1:nc rr{m}(jj2) = rate{m}(jj1); end  
           else
               if (numel(jj1)>=2)
                   xnow = refxx(jj2) - refjoint(i-1);               
                   for m = 1:nc 
                       rr{m}(jj2) = interpn(xx(jj1)-xjoint(i-1), rate{m}(jj1), xnow); 
                   end 
                   o1now = interpn(xx(jj1)-xjoint(i-1), occutime(jj1), xnow); o1now = o1now*sum(occutime(jj1))/sum(o1now); %%%to keep overall time the same
                   occut(jj2) = o1now;
               else
                   if numel(jj1) == 1
                      for m = 1:nc rr{m}(jj2) = rate{m}(jj1)*ones(1, numel(jj2)); end 
                      occut(jj2) = occutime(jj1)*ones(1, numel(jj2))/numel(jj2); 
                   end
               end
           end
       end
   end
end
function flipsel = settrajdir(nameword)%%%%have to do a namesearch
flipsel = false(numel(nameword));
for i = 1:numel(nameword)
    if strncmp(nameword{i}, '_RL', 3) flipsel(i) = 1; end  %%for _RL
    if strncmpi(nameword{i}, '_rightleft', 7) flipsel(i) = 1; end  %%for _rightleft
    if strncmpi(nameword{i}, '_NCW', 3) flipsel(i) = 1; end  %%for _NCW,
    if strncmpi(nameword{i}, '_center', 7) flipsel(i) = 1; end  %%for _centerleft, _centerright
    if strncmp(nameword{i}, '_CL', 2) flipsel(i) = 1; end  %%for _CL
    if strncmp(nameword{i}, '_CR', 2) flipsel(i) = 1; end  %%for _CR
end

function refjoint = getrefjoint(pinfo, data, behav, bhdata, cellind, refevName)
refjoint = []; finaldirnow = pinfo.general.finaldir{cellind}; 
evName = pinfo.general.eventname{cellind}; evTime = data.events.eventtimes{cellind}; evType = pinfo.parm.eventtype{cellind};    
iii = find(strcmp(evName, refevName)); 
if (numel(iii) == 1)
       j = iii(1); 
       %%%%locate event position data
       evSess = identifysession(evTime{j}, pinfo.general.sessionname{cellind}, pinfo.general.sessionstartT{cellind}, pinfo.general.sessionendT{cellind});
       posid = []; evid = [];
       if (~isempty(evSess))
           posid = find( strcmp(behav.general.finaldir, finaldirnow) & strcmp(behav.general.sessname, evSess) );
       else
           disp(['--------> warning: no session identified for the event: ', finaldirnow, '___', evName{j}]);  
       end
       if numel(posid) == 1
            if (isfield(behav.general, 'eventname'))
                evid = find(strcmp(behav.general.eventname{posid}, evName{j}));
            else
                evid = find(strcmp(behav.behavior.eventname{posid}, evName{j}));
            end
       else
            disp(['--------> warning: no or more than 1 positon files match the session: ', finaldirnow, '___', evSess]); 
       end
       if (numel(posid)~=1)||(numel(evid)~=1)
            disp(['--------> warning: no or more than 1 position event files match the event: ', finaldirnow, '___', evName{j}]);
       else
            refjoint = bhdata.event.Xjoint{posid}{evid}; 
       end
else
       disp(['--------> warning: no or multiple event files found match with ', refevName]);
end

function evSess = identifysession(evTime, sessionname, startT, endT)
evSess = []; minevstart = min(evTime.start); maxevend = max(evTime.ent);
iii = find( (startT <= minevstart) & (endT >= maxevend) );
if (numel(iii) == 1)
    evSess = sessionname{iii};
end

function evsel = checkifcorrectevents(evName, evType, evkeyword, evkeytype, evkeynoword, evkeynotype)
evsel = true(1,numel(evName));
for (i = 1:numel(evName))
    if (~isempty(evkeytype))||(~isempty(evkeyword)) %%%for inclusion
      if isempty(evkeytype)
         if ~contains(lower(evName{i}), lower(evkeyword)) evsel(i) = 0; end 
      elseif isempty(evkeyword)
         if ~strncmpi(evType{i}, evkeytype, 3) evsel(i) = 0; end 
      else
         if ~(contains(lower(evName{i}), lower(evkeyword)) && strncmpi(evType{i}, evkeytype, 3))
             evsel(i) = 0; 
         end 
      end
    end
    if (~isempty(evkeynotype))||(~isempty(evkeynoword)) %%%for exclusion
       if isempty(evkeynotype)
          if contains(lower(evName{i}), lower(evkeynoword)) evsel(i) = 0; end 
       elseif isempty(evkeynoword)
          if strncmpi(evType{i}, evkeynotype, 3) evsel(i) = 0; end
       else
          if contains(lower(evName{i}), lower(evkeynoword)) || strncmpi(evType{i}, evkeynotype, 3)
              evsel(i) = 0; 
          end   
       end
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

function PT = findevtpairtype(evtname1, evtname2)
PT = 'ND';
if (~isempty(evtname1)) && (~isempty(evtname2))
if strcmp(evtname1, evtname2)
    PT = 'self';
else
    PT = 'crossSess'; P2 = []; P3 = [];
    [s1, t1] = strtok(evtname1, '_'); [s2, t2] = strtok(evtname2, '_');
    if strcmp(s1, s2) PT = 'sameSess'; end
    [s11, t11] = strtok(t1(2:numel(t1)), '_'); [s21, t22] = strtok(t2(2:numel(t2)), '_'); 
    if ~isempty(s11) || ~isempty(s21) 
       if strcmp(s11, s21) 
          P2 = 'sameTraj';
       else
          P2 = 'crossTraj';
       end
    end
    if ~isempty(t11) || ~isempty(t22) 
       if ~strcmp(t11, t22) 
          P3 = 'crossGrp';
       end
    end
    %if evtname1(numel(evtname1)) ~= evtname2(numel(evtname2)) P3 = 'crossGrp'; end
    if (~isempty(P3)) && (~isempty(P2))
       PT = strcat(PT, '_', P2, '_', P3);
    elseif ~isempty(P2)
       PT = strcat(PT, '_', P2);
    elseif ~isempty(P3)
       PT = strcat(PT, '_', P3);
    end
end
end

function [grplapind, grpname] = findlapgroups(sessevname, evT, nlap, parm)
%%%%%%now identify early/late groups within each session
lapind = 1:nlap; grplapind = cell(1,2); grpname = cell(1,2);
if strncmpi(parm.splitmethod, 'lap', 2) %%%group by laps
   iii = find(lapind <= parm.grouplength);
   if numel(iii)>= parm.grpminlapnum
       earlyind = lapind(iii); 
       grpname{1} = sessevname; grplapind{1} = earlyind; 
   end
   iii = find(max(lapind) - lapind <= parm.grouplength);
   if numel(iii)>= parm.grpminlapnum
      lateind = lapind(iii); 
      grpname{2} = sessevname; grplapind{2} = lateind; 
   end
else %%%group by time
   tstartnow = evT.start(lapind); tendnow = evT.ent(lapind);
   ts = min(tstartnow); te = max(tendnow);
   iii = find(tendnow-ts <= parm.grouplength); 
   if numel(iii)>= parm.grpminlapnum
      earlyind = lapind(iii); 
      grpname{1} = sessevname; grplapind{2} = earlyind; 
   end
   iii = find(te - tstartnow <= parm.grouplength); 
   if numel(iii)>= parm.grpminlapnum
       lateind = lapind(iii); 
       grpname{2} = sessevname; grplapind{2} = lateind; 
   end 
end

function [D1rate, nact] = getlinearmap(spiketime, evTime, lappostime, lapx, xbin, occutime, smParm, frametime)
nx = numel(xbin); D1rate = NaN*ones(nx,1); nlap = numel(lappostime); count = zeros(nx,1); 
frametime = 1.03*frametime; %%%slightly expand frametime to assure detection, due to slightly variable frame-to-frame time
if (nx > 1)
    binsize = xbin(2) - xbin(1); sigma = round(smParm.d1sigma/binsize); %%sigma now in number of bins
%%%%%count spikes in xbin lap by lap
for (i = 1:nlap)
    counti = zeros(nx,1);
    spikenow = sort( spiketime( (spiketime>=evTime.start(i)) & (spiketime<=evTime.ent(i)) ) );
    spikex = NaN*ones(size(spikenow)); [laptimenow, iii] = sort(lappostime{i}); allxnow = lapx{i}(iii);
    lastpoint = 1;  %this only for saving time
    for (j = 1:numel(spikenow))
         for (k = lastpoint:numel(laptimenow)-1) %find corresponding time in position data
             if abs(laptimenow(k) - spikenow(j)) <= frametime 
             %if (laptimenow(k) <= spikenow(j)) && (laptimenow(k+1) > spikenow(j)) 
                 spikex(j) = allxnow(k); lastpoint = k; break; 
             end
         end
    end
    for (j = 1:nx-1) counti(j) = numel(find((spikex>=xbin(j)) & (spikex<xbin(j+1)))); end
    counti(nx) = numel(find(spikex>=xbin(nx)));
    count = count + counti;
end
%%%%%compute firing rate
for (i = 1:nx)
     if (occutime(i) ~= 0) D1rate(i) = count(i)/occutime(i); end
end
nact = numel(find(count>0))/numel(find(occutime>0));
%%%%%%%%%%The issue here is that, for single lap and after stopping removal, rates have NaN values at some position points
%%now do a 1d smoothing
if (strcmpi(smParm.d1sm, 'yes'))
    otime = ones(size(occutime));
    D1rate = OneDSmooth_new(D1rate, otime, sigma, smParm.d1Nsig); %%%smooth first - NaNs (occutime = 0) are not removed by smoothing
    %%%need to do a 1D interpolation of the NaN data points
    X = (1:numel(D1rate))'; iii = find(D1rate > -10); 
    if numel(iii)>1
       R = D1rate(iii); Y=X(iii);
       D1rate = interpn(Y, R, X);
       %%%%The following has the problem of reducing rate at the end of traj
       %XX = [-smParm.d1Nsig*sigma:1:smParm.d1Nsig*sigma]; 
       %wind = normpdf(XX, 0, sigma); wind = wind/sum(wind); nw = numel(wind); %will have a delay = (n-1)/2 bins
       %D1rate = conv(D1rate, wind); D1rate = D1rate((nw-1)/2+1:numel(D1rate)-(nw-1)/2);
    end
end
end
