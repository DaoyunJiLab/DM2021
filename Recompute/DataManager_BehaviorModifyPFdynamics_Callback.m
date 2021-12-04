function DataManager_BehaviorModifyPFdynamics_Callback
%%%%%%%%%%%Plot lap by lap changes of place field/firing features after correcting with speed
%%%%%%%%%%%    regression.
disp('Computing speed/headdir correction of place field dynamics');

hf = gcbf; hgroup = getappdata(hf, 'hgroup'); hfield = getappdata(hf, 'hfield');
plotparm = getappdata(hf, 'plotparm'); vv = plotparm.showdetail; %if vv=1, shoe details of computation
%cc = plotparm.compute; %if cc=0, assign parameters for computation only; if 1, do the real computation
ow = plotparm.overwrite; %if ow=1, plot into the current database
fname = get(hf, 'Name'); ftt = strfind(fname, '__'); currentfilename = fname(ftt+2:numel(fname));
[~, ~, ee] = fileparts(currentfilename);
[writefilename, ok] = getoutputfile(ee, currentfilename, ow);

if ok && (~strcmp(ee, '.spikedb'))
    disp('-----> not a spikedb; aborted'); ok = 0;
else
    pinfo = getappdata(hf, 'pinfo'); data = getappdata(hf, 'data'); 
    %get selected cellind
    groupselection = getappdata(hgroup, 'selection'); cellind = []; grpind = find(groupselection == 1);
    for (kk = 1:numel(grpind)) cellind = union(cellind, data.grouplist.groupindex{grpind(kk)}); end 
    if numel(cellind) ==0 
        disp('-----> no cells selected; aborted'); ok = 0;
    else  %%%check if field dynam already computed
        if ~isfield(pinfo, 'fielddynam') disp('-----> field dynamics not computed yet; aborted'); ok = 0; end
    end
end
if ok 
    if (plotparm.linkbehav == 0)
        disp(['--------> no behavioral data linked']); ok = 0;
    else
        behav = getappdata(hf, 'behav'); bhdata = getappdata(hf, 'bhdata');
    end
end
if ok
    S = {'Speed'; 'Head direction'};
    [param, ok] = listdlg('ListString', S, 'PromptString', 'Correct for what parameters?'); 
    %param = S(sel);
end
% if ok
%     T = {'Speed'; 'Head direction'};
%     [param, ok] = listdlg('ListString', S, 'PromptString', 'Correct for what parameters?'); 
%     %param = S(sel);
% end
if ok
   [pinfo, data] = initialnewfield(pinfo, data);
   %%% Strategy: get all unique dates; For each day, determine raw session
   %%% linear speeds; For each cell, determine speed within fields for the
   %%% right traj/laps; regress and subtract but maintain ref values
   for (ttt = 1:numel(cellind))
       i = cellind(ttt);
       disp(strcat('-----> correcting field dynamics ---', pinfo.general.parmfile{i}));
       %%%get field boundaries
       fstart = pinfo.fielddynam.runFstart{i}; fend = pinfo.fielddynam.runFend{i}; 
       %%%get speed/headdir from behavioral database
       finaldirnow = pinfo.general.finaldir{i}; 
       evName = pinfo.general.eventname{i}; evTime = data.events.eventtimes{i}; %evType = pinfo.parm.eventtype{i}; 
       %%%%identify behavioral parameters lap by lap
       lapevname = pinfo.fielddynam.run1DEvtname{i}; lapevnum = pinfo.fielddynam.runWithinLapNum{i};
       nlap = numel(lapevname);
       trajspeed = NaN*ones(1, nlap); fspeed = NaN*ones(1, nlap); fheaddir = NaN*ones(1, nlap);
       for (kkk = 1:nlap)
           spikeevid = find(strcmp(evName, lapevname{kkk}));
           spikeevSessname = identifysession(evTime{spikeevid}, pinfo.general.sessionname{i}, pinfo.general.sessionstartT{i}, pinfo.general.sessionendT{i});
           posid = []; evid = []; %%%in behav database
           if (~isempty(spikeevSessname))
                posid = find( strcmp(behav.general.finaldir, finaldirnow) & strcmp(behav.general.sessname, spikeevSessname) );
           end
           if numel(posid == 1)
                if (isfield(behav.general, 'eventname'))
                    evid = find(strcmp(behav.general.eventname{posid}, evName{spikeevid}));
                else
                    evid = find(strcmp(behav.behavior.eventname{posid}, evName{spikeevid}));
                end
           end
           if (numel(posid)~=1)||(numel(evid)~=1)
               disp(['------------------> 1D field dynamics correction not computed for this event: no or more than 1 positon/event files match the session: ', finaldirnow, '___', evName{j}]);
           else
               xbin = bhdata.event.Xbin{posid}{evid}; framerate = behav.parm.framerate(posid);
               lappostime = bhdata.event.LapAllPostimestamp{posid}{evid}{lapevnum(kkk)}; %%% Be careful: some (stopping) data points may be removed
               lapx = bhdata.event.LapAllX{posid}{evid}{lapevnum(kkk)}; lapspeed = bhdata.event.LapAll1DSpeed{posid}{evid}{lapevnum(kkk)}; 
               lapdir = bhdata.event.LapAllHeadDir{posid}{evid}{lapevnum(kkk)};
               trajspeed(kkk) = bhdata.event.LapMeanSpeed{posid}{evid}(lapevnum(kkk));
               if ~isempty(fstart)
                  [fspeed(kkk), fheaddir(kkk)] = findlapparm(lappostime, lapx, lapspeed, lapdir, xbin, fstart, fend, data, i, framerate);
               end
           end
       end
       %%%%%assign pinfo variables below
       reflap = pinfo.fielddynam.runRefLapNum{i}; fheaddir = anglerenorm(fheaddir); %%%head dir now in [-180 180]
       pinfo = correctandputvarvalues(pinfo, i, reflap, trajspeed, fspeed, fheaddir, param, fstart);
   end
end
if ok 
    if (ow)
        iii = get(hf, 'Children'); delete(iii(~strcmp(get(iii, 'Type'), 'uimenu'))); 
        fname = get(hf, 'Name'); ftt = strfind(fname, '__'); 
        newname = fname(1:ftt-1); set(hf, 'Name', newname); hmain = hf; 
    else
        hmain = DataManager_DataManager_Callback; plotparm.linkspike = 1; plotparm.linkeeg = 0; plotparm.linkbehav = 0;
    end
    setappdata(hmain,'plotparm', plotparm);
    set(hmain, 'Name', strcat(get(hmain, 'Name'), '__', writefilename));
    s = whos('data'); Gb1 = s.bytes/(1024^3);
    s = whos('pinfo'); Gb2 = s.bytes/(1024^3);
    if Gb1 + Gb2 < 2
       save(writefilename, 'pinfo', 'data'); %disp('Small one');
    else
       save(writefilename, 'pinfo', 'data', '-mat', '-v7.3'); %disp('Big one');
    end
    %save(writefilename, 'pinfo', 'data', '-mat', '-v7.3');
    setappdata(hmain, 'pinfo', pinfo); setappdata(hmain, 'data', data);
    DataManager_PlotSpikeDatabase(hmain, pinfo, data);
    pinfo = []; data = [];
end
disp('**********************');

function [fspeed, fheaddir] = findlapparm(lappostime, lapx, lapspeed, lapdir, xbin, fstart, fend, data, cellind, framerate)
fspeed = NaN; fheaddir = NaN;
xxbin = data.fielddynam.evt1Dxbin{cellind}; %%this is the average, aligned trajectory
%%%determine which bins in xbin belong to field [fstart fend] of xxbin, assuming binsizes are the same
%%%   use center alignment schemes: center of xbin is the same as center of xxbin
mid = median(xxbin); 
fstart = fstart - mid; fend = fend - mid; lapx = lapx - mid;
%disp([fstart fend]); disp(['***---']); disp(size(lapx));
iii = find( (lapx>=fstart) & (lapx<=fend) );
if ~isempty(iii)
    %x = lapx(iii); 
    fdir = lapdir(iii); fdir = fdir(~isnan(fdir)); fsp = lapspeed(iii); fsp = fsp(~isnan(fsp)); 
         %%%be careful, there are gaps - use framerate instead of duration
         %%%although there could be stopping points within the field, but they were removed from both x, fdir, lappostime and from spikes
         %%% What if a field is traveled for multiple (N) times -- speed (not headdir) needs to be compensated for this 
              % this not work in this case  fspeed = (x(numel(x)) - x(1))/(numel(iii)/framerate); 
    if (~isempty(fsp))
        fspeed = mean(fsp); %%%decide to do mean instead of median  since stopping already removed in most cases
    end
    if (~isempty(fdir))
        [aa,~] = size(fdir); if (aa==1) fdir = fdir'; end %%fdir needs to be column vector for circ_median in radians
        fheaddir = circ_mean(fdir*2*pi/360); fheaddir = fheaddir*360/2/pi;
    end
end

function pinfo = correctandputvarvalues(pinfo, i, reflap, trajspeed, fspeed, fheaddir, param, fstart)
%%%re-compute and correct variables in the field cfielddynam, which mirrors everything in fielddynam
%%%%%%%%%%%%%%% evt trajectory variables - use trajspeed for 
pinfo.cfielddynam.runLapMeanspeed{i} = trajspeed;
pinfo.cfielddynam.runLapMeanrate{i} = regresandcorrect_traj(pinfo.fielddynam.runLapMeanrate{i}, trajspeed);
pinfo.cfielddynam.runLapSparsity{i} = regresandcorrect_traj(pinfo.fielddynam.runLapSparsity{i}, trajspeed);
pinfo.cfielddynam.runLapSpatialInfo{i} = regresandcorrect_traj(pinfo.fielddynam.runLapSpatialInfo{i}, trajspeed);
pinfo.cfielddynam.runLapSpatialInfotime{i} = regresandcorrect_traj(pinfo.fielddynam.runLapSpatialInfotime{i}, trajspeed);
pinfo.cfielddynam.runLapActarea{i} = regresandcorrect_traj(pinfo.fielddynam.runLapActarea{i}, trajspeed);
pinfo.cfielddynam.runLapCrr{i} = regresandcorrect_traj(pinfo.fielddynam.runLapCrr{i}, trajspeed);
%%%%%%%%%%%%%%% evt field variables - use fspeed, fheaddir(if specified in param) for 
if ~isempty(fstart) %%if field exists
%[pinfo.cfielddynam.runLapfComX{i},pinfo.cfielddynam.runRefLapfComX{i}] = regresandcorrect_field(pinfo.fielddynam.runLapfComX{i}, fspeed, fheaddir, param, reflap);
%[pinfo.cfielddynam.runLapfPeakX{i},pinfo.cfielddynam.runRefLapfPeakX{i}] = regresandcorrect_field(pinfo.fielddynam.runLapfPeakX{i}, fspeed, fheaddir, param, reflap);
%[pinfo.cfielddynam.runLapfMeanRate{i},pinfo.cfielddynam.runRefLapfMeanRate{i}] = regresandcorrect_field(pinfo.fielddynam.runLapfMeanRate{i}, fspeed, fheaddir, param, reflap);
%[pinfo.cfielddynam.runLapfPeakRate{i},pinfo.cfielddynam.runRefLapfPeakRate{i}] = regresandcorrect_field(pinfo.fielddynam.runLapfPeakRate{i}, fspeed, fheaddir, param, reflap);
%[pinfo.cfielddynam.runLapfFieldSize{i},pinfo.cfielddynam.runRefLapfFieldSize{i}] = regresandcorrect_field(pinfo.fielddynam.runLapfFieldSize{i}, fspeed, fheaddir, param, reflap);
%[pinfo.cfielddynam.runLapfSkew{i},pinfo.cfielddynam.runRefLapfSkew{i}] = regresandcorrect_field(pinfo.fielddynam.runLapfSkew{i}, fspeed, fheaddir, param, reflap);
%[pinfo.cfielddynam.runLapfKurt{i},pinfo.cfielddynam.runRefLapfKurt{i}] = regresandcorrect_field(pinfo.fielddynam.runLapfKurt{i}, fspeed, fheaddir, param, reflap);
pinfo.cfielddynam.runLapfSpeed{i} = fspeed; pinfo.cfielddynam.runLapfHeaddir{i} = fheaddir;
pinfo.cfielddynam.runLapfComX{i} = regresandcorrect_field(pinfo.fielddynam.runLapfComX{i}, fspeed, fheaddir, param, reflap);
pinfo.cfielddynam.runLapfPeakX{i} = regresandcorrect_field(pinfo.fielddynam.runLapfPeakX{i}, fspeed, fheaddir, param, reflap);
pinfo.cfielddynam.runLapfMeanRate{i} = regresandcorrect_field(pinfo.fielddynam.runLapfMeanRate{i}, fspeed, fheaddir, param, reflap);
pinfo.cfielddynam.runLapfPeakRate{i} = regresandcorrect_field(pinfo.fielddynam.runLapfPeakRate{i}, fspeed, fheaddir, param, reflap);
pinfo.cfielddynam.runLapfFieldSize{i} = regresandcorrect_field(pinfo.fielddynam.runLapfFieldSize{i}, fspeed, fheaddir, param, reflap);
pinfo.cfielddynam.runLapfSkew{i} = regresandcorrect_field(pinfo.fielddynam.runLapfSkew{i}, fspeed, fheaddir, param, reflap);
pinfo.cfielddynam.runLapfKurt{i} = regresandcorrect_field(pinfo.fielddynam.runLapfKurt{i}, fspeed, fheaddir, param, reflap);
end

function evSess = identifysession(evTime, sessionname, startT, endT)
evSess = []; minevstart = min(evTime.start); maxevend = max(evTime.ent);
iii = find( (startT <= minevstart) & (endT >= maxevend) );
if (numel(iii) == 1)
    evSess = sessionname{iii};
end

function [pinfo, data] = initialnewfield(pinfo, data)
nspike = numel(pinfo.general.parmfile);
if (~isfield(pinfo, 'cfielddynam')) pinfo.cfielddynam = []; end
   %%%session field dynam -- do not compute the corrected version of 2D variables at this time
           %%%if (~isfield(pinfo.fielddynam, 'sessRefSegNum')) pinfo.fielddynam.sessRefSegNum = cell(1, nspike); end
   %%%evt ref laps - no need to repeat here
%if (~isfield(pinfo.cfielddynam, 'runRefLapNum')) pinfo.cfielddynam.runRefLapNum = cell(1, nspike); end
%if (~isfield(pinfo.cfielddynam, 'run1DEvtname')) pinfo.cfielddynam.run1DEvtname = cell(1, nspike); end 
%if (~isfield(pinfo.cfielddynam, 'runLapNum')) pinfo.cfielddynam.runLapNum = cell(1, nspike); end
%if (~isfield(pinfo.cfielddynam, 'runWithinLapNum')) pinfo.cfielddynam.runWithinLapNum = cell(1, nspike); end
   %%% evt trajectory properties
if (~isfield(pinfo.cfielddynam, 'runLapMeanspeed')) pinfo.cfielddynam.runLapMeanspeed = cell(1, nspike); end   
if (~isfield(pinfo.cfielddynam, 'runLapMeanrate')) pinfo.cfielddynam.runLapMeanrate = cell(1, nspike); end   
if (~isfield(pinfo.cfielddynam, 'runLapSparsity')) pinfo.cfielddynam.runLapSparsity = cell(1, nspike); end
if (~isfield(pinfo.cfielddynam, 'runLapSpatialInfo')) pinfo.cfielddynam.runLapSpatialInfo = cell(1, nspike); end
if (~isfield(pinfo.cfielddynam, 'runLapSpatialInfotime')) pinfo.cfielddynam.runLapSpatialInfotime = cell(1, nspike); end
if (~isfield(pinfo.cfielddynam, 'runLapActarea')) pinfo.cfielddynam.runLapActarea = cell(1, nspike); end
if (~isfield(pinfo.cfielddynam, 'runLapCrr')) pinfo.cfielddynam.runLapCrr = cell(1, nspike); end
   %%% evt field properties
if (~isfield(pinfo.cfielddynam, 'runLapfSpeed')) pinfo.cfielddynam.runLapfSpeed = cell(1, nspike); end
if (~isfield(pinfo.cfielddynam, 'runLapfHeaddir')) pinfo.cfielddynam.runLapfHeaddir = cell(1, nspike); end
if (~isfield(pinfo.cfielddynam, 'runLapfComX')) pinfo.cfielddynam.runLapfComX = cell(1, nspike); end
if (~isfield(pinfo.cfielddynam, 'runLapfPeakX')) pinfo.cfielddynam.runLapfPeakX = cell(1, nspike); end
if (~isfield(pinfo.cfielddynam, 'runLapfMeanRate')) pinfo.cfielddynam.runLapfMeanRate = cell(1, nspike); end
if (~isfield(pinfo.cfielddynam, 'runLapfPeakRate')) pinfo.cfielddynam.runLapfPeakRate = cell(1, nspike); end
if (~isfield(pinfo.cfielddynam, 'runLapfFieldSize')) pinfo.cfielddynam.runLapfFieldSize = cell(1, nspike); end
if (~isfield(pinfo.cfielddynam, 'runLapfSkew')) pinfo.cfielddynam.runLapfSkew = cell(1, nspike); end
if (~isfield(pinfo.cfielddynam, 'runLapfKurt')) pinfo.cfielddynam.runLapfKurt = cell(1, nspike); end
   %%%output reference values as well - these maintain the same after correction
%if (~isfield(pinfo.cfielddynam, 'runRefLapfComX')) pinfo.cfielddynam.runRefLapfComX = cell(1, nspike); end
%if (~isfield(pinfo.cfielddynam, 'runRefLapfPeakX')) pinfo.cfielddynam.runRefLapfPeakX = cell(1, nspike); end
%if (~isfield(pinfo.cfielddynam, 'runRefLapfMeanRate')) pinfo.cfielddynam.runRefLapfMeanRate = cell(1, nspike); end
%if (~isfield(pinfo.cfielddynam, 'runRefLapfPeakRate')) pinfo.cfielddynam.runRefLapfPeakRate = cell(1, nspike); end
%if (~isfield(pinfo.cfielddynam, 'runRefLapfFieldSize')) pinfo.cfielddynam.runRefLapfFieldSize = cell(1, nspike); end
%if (~isfield(pinfo.cfielddynam, 'runRefLapfSkew')) pinfo.cfielddynam.runRefLapfSkew = cell(1, nspike); end
%if (~isfield(pinfo.cfielddynam, 'runRefLapfKurt')) pinfo.cfielddynam.runRefLapfKurt = cell(1, nspike); end

function YY = regresandcorrect_traj(YY, XX)
if isempty(XX)
    YY = NaN*ones(size(YY));
else
    ind = find( (~isnan(XX)) & (~isnan(YY)) ); 
    if ~isempty(ind) 
          nind = setdiff( (1:numel(YY)), ind);
          xx= XX(ind); yy = YY(ind); %%%get rid of the non-numbers
          aaa = size(xx); if (aaa(1)==1) xx = xx'; end %%%all column vectors
          aaa = size(yy); if (aaa(1)==1) yy = yy'; end
          n = numel(xx); Xparm = [ones(n,1) xx]; %disp('****'); disp([yy Xparm]);
          [BBB, BBBint, R, Rint, stat] =  regress(yy, Xparm, 0.05);
          Y = Xparm * BBB; %regression result: line in (X, Y) = (xx, L4)
          YY(ind) = YY(ind)-Y'; YY(nind) = NaN*ones(1, numel(nind));
          %%%%maintain the same median value after the regression
          %YY = YY - median(YY) + median(yy);
    end
end

function YY = regresandcorrect_field(YY, XX, hd, param, reflap)
%param =1 if speed only; 2 = if hd only; [1 2] = if both
refnow = YY(reflap); mref = median(refnow(~isnan(refnow)));
if (numel(param) ==1) && (param == 2) %%%if correct by hd only
   XX = hd;
end
%%%flip all to column vectors
aaa = size(XX); if (aaa(1)==1) XX = XX'; end %%%all column vectors
aaa = size(YY); if (aaa(1)==1) YY = YY'; end
aaa = size(hd); if (aaa(1)==1) hd = hd'; end
if isempty(XX)
    YY = NaN*ones(size(YY));
else
    ind = find( (~isnan(XX)) & (~isnan(YY)) );
    ik = find( (~isnan(XX)) & (~isnan(YY)) & (~isnan(hd)) );
    if ~isempty(ind) 
       if (numel(param)==2) && (numel(ik)>=numel(ind)/2) %%%if sufficeint head direction data points
          nind = setdiff( (1:numel(YY)), ik);
          xx= XX(ik); yy = YY(ik); hd = hd(ik); %%%get rid of the non-numbers
          n = numel(xx); Xparm = [ones(n,1) xx hd]; %disp('****'); disp([yy Xparm]);
          [BBB, BBBint, R, Rint, stat] =  regress(yy, Xparm, 0.05);
          Y = Xparm * BBB; %regression result: line in [(xx,hd) Y] = (xx, L4)
          YY(ind) = YY(ind)-Y; YY(nind) = NaN*ones(numel(nind),1); 
        else 
          nind = setdiff( (1:numel(YY)), ind);
          xx= XX(ind); yy = YY(ind); %%%get rid of the non-numbers
          n = numel(xx); Xparm = [ones(n,1) xx]; %disp('****'); disp([yy Xparm]);
          [BBB, BBBint, R, Rint, stat] =  regress(yy, Xparm, 0.05);
          Y = Xparm * BBB; %regression result: line in (X, Y) = (xx, L4)
          YY(ind) = YY(ind)-Y; YY(nind) = NaN*ones(numel(nind),1);
       end
       %%%%maintain the same median value for reference laps after the regression
       %refthen = YY(reflap); mrefthen = median(refthen(~isnan(refthen)));
       %YY = YY + mref - mrefthen;
    end
end

function hdir = anglerenorm(hdir) %%%compute angle variation [-180 180]: input in [0 360]
xx = hdir; xx = xx*2*pi/360; aaa = size(xx); if (aaa(1)==1) xx = xx'; end %%%all column vectors in radians for circ_mean, circ_median
xx = xx-circ_mean(xx); xx = mod(xx*360/2/pi, 360); ii = find(xx>180); xx(ii) = xx(ii)-360; hdir = xx;

function [writefilename, okk] = getoutputfile(ext, currentfilename, ow)
okk = 1; writefilename = [];
if (ow == 0)
   [fname, pname] = uiputfile(fullfile(cd, ext), 'Write the new database to:');
   if (numel(fname)>1)
      writefilename = fullfile(pname, fname);
   else
      okk = 0;
   end
else
   writefilename = currentfilename;      
end

