function DataManager_FindLapConsistency_GoodLaps
%%% (pinfo,data,cellind,behav,bhdata,vv)
MinLapMeanSpeed = 50; %%%MinInstantSpeed = 20; %%% not used minimum instantaneous speed
MinLapNum = 4; %if number of good laps smaller than this, do not analyze 1d rate maps
MinLapMeanRate = 0.25; MaxLapMeanRate = 4; %%%%select active cells, but not interneurons

hf = gcbf; pinfo = getappdata(hf, 'pinfo'); data = getappdata(hf, 'data'); tagmark = get(gcbo, 'Tag');
hgroup = getappdata(hf, 'hgroup');hfield = getappdata(hf, 'hfield');
ok = 1; plotparm = getappdata(hf, 'plotparm'); vv = plotparm.showdetail; %if vv=1, shoe details of computation
cc = plotparm.compute; %if cc=0, assign parameters for computation only; if 1, do the real computation
ow = plotparm.overwrite; %if ow=1, plot into the current database
[writefilename, okk] = getoutputfile(hf, ow);

%get selected cellind
groupselection = getappdata(hgroup, 'selection'); cellind = []; grpind = find(groupselection == 1); 
for (kk = 1:numel(grpind)) cellind = union(cellind, data.grouplist.groupindex{grpind(kk)}); end

if (plotparm.linkbehav == 0)
    disp(['--------> no behavioral data linked; aborted.']);
else
behav = getappdata(hf, 'behav'); bhdata = getappdata(hf, 'bhdata');
if ~isfield(pinfo, 'GRfield')
    disp(['--------> refined (lap filtered) field computation not done; aborted.']);
else
if (~isempty(cellind)) && okk
    [pinfo,data] = DoItNow(pinfo,data,behav,bhdata,  cellind, MinLapMeanSpeed, MinLapNum, MinLapMeanRate, MaxLapMeanRate, vv);
    %%%%%%%%save and plot the new database
    if (~isempty(pinfo.general.parmfile))
        if (ok)
            save(writefilename, 'pinfo', 'data');
            if (ow)
               iii = get(hf, 'Children'); delete(iii(~strcmp(get(iii, 'Type'), 'uimenu'))); 
               fname = get(hf, 'Name'); ftt = strfind(fname, '__'); 
               newname = fname(1:ftt-1); set(hf, 'Name', newname); hmain = hf; 
            else
               hmain = DataManager_DataManager_Callback; plotparm.linkspike = 1; plotparm.linkeeg = 0; plotparm.linkbehav = 0;
            end
            setappdata(hmain,'plotparm', plotparm);
            DataManager_PlotSpikeDatabase(hmain, pinfo, data);
            set(hmain, 'Name', strcat(get(hmain, 'Name'), '__', writefilename));
            setappdata(hmain, 'pinfo', pinfo); setappdata(hmain, 'data', data); pinfo = []; data = [];
        end
    else
        disp('-------------> no cells in the database!');
    end
else
    disp('--------------> no groups selected, groups do not contain any cells');
end
end
end
disp('**********************');


%%%%%%%%Changes from the standard PF calculation:
%%%%%%%%%%%%%% 1. Filter out laps with obvious stopping points on linear track
%%%%%%%%%%%%%%  *********** 2. Filter out stopping periods on open platform and linear track ******* Not done yet, need to recompute behavdb 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pinfo,data] = DoItNow(pinfo,data,behav,bhdata, cellind, MinLapMeanSpeed, MinLapNum, MinLapMeanRate, MaxLapMeanRate, vv)
%variable to assign
if (~isempty(cellind)) && (~isempty(behav))
   nspike = numel(pinfo.general.parmfile);
   if (~isfield(pinfo.parm, 'MinLapMeanSpeed')) pinfo.parm.MinLapMeanSpeed = MinLapMeanSpeed*ones(1, nspike); end%%%minimum mean speed of the lap
   if (~isfield(pinfo.parm, 'MinLapNum')) pinfo.parm.MinLapNum = MinLapNum*ones(1, nspike); end %%%minimum number of laps (if good laps smaller than, do not analyze)
   if (~isfield(pinfo.parm, 'MinLapMeanRate')) pinfo.parm.MinLapMeanRate = MinLapMeanRate*ones(1, nspike); end
   if (~isfield(pinfo.parm, 'MaxLapMeanRate')) pinfo.parm.MaxLapMeanRate = MaxLapMeanRate*ones(1, nspike); end
   if (~isfield(pinfo, 'GRfielddynam')) pinfo.GRfielddynam = []; end
   %%%evt field dynam
   if (~isfield(pinfo.GRfielddynam, 'run1DEvtname')) pinfo.GRfielddynam.run1DEvtname = cell(1, nspike); end 
   if (~isfield(pinfo.GRfielddynam, 'runTrajMeanrate')) pinfo.GRfielddynam.runTrajMeanrate = cell(1, nspike); end 
   if (~isfield(pinfo.GRfielddynam, 'runTrajSpatInfo')) pinfo.GRfielddynam.runTrajSpatInfo = cell(1, nspike); end
   if (~isfield(pinfo.GRfielddynam, 'runTrajSpatInfotime')) pinfo.GRfielddynam.runTrajSpatInfotime = cell(1, nspike); end
   if (~isfield(pinfo.GRfielddynam, 'runTrajSparsity')) pinfo.GRfielddynam.runTrajSparsity = cell(1, nspike); end
   if (~isfield(pinfo.GRfielddynam, 'runMeanLapCrr')) pinfo.GRfielddynam.runMeanLapCrr = cell(1, nspike); end
   if (~isfield(pinfo.GRfielddynam, 'runMeanLapActarea')) pinfo.GRfielddynam.runMeanLapActarea = cell(1, nspike); end
   if (~isfield(pinfo.GRfielddynam, 'runMeanLapSparsity')) pinfo.GRfielddynam.runMeanLapSparsity = cell(1, nspike); end
   if (~isfield(pinfo.GRfielddynam, 'runMeanLapSpatInfo')) pinfo.GRfielddynam.runMeanLapSpatInfo = cell(1, nspike); end
   if (~isfield(pinfo.GRfielddynam, 'runMeanLapSpatInfotime')) pinfo.GRfielddynam.runMeanLapSpatInfotime = cell(1, nspike); end
   if (~isfield(pinfo.GRfielddynam, 'runMeanLapActarea')) pinfo.GRfielddynam.runMeanLapActarea = cell(1, nspike); end
   if (~isfield(pinfo.GRfielddynam, 'runMeanLapComX')) pinfo.GRfielddynam.runMeanLapComX = cell(1, nspike); end
   if (~isfield(pinfo.GRfielddynam, 'runMeanLapPeakX')) pinfo.GRfielddynam.runMeanLapPeakX = cell(1, nspike); end
   if (~isfield(pinfo.GRfielddynam, 'runMeanLapMeanRate')) pinfo.GRfielddynam.runMeanLapMeanRate = cell(1, nspike); end
   if (~isfield(pinfo.GRfielddynam, 'runMeanLapPeakRate')) pinfo.GRfielddynam.runMeanLapPeakRate = cell(1, nspike); end

   if (~isfield(pinfo.GRfielddynam, 'runVarLapCrr')) pinfo.GRfielddynam.runVarLapCrr = cell(1, nspike); end
   if (~isfield(pinfo.GRfielddynam, 'runVarLapActarea')) pinfo.GRfielddynam.runVarLapActarea = cell(1, nspike); end
   if (~isfield(pinfo.GRfielddynam, 'runVarLapSparsity')) pinfo.GRfielddynam.runVarLapSparsity = cell(1, nspike); end
   if (~isfield(pinfo.GRfielddynam, 'runVarLapSpatInfo')) pinfo.GRfielddynam.runVarLapSpatInfo = cell(1, nspike); end
   if (~isfield(pinfo.GRfielddynam, 'runVarLapSpatInfotime')) pinfo.GRfielddynam.runVarLapSpatInfotime = cell(1, nspike); end
   if (~isfield(pinfo.GRfielddynam, 'runVarLapActarea')) pinfo.GRfielddynam.runVarLapActarea = cell(1, nspike); end
   if (~isfield(pinfo.GRfielddynam, 'runVarLapComX')) pinfo.GRfielddynam.runVarLapComX = cell(1, nspike); end
   if (~isfield(pinfo.GRfielddynam, 'runVarLapPeakX')) pinfo.GRfielddynam.runVarLapPeakX = cell(1, nspike); end
   if (~isfield(pinfo.GRfielddynam, 'runVarLapMeanRate')) pinfo.GRfielddynam.runVarLapMeanRate = cell(1, nspike); end
   if (~isfield(pinfo.GRfielddynam, 'runVarLapPeakRate')) pinfo.GRfielddynam.runVarLapPeakRate = cell(1, nspike); end

   %%%data variables to assign
   %if (~isfield(data, 'GRfielddynam')) data.GRfielddynam = []; end
   %if (~isfield(data.GRfielddynam, 'evt1DlapRateMaps')) data.field.evt1DlapRateMaps = cell(1, nspike); end
end
for (jjjk = 1:numel(cellind))
    i = cellind(jjjk); 
    %%%%load computing parameters
    %fd1DBaseSess = pinfo.parm.fd1DBaseSess{i}; fd1DBaseLap = pinfo.parm.fd1DBaseLap{i}; 
    smParm.d1sigma = pinfo.parm.fSmooth1DSigma(i); smParm.d1Nsig = pinfo.parm.fSmooth1DNSigma(i); smParm.d1sm = pinfo.parm.fSmooth1D{i};
    smParm.MinLapMeanSpeed = pinfo.parm.MinLapMeanSpeed(i); smParm.MinLapNum = pinfo.parm.MinLapNum(i); 
    smParm.MinLapMeanRate = pinfo.parm.MinLapMeanRate(i); smParm.MaxLapMeanRate = pinfo.parm.MaxLapMeanRate(i); 
    %threshold1Drate = pinfo.parm.f1DThresholdRate(i); max1Dgap = pinfo.parm.f1DMaxGap(i); min1Dpeakrate = pinfo.parm.f1DMinPeakRate(i); 
    timeunit = pinfo.parm.timeunit(i); %base1DrateN = pinfo.parm.f1DBaseRateGridN(i);  
    disp(strcat('-----> compute field dynamics ---', pinfo.general.parmfile{i}));
    spiketime = timeunit*data.spike.spiketime{i};
    %%%%%%%%compute lap-by-lap 1D rate maps on the trajectory the cell's most prominent field (maximum peak rate) resides
    fprate = pinfo.GRfield.PF1DInPeakrate{i}; 
    if (~isempty(fprate)) %%%if fields exist on linear track
       fieldevName = pinfo.GRfield.PF1Devt{i}; 
       evName = pinfo.general.eventname{i}; evTime = data.events.eventtimes{i}; evType = pinfo.parm.eventtype{i}; 
       finaldirnow = pinfo.general.finaldir{i}; nev = 0;
       for (j = 1:numel(evTime))
          trajmeanrate = pinfo.firing.evtmeanrate{i}(j); 
         if (strcmp(evType{j}, 'run')) && (~isempty(find(strcmp(fieldevName,evName{j}))))...
                 && (trajmeanrate >= smParm.MinLapMeanRate) && (trajmeanrate <= smParm.MaxLapMeanRate) %%%if this is a run event on the selected traj 
            %%%%locate event position data
            evSess = identifysession(evTime{j}, pinfo.general.sessionname{i}, pinfo.general.sessionstartT{i}, pinfo.general.sessionendT{i});
            posid = []; evid = [];
            if (~isempty(evSess))
                posid = find( strcmp(behav.general.finaldir, finaldirnow) & strcmp(behav.general.sessname, evSess) );
            end
            if numel(posid == 1)
                if (isfield(behav.general, 'eventname'))
                    evid = find(strcmp(behav.general.eventname{posid}, evName{j}));
                else
                    evid = find(strcmp(behav.behavior.eventname{posid}, evName{j}));
                end
            end
            if (numel(posid)~=1)||(numel(evid)~=1)
               disp(['------------------> Warning: 1D field dynamics not computed for this event: no or more than 1 positon/event files match the session: ', finaldirnow, '___', evName{j}]);
            else
               xbin = bhdata.event.Xbin{posid}{evid}; 
               lapmeanspeed = bhdata.event.LapMeanSpeed{posid}{evid};
               %%%%need to compute lap-by-lap meanrate first
               alllap = numel(lapmeanspeed); lapmeanrate = zeros(1, alllap);
               for (tj = 1:alllap)
                   spikenow = sort( spiketime( (spiketime>=evTime{j}.start(tj)) & (spiketime<=evTime{j}.ent(tj)) ) );
                   lapmeanrate(tj) = numel(spikenow)/ (evTime{j}.ent(tj) - evTime{j}.start(tj) );
               end
               %%%%select the right lap to analyze              
               ttjj = find( (lapmeanspeed >= smParm.MinLapMeanSpeed) & (lapmeanrate >= smParm.MinLapMeanRate) & (lapmeanrate <= smParm.MaxLapMeanRate) );
               nlap = numel(ttjj); 
               if (nlap >= smParm.MinLapNum)
                   nev = nev + 1; 
                   pinfo.GRfielddynam.runTrajSpatInfo{i}(nev) = pinfo.GRfield.runSptlInfo{i}(j);
                   pinfo.GRfielddynam.runTrajSpatInfotime{i}(nev) = pinfo.GRfield.runSptlInfotime{i}(j);
                   pinfo.GRfielddynam.runTrajSparsity{i}(nev) = pinfo.GRfield.runSparsty{i}(j);
                   pinfo.GRfielddynam.runTrajMeanrate{i}(nev) = pinfo.firing.evtmeanrate{i}(j);
                   D1ratenow = cell(1, nlap);
                   spinfo = NaN*ones(1, nlap); spinfotime = NaN*ones(1, nlap); 
                   sparsity = NaN*ones(1, nlap); acarea = NaN*ones(1, nlap);
                   mrate = NaN*ones(1, nlap); prate = NaN*ones(1, nlap); comX = NaN*ones(1, nlap); ploc = NaN*ones(1, nlap); 
                   crr = NaN*ones(1, nlap*(nlap-1));
                   for (tj = 1:nlap)
                       jnow = ttjj(tj);
                       occutime = bhdata.event.LapOccuptime{posid}{evid}{jnow}; framerate = behav.parm.framerate(posid);
                       lappostime{1} = bhdata.event.LapAllPostimestamp{posid}{evid}{jnow}; lapx{1} = bhdata.event.LapAllX{posid}{evid}{jnow};
                       evok.start = evTime{j}.start(jnow); evok.ent = evTime{j}.ent(jnow); 
                       [D1ratenow{tj}, acarea(tj)] = getlinearmap(spiketime, evok, lappostime, lapx, xbin, occutime, smParm, 1/framerate);
                       [spinfo(tj), spinfotime(tj), sparsity(tj)] = findratemapprop(D1ratenow{tj}, occutime);
                       [mrate(tj), prate(tj), ~, comX(tj), ploc(tj), ~, ~] = find1DFFprop(D1ratenow{tj}, xbin); %, min(xbin), max(xbin));
                   end
                   for (ti = 1:nlap)
                       for (tj = ti+1:nlap)
                           [crr((ti-1)*nlap+tj), ~] = correlateratemaps(D1ratenow{ti}, D1ratenow{tj}, xbin, xbin, 0, 0);
                       end
                   end
                   pinfo.GRfielddynam.run1DEvtname{i}{nev} = evName{j};
                   [mcrr, vcrr] = findMV(crr); pinfo.GRfielddynam.runMeanLapCrr{i}(nev) = mcrr; pinfo.GRfielddynam.runVarLapCrr{i}(nev) = vcrr; 
                   [mcrr, vcrr] = findMV(spinfo); pinfo.GRfielddynam.runMeanLapSpatInfo{i}(nev) = mcrr; pinfo.GRfielddynam.runVarLapSpatInfo{i}(nev) = vcrr; 
                   [mcrr, vcrr] = findMV(spinfotime); pinfo.GRfielddynam.runMeanLapSpatInfotime{i}(nev) = mcrr; pinfo.GRfielddynam.runVarLapSpatInfotime{i}(nev) = vcrr; 
                   [mcrr, vcrr] = findMV(acarea); pinfo.GRfielddynam.runMeanLapActarea{i}(nev) = mcrr; pinfo.GRfielddynam.runVarLapActarea{i}(nev) = vcrr; 
                   [mcrr, vcrr] = findMV(sparsity); pinfo.GRfielddynam.runMeanLapSparsity{i}(nev) = mcrr; pinfo.GRfielddynam.runVarLapSparsity{i}(nev) = vcrr; 
                   [mcrr, vcrr] = findMV(mrate); pinfo.GRfielddynam.runMeanLapMeanRate{i}(nev) = mcrr; pinfo.GRfielddynam.runVarLapMeanRate{i}(nev) = vcrr; 
                   [mcrr, vcrr] = findMV(prate); pinfo.GRfielddynam.runMeanLapPeakRate{i}(nev) = mcrr; pinfo.GRfielddynam.runVarLapPeakRate{i}(nev) = vcrr; 
                   [mcrr, vcrr] = findMV(mrate); pinfo.GRfielddynam.runMeanLapMeanRate{i}(nev) = mcrr; pinfo.GRfielddynam.runVarLapMeanRate{i}(nev) = vcrr; 
                   [mcrr, vcrr] = findMV(comX); pinfo.GRfielddynam.runMeanLapComX{i}(nev) = mcrr; pinfo.GRfielddynam.runVarLapComX{i}(nev) = vcrr; 
                   [mcrr, vcrr] = findMV(ploc); pinfo.GRfielddynam.runMeanLapPeakX{i}(nev) = mcrr; pinfo.GRfielddynam.runVarLapPeakX{i}(nev) = vcrr; 
               end
            end
        end
       end
    else
        disp('----------------> 1D dynamics not computed: no place fields found');
    end
end
%%%%%%%%%%%%%%%END OF THE MAIN ROUTINE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mcrr, vcrr] = findMV(crr)
mcrr = NaN; vcrr = NaN;
crr = crr(~isnan(crr)); 
if (~isempty(crr))
    mcrr = mean(crr); vcrr = std(crr);
end

function [mrate, prate, fsize, comX, ploc, skew, kurt] = find1DFFprop(D1map, xxbin) %, fstart, fend)
mrate = NaN; prate = NaN; fsize=NaN; comX=NaN; ploc=NaN; skew=NaN; kurt = NaN; indstart = 1; indend = numel(xxbin);
% for (i = 2:numel(xxbin))
%     if ((xxbin(i)>fstart) && (xxbin(i-1)<=fstart)) 
%        indstart = i; break
%     end
% end
% for (i = 2:numel(xxbin))
%     if ((xxbin(i)>fend) && (xxbin(i-1)<=fend)) 
%        indend = i; break
%     end
% end
if (indend > indstart)
   ratenow = D1map(indstart:indend); xnow = xxbin(indstart:indend);
   iii = find(ratenow > -10); 
   if (numel(iii) > 2)
       ratenow = ratenow(iii); xnow = xnow (iii); mrate = mean(ratenow); [prate,ii] = max(ratenow); ploc = xnow(ii); 
       fsize = sum(ratenow)*mean(diff(xxbin));
       densfun = ratenow/sum(ratenow); comX = xnow * densfun; %locvar = sqrt( ((xnow-comX).^2) * densfun );
       %if (locvar > 0)
       %    X = (xnow-comX)/locvar; skew = (X.^3) * densfun; kurt = (X.^4)*densfun - 3;
       %end
   end
end

function [crr,pp] = correlateratemaps(rate1, rate2, XY1, XY2, RC, theta)
%%%%here RC and theta are the rotation center and angle: RC[x y](coordinates in pixels), 
%%%% theta (angles in degree): 0, 90, -90, 180, 1 (x mirror), -1 (y mirror)  
%%%% assume XYgrid1 and XYgrid2 (could be 1D or 2D) have the same bin sizes
crr = NaN; pp = NaN;
%%%%first transform the rotated rate maps for 2D maps; for 1D maps: only do the x mirror (only choice 1 is valid... wait... look at directionaly algorithm first)
%%%%for 2D maps: transform XY grid coordinates, then choose common grids and their rates
%%%%%%%%%%determine 1D or 2D rate maps
if (~iscell(XY1)) %%if a 1D rate map: add the cell dimension and RC input replaced by traj2 center
    nd = 1; binsize = XY2(2) - XY2(1);
else
    nd = 2; binsize = XY2{1}(2) - XY2{1}(1);
end
%%%%%first, this is the first map rates and coordinates (x,y), & the second map without rotation/mirror, shifting, etc.
if (nd == 1)
    r1 = rate1; x1 = XY1; r2 = rate2; x2 = XY2; RC = (min(XY2) + max(XY2))/2; 
elseif (nd == 2)
    nx1 = numel(XY1{1}); ny1 = numel(XY1{2});
    r1 = reshape(rate1, nx1*ny1, 1); x1 = zeros(nx1*ny1,1); y1 = zeros(nx1*ny1,1); 
    for (i = 1:nx1)
        for (j = 1:ny1)
            x1((i-1)*ny1+j) = XY1{1}(i); y1((i-1)*ny1+j) = XY1{2}(j);
        end
    end
    nx2 = numel(XY2{1}); ny2 = numel(XY2{2});
    r2 = reshape(rate2, nx2*ny2, 1); x2 = zeros(nx2*ny2,1); y2 = zeros(nx2*ny2,1); 
    for (i = 1:nx2)
        for (j = 1:ny2)
            x2((i-1)*ny2+j) = XY2{1}(i); y2((i-1)*ny2+j) = XY2{2}(j);
        end
    end
end
%%%%%%second: work on the second coordinates (rotated, mirrored, shifted, etc)
%%%%rotation or mirror transformation: shift cartesian origin to RC, then rotate in polar coordinates, then shift cartesian orgin back
%%%%%%%%%%%%shift origin
x2 = x2 - RC(1); if (nd == 2) y2 = y2-RC(2); end
%%%%%%%%%%%%rotate or mirror
if (abs(theta) > 2) %if rotate 2D maps
    [ang, rr] = cart2pol(x2, y2); ang = ang + theta*2*pi/360; [x2, y2] = pol2cart(ang, rr); %coordinates for all rate2 points
elseif (theta == 1) %if mirror X
    x2 = -x2;
elseif (theta == -1) %if mirror y: only available for 2D maps
    if (nd ==2) y2 = -y2; end
end
%%%%%%%%%%%%shift origin back
if (nd ==1) x2 = x2 + RC; end
if (nd == 2) x2 = x2+RC(1); y2 = y2+RC(2); end
%%%%match up (r1, x1, y1) with (r2, x2, y2): coordinate (x y) need to be within binsize 
iii = find(r1>-10); r1 = r1(iii); x1 = x1(iii); if (nd ==2) y1 = y1(iii); end
jjj = find(r2>-10); r2 = r2(jjj); x2 = x2(jjj); if (nd ==2) y2 = y2(jjj); end
nmatch = 0; rr1 = []; rr2 = []; %%%%have to work point-by-point
if (nd ==1)
    for (i = 1:numel(x1))
        [mmm, iii] = min( abs(x2-x1(i)) );
        if (mmm <= binsize)
            nmatch = nmatch + 1; rr1(nmatch) = r1(i); rr2(nmatch) = r2(iii);
        end
    end
elseif (nd == 2)
    for (i = 1:numel(x1))
        [mmm, iii] = min( abs(x2-x1(i)) + abs(y2-y1(i)) );
        if (mmm <= 2*binsize)
            nmatch = nmatch + 1; rr1(nmatch) = r1(i); rr2(nmatch) = r2(iii);
        end
    end
end
if (nmatch>=10) %%%%if overlap more than 10 bins
    m1 = mean(rr1); m2 = mean(rr2); s1 = std(rr1); s2 = std(rr2);
    if (s1~=0)&&(s2~=0)
        [R,P] = corrcoef(rr1',rr2'); crr = R(1,2); pp = P(1,2);
        %rr1 = (rr1-m1)/s1; rr2 = (rr2-m2)/s2; crr = rr1*rr2'/nmatch; 
    end
end

function [D1rate, nact] = getlinearmap(spiketime, evTime, lappostime, lapx, xbin, occutime, smParm, frametime)
nx = numel(xbin); D1rate = zeros(nx,1); %D1rate = NaN*ones(nx,1); 
nlap = numel(lappostime); count = zeros(nx,1); spikex = []; nact = NaN;
frametime = 1.03*frametime; %%%slightly expand frametime to assure detection, due to slightly variable frame-to-frame time
if (nx > 1)
    binsize = xbin(2) - xbin(1); sigma = round(smParm.d1sigma/binsize); %%sigma now in number of bins
%%%%%count spikes in xbin lap bu lap
for (i = 1:nlap)
    counti = zeros(nx,1);
    spikenow = sort( spiketime( (spiketime>=evTime.start(i)) & (spiketime<=evTime.ent(i)) ) );
    spikex = NaN*ones(size(spikenow)); [laptimenow, iii] = sort(lappostime{i}); allxnow = lapx{i}(iii);
    lastpoint = 1;  %this only for saving time
    for (j = 1:numel(spikenow))
         for (k = lastpoint:numel(laptimenow)) %find corresponding time in position data
             if abs(laptimenow(k) - spikenow(j)) <= frametime 
             %if (laptimenow(k) <= spikenow(j)) && (laptimenow(k+1) > spikenow(j)) 
                 spikex(j) = allxnow(k); lastpoint = k; break; 
             end
         end
    end
    for (j = 1:nx) counti(j) = numel(find((spikex>=xbin(j)) & (spikex<xbin(j)+binsize))); end
    %counti(nx) = numel(find(spikex>=xbin(nx)));
    count = count + counti;
end
%%%%%compute firing rate
for (i = 1:nx)
     if (occutime(i) ~= 0) D1rate(i) = count(i)/occutime(i); end
end
nact = numel(find(count>0))/numel(find(occutime>0));
%%now do a 1d smoothing
if (strcmpi(smParm.d1sm, 'yes'))
    otime = ones(size(occutime));
    D1rate = OneDSmooth_new(D1rate, otime, sigma, smParm.d1Nsig); %%%smooth first: smoothing does not remove NaN's (occutime = 0)
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

function [spinfo, spinfotime, sparsty] = findratemapprop(rate, occu)
%%%%%%another variable: mutual information is not computed here
spinfo = NaN; spinfotime = NaN; sparsty = NaN; %size(rate)
%size(occu)
[ny, nx] = size(rate); rate = reshape(rate, nx*ny, 1); occu = reshape(occu, nx*ny, 1); mr = max(rate);
if (sum(occu) > 0)
   occu = occu / sum(occu); jjj = find(rate>-1); meanrate = sum(occu(jjj).*rate(jjj));
   spinfo = 0; spinfotime = 0; sparsty = 0;  
   if (meanrate > 0)
       rateok = rate / meanrate; iii = find(rateok > 0); 
       spinfo = sum(occu(iii).*rateok(iii).*log2(rateok(iii))); 
       spinfotime = spinfo*meanrate; sparsty = 1/sum(occu(iii).*rateok(iii).*rateok(iii));
   end
end

function baserate = computebaserate(ratenow, baserategridN)
baserate = 0;
ratenow = sort(ratenow(ratenow>-10)); nn = round(baserategridN*numel(ratenow));
if (nn > 0) baserate = mean(ratenow(1:nn)); end

function [D2refmap, XYrefbin, refsegNumout] = resolve2Drefmap(segSess, segNum, withinsegNum, D2ratenow, XYsess, XYbin, BaseSess, BaseSeg)
D2refmap = []; XYrefbin = []; refsegNumout = NaN;
ises = find( strcmpi(segSess,BaseSess) ); ibin = find( strcmpi(XYsess, BaseSess) );
if isempty(ises) || isempty(ibin)
    disp(['----------------> 2D dynamics not computed: no reference sessions found: ', BaseSess]);
else
    mapnow = D2ratenow(ises); withinsegnow = withinsegNum(ises); segNum = segNum(ises);
    refsegnum = str2num(BaseSeg);  
    if (isempty(refsegnum))
        if (strncmpi(BaseSeg, 'last', 4))
            [str, tok] = strtok(lower(BaseSeg), 'last');
            if (~isempty(str2num(str))) 
                refsegnum = (max(withinsegnow)-(str2num(str)-1)):max(withinsegnow); 
            else
                [str, tok] = strtok(str, '_');
                if (~isempty(str2num(str))) refsegnum = (max(withinsegnow)-(str2num(str)-1)):max(withinsegnow);  end
            end
        elseif (strncmpi(BaseSeg, 'first', 5))
            [str, tok] = strtok(lower(BaseSeg), 'first');
            if (~isempty(str2num(str))) 
                refsegnum = min(withinsegnow):(min(withinsegnow)+(str2num(str)-1)); 
            else
                [str, tok] = strtok(str, '_');
                if (~isempty(str2num(str))) refsegnum = min(withinsegnow):(min(withinsegnow)+(str2num(str)-1));  end
            end 
        end
    end
    if (~isempty(refsegnum)) || (strncmpi(BaseSeg, 'average', 4)) %%%if average or segments specified, do an average map
        if (~isempty(refsegnum)) 
            [~, iw, ir] = intersect( withinsegnow, refsegnum ); mapnow = mapnow(iw); segNum = segNum(iw);
        end
        refsegNumout = segNum;
        if (~isempty(mapnow))
            XYrefbin = XYbin{ibin(1)};
            D2refmap = mapnow{1}; [mm,nn] = size(D2refmap);
            for (i = 1:mm)
                for (j = 1:nn)
                    for (k = 2:numel(mapnow));
                        D2refmap(i,j) = D2refmap(i,j) + mapnow{k}(i,j);
                    end
                    D2refmap(i,j) = D2refmap(i,j)/numel(mapnow);
                end
            end
        else
           disp(['----------------> 2D dynamics not computed: no reference session segments found: ', BaseSess, '__', BaseSeg]);
        end
    elseif (strcmpi(BaseSeg, 'last')) %%if last segment of a session
        XYrefbin = XYbin{ibin(1)}; D2refmap = mapnow{numel(mapnow)}; refsegNumout = segNum(numel(mapnow));
    elseif (strcmpi(BaseSeg, 'first')) %%if first segment of a session
        XYrefbin = XYbin{ibin(1)}; D2refmap = mapnow{1}; refsegNumout = segNum(1);
    else
        disp(['----------------> 2D dynamics not computed: no reference session segments found: ', BaseSess, '__', BaseSeg]);
    end
end

function [D1refmap, xrefbin, reflapNum] = resolve1Drefmap(lapevName, lapNum, withinlapNum, D1ratenow, xev, xbin, BaseSess, BaseSeg)
D1refmap = []; xrefbin = []; reflapNum = NaN; nev = numel(xev); nlap = numel(lapNum); %%%BaseSeg = BaseLap; use the name seg for lap
evid = zeros(1, nev); lapid = zeros(1, nlap);
for (i = 1:nev)
    if (~isempty(strfind(lower(xev{i}), lower(BaseSess)))) evid(i) = 1; end
end
for (i = 1:nlap)
    if (~isempty(strfind(lower(lapevName{i}), lower(BaseSess)))) lapid(i) = 1; end
end
ises = find( lapid ); ibin = find( evid );
if (isempty(ises) || isempty(ibin))
    disp(['----------------> 1D dynamics not computed: no reference sessions found: ', BaseSess]);
else
    mapnow = D1ratenow(ises); withinlapnow = withinlapNum(ises); lapNum = lapNum(ises);
    refsegnum = str2num(BaseSeg);  
    if (isempty(refsegnum))
        if (strncmpi(BaseSeg, 'last', 4))
            [str, tok] = strtok(lower(BaseSeg), 'last');
            if (~isempty(str2num(str))) 
                refsegnum = (max(withinlapnow)-(str2num(str)-1)):max(withinlapnow); 
            else
                [str, tok] = strtok(str, '_');
                if (~isempty(str2num(str))) refsegnum = (max(withinlapnow)-(str2num(str)-1)):max(withinlapnow);  end
            end
        elseif (strncmpi(BaseSeg, 'first', 5))
            [str, tok] = strtok(lower(BaseSeg), 'first');
            if (~isempty(str2num(str))) 
                refsegnum = min(withinlapnow):(min(withinlapnow)+(str2num(str)-1)); 
            else
                [str, tok] = strtok(str, '_');
                if (~isempty(str2num(str))) refsegnum = min(withinlapnow):(min(withinlapnow)+(str2num(str)-1));  end
            end 
        end
    end
    if (~isempty(refsegnum)) || (strncmpi(BaseSeg, 'average', 4)) %%%if average or segments specified, do an average map
        if (~isempty(refsegnum)) 
           [~, iw, ir] = intersect( withinlapnow, refsegnum ); mapnow = mapnow(iw); lapNum = lapNum(iw);
        end
        reflapNum = lapNum;
        if (~isempty(mapnow))
            xrefbin = xbin{ibin(1)}; D1refmap = mapnow{1}; nmap = numel(mapnow);
            for (k = 2:nmap);
                 D1refmap = D1refmap + mapnow{k};
            end
            D1refmap = D1refmap/nmap;
        else
            disp(['----------------> 1D dynamics not computed: no reference session segments found: ', BaseSess, '__', BaseSeg]);
        end
    elseif (strcmpi(BaseSeg, 'last')) %%if last segment of a session
        xrefbin = xbin{ibin(1)}; D1refmap = mapnow{numel(mapnow)}; reflapNum = lapNum(numel(mapnow));
    elseif (strcmpi(BaseSeg, 'first')) %%if first segment of a session
        xrefbin = xbin{ibin(1)}; D1refmap = mapnow{1}; reflapNum = lapNum(1);
    else
        disp(['----------------> 1D dynamics not computed: no reference session segments found: ', BaseSess, '__', BaseSeg]);
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
   
function evSess = identifysession(evTime, sessionname, startT, endT)
evSess = []; minevstart = min(evTime.start); maxevend = max(evTime.ent);
iii = find( (startT <= minevstart) & (endT >= maxevend) );
if (numel(iii) == 1)
    evSess = sessionname{iii};
end   
   
   