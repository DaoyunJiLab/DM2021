function DM_Seq_GenerateTemplate_2DActiveCells_Callback
%------------This is for group-generation of templates in 2D spaces %%%%%%%%%%%%%%%%%%
%------------Generate templates on all cells active (removed if not within rate parameters) in selected events
%------------   redo 2D rate maps here to make sure that all cells are included
%NOT DONE YET------------ If the event can be linearized (behav.parm.eventPosltr exist), linearization joints are loaded, but not used in generating the template (rate maps)
%%%Template defined as spatially binned and smoothed rate maps 
%%%%           spatial bin size is pre-defined in the figure-associated data sructure
%%%output as a .tmp file in the template directory with filename indicating different types of template

hf = gcbf; hgroup = getappdata(hf, 'hgroup'); ok = 1;
plotparm = getappdata(hf, 'plotparm'); vv = plotparm.showdetail; %if vv=1, shoe details of computation
cc = plotparm.compute; %if cc=0, assign parameters for computation only; if 1, do the real computation
ow = plotparm.overwrite; %if ow=1, plot into the current database
fname = get(hf, 'Name'); ftt = strfind(fname, '__'); currentfilename = fname(ftt+2:numel(fname));
[~, ~, ee] = fileparts(currentfilename);
if ~strcmp(ee, '.spikedb')
    disp('-----> not a spikedb; aborted'); ok = 0;
else
    pinfo = getappdata(hf, 'pinfo'); data = getappdata(hf, 'data');
    hgroup = getappdata(hf, 'hgroup'); groupselection = getappdata(hgroup, 'selection'); cellind = []; grpind = find(groupselection == 1); 
    for (kk = 1:numel(grpind)) cellind = union(cellind, data.grouplist.groupindex{grpind(kk)}); end
    if (plotparm.linkbehav == 0)
        disp(['--------> no behav data linked']); ok = 0;
    else
        behav = getappdata(hf, 'behav'); bhdata = getappdata(hf, 'bhdata');
    end
end
if ok
    input = inputdlg({'Enter event keyword'; 'Enter event type'; 'Enter event keyNOword'; 'Enter event NOtype'}, 'Events selection', 4, {'track'; 'run'; 'first'; ''}); 
    if (~isempty(input))
        tmpparm.evkeyword = input{1}; tmpparm.evkeytype = input{2}; tmpparm.evkeynoword = input{3}; tmpparm.evkeynotype = input{4}; 
    else
        ok = 0;
    end
end
if ok
    input = inputdlg({'Min rate (Hz)'; 'Max rate'}, 'Cell parameters', 2, {'0.5'; '10'}); 
    if (~isempty(input))
        tmpparm.minrate = str2num(input{1}); tmpparm.maxrate = str2num(input{2});   
    else
        ok = 0;
    end
end
if ok
    input = inputdlg({'Smooth 2D rate maps?'; 'Sigma (bins)'; '# of sigmas'}, 'Smoothing parameters', 3, {'Yes'; '3'; '5'}); 
    if (~isempty(input))
        tmpparm.smoothparm.smmode = input{1}; tmpparm.smoothparm.sigma = str2num(input{2}); tmpparm.smoothparm.Nsigma = str2num(input{3});  
    else
        ok = 0;
    end
end
if ok 
    [fname, pname] = uiputfile(fullfile(cd, '.tmp'), 'Write templates to:');
    if (numel(fname)>1)
       writefilename = fullfile(pname, fname);
    else
       ok = 0;
    end
end
if ok
    tmp = findalltemplates(pinfo, data, behav, bhdata, cellind, tmpparm); 
    if (~isempty(tmp))
       pinfo = []; pinfo = tmp; save(writefilename, 'pinfo', '-mat');
    end
    disp(['-----> number of templates generated: ', num2str(numel(tmp))]);
end
disp('**********************');

function tmpout = findalltemplates(pinfo, data, behav, bhdata, cellind, tmpparm)
%get selected cellind
tmpout = []; nt = 0; 
%%%%search for all the possible finaldir
allfinaldir = unique(pinfo.general.finaldir(cellind));
disp(['-----> number of final dirs: ', num2str(numel(allfinaldir))]);
for (i = 1:numel(allfinaldir))
    fdirnow = allfinaldir{i};
    disp(['-------> final dir now: ', fdirnow]);
    ind = find(strcmp(pinfo.general.finaldir(cellind), allfinaldir{i}));
    cellid = cellind(ind); %%%%all cells in this final directory
    cellletter = DE_Seq_AssignLetter(numel(cellid));  %%This is a temporary assignment so that all cells within a day (for all events) remain the same
    evType = pinfo.parm.eventtype{cellid(1)}; evName = pinfo.general.eventname{cellid(1)};
    evTimes = data.events.eventtimes{cellid(1)};
    indok = find( ifmatchevent(evName, evType, tmpparm.evkeyword, tmpparm.evkeytype, tmpparm.evkeynoword, tmpparm.evkeynotype) );
    for (j = 1:numel(indok))
         evnamenow = evName{indok(j)}; evTimenow = evTimes{indok(j)}; evindnow = indok(j);
         disp(['----------> event now: ', evnamenow]);
         tmpnow = findtmp_singleevt(pinfo, data, behav, bhdata, cellid, cellletter, fdirnow, evnamenow, evTimenow, evindnow, tmpparm);
         if ~isempty(tmpnow)
            nt = nt + 1; tmp(nt) = tmpnow;
         end
    end
end
if (nt > 0) tmpout = tmp; end

function tmp = findtmp_singleevt(pinfo, data, behav, bhdata, cellid, cellletter, fdirnow, evnamenow, evTimes, evindnow, tmpparm)
tmp = []; ncell = 0;
for (ttjjk = 1:numel(cellid))
    i = cellid(ttjjk); meanrate = pinfo.firing.evtmeanrate{i}(evindnow); 
    if (meanrate>=tmpparm.minrate) && (meanrate<=tmpparm.maxrate) %%% select cells active in this event
        ncell = ncell + 1; 
        cellokind(ncell) = i; letterok(ncell) = cellletter(ttjjk);
    end
end
if (ncell<4)
    disp('------------> not enough active cells in the event');
else
    filename = pinfo.general.parmfile(cellokind); fileletter = letterok;
    %%%%compute spike cnt for template/spatial rates
    spikedata = data.spike.spiketime(cellokind); 
    for (i = 1:numel(cellokind))  spikedata{i} = spikedata{i}*pinfo.parm.timeunit(cellokind(i)); end
    [binsize, framerate, postimestamp, xpos, ypos] = findposdata(pinfo, behav, bhdata, fdirnow, evnamenow, evTimes, cellokind(1));
    [xgrid, ygrid, ratemap] = generate2Dspacetmp(spikedata, evTimes, binsize, postimestamp, xpos, ypos, framerate, tmpparm.smoothparm);
    %%%%%%%%%% Do a 2D sorting: peak location x + y along diagonal line
    peakx = zeros(1, ncell); peaky = zeros(1, ncell);
    for (i = 1:numel(ratemap)) [peakx(i), peaky(i)] = findpeakxyind(ratemap{i}); end
    [~, IX] = sort(peakx + peaky); fileletter = fileletter(IX);
    filename = filename(IX); ratemap = ratemap(IX);  
    %%%%%%%assign the output structure
    tmp.tmp.tmptype = '2Dratemap'; tmp.tmp.tmplength = NaN;
    tmp.tmp.tmpfilename = filename; tmp.tmp.tmpfileletter = fileletter; tmp.tmp.tmpfinalseq{1} = fileletter; 
    tmp.tmp.xgrid = xgrid; tmp.tmp.ygrid = ygrid; tmp.tmp.ratemap = ratemap;
    %%%others
    tmp.tmp.tmpscore = []; tmp.tmp.tmporderpair = [];  tmp.tmp.crrT = [];
    tmp.parm.animaldate = strcat(pinfo.general.animalname{cellokind(1)}, '_', pinfo.general.datedir{cellokind(1)}); 
    tmp.evtfile = evnamenow; tmp.parm.option = 'spatial'; 
    tmp.parm.smooth = tmpparm.smoothparm.smmode; tmp.parm.smoothsig = tmpparm.smoothparm.sigma; tmp.parm.smoothNsig = tmpparm.smoothparm.Nsigma;
    tmp.parm.spatialbinsize = binsize;
    tmp.parm.minrate = tmpparm.minrate; tmp.parm.maxrate = tmpparm.maxrate;
    tmp.parm.evkeyword = tmpparm.evkeyword; tmp.parm.evkeytype = tmpparm.evkeytype;
    tmp.parm.spatialunit = 'cm'; tmp.parm.epindex = []; 
    %%%junk parameters
    tmp.parm.peaksearchhalfbinnum = []; tmp.parm.epindex = NaN; tmp.parm.starttime = 0; tmp.parm.endtime = 0;
    tmp.parm.lapconsisthres = []; tmp.parm.kurtthres = []; tmp.parm.zscorethres = [];
    tmp.parm.maxlag = []; tmp.parm.crrmode = []; tmp.tmp.tmpscore = []; tmp.tmp.tmporderpair = [];  tmp.tmp.crrT = [];
end

function [binsize, framerate, postimestamp, xpos, ypos] = findposdata(pinfo, behav, bhdata, fdirnow, evName, evTimes, i)
posid = []; evid = []; binsize = NaN; postimestamp = []; xpos = []; ypos = []; joint = [];
evSess = identifysession(evTimes, pinfo.general.sessionname{i}, pinfo.general.sessionstartT{i}, pinfo.general.sessionendT{i});
posid = find( strcmp(behav.general.finaldir, fdirnow) & strcmp(behav.general.sessname, evSess) );
if numel(posid) == 1
   if (isfield(behav.general, 'eventname'))
       evid = find(strcmp(behav.general.eventname{posid}, evName));
   else
       evid = find(strcmp(behav.behavior.eventname{posid}, evName));
   end
end
if (numel(posid)~=1)||(numel(evid)~=1)
    disp(['------------------> Warning: lap rate map not computed for this event: no or more than 1 positon/event files match the session: ', finaldirnow, '___', evName]);
else
    postimestamp = bhdata.pos.postimestamp{posid}*behav.parm.timeunit(posid);
    Pmarker = behav.parm.sessPmarker{posid}; allposmarker = behav.general.posMarker{posid}; 
    xpos = []; ypos = []; ik = find(strcmp(allposmarker, Pmarker)); 
    if (numel(ik) == 1) xpos = bhdata.pos.XX{posid}{ik}*behav.parm.pixelXSize(posid); ypos = bhdata.pos.YY{posid}{ik}*behav.parm.pixelYSize(posid); end
    binsize = behav.parm.s2dbinsize(posid); framerate = behav.parm.framerate(posid);
end
function evSess = identifysession(evTime, sessionname, startT, endT)
evSess = []; minevstart = min(evTime.start); maxevend = max(evTime.ent);
iii = find( (startT <= minevstart) & (endT >= maxevend) );
if (numel(iii) == 1)
    evSess = sessionname{iii};
end   

function [xp, yp, ratemap] = generate2Dspacetmp(data, ep, gridbin, postimestamp, xpos, ypos, framerate, smoothparm)
ncell = numel(data); ratemap = cell(1, ncell); sigma = round(smoothparm.sigma); Nsigma = round(smoothparm.Nsigma); 
%%%%spatial grids
minx = 0.9*min(xpos); maxx = 1.1*max(xpos); miny = 0.9*min(ypos); maxy = 1.1*max(ypos);
xp = minx:gridbin:maxx; yp = miny:gridbin:maxy;
gridnumx = numel(xp); gridnumy = numel(yp);
%%%get all position points in events
indall = [];
for (m = 1:numel(ep.start))
    tindex = find( (postimestamp>=ep.start(m)) & (postimestamp<=ep.ent(m)) );
    indall = union(indall, tindex);
end
postime = postimestamp(indall); xeppos = xpos(indall); yeppos = ypos(indall);
[postime, iii] = sort(postime); xeppos = xeppos(iii); yeppos = yeppos(iii);
%%%occupancy in spatial grids
v= 0; occupoint = OccupancyPoint(xp, yp, xeppos, yeppos,v); %output occupoint(1:gridnumy,1:gridnumx), how many position points in a grid
occutime = occupoint/framerate; % * timeunit(3) * (postimestamp(2)-postimestamp(1));
%%%%smooth occupancymaps - use 1/4 sigma of ratemaps
if strncmpi(smoothparm.smmode, 'y',1)
   sigmanow = round(sigma/4); if (sigmanow<1) sigmanow =1; end
   smoccutime = TwoDSmooth_new(occutime, ones(size(occutime)), sigmanow, sigmanow, 1, 1); %%%%%%occutime smoothing less than rates
end
%%%%%%for each cell, compute rate within grids
for (k = 1:ncell)
    rate = NaN*ones(gridnumy, gridnumx); %%%for unvisited grids, rate = NaN;
    %%%%filter through events
    indall = [];
    for (m = 1:numel(ep.start))
        tindex = find( (data{k}>=ep.start(m)) & (data{k}<=ep.ent(m)) );
        indall = union(indall, tindex);
    end
    spiketimenow = data{k}(indall); spiketimenow = sort(spiketimenow);
    %%%%work out spike positions
    [spikex, spikey] = findspikepos(postime, xeppos, yeppos, spiketimenow);
    %%%%number of spikes in grids
    v= 0; spikecount = OccupancyPoint(xp, yp, spikex, spikey,v); %output occuspikes(1:gridnumy,1:gridnumx), how many spikes in a grid
    %%%%compute rates
    for (j = 1:gridnumy)
    for (i = 1:gridnumx)
        if (occutime(j,i) ~= 0) 
           rate(j,i) = spikecount(j,i)/occutime(j,i); %%%This makes sure true rate before smoothing is accurate
        elseif (smoccutime(j,i) >= 0) %%%%reset rate = 0 if smoocutime ~= 0 but occutime = 0
           rate(j,i) = 0; 
        end
    end
    end
    %%%%smooth ratemaps
    if strncmpi(smoothparm.smmode, 'y',1)
       rate = TwoDSmooth_new(rate, smoccutime, sigma, sigma, Nsigma, Nsigma); %%%%%%this is the main smoothing engine
    end
    ratemap{k} = rate;
end
function [x,y] = findspikepos(postime, xpos, ypos, spiketime) 
nspike = numel(spiketime); lastpoint = 1;  %this only for saving time
x = NaN*ones(1,nspike); y = NaN*ones(1,nspike);
for (j = 1:nspike)
    for (k = lastpoint: numel(xpos)-1) %find corresponding time in position data
         if (postime(k) <= spiketime(j)) && (postime(k+1) > spiketime(j)) 
             x(j) = xpos(k); y(j) = ypos(k); 
             lastpoint = k;
             break; 
         end
    end
end

function [peakx, peaky] = findpeakxyind(ratemap)
[pprate, IY] = max(ratemap); [~, IX] = max(pprate); 
peakx = IX; peaky = IY(IX); %%%For ratemaps, x = column, y = row 

function evsel = ifmatchevent(evName, evType, evkeyword, evkeytype, evkeynoword, evkeynotype)
%%%decide if a bunch of event evName{i}/evType{i} match evkeyword/evkeytype
evsel = ones(1,numel(evName));
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



