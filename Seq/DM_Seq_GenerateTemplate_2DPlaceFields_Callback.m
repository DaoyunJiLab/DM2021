function DM_Seq_GenerateTemplate_2DPlaceFields_Callback
%------------This is for group-generation of templates in 2D spaces %%%%%%%%%%%%%%%%%%
%------------Generate templates on all cells with 2D fields already defined in selected events
%%%%%%%%%%%%%%%%%%%%%%%% Does not make sense for this option - 2D fields are defined for the entire session, not for individual events %%%%%%%%%%%%%%%%%%%%%

%----------     no need to redo the 2D rate maps, just pick from the spikedb database
%------------ If the event can be linearized (behav.parm.eventPosltr exist), linearization joints are loaded, but not used in generating the template (rate maps)
%%%Template defined as spatially binned and smoothed rate maps 
%%%%           spatial bin size is pre-defined in the figure-associated data sructure
%%%output as a .tmp file in the template directory with filename indicating different types of template

disp('-----> This option for template generation is unavailable; aborted'); ok = 0;
if ok
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
    if ~isfield(pinfo, 'field')
        disp('-----> place fields not defined in the database; aborted'); ok = 0;
    else
        if (plotparm.linkbehav == 0)
           disp(['--------> no behav data linked']); ok = 0;
        else
           behav = getappdata(hf, 'behav'); bhdata = getappdata(hf, 'bhdata');
        end
    end
end
end
if ok
    input = inputdlg({'Enter event keyword'; 'Enter event type'}, 'Events for template', 2, {'track'; 'run'}); 
    if (~isempty(input))
        tmpparm.evkeyword = input{1}; tmpparm.evkeytype = input{2};
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
    cellletter = DE_Seq_AssignLetter(numel(cellid)); 
    evType = pinfo.parm.eventtype{cellid(1)}; evName = pinfo.general.eventname{cellid(1)};
    evTimes = data.events.eventtimes{cellid(1)};
    cii = strfind(fdirnow, 'final'); animaldate = fdirnow(1:cii-2);
    ind = find( fmatchevent(evName, evType, tmpparm.evkeyword, tmpparm.evkeytype) );
    for (j = 1:numel(ind))
         evnamenow = evName{ind(j)}; 
         disp(['----------> event now: ', evnamenow]);
         tmpnow = findtmp_singleevt(pinfo, data, behav, bhdata, cellid, cellletter, fdirnow, evnamenow, evTimes{ind(j)}, animaldate, evnamenow);
         if ~isempty(tmpnow)
            nt = nt + 1; tmp(nt) = tmpnow;
         end
    end
end
if (nt > 0) tmpout = tmp; end

function tmp = findtmp_singleevt(pinfo, data, behav, bhdata, cellid, cellletter, fdirnow, evnamenow, evTimes)
tmp = []; ncell = 0; timebin = 0.1; %%%default time binsize for template generation
for (ttjjk = 1:numel(cellid))
    i = cellid(ttjjk);
    allpfid = find(strcmp(pinfo.field.PF1Devt{i}, evnamenow));
    if ~isempty(allpfid)
        %pr = pinfo.field.PF1DInPeakrate{i}(allpfid); %loc = pinfo.field.PF1DLocPeakX{i}(allpfid);
        %[~, iii] = max(pr); 
        ncell = ncell + 1; 
        %peakloc(ncell) = loc(iii); 
        cellokind(ncell) = i; letterok(ncell) = cellletter(ttjjk);
    end
end
if (ncell<4)
    disp('------------> not enough cells with fields on the trajectory.');
else
    %[~, iii] = sort(peakloc); cellokind = cellokind(iii); tmp.tmp.tmpfinalseq{1} = letterok(iii); tmp.tmp.tmpfileletter = letterok(iii);
    filename = pinfo.general.parmfile(cellokind); fileletter = letterok;
    
    %%%%compute spike cnt for template/spatial rates
    spikedata = data.spike.spiketime(cellokind); 
    for (i = 1:numel(cellokind))  spikedata{i} = spikedata{i}*pinfo.parm.timeunit(cellokind(i)); end
    [binsize, framerate, postimestamp, xpos, ypos, joint] = findposdata(pinfo, behav, bhdata, fdirnow, evnamenow, evTimes, cellokind(1));
    [spikecnt, spiketime, tmplength, Xgrid, spikespatialrate, spikeX] = generatespacetmp(spikedata, evTimes, timebin, binsize, postimestamp, xpos, ypos, joint, framerate);
    %for (i = 1:numel(spikecnt)) [~, peakbin(i)] = max(spikecnt{i}); end
    for (i = 1:numel(spikecnt)) [~, peakbin(i)] = max(spikespatialrate{i}); end
    [~, IX] = sort(peakbin); fileletter = fileletter(IX);
    filename = filename(IX); spikecnt = spikecnt(IX); spiketime = spiketime(IX); %plotcolor = plotcolor(IX);
    spikeX = spikeX(IX); spikespatialrate = spikespatialrate(IX); 
    
    %%%%%%%assign the output structure
    tmp.tmp.tmptype = '2Dratemap'; tmp.tmp.tmplength = tmplength;
    tmp.tmp.tmpfilename = filename; tmp.tmp.tmpfileletter = fileletter; tmp.tmp.tmpfinalseq{1} = fileletter; 
    tmp.tmp.tmpspiketime = spiketime; tmp.tmp.tmpspikecnt = spikecnt;
    tmp.tmp.Xgrid = Xgrid; tmp.tmp.tmpspacerate = spikespatialrate; tmp.tmp.spikeX = spikeX;
    tmp.tmp.tmpscore = []; tmp.tmp.tmporderpair = [];  tmp.tmp.crrT = [];
    tmp.parm.animaldate = strcat(pinfo.general.animalname{cellokind(1)}, '_', pinfo.general.datedir{cellokind(1)}); 
    tmp.evtfile = evnamenow;
    tmp.parm.option = 'spatial'; tmp.parm.smooth = pinfo.parm.fSmooth1D{cellokind(1)}; 
    %     %%%%parameters
    tmp.parm.smoothsig = pinfo.parm.fSmooth1DSigma(cellokind(1));
    tmp.parm.epindex = NaN; pinfo.parm.nepdisplay = 1; tmp.parm.starttime = 0; tmp.parm.endtime = 0;
    tmp.parm.timebinsize = timebin; tmp.parm.spatialbinsize = binsize;
    %%%junk parameters
    tmp.parm.peaksearchhalfbinnum = [];
    tmp.parm.lapconsisthres = []; tmp.parm.kurtthres = []; tmp.parm.zscorethres = [];
    tmp.parm.maxlag = []; tmp.parm.crrmode = [];
    tmp.tmp.tmpscore = []; tmp.tmp.tmporderpair = [];  tmp.tmp.crrT = [];
end

function [binsize, framerate, postimestamp, xpos, ypos, joint] = findposdata(pinfo, behav, bhdata, fdirnow, evName, evTimes, i)
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
    disp(['------------------> Warning: lap rate map not computed for this event: no or more than 1 positon/event files match the session: ', finaldirnow, '___', evName{j}]);
else
    postimestamp = bhdata.pos.postimestamp{posid}*behav.parm.timeunit(posid);
    Pmarker = behav.parm.sessPmarker{posid}; allposmarker = behav.general.posMarker{posid}; 
    xpos = []; ypos = []; ik = find(strcmp(allposmarker, Pmarker)); 
    if (numel(ik) == 1) xpos = bhdata.pos.XX{posid}{ik}*behav.parm.pixelXSize(posid); ypos = bhdata.pos.YY{posid}{ik}*behav.parm.pixelYSize(posid); end
    binsize = behav.parm.s1dbinsize(posid); framerate = behav.parm.framerate(posid);
    posltrfile = bhdata.pos.ltrfilename{posid}; posltr = bhdata.pos.posltr{posid}; evPosltr = behav.parm.eventPosltr{posid};
    joint = [];
    for (tt = 1:numel(posltrfile))
        if (strcmp(posltrfile{tt}, evPosltr{evid}))
             joint = posltr{tt}; 
             joint(:,1) = joint(:,1)*behav.parm.pixelXSize(posid); joint(:,2) = joint(:,2)*behav.parm.pixelYSize(posid);
             break
        end
    end
end
function evSess = identifysession(evTime, sessionname, startT, endT)
evSess = []; minevstart = min(evTime.start); maxevend = max(evTime.ent);
iii = find( (startT <= minevstart) & (endT >= maxevend) );
if (numel(iii) == 1)
    evSess = sessionname{iii};
end   

function [spikecnt, spiketime, timelength, grid1D, stspikecnt, xout] = generatespacetmp(data, ep, bin, binsize, postimestamp, xpos, ypos, joint, framerate)
%minimum time resolution for constructing temporal template, should always use this value
smoothparm.smmode = 'yes'; smoothparm.sig = 2; %%%smooth 2 bins away
nspike = numel(data); nep = numel(ep.start);
   %%%compute template spatial profile in a binned grid
   %%%%%%%%%1. compute occupancy time at each spatial grid
   %%%%%%%%%2. get all position points in each grid for each episode
   for (i = 1:nep) %for each episode
          tindex = find( (postimestamp>=ep.start(i)) & (postimestamp<=ep.ent(i)) );
          timestamp{i} = postimestamp(tindex); xeppos{i} = xpos(tindex); yeppos{i} = ypos(tindex);
   end
   for (i = 1:nep) %for each episode
       [posxout{i}, pout] = LinearizePath(joint, xeppos{i}, yeppos{i}, 0);
   end
   %count position points in each grid
   leng = pout(numel(pout)); grid1d = ceil(leng/binsize); %space length and bins
   Xgrid = ((1:grid1d+1)'-1)*binsize;
   Yposcnt = zeros(grid1d+1, 1); %initial assignment of the position point count
   for (k = 1:nep)
       for (i = 1:grid1d)
            cnt = find( (posxout{k}>=Xgrid(i)) & (posxout{k}<Xgrid(i+1)) ); %count points in grid
            if (~isempty(cnt))
                Yposcnt(i)=Yposcnt(i)+numel(cnt); 
            end
            gridepcnt{k}{i} = cnt;
            cnt = [];
       end
       cnt = find(posxout{k}>=grid1d*binsize);
       if (~isempty(cnt))
           Yposcnt(grid1d+1)=Yposcnt(grid1d+1)+numel(cnt); 
       end
       gridepcnt{k}{grid1d+1} = cnt;
       cnt = [];
   end
   %occupancy time
   occuptime = Yposcnt / framerate;  
   %linearize spikes
   for (i = 1:nspike) 
       for (j = 1:nep)
           index = find( (data{i}>=ep.start(j)) & (data{i}<=ep.ent(j)) );
           xout{i}{j} = []; spiketimenow{i}{j} = [];
           if (~isempty(index))
              spiketimenow{i}{j} = data{i}(index); 
                 xout{i}{j} = zeros(1, numel(index));
                 lastpoint = 1;
                 for (k = 1:numel(spiketimenow{i}{j}))
                     if (spiketimenow{i}{j}(k) >= timestamp{j}(numel(timestamp{j})))
                         m = numel(timestamp{j});
                         xout{i}{j}(k) = posxout{j}(m);
                         [ingrid{i}{j}(k), correctratio{i}{j}(k)] = findcorrectratio(m, gridepcnt{j}, posxout{j}(m), Xgrid);
                     elseif (spiketimenow{i}{j}(k) < timestamp{j}(1))
                         m = 1;
                         xout{i}{j}(k) = posxout{j}(m);
                         [ingrid{i}{j}(k), correctratio{i}{j}(k)] = findcorrectratio(m, gridepcnt{j}, posxout{j}(m), Xgrid);
                     else
                         for (m = lastpoint:numel(timestamp{j})-1)
                             if ((timestamp{j}(m) <= spiketimenow{i}{j}(k)) & (timestamp{j}(m+1) > spiketimenow{i}{j}(k))) 
                                xout{i}{j}(k) = posxout{j}(m); 
                                [ingrid{i}{j}(k), correctratio{i}{j}(k)] = findcorrectratio(m, gridepcnt{j}, posxout{j}(m), Xgrid);
                                lastpoint = m;
                                break; 
                             end
                         end
                     end
                     if (xout{i}{j}(k) == 0)
                         disp(['--------> warning: zero position for spike at: ', num2str(spiketimenow{i}{j}(k))]);
                     end
                 end
           end
       end
   end
   %get spatial firing rate and spatial raster
   for (i = 1:nspike)
       stspikecnt{i} = zeros(grid1d+1, 1);
       for (k = 1:grid1d)
           for (j = 1:nep)
               stspikecnt{i}(k) = stspikecnt{i}(k) + numel( find( (xout{i}{j}>=Xgrid(k)) & (xout{i}{j}<Xgrid(k+1)) ) );
           end
           if (occuptime(k) > 0)
              stspikecnt{i}(k) = stspikecnt{i}(k)/occuptime(k); %spikecnt already firing rate
           end
       end
       for (j = 1:nep)
           stspikecnt{i}(grid1d+1) = stspikecnt{i}(grid1d+1) + numel( find(xout{i}{j}>=grid1d*binsize) ); %if spikes at last bin
       end
       if (occuptime(grid1d+1) > 0)
           stspikecnt{i}(grid1d+1) = stspikecnt{i}(grid1d+1)/occuptime(grid1d+1); %spikecnt already firing rate
       end
   end
   if (strncmpi(smoothparm.smmode, 'yes', 1))
       sigma = smoothparm.sig; XX = [-5*sigma:1:5*sigma]; wind = normpdf(XX, 0, sigma); wind = wind/sum(wind); nw = numel(wind); %will have a delay = (n-1)/2 bins
       for (i = 1:numel(stspikecnt)) %smoothing
           Yrate = stspikecnt{i};
           Yrate = conv(Yrate, wind); Yrate = Yrate((nw-1)/2+1:numel(Yrate)-(nw-1)/2);
           stspikecnt{i} = Yrate; %conv(stspikecnt{i}, wind)/sum(wind); stspikecnt{i} = stspikecnt{i}(3:numel(stspikecnt{i})-2);
       end
   end
   %now get spiketimes with respect to average trajectory time
   ctime = occuptime / nep; %occupancy time per run per grid (1 pixel) in second
   timelength = sum(ctime);
   for (i = 1:nspike) %convert from each spike xout{i}{j} to time 
       for (j = 1:nep)
           nnn = numel(xout{i}{j}); spiketime{i}{j} = zeros(1,nnn);
           for (k = 1:nnn)
               ingridnow = ingrid{i}{j}(k);
               spiketime{i}{j}(k) = sum(ctime(1:ingridnow)) - ctime(ingridnow) * correctratio{i}{j}(k);
           end
       end
   end
   %%bin spikes in the time resolution specified
   nbin = ceil(timelength/bin); bintimenow =  0 + ((1:nbin)-1) * bin;
   for (i = 1:nspike)
       spikecnt{i} = zeros(1,nbin+1); %spike count row vector
       for (j = 1:nep)
           for ( kk = 1:nbin)
              index = find( (spiketime{i}{j} >= bintimenow(kk)) & (spiketime{i}{j} < bintimenow(kk)+bin) );
              spikecnt{i}(kk) = spikecnt{i}(kk) + numel(index);
           end
           spikecnt{i}(nbin+1) = spikecnt{i}(nbin+1) + numel(find(spiketime{i}{j} >= nbin*bin));
       end
       spikecnt{i} = spikecnt{i} / (nep*bin);
   end
   %%smooth spike count
   if (strncmpi(smoothparm.smmode, 'yes', 1))
      sigma = smoothparm.sig; XX = [-5*sigma:1:5*sigma]; wind = normpdf(XX, 0, sigma); wind = wind/sum(wind); nw = numel(wind); %will have a delay = (n-1)/2 bins
      for (i = 1:numel(spikecnt)) %smoothing
          Yrate = spikecnt{i};
          Yrate = conv(Yrate, wind); Yrate = Yrate((nw-1)/2+1:numel(Yrate)-(nw-1)/2);
          spikecnt{i} = Yrate; 
      end
   end
grid1D{1} = Xgrid;

function [ingrid, correctratio] = findcorrectratio(m, gridepcnt, posout, Xgrid)
%%find a timing correction ratio for a spike closest to timestamp(m)
%%gridepcnt{ngrid} contains indices in timestamp that all locate in a grid
ingrid = numel(Xgrid); %find which grid in Xgrid posout belongs to
for (i = 2:numel(Xgrid))
    if (posout < Xgrid(i))
        ingrid = i - 1; break
    end
end
if (isempty(gridepcnt{ingrid}))
    %disp(['--------> warning: spike grid index in an empty episode grid']);
    correctratio = 0;
else
    totpoint = numel(gridepcnt{ingrid}); index = find( gridepcnt{ingrid} > m);
    correctratio = numel(index)/totpoint;
end

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



