function DM_Seq_GenerateTemplate_PlaceFields_Callback
%------------This is for group-generation of linear track place field templates%%%%%%%%%%%%%%%%%%
%------------   redo the 1D rate maps here to extract temporal length of the templates
%------------Now: generate templates on all selected cells

%%%Generate template based on (lap-averaged) place fields (on linear track) on a spikedb 
%%%template defined as binned spike count of multiple spike trains
%%%%bin size is pre-defined in the figure-associated data sructure
%%%output as a .tmp file in the template directory with filename indicating
%%%   different types of template

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
if ok
    input = inputdlg({'Enter event keyword'; 'Enter event type'; 'Enter event keyNOword'; 'Enter event NOtype'; 'Time binsize (s)'}, 'Binsize/events for template', 5, {'track'; 'run'; 'first'; ''; '0.1'}); 
    if (~isempty(input))
        tmpparm.evkeyword = input{1}; tmpparm.evkeytype = input{2}; tmpparm.evkeynoword = input{3}; tmpparm.evkeynotype = input{4};
        tmpparm.timebin = str2num(input{5}); %%%time bin for translating spatial templates to temporal templates
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
    indok = find( ifmatchevent(evName, evType, tmpparm.evkeyword, tmpparm.evkeytype, tmpparm.evkeynoword, tmpparm.evkeynotype) );
    for (j = 1:numel(indok))
            evnamenow = evName{indok(j)}; 
            disp(['----------> event now: ', evnamenow]);
            tmpnow = findtmp_singleevt(pinfo, data, behav, bhdata, cellid, cellletter, fdirnow, evnamenow, evTimes{indok(j)}, tmpparm);
            if ~isempty(tmpnow)
               nt = nt + 1; tmp(nt) = tmpnow;
            end
   end
end
if (nt > 0) tmpout = tmp; end

function tmp = findtmp_singleevt(pinfo, data, behav, bhdata, cellid, cellletter, fdirnow, evnamenow, evTimes, tmpparm)
tmp = []; ncell = 0; timebin = tmpparm.timebin; %%%default time binsize for template generation
for (ttjjk = 1:numel(cellid))
    i = cellid(ttjjk);
    allpfid = find(strcmp(pinfo.field.PF1Devt{i}, evnamenow)); %%% select cells that have fields in this event
    if ~isempty(allpfid)
        ncell = ncell + 1; 
        cellokind(ncell) = i; letterok(ncell) = cellletter(ttjjk);
    end
end
if (ncell<4)
    disp('------------> not enough cells with fields on the trajectory.');
else
    filename = pinfo.general.parmfile(cellokind); fileletter = letterok;
    %%%%compute spike cnt for template/spatial rates
    spikedata = data.spike.spiketime(cellokind); 
    for (i = 1:numel(cellokind))  spikedata{i} = spikedata{i}*pinfo.parm.timeunit(cellokind(i)); end
    [binsize, Xgrid, occuptime, framerate, lappostimestamp, lapx] = findposdata(pinfo, behav, bhdata, fdirnow, evnamenow, evTimes, cellokind(1));
    smoothparm.smmode = pinfo.parm.fSmooth1D{cellokind(1)}; smoothparm.sigma = pinfo.parm.fSmooth1DSigma(cellokind(1)); smoothparm.Nsigma = pinfo.parm.fSmooth1DNSigma(cellokind(1));
    [spikecnt, spiketime, tmplength, Xgrid, spikespatialrate, spikeX] = generatespacetmp(spikedata, evTimes, timebin, binsize, Xgrid, occuptime, lappostimestamp, lapx, framerate, smoothparm);
    for (i = 1:numel(spikecnt)) [~, peakbin(i)] = max(spikespatialrate{i}); end
    [~, IX] = sort(peakbin); fileletter = fileletter(IX);
    filename = filename(IX); spikecnt = spikecnt(IX); spiketime = spiketime(IX); 
    spikeX = spikeX(IX); spikespatialrate = spikespatialrate(IX); 
    
    %%%%%%%assign the output structure
    tmp.tmp.tmptype = 'pfratetmp'; tmp.tmp.tmplength = tmplength; %%% This is length in time
    tmp.tmp.tmpfilename = filename; tmp.tmp.tmpfileletter = fileletter; tmp.tmp.tmpfinalseq{1} = fileletter; 
    tmp.tmp.tmpspiketime = spiketime; tmp.tmp.tmpspikecnt = spikecnt;
    tmp.tmp.Xgrid = Xgrid; tmp.tmp.tmpspacerate = spikespatialrate; tmp.tmp.spikeX = spikeX;
    tmp.tmp.tmpscore = []; tmp.tmp.tmporderpair = [];  tmp.tmp.crrT = [];
    tmp.parm.animaldate = strcat(pinfo.general.animalname{cellokind(1)}, '_', pinfo.general.datedir{cellokind(1)}); 
    tmp.evtfile = evnamenow;
    tmp.parm.option = 'spatial';  
    %     %%%%parameters
    tmp.parm.smooth = smoothparm.smmode; tmp.parm.smoothsig = smoothparm.sigma; tmp.parm.smoothNsig = smoothparm.Nsigma;
    tmp.parm.epindex = NaN; pinfo.parm.nepdisplay = 1; tmp.parm.starttime = 0; tmp.parm.endtime = 0;
    tmp.parm.timebinsize = timebin; tmp.parm.spatialbinsize = binsize; tmp.parm.evkeyword = tmpparm.evkeyword; tmp.parm.evkeytype = tmpparm.evkeytype;
    %%%junk parameters
    tmp.parm.peaksearchhalfbinnum = [];
    tmp.parm.lapconsisthres = []; tmp.parm.kurtthres = []; tmp.parm.zscorethres = [];
    tmp.parm.maxlag = []; tmp.parm.crrmode = [];
    tmp.tmp.tmpscore = []; tmp.tmp.tmporderpair = [];  tmp.tmp.crrT = [];
end

function [binsize, Xgrid, occuptime, framerate, lappostimestamp, lapx] = findposdata(pinfo, behav, bhdata, fdirnow, evName, evTimes, i)
posid = []; evid = []; binsize = NaN; postimestamp = []; xpos = []; ypos = []; joint = []; Xgrid = [];
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
    lappostimestamp = bhdata.event.LapAllPostimestamp{posid}{evid}; lapx = bhdata.event.LapAllX{posid}{evid};
    binsize = behav.parm.s1dbinsize(posid); framerate = behav.parm.framerate(posid); 
    Xgrid = bhdata.event.Xbin{posid}{evid}; occuptime = bhdata.event.Occuptime{posid}{evid}; 
end
function evSess = identifysession(evTime, sessionname, startT, endT)
evSess = []; minevstart = min(evTime.start); maxevend = max(evTime.ent);
iii = find( (startT <= minevstart) & (endT >= maxevend) );
if (numel(iii) == 1)
    evSess = sessionname{iii};
end   

function [spikecnt, spiketime, timelength, grid1D, stspikecnt, xout] = generatespacetmp(data, ep, bin, binsize, Xgrid, occuptime, lappostimestamp, lapx, framerate, smoothparm)
%minimum time resolution for constructing temporal template, should always use this value
nspike = numel(data); nep = numel(ep.start); timestamp = lappostimestamp; posxout = lapx; 
frametime = 1/framerate; frametime = 1.03*frametime; %%% to make the detection more reliable
   %count position points in each grid
   grid1d = numel(Xgrid);
   Yposcnt = zeros(grid1d, 1); %initial assignment of the position point count
   for (k = 1:nep)
       for (i = 1:grid1d)
            cnt = find( (posxout{k}>=Xgrid(i)) & (posxout{k}<Xgrid(i)+binsize) ); %count points in grid
            if (~isempty(cnt))
                Yposcnt(i)=Yposcnt(i)+numel(cnt); 
            end
            gridepcnt{k}{i} = cnt;
            cnt = [];
       end
   end
%%%%%%%%%%%Use the behavdb occutime to trunk trajectories
%    occuptime = Yposcnt / framerate;  
   %linearize spikes
   for (i = 1:nspike) 
       for (j = 1:nep)
           index = find( (data{i}>=ep.start(j)) & (data{i}<=ep.ent(j)) );
            spiketimenow{i}{j} = sort(data{i}(index)); xout{i}{j} = NaN*ones(1, numel(index));
            ingrid{i}{j} = zeros(1, numel(index)); correctratio{i}{j} = zeros(1, numel(index));
           if (~isempty(index))
                 lastpoint = 1;
                 for (k = 1:numel(spiketimenow{i}{j}))
                     if (spiketimenow{i}{j}(k) > timestamp{j}(numel(timestamp{j}))+frametime)
                         m = numel(timestamp{j});
                         xout{i}{j}(k) = NaN; %posxout{j}(m);
                         [ingrid{i}{j}(k), correctratio{i}{j}(k)] = findcorrectratio(m, gridepcnt{j}, posxout{j}(m), Xgrid);
                     elseif (spiketimenow{i}{j}(k) < timestamp{j}(1)-frametime)
                         m = 1;
                         xout{i}{j}(k) = NaN; %posxout{j}(m);
                         [ingrid{i}{j}(k), correctratio{i}{j}(k)] = findcorrectratio(m, gridepcnt{j}, posxout{j}(m), Xgrid);
                     else
                         for (m = lastpoint:numel(timestamp{j}))
                             if abs(timestamp{j}(m) - spiketimenow{i}{j}(k)) <= frametime 
                             %elseif ((timestamp{j}(m) <= spiketimenow{i}{j}(k)) && (timestamp{j}(m+1) > spiketimenow{i}{j}(k))) 
                                xout{i}{j}(k) = posxout{j}(m); 
                                [ingrid{i}{j}(k), correctratio{i}{j}(k)] = findcorrectratio(m, gridepcnt{j}, posxout{j}(m), Xgrid);
                                lastpoint = m;
                                break; 
                             end
                         end
                     end
                     if isnan(xout{i}{j}(k))
                         disp(['--------> warning: undetermined position for spike at: ', num2str(spiketimenow{i}{j}(k))]);
                     end
                 end
           end
       end
   end
   %get spatial firing rate and spatial raster
   for (i = 1:nspike)
       stspikecnt{i} = zeros(grid1d, 1);
       for (k = 1:grid1d)
           for (j = 1:nep)
               stspikecnt{i}(k) = stspikecnt{i}(k) + numel( find( (xout{i}{j}>=Xgrid(k)) & (xout{i}{j}<Xgrid(k)+binsize) ) );
           end
           if (occuptime(k) > 0)
              stspikecnt{i}(k) = stspikecnt{i}(k)/occuptime(k); %spikecnt already firing rate
           end
       end
   end
   if (strncmpi(smoothparm.smmode, 'yes', 1))
       sigma = round(smoothparm.sigma/binsize); Nsig = smoothparm.Nsigma; 
       for (i = 1:nspike)
           stspikecnt{i} = smoothinterp(stspikecnt{i}, occuptime, sigma, Nsig);
       end
   end
   %now get spiketimes with respect to average trajectory time
   ctime = occuptime / nep; %occupancy time per run per grid (1 pixel) in second
   timelength = sum(ctime); %%%this is length in time
   for (i = 1:nspike) %convert from each spike xout{i}{j} to time 
       for (j = 1:nep)
           nnn = numel(xout{i}{j}); spiketime{i}{j} = NaN*ones(1,nnn);
           for (k = 1:nnn)
               ingridnow = ingrid{i}{j}(k); 
               if ingridnow>0
                  spiketime{i}{j}(k) = sum(ctime(1:ingridnow)) - ctime(ingridnow) * correctratio{i}{j}(k);
               end
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
      sigma = round(smoothparm.sigma/binsize); Nsig = smoothparm.Nsigma; 
        for (i = 1:nspike)
      spikecnt{i} = smoothinterp(spikecnt{i}, ones(size(spikecnt{i})), sigma, Nsig);
        end
   end
grid1D{1} = Xgrid;

function D1rate = smoothinterp(D1rate, occutime, sigma, Nsig)
otime = ones(size(occutime));
D1rate = OneDSmooth_new(D1rate, otime, sigma, Nsig); %%%smooth first: smoothing does not remove NaN's (occutime = 0)
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
evsel = zeros(1, numel(evName));
for (i = 1:numel(evName))
    if (~isempty(evkeytype))||(~isempty(evkeyword)) %%%for inclusion
      if isempty(evkeytype)
         if (~isempty(strfind(lower(evName{i}), lower(evkeyword)))) evsel(i) = 1; end 
      elseif isempty(evkeyword)
         if strncmpi(evType{i}, evkeytype, 3) evsel(i) = 1; end 
      else
         if (~isempty(strfind(lower(evName{i}), lower(evkeyword)))) && strncmpi(evType{i}, evkeytype, 3)
             evsel(i) = 1; 
         end 
      end
    end
    if (~isempty(evkeynotype))||(~isempty(evkeynoword)) %%%for exclusion
       if isempty(evkeynotype)
          if (~isempty(strfind(lower(evName{i}), lower(evkeynoword)))) evsel(i) = 0; end 
       elseif isempty(evkeynoword)
          if strncmpi(evType{i}, evkeynotype, 3) evsel(i) = 0; end
       else
          if (~isempty(strfind(lower(evName{i}), lower(evkeynoword)))) && strncmpi(evType{i}, evkeynotype, 3)
              evsel(i) = 0; 
          end   
       end
    end
end



