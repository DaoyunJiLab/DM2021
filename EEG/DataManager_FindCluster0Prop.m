function [eeg, eegdata] = DataManager_FindCluster0Prop(eeg, eegdata, eegind, vv)

%%%%first add filtered cluster0-band EEG traces
iii = find( strcmp(eeg.parm.band(eegind), 'cluster0') );
SS = questdlg(['Found ', num2str(numel(iii)), ' cluster0 files in the database. Add/modify session cluster0?']);
ok = 1;
if (strcmp(SS, 'Yes'))
    %%%%%search cluster0 files: collapse all cluster0s of a sessions and an area in to single cluster0 file;
    %%%%%                       then normalize (according to the selected mode)
    iii = find( (strcmp(eeg.parm.band, 'broad') & ~strcmp(eeg.general.recarea, 'MSC')) | strcmp(eeg.parm.band, 'cluster0') ); %%%for each broad band eeg file, create a cluster0.spm file
    iii = intersect(iii, eegind);
    for (i = 1:numel(iii))
        sessid = iii(i); infix = 'cluster0';
        %recareanow = eeg.general.recarea{sessid}(1:2); sessnow = eeg.general.sessname{sessid};
        fdirnow = eeg.general.finaldir{sessid}; animnamenow = eeg.general.animalname{sessid};
        %cl0filenow = fullfile(fdirnow, 'spikedata', ['cluster0_', animnamenow, '_', recareanow, '_', sessnow, '.spm']);
        pp1{1} = strcat(fdirnow, filesep, 'spikedata');
        [fff, nnn] =GetAllFile(pp1, 'cluster0', '.spm');
        if isempty(fff) %%%if cluster0 files not computed yet, create them now
            disp(['--------> create cluster0 files in directory: ', fdirnow]);
            Spike_CreateTotalCluster0(fdirnow, 1); %%%%this create cluster0 files for all the animals, areas and sessions
        end
    end
    nnfile = numel(eeg.general.eegfile);
    extra.ss = cell(1, nnfile); extra.mm = cell(1, nnfile); 
    extra.ma = cell(1, nnfile); extra.mi= cell(1, nnfile); %%%extra parameters to carry over from here to later cluster0 analysis
    if isfield(eeg, 'cluster0')
       if isfield(eeg.cluster0, 'sessCntSTD') extra.ss = eeg.cluster0.sessCntSTD; end
       if isfield(eeg.cluster0, 'sessCntMean') extra.mm = eeg.cluster0.sessCntMean; end
       if isfield(eeg.cluster0, 'sessCntMax') extra.ma = eeg.cluster0.sessCntMax; end
       if isfield(eeg.cluster0, 'sessCntMin') extra.mi = eeg.cluster0.sessCntMin; end
    end
    for (i = 1:numel(iii))
        sessid = iii(i); infix = 'cluster0';
        recareanow = eeg.general.recarea{sessid}; sessnow = eeg.general.sessname{sessid};
        fdirnow = eeg.general.finaldir{sessid}; animnamenow = eeg.general.animalname{sessid};
        cl0filenow = fullfile(fdirnow, 'spikedata', ['cluster0_', animnamenow, '_', recareanow(1:2), '_', sessnow, '.spm']);
        if (exist(cl0filenow, 'file') ~= 2)
            cl0filenow = fullfile(fdirnow, 'spikedata', ['cluster0_', animnamenow, '_', recareanow, '_', sessnow, '.spm']);
            if (exist(cl0filenow, 'file') ~= 2)
                if strncmpi(recareanow, 'ctx', 2) || strncmpi(recareanow, 'V', 1)
                    cl0filenow = fullfile(fdirnow, 'spikedata', ['cluster0_', animnamenow, '_', 'V1', '_', sessnow, '.spm']);
                    if (exist(cl0filenow, 'file') ~= 2)
                        cl0filenow = fullfile(fdirnow, 'spikedata', ['cluster0_', animnamenow, '_', 'V2', '_', sessnow, '.spm']);
                        if (exist(cl0filenow, 'file') ~= 2)
                            cl0filenow = fullfile(fdirnow, 'spikedata', ['cluster0_', animnamenow, '_', 'V', '_', sessnow, '.spm']);
                        end
                    end
                end
            end
        end
        if (exist(cl0filenow, 'file') == 2)
           sptime = ReadSpikeTime(cl0filenow);
           %%%%%compute power file 
           disp(['--------> compute cluster0 power of file: ', cl0filenow]);
           starttime = eeg.general.sessstartT{sessid}; endtime = eeg.general.sessendT{sessid};
           parm.smooth = eeg.parm.cl0Smooth{sessid}; parm.normmode = eeg.parm.cl0NormMode{sessid}; parm.normmeth = eeg.parm.cl0NormMeth{sessid};
           parm.binsize = eeg.parm.cl0Binsize(sessid); parm.sigma = eeg.parm.cl0Sigma(sessid); parm.Nsigma = eeg.parm.cl0NSigma(sessid);
           [timepoint, cnt, newmode, ss, mm, ma, mi] = computecl0pwr(sptime, starttime, endtime, parm, eeg, sessid);
           [eeg, eegdata, extra, eegind] = appendcl0file(eeg, eegdata, sessid, cl0filenow, timepoint, cnt, newmode, extra, eegind, ss, mm, ma, mi);
        else
            disp(['--------> Warning: cluster0 .spm file not found: ', cl0filenow]);
        end
    end
    %eegind = 1:numel(eeg.general.eegfile);
elseif (strcmp(SS, 'Cancel'))
    ok = 0; disp(['--------------> action cancelles']);
else
    %%%need to re-asign extra
    nnfile = numel(eeg.general.eegfile);
    extra.ss = cell(1, nnfile); extra.mm = cell(1, nnfile); extra.ma = cell(1, nnfile); %%%extra parameters to carry over from here to later cluster0 analysis
    if isfield(eeg, 'cluster0')
       if isfield(eeg.cluster0, 'sessCntSTD') extra.ss = eeg.cluster0.sessCntSTD; end
       if isfield(eeg.cluster0, 'sessCntMean') extra.mm = eeg.cluster0.sessCntMean; end
       if isfield(eeg.cluster0, 'sessCntMax') extra.ma = eeg.cluster0.sessCntMax; end
    end
end
if ok
   [eeg, eegdata] = DoCluster0Analysis(eeg, eegdata, eegind, extra, vv);
end

function [eeg, eegdata] = DoCluster0Analysis(eeg, eegdata, eegind, extra, vv)
%%%%%%%%%% Do simple analysis on cluster0 files: detect frame start and end times, compute basic frame properties
neeg = numel(eeg.general.eegfile);
%%%the following variables are assigned
% if (~isfield(eegdata, 'cluster0')) eegdata.cluster0 = []; end
% if (~isfield(eegdata.cluster0, 'timepoint')) eegdata.cluster0.tiempoint = cell(1, neeg); end %%all computed results are here
%%%%displayed variables
if (~isfield(eeg, 'cluster0')) eeg.cluster0 = []; end
if (~isfield(eeg.cluster0, 'startThreshold')) eeg.cluster0.startThreshold = cell(1, neeg); end
if (~isfield(eeg.cluster0, 'peakThreshold')) eeg.cluster0.peakThreshold = cell(1, neeg); end
if (~isfield(eeg.cluster0, 'sessCntSTD')) eeg.cluster0.sessCntSTD = cell(1, neeg); end
if (~isfield(eeg.cluster0, 'sessCntMean')) eeg.cluster0.sessCntMean = cell(1, neeg); end
if (~isfield(eeg.cluster0, 'sessCntMax')) eeg.cluster0.sessCntMax = cell(1, neeg); end

if (~isfield(eeg.cluster0, 'sessNum')) eeg.cluster0.sessNum = cell(1, neeg); end
if (~isfield(eeg.cluster0, 'sessRate')) eeg.cluster0.sessRate = cell(1, neeg); end
if (~isfield(eeg.cluster0, 'sessMeanAmp')) eeg.cluster0.sessMeanAmp = cell(1, neeg); end
if (~isfield(eeg.cluster0, 'sessMeanDur')) eeg.cluster0.sessMeanDur = cell(1, neeg); end
if (~isfield(eeg.cluster0, 'sessMeanEng')) eeg.cluster0.sessMeanEng = cell(1, neeg); end
if (~isfield(eeg.cluster0, 'sessMeanSDur')) eeg.cluster0.sessMeanSDur = cell(1, neeg); end

if (~isfield(eeg.cluster0, 'evtNum')) eeg.cluster0.evtNum = cell(1, neeg); end
if (~isfield(eeg.cluster0, 'evtRate')) eeg.cluster0.evtRate = cell(1, neeg); end
if (~isfield(eeg.cluster0, 'evtMeanAmp')) eeg.cluster0.evtMeanAmp = cell(1, neeg); end
if (~isfield(eeg.cluster0, 'evtMeanDur')) eeg.cluster0.evtMeanDur = cell(1, neeg); end
if (~isfield(eeg.cluster0, 'evtMeanEng')) eeg.cluster0.evtMeanEng = cell(1, neeg); end
if (~isfield(eeg.cluster0, 'evtMeanSDur')) eeg.cluster0.evtMeanSDur = cell(1, neeg); end

if (~isfield(eeg.cluster0, 'sessStartT')) eeg.cluster0.sessStartT = cell(1, neeg); end
if (~isfield(eeg.cluster0, 'sessEndT')) eeg.cluster0.sessEndT = cell(1, neeg); end
if (~isfield(eeg.cluster0, 'sessPeakT')) eeg.cluster0.sessPeakT = cell(1, neeg); end
if (~isfield(eeg.cluster0, 'sessDur')) eeg.cluster0.sessDur = cell(1, neeg); end
if (~isfield(eeg.cluster0, 'sessSDur')) eeg.cluster0.sessSDur = cell(1, neeg); end
if (~isfield(eeg.cluster0, 'sessAmp')) eeg.cluster0.sessAmp = cell(1, neeg); end
if (~isfield(eeg.cluster0, 'sessEng')) eeg.cluster0.sessEng = cell(1, neeg); end

for (iiik = 1:numel(eegind))
i = eegind(iiik);
%if strcmp(eeg.parm.band{i}, 'cluster0') %%only do this for filtered EEG traces
    disp(['--------> cluster0 analysis: ', eeg.general.eegfile{i}]);
    sth = eeg.parm.cl0StartThres(i); pth = eeg.parm.cl0PeakThres(i); 
    maxgap = eeg.parm.cl0MaxGap(i); binsize = eeg.parm.cl0Binsize(i);
    maxgap = round(maxgap/binsize); %%%maxgap now in number of data points
    mindur = eeg.parm.cl0MinDur(i); maxdur = eeg.parm.cl0MaxDur(i); 
    %%%here are data to analyze
    if strcmp(eeg.parm.band{i}, 'cluster0')
        timepoint = eegdata.cluster0.timepoint{i}; cnt = eegdata.cluster0.cnt{i};
    else
        [timepoint, cnt, ~, ~] = ReadEEGFile(eeg.general.eegfile{i}, eeg.parm.timeunit(i), eeg.parm.buffersize(i));% timestamp/point in second, dat/cnt in mV
    end
    %%%thresholds ans spike cnt parameters
    eeg.cluster0.startThreshold{i} = sth; eeg.cluster0.peakThreshold{i} = pth; 
    eeg.cluster0.sessCntSTD{i} = extra.ss{i}; eeg.cluster0.sessCntMean{i} = extra.mm{i}; eeg.cluster0.sessCntMax{i} = extra.ma{i};
    %%%detect cluster0 frames (UP states)
    cl0 = findallcluster0s(timepoint, cnt, sth, pth, maxgap, mindur, maxdur); %%cl0.startT(i), .peakT, .endT, .dur, .amp, .eng, .sdur; 
    %%%session variables
    if (~isempty(eeg.general.sesslength{i})) && (eeg.general.sesslength{i}>0)
        eeg.cluster0.sessNum{i} = numel(cl0.startT); eeg.cluster0.sessRate{i} = numel(cl0.startT)/eeg.general.sesslength{i};
        eeg.cluster0.sessMeanAmp{i} = mean(cl0.amp); eeg.cluster0.sessMeanDur{i} = mean(cl0.dur); 
        eeg.cluster0.sessMeanEng{i} = mean(cl0.eng);
        eeg.cluster0.sessStartT{i} = cl0.startT; eeg.cluster0.sessPeakT{i} = cl0.peakT; eeg.cluster0.sessEndT{i} = cl0.endT; 
        eeg.cluster0.sessAmp{i} = cl0.amp; eeg.cluster0.sessDur{i} = cl0.dur; eeg.cluster0.sessSDur{i} = cl0.sdur; eeg.cluster0.sessEng{i} = cl0.eng;
        eeg.cluster0.sessMeanSDur{i} = mean(cl0.sdur);
        %%%event variables
        evTime = eegdata.event.eventtimes{i};
        for (j = 1:numel(evTime))
            eeg.cluster0.evtNum{i}(j) = 0; eeg.cluster0.evtRate{i}(j) = 0; eeg.cluster0.evtMeanAmp{i}(j) = NaN; 
            %eeg.cluster0.evtMeanFreq{i}(j) = NaN; 
            eeg.cluster0.evtMeanDur{i}(j) = NaN; eeg.cluster0.evtMeanSDur{i}(j) = NaN;
            eeg.cluster0.evtMeanEng{i}(j) = NaN;
            st = evTime{j}.start; et = evTime{j}.ent; nev = numel(st); evtid = []; leng = 0;
            for (ttt = 1:nev)
                ik = find( (cl0.peakT>=st(ttt)) & (cl0.peakT<et(ttt)) ); evtid = union(evtid, ik); ik = []; leng = leng + et(ttt)-st(ttt);
            end
            if (leng > 0)
                eeg.cluster0.evtNum{i}(j) = numel(cl0.startT(evtid)); eeg.cluster0.evtRate{i}(j) = numel(cl0.startT(evtid))/leng;
                eeg.cluster0.evtMeanAmp{i}(j) = mean(cl0.amp(evtid)); eeg.cluster0.evtMeanDur{i}(j) = mean(cl0.dur(evtid));
                eeg.cluster0.evtMeanEng{i}(j) = mean(cl0.eng(evtid));
                eeg.cluster0.evtMeanSDur{i}(j) = mean(cl0.sdur(evtid));
            end
        end  
    end
%end
end

function [timepoint, cnt, newmode, ss, mm, ma, mi] = computecl0pwr(sptime, starttime, endtime, parm, eeg, sessid)
[cnt, timepoint] = gethistcnt(sptime, starttime, endtime, parm); %%%counted and smoothed if specified
newmode = parm.normmode;
%%%%%%%%%%now need to locate the normalization factors: run, self, or none
if (strcmp(parm.normmode, 'self'))
   ss = std(cnt); mm = mean(cnt); ma = max(cnt); mi = min(cnt);
   disp(['---------------> self std: ', num2str(ss), '; mean:', num2str(mm), '; max:', num2str(ma), '; min:', num2str(mi)]);
elseif (strcmp(parm.normmode, 'run'))
   runcnt = findrunsessstd(eeg, eeg.general.animalname{sessid}, eeg.general.finaldir{sessid}, eeg.general.recarea{sessid}, parm); 
   if (~isempty(runcnt))
       ss = std(runcnt); mm = mean(runcnt); ma = max(runcnt); mi = min(cnt);
       disp(['---------------> run std: ', num2str(ss), '; mean:', num2str(mm), '; max:', num2str(ma), '; min:', num2str(mi)]);
   else
       disp('---------------> Warning: run cluster0 data not found; use self normalization mode instead');
       ss = std(cnt); mm = mean(cnt); ma = max(cnt); mi = min(cnt);
       disp(['---------------> self std: ', num2str(ss), '; mean:', num2str(mm), '; max:', num2str(ma), '; min:', num2str(mi)]);
       newmode = 'self';
   end
else
   mm = mean(cnt); ma = max(cnt); mi = min(cnt); ss = std(cnt);
end
%%%%%%%%%%%%% normalize: smooth first to avoid singular normalization factors
if (strncmpi(parm.normmeth, 'Z-transform', 3))
    ssss = std(cnt); cnt = cnt - mm; if (ssss~=0) cnt = cnt/ssss; end
elseif (strncmpi(parm.normmeth, '0-1', 3)) 
    if (ma-mi > 0) cnt = (cnt-mi)/(ma-mi);
    else
        cnt = zeros(size(cnt));
    end
elseif (strncmpi(parm.normmeth, 'mean-subtraction', 3)) 
    cnt = cnt-mm;
end

function [eeg, eegdata, extra, eegind] = appendcl0file(eeg, eegdata, index, eegfilenow, timepoint, cnt, newmode, extra, eegind, ss, mm, ma, mi)
if strcmp(eeg.parm.band{index}, 'cluster0') %%%if the analyzed file itself is already a cluster0 files, just replace the parameters
    extra.ss{index} = ss; extra.mm{index} = mm; extra.ma{index} = ma; extra.mi{index} = mi;
    eegdata.cluster0.timepoint{index} = timepoint; eegdata.cluster0.cnt{index} = cnt;
else
nfile = numel(eeg.general.eegfile); [~, ID, ~] = fileparts(eegfilenow);
newind = nfile + 1; eegind = union(eegind, newind);
%%%%%%%%% reassign general variables
eeg.general.eegfile{nfile+1} = eegfilenow; eeg.general.eegID{nfile+1} = ID;
eeg.general.eeggain{nfile+1} = 1; eeg.general.freq{nfile+1} = NaN;
eeg.general.eegTT{nfile+1} = 'multiple'; 
eeg.general.recarea{nfile+1} = eeg.general.recarea{index}(1:2);
%%%%%%%%% re-assign parm variables 
eeg.parm.band{nfile+1} = 'cluster0'; eeg.parm.buffersize(nfile+1) = NaN; eeg.parm.cl0NormMode{nfile+1} = newmode;
%%%%%%%%%keep old general and parm variables
genvar = fieldnames(eeg.general);
for (i = 1:numel(genvar))
    if (~strcmp(genvar{i}, 'eegfile')) && (~strcmp(genvar{i}, 'eegID')) && (~strcmp(genvar{i}, 'eeggain'))...
            && (~strcmp(genvar{i}, 'freq')) && (~strcmp(genvar{i}, 'eegTT')) && (~strcmp(genvar{i}, 'recarea'))
        eeg.general.(genvar{i}){nfile+1} = eeg.general.(genvar{i}){index};
    end
end
parmvar = fieldnames(eeg.parm);
for (i = 1:numel(parmvar))
    if (~strcmp(parmvar{i}, 'band')) && (~strcmp(parmvar{i}, 'cl0NormMode'))
        if (iscell(eeg.parm.(parmvar{i})))
           eeg.parm.(parmvar{i}){nfile+1} = eeg.parm.(parmvar{i}){index};
        elseif (isnumeric(eeg.parm.(parmvar{i})))
           eeg.parm.(parmvar{i})(nfile+1) = eeg.parm.(parmvar{i})(index);
        end
    end
end
%%%%pass the normalization parameters to extra
extra.ss{nfile+1} = ss; extra.mm{nfile+1} = mm; extra.ma{nfile+1} = ma; extra.mi{nfile+1} = mi;
%%%%%%%%%% add new data to eegdata
eegdata.cluster0.timepoint{nfile+1} = timepoint; eegdata.cluster0.cnt{nfile+1} = cnt;
eegdata.event.eventtimes{nfile+1} = eegdata.event.eventtimes{index}; 
grpname = eegdata.grouplist.groupname; zeroind = strcmp(grpname, 'List0'); 
eegdata.grouplist.groupindex{zeroind} = 1:numel(eeg.general.eegfile);
%%%%%Also need to add empty variables for other types of analysis: ripple/theta etc.
catnames = fieldnames(eeg);
for (i = 1:numel(catnames))
    if (~strcmp(catnames{i}, 'general')) && (~strcmp(catnames{i}, 'parm'))
        if (~isempty(eeg.(catnames{i})))
            subvar = fieldnames(eeg.(catnames{i}));
            for (j = 1:numel(subvar))
                if (~isempty(eeg.(catnames{i}).(subvar{j})))
                    %disp([catnames{i}, '.', subvar{j}]);
                    eeg.(catnames{i}).(subvar{j}){nfile+1} = []; 
                end
            end
        end
    end
end
end

function cl0 = findallcluster0s(timestamp, dat, sth, pth, maxgap, mindur, maxdur) %%cl0.startT(i), .peakT, .endT, .dur, .amp, .eng, .sdur; 
cl0.startT = []; cl0.peakT = []; cl0.endT = []; cl0.dur = []; cl0.amp = []; cl0.eng = []; cl0.sdur = [];
%index = find(abs(dat) >= sth); %all threshold crossing for both positive crossing and negative crossing
index = find(dat >= sth);
if (~isempty(index))
    [startindex, endindex] = findeventindex(index, maxgap); %%start/end indices for start threshold crossing
    cl0 = findamp(dat, timestamp, startindex, endindex, pth, mindur, maxdur); %%get peak values and filter through peak/duration
    %disp([startindex' endindex']);
end

function [startindex, endindex] = findeventindex(index, groupspan)
%%index = a vector of integer indices of EEG file, in ascending order 
%%groupspan = any points with gap smaller than group span belong to same event
%%[startindex endindex] = detected event start index and end index, real numbers in index
startindex = []; endindex = [];
if (numel(index) == 1)
   startindex = index(1); endindex = index(1); 
else
   index = sort(index); %% with indexing [1,NN]
   indexdiff = diff(index); %get difference with indexing [1, NN-1] = original [2, NN]
   spanindex = find(indexdiff > groupspan); %%get all gaps bigger than groupspan, spanindex are indices in indexdiff
   if (isempty(spanindex)) %%if all within groupspan = only one event
       startindex = index(1); endindex = index(numel(index));
   else
       nindex = numel(spanindex); startindex = ones(1,nindex+1); endindex = ones(1,nindex+1);
       startindex(1) = index(1);
       endindex(1) = index(spanindex(1));
       for (i = 2:nindex)
           startindex(i) = index(spanindex(i-1)+1);
           endindex(i) = index(spanindex(i));
       end
       startindex(nindex+1) = index(spanindex(nindex)+1);
       endindex(nindex+1) = index(numel(index));
   end
end

function cl0 = findamp(dat, timestamp, startindex, endindex, pth, mindur, maxdur)
%%individual event indicated by [startindex, endindex], indices across startthreshold
%%peaktime, amplitude = maximum values within [startindex endindex] (minimum for slow wave detection)
%%starttime = index in [startindex-groupspan startindex] but cross the startthreshold
%%endtime = index in [endindex endindex+groupspan]but cross the startthreshold
%%area = integrate through [tarttime endtime]
cl0.startT = NaN*ones(size(startindex)); cl0.peakT = NaN*ones(size(startindex)); cl0.endT = NaN*ones(size(startindex)); 
cl0.dur = NaN*ones(size(startindex)); cl0.amp = NaN*ones(size(startindex)); cl0.eng = NaN*ones(size(startindex)); 
cl0.sdur = NaN*ones(size(startindex)); ntime = numel(timestamp);
for (i = 1:numel(startindex))
     datnow = dat(startindex(i):endindex(i)); timenow = timestamp(startindex(i):endindex(i)); 
     [abspeak, peakindex] = max(abs(datnow)); %%if either positive crossing or negative crossing
     cl0.amp(i) = abspeak; cl0.peakT(i) = timenow(peakindex);
     cl0.startT(i) = timestamp(min([startindex(i) ntime])); 
     cl0.endT(i) = timestamp(endindex(i));
     cl0.dur(i) = cl0.endT(i) - cl0.startT(i);
     cl0.eng(i) = mean(abs(datnow))* cl0.dur(i);
end
for (i = 1:numel(startindex)-1)
     cl0.sdur(i) = cl0.startT(i+1) - cl0.endT(i);
end
cl0.sdur(numel(startindex)) = mean(cl0.sdur(1:numel(startindex)-1)); %%%the last value takes the mean of all the previous ones
%%%%filter through peak threshold and duration
ind = find((cl0.amp>=pth) & ((cl0.dur>=mindur)&(cl0.dur<=maxdur))); 
fcl0 = fieldnames(cl0);
for (i = 1:numel(fcl0))
    cl0.(fcl0{i}) = cl0.(fcl0{i})(ind);
end

function runcnt = findrunsessstd(eeg, animname, fdir, recarea, parm)
runcnt = []; 
nfile = 0; sptime = []; sT = []; eT = [];
%%%look for all the run cluster0.spm data
pp{1} = strcat(fdir, filesep, 'spikedata');
[fpath, fname] = GetAllFile(pp, 'cluster0_', '.spm');
for (i = 1:numel(fname))
    if ~isempty(strfind(fname{i}, animname)) && ~isempty(strfind(fname{i}, recarea))
       ki = strfind(fname{i}, '_'); [str, tok] = strtok(fname{i});
       sessthen = str([ki(numel(ki))+1:numel(str)-4]);
       iii = find( strcmp(eeg.general.finaldir, fdir) & strcmp(eeg.general.sessname, sessthen) );
       if (~isempty(iii))
          sesstype = eeg.parm.sessType{iii(1)};
          if strncmpi(sesstype, 'linear', 4) || strncmpi(sesstype, 'open', 4) %%%%if linear or open session
             nfile = nfile + 1;
             spnow = ReadSpikeTime(fullfile(fpath{i}, fname{i}));
             [runstart, runend] = findrunepisodes(fdir, sessthen);
             if (~isempty(runstart))
             for (k = 1:numel(runstart))
                 if (runend(k)-runstart(k)) > 100*parm.binsize  %%if the run episode is not so short
                     spok = spnow( (spnow>=runstart(k)) & (spnow<=runend(k)) );
                     [cnt, ~] = gethistcnt(spok, runstart(k), runend(k), parm);
                     runcnt = [runcnt; cnt];
                 end
             end
             end
             sptime{nfile} = spnow; sT(nfile) = eeg.general.sessstartT{iii(1)}; eT(nfile) = eeg.general.sessendT{iii(1)};
          end
       end
    end
end
%numel(runcnt)
if (numel(runcnt)<100) %%%%if less than 100 samples for running, use the entire run session
    for (i = 1:nfile)
        [cnt, ~] = gethistcnt(sptime{i}, sT(i), eT(i), parm);
        runcnt = [runcnt; cnt];
    end
    if (~isempty(runcnt))
       disp('---------------> Warning: motion power file not found or no running detected; use whole run session for computing thresholds');
    end
end

function [cnt, timepoint] = gethistcnt(sptime, starttime, endtime, parm)
nbin = floor( (endtime-starttime)/parm.binsize ); %number of bins for this session
timepoint = starttime + ((1:nbin+1)-1)*parm.binsize; %timebin in second = row vector
cnt = hist(sptime, timepoint); %output column or row vector, centerred time bins, same as spiketime = column vector
cnt = cnt(1:nbin); timepoint = timepoint(1:nbin); %column vector
if (strncmpi(parm.smooth, 'yes', 1))  %%if smooth
    XX = [-parm.Nsigma*parm.sigma:1:parm.Nsigma*parm.sigma]; 
    wind = normpdf(XX, 0, parm.sigma); wind = wind/sum(wind); nw = numel(wind); %will have a delay = (n-1)/2 bins
    cnt = conv(cnt, wind); cnt = cnt((nw-1)/2+1:numel(cnt)-(nw-1)/2); 
end

function [runstart, runend] = findrunepisodes(finaldir, sessname)
runstart = []; runend = []; powerfile = [];
posdir{1} = strcat(finaldir, filesep, 'pos');
[filepath, filename] = GetAllFile(posdir, '', '.pwr'); 
for  (j = 1:numel(filename)) %%%green diodes first: better signal in most of recordings
     if (~isempty(strfind(filename{j}, 'green'))) && (~isempty(strfind(filename{j}, sessname)))
        powerfile = fullfile(filepath{j}, filename{j}); break
     end
end
if isempty(powerfile)
   for (j = 1:numel(filename))
        if (~isempty(strfind(filename{j}, sessname)))
            powerfile = fullfile(filepath{j}, filename{j}); break
        end
   end
end
if (~isempty(powerfile))
    [data, ~, ~, ~, ~] = ReadPowergramData(powerfile);
    timenow = 0.0001*data(:,1); pwrnow = data(:,2); %%%% timestamps in power data are in 0.1ms
    state = ones(size(pwrnow)); %%%initial assignment: intermediate state
    iii = find(pwrnow >=40); state(iii) = 2*ones(size(iii)); %%%%running state: >50 fixed threshold in pixels per second
    %jjj = find(pwrnow <=30); state(jjj) = zeros(size(jjj)); %%%%stopp state: <30 fixed threshold in pixels per second
    [stateStart, stateEnd] = findstateepisode(state, 1, numel(state), 2);
    runstart = timenow(stateStart); runend = timenow(stateEnd);
    %iii = find(runend-runstart>2); runstart = runstart(iii); runend = runend(iii); %%%select events longer than 2s
end
%disp([runstart runend]);
    
function [stateStart, stateEnd] = findstateepisode(state, starttime, endtime, value)
%find continous sleep episode with the same sleep state value (0=wake, 1=tran, 2=sws, 3=rem)
%state = classified point by point sleep state, 
%starttime/endtime = start/end time point of sleep classification
stateStart = []; stateEnd = [];
istate= 0; %number of episode
i = starttime;
while (i < endtime)
    if (i == 1) 
        if (state(i) == value)&&(state(i+1)==value)
           istate = istate + 1;
           stateStart(istate) = i; %%assign the starting point
           j = i + 1;
           endflag = 0;  %%endflag = 0 continue to search for an ending point
           while ((j <= endtime) && (endflag == 0))
             if (j == endtime)   %%if this is the last point of the whole data
                stateEnd(istate) = j;
                i = j + 1;
                endflag = 1;
             else
                if ((state(j-1) == value) && (state(j) == value) && (state(j+1) ~= value))   %%detect state ending point
                    stateEnd(istate) = j; %%assign the ending point
                    endflag = 1;
                    i = j + 1;
                else
                    j = j +1;
                end
             end
           end
        else
           i = i + 1;
        end
    else
        if ((state(i) == value) && (state(i-1) ~=value) && (state(i+1) == value))   %%detect state starting point
           istate = istate + 1;
           stateStart(istate) = i; %%assign the starting point
           j = i + 1;
           endflag = 0;  %%endflag = 0 continue to search for an ending point
           while ((j <= endtime) && (endflag == 0))
             if (j == endtime)   %%if this is the last point of the whole data
                stateEnd(istate) = j;
                i = j + 1;
                endflag = 1;
             else
                if ((state(j-1) == value) && (state(j) == value) && (state(j+1) ~= value))   %%detect state ending point
                    stateEnd(istate) = j; %%assign the ending point
                    endflag = 1;
                    i = j + 1;
                else
                    j = j +1;
                end
             end
           end
        else
           i = i + 1;
        end
    end
end