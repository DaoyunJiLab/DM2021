function DataManager_EEG_speedThetaPeakPowerRun_Callback
%%% Compute theta power/frequency in segmented time periods during running events
%%% Also estimated the 2D speed within the same time periods
%%% Use Hilbert-Huang transform to analyze filtered theta band LFPs
%%%       Or use pmtm to analyze broad band LFPs

Q = questdlg('Read (or compute) from a file?');
if strcmp(Q, 'Yes') %%%read from a file
    [p,n] = uigetfile(cd);
    S = load(fullfile(p,n)); OutD = S.OutD; S = [];
else
    
hf = gcbf; pinfo = getappdata(hf, 'eeg'); data = getappdata(hf, 'eegdata'); plotparm = getappdata(hf, 'plotparm');
hgroup = getappdata(hf, 'hgroup'); 
[writefilename, okk] = getoutputfile;
if okk
    if (plotparm.linkbehav == 0)
       disp(['--------> no behav data linked']); okk = 0;
    else
       behav = getappdata(hf, 'behav'); bhdata = getappdata(hf, 'bhdata');
    end
end
if okk
   groupselection = getappdata(hgroup, 'selection'); grpind = find(groupselection == 1); ngrp = numel(grpind); eegind = cell(ngrp);
   for (kk = 1:ngrp) eegind{kk} = data.grouplist.groupindex{grpind(kk)}; end
   OutD.groupname = data.grouplist.groupname{grpind};
   OutD.bandPeakF = cell(ngrp); OutD.bandTotalPower = cell(ngrp);
   OutD.bandType = cell(ngrp); OutD.speed2D = cell(ngrp); OutD.speed1D = cell(ngrp); %%%%2D speed is not computed yet: empty for now
   %%%%%%%%%%%%below are the parameters for band specifications
   parm.bandF = [6 10]; %%%%frequency range in Hz to focus
   for i = 1:ngrp parm.targetLFPtype{i} = 'broad'; end  %%% default broad LFPs
   parm.computeMethod = 'multi-taper'; %%% default method is multitaper
   parm.timeparse.win = 0.5; %%%default 0.5 s time window
   parm.normrange = [4 200]; parm.normalize = 'no'; %%normalization frequency range 
   %%%%%%% event selection critetia:
   parm.event.evkeyword = 'Track'; parm.event.evkeytype = 'run'; 
   parm.event.evkeynoword = 'active'; parm.event.evkeynotype = 'others'; 

   %%%%%now re-detect LFP type
   for i = 1:ngrp
       if ~isempty(eegind{i})
           parm.targetLFPtype{i} = pinfo.parm.band{eegind{i}(1)}; OutD.bandType{i} =  pinfo.parm.band(eegind{i});
           OutD.speed2D{i} = cell(1, numel(eegind{i})); OutD.speed1D{i} = cell(1, numel(eegind{i})); 
           OutD.bandPeakF{i} = cell(1, numel(eegind{i})); OutD.bandTotalPower{i} = cell(1, numel(eegind{i})); 
           if ~isempty(find(~strcmp(pinfo.parm.band(eegind{i}), parm.targetLFPtype{i}))) parm.targetLFPtype{i} = 'broad'; end
       end
   end

   if isempty(find(~strcmp(parm.targetLFPtype, 'broad'))) 
       parm.computeMethod = 'Hiltbert-Huang'; %%%Two choices here: 'Hilbert-Huang' or 'multi-taper' 
   end
   OutD.parm = parm;

   %%%%perform the spectral analysis
   for i = 1:ngrp
   for j = 1:numel(eegind{i})
        eegidnow = eegind{i}(j);
        %%%%% get eegdata
        [T, dat, gain, fs] = ReadEEGFile(eeg.general.eegfile{i}, eeg.parm.timeunit(i), eeg.parm.buffersize(i)); % timestamp in second, dat in mV 
        dat = dat-mean(dat);
        %%%%% get event times
        [lapevT, lapposT, lapS1D] = geteventdata(pinfo, eegidnow, behav, bhdata, parm.event); 
        if ~isempty(evT.start)
           [F,P,S1D] = getspectralspeednumbers(T, dat, fs, lapevT, lapposT, lapS1D, parm);
        end
        OutD.bandPeakF{i}{j} = F; OutD.bandTotalPower{i}{j} = P; OutD.speed1D{i}{j} = S1D;
   end
   end
   save(writefilename, 'OutD');

else
   disp(['--------> aborted']);
end
plotOutD(OutD);
end

function [F,P,S1D] = getspectralspeednumbers(T, dat, fs, lapevT, lapposT, lapS1D, parm)
%%%% compute for single session, still multiple events - lapevT{i} with multiple laps - lapposT{i}{j}



function [evTout, lapTout, lapS1Dout] = geteventdata(pinfo, eegidnow, behav, bhdata, parm)
evTout = []; lapTout = []; lapS1Dout = []; %%%evTout{i}.start(j)/ent' lapT{i}{j} = i - events; j = lap number
%%%%identify behav event ID: posid
finaldirnow = pinfo.general.finaldir{eegidnow}; sessnamenow = pinfo.general.sessname{eegidnow}; 
posid = find( strcmp(behav.general.finaldir, finaldirnow) & strcmp(behav.general.sessname, sessnamenow) );
%%%%select events
if ~isempty(posid)
   evName = behav.general.eventname{posid}; evType = behav.parm.eventtype{posid};
   evT = bhdata.event.eventtimes{posid}; 
   lapT = bhdata.event.LapAllPostimestamp{posid}; lapS1D = bhdata.event.LapAll1DSpeed{posid};
   evsel = checkifcorrectevents(evName, evType, parm.evkeyword, parm.evkeytype, parm.evkeynoword, parm.evkeynotype);
   evpos = find(evsel == 1); 
   if ~isempty(evpos)
      evTout = evT(evpos); lapTout = lapT(evpos); lapS1Dout = lapS1D(evpos);
   else
      disp(['------------> warning: no events found in ' finaldirnow ': ' sessnamenow]);
   end
else
   disp(['------------> warning: no matching behavior data found for ' finaldirnow ': ' sessnamenow]);
end
function evsel = checkifcorrectevents(evName, evType, evkeyword, evkeytype, evkeynoword, evkeynotype)
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

function [writefilename, okk] = getoutputfile
okk = 1; writefilename = [];
[fname, pname] = uiputfile(fullfile(cd, '*.eegdb'), 'Write the new spike database to:');
if (numel(fname)>1)
   writefilename = fullfile(pname, fname);
else
   okk = 0;
end
