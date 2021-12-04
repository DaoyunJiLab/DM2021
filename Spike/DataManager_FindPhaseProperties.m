function [pinfo,data] = DataManager_FindPhaseProperties(pinfo,data,eeg,eegdata,behav,bhdata,cellind)
%%%Compute phase properties of individual cells
%%%mean phase/tuning significance for each session and each event
%%%phase precession within every place field (start phase, end phase,
%%%   phase_postion correlatin and significance

%variable to assign
if (~isempty(cellind))
   nspike = numel(pinfo.general.parmfile);
   if (~isfield(pinfo, 'phase')) pinfo.phase = []; end
   %%%session field dynam
   if (~isfield(pinfo.phase, 'sessPeakPhase')) pinfo.phase.sessPeakPhase = cell(1, nspike); end %%%Theta peak phase
   if (~isfield(pinfo.phase, 'sessMeanPhase')) pinfo.phase.sessMeanPhase = cell(1, nspike); end
   if (~isfield(pinfo.phase, 'sessPhaseVar')) pinfo.phase.sessPhaseVar = cell(1, nspike); end
   if (~isfield(pinfo.phase, 'sessRayLength')) pinfo.phase.sessRayLength = cell(1, nspike); end  %%%Rayleigh length of all spikes' phases
   if (~isfield(pinfo.phase, 'sessRayleighP')) pinfo.phase.sessRayleighP = cell(1, nspike); end  %%%Rayleigh P-value of all spikes' phases
   %%%evt field dynam
   if (~isfield(pinfo.phase, 'evtPeakPhase')) pinfo.phase.evtPeakPhase = cell(1, nspike); end
   if (~isfield(pinfo.phase, 'evtMeanPhase')) pinfo.phase.evtMeanPhase = cell(1, nspike); end
   if (~isfield(pinfo.phase, 'evtPhaseVar')) pinfo.phase.evtPhaseVar = cell(1, nspike); end
   if (~isfield(pinfo.phase, 'evtRayLength')) pinfo.phase.evtRayLength = cell(1, nspike); end  %%%Rayleigh length of all spikes' phases
   if (~isfield(pinfo.phase, 'evtRayleighP')) pinfo.phase.evtRayleighP = cell(1, nspike); end  %%%Rayleigh P-value of all spikes' phases
   %%%%%phase precession within place/time fields of 1D track
   if (~isfield(pinfo.phase, 'fieldPeakPhase')) pinfo.phase.fieldPeakPhase = cell(1, nspike); end
   if (~isfield(pinfo.phase, 'fieldMeanPhase')) pinfo.phase.fieldMeanPhase = cell(1, nspike); end
   if (~isfield(pinfo.phase, 'fieldPhaseVar')) pinfo.phase.fieldPhaseVar = cell(1, nspike); end
   if (~isfield(pinfo.phase, 'fieldRayLength')) pinfo.phase.fieldRayLength = cell(1, nspike); end
   if (~isfield(pinfo.phase, 'fieldRayleighP')) pinfo.phase.fieldRayleighP = cell(1, nspike); end
   if (~isfield(pinfo.phase, 'fieldStartPhase')) pinfo.phase.fieldStartPhase = cell(1, nspike); end
   if (~isfield(pinfo.phase, 'fieldEndPhase')) pinfo.phase.fieldEndPhase = cell(1, nspike); end
   if (~isfield(pinfo.phase, 'fieldPPcircCrr')) pinfo.phase.fieldPPcircCrr = cell(1, nspike); end
   if (~isfield(pinfo.phase, 'fieldPPCorr')) pinfo.phase.fieldPPCorr = cell(1, nspike); end
   if (~isfield(pinfo.phase, 'fieldPPCorrP')) pinfo.phase.fieldPPCorrP = cell(1, nspike); end
   if (~isfield(pinfo.phase, 'fieldPPR2')) pinfo.phase.fieldPPR2 = cell(1, nspike); end
   if (~isfield(pinfo.phase, 'fieldPPSlope')) pinfo.phase.fieldPPSlope = cell(1, nspike); end
   if (~isfield(pinfo.phase, 'fieldPPInt')) pinfo.phase.fieldPPInt = cell(1, nspike); end
   if (~isfield(pinfo.phase, 'fieldPPOffset')) pinfo.phase.fieldPPOffset = cell(1, nspike); end    
   
   %%%data variables to assign
   if (~isfield(data, 'phase')) data.phase = []; end
   if (~isfield(data.phase, 'SpikePhases')) data.phase.SpikePhases = cell(1, nspike); end
end
currentfinaldir = []; currentsessname = []; currentEEGphase = []; currenttimepont = [];
for (jjjk = 1:numel(cellind))
    i = cellind(jjjk); finaldir = pinfo.general.finaldir{i};
    %%%%load computing parameters 
    parm.recarea = pinfo.parm.phaseEEGarea{i}; parm.band = pinfo.parm.phaseEEGband{i}; parm.EEGkeyword = pinfo.parm.phaseEEGkeyword{i};
    parm.sesskeyword = pinfo.parm.phaseSessKeyword{i}; parm.sesskeytype = pinfo.parm.phaseSessKeytype{i};
    parm.evtkeyword = pinfo.parm.phaseEvtKeyword{i}; parm.evtkeytype = pinfo.parm.phaseEvtKeytype{i};
    %parm.allslope = pinfo.parm.phaseMinSlope(i):pinfo.parm.phaseSlopeStep(i):pinfo.parm.phaseMaxSlope(i);
    %parm.alloffset = pinfo.parm.phaseMinOffset(i):pinfo.parm.phaseOffsetStep(i):pinfo.parm.phaseMaxOffset(i);
    disp(['---> compute phase properties and phase precession variables (', num2str(jjjk), ' out of ', num2str(numel(cellind)), '): ', pinfo.general.parmfile{i}]);
    spiketime = sort(pinfo.parm.timeunit(i)*data.spike.spiketime{i}); [a,b] =size(spiketime); if a==1 spiketime = spiketime'; end %%%all column vectors
    %%%%locate EEG phases
    if ~strcmp(finaldir, currentfinaldir)
        currentfinaldir = finaldir;
        [currentsessname, currentEEGphase, currenttimepoint] = locateEEGphases(eeg, eegdata, finaldir, parm);   %%% timepoint already all column vectors
    end
    %%%%identify phases for all spikes of the current cell within each session
    spikephase = NaN*ones(size(spiketime));
    sessionname = pinfo.general.sessionname{i}; sessstartT = pinfo.general.sessionstartT{i}; sessendT = pinfo.general.sessionendT{i};
    nsess = numel(sessionname); sesstype = pinfo.parm.sessType{i}; pinfo.phase.sessPeakPhase{i} = NaN*ones(1, nsess);
    pinfo.phase.sessMeanPhase{i} = NaN*ones(1, nsess); pinfo.phase.sessPhaseVar{i} = NaN*ones(1, nsess);
    pinfo.phase.sessRayLength{i} = NaN*ones(1, nsess); pinfo.phase.sessRayleighP{i} = NaN*ones(1, nsess);
    for (tt = 1:nsess)
        iii = find(strcmp(currentsessname, sessionname{tt}));
        if (numel(iii)~=1) || isempty(currenttimepoint{iii})
            disp(['-------> warning: no EEG found in the session: ', sessionname{tt}]);
        else
            jjj = find( (spiketime>=sessstartT(tt)) & (spiketime<=sessendT(tt)) );
            ind = knnsearch(currenttimepoint{iii}, spiketime(jjj)); spikephase(jjj) = currentEEGphase{iii}(ind);
            if filterevents(sessionname{tt}, sesstype{tt}, parm.sesskeyword, parm.sesskeytype)
               if ~isempty(jjj)
                  P = findphaseprop(spikephase(jjj)); pinfo.phase.sessPeakPhase{i}(tt) = P.peakphase; 
                  pinfo.phase.sessMeanPhase{i}(tt) = P.meanphase; pinfo.phase.sessPhaseVar{i}(tt) = P.phasevar;
                  pinfo.phase.sessRayLength{i}(tt) = P.raylength; pinfo.phase.sessRayleighP{i}(tt) = P.rayp;
               end
            end
        end
    end
    [a,b] =size(spikephase); if a==1 spikephase = spikephase'; end %%%all column vectors
    data.phase.SpikePhases{i} = spikephase;
    %%%%%compute event phase properties
    evName = pinfo.general.eventname{i}; evTime = data.events.eventtimes{i}; evType = pinfo.parm.eventtype{i};       
    nev = numel(evName); pinfo.phase.evtPeakPhase{i} = NaN*ones(1, nev);
    pinfo.phase.evtMeanPhase{i} = NaN*ones(1, nev); pinfo.phase.evtPhaseVar{i} = NaN*ones(1, nev);
    pinfo.phase.evtRayLength{i} = NaN*ones(1, nev); pinfo.phase.evtRayleighP{i} = NaN*ones(1, nev);
    data.phase.evtSpikeTimes{i} = cell(1, numel(evName)); data.phase.evtSpikePhases{i} = cell(1, numel(evName));
    if ~isempty(find(~isnan(spikephase)))
    for (tt = 1:numel(evName))
        if filterevents(evName{tt}, evType{tt}, parm.evtkeyword, parm.evtkeytype)
            [~, ~, spikeid] = SpikeEventFilter(spiketime, evTime{tt});
            phasenow = spikephase(spikeid); P = findphaseprop(phasenow); pinfo.phase.evtPeakPhase{i}(tt) = P.peakphase;
            pinfo.phase.evtMeanPhase{i}(tt) = P.meanphase; pinfo.phase.evtPhaseVar{i}(tt) = P.phasevar;
            pinfo.phase.evtRayLength{i}(tt) = P.raylength; pinfo.phase.evtRayleighP{i}(tt) = P.rayp;
        end
    end
    end
    %%%%%compute phase prcession within 1D place/time fields
    fprate = pinfo.field.PF1DInPeakrate{i}; pinfo.phase.fieldPeakPhase{i} = NaN*ones(size(fprate));
    pinfo.phase.fieldMeanPhase{i} = NaN*ones(size(fprate)); pinfo.phase.fieldPhaseVar{i} = NaN*ones(size(fprate));
    pinfo.phase.fieldRayLength{i} = NaN*ones(size(fprate)); pinfo.phase.fieldRayleighP{i} = NaN*ones(size(fprate));
    pinfo.phase.fieldStartPhase{i} = NaN*ones(size(fprate)); pinfo.phase.fieldEndPhase{i} = NaN*ones(size(fprate));
    pinfo.phase.fieldPPCorr{i} = NaN*ones(size(fprate)); pinfo.phase.fieldPPCorrP{i} = NaN*ones(size(fprate));
    pinfo.phase.fieldPPSlope{i} = NaN*ones(size(fprate)); pinfo.phase.fieldPPOffset{i} = NaN*ones(size(fprate));
    pinfo.phase.fieldPPInt{i} = NaN*ones(size(fprate)); pinfo.phase.fieldPPR2{i} = NaN*ones(size(fprate));   
    pinfo.phase.fieldPPcircCrr{i} = NaN*ones(size(fprate));
    if ~isempty(find(~isnan(spikephase)))
    for tt = 1:numel(fprate) %%%for every field exist on linear track
        evNamenow = pinfo.field.PF1Devt{i}{tt}; [~,evTrajnow] = strtok(evNamenow, '_');
        j = find(strcmp(evName, evNamenow)); 
        if (numel(j)==1)&& filterevents(evName{j}, evType{j}, parm.evtkeyword, parm.evtkeytype) %%%if this is a run event on the selected traj 
            %%%%locate event position data
            evSess = identifysession(evTime{j}, sessionname, sessstartT, sessendT);
            posid = []; evid = [];
            if (~isempty(evSess))
                posid = find( strcmp(behav.general.finaldir, finaldir) & strcmp(behav.general.sessname, evSess) );
            end
            if numel(posid == 1)
                if (isfield(behav.general, 'eventname'))
                    evid = find(strcmp(behav.general.eventname{posid}, evName{j}));
                else
                    evid = find(strcmp(behav.behavior.eventname{posid}, evName{j}));
                end
            end
            if (numel(posid)~=1)||(numel(evid)~=1)
               disp(['------> phase precession not computed for this event: no or more than 1 positon/event files match the session: ', finaldir, '___', evName{j}]);
            else
               %fstart = pinfo.field.PF1DBoundStart{i}(tt); fend = pinfo.field.PF1DBoundEnd{i}(tt); fstart = fstart+1; fend = fend-1;
               fstart = pinfo.field.PF1DLocStartX{i}(tt); fend = pinfo.field.PF1DLocEndX{i}(tt);
               nlapnow = numel(evTime{j}.start);
               spikeid = [];  pos = []; %startid = []; endid = [];
               for (tj = 1:nlapnow)
                   lappostime = bhdata.event.LapAllPostimestamp{posid}{evid}{tj}; lapx = bhdata.event.LapAllX{posid}{evid}{tj};
                   lapstart = evTime{j}.start(tj); lapent = evTime{j}.ent(tj);
                   spikeidnow = find((spiketime>=lapstart) & (spiketime<=lapent)); spikenow = spiketime(spikeidnow); 
                   ind = knnsearch(lappostime, spikenow); posnow = lapx(ind);
                   jjj = find( (posnow>=fstart) & (posnow<=fend) );
                   if ~isempty(jjj)
                      spikeid = [spikeid; spikeidnow(jjj)]; pos = [pos; posnow(jjj)'];
                      %startid = [startid spikeidnow(jjj(1))]; endid =[endid spikeidnow(jjj(numel(jjj)))];
                   end
               end
               %%%%%data.phase.evtSpikeTimes{i}{tt} = spiketimeout;  data.phase.evtSpikePhases{i}{tt} = phasenow; 
               phase = spikephase(spikeid);
               P = findphaseprop(phasenow); %Q = findphaseprop(spikephase(startid)); R = findphaseprop(spikephase(endid));
               pinfo.phase.fieldMeanPhase{i}(tt) = P.meanphase; pinfo.phase.fieldPhaseVar{i}(tt) = P.phasevar; pinfo.phase.fieldPeakPhase{i}(tt) = P.peakphase;
               pinfo.phase.fieldRayLength{i}(tt) = P.raylength; pinfo.phase.fieldRayleighP{i}(tt) = P.rayp;
               %pinfo.phase.fieldStartPhase{i}(tt) = Q.meanphase; pinfo.phase.fieldEndPhase{i}(tt) = R.meanphase;
               %%%%computing phase precession below
               [circR, offset, opslope, opint, oprr, oppp, opr2, sphase, ephase] = findphaseprec(phase, pos);
               pinfo.phase.fieldPPCorr{i}(tt) = oprr; pinfo.phase.fieldPPCorrP{i}(tt) = oppp;
               pinfo.phase.fieldPPSlope{i}(tt) = opslope; pinfo.phase.fieldPPOffset{i}(tt) = offset;
               pinfo.phase.fieldPPInt{i}(tt) = opint; pinfo.phase.fieldPPR2{i}(tt) = opr2; 
               pinfo.phase.fieldPPcircCrr{i}(tt) = circR; pinfo.phase.fieldStartPhase{i}(tt) = sphase; pinfo.phase.fieldEndPhase{i}(tt) = ephase;             
            end
        end
    end
    end
end
disp('*****************');
%%%%%%%%%%%%%%%END OF THE MAIN ROUTINE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [circR, offset, opslope, opint, oprr, oppp, opr2, sphase, ephase] = findphaseprec(phasenow, posnow)
%allslope = -20:0.1:10; nslope = numel(allslope); 
%alloffset = -1000:10:1000; noffset = numel(alloffset); 
offset = NaN; opslope = NaN; opint = NaN; oprr = NaN; oppp = NaN; opr2 =NaN; circR = NaN; sphase = NaN; ephase = NaN;
nspike = numel(phasenow);
if (nspike > 2)
   circR = findcirccorrelation(phasenow', posnow');
   XX = [ones(nspike,1) posnow]; %nspike x 2 for two linear parameter regression
   r2 = NaN*ones(1, 359); slope = NaN*ones(1, 359); 
   for (i =1:359)
        phaseok = shiftphase(phasenow, i);  %%%%This is the key transformation to make sure the offset is done correctly
        %d(i) = Utilities_FindCorrCoef_self(phaseok',posnow'); %%%correlation as a metrics %%%a, b are row vectors
        [B1, Bint1, R, Rint, stat] = regress(phaseok, XX, 0.05); %%%regression fitness (R2) as meric 
        r2(i) = stat(1); slope(i) = B1(2);
   end 
   iii = find(slope<0); 
   if ~isempty(iii)
      [~,ii] = max(r2(iii)); offset = iii(ii); phasenow = shiftphase(phasenow, offset);
      oprr = Utilities_FindCorrCoef_self(phasenow',posnow'); %%%a, b are row vectors
      [B1, Bint1, R, Rint, stat] = regress(phasenow, XX, 0.05);
      opint = B1(1); opslope = B1(2); oppp = stat(3); opr2 = stat(1); %R2 = 1-sum(residuls) = percentatge of variability accounted by the regression
      Xnow = [ones(2,1) [min(posnow); max(posnow)]];
      YY = Xnow * B1; sphase = YY(1); ephase = YY(2);%regression result: line in (X, Y) = (xx, L4)
   end
end
function phasenow = shiftphase(phasenow, offset)  %%%%This is the key transformation to make sure the offset is done correctly
ind = find(phasenow<offset); phasenow(ind) = phasenow(ind) + 360;
function circR = findcirccorrelation(phasenow, posnow) %%%circular correlation between phase and position, following M Metha Science paper
S = sin(phasenow); C = cos(phasenow);
rcs = Utilities_FindCorrCoef_self(S,C); %%%a, b are row vectors
rcx = Utilities_FindCorrCoef_self(posnow,C); rsx = Utilities_FindCorrCoef_self(posnow,S);
circR = sqrt( (rcx^2 + rsx^2 - 2*rcx*rsx*rcs) / (1-rcs^2) );



% function [offset, opslope, opint, oprr, oppp, opr2] = findphaseprec(phasenow, posnow, allslope, alloffset)
% offset = NaN; opslope = NaN; opint = NaN; oprr = NaN; oppp = NaN; opr2 =NaN;
% nspike = numel(phasenow);
% if (nspike > 2)
%    XX = [ones(nspike,1) posnow]; %nspike x 2 for two linear parameter regression
%    r = zeros(1, 359); p = ones(1, 359);
%    for (i =1:359)
%        phaseok = mod(phasenow+i, 360);
%        [R, P] = corrcoef(phaseok, posnow); r(i) = R(1,2); p(i) = P(1,2);
%        %r(i) = Utilities_FindCorrCoef_self(phaseok',posnow'); %%%a, b are row vectors
%    end
%    %[~, ii] = max(abs(r)); 
%    [~,ii] = min(p);
%    offset = ii; phasenow = mod(phasenow+offset, 360);
%    oprr = Utilities_FindCorrCoef_self(phasenow',posnow'); %%%a, b are row vectors
%    [B1, Bint1, R, Rint, stat] = regress(phasenow, XX, 0.05);
%    opint = B1(1); opslope = B1(2); oppp = stat(3); opr2 = stat(1); %R2 = 1-sum(residuls) = percentatge of variability accounted by the regression
% end

% function [offset, opslope, opint, oprr, oppp, opr2] = findphaseprec(phasenow, posnow, allslope, alloffset)
% %allslope = -20:0.1:10; nslope = numel(allslope); 
% %alloffset = -1000:10:1000; noffset = numel(alloffset); 
% offset = NaN; opslope = NaN; opint = NaN; oprr = NaN; oppp = NaN; opr2 =NaN;
% nspike = numel(phasenow);
% if (nspike > 2)
%    XX = [ones(nspike,1) posnow]; %nspike x 2 for two linear parameter regression
% 
%    nslope = numel(allslope); noffset = numel(alloffset); 
%    errval = zeros(nslope*noffset,1);
%    %%pretreat the data: move the mean position to the center:
%    meanpos = mean(posnow); posok = posnow - meanpos;
%    for (i = 1:nslope)
%         for (j = 1:noffset)
%             ddd = mod(phasenow-(allslope(i)*posok+alloffset(j)), 360);
%             dd1 = ddd(find(ddd<=180)); dd2 = 360-ddd(find(ddd>180));
%             errval((i-1)*noffset+j) = errval((i-1)*noffset+j) + sum(dd1) + sum(dd2);  %sum of error: careful! has to be circular difference
%         end
%    end
%    [CCC,III] = min(errval); slopeind = ceil(III/noffset); offsetind = III-(slopeind-1)*noffset;
%    opslope = allslope(slopeind); offset = mod(alloffset(offsetind), 360) - opslope*meanpos;
%    %now rectify phase values to get correlation values and regression parameters
%    ddd = mod(phasenow-(opslope*posnow+offset), 360); ii1 = find(ddd<=180); ii2 = find(ddd>180);
%    phasenow(ii1) = opslope*posnow(ii1)+offset+ddd(ii1); phasenow(ii2) = opslope*posnow(ii2)+offset+ddd(ii2)-360;
%    
%    oprr = Utilities_FindCorrCoef_self(phasenow',posnow'); %%%a, b are row vectors
%    [B1, Bint1, R, Rint, stat] = regress(phasenow, XX, 0.05);
%    opint = B1(1); opslope = B1(2); oppp = stat(3); opr2 = stat(1); %R2 = 1-sum(residuls) = percentatge of variability accounted by the regression
% end


function P = findphaseprop(phasenow)
P.peakphase = NaN; P.meanphase = NaN; P.phasevar = NaN; P.raylength = NaN; P.rayp = NaN;
nspike = numel(phasenow); %phase in [0 360)
if (nspike > 1)
   P.peakphase = getpeakphase(phasenow);
   CC = exp(1i*phasenow*2*pi/360); %%now a complex number
   %first do the corcular statistics 
   finalV = sum(CC);
   Rayleng = abs(finalV); meanphase = angle(finalV); %%here phase is in [-pi pi), need to change to [0 360)
   if (meanphase >= 0)
       P.meanphase = 180 * meanphase / pi;
   else
       P.meanphase = 360 + (180 * meanphase / pi);
   end
   P.phasevar = 1 - Rayleng/nspike;
   %next do the Rayleigh test
   dist = zeros(1, nspike); %distance from origin
   for (k = 1:nspike)
       dist(k) = abs( sum(CC(1:k)) ); %samples of distance
   end
   b = raylfit(dist); %estimate of Rayleigh distribution parameter
   RayP = 1 - raylcdf(Rayleng, b); %final vector length's probability from Rayleigh distribution
   P.raylength = Rayleng; P.rayp = RayP;
end
function peak = getpeakphase(phase)
phasebin = [0:6:720]; phasenow = [phase phase+360]; phasecnt = histcounts(phasenow, phasebin); [~,ii] = max(phasecnt); peak=3+phasebin(ii); peak = mod(peak, 360);

function [currentsessname, currentEEGphase, currenttimepoint] = locateEEGphases(eeg, eegdata, finaldir, parm) 
currentsessname = []; currentEEGphase = []; currenttimepoint = [];
%%%%first select the eeg files in the same finaldir with the same recarea and band
iii = find( strcmp(eeg.general.recarea, parm.recarea) & strcmp(eeg.general.finaldir, finaldir) & strcmp(eeg.parm.band, parm.band) );
if ~isempty(iii)
   sessname = eeg.general.sessname(iii); sesstype = eeg.parm.sessType(iii); filename = eeg.general.eegfile(iii);
   timeunit = eeg.parm.timeunit(iii); buffersize = eeg.parm.buffersize(iii);
   ok2 = zeros(1, numel(iii));
   for (i = 1:numel(iii)) ok2(i) = (~isempty(strfind(lower(filename{i}), lower(parm.EEGkeyword)))); end
   jjj = find(ok2);
   if ~isempty(jjj)
      currentsessname = sessname(jjj); filename = filename(jjj); timeunit = timeunit(jjj); buffersize = buffersize(jjj); 
      currentEEGphase = cell(1, numel(jjj)); currenttimepoint = cell(1, numel(jjj));
      for (i = 1:numel(jjj))
       if exist(filename{i}, 'file')~=2
           filename{i} = strrep(filename{i}, 'I:', 'J:'); 
       end
       disp(['----------> EEG file now: ', filename{i}]);
       [currenttimepoint{i}, dat] = ReadEEGFile(filename{i}, timeunit(i), buffersize(i));
       [a,b] =size(currenttimepoint{i}); if a==1 currenttimepoint{i} = currenttimepoint{i}'; end 
       %%%%%%%%%%use hilbert transform to compute instantaneous phases (output x is a complex number)
       x = hilbert(dat); currentEEGphase{i} = mod(angle(x)*360/(2*pi), 360); %phase in (0 360)
      end
   else
       disp(['----------> no EEG files match with keyword within selected sessions: ', finaldir]); 
   end
else
   disp(['----------> no EEG file match with the date or EEG area/band: ', finaldir]); 
end

% function [currentsessname, currentEEGphase, currenttimepoint] =
% locateEEGphases(eeg, eegdata, finaldir, parm) %%%%%This one restricted to
% the selected sessions with keyword
% currentsessname = []; currentEEGphase = []; currenttimepoint = [];
% %%%%first select the eeg files in the same finaldir with the same recarea and band
% iii = find( strcmp(eeg.general.recarea, parm.recarea) & strcmp(eeg.general.finaldir, finaldir) & strcmp(eeg.parm.band, parm.band) );
% if ~isempty(iii)
%    sessname = eeg.general.sessname(iii); sesstype = eeg.parm.sessType(iii); filename = eeg.general.eegfile(iii);
%    timeunit = eeg.parm.timeunit(iii); buffersize = eeg.parm.buffersize(iii);
%    ok1 = zeros(1, numel(iii)); ok2 = zeros(1, numel(iii));
%    for (i = 1:numel(iii)) ok1(i) = filterevents(sessname{i}, sesstype{i}, parm.sesskeyword, parm.sesskeytype); end
%    for (i = 1:numel(iii)) ok2(i) = (~isempty(strfind(lower(filename{i}), lower(parm.EEGkeyword)))); end
%    jjj = find(ok1.*ok2);
%    if ~isempty(jjj)
%       currentsessname = sessname(jjj); filename = filename(jjj); timeunit = timeunit(jjj); buffersize = buffersize(jjj); 
%       currentEEGphase = cell(1, numel(jjj)); currenttimepoint = cell(1, numel(jjj));
%       for (i = 1:numel(jjj))
%        if exist(filename{i}, 'file')~=2
%            filename{i} = strrep(filename{i}, 'I:', 'J:'); 
%        end
%        disp(['----------> EEG file now: ', filename{i}]);
%        [currenttimepoint{i}, dat] = ReadEEGFile(filename{i}, timeunit(i), buffersize(i));
%        [a,b] =size(currenttimepoint{i}); if a==1 currenttimepoint{i} = currenttimepoint{i}'; end 
%        %%%%%%%%%%use hilbert transform to compute instantaneous phases (output x is a complex number)
%        x = hilbert(dat); currentEEGphase{i} = mod(angle(x)*360/(2*pi), 360); %phase in (0 360)
%       end
%    else
%        disp(['----------> no EEG files match with keyword within selected sessions: ', finaldir]); 
%    end
% else
%    disp(['----------> no EEG file match with EEG criteria: ', finaldir]); 
% end

function ok = filterevents(evName, evType, evkeyword, evkeytype)
ok = 1;
if (~isempty(evkeytype))||(~isempty(evkeyword))
    if isempty(evkeytype)
       if isempty(strfind(lower(evName), lower(evkeyword))) ok = 0; end 
    elseif isempty(evkeyword)
       if ~strncmpi(evType, evkeytype, 3) ok = 0; end 
    else
       if isempty(strfind(lower(evName), lower(evkeyword))) || (~strncmpi(evType, evkeytype, 3))
           ok = 0; 
       end 
    end
end

function evSess = identifysession(evTime, sessionname, startT, endT)
evSess = []; minevstart = min(evTime.start); maxevend = max(evTime.ent);
iii = find( (startT <= minevstart) & (endT >= maxevend) );
if (numel(iii) == 1)
    evSess = sessionname{iii};
end

