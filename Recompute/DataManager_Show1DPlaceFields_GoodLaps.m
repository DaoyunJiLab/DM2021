function DataManager_Show1DPlaceFields_GoodLaps
%%Plot the computed bhdata
hf = gcbf; pinfo = getappdata(hf, 'pinfo'); data = getappdata(hf, 'data');
hgroup = getappdata(hf, 'hgroup'); groupselection = getappdata(hgroup, 'selection');
hfield = getappdata(hf, 'hfield'); plotparm = getappdata(hf, 'plotparm');
hspike = getappdata(hf, 'hspike'); spikeselection = getappdata(hspike, 'selection'); 
spikeind = find(spikeselection == 1); %%%%%spikeind contains what will be plotted
behav = []; bhdata = [];
if (plotparm.linkbehav == 0);
    disp(['--------> no behavioral data linked']);
else
    behav = getappdata(hf, 'behav'); bhdata = getappdata(hf, 'bhdata');
    plot1dGRfields(pinfo, data, behav, bhdata, spikeind, 'Show1Dfields');
end
disp('************************');

function plot1dGRfields(pinfo, data, behav, bhdata, spikeind, tag)
if ~isfield(pinfo, 'field')
    msgbox('Fields not defined in the data base');
else
    if (numel(spikeind)>10)
        disp('----------> Too many cells selected; Only plot fields for the first 10 cells');
        spikeind = spikeind(1:10);
    end
    for (tt = 1:numel(spikeind))
        i = spikeind(tt); %%for the current cell, get event rate maps and defined fields
        MinLapMeanSpeed = 0; if isfield(pinfo.parm, 'minLapSpeed') MinLapMeanSpeed = pinfo.parm.minLapSpeed(i); end
        clname = pinfo.general.clname{i};
        evName = pinfo.general.eventname{i}; evTime = data.events.eventtimes{i}; evType = pinfo.parm.eventtype{i}; 
        D1rate = data.GRfield.evt1DRateMaps{i}; finaldirnow = pinfo.general.finaldir{i};
        spiketime = pinfo.parm.timeunit(i)*data.spike.spiketime{i};
        for (j = 1:numel(evTime))
             xbin = 5*(1:numel(D1rate{j})); xjoint = []; nlap = 0; %default xbin
             if (strcmp(evType{j}, 'run')) && (~isempty(D1rate{j}))
                if (~isempty(behav)) & (~isempty(bhdata)) %%%get xbin from behavioral data
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
                   if (numel(posid)~=1)|(numel(evid)~=1)
                      disp(['-------------> field not computed: no or more than 1 positon/event files match the session: ', finaldirnow, '___', evName{j}]);
                   else
                      xbin = bhdata.event.Xbin{posid}{evid}; xjoint = bhdata.event.Xjoint{posid}{evid};
                   end
                   %%%%get lap-by-lap spike rasters
                   lapmeanspeed = bhdata.event.LapMeanSpeed{posid}{evid}; ttjj = find(lapmeanspeed >= MinLapMeanSpeed);
                   nlap = numel(ttjj); spikeX = cell(1, nlap);
                   for tj = 1:nlap
                        jnow = ttjj(tj);
                        lappostime = bhdata.event.LapAllPostimestamp{posid}{evid}{jnow}; 
                        lapx = bhdata.event.LapAllX{posid}{evid}{jnow};
                        evok.start = evTime{j}.start(jnow); evok.ent = evTime{j}.ent(jnow);
                        [spiketimenow, epid] = SpikeEventFilter(spiketime, evok);
                        spikeX{tj} = findspikex(lappostime, lapx, spiketimenow);
                   end
                end
                %%%%%%%plot spike rasters
                hg = figure('Name', strcat(tag, '---', clname)); mr = max(D1rate{j}); if isnan(mr) || (mr ==0 ) mr = 1; end
                posvec = [0.08 0.05 0.90 0.90];
                for (tti=1:nlap) %first tick trace is on the top.
                     mtickposvec{tti} = [posvec(1) posvec(2)+posvec(4)/3+(nlap-tti)*posvec(4)*2/3/nlap posvec(3) posvec(4)*2/3/nlap];    % multiple tick traces
                     hmtickaxes(tti) = axes('Parent', hg, 'Units', 'normalized', 'Position', mtickposvec{tti}, 'Visible', 'off', 'NextPlot', 'add');
                     xlim ([min(xbin)-10 max(xbin)+10]); %plot x-axis as linearized track
                     ylim ([0 1.2]);
                end
                for (tti=1:nlap) %for each episode
                    if (~isempty(spikeX{tti}))
                    for (ttj=1:numel(spikeX{tti}))
                        line([spikeX{tti}(ttj) spikeX{tti}(ttj)], [0 1], 'Parent', hmtickaxes(tti), 'LineWidth', 2, 'Color', [0 0 1]) % plot spikes as ticks
                    end
                    end
                end
                %%%plot rate maps and xjoint
                rateposvec = [posvec(1) posvec(2) posvec(3) posvec(4)/3]; %mr
                hax = axes('Parent', hg, 'Units', 'normalized', 'Position', rateposvec, 'NextPlot', 'add', 'XLim', [min(xbin)-10 max(xbin)+10], 'YLim', [-0.2 1.2]*mr);
                xlabel ('X (pixel)'); ylabel('Firing rate (Hz)'); 
                line(xbin, D1rate{j}, 'Parent', hax, 'LineWidth', 4, 'Color', [0 0 0]);
                for (k = 1:numel(xjoint))
                    line([xjoint(k) xjoint(k)], [0 max(D1rate{j})], 'Parent', hax, 'LineWidth', 0.5, 'Color', [0.5 0.5 0.5]);
                end
                %%%find and plot fields belong to this event (trajectory)
                pfev = pinfo.GRfield.PF1Devt{i}; pfid = find(strcmp(pfev, evName{j})); baserate = pinfo.GRfield.PF1DBaseRate{i}(j);
                for (k = 1:numel(pfid))
                    sbound = pinfo.GRfield.PF1DBoundStart{i}(pfid(k)); ebound = pinfo.GRfield.PF1DBoundEnd{i}(pfid(k)); 
                    sloc = pinfo.GRfield.PF1DLocStartX{i}(pfid(k)); eloc = pinfo.GRfield.PF1DLocEndX{i}(pfid(k));
                    ploc = pinfo.GRfield.PF1DLocPeakX{i}(pfid(k)); cloc = pinfo.GRfield.PF1DLocComX{i}(pfid(k));
                    prate = pinfo.GRfield.PF1DInPeakrate{i}(pfid(k));
                    line([sbound sbound], [0 prate], 'Parent', hax, 'LineWidth', 1, 'Color', [0 0 1]);
                    line([ebound ebound], [0 prate], 'Parent', hax, 'LineWidth', 1, 'Color', [0 0 1]);
                    line([sbound ebound], [prate prate], 'Parent', hax, 'LineWidth', 1, 'Color', [1 0 0]);
                    line([sloc sloc], [0 prate], 'Parent', hax, 'LineWidth', 1, 'Color', [1 0 0]);
                    line([eloc eloc], [0 prate], 'Parent', hax, 'LineWidth', 1, 'Color', [1 0 0]);
                    line([ploc ploc], [0 prate], 'Parent', hax, 'LineWidth', 1, 'Color', [0 1 0]);
                    line([cloc cloc], [0 prate], 'Parent', hax, 'LineWidth', 1, 'Color', [0 1 1]);
                end
                line([xbin(1) xbin(numel(xbin))], [baserate baserate], 'Parent', hax, 'LineWidth', 1, 'Color', [0 0 0]);
                str = [evName{j} '; Number of fields identified: ' num2str(numel(pfid))];
                text('Interpreter', 'none', 'Parent', hax, 'String', str, 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.02 0.94]);
                %text('Interpreter', 'none', 'Parent', hax, 'String', str, 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.05 0.92]);    
            end
        end
    end
end

function evSess = identifysession(evTime, sessionname, startT, endT)
evSess = []; minevstart = min(evTime.start); maxevend = max(evTime.ent);
iii = find( (startT <= minevstart) & (endT >= maxevend) );
if (numel(iii) == 1)
    evSess = sessionname{iii};
end

function spikeX = findspikex(postime, posx, spiketime)
spikeX = zeros(size(spiketime));
[postime, iii] = sort(postime); posx = posx(iii);
lastpoint = 1;  %this only for saving time
for (j = 1:numel(spiketime))
    for (k = lastpoint:numel(postime)-1) %find corresponding time in position data
         if (postime(k) <= spiketime(j)) && (postime(k+1) > spiketime(j)) 
             spikeX(j) = posx(k); lastpoint = k; break; 
         end
    end
end


