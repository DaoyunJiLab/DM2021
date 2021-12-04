function [behav, bhdata] = DataManager_FindBehavParm(behav, bhdata)
%%set up all parameters for computing behavioral bhdata
%%%behav.parm.sessFmarker, behav.sessBmarker, behav.parm.sessPmarker, behav.parm.sessType, behav.parm.eventType
%%%        will be assigned in FindBehavParm

%%%%variables to assign
defaultbin = 1; %2.5; %%default spatial binsize in cm;
defaultseg = 180; %%default segment time in seconds 
nsess = numel(behav.general.datedir);
behav.parm.timeunit = 0.0001*ones(1, nsess);   %timeunit = 100us
behav.parm.framerate = NaN*ones(1,nsess);
behav.parm.sessFmarker = cell(1,nsess); behav.parm.sessBmarker = cell(1,nsess); behav.parm.sessPmarker = cell(1,nsess);
behav.parm.sessType = cell(1,nsess); 
behav.parm.eventType = cell(1,nsess); behav.parm.eventPosltr = cell(1,nsess);
behav.parm.pixelXSize = 0.25*ones(1, nsess);  %number of cm per pixel: horizontal direction 
behav.parm.pixelYSize = 0.25*ones(1, nsess); %number of cm per pixel: vertical direction
behav.parm.s1dbinsize = NaN*ones(1, nsess);  %1d grid binsize in cm: linearized track spatial binsize 
behav.parm.s2dbinsize = NaN*ones(1, nsess);   %2d track spatial binsize in pixels
behav.parm.s2dminX = zeros(1, nsess); behav.parm.s2dmaxX = 200*ones(1, nsess);
behav.parm.s2dminY = zeros(1, nsess); behav.parm.s2dmaxY = 200*ones(1, nsess); %%%hard assignments of 2D space boundaries in cm; for better alignments across sessions
behav.parm.sessSegTime = NaN*ones(1, nsess);   %duration of a segment in a session in second (3 min, 5 min ect.)
behav.parm.trajstart = 0*ones(1, nsess);   %in cm now; starting point of 1D linearized trajectory: traj outside the [start, end] points are excluded for 1D place field identtification/correlation
behav.parm.trajend = 10000*ones(1, nsess);   %in cm now; ending point of 1D linearized trajectory
%%%assign front, back, position files
frate = NaN;
for (i = 1:nsess)
     sessnow = behav.general.sessname{i};
     behav.parm.sessType{i} = findsesstype(sessnow);
     postimenow = bhdata.pos.postimestamp{i}; 
     nframe = numel(postimenow); 
     if (nframe > 0)
         frate = (nframe-1)/((postimenow(nframe) - postimenow(1))*behav.parm.timeunit(i));
     else
         frate = NaN;
     end
     markernow = behav.general.posMarker{i};
     behav.parm.sessFmarker{i} = []; behav.parm.sessBmarker{i} = []; behav.parm.sessPmarker{i} = []; 
     
     if numel(markernow) == 1 %%%if only one color
         behav.parm.sessPmarker{i} = markernow{1}; behav.parm.sessFmarker{i} = markernow{1}; behav.parm.sessBmarker{i} = markernow{1};
     elseif numel(markernow)>1
         locnow = behav.general.posLocation{i}; animnow = behav.general.posAnimalname{i};
         if (~isempty(strfind(lower(sessnow), 'obs'))) %%%if observation session
             for (j = 1:numel(markernow))
                 if (~isempty(strfind(lower(animnow{j}), 'dem'))) && (~isempty(strfind(markernow{j}, 'red'))) %%assign to green or back of demo
                     behav.parm.sessPmarker{i} = markernow{j}; break
                 elseif (~isempty(strfind(lower(animnow{j}), 'dem'))) && (~isempty(strfind(locnow{j}, 'front')))
                     behav.parm.sessPmarker{i} = markernow{j}; break
                 end
             end
             if isempty(behav.parm.sessPmarker{i}) %%if still unassigned, assign the animal name
                 for (j = 1:numel(markernow))
                 if (~isempty(strfind(lower(animnow{j}), 'dem')))
                     behav.parm.sessPmarker{i} = markernow{j}; break
                 end
                 end
             end
             for (j = 1:numel(locnow))
                 if (isempty(strfind(lower(animnow{j}), 'dem'))) && (~isempty(strfind(locnow{j}, 'front'))) %%front assign to the observer's front 
                    behav.parm.sessFmarker{i} = markernow{j}; break
                 end 
             end
             for (j = 1:numel(locnow))
                 if (isempty(strfind(lower(animnow{j}), 'dem'))) && (~isempty(strfind(locnow{j}, 'back'))) %%back assign to the observer's back 
                    behav.parm.sessBmarker{i} = markernow{j}; break
                 end 
             end
         else %%%if non-observation session
             for (j = 1:numel(markernow))
                 if (~isempty(strfind(markernow{j}, 'red'))) && (j<=numel(animnow)) && isempty(strfind(lower(animnow{j}), 'dem')) && (~isempty(animnow{j}))
                     behav.parm.sessPmarker{i} = markernow{j}; break
                 end
             end
             if isempty(behav.parm.sessPmarker{i}) %%if still unassigned, assign the animal name
                 for (j = 1:numel(markernow))
                 if (~isempty(animnow{j})) && (isempty(strfind(lower(animnow{j}), 'dem')))
                     behav.parm.sessPmarker{i} = markernow{j}; break
                 end
                 end
             end
             if isempty(behav.parm.sessPmarker{i}) %%if still unassigned, assign the animal name
                 for (j = 1:numel(markernow))
                 if (~isempty(strfind(animnow{j}, 'DemoS'))) && (~isempty(strfind(lower(markernow{j}), 'red')))
                     behav.parm.sessPmarker{i} = markernow{j}; break
                 end
                 end
             end
             for (j = 1:numel(locnow))
                 if (~isempty(strfind(locnow{j}, 'front'))) && (isempty(strfind(lower(animnow{j}), 'dem'))) && (~isempty(animnow{j}))
                     behav.parm.sessFmarker{i} = markernow{j}; break
                 end
             end
             for (j = 1:numel(locnow))
                 if (~isempty(strfind(locnow{j}, 'back'))) && (isempty(strfind(lower(animnow{j}), 'dem'))) && (~isempty(animnow{j}))
                     behav.parm.sessBmarker{i} = markernow{j}; break
                 end
             end
         end
     end
          
     behav.parm.s1dbinsize(i) = defaultbin; behav.parm.s2dbinsize(i) = defaultbin;
     behav.parm.sessSegTime(i) = defaultseg; behav.parm.framerate(i) = frate;
     nev = numel(behav.general.eventname{i});
     for (j = 1:nev)
         evnamenow = behav.general.eventname{i}{j};
         [eTnow, ePnow] = findeventtypesess(evnamenow, bhdata.pos.ltrfilename{i});
         behav.parm.eventType{i}{j} = eTnow; behav.parm.eventPosltr{i}{j} = ePnow; 
     end
end

function sessType = findsesstype(sessnow)
sessType = 'others';
if (~isempty(strfind(sessnow, 'open'))) | (~isempty(strfind(sessnow, 'field'))) | (~isempty(strfind(sessnow, 'OF'))) ...
        | (~isempty(strfind(sessnow, 'OP'))) | (~isempty(strfind(sessnow, 'platform'))) | (~isempty(strfind(sessnow, 'circ')))...
        | (~isempty(strfind(sessnow, 'square'))) | (~isempty(strfind(sessnow, 'Open')))...
        | (~isempty(strfind(sessnow, 'QS'))) | (~isempty(strfind(sessnow, 'QW'))) | (~isempty(strfind(sessnow, 'QN')))...
        | (~isempty(strfind(sessnow, 'QE'))) | (~isempty(strfind(sessnow, 'QX'))) | (~isempty(strfind(lower(sessnow), 'box')))
    sessType = 'open';
elseif (~isempty(strfind(sessnow, 'track'))) | (~isempty(strfind(sessnow, 'linear'))) | (~isempty(strfind(sessnow, 'LT')))...
        | (~isempty(strfind(sessnow, 'Linear'))) | (~isempty(strfind(sessnow, 'Track'))) | (~isempty(strfind(sessnow, 'TRK')))...
        | (~isempty(strfind(lower(sessnow), 'run'))) | (~isempty(strfind(lower(sessnow), 'obs'))) 
    sessType = 'linear';
elseif (~isempty(strfind(sessnow, 'sleep'))) | (~isempty(strfind(sessnow, 'Sleep'))) | (~isempty(strfind(sessnow, 'rest')))...
        | (~isempty(strfind(sessnow, 'Rest'))) | (~isempty(strfind(sessnow, 'SLP')))
    sessType = 'sleep';
end
    
function [evType, evltr] = findeventtypesess(evn, ltrfiles)
evType = 'others'; evSess = 'notfound'; evltr = 'notfound';
if ( (~isempty(strfind(evn, 'cw'))) | (~isempty(strfind(evn, 'acw'))) | (~isempty(strfind(evn, 'CW'))) | (~isempty(strfind(evn, 'ACW')))...
        | (~isempty(strfind(evn, 'leftright'))) | (~isempty(strfind(evn, 'rightleft'))) | (~isempty(strfind(evn, 'Leftright')))...
        | (~isempty(strfind(evn, 'Rightleft'))) | (~isempty(strfind(evn, 'ncw'))) | (~isempty(strfind(evn, 'NCW')))...
        | (~isempty(strfind(evn, 'RL'))) | (~isempty(strfind(evn, 'LR'))) | (~isempty(strfind(lower(evn), 'run')))...
        | (~isempty(strfind(evn, 'lefttoright'))) | (~isempty(strfind(evn, 'righttoleft'))) ...
        | (~isempty(strfind(evn, 'CenterRight'))) | (~isempty(strfind(evn, 'RightCenter'))) ...
        | (~isempty(strfind(evn, 'CenterLeftt'))) | (~isempty(strfind(evn, 'LeftCenter'))) )...
    & ( (isempty(strfind(evn, 'ripple'))) & (isempty(strfind(evn, 'spindle'))) & (isempty(strfind(evn, 'full'))) ...
        & (isempty(strfind(evn, 'session'))) & (isempty(strfind(evn, 'allep'))) & (isempty(strfind(evn, 'Wrong')))  )
   evType = 'run';
elseif (~isempty(strfind(evn, 'stop'))) | (~isempty(strfind(evn, 'Stop'))) | (~isempty(strfind(evn, 'eat'))) | (~isempty(strfind(evn, 'Eat')))...
        | (~isempty(strfind(evn, 'leftleft'))) | (~isempty(strfind(evn, 'rightright'))) | (~isempty(strfind(evn, 'Leftleft')))...
        | (~isempty(strfind(evn, 'Rightright')))
    evType = 'stop';
% elseif (~isempty(strfind(evn, 'sleep'))) | (~isempty(strfind(evn, 'Sleep')))...
%         | (~isempty(strfind(evn, 'rest'))) | (~isempty(strfind(evn, 'Rest')))
%     evType = 'sleep';
elseif (~isempty(strfind(evn, 'REM'))) | (~isempty(strfind(evn, 'rem')))
    evType = 'rem';
elseif (~isempty(strfind(evn, 'SWS'))) | (~isempty(strfind(evn, 'sws')))
    evType = 'sws';    
end
if (~isempty(strfind(evn, 'ripp')))
    evType = 'ripple';
end
if (~isempty(strfind(evn, 'theta')))
    evType = 'theta';
end
if (~isempty(strfind(lower(evn), 'stop'))) | (~isempty(strfind(evn, 'leftleft'))) | (~isempty(strfind(evn, 'rightright'))) %%%stop keyword override run keyword
    evType = 'stop';
end


ntotstr = 0;
for (i = 1:numel(ltrfiles))
    [str, ~] = strtok(ltrfiles{i}, '.');
    if contains(lower(evn), lower(str))
        nstrnow = numel(str); 
        if (nstrnow > ntotstr)
            evltr = ltrfiles{i}; ntotstr = nstrnow;
        end
    end
end
if ntotstr==0
   evnnow = strrep(evn, 'demo', 'obsv');
   for (i = 1:numel(ltrfiles))
       [str, ~] = strtok(ltrfiles{i}, '.');
       if contains(lower(evnnow), lower(str))
        nstrnow = numel(str); 
        if (nstrnow > ntotstr)
            evltr = ltrfiles{i}; ntotstr = nstrnow;
        end
       end
   end
end
    
    
    