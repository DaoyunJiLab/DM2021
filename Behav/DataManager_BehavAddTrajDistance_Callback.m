function DataManager_BehavAddTrajDistance_Callback
%compute common behavior parameters
hf = gcbf; behav = getappdata(hf, 'behav'); bhdata = getappdata(hf, 'bhdata');
hgroup = getappdata(hf, 'hgroup'); groupselection = getappdata(hgroup, 'selection');
grpind = find(groupselection == 1); cellind = bhdata.grouplist.groupindex(grpind);
sessind = [];
for (i = 1:numel(grpind))
    sessind = union(sessind, cellind{i});
end

%field to assign
nsess = numel(behav.general.datedir);
if (~isfield(bhdata, 'sess')) bhdata.sess = []; end
if (~isfield(bhdata.sess, 'AllTrajDistance')) bhdata.sess.AllTrajDistance = cell(1, nsess); end
%if (~isfield(bhdata.sess, 'AllVel')) bhdata.sess.AllVel = cell(1, nsess); end %%%velocity at each time points of each session

%%%%%%added later to correct event traj length
if (~isfield(bhdata.event, 'LapTrajleng')) bhdata.event.LapTrajleng = cell(1, nsess); end %%%{ntraj}: lap-by-lap mean trajectory lenth in pixels
if (~isfield(behav.behavior, 'evtMeanTrajleng')) behav.behavior.evtMeanTrajleng = cell(1, nsess); end %number: mean traj length (pixel), averaged across all laps and all trajectories
if (~isfield(bhdata.event, 'LapTravelleng')) bhdata.event.LapTravelleng = cell(1, nsess); end %%%{ntraj}: lap-by-lap mean travel lenth in pixels
if (~isfield(behav.behavior, 'evtMeanTravelleng')) behav.behavior.evtMeanTravelleng = cell(1, nsess); end
%now work on each day
for ttp = 1:numel(sessind)
    i = sessind(ttp); %%%current session
    disp(strcat('-----> compute trajectory distance ---', behav.general.sessID{i}));
    %%%common info
    postimestamp = bhdata.pos.postimestamp{i}*behav.parm.timeunit(i); 
    Pmarker = behav.parm.sessPmarker{i}; 
    allposmarker = behav.general.posMarker{i}; 
    posXX = []; posYY = []; ik = find(strcmp(allposmarker, Pmarker)); 
    if (numel(ik) == 1) posXX = bhdata.pos.XX{i}{ik}; posYY = bhdata.pos.YY{i}{ik}; end
    %%compute session trajectory
    if (~isempty(posXX))
       D = computetrajdistance(postimestamp, posXX, posYY); 
       bhdata.sess.AllTrajDistance{i} = D; 
       evname = behav.general.eventname{i}; evTimes = bhdata.event.eventtimes{i};  
       evType = behav.parm.eventType{i}; %evPosltr = behav.parm.eventPosltr{i}; posltrfile = bhdata.pos.ltrfilename{i};
       %%compute event traj length
       for (j = 1:numel(evname))
           evstart = evTimes{j}.start; evend = evTimes{j}.ent; nev = numel(evstart);
           behav.behavior.evtMeanTrajleng{i}(j) = NaN;
           if (strcmp(evType{j}, 'run'))  %if this is run event
               evposind = cell(1, nev); 
               for (tt = 1:nev)
                   evposind{tt} = find( (postimestamp>=evstart(tt)) & (postimestamp<=evend(tt)) );
               end
               xnow = posXX(evposind{tt}); ynow = posYY(evposind{tt});
               for (tt = 1:nev)
                   bhdata.event.LapTrajleng{i}{j}(tt) = findtrajleng(xnow, ynow);
                   bhdata.event.LapTravelleng{i}{j}(tt) = findtravelleng(bhdata.event.LapAllX{i}{j}{tt});
               end
               behav.behavior.evtMeanTrajleng{i}(j) = mean(bhdata.event.LapTrajleng{i}{j});
               behav.behavior.evtMeanTravelleng{i}(j) = mean(bhdata.event.LapTravelleng{i}{j});
           end
       end
       
    else
        disp('----------> warning: position data not found');
    end
end
setappdata(hf, 'behav', behav); setappdata(hf, 'bhdata', bhdata);
disp('***************');

function D = computetrajdistance(postimestamp, xnow, ynow)
D = zeros(size(postimestamp)); D(1) = 0; prevs = 0;
for (i = 2:numel(postimestamp))
    dis = sqrt( (xnow(i)-xnow(i-1))^2 + (ynow(i)-ynow(i-1))^2 );
    if (dis < 20)
        D(i) = D(i-1) + dis; prevs = dis;
    else
        D(i) = D(i-1) + prevs;
    end
end

function trajleng = findtrajleng(xnow, ynow)
trajleng = 0; 
prevs = 0; %%%%if there is a jump between two time points, use the previous distance
for (i = 2:numel(xnow))
    dis = sqrt( (xnow(i)-xnow(i-1))^2 + (ynow(i)-ynow(i-1))^2 );
    if (dis < 20);
       trajleng = trajleng + dis; prevs = dis;
    else
       trajleng = trajleng + prevs;
    end
end

function travleng = findtravelleng(xnow)
travleng = 0; 
prevs = 0; %%%%if there is a jump between two time points, use the previous distance
for (i = 2:numel(xnow))
    dis = abs(xnow(i) - xnow(i-1));
    if (dis < 20);
       travleng = travleng + dis; prevs = dis;
    else
       travleng = travleng + prevs;
    end
end
%disp(['--------------> travel length is now: ', num2str(travleng)]); 
