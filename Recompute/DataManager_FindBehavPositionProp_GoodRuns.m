function [behav, bhdata] = DataManager_FindBehavPositionProp_GoodRuns(behav, bhdata, sessind, minRunspeed, meanLapSpeed, vv)
%compute common behavior parameters

%field to assign
nsess = numel(behav.general.datedir);
if (~isfield(bhdata, 'sess')) bhdata.sess = []; end
if (~isfield(bhdata.sess, 'AllVel')) bhdata.sess.AllVel = cell(1, nsess); end %%%velocity at each time points of each session
if (~isfield(bhdata.sess, 'AllHeadDir')) bhdata.sess.AllHeadDir = cell(1, nsess); end %%%Head dir at each time points of each session
if (~isfield(bhdata.sess, 'gridXYbin')) bhdata.sess.gridXYbin = cell(1, nsess); end %%%(nsess} spatial 2D bins; {1} = xgrid, {2} = ygrid
%                            %%%but content shound be arranged in [j,i] while j index ygrid and i index xgrid (for pcolor plot requirement).
if (~isfield(bhdata.sess, 'gridOccuptime')) bhdata.sess.gridOccuptime = cell(1, nsess); end %%%occupancy time at each 2D grid of the space, computed for each session
if (~isfield(bhdata.sess, 'gridSegOccuptime')) bhdata.sess.gridSegOccuptime = cell(1, nsess); end %%%occupancy time at each 2D grid of the space, computed for each segment of each session
%bhdata.sess.gridVel = cell(1, nsess); %%%Velocity at each 2D grid of the space, computed for each session
%bhdata.sess.gridSegVel = cell(1, nsess); %%%Velocity at each 2D grid of the space, computed for each segment of each session
%bhdata.sess.gridHeadDir = cell(1, nsess); %%%HeadDir at each 2D grid of the space, computed for each session
%bhdata.sess.gridSegHeadDir = cell(1, nsess); %%%HeadDir at each 2D grid of the space, computed for each segment of each session
if (~isfield(bhdata, 'event')) bhdata.event = []; end
if (~isfield(bhdata.event, 'LapAllPostimestamp')) bhdata.event.LapAllPostimestamp = cell(1,nsess); end %{ntraj}: all timestamps in a lap 
if (~isfield(bhdata.event, 'LapAllX')) bhdata.event.LapAllX = cell(1,nsess); end %{ntraj}: linearized X at each time point of a lap
if (~isfield(bhdata.event, 'LapAllY')) bhdata.event.LapAllY = cell(1,nsess); end %{ntraj}: linearized lateral position Y at each time point of a lap
if (~isfield(bhdata.event, 'LapAllVel')) bhdata.event.LapAllVel = cell(1,nsess); end %{ntraj}: velocity at each time point of a lap
if (~isfield(bhdata.event, 'LapAllHeadDir')) bhdata.event.LapAllHeadDir = cell(1,nsess); end %{ntraj}: headdir at each time point of a lap
if (~isfield(bhdata.event, 'LapOccuptime')) bhdata.event.LapOccuptime = cell(1, nsess); end %%%occupancy time at each 1D grid of the space, computed for each lap of each event
if (~isfield(bhdata.event, 'LapMeanSpeed')) bhdata.event.LapMeanSpeed = cell(1, nsess); end %%%{ntraj}: lap-by-lap mean speed
if (~isfield(bhdata.event, 'LapTrajleng')) bhdata.event.LapTrajleng = cell(1, nsess); end %%%{ntraj}: lap-by-lap mean trajectory lenth in pixels
if (~isfield(bhdata.event, 'LapDur')) bhdata.event.LapDur = cell(1, nsess); end %%%{ntraj}: lap-by-lap mean trajectory duration in seconds
if (~isfield(bhdata.event, 'Xbin')) bhdata.event.Xbin = cell(1, nsess); end %%%{ntraj}: linearized spatial 1D bins
if (~isfield(bhdata.event, 'Xjoint')) bhdata.event.Xjoint = cell(1, nsess); end %%%{ntraj}: linearized spatial 1D joint coordinate
if (~isfield(bhdata.event, 'Occuptime')) bhdata.event.Occuptime = cell(1, nsess); end %%%occupancy time at each 1D grid of the space, computed for each event

if (~isfield(behav, 'behavior')) behav.behavior = []; end
%behav.behavior.eventname = bhdata.event.eventname;
if (~isfield(behav.behavior, 'runlapnum')) behav.behavior.runlapnum = cell(1, nsess); end %%number of laps during each run session
if (~isfield(behav.behavior, 'totalrunlap')) behav.behavior.totalrunlap = cell(1, nsess); end %%number of total laps durint the run session
if (~isfield(behav.behavior, 'meanRunSpeed')) behav.behavior.meanRunSpeed = cell(1, nsess); end%number: mean speed (pixel/s), averaged across all laps and all trajectories
if (~isfield(behav.behavior, 'meanRunTrajleng')) behav.behavior.meanRunTrajleng = cell(1, nsess); end %number: mean traj length (pixel), averaged across all laps and all trajectories
if (~isfield(behav.behavior, 'meanRunDur')) behav.behavior.meanRunDur = cell(1, nsess); end %number: mean traj length (s), averaged across all laps and all trajectories
if (~isfield(behav.behavior, 'meanStopDur')) behav.behavior.meanStopDur = cell(1, nsess); end %number: mean traj length (s), averaged across all laps and all trajectories
if (~isfield(behav.behavior, 'meanSleepDur')) behav.behavior.meanSleepDur = cell(1, nsess); end %number: mean traj length (s), averaged across all laps and all trajectories
if (~isfield(behav.behavior, 'meanSWSDur')) behav.behavior.meanSWSDur = cell(1, nsess); end %number: mean traj length (s), averaged across all laps and all trajectories
if (~isfield(behav.behavior, 'meanREMDur')) behav.behavior.meanREMDur = cell(1, nsess); end %number: mean traj length (s), averaged across all laps and all trajectories

%now work on each day
for ttp = 1:numel(sessind)
    i = sessind(ttp); %%%current session
    disp(strcat('-----> compute position data properties ---', behav.general.sessID{i}));
    %%%common info
    posltrfile = bhdata.pos.ltrfilename{i}; posltr = bhdata.pos.posltr{i};
    %%%check in the data
    evname = behav.general.eventname{i}; evTimes = bhdata.event.eventtimes{i};  
    evType = behav.parm.eventType{i}; evPosltr = behav.parm.eventPosltr{i}; 
    framerate = behav.parm.framerate(i); s2dbinsize = behav.parm.s2dbinsize(i); s1dbinsize = behav.parm.s1dbinsize(i);
   
    postimestamp = bhdata.pos.postimestamp{i}*behav.parm.timeunit(i); 
    stime = behav.general.sessstartT{i}; etime = behav.general.sessendT{i}; 
    segtime = behav.parm.sessSegTime(i); nseg = floor((etime-stime)/segtime);
    Pmarker = behav.parm.sessPmarker{i}; Fmarker = behav.parm.sessFmarker{i}; Bmarker = behav.parm.sessBmarker{i}; 
    allposmarker = behav.general.posMarker{i}; 
    posXX = []; posYY = []; ik = find(strcmp(allposmarker, Pmarker)); 
    if (numel(ik) == 1) posXX = bhdata.pos.XX{i}{ik}; posYY = bhdata.pos.YY{i}{ik}; end
    fXX = []; fYY = []; ik = find(strcmp(allposmarker, Fmarker)); 
    if (numel(ik) == 1) fXX = bhdata.pos.XX{i}{ik}; fYY = bhdata.pos.YY{i}{ik}; end
    bXX = []; bYY = []; ik = find(strcmp(allposmarker, Bmarker)); 
    if (numel(ik) == 1) bXX = bhdata.pos.XX{i}{ik}; bYY = bhdata.pos.YY{i}{ik}; end
    %%compute session parameters
    AllVel = ComputeV(postimestamp, posXX, posYY, 1); AllHeadDir = ComputeDir(postimestamp, fXX, bXX, fYY, bYY);
    bhdata.sess.AllVel{i} = AllVel; bhdata.sess.AllHeadDir{i} = AllHeadDir; 
    nxbin = ceil((max(posXX)-min(posXX))/s2dbinsize); nybin = ceil((max(posYY)-min(posYY))/s2dbinsize);
    if (~isempty(posXX))
        xybin{1} = min(posXX) + ((1:nxbin)-1)*s2dbinsize; 
    else
        xybin{1} = [];
    end
    if (~isempty(posYY)) 
        xybin{2} = min(posYY) + ((1:nybin)-1)*s2dbinsize;
    else
        xybin{2} = [];
    end
    %[occutime, velnow, dirnow] = find2Doccutimeetc(postimestamp, posXX, posYY, AllVel, AllHeadDir, xybin);
    occutime = OccupancyPoint(xybin{1}, xybin{2}, posXX, posYY,0)/framerate;
    bhdata.sess.gridXYbin{i} = xybin; bhdata.sess.gridOccuptime{i} = occutime;
    %bhdata.sess.gridVel{i} = velnow; bhdata.sess.gridHeadDir{i} = dirnow;
    for (j = 1:nseg)
            stnow = stime + (j-1)*segtime; etnow = stnow + segtime; if (etnow>etime) etnow = etime; end
            tii = find( (postimestamp>=stnow)&(postimestamp<etnow) );
            %[occutime, velnow, dirnow] = find2Doccutimeetc(postimestamp(tii), posXX(tii), posYY(tii), AllVel(tii), AllHeadDir(tii), xybin);
            occutime = OccupancyPoint(xybin{1}, xybin{2}, posXX(tii), posYY(tii), 0)/framerate;
            bhdata.sess.gridSegOccuptime{i}{j} = occutime; %bhdata.sess.gridVel{k}{i}{j} = velnow; bhdata.sess.gridHeadDir{k}{i}{j} = dirnow;
    end
    %%compute event variables
    runlapnow = [];
    for (j = 1:numel(evname))
            if (strcmp(evType{j}, 'run'))  %if this is run event
                joint = [];
                for (tt = 1:numel(posltrfile))
                    if (strcmp(posltrfile{tt}, evPosltr{j}))
                        joint = posltr{tt}; break
                    end
                end
                evstart = evTimes{j}.start; evend = evTimes{j}.ent; runlapnow(j) = numel(evstart);
                occunow = cell(1, numel(evstart)); %disp(['------- get in now: ', evname{j}]);
                for (tt = 1:numel(evstart))
                    tii = find( (postimestamp>=evstart(tt)) & (postimestamp<=evend(tt)) );
                    postimenow = postimestamp(tii); xnow = posXX(tii); ynow = posYY(tii);
                    bhdata.event.LapAllPostimestamp{i}{j}{tt} = postimenow;
                    bhdata.event.LapAllVel{i}{j}{tt} = AllVel(tii); bhdata.event.LapAllHeadDir{i}{j}{tt} = AllHeadDir(tii);
                    bhdata.event.LapDur{i}{j}(tt) = evend(tt)-evstart(tt);
                    bhdata.event.LapTrajleng{i}{j}(tt) = findtrajleng(xnow, ynow);
                    %[occunow{tt}, allX, allY, xbin] = findlapoccuetc(postimenow, xnow, ynow, joint, s1dbinsize); 
                    [occunow{tt}, xbin, allX, allY, pout] = OccupancyTimeEps_binsize(joint, xnow, ynow, framerate, s1dbinsize, 0);
                    bhdata.event.LapOccuptime{i}{j}{tt} = occunow{tt}; bhdata.event.LapAllX{i}{j}{tt} = allX; 
                    bhdata.event.LapAllY{i}{j}{tt} = allY; bhdata.event.Xbin{i}{j} = xbin; bhdata.event.Xjoint{i}{j} = pout;
                    bhdata.event.LapMeanSpeed{i}{j}(tt) = pout(numel(pout))/(evend(tt)-evstart(tt));
                end
                alloccu = zeros(size(occunow{1}));
                for (tt = 1:numel(evstart))
                    alloccu = alloccu + occunow{tt};
                end
                bhdata.event.Occuptime{i}{j} = alloccu; behav.behavior.runlapnum{i}{j} = runlapnow(j);
            end
    end
    if (~isempty(runlapnow)) behav.behavior.totalrunlap{i} = sum(runlapnow); end

    runevid = find( strcmp(evType, 'run') ); 
    if (~isempty(runevid))
       runsp = []; rundur = []; runleng = [];
       for (j = 1:numel(runevid))
           runsp = [runsp bhdata.event.LapMeanSpeed{i}{runevid(j)}]; rundur = [rundur bhdata.event.LapDur{i}{runevid(j)}];
           runleng = [runleng bhdata.event.LapTrajleng{i}{runevid(j)}];
       end
       behav.behavior.meanRunSpeed{i} = mean(runsp); behav.behavior.meanRunTrajleng{i} = mean(runleng);
       behav.behavior.meanRunDur{i} = mean(rundur); 
    end
    
    stopevid = find( strcmp(evType, 'stop') ); 
    if (~isempty(stopevid))
       stopdur = [];
       for (j = 1:numel(stopevid))
           stopdur = [stopdur evTimes{stopevid(j)}.ent-evTimes{stopevid(j)}.start];
       end
       behav.behavior.meanStopDur{i} = mean(stopdur);
    end
    
%     stopevid = find( strcmp(evType, 'sleep') );
%     if (~isempty(stopevid))
%        stopdur = [];
%        for (j = 1:numel(stopevid))
%            stopdur = [stopdur evTimes{stopevid(j)}.ent-evTimes{stopevid(j)}.start];
%        end
%        behav.behavior.meanSleepDur{i} = mean(stopdur);
%     end
    
    stopevid = find( strcmp(evType, 'sws') ); 
    if (~isempty(stopevid))
       stopdur = [];
       for (j = 1:numel(stopevid))
           stopdur = [stopdur evTimes{stopevid(j)}.ent-evTimes{stopevid(j)}.start];
       end
       behav.behavior.meanSWSDur{i} = mean(stopdur);
    end
    
    stopevid = find( strcmp(evType, 'rem') ); 
    if (~isempty(stopevid))
       stopdur = [];
       for (j = 1:numel(stopevid))
           stopdur = [stopdur evTimes{stopevid(j)}.ent-evTimes{stopevid(j)}.start];
       end
       behav.behavior.meanREMDur{i} = mean(stopdur);
    end
end

function trajleng = findtrajleng(xnow, ynow)
trajleng = 0;
for (i = 2:numel(xnow))
    trajleng = trajleng + sqrt( (xnow(i)-xnow(i-1))^2 + (ynow(i)-ynow(i-1))^2 );
end



