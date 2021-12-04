function [behav, bhdata] = DataManager_FindBehavPositionProp(behav, bhdata, sessind, vv)
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
if (~isfield(bhdata.event, 'LapAll1DSpeed')) bhdata.event.LapAll1DSpeed = cell(1, nsess); end %%%linearized speed at every time point of every lap --- could be negative - running backwards
if (~isfield(bhdata.event, 'LapAllVel')) bhdata.event.LapAllVel = cell(1,nsess); end %{ntraj}: velocity at each time point of a lap
if (~isfield(bhdata.event, 'LapAllHeadDir')) bhdata.event.LapAllHeadDir = cell(1,nsess); end %{ntraj}: headdir at each time point of a lap
if (~isfield(bhdata.event, 'LapOccuptime')) bhdata.event.LapOccuptime = cell(1, nsess); end %%%occupancy time at each 1D grid of the space, computed for each lap of each event
if (~isfield(bhdata.event, 'LapMed1DSpeed')) bhdata.event.LapMed1DSpeed = cell(1, nsess); end %%%{ntraj}: lap-by-lap median linearized speed 
if (~isfield(bhdata.event, 'LapMeanSpeed')) bhdata.event.LapMeanSpeed = cell(1, nsess); end %%%{ntraj}: lap-by-lap mean speed
if (~isfield(bhdata.event, 'LapTrajleng')) bhdata.event.LapTrajleng = cell(1, nsess); end %%%{ntraj}: lap-by-lap mean trajectory lenth in pixels
if (~isfield(bhdata.event, 'LapDur')) bhdata.event.LapDur = cell(1, nsess); end %%%{ntraj}: lap-by-lap mean trajectory duration in seconds
if (~isfield(bhdata.event, 'Xbin')) bhdata.event.Xbin = cell(1, nsess); end %%%{ntraj}: linearized spatial 1D bins
if (~isfield(bhdata.event, 'Xjoint')) bhdata.event.Xjoint = cell(1, nsess); end %%%{ntraj}: linearized spatial 1D joint coordinate
if (~isfield(bhdata.event, 'Occuptime')) bhdata.event.Occuptime = cell(1, nsess); end %%%occupancy time at each 1D grid of the space, computed for each event

if (~isfield(behav, 'behavior')) behav.behavior = []; end
%behav.behavior.eventname = bhdata.event.eventname;
if (~isfield(behav.behavior, 'sessDur')) behav.behavior.sessDur = cell(1, nsess); end %number: mean traj length (s), averaged across all laps and all trajectories
if (~isfield(behav.behavior, 'sessMeanSpeed')) behav.behavior.sessMeanSpeed = cell(1, nsess); end%number: mean speed (pixel/s), averaged across all laps and all trajectories
if (~isfield(behav.behavior, 'sessMedSpeed')) behav.behavior.sessMedSpeed = cell(1, nsess); end%number: median speed (pixel/s), averaged across all laps and all trajectories
if (~isfield(behav.behavior, 'evtNum')) behav.behavior.evtNum = cell(1, nsess); end %%number of laps during each run session
if (~isfield(behav.behavior, 'evtTotDur')) behav.behavior.evtTotDur = cell(1, nsess); end
if (~isfield(behav.behavior, 'evtMeanDur')) behav.behavior.evtMeanDur = cell(1, nsess); end %number: mean duration (s), averaged across all laps and all trajectories
if (~isfield(behav.behavior, 'evtMedDur')) behav.behavior.evtMedDur = cell(1, nsess); end
if (~isfield(behav.behavior, 'evtMeanSpeed')) behav.behavior.evtMeanSpeed = cell(1, nsess); end%number: mean speed (pixel or cm/s), averaged across all laps and all trajectories
if (~isfield(behav.behavior, 'evtMedSpeed')) behav.behavior.evtMedSpeed = cell(1, nsess); end
if (~isfield(behav.behavior, 'evtMean1DSpeed')) behav.behavior.evtMean1DSpeed = cell(1, nsess); end%number: median speed (pixel or cm/s) of lap-lap median 1D speeds, across all laps and all trajectories
if (~isfield(behav.behavior, 'evtMed1DSpeed')) behav.behavior.evtMed1DSpeed = cell(1, nsess); end
if (~isfield(behav.behavior, 'evtMeanTrajleng')) behav.behavior.evtMeanTrajleng = cell(1, nsess); end %number: mean traj length (pixel), averaged across all laps and all trajectories
if (~isfield(behav.behavior, 'evtMedTrajleng')) behav.behavior.evtMedTrajleng = cell(1, nsess); end 

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
    trajstart = -1; trajend = -1;
    if isfield(behav.parm, 'trajstart') trajstart = behav.parm.trajstart(i); end
    if isfield(behav.parm, 'trajend') trajend = behav.parm.trajend(i); end
   
    postimestamp = bhdata.pos.postimestamp{i}*behav.parm.timeunit(i); 
    stime = behav.general.sessstartT{i}; etime = behav.general.sessendT{i}; 
    segtime = behav.parm.sessSegTime(i); nseg = floor((etime-stime)/segtime);
    Pmarker = behav.parm.sessPmarker{i}; Fmarker = behav.parm.sessFmarker{i}; Bmarker = behav.parm.sessBmarker{i}; 
    allposmarker = behav.general.posMarker{i}; 
    posXX = []; posYY = []; ik = find(strcmp(allposmarker, Pmarker)); 
    if (numel(ik) == 1) posXX = bhdata.pos.XX{i}{ik}*behav.parm.pixelXSize(i); posYY = bhdata.pos.YY{i}{ik}*behav.parm.pixelYSize(i); end
    fXX = []; fYY = []; ik = find(strcmp(allposmarker, Fmarker)); 
    if (numel(ik) == 1) fXX = bhdata.pos.XX{i}{ik}*behav.parm.pixelXSize(i); fYY = bhdata.pos.YY{i}{ik}*behav.parm.pixelYSize(i); end
    bXX = []; bYY = []; ik = find(strcmp(allposmarker, Bmarker)); 
    if (numel(ik) == 1) bXX = bhdata.pos.XX{i}{ik}*behav.parm.pixelXSize(i); bYY = bhdata.pos.YY{i}{ik}*behav.parm.pixelYSize(i); end
    %%compute session parameters
    AllVel = ComputeV(postimestamp, posXX, posYY, 1); AllHeadDir = ComputeDir(postimestamp, fXX, bXX, fYY, bYY);
    bhdata.sess.AllVel{i} = AllVel; bhdata.sess.AllHeadDir{i} = AllHeadDir; 
    minx = min(posXX); miny = min(posYY); maxx = max(posXX); maxy = max(posYY); %%% ---This is the old way, not easy for alignment
    %%%minx = 0; miny = 0; maxx = max(posXX); maxy = max(posYY); 
    if isfield(behav.parm, 's2dminX') %%%% this is new way for better assignment
       minx = behav.parm.s2dminX(i); miny = behav.parm.s2dminY(i); maxx = behav.parm.s2dmaxX(i); maxy = behav.parm.s2dmaxY(i);
    end
    nxbin = ceil((maxx-minx)/s2dbinsize); nybin = ceil((maxy-miny)/s2dbinsize);
    if (~isempty(posXX))
        xybin{1} = minx + ((1:nxbin)-1)*s2dbinsize; 
    else
        xybin{1} = [];
    end
    if (~isempty(posYY)) 
        xybin{2} = miny + ((1:nybin)-1)*s2dbinsize;
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
            if ~isempty(posXX)
                pxnow = posXX(tii); pynow = posYY(tii);
            else
                pxnow = []; pynow = [];
            end
            occutime = OccupancyPoint(xybin{1}, xybin{2}, pxnow, pynow, 0)/framerate;
            bhdata.sess.gridSegOccuptime{i}{j} = occutime; %bhdata.sess.gridVel{k}{i}{j} = velnow; bhdata.sess.gridHeadDir{k}{i}{j} = dirnow;
    end
    behav.behavior.sessDur{i} = etime-stime; aa = AllVel(~isnan(AllVel));
    behav.behavior.sessMeanSpeed{i} = mean(aa); behav.behavior.sessMedSpeed{i} = median(aa);
    %%compute event variables
    for (j = 1:numel(evname))
        evstart = evTimes{j}.start; evend = evTimes{j}.ent; nev = numel(evstart);
        behav.behavior.evtNum{i}(j) = nev; behav.behavior.evtTotDur{i}(j) = sum(evend-evstart);
        behav.behavior.evtMeanDur{i}(j) = mean(evend-evstart); behav.behavior.evtMedDur{i}(j) = median(evend-evstart); 
        evposind = cell(1, nev); tii = evposind; allposind = [];
        for (tt = 1:nev)
             evposind{tt} = find( (postimestamp>=evstart(tt)) & (postimestamp<=evend(tt)) );
             allposind = union(allposind, evposind{tt});
        end
        if ~isempty(posXX)%%%mean speed is replaced by overall 1D speed if can be linearized
            aa = AllVel(allposind); aa = aa(~isnan(aa));
            behav.behavior.evtMeanSpeed{i}(j) = mean(aa); behav.behavior.evtMedSpeed{i}(j) = median(aa); 
        else
            behav.behavior.evtMeanSpeed{i}(j) = NaN; behav.behavior.evtMedSpeed{i}(j) = NaN; 
        end
        behav.behavior.evtMeanTrajleng{i}(j) = NaN; behav.behavior.evtMedTrajleng{i}(j) = NaN;
        behav.behavior.evtMean1DSpeed{i}(j) = NaN; behav.behavior.evtMed1DSpeed{i}(j) = NaN;
        %%%%%%%common event variables
        bhdata.event.LapDur{i}{j} = evend-evstart;
        for (tt = 1:nev)
            bhdata.event.LapAllPostimestamp{i}{j}{tt} = postimestamp(evposind{tt});
            if ~isempty(posXX)
               bhdata.event.LapAllVel{i}{j}{tt} = AllVel(evposind{tt}); 
               bhdata.event.LapAllHeadDir{i}{j}{tt} = AllHeadDir(evposind{tt}); 
            else
               bhdata.event.LapAllVel{i}{j}{tt} = NaN; 
               bhdata.event.LapAllHeadDir{i}{j}{tt} = NaN; 
            end
        end
        %%%%%%%run-type event variables
        bhdata.event.LapTrajleng{i}{j} = []; bhdata.event.LapAll1DSpeed{i}{j} = []; bhdata.event.LapMed1DSpeed{i}{j} = [];
        bhdata.event.LapOccuptime{i}{j} = []; bhdata.event.LapAllX{i}{j} = []; 
        bhdata.event.LapAllY{i}{j} = []; bhdata.event.Xbin{i}{j} = []; bhdata.event.Xjoint{i}{j} = [];
        bhdata.event.LapMeanSpeed{i}{j} = [];
        bhdata.event.Occuptime{i}{j} = []; spaceunit = (behav.parm.pixelXSize(i)+ behav.parm.pixelYSize(i))/2;
        if (strcmp(evType{j}, 'run'))  %if this is run event
                joint = [];
                for (tt = 1:numel(posltrfile))
                    if (strcmp(posltrfile{tt}, evPosltr{j}))
                        joint = posltr{tt}; 
                        joint(:,1) = joint(:,1)*behav.parm.pixelXSize(i); joint(:,2) = joint(:,2)*behav.parm.pixelYSize(i);
                        break
                    end
                end
                if (~isempty(joint)) && (~isempty(posXX))
                alloccu = []; occunow = cell(1, numel(evstart)); xnowexit = 0; lapsspeed = []; %disp(['------- get in now: ', evname{j}]);
                for (tt = 1:nev)
                    xnow = posXX(evposind{tt}); ynow = posYY(evposind{tt});
                    %[occunow{tt}, allX, allY, xbin] = findlapoccuetc(postimenow, xnow, ynow, joint, s1dbinsize); 
                    if ~isempty(xnow)
                       bhdata.event.LapTrajleng{i}{j}(tt) = findtrajleng(xnow, ynow);
                       [occunow{tt}, xbin, allX, allY, pout] = OccupancyTimeEps_binsize(joint, xnow, ynow, framerate, s1dbinsize, 0, trajstart, trajend, spaceunit);
                       bhdata.event.LapOccuptime{i}{j}{tt} = occunow{tt}; bhdata.event.LapAllX{i}{j}{tt} = allX; 
                       bhdata.event.LapAllY{i}{j}{tt} = allY; bhdata.event.Xbin{i}{j} = xbin; bhdata.event.Xjoint{i}{j} = pout;
                       timenow = postimestamp(evposind{tt}); aa = diff(allX') ./ diff(timenow); lapsspeed = [lapsspeed; aa]; %%%aa is a column vector
                       bhdata.event.LapAll1DSpeed{i}{j}{tt} = [NaN; aa]'; bhdata.event.LapMed1DSpeed{i}{j}(tt) = median(aa(~isnan(aa)));
                       bhdata.event.LapMeanSpeed{i}{j}(tt) = pout(numel(pout))/(evend(tt)-evstart(tt));
                       if (xnowexit == 0) 
                           alloccu = occunow{tt}; xnowexit = 1;
                       else
                           alloccu = alloccu + occunow{tt};
                       end
                    else
                        bhdata.event.LapOccuptime{i}{j}{tt} = []; bhdata.event.LapAllX{i}{j}{tt} = []; 
                        bhdata.event.LapAllY{i}{j}{tt} = []; bhdata.event.Xbin{i}{j} = []; bhdata.event.Xjoint{i}{j} = [];
                        bhdata.event.LapMeanSpeed{i}{j}(tt) = NaN; bhdata.event.LapTrajleng{i}{j}(tt) = NaN;
                        bhdata.event.LapAll1DSpeed{i}{j}{tt} = []; bhdata.event.LapMed1DSpeed{i}{j}(tt) = NaN;
                    end
                end
                bhdata.event.Occuptime{i}{j} = alloccu; 
                %aa = bhdata.event.LapMeanSpeed{i}{j}; behav.behavior.evtMeanSpeed{i}(j) = mean(aa(~isnan(aa)));
                aa = bhdata.event.LapTrajleng{i}{j}; aa = aa(~isnan(aa)); 
                behav.behavior.evtMeanTrajleng{i}(j) = mean(aa); behav.behavior.evtMedTrajleng{i}(j) = median(aa);
                aa = lapsspeed(~isnan(lapsspeed));
                behav.behavior.evtMean1DSpeed{i}(j) = mean(aa); behav.behavior.evtMed1DSpeed{i}(j) = median(aa);
                end
         end
    end
end

function trajleng = findtrajleng(xnow, ynow)
trajleng = 0; 
prevs = 0; %%%%if there is a jump between two time points, use the previous distance
for (i = 2:numel(xnow))
    dis = sqrt( (xnow(i)-xnow(i-1))^2 + (ynow(i)-ynow(i-1))^2 );
    if (dis < 20)
       trajleng = trajleng + dis; prevs = dis;
    else
       trajleng = trajleng + prevs;
    end
end



