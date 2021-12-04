function DataManager_BehavShowBhdata_Callback
%%Plot the computed bhdata
hf = gcbf; behav = getappdata(hf, 'behav'); bhdata = getappdata(hf, 'bhdata');
hgroup = getappdata(hf, 'hgroup'); groupselection = getappdata(hgroup, 'selection');
hfield = getappdata(hf, 'hfield'); plotparm = getappdata(hf, 'plotparm');
hspike = getappdata(hf, 'hspike'); spikeselection = getappdata(hspike, 'selection'); 
dateind = find(spikeselection == 1); %%%%%dateind here is actually sessind (session IDs)

tagmark = get(gcbo, 'Tag'); 
if (strcmp(tagmark, 'evlinearx'))
    plotevlinearx(behav, bhdata, dateind, 'evLinearX');
elseif (strcmp(tagmark, 'evlineary'))
    plotevlinearx(behav, bhdata, dateind, 'evLinearY');
elseif (strcmp(tagmark, 'evvel'))
    plotevVelDir(behav, bhdata, dateind, 'evVel');
elseif (strcmp(tagmark, 'evheaddir'))
    plotevVelDir(behav, bhdata, dateind, 'evHeaddir');
elseif (strcmp(tagmark, 'sesstraj'))
    plotsesstraj(behav, bhdata, dateind, 'sessTraj');
elseif (strcmp(tagmark, 'evtraj'))
    plotevVelDir(behav, bhdata, dateind, 'evTraj');
elseif (strcmp(tagmark, 'sessoccup'))
    plotsesoccup(behav, bhdata, dateind, 'sessOccupancy');
elseif (strcmp(tagmark, 'evoccup'))
    plotevoccup(behav, bhdata, dateind, 'evOccupancy');
elseif (strcmp(tagmark, 'showposdata')) %%%for all other options to plot data
    showalldata(behav, bhdata, dateind, 'allposdata');
end
disp('************************');

function plotevoccup(behav, bhdata, dateind, tmark)
for (i = 1:numel(dateind))
    sess = behav.general.sessID{dateind(i)};
    evname = behav.general.eventname{dateind(i)}; evType = behav.parm.eventType{dateind(i)}; 
    lapX = bhdata.event.LapAllX{dateind(i)};  xj = bhdata.event.Xjoint{dateind(i)};
    laptime = bhdata.event.LapAllPostimestamp{dateind(i)};
    Xbin = bhdata.event.Xbin{dateind(i)}; occup = bhdata.event.Occuptime{dateind(i)};  
    iiev = find(strcmp(evType, 'run')); evname = evname(iiev); lapX = lapX(iiev); 
    xj = xj(iiev); laptime = laptime(iiev); Xbin = Xbin(iiev); occup = occup(iiev);
    if (~isempty(lapX))
       plotlinearnow(Xbin, occup, laptime, xj, sess, evname, 'Linearized X(cm)', tmark);
    else
       disp(['-------------> no run events found: ', sess]);
    end
end

function plotsesoccup(behav, bhdata, dateind, tmark)
ok = 1; evname = 'None'; nsess = numel(dateind);
eventlist = geteventlist(behav, dateind); choice = listdlg('ListString', eventlist);
if choice > 0 
    evname = eventlist{choice};
else 
    ok =0;
end
if ok
    %%%%ask if smoothing
    pp = {'Smooth?'; 'Smoothing sigma (bin)'; 'Peak threshold'}; 
    def = {'Yes'; '2'; '3'};
    III=inputdlg(pp, 'Smoothing parameters', 3, def, 'on');
    ifsmooth = III{1}; sigma = str2num(III{2}); pth = str2num(III{3});
end
if ok
for (j = 1:nsess)
    i = dateind(j);
        sess = behav.general.sessID{i}; sessname = behav.general.sessname{i}; 
        XYbin = bhdata.sess.gridXYbin{i};
        if ~strcmp(evname, 'None') %%%filter through events if selected
            evid = find(strcmp(behav.general.eventname{i}, evname));
            if ~isempty(evid)
               evtimes = bhdata.event.eventtimes{i}{evid}; evstart = evtimes.start; evend = evtimes.ent;
            else
               evstart = 0; evend = 0;
            end
            Pmarker = behav.parm.sessPmarker{i}; allposmarker = behav.general.posMarker{i}; 
            ik = find(strcmp(allposmarker, Pmarker));
            posXX = bhdata.pos.XX{i}{ik}*behav.parm.pixelXSize(i); 
            posYY = bhdata.pos.YY{i}{ik}*behav.parm.pixelYSize(i);
            postimestamp = bhdata.pos.postimestamp{i}*behav.parm.timeunit(i); allposind = [];
            for (tt = 1:numel(evstart))
                 iii = find( (postimestamp>=evstart(tt)) & (postimestamp<=evend(tt)) );
                 allposind = union(allposind, iii);
            end
            posXX = posXX(allposind); posYY = posYY(allposind); framerate = behav.parm.framerate(i);
            occup = OccupancyPoint(XYbin{1}, XYbin{2}, posXX, posYY,0)/framerate;
        else
            occup = bhdata.sess.gridOccuptime{i}; %segoccup = bhdata.sess.gridSegOccuptime{i}{j};
        end
        if (~isempty(occup))
            if strncmpi(ifsmooth, 'Yes', 1) occup = TwoDSmooth_fast(occup, sigma, 5, sigma, 5); end
            aa = sum(sum(occup)); if aa>0 occup = occup/aa; end
            minc = 0; occup = occup * numel(XYbin{1}) * numel(XYbin{2}); %%%normalize to 1/nbins
            %maxc = max(max(occup)); 
            [mm, nn] = size(occup); maxc = prctile(reshape(occup, mm*nn, 1), 99); 
            if maxc==minc maxc = minc+1; end
            [xp,yp] = findmappeaks(occup, XYbin{1}, XYbin{2}, pth);
            hf = figure('Name', strcat(sess, '---', tmark));
            xleng = max(XYbin{1})-min(XYbin{1}); xL = [min(XYbin{1})-xleng*0.1 max(XYbin{1})+xleng*0.1];
            yleng = max(XYbin{2})-min(XYbin{2}); yL = [min(XYbin{2})-yleng*0.1 max(XYbin{2})+yleng*0.1];
            hax = axes('Parent', hf, 'NextPlot', 'add', 'XLim', xL, 'YLim', yL); 
            str = sessname;
            cmap = colormap(jet); %cmap(1,:) = [1 1 1]; %change background to white
            hp = pcolor(XYbin{1}, XYbin{2}, occup); colormap(cmap); colorbar('vert', 'peer', hax);
            for tt = 1:numel(xp)
                line(xp(tt), yp(tt), 'Parent', hax, 'LineStyle', 'none',...
                            'Marker', 'x', 'MarkerSize', 10, 'MarkerEdgeColor', [1 1 1], 'MarkerFaceColor', [1 1 1]);
            end
            xlabel('X (cm)'); ylabel('Y (cm)'); 
            set(hp, 'EdgeColor', 'none', 'LineStyle', 'none', 'Parent', hax); set(hax, 'CLim', [minc maxc]);
            text('Interpreter', 'none', 'Parent', hax, 'String', str, 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.05 0.96]);
        end
end
end
function [peakx, peaky] = findmappeaks(D2rate, Xgrid, Ygrid, pth, N)
%%%define peaks as the ones with values larger than its neighbors
%%%do not count the edges
i = 2:(numel(Ygrid)-1); j = 2:(numel(Xgrid)-1); D2ratenow = D2rate(i,j);
N1 = sign(D2rate(i,j) - D2rate(i-1, j-1)); N2 = sign(D2rate(i,j) - D2rate(i, j-1)); N3 = sign(D2rate(i,j) - D2rate(i+1, j-1));
N4 = sign(D2rate(i,j) - D2rate(i-1, j)); N5 = sign(D2rate(i,j) - D2rate(i+1, j));
N6 = sign(D2rate(i,j) - D2rate(i-1, j+1)); N7 = sign(D2rate(i,j) - D2rate(i, j+1)); N8 = sign(D2rate(i,j) - D2rate(i+1, j+1));
%%%%identify the peaks: how to deal with plateau peaks - inclued first, then excluded repeated ones
allN = N1+N2+N3+N4+N5+N6+N7+N8; allabsN = abs(N1)+abs(N2)+abs(N3)+abs(N4)+abs(N5)+abs(N6)+abs(N7)+abs(N8);
iii = find( (allN>0) & (allN==allabsN) );
peakV = D2ratenow(iii); [peakV, it] = unique(peakV); iii = iii(it);
[peakV, jj] = sort(peakV, 'descend'); iii = iii(jj);
tk = find(peakV>pth); iii = iii(tk);
ny = numel(Ygrid)-2; peakx = Xgrid(j(ceil(iii/ny))); 
t = mod(iii, ny); tt = find(t ==0); t(tt) = ny*ones(size(tt)); peaky = Ygrid(i(t));

function plotsesstraj(behav, bhdata, dateind, tmark)
ok = 1; evname = 'None'; nsess = numel(dateind);
eventlist = geteventlist(behav, dateind); choice = listdlg('ListString', eventlist);
if choice > 0 
    evname = eventlist{choice};
else 
    ok =0;
end
if ok
for (i = 1:nsess)
    j = dateind(i);
    sess = behav.general.sessID{j}; sessname = behav.general.sessname{j};
    posXX = []; posYY = []; %postime = bhdata.pos.postimestamp{i}{j}*behav.parm.timeunit(i); 
    Pmarker = behav.parm.sessPmarker{j}; allposmarker = behav.general.posMarker{j}; 
    ik = find(strcmp(allposmarker, Pmarker)); 
    if (numel(ik) == 1) 
        posXX = bhdata.pos.XX{j}{ik}*behav.parm.pixelXSize(j); 
        posYY = bhdata.pos.YY{j}{ik}*behav.parm.pixelYSize(j);
        if ~strcmp(evname, 'None') %%%filter through events if selected
            evid = find(strcmp(behav.general.eventname{j}, evname)); 
            if ~isempty(evid)
               evtimes = bhdata.event.eventtimes{j}{evid}; evstart = evtimes.start; evend = evtimes.ent;
            else
               evstart = 0; evend = 0;
            end
            postimestamp = bhdata.pos.postimestamp{j}*behav.parm.timeunit(j); allposind = [];
            for (tt = 1:numel(evstart))
                 iii = find( (postimestamp>=evstart(tt)) & (postimestamp<=evend(tt)) );
                 allposind = union(allposind, iii);
            end
            posXX = posXX(allposind); posYY = posYY(allposind);
        end
    end
    if (mean(posXX)>0)
        XYbin = bhdata.sess.gridXYbin{j};
        xleng = max(XYbin{1})-min(XYbin{1}); xL = [min(XYbin{1})-xleng*0.1 max(XYbin{1})+xleng*0.1];
        yleng = max(XYbin{2})-min(XYbin{2}); yL = [min(XYbin{2})-yleng*0.1 max(XYbin{2})+yleng*0.1];
        hf = figure('Name', strcat(sess, '---', tmark)); 
        hax = axes('Parent', hf, 'NextPlot', 'add', 'XLim', xL, 'YLim', yL); 
        str = sessname;
        line(posXX, posYY, 'Parent', hax, 'LineWidth', 1, 'Color', [0 0 0]);
        text('Interpreter', 'none', 'Parent', hax, 'String', str, 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.05 0.96]);
        xlabel('X (cm)'); ylabel('Y (cm)'); 
    end
end
end
function eventlist = geteventlist(behav, dateind)
eventlist{1} = 'None';
for i = 1:numel(dateind)
    j = dateind(i); eventnames = behav.general.eventname{i};
    eventlist = union(eventlist, eventnames);
end
eventlist = unique(eventlist);

function plotevVelDir(behav, bhdata, dateind, tmark)
ok = 1; ev2dname = 'None'; 
eventlist = geteventlist(behav, dateind); choice = listdlg('ListString', eventlist);
if choice > 0 
    ev2dname = eventlist{choice};
else 
    ok =0;
end
if ok
   if strcmp(tmark, 'evVel')
      %%%%ask if smoothing
      pp = {'Smooth speed?'; 'Speed smoothing sigma (s)'}; 
      def = {'Yes'; '0.5'};
      III=inputdlg(pp, 'Parameters for identifying stopping points', 2, def, 'on');
      if ~isempty(III)
          ifsmooth = III{1}; sigma = str2num(III{2});
      else
          ok = 0;
      end
   end
end
for (i = 1:numel(dateind))
    finaldir = behav.general.sessID{dateind(i)}; framerate = behav.parm.framerate(dateind(i));
    evname = behav.general.eventname{dateind(i)}; evType = behav.parm.eventType{dateind(i)}; 
    lapX = bhdata.event.LapAllX{dateind(i)};  xj = bhdata.event.Xjoint{dateind(i)}; laptime = bhdata.event.LapAllPostimestamp{dateind(i)};
    iiev = find(strcmp(evType, 'run')); 
    if ~isempty(iiev)
       evname = evname(iiev); lapX = lapX(iiev); xj = xj(iiev); laptime = laptime(iiev);
       if strcmp(tmark, 'evVel') 
          lapV = bhdata.event.LapAllVel{dateind(i)}(iiev); lapS = [];
          lapV = smoothspeed(lapV, ifsmooth, sigma, framerate);
          if isfield(bhdata.event, 'LapAll1DSpeed')
             lapS = bhdata.event.LapAll1DSpeed{dateind(i)}(iiev);
             lapS = smoothspeed(lapS, ifsmooth, sigma, framerate);
          end
       elseif strcmp(tmark, 'evHeaddir')
          lapV = bhdata.event.LapAllHeadDir{dateind(i)}(iiev); 
       elseif strcmp(tmark, 'evTraj')
          lapV = bhdata.event.LapAllY{dateind(i)}(iiev);
       end
       if (~isempty(lapV))
           plotlinearnow(lapX, lapV, laptime, xj, finaldir, evname, 'Linearized X (cm)', tmark);
           if strcmp(tmark, 'evVel') && (~isempty(lapS)) 
              plotlinearnow(lapX, lapS, laptime, xj, finaldir, evname, 'Linearized X (cm)', 'ev1dSpeed (cm/s)');
           end
       end
    else
        disp(['-------------> no run events found: ', finaldir]);
    end
    %%%%plot 2D speed as a simple line
    allV = bhdata.sess.AllVel{dateind(i)}; allV = smoothsimplespeed(allV, ifsmooth, sigma, framerate);
    j = dateind(i); postimestamp = bhdata.pos.postimestamp{j}*behav.parm.timeunit(j);
    if ~strcmp(ev2dname, 'None') %%%filter through events if selected
            evid = find(strcmp(behav.general.eventname{j}, ev2dname)); 
            if ~isempty(evid)
               evtimes = bhdata.event.eventtimes{j}{evid}; evstart = evtimes.start; evend = evtimes.ent;
            else
               evstart = 0; evend = 0;
            end
             allposind = [];
            for (tt = 1:numel(evstart))
                 iii = find( (postimestamp>=evstart(tt)) & (postimestamp<=evend(tt)) );
                 allposind = union(allposind, iii);
            end
            allposind = sort(allposind);
            postimestamp = postimestamp(allposind); allV = allV(allposind);
    end
    plotsimpleline(postimestamp, allV, finaldir, ev2dname, 'Time (s)', '2d velocity (cm/s)');
end
function lapV = smoothspeed(lapV, ifsmooth, sigma, framerate)%%%lapV{ntraj}{nlap}
if strncmpi(ifsmooth, 'yes', 1)
   sigma = round(sigma*framerate); %%%sigma now if number of frames
   XX = [-5*sigma:1:5*sigma]; wind = normpdf(XX, 0, sigma); wind = wind/sum(wind); nw = numel(wind); %will have a delay = (n-1)/2 bins
   for i = 1:numel(lapV)
       for j = 1:numel(lapV{i})
           vel = lapV{i}{j}; vel = conv(vel, wind); vel = vel( ((nw-1)/2+1):(numel(vel)-(nw-1)/2) );
           lapV{i}{j} = vel; %disp([num2str(i) '; ' num2str(j) '; doing now']);
       end
   end
end
function vel = smoothsimplespeed(vel, ifsmooth, sigma, framerate)%%%lapV()
if strncmpi(ifsmooth, 'yes', 1)
   sigma = round(sigma*framerate); %%%sigma now if number of frames
   XX = [-5*sigma:1:5*sigma]; wind = normpdf(XX, 0, sigma); wind = wind/sum(wind); nw = numel(wind); %will have a delay = (n-1)/2 bins
   vel = conv(vel, wind); vel = vel( ((nw-1)/2+1):(numel(vel)-(nw-1)/2) );
end

function plotevlinearx(behav, bhdata, dateind, tagmark)
for (i = 1:numel(dateind))
    finaldir = behav.general.finaldir{dateind(i)};
    evname = behav.general.eventname{dateind(i)}; evType = behav.parm.eventType{dateind(i)}; 
    laptime = bhdata.event.LapAllPostimestamp{dateind(i)}; %%{ntraj}{nlap}
    xj = bhdata.event.Xjoint{dateind(i)};
    iiev = find(strcmp(evType, 'run')); evname = evname(iiev); 
    if strcmp(tagmark, 'evLinearX') 
        lapX = bhdata.event.LapAllX{dateind(i)}(iiev); 
    elseif strcmp(tagmark, 'evLinearY')
        lapX = bhdata.event.LapAllY{dateind(i)}(iiev);
    end
    laptime = laptime(iiev); xj = xj(iiev); %%plot only those 'run' events
    %lapX = patchupNaN(lapX, laptime, xj); %%%%%interpolate NaN values ---not necessary, already dezeroed
    plotlinearnow(laptime, lapX, laptime, xj, finaldir, evname, 'Time(s)', tagmark);
end

function plotsimpleline(postimestamp, allV, sess, evname, xla, tagmark)
hf = figure('Name', strcat(sess, '---', tagmark)); hax = axes('Parent', hf, 'NextPlot', 'add'); str = [];
line(postimestamp, allV, 'Parent', hax, 'Color', [0 0 0], 'LineWidth', 1); 
str{1} = ['Event: ', evname]; yla = determineylabel(tagmark);
for (i = 1:numel(str))
    text('Interpreter', 'none', 'Parent', hax, 'String', str{i}, 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.05 0.96-(i-1)*0.04]);
end 
xlabel(xla); ylabel(yla);

function plotlinearnow(lapX, lapY, lapT, xj, sess, evname, xla, tagmark)
nev = numel(lapX); col{1} = [0 0 1]; col{2} = [0 1 0]; col{3} = [1 0 0]; col{4} = [0.5 0.5 0.5]; col{5} = [0.5 1 1]; col{6} = [0 0 0];
col{7} = [0.1 0 0]; col{8} = [0 0 0.8]; col{9} = [0.2 1 0]; col{10} = [0.5 0.9 0.5]; col{11} = [0.5 0.2 1]; col{12} = [0.3 0.5 0.1];
ok = 1;
if nev>=3
    [evname, lapX, lapY, lapT, xj, ok] = selectevstoshow(lapX, lapY, lapT, xj, evname);
end
nev = numel(lapX);
if ok
if (strcmp(tagmark, 'evLinearX')) && (nev == 2) % if 2 trajs, plot up and downs
    %if (lapT{1}{1}(1) > lapT{2}{1}(1)) %%%reverse the values of the second trajectory -- this line crashes when first event contains no points
    if min(gettime(lapT{1})) > min(gettime(lapT{2}))    
        mV = xj{1}(numel(xj{1})); xj{1} = mV - xj{1}; for (ki = 1:numel(lapY{1})) lapY{1}{ki} = mV - lapY{1}{ki}; end 
    else
        mV = xj{2}(numel(xj{2})); xj{2} = mV - xj{2}; for (ki = 1:numel(lapY{2})) lapY{2}{ki} = mV - lapY{2}{ki}; end 
    end
elseif (strcmp(tagmark, 'evVel') || strcmp(tagmark, 'evHeaddir') || strcmp(tagmark, 'evTraj'))  && (nev == 2)
    %if (lapT{1}{1}(1) > lapT{2}{1}(1)) %%%reverse the values of the second trajectory -- this line crashes when first event contains no points
    if min(gettime(lapT{1})) > min(gettime(lapT{2}))        
        mV = xj{1}(numel(xj{1})); xj{1} = mV - xj{1}; for (ki = 1:numel(lapX{1})) lapX{1}{ki} = mV - lapX{1}{ki}; end 
    else
        mV = xj{2}(numel(xj{2})); xj{2} = mV - xj{2}; for (ki = 1:numel(lapX{2})) lapX{2}{ki} = mV - lapX{2}{ki}; end 
    end
elseif (strcmp(tagmark, 'evOccupancy')) && (nev == 2)
    %if (lapT{1}{1}(1) > lapT{2}{1}(1)) %%%reverse the values of the second trajectory -- this line crashes when first event contains no points
    if min(gettime(lapT{1})) > min(gettime(lapT{2}))        
        mV = xj{1}(numel(xj{1})); xj{1} = mV - xj{1}; lapX{1} = mV - lapX{1};
    else
        mV = xj{2}(numel(xj{2})); xj{2} = mV - xj{2}; lapX{2} = mV - lapX{2};
    end
end
minY = 0; maxY = 0;
for (i = 1:numel(lapY))
    if (strcmp(tagmark, 'evOccupancy'))
        minY = min([minY min(lapY{i})]); maxY = max([maxY max(lapY{i})]);
    else
        for (j = 1:numel(lapY{i}))  
           if (minY > min(lapY{i}{j})) minY = min(lapY{i}{j}); end
           if (maxY < max(lapY{i}{j})) maxY = max(lapY{i}{j}); end
        end
    end
end
if (nev>0)
    hf = figure('Name', strcat(sess, '---', tagmark)); hax = axes('Parent', hf, 'NextPlot', 'add'); str = [];
    for (i = 1:nev)
        if (strcmp(tagmark, 'evOccupancy'))
             line(lapX{i}, lapY{i}, 'Parent', hax, 'Color', col{i}, 'LineWidth', 2); 
        else
            for (j = 1:numel(lapY{i}))
                 line(lapX{i}{j}, lapY{i}{j}, 'Parent', hax, 'LineStyle', 'none',...
                            'Marker', '.', 'MarkerSize', 5, 'MarkerEdgeColor', col{i}, 'MarkerFaceColor', col{i}); % 'Color', col{i}, 'LineWidth', 2);  
                 if (strcmp(tagmark, 'evLinearX'))
                    set(hax, 'YLim', [-10 max(xj{i})]); %set(hax, 'YLim', [-10 500]);
                    for (k = 1:numel(xj{i}))
                        if ~isempty(lapX{i}{j})
                        line([min(lapX{i}{j}) max(lapX{i}{j})], [xj{i}(k) xj{i}(k)], 'Parent', hax, 'Color', col{i},'LineWidth', 1); 
                        end
                    end
                 end
            end
        end
        yla = determineylabel(tagmark);
        str{i} = ['Event: ', evname{i}];
        if (strcmp(tagmark, 'evVel')) || (strcmp(tagmark, 'evHeaddir')) || (strcmp(tagmark, 'evTraj')) || (strcmp(tagmark, 'evOccupancy')) 
           if (minY < maxY)
           for (k = 1:numel(xj{i}))
                line([xj{i}(k) xj{i}(k)], [minY 1.2*maxY], 'Parent', hax, 'Color', col{i}, 'LineWidth', 1);
           end
           end
        end
    end
    for (i = 1:nev)
       text('Interpreter', 'none', 'Parent', hax, 'String', str{i}, 'Color', col{i}, 'Units', 'normalized', 'Position', [0.05 0.96-(i-1)*0.04]);
    end 
    xlabel(xla); ylabel(yla);
end
end

function yla = determineylabel(tagmark)
yla = tagmark;
if (strcmp(tagmark, 'evOccupancy'))
    yla = strcat(tagmark, '(s)');
elseif (strcmp(tagmark, 'evTraj'))
    yla = 'Linearized Y (cm)';
elseif (strcmp(tagmark, 'evHeaddir'))
    yla = strcat(tagmark, '(o)');
elseif (strcmp(tagmark, 'evVel'))
    yla = strcat(tagmark, '(cm/s)'); 
elseif (strcmp(tagmark, 'evLinearY'))
    yla = strcat(tagmark, '(cm)');
elseif (strcmp(tagmark, 'evLinearX'))
    yla = strcat(tagmark, '(cm)');
end

function lapX = patchupNaN(lapX, laptime, xj) %%%%%interpolate NaN values
for (i = 1:numel(lapX)) %%for each trajectory
    for (j = 1:numel(lapX{i})) %%for each lap
%         if (j == 6)
%             disp([laptime{i}{j} lapX{i}{j}']);
%             disp([]);
%         end
         xnow = lapX{i}{j}; intpos = xnow(1); %xj{i}(1);
         for (k = 2:numel(xnow))
             if (isnan(xnow(k))) || (xnow(k) == 0)
                 xnow(k) = intpos;
             else
                 intpos = xnow(k);
             end
         end
         lapX{i}{j} = xnow;
    end
end

function T = gettime(lapT)
T = [];
for i = 1:numel(lapT)
    if ~isempty(lapT{i})
        T = [T; lapT{i}];
    end
end

function [evname, lapX, lapY, lapT, xj, ok] = selectevstoshow(lapX, lapY, lapT, xj, evname)
[sss, ok] = listdlg('ListString', evname, 'PromptString', 'Choose one or two events to show', 'SelectionMode', 'multiple');
if ok
   if (numel(sss)>2) sss = sss(1:2); end
   evname = evname(sss); lapX = lapX(sss); lapY = lapY(sss); lapT = lapT(sss); xj = xj(sss);
end

function showalldata(behav, bhdata, dateind, tmark)
for (j = 1:numel(dateind))
    k = dateind(j); sess = behav.general.sessID{k}; finaldir = behav.general.finaldir{k}; timeunit = behav.parm.timeunit(k);
    allmarker = behav.general.posMarker{k}; pmarker = behav.parm.sessPmarker{k};
    nm = numel(allmarker); posXX = cell(1, nm); posYY = cell(1, nm); 
    for (ttk = 1:nm)
         posXX{ttk} = bhdata.pos.XX{k}{ttk}*behav.parm.pixelXSize(k); posYY{ttk} = bhdata.pos.YY{k}{ttk}*behav.parm.pixelYSize(k);
    end
    postime = bhdata.pos.postimestamp{k};
if (~isempty(postime))    
    [MCroot, MCname, DAname, DEname] = CurrentVersion;
    hmain = figure('Name', DAname, 'NumberTitle', 'off', 'NextPlot', 'add', 'MenuBar', 'figure', 'Unit', 'normalized',...
       'OuterPosition', [0.204 0.425 0.779 0.576], 'Position', [0.207 0.429 0.773 0.496]);
    hf = uimenu(hmain, 'Label','Add Plot');
    uimenu(hf,'Label','Add New Spike','Callback','DataAnimator_AddSpike_Callback');
    uimenu(hf,'Label','Add New Position','Callback','DataAnimator_AddPosition_Callback');
    uimenu(hf,'Label','Add New EEG','Callback','DataAnimator_AddEEG_Callback');
    uimenu(hf,'Label','Add New Powergram/ratio','Callback','DataAnimator_AddPower_Callback', 'Separator', 'on');
    uimenu(hf,'Label','Add New PSD','Callback','DataAnimator_AddPSD_Callback');
    uimenu(hf,'Label','Add New SleepClass','Callback','DataAnimator_AddSleepClass_Callback');
    ndatafile = 0; displaysetting = []; haxes = {}; data = {};
    setappdata(hmain, 'ndatafile', ndatafile); setappdata(hmain, 'displaysetting', displaysetting);
    setappdata(hmain, 'haxes', haxes); setappdata(hmain, 'data', data);
    ndatafile = ndatafile + 1;
    displaysetting{ndatafile}{1} = ndatafile;  
    for (tt = 1:numel(allmarker))
         displaysetting{ndatafile}{2}{tt} = strcat(finaldir, '__', sess, '__', allmarker{tt});
    end
    displaysetting{ndatafile}{3} = 'p'; %flag for position data type
    displaysetting{ndatafile}{4} = [0 0 1]; %default plot, continuous mode (0), default diode color for plotting track shape (1= red)
    finitpos = 713; %default vaule for where binary data starts in the file
    binarysize = 0;%not used anymore
    ndisplay = 0;   %indicator for how many times file has been read, push scroll bar increase or decrease this number by 1
    ndisplayframe = 1;  %default number of frames being plot at once
    maxdisplay = 10000;
    
    displaysetting{ndatafile}{5} = [ndisplay, ndisplayframe, finitpos, binarysize, 0, 0, 0, 0, maxdisplay]; %more parameters read in PositionDisplay.m
    displaysetting{ndatafile}{6} = [0.5 0.5 0.5]; %default backgrnd
    displaysetting{ndatafile}{7} = [1 0 0]; %default front color
    displaysetting{ndatafile}{8} = [3 0.5]; %default dot line width
    displaysetting{ndatafile}{9} = [0 1 0]; %default back color
    displaysetting{ndatafile}{10}= [0 0 1]; %default line color
    displaysetting{ndatafile}{11}= [-1 -1 -1; -1 -1 -1; -1 -1 -1; -1 -1 -1]; %dumy assignments

    set(0, 'ShowHiddenHandles', 'on');
    H = get(hmain, 'Children');
    TypeH = get(H, 'Type');
    for (i = 1:size(TypeH))
        if (strcmp(TypeH{i}, 'uicontrol')) delete(H(i)); end
    end
    delete(findobj(hmain, 'Tag', 'colorbar'));
    set(0, 'ShowHiddenHandles', 'off');
    
    setappdata(hmain, 'ndatafile', ndatafile);
    setappdata(hmain, 'displaysetting', displaysetting);  
    %%re-allocacte spaces for all the plots
    [axposvector, uiposvector] = DataAnimator_PlotSpaceAuto(ndatafile,0.8);

    %%re-plot all axes and uicontrols
    %DataAnimator_Replot(haxes, ndatafile-1, axposvector, uiposvector);
    %%add new position(raw) plot
    %%[haxes{ndatafile}, hbar] = DataAnimator_PositionDisplay(ndatafile, hmain, axposvector{ndatafile});
    idnum = ndatafile;
    ndisplay = displaysetting{idnum}{5}(1); %start display #
    ndisplayframe = displaysetting{idnum}{5}(2); %how many frame to display at a time
    finitpos = displaysetting{idnum}{5}(3); %starting binary position
    binarysize = displaysetting{idnum}{5}(4); %total binary size in byte

    
    data{idnum}.timestamp = []; 
    for (i = 1:4)
        data{idnum}.x{i} = []; data{idnum}.y{i} = []; 
    end
   nfchoice = 0; nframe = numel(postime); data{idnum}.timestamp = postime;
   for i = 1:numel(displaysetting{idnum}{2})
       fname = displaysetting{idnum}{2}{i};
       data{idnum}.x{i} = posXX{i}; 
       data{idnum}.y{i} = posYY{i};
       if (~isempty(findstr(fname, 'red'))) %if red diode, 
         displaysetting{idnum}{11}(i,:) = [1 0 0]; %displaysetting{idnum}{4}(3) = i; %%%red is the default pos color
       elseif (~isempty(findstr(fname, 'green'))) %if green diode, 
         displaysetting{idnum}{11}(i,:) = [0 1 0]; 
       elseif (~isempty(findstr(fname, 'blue'))) %if blue diode, 
         displaysetting{idnum}{11}(i,:) = [0 0 1]; 
       elseif (~isempty(findstr(fname, 'int'))) %if intensity threshold, 
         displaysetting{idnum}{11}(i,:) = [0.8 0.8 0.8]; 
       end
       if contains(fname, pmarker)
         displaysetting{idnum}{4}(3) = i; %%%This is the pos color for plotting whole track
         %disp(fname)
       end
       XX = []; YY = [];
   end

%disp(strcat('-----> total number of frames:', num2str(nframe)));
%disp(strcat('-----> start time stamp:', num2str(data{idnum}.timestamp(1))));
%disp(strcat('-----> end time stamp:', num2str(data{idnum}.timestamp(nframe))));
minx = min([min(data{idnum}.x{1}) min(data{idnum}.x{2}) min(data{idnum}.x{3}) min(data{idnum}.x{4})]);
maxx = max([max(data{idnum}.x{1}) max(data{idnum}.x{2}) max(data{idnum}.x{3}) max(data{idnum}.x{4})]);
miny = min([min(data{idnum}.y{1}) min(data{idnum}.y{2}) min(data{idnum}.y{3}) min(data{idnum}.y{4})]);
maxy = max([max(data{idnum}.y{1}) max(data{idnum}.y{2}) max(data{idnum}.y{3}) max(data{idnum}.y{4})]);
maxdisplay = ceil(nframe / ndisplayframe);

%%write more setting parameters
displaysetting{idnum}{5}(5) = minx; displaysetting{idnum}{5}(6) = maxx; displaysetting{idnum}{5}(7) = miny;
displaysetting{idnum}{5}(8) = maxy; displaysetting{idnum}{5}(9) = maxdisplay;

hpos = axes('Parent', hmain, 'Units', 'normalized', 'Position', axposvector{idnum}, 'FontSize', 8, 'NextPlot', 'add'); 
xlabel('X'); ylabel('Y'); ylim([0.8*miny 1.2*maxy]); xlim([0.8*minx 1.2*maxx]);

%%% plot whole data to show up the track
g1 = displaysetting{idnum}{4}(3);
       plot(data{idnum}.x{g1}, data{idnum}.y{g1}, 'Parent', hpos, 'LineStyle', 'none', 'Marker', 'o',...
           'MarkerEdgeColor', [0.2 0.2 0.2], 'MarkerFaceColor', [0.2 0.2 0.2], 'MarkerSize', 3); %ground: light gray circle 5 points
%%% plot current position
for (i = 1:numel(displaysetting{idnum}{2}))
    if (~isempty(data{idnum}.x{i}))
        plotcolornow = displaysetting{idnum}{11}(i,:);
        plot(data{idnum}.x{i}(ndisplay+1), data{idnum}.y{i}(ndisplay+1), 'Parent', hpos, 'LineStyle', 'none', 'Marker', 'o',...
           'MarkerEdgeColor', plotcolornow, 'MarkerFaceColor', plotcolornow, 'MarkerSize', 12);     %front: red circle 3 points
    end
end
infotext = ['Current time(s): ',num2str(timeunit*data{idnum}.timestamp(ndisplay+1), '%10.5f')];
%infotext = strcat('Timestamps:',num2str(data{idnum}.timestamp(ndisplay+1)));
text('String', infotext, 'Interpreter', 'none', 'Parent', hpos, 'Tag', 'info',...
                      'Position', [0.01 0.95], 'Units', 'normalized', 'Color', [0 0 0], 'FontSize', 8); % add buffer and data info
for (i = 1:numel(displaysetting{idnum}{2}))
    fname = displaysetting{idnum}{2}{i}; posnow = 0.9 - (i-1)*0.05; plotcolornow = displaysetting{idnum}{11}(i,:);
    titletext = strcat(num2str(idnum), ' = ', fname,' -Position');
    text('String', titletext, 'Interpreter', 'none', 'Parent', hpos, 'Tag', 'title',...
              'Position', [0.01 posnow], 'Units', 'normalized', 'Color', plotcolornow, 'FontSize', 8); % add a title
end
%%plot a control bar with the axes (%%%main figure (hmain))
barpos = [axposvector{idnum}(1) axposvector{idnum}(2)-0.10*axposvector{idnum}(4)/0.80 axposvector{idnum}(3) axposvector{idnum}(4)*0.06];
hbar = uicontrol('Style', 'slider', 'Parent', hmain, 'Units', 'normalized', 'Position', barpos,...
                 'Min', 0, 'Max', maxdisplay-1, 'SliderStep', [1/maxdisplay, 10/maxdisplay], 'Callback', 'DataAnimator_PositionBar_Callback');
set(hbar, 'Value', ndisplay); %%initial bar position
setappdata(hbar, 'idnum', idnum);
setappdata(hmain, 'displaysetting', displaysetting);
setappdata(hmain, 'data', data); haxes{ndatafile} = hpos;
    %%add new control panel
    DataAnimator_PositionControl(hmain, haxes{ndatafile}, ndatafile, uiposvector{idnum});
    setappdata(hmain, 'haxes', haxes);

end
end