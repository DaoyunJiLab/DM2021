function DM_Seq_ShowSequences_Callback
%%%display sequencing results of a .seqdb

hf = gcbf; pinfonow = getappdata(hf, 'pinfo'); datanow = getappdata(hf, 'data'); tagnow = get(gcbo, 'Tag');
plotparm = getappdata(hf, 'plotparm');
hspike = getappdata(hf, 'hspike'); spikeselection = getappdata(hspike, 'selection'); cellind = find(spikeselection==1);
%hgroup = getappdata(hf, 'hgroup'); groupselection = getappdata(hgroup, 'selection'); cellind = []; grpind = find(groupselection == 1); 
%for (kk = 1:numel(grpind)) cellind = union(cellind, datanow.grouplist.groupindex{grpind(kk)}); end
   ok = 1;
disp('-----> Displaying decoding/sequencing results ......');
if (plotparm.evselect == 0) %%%select sessions
      disp('-----------> Sequences cannot be displayed: need to select the event option'); ok = 0;
elseif ~isfield(datanow, 'data')
      ok = 0; disp('-----------> aborted: data not available for stripped databases');
end
if ok
   if (numel(cellind)>1) disp('-----------> more than 1 templates selected; only display the first selection'); end
   cellind = cellind(1); %%%%for sequence disiplay: only display the  first selected template
   if (plotparm.evselect == 0) %%%select sessions
       seqdata = datanow.data.sessseq{cellind}; 
       evnames = pinfonow.general.sessionname{cellind}; infstr = 'sess';
       eventoption = pinfonow.parm.sessionoption{cellind}; 
       timingoption = pinfonow.parm.sessiontimingoption{cellind};
   else
    seqdata = []; someparm.ifshuf = 0;
    if strcmp(tagnow, 'showsequences')%%%select either actual copy or shuffled copy
       seqdata = datanow.data.evseq{cellind}; 
    elseif strcmp(tagnow, 'showshufflecopy')
        if isfield(datanow.data, 'shufevseq')
           seqdata = findshuffledatacopy(datanow.data.shufevseq{cellind}, pinfonow.parm, cellind); someparm.ifshuf = 1;
        else
           disp('-----------> shuffled data cannot be dispalyed: not exit or saved');
        end
    end
    someparm.pthre = []; someparm.sthre = []; someparm.maxstepdis = []; someparm.minstepdis = [];
    someparm.minmedpprob = []; someparm.maxNstep = []; someparm.minNstep = []; someparm.maxleng = []; someparm.minleng = [];
    if isfield(pinfonow.parm, 'probthres') someparm.pthre = pinfonow.parm.probthres{cellind}; end
    if isfield(pinfonow.parm, 'speedthres') someparm.sthre = pinfonow.parm.speedthres{cellind}; end
    if isfield(pinfonow.parm, 'maxstepdis') someparm.maxstepdis = pinfonow.parm.maxstepdis{cellind};someparm.minstepdis = pinfonow.parm.minstepdis{cellind}; end
    if isfield(pinfonow.parm, 'minmedpprob') someparm.minmedpprob = pinfonow.parm.minmedpprob{cellind}; end
    if isfield(pinfonow.parm, 'maxNstep') someparm.maxNstep = pinfonow.parm.maxNstep{cellind};someparm.minNstep = pinfonow.parm.minNstep{cellind}; end
    if isfield(pinfonow.parm, 'maxleng') someparm.maxleng = pinfonow.parm.maxleng{cellind};someparm.minleng = pinfonow.parm.minleng{cellind}; end
    tmpID = pinfonow.general.tmpID{cellind}; %tmpType = pinfonow.parm.tmpType{cellind}; 
    animaldate = strcat(pinfonow.general.animalname{cellind}, '_', pinfonow.general.datedir{cellind});
    evnames = pinfonow.general.eventname{cellind}; %evType = pinfonow.parm.eventType{cellind}; 
    infstr = 'evt';
    someparm.eventoption = pinfonow.parm.eventoption{cellind}; 
    someparm.timingoption = [];
    if isfield(pinfonow.parm, 'eventtimingoption')
       someparm.timingoption = pinfonow.parm.eventtimingoption{cellind};
    end
   end
end
if ok
   rrmode = []; if isfield(pinfonow.parm, 'rankmode') rrmode = pinfonow.parm.rankmode{cellind}; end 
   parentfilename = []; if isfield(datanow, 'parentfile') parentfilename = datanow.parentfile; end
   tmpevname = pinfonow.tmp.sessevname{cellind}; tmpfilename = pinfonow.tmp.filename{cellind};
   rrank = []; rfile = [];
   [MCroot, ~, ~, ~, ~, ~] = CurrentVersion;
   if strcmp(rrmode, 'PairRank')
     rfile = fullfile(MCroot, 'DataExplorer', 'Sequence', 'rrankprob.mat'); 
   elseif strcmp(rrmode, 'OrderRank')
     rfile = fullfile(MCroot, 'DataExplorer', 'Sequence', 'orderrankprob.mat');
   end
   if (exist(rfile, 'file') == 2)
      S = load(rfile); rrank = S.rrank; S = [];
   elseif strcmp(rrmode, 'PairRank')
      ok = 0; disp('-----------> pair ranking probability file for selected rank mode not found; aborted');
   end
end
if ok
    seqtype = [];
    if (isfield(pinfonow.parm, 'seqtype')) someparm.seqtype = pinfonow.parm.seqtype{cellind}; end
    if contains(someparm.seqtype, 'evtitemized') 
        nev = numel(tmpID); evnames = erase(evnames, [animaldate '_']); 
        title = strcat('Event-', evnames, '___Tmp-', tmpID);
    else
        nev = numel(evnames); title = strcat('Tmp-', tmpID, '___Event-', evnames);
    end
end
if ok
  for (i = 1:nev)
   %disp(seqdata)
   seqnow = seqdata{i}; 
   if (isempty(someparm.seqtype)) || (contains(someparm.seqtype, 'Bayesian'))  %%%if bayesian decoding
       evT = datanow.events.eventimes{cellind}{i};
       if ~(someparm.ifshuf)
           if contains(someparm.seqtype, '2D')
              showevind = datanow.events.targevind{cellind}{i};
           else
              showevind = datanow.events.matchevind{cellind}{i};
           end
       else
           showevind = [];
       end
       if contains(someparm.seqtype, 'evtitemized')
          plotcurrentdecodingevents(hf, seqnow, evnames, showevind, evT, cellind, pinfonow, datanow, title{i}, tmpevname{i}, tmpfilename{i}, someparm);
       else %if template itemized
          plotcurrentdecodingevents(hf, seqnow, evnames{i}, showevind, evT, cellind, pinfonow, datanow, title{i}, tmpevname, tmpfilename, someparm);
       end
   else %if a regular sequence matching method
   [mm, nn] = size(seqnow.ep);
   seq = seqnow.seq; seqstart = seqnow.seqstart; seqend = seqnow.seqend; seqmarker = seqnow.seqmarker;
   matchscore = seqnow.matchscore; posmatchprob = seqnow.posmatchprob; negmatchprob = seqnow.negmatchprob;
   %%%wrap the results into a standard database form and save results
   SSS = []; StT = []; EdT = []; MMM = []; SCr = []; posP = []; negP = [];
   grpname = cell(mm*nn,1); grpind = cell(mm*nn,1); grpcrit = cell(mm*nn,1); grptype = cell(mm*nn,1); 
   for (ii = 1:mm) %%%event index
       for (jj = 1:nn) %%%windowsize
           indstart = numel(SSS)+1;
           SSS = [SSS; seq{ii,jj}]; 
           SCr = [SCr; matchscore{ii,jj}]; posP =[posP; posmatchprob{ii,jj}]; negP =[negP; negmatchprob{ii,jj}];
           StT = [StT; seqstart{ii,jj}]; EdT = [EdT; seqend{ii,jj}]; MMM = [MMM; seqmarker{ii,jj}];
           indend = numel(SSS);
           if ~strcmp(eventoption, 'FreeSequences')
              grpname{(ii-1)*nn + jj} = ['ev',num2str(ii),'_ws',num2str(jj)]; 
           else
              grpname{(ii-1)*nn + jj} = ['ev',num2str(ii)]; 
           end
           grpind{(ii-1)*nn + jj} = indstart:indend; 
           grpcrit{(ii-1)*nn + jj}{1} = 'Manual'; grptype{(ii-1)*nn + jj} = 'program';
       end
   end
   
   pinfo = []; data = []; pinfo.work = struct([]); 
   if isfield(seqnow, 'seqdata')
      data.seq.matchspikedata = seqnow.seqdata; 
   elseif isfield(seqnow, 'cellpeaktime')
      data.seq.matchspikedata = seqnow.cellpeaktime;
   else
      data.seq.matchspikedata = [];
   end
   if isfield(seqnow, 'timepoint')
      data.seq.timepoint = seqnow.timepoint; 
   else
      data.seq.timepoint = [];
   end
   if isfield(seqnow, 'allletter')
      data.seq.allletter = seqnow.allletter; 
   else
      data.seq.allletter = [];
   end
   if isfield(seqnow, 'alltime')
      data.seq.alltime = seqnow.alltime; 
   else
      data.seq.alltime = [];
   end
   [mm,nn] = size(seqnow.ep);
   if (mm>1) && (nn == 1)
      data.seq.ep = seqnow.ep';
   else
      data.seq.ep = seqnow.ep;
   end
   data.seq.matchletter = cell2mat(pinfonow.tmp.fileletter{cellind});
   data.seq.timingoption = timingoption; data.seq.tmprank = pinfonow.tmp.tmprank{cellind};
   data.seq.rrank = rrank; data.seq.parentfilename = parentfilename;
   
   %data.seq.peakth = PeakTh;
   %data.seq.laptime = laptime; data.seq.lappos = lappos; data.seq.trajvert = pout;
   pinfo.seq.seqID = []; pinfo.seq.sequence = []; pinfo.seq.matchscore = []; pinfo.seq.posmatchprob = []; pinfo.seq.negmatchprob = [];
   pinfo.seq.seqstart = []; pinfo.seq.seqend = []; pinfo.seq.seqmarker = [];   
   for (j = 1:numel(SSS))
       pinfo.seq.seqID{j} = strcat(num2str(j), '_', SSS{j});
       pinfo.seq.sequence{j} = SSS{j}; 
       pinfo.seq.matchscore{j} = SCr(j); pinfo.seq.posmatchprob{j} = posP(j); pinfo.seq.negmatchprob{j} = negP(j);
       pinfo.seq.seqstart{j} = StT(j); pinfo.seq.seqend{j} = EdT(j); pinfo.seq.seqmarker{j} = MMM{j};
   end
   
   if isempty(pinfo.seq.sequence)
       disp(['---------> warning: no sequences:  ', title]);
   else
   
   grpname = ['List0'; grpname]; grpind = [1:numel(SSS); grpind]; grpcrit = ['all'; grpcrit]; grptype = ['Manual'; grptype];
   data.grouplist.groupname = grpname; data.grouplist.groupindex = grpind;
   data.grouplist.grouptype = grptype; data.grouplist.groupcrit = grpcrit;
   
   pinfo.general.tmpfile{1} = pinfonow.tmp.tmpfile{cellind}; 
   pinfo.general.epfile{1} = evnames{i}; 
   pinfo.general.sessionflag{1} = infstr; 
   pinfo.general.animaldate{1} = strcat(pinfonow.general.animalname{cellind}, '_', pinfonow.general.datedir{cellind}); 
   pinfo.general.cellname{1} = pinfonow.tmp.filename{cellind}; 
   pinfo.general.cellletter{1} = pinfonow.tmp.fileletter{cellind};
   %for (i = 1:numel(matchletter)) pinfo.general.cellletter{1}{i} = char(matchletter(i)); end 
   pinfo.general.tmpseq{1} = pinfonow.tmp.tmpseq{cellind}; pinfo.general.matchtmpseq{1} = pinfonow.tmp.tmpseq{cellind};
   %pinfo.general.matchtmpseq{1} = seltmp;
   
   %pinfo.general.cellletter{1} = matchletter';
   if (plotparm.evselect == 0) %%if sessions
      pinfo.parm.timingoption{1} = pinfonow.parm.sessiontimingoption{cellind}; 
   else
      pinfo.parm.timingoption{1} = pinfonow.parm.eventtimingoption{cellind};
   end
   pinfo.parm.eventoption{1} = pinfonow.parm.eventoption{cellind}; 
   %pinfo.parm.minwindowsize = minwindowsize; pinfo.parm.maxwindowsize = maxwindowsize; 
   %pinfo.parm.windowsizestep = windowsizestep; pinfo.parm.windowshifttime = windowshifttime;
   %pinfo.parm.resolution = res; 
   pinfo.parm.sigma(1) = pinfonow.parm.timesigma{cellind};
   %pinfo.parm.shuffleoption = shuffleoption; pinfo.parm.shuffletime = shuffletime; pinfo.parm.shufflemode = shufflemode; 
   %save(seqfile, 'pinfo', 'data', '-mat'); data.seq.selffilename = seqfile;
   [MCroot, MCname, DAname, DEname, ANname, MAname] = CurrentVersion;
   hmain = figure('Name', MAname, 'NumberTitle', 'off', 'NextPlot', 'add',...
    'MenuBar', 'figure', 'Units', 'normalized', 'Position', [0.05 0.2 0.9 0.7]);
   DE_Seq_DisplaySequences(hmain, pinfo, data); 
   namenow = strcat(get(hmain, 'Name'), '__', title);
   set(hmain, 'Name', namenow);
   setappdata(hmain, 'pinfo', pinfo); setappdata(hmain, 'data', data);
   end
   end
  end
end
disp('*******************************');

function plotcurrentdecodingevents(hf, seqnow, evName, showevind, evT, cellid, pinfonow, datanow, selffilename, tmpevname, tmpfilename, someparm)
%%%work out the behavioral data
ok = 1;
plotparm = getappdata(hf, 'plotparm');
if (plotparm.linkbehav == 0)
       disp(['-----------> no behav data linked']); ok = 0;
else
       behav = getappdata(hf, 'behav'); bhdata = getappdata(hf, 'bhdata');
end
if ok
   finaldirnow = pinfonow.general.finaldir{cellid};
   %%%idenity animal's 2D and 1D (if possible) positions
   posid = []; evid = []; pout = []; joint = [];
   postimestamp = []; fXX = []; fYY = []; bXX = []; bYY = []; posXX = []; posYY = []; XX = []; YY = []; 
   [evsessname, ~] = identifysession(evT, pinfonow.general.sessionname{cellid}, pinfonow.general.sessionstartT{cellid}, pinfonow.general.sessionendT{cellid});
   posid = find( strcmp(behav.general.finaldir, finaldirnow) & strcmp(behav.general.sessname, evsessname) );
   if numel(posid) == 1
      %%%%2D position data 
      framerate = behav.parm.framerate(posid); 
      postimestamp = bhdata.pos.postimestamp{posid}; 
      stime = behav.general.sessstartT{posid}; etime = behav.general.sessendT{posid}; 
      Pmarker = behav.parm.sessPmarker{posid}; Fmarker = behav.parm.sessFmarker{posid}; Bmarker = behav.parm.sessBmarker{posid}; 
      allposmarker = behav.general.posMarker{posid}; 
      posXX = []; posYY = []; ik = find(strcmp(allposmarker, Pmarker)); 
      if (numel(ik) == 1) posXX = bhdata.pos.XX{posid}{ik}*behav.parm.pixelXSize(posid); posYY = bhdata.pos.YY{posid}{ik}*behav.parm.pixelYSize(posid); end
      fXX = []; fYY = []; ik = find(strcmp(allposmarker, Fmarker)); 
      if (numel(ik) == 1) fXX = bhdata.pos.XX{posid}{ik}*behav.parm.pixelXSize(posid); fYY = bhdata.pos.YY{posid}{ik}*behav.parm.pixelYSize(posid); end
      bXX = []; bYY = []; ik = find(strcmp(allposmarker, Bmarker)); 
      if (numel(ik) == 1) bXX = bhdata.pos.XX{posid}{ik}*behav.parm.pixelXSize(posid); bYY = bhdata.pos.YY{posid}{ik}*behav.parm.pixelYSize(posid); end
      XX = bhdata.pos.XX{posid}; YY = bhdata.pos.YY{posid};
      for (ttkk = 1:numel(XX))
          XX{ttkk} = XX{ttkk}*behav.parm.pixelXSize(posid); YY{ttkk} = YY{ttkk}*behav.parm.pixelYSize(posid); 
      end
      
      if (isfield(behav.general, 'eventname'))
          evid = find(strcmp(behav.general.eventname{posid}, evName));
      else
          evid = find(strcmp(behav.behavior.eventname{posid}, evName));
      end 
   else 
      disp(['-------------> position data for the event not found: ', finaldirnow, '___', evName]);
   end
   if (numel(evid)==1)
       realtime = bhdata.event.LapAllPostimestamp{posid}{evid}; 
       reallocx = cell(size(realtime));
       %if  numel(seqnow.Xgrid) == 1 %%%if id decoding: linearize according to the template trajectory
                                    %%%for 2d decoding, 1d linearization is done at the time of display
           if strcmp(evName, tmpevname)%%%if current event is the same as template event
              pout = bhdata.event.Xjoint{posid}{evid}; joint = bhdata.pos.posltr{posid}{evid};
              if numel(seqnow.Xgrid) == 1 reallocx = bhdata.event.LapAllX{posid}{evid}; end  
           else   %if current event not the same as template event: linearize the current postion on the template trajectory(why?) 
              [joint, spaceunit] = locatejoint(behav, bhdata, tmpevname, tmpfilename);
              if ~isempty(joint)
                 joint(:,1) = joint(:,1)*behav.parm.pixelXSize(posid); joint(:,2) = joint(:,2)*behav.parm.pixelYSize(posid);
                 if (~isempty(joint)) && (numel(seqnow.Xgrid) == 1)
                    [realtime, reallocx, pout] = linearizenow(joint, realtime, postimestamp*behav.parm.timeunit(posid), posXX, posYY, spaceunit);
                 end
              end
           end
       %end
   else
       disp(['-------------> 1D positions not computed: ', finaldirnow, '___', evName]); 
   end
   %%%%%%save results as a session file
   if ok 
      pinfo = []; pinfo.sesstype = 'decodedpos'; pinfo.plottype = 'colorplots';
      pinfo.data = []; pinfo.data{1} = []; pinfo.plotid = 1;
      [mm,nn] = size(seqnow.timepoint); timepoint = []; prob = []; seqstart = []; seqend = []; seqtime = []; seqloc = cell(1,2);
      for (ii = 1:mm)
       for (jj = 1:nn)
           %disp([size(seqnow.timepoint) size(seqnow.seqdata)]);
           %disp([size(seqnow.timepoint{ii,jj}) size(seqnow.seqdata{ii,jj})]);
           timepoint = [timepoint seqnow.timepoint{ii,jj}{1}]; prob = [prob seqnow.seqdata{ii,jj}{1}];
           seqstart = [seqstart seqnow.seqstart{ii,jj}]; seqend = [seqend seqnow.seqend{ii,jj}];
           if ~contains(someparm.seqtype, '2D') %%%for 1D decoding
              for (k = 1:numel(seqnow.quantdata{ii,jj}.timepoint))
                  seqtime = [seqtime seqnow.quantdata{ii,jj}.timepoint{k}]; seqloc{1} = [seqloc{1} seqnow.quantdata{ii,jj}.peakloc{k}];
              end
           else
              %for (k = 1:numel(seqnow.timepoint{ii,jj}{1}))
                  seqloc{1} = [seqloc{1} seqnow.quantdata{ii,jj}.peakx{1}]; seqloc{2} = [seqloc{2} seqnow.quantdata{ii,jj}.peaky{1}];
              %end
           end
              
       end
      end
      if ~contains(someparm.eventoption, 'FreeSeq') %%%% if not Free sequencing
         seqstart = seqstart(showevind); seqend = seqend(showevind);
      end
      if contains(someparm.seqtype, '2D') seqtime = timepoint; end
      alltime = []; allposx = []; allposy = [];
      for (i = 1:numel(realtime)) alltime = [alltime; realtime{i}]; end
      if numel(seqnow.Xgrid) == 1 %%%if id decoding
         for (i = 1:numel(realtime))
              allposx = [allposx reallocx{i}];
         end
      elseif numel(seqnow.Xgrid) == 2 %%%if 2d decoding
         alltime = postimestamp*behav.parm.timeunit(posid); allposx = posXX; allposy = posYY;
      end   
      pinfo.parm.sessionstart = stime; pinfo.parm.sessionend = etime;
      pinfo.data{1}.timepoint = timepoint; pinfo.data{1}.prob = prob; pinfo.data{1}.Xgrid = seqnow.Xgrid;    
      pinfo.data{1}.alltime = alltime; pinfo.data{1}.allposx= allposx; pinfo.data{1}.allposy = allposy; 
      pinfo.data{1}.trajvert = pout; pinfo.data{1}.joint = joint;
      if ~isempty(someparm.pthre) pinfo.data{1}.peakthres = someparm.pthre; end
      if ~isempty(someparm.sthre) pinfo.data{1}.speedthres = someparm.sthre; end
      if ~isempty(someparm.maxstepdis) pinfo.data{1}.maxstepdis = someparm.maxstepdis; pinfo.data{1}.minstepdis = someparm.minstepdis; end
      if ~isempty(someparm.minmedpprob) pinfo.data{1}.minmedpprob = someparm.minmedpprob; end
      if ~isempty(someparm.maxNstep) pinfo.data{1}.maxNstep = someparm.maxNstep; pinfo.data{1}.minNstep = someparm.minNstep; end
      if ~isempty(someparm.maxleng) pinfo.data{1}.maxleng = someparm.maxleng; pinfo.data{1}.minleng = someparm.minleng; end      
      pinfo.selffilename = []; pinfo.seqstart = seqstart; pinfo.seqend = seqend; pinfo.seqtime = seqtime; pinfo.seqloc = seqloc; %disp(['---> # of sequences', numel(seqstart)]);
      %%%position data
      pinfo.posfile = 'Multiple (unspecified)'; pinfo.posdata.postimestamp = postimestamp; pinfo.posdata.xfront = fXX;
      pinfo.posdata.yfront = fYY; pinfo.posdata.xback = bXX; pinfo.posdata.yback = bYY; pinfo.posdata.xpos = posXX;
      pinfo.posdata.ypos = posYY; pinfo.posdata.XX = XX; pinfo.posdata.YY = YY; pinfo.diodepos = [];
      %%%%need to plugin more session (dumy) parameters
      pinfo = assigndumyvariables(pinfo, seqnow, pinfonow, datanow, evName, evT, cellid); 
      %%%display results as a rate session file
      hmain = DataExplorer_DataExplorer_Callback;
      %[MCroot, MCname, DAname, DEname, ANname, MAname] = CurrentVersion;
      fname = strcat(get(hmain,'Name'), '__', selffilename); if someparm.ifshuf fname = strcat(fname, '__Shuffle'); end; set(hmain, 'Name', fname); 
      %hmain = figure('Name', fname, 'NextPlot', 'add', 'NumberTitle', 'off', 'Color', [0 0 0], 'Unit', 'normalized',...
      %      'OuterPosition', [0.132 0.248 0.848 0.697], 'Position', [0.135 0.252 0.842 0.617]);
      setappdata(hmain, 'pinfo', pinfo);
      pinfo = DE_PlotRasterAndAll(hmain, pinfo);
      setappdata(hmain, 'pinfo', pinfo); 
   end
end

function pinfo = assigndumyvariables(pinfo, seqnow, pinfonow, datanow, evName, ev, cellid)
%%%%%%%%%%%%%%%%%%%%%%%here is the big variable assignment
pinfo.parm.bin = 0.1; pinfo.parm.maxlag = 10;
pinfo.parm.grid1d = 10; pinfo.parm.grid2d = 10; pinfo.parm.frametime = 0.033; pinfo.parm.wavepoint = 32;
pinfo.parm.nepdisplay = 1; pinfo.parm.plotnumber = 1; pinfo.parm.timeunit = [1 0.0001 0.0001]; pinfo.parm.plotduration = 10;
pinfo.parm.sessionflag = ev.marker; pinfo.parm.setbacktime = 2;
pinfo.parm.searchtime = 0.5; pinfo.parm.align = 'start';
pinfo.parm.buffersize = 512; pinfo.parm.freq = 2002; pinfo.parm.filegain = []; freq = 2002;
pinfo.parm.windowpoint = round(2*freq); pinfo.parm.shiftpoint = round(1*freq);
leftcontrolwidth = 0.05; rightcontrolwidth = 0.2; plotwidth = 1-rightcontrolwidth-leftcontrolwidth;
scrollbarheight = 0.05; addplotheight = 0.30; spikeplotheight = 1-addplotheight-scrollbarheight;
pinfo.parm.leftcontrolwidth = leftcontrolwidth; pinfo.parm.rightcontrolwidth = rightcontrolwidth;
pinfo.parm.scrollbarheight = scrollbarheight; pinfo.parm.addplotheight = addplotheight;
pinfo.filename = pinfonow.general.tmpID{1}; pinfo.filearea{1} = pinfonow.general.recarea{1}; 
pinfo.datatype{1} = 'spike';  pinfo.filegain = [];
ndata = numel(pinfo.data); id = 1:ndata; currentid = 1; plotid = id;
pinfo.id = id; pinfo.plotid = plotid; 
plotselect = zeros(1, ndata); plotselect(1) = 1; pinfo.plotselect = plotselect; %default all plot deselected = 0 except for 1st one (0=not selected)            
pinfo.epfile = evName; pinfo.ep = ev; pinfo.animaldate = strcat(pinfonow.general.animalname{cellid}, '_', pinfonow.general.datedir{cellid}); 
pinfo.epselect = zeros(1, numel(ev.start));
pinfo.jicolor = []; pinfo.plotcolor = [];

function [evSess, evsessid] = identifysession(evTime, sessionname, startT, endT)
evSess = []; evsessid = []; minevstart = min(evTime.start); maxevend = max(evTime.ent);
iii = find( (startT <= minevstart) & (endT >= maxevend) );
if (numel(iii) == 1)
    evSess = sessionname{iii}; evsessid = iii;
end
  
function [joint, spaceunit] = locatejoint(behav, bhdata, tmpevname, tmpfilename)
joint = [];
%%%get the final directory
tind = strfind(tmpfilename{1}, 'final'); tmpfdir = tmpfilename{1}(1:tind+4);
sessid = []; evid = []; nsess = numel(behav.general.finaldir);
for (i = 1:nsess)
    evnamenow = behav.general.eventname{i}; finaldirnow = behav.general.finaldir{i};
    jjj = find(strcmp(evnamenow, tmpevname));
    if (numel(jjj)==1) && strcmp(finaldirnow, tmpfdir)
        sessid = i; evid = jjj; break
    end
end
if ~isempty(sessid)
    posltrfile = bhdata.pos.ltrfilename{sessid}; posltr = bhdata.pos.posltr{sessid};
    evPosltr = behav.parm.eventPosltr{sessid}; spaceunit = (behav.parm.pixelXSize(sessid) + behav.parm.pixelYSize(sessid))/2;
    for (tt = 1:numel(posltrfile))
        if (strcmp(posltrfile{tt}, evPosltr{evid}))
           joint = posltr{tt}; break
        end
    end
end

function [realtime, reallocx, pout] = linearizenow(joint, realtime, postimestamp, xpos, ypos, spaceunit)
reallocx = cell(size(realtime)); pout = [];
for (tt = 1:numel(realtime))
    if ~isempty(realtime{tt})
       evposind = find( (postimestamp>=min(realtime{tt})) & (postimestamp<=max(realtime{tt})) );
       realtime{tt} = postimestamp(evposind); reallocx{tt} = NaN*ones(size(realtime{tt}));
       if ~isempty(evposind)
          [xout, yout, pout] = LinearizePathXY_new(joint, xpos(evposind), ypos(evposind), 0.9, 0, spaceunit);
          reallocx{tt} = xout; 
       end
    end
end

function seqdata = findshuffledatacopy(shufevseq, parm, cellind)
%%%%shufevseq{nushuffle, nev}, but not all shuffled copies are saved - randomly choose one every time
seqdata = []; savecopynum = []; 
if isfield(parm, 'Nshuffle') && isfield(parm, 'saveshufrawdatarate')
   nshuffle = parm.Nshuffle{cellind}; saverate = parm.saveshufrawdatarate{cellind};
   savecopynum = 1:saverate:nshuffle;
end
if isempty(savecopynum)
    disp('-----------> shuffled data cannot be dispalyed: not exist or saved');
else
    kk = randperm(numel(savecopynum)); copynumnow = savecopynum(kk(1));
    seqdata = shufevseq(copynumnow, :); %this way () means output remains a cell array of {1, nev}    
end

