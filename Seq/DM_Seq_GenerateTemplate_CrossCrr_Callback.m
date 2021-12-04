function DM_Seq_GenerateTemplate_CrossCrr_Callback
%%Generate template based on (lap-averaged if linear) cross-crr peak times (on linear tracks or open) on a crrdb 
%%template defined as a sequence
%%output as a .tmp file 

hf = gcbf; hgroup = getappdata(hf, 'hgroup'); ok = 1;
plotparm = getappdata(hf, 'plotparm'); vv = plotparm.showdetail; %if vv=1, shoe details of computation
cc = plotparm.compute; %if cc=0, assign parameters for computation only; if 1, do the real computation
ow = plotparm.overwrite; %if ow=1, plot into the current database
fname = get(hf, 'Name'); ftt = strfind(fname, '__'); currentfilename = fname(ftt+2:numel(fname));
[pp, nn, ee] = fileparts(currentfilename);
if ~strcmp(ee, '.crrdb')
    disp('-----> not a crrdb; aborted'); ok = 0;
else
    pinfo = getappdata(hf, 'pinfo'); data = getappdata(hf, 'data');
    hgroup = getappdata(hf, 'hgroup'); groupselection = getappdata(hgroup, 'selection'); cellind = []; grpind = find(groupselection == 1); 
    for (kk = 1:numel(grpind)) cellind = union(cellind, data.grouplist.groupindex{grpind(kk)}); end
    se = {'sessions'; 'events'};
    [ttt, ok] = listdlg('ListString', se, 'PromptString', 'Which time periods');
    if ok 
        sessev = se{ttt}; crrmode = 'crr';
        if strcmp(sessev, 'events') 
           cat = {'crr'; 'lapcrr'; 'lapspatialcrr'};
           [sss, ok] = listdlg('ListString', cat, 'PromptString', 'Which type of crr?');
           if ok
              crrmode =  cat{sss}; 
           end
        end
        if ~isfield(pinfo, crrmode) ok = 0; end
    end
end
if ok
    ccparm.pzsc = []; ccconsis = []; ccparm.kurt = []; ccparm.ptime = []; ccparm.consis = [];
    if ~strcmp(crrmode, 'crr') %%for lapcrr or spatial lapcrr
       input = inputdlg({'Enter min crr stability'; 'Enter min crr peak zscore'; 'Enter min crr kurtosis'; 'Enter max crr peak time/loc (s/pixels)'},...
          'Pair filtering thresholds', 4, {'0.2'; '1'; '-5'; '2.5'}); 
       if (~isempty(input))
          ccparm.consis = str2num(input{1}); ccparm.pzsc = str2num(input{2}); ccparm.kurt = str2num(input{3}); ccparm.ptime = str2num(input{4}); 
       else
          ok = 0;
       end
    else %if strcmp(crrmode, 'crr')
       input = inputdlg({'Enter min peak crr threshold (Z)';'Enter max crr peaktime threshold';'Enter min crr kurtosis threshold'},...
           'Crr selection Parameters', 3, {'2.5'; '1.501'; '-5'}); 
       if (~isempty(input))
          ccparm.pzsc = str2num(input{1}); ccparm.ptime = str2num(input{2}); ccparm.kurt = str2num(input{3});
       else
          ok = 0;
       end
    end
end
if ok
    seqparm.maxlet = 10; seqparm.nrandperm = 500000; seqparm.orderpairthre = 0.9; seqparm.allpairthre = 0.6;
    input = inputdlg({'Max computing letters'; 'Random permuation number'; 'Ordered pair threshold';...
           'All pair threshold'; 'Neighbor pair threshold'}, 'Sequence finalization parameters', 5, {'10'; '500000'; '0.8'; '0.9'; '0.9'}); 
    if (~isempty(input))
        seqparm.maxlet = str2num(input{1}); seqparm.nrandperm = str2num(input{2}); seqparm.orderpairthre = str2num(input{3}); 
        seqparm.allpairthre = str2num(input{4}); seqparm.neibpairthre = str2num(input{5});
    else
        ok = 0;
    end
end
if ok 
    [fname, pname] = uiputfile(fullfile(cd, '.tmp'), 'Write templates to:');
    if (numel(fname)>1)
       writefilename = fullfile(pname, fname);
    else
       ok = 0;
    end
end
if ok
    filenamenow = writefilename;
    tmp = findalltemplates(pinfo, cellind, crrmode, sessev, ccparm, seqparm); 
    if (~isempty(tmp))
        pppinfo = pinfo;
        pinfo = []; pinfo = tmp; save(filenamenow, 'pinfo', '-mat');
        pinfo = []; pinfo = pppinfo;
    end
    disp(['-----> number of templates generated: ', num2str(numel(tmp))]);
end
disp('**********************');

function tmpout = findalltemplates(pinfo, cellind, crrmode, sessev, ccparm, seqparm)
%get selected cellind (this is actually pairid)
tmpout = []; nt = 0;
%%%%search for all the possible finaldir
allfinaldir = cell(1, numel(cellind));
for (i = 1:numel(cellind))
    aaa = strfind(pinfo.general.finaldir{cellind(i)}, '__'); 
    allfinaldir{i} = pinfo.general.finaldir{cellind(i)}(1:aaa(1)-1);
    %allfinaldir{i} = strtok(pinfo.general.finaldir{cellind(i)}, '_');
end
unifinaldir = unique(allfinaldir);
disp(['-----> number of final dirs: ', num2str(numel(unifinaldir))]);
for (i = 1:numel(unifinaldir))
    fdirnow = unifinaldir{i};
    disp(['-------> final dir now: ', fdirnow]);
    ind = find(strcmp(pinfo.general.crrtype(cellind), 'cross') & strcmp(allfinaldir, unifinaldir{i}));
    pairid = cellind(ind); pairid = sort(pairid); %%%%all pairs in this final directory - sort o make sure cell letter assignment is consistent.
    pairparmfile = pinfo.general.parmfile(pairid); cellparmfile = listcellsinpairs(pairparmfile);
    if (numel(cellparmfile)>62) %%%
        disp('----------> number of cells exceeds template capacity; try to be more selective; aborted');
    elseif (numel(cellparmfile)<4)
        disp('----------> too few cells; aborted');
    else
        cellletter = DE_Seq_AssignLetter(numel(cellparmfile), 'allletters'); paircellid = findcellid(pairparmfile, cellparmfile);
        if strcmp(sessev, 'events')
           evType = pinfo.parm.eventtype{pairid(1)}; evName = pinfo.general.eventname{pairid(1)};
           ind = find( strcmp(evType, 'run') );
           for (j = 1:numel(ind))
               evnamenow = evName{ind(j)}; disp(['----------> event now: ', evnamenow]);
               [ptime, pzsc, consis, kurt] = findcrrparametervalues(pinfo, pairid, crrmode, sessev, evnamenow);
               tmpnow = findtmp_singleevtss(pinfo, ptime, pzsc, consis, kurt, paircellid, cellparmfile, cellletter, fdirnow, evnamenow, ccparm, seqparm, crrmode, sessev);
               if ~isempty(tmpnow)
                  nt = nt + 1; tmp(nt) = tmpnow;
               end
           end
        elseif strcmp(sessev, 'sessions')
           sessName = pinfo.general.sessionname{pairid(1)};
           for (j = 1:numel(sessName))
               sessnow = sessName{j}; disp(['----------> session now: ', sessnow]);
               [ptime, pzsc, consis, kurt] = findcrrparametervalues(pinfo, pairid, crrmode, sessev, sessnow);
               %disp(ptime); disp(pzsc);
               tmpnow = findtmp_singleevtss(pinfo, ptime, pzsc, consis, kurt, paircellid, cellparmfile, cellletter, fdirnow, sessnow, ccparm, seqparm, crrmode, sessev);
               if ~isempty(tmpnow)
                  nt = nt + 1; tmp(nt) = tmpnow;
               end
           end
        end
    end
end
if (nt > 0) tmpout = tmp; end

function tmp = findtmp_singleevtss(pinfo, ptime, pzsc, consis, kurt, paircellid, cellparmfile,  cellletter, fdirnow, evnamenow, ccparm, seqparm, crrmode, sessev)
tmp = [];
maxlet = seqparm.maxlet;
nrandperm = seqparm.nrandperm; %%%number of random sequences to search for the best fit if number of letters > maxlet
orderpairthre = seqparm.orderpairthre; %%%%at least this fraction of pairs are verified by their cross-correlation orders
neibpairthre = seqparm.neibpairthre; 
allpairthre = seqparm.allpairthre; %%%%at least this fraction of pairs among all pairs are verified by their cross-correlation orders
searchoptions = optimset('TolX', 0.001, 'MaxIter', 100000);
tmpfilename = []; bestseq = []; bestscore = 0; finalseq = []; finalscore = []; bestrank = []; finalrank = [];
ok = 1; ind1 = 1:numel(ptime); ind2 = 1:numel(ptime); ind3 = 1:numel(ptime); %ind4 = 1:numel(ptime);
if ~isempty(ccparm.pzsc) ind1 = find(pzsc>=ccparm.pzsc); end
if ~isempty(ccparm.consis) ind2 = find(consis>=ccparm.consis); end
if ~isempty(ccparm.ptime) ind3 = find(abs(ptime)<=ccparm.ptime); end
%disp(ind1); disp(ind2); disp(ind3);
indpass = intersect(intersect(ind1, ind2), ind3);
iii = ones(1, numel(ptime)); %%%just in case cells in some pairs are not identified in paircellid
for (i = 1:numel(ptime)) 
    if (isempty(find(paircellid{i} == 0))) iii(i) = 1; end
end
ind4 = find( iii == 1);
pairind = intersect(indpass, ind4); npair = numel(pairind);  
indorder = cell(1, npair); orderpairletter= cell(1, npair); indrank = cell(1, npair); pairPtime = NaN*ones(1, npair);
allind = []; outfileletter = []; finalorderpair = []; finalpairPtime = [];
for (i = 1:npair)
    if ptime(pairind(i)) >0
        IX = [1 2]; pairPtime(i) = ptime(pairind(i));
    elseif ptime(pairind(i)) <0
        IX = [2 1]; pairPtime(i) = -ptime(pairind(i));
    else
        if rand>0.5 %%%for avoid bias in order assigning
           IX = [1 2]; pairPtime(i) = 0; 
        else
           IX = [2 1]; pairPtime(i) = 0;
        end
    end
    indorder{i} = paircellid{pairind(i)}(IX); allind = [allind paircellid{pairind(i)}];
    orderpairletter{i} = cellletter(indorder{i});    
end
%%%%%%now how to link the indices in indorder into sequences and find th longest
%%%%%%%%%%%%%use a brutal complet search of all the possible sequences and then choose the best match
allind = unique(allind); tmpfilename = cellparmfile(allind); 
nletter = numel(allind); nallpair = nletter*(nletter-1)/2;
if (nletter<4)
    disp(['------------> not enough correlated cells: N = ', num2str(nletter)]);
else
    disp(['------------> number of correlated cells: N = ', num2str(nletter)]);
    if (numel(allind)<=maxlet) %%%%if less than 12 letter, complete permutation
        allseq = perms(allind);
    else %%%is more than 12 letters, random choose 100000
        allseq = zeros(nrandperm, numel(allind));
        for (i = 1:nrandperm)
            allseq(i,:) = allind(randperm(numel(allind)));
        end
    end
    [bestseqind, bestscore, bestrank] = findbestseq(allseq, indorder, pairPtime); %%%bestseqind = [nseq, nletter] contains only sequences of filename indices 
    if bestscore <= 2 %%%if less than 3 pairs matched
        disp('----------------> best matching scores are too low');
    else
    [nseq, nletter] = size(bestseqind); bestseq = cell(nseq, 1);
    for (ij = 1:nseq) bestseq{ij} = cellletter(bestseqind(ij,:)); end
    disp(['----------------> best score is: ', num2str(bestscore), ' out of ', num2str(npair), ' (total pairs: ', num2str(nallpair), ')']);
    outfileletter = cellletter(allind);
    [finalseq, finalscore, finalorderpair, finalpairPtime] = findfinalseq(bestseq, bestscore, bestrank, pairPtime, orderpairletter, outfileletter, orderpairthre, neibpairthre, allpairthre);
    if isempty(finalseq)
        finalrank = []; finalfiterr = [];
        disp('----------------> final sequences did not pass the pair check');
    else
        nfinalletter = numel(finalseq{1});
        [finalrank, finalfiterr] = findfinalrank(finalseq, finalorderpair, finalpairPtime, searchoptions);
        disp(['----------------> final score is: ', num2str(finalscore), ' out of ', num2str(numel(finalorderpair)), ' (total pairs: ', num2str(nfinalletter*(nfinalletter-1)/2), ')']);
       
    %%%%%%output template
    tmp.tmp.tmpfinalseq = finalseq; tmp.tmp.finalrank = finalrank; tmp.tmp.tmpfileletter = outfileletter;
    tmp.tmp.tmpfilename = tmpfilename; tmp.tmp.tmptype = 'crrseq'; tmp.tmp.tmplength = numel(tmpfilename); 
    tmp.tmp.tmpseq = bestseq; tmp.tmp.tmpscore = bestscore; tmp.tmp.tmprank = bestrank;
    tmp.tmp.tmporderpair = orderpairletter; tmp.tmp.tmpfinalscore = finalscore;
    tmp.tmp.tmpfinalorderpair = finalorderpair; tmp.tmp.finalfiterr = finalfiterr;
    tmp.tmp.tmppairPtime = pairPtime; tmp.tmp.tmpfinalpairPtime = finalpairPtime;
    %%%%parameters %%%first aninmaldate
    animaldate = findanimaldate(pinfo, fdirnow);
    tmp.parm.animaldate = animaldate; tmp.parm.sessev = sessev;       
    tmp.evtfile = evnamenow; tmp.parm.crrmode = crrmode; 
    tmp.parm.lapconsisthres = ccparm.consis; tmp.parm.kurtthres = ccparm.kurt; tmp.parm.zscorethres = ccparm.pzsc;
    tmp.parm.maxlet = seqparm.maxlet; tmp.parm.nrandperm = seqparm.nrandperm; 
    tmp.parm.orderpairthre = seqparm.orderpairthre; tmp.parm.allorderthre = seqparm.allpairthre;
    tmp.parm.smooth = []; tmp.parm.smoothsig = [];
    tmp.parm.epindex = NaN; tmp.parm.nepdisplay = 1; tmp.parm.starttime = 0; tmp.parm.endtime = 0;
    tmp.parm.timebinsize = []; tmp.parm.spatialbinsize = [];
    %%%junk parameters
    tmp.parm.peaksearchhalfbinnum = []; tmp.parm.maxlag = []; tmp.parm.option = []; tmp.tmp.crrT = [];
    %%%%dummy variables for rate template
    tmp.tmp.tmpspikecnt = []; tmp.tmp.tmpspiketime = [];
    end
    end
end

function [finalseq, finalscore, finalorderpair, finalpairPtime] = findfinalseq(bestseq, bestscore, bestrank, pairPtime, orderpair, fileletter, orderpairthre, neibpairthre, allpairthre)
remlet = findremlet(bestseq, orderpair, bestscore, bestrank, orderpairthre, neibpairthre, allpairthre);
while ~isempty(remlet) 
     nseq = numel(bestseq); npair = numel(orderpair); nletter = numel(bestseq{1}); remind = [];
     for (i = 1:npair)
         if ~isempty(strfind(orderpair{i}, remlet)) remind = union(remind, i); end
     end
     IX = setdiff( (1:npair), remind ); orderpair = orderpair(IX); pairPtime = pairPtime(IX);
     indorder = []; for (i = 1:numel(orderpair)) indorder{i} = fileind(orderpair{i}, fileletter); end
     allseq = zeros(nseq, nletter-1);
     for (i = 1:nseq)
         ii = find(bestseq{i} == remlet); IY = setdiff((1:nletter), ii);
         allseq(i,:) = fileind(bestseq{i}(IY), fileletter); 
     end
     if (isempty(allseq))
         break
     else
         allseq = unique(allseq, 'rows');
         [bestseqind, bestscore, bestrank] = findbestseq(allseq, indorder, pairPtime); %%%bestseqind = [nseq, nletter] contains only sequences of filename indices 
         [nseq, ~] = size(bestseqind); bestseq = cell(nseq, 1);
         for (ij = 1:nseq) bestseq{ij} = fileletter(bestseqind(ij,:)); end
         [bestseq, iii] = unique(bestseq); bestrank = bestrank(iii);
         remlet = findremlet(bestseq, orderpair, bestscore, bestrank, orderpairthre, neibpairthre, allpairthre);
     end
end
finalseq = bestseq; finalscore = bestscore; finalorderpair = orderpair; finalpairPtime = pairPtime;
%%%%final check: how tight finalseq fits the ordered pairs
checkscore = zeros(1, numel(finalseq));
for (i = 1:numel(finalseq))
    checkscore(i) = checkallpairs(finalseq{i}, orderpair);
end
%iii = find(checkscore>=allpairthre); finalseq = finalseq(iii);
iii = find(checkscore ==max(checkscore)); finalseq = finalseq(iii); 
disp(['------------> (exp.) final sequence: ', finalseq{1}, '(neighbor pair score: ', num2str(max(checkscore)), ')']);
function remlet = findremlet(bestseq, orderpair, bestscore, bestrank, orderpairthre, neibpairthre, allpairthre)
nletter = numel(bestseq{1}); qua1 = bestscore/numel(orderpair); qua2 = bestscore/(nletter*(nletter-1)/2);
remlet =[];
letterid = sort(bestseq{1}); nletter = numel(letterid); 
floatdist = zeros(nletter,1);
for (i = 1:nletter) floatdist(i) = findfloatdist(letterid(i),bestseq); end
checkscore = zeros(1, numel(bestseq));
for (i = 1:numel(bestseq)) checkscore(i) = checkallpairs(bestseq{i}, orderpair); end
if ~(uniquerank(bestrank, bestseq, letterid, floatdist) && (qua1>=orderpairthre) && ((qua2>=allpairthre)|| (max(checkscore)>=neibpairthre))) %%%if not satisfied, stop 
%if (qua1<orderpairthre) || (numel(bestseq)>1)     
%    [~, ii] = max(floatdist); 
    nincludepair = zeros(nletter,1);
    for (i = 1:nletter)
        nincludepair(i) = findincludepair(letterid(i), orderpair);
    end
    [~, ii] = min(nincludepair); 
    if (min(nincludepair) < max(nincludepair)) 
        remlet = letterid(ii(1));
%     if (min(floatdist) < max(floatdist))
%         remlet = letterid(ii(1));
    else %if all letters have equal number of pairs: remove the most floating letters
%         nincludepair = zeros(nletter,1);
%         for (i = 1:nletter)
%             nincludepair(i) = findincludepair(letterid(i), orderpair);
%         end
%         [~, ii] = min(nincludepair); 
%         if (min(nincludepair) < max(nincludepair)) 
%             remlet = letterid(ii(1));
%         end
        [~, ii] = max(floatdist); 
        if (min(floatdist) < max(floatdist)) 
            remlet = letterid(ii(1));
        end
    end
end
function out = uniquerank(bestrank, bestseq, letterid, floatdist) %%%if the only floaters are the 0-peak pairs: they shold next to each other
out = 1;
if (numel(bestrank) < 2)
   nextrankcheck = 1; floatcheck = 1;
   flet = letterid(floatdist>0); 
   for (i = 1:numel(bestrank))
       ur = unique(bestrank{i}); nextornot = zeros(1, numel(ur)); eletnow = [];
       for (j = 1:numel(ur))
        iii = find(bestrank{i} == ur(j));
        if (numel(iii)>=2)
            eletnow = union(eletnow, bestseq{i}(iii));
            if ~isempty(find( diff(iii) ~= 1 )) nextornot(j) = 1; end
        end
       end
       if ~isempty(find(nextornot == 1))
        nextrankcheck = 0; break
       end
       if ~isempty(find(~ismember(flet, eletnow)))
        floatcheck = 0; break
       end
   end
   out = nextrankcheck & floatcheck; 
end

% yes = 1;
% for (i = 1:numel(bestrank))
%     for (j = i+1:numel(bestrank))
%         if ~isempty(find((bestrank{i}-bestrank{j}) ~= 0))
%             yes = 0; break
%         end
%     end
% end
function frac = checkallpairs(bestseq, orderpair)
frac = 0; nlet = numel(bestseq); 
if (nlet > 1)
   mm = zeros(1,nlet-1);
   for (i = 1:nlet-1)
        pairnow = bestseq(i:i+1); if (~isempty(find(strcmp(orderpair, pairnow)))) mm(i)=1; end
   end
   frac = numel(find(mm==1))/(nlet-1);
end
function floatdist = findfloatdist(letterid, bestseq)
ind = ones(1,numel(bestseq));
for (i=1:numel(bestseq)) ind(i) = find(bestseq{i} == letterid); end
floatdist = max(ind) - min(ind);
function ncheck = checkpairs(orderpair, seq)
nincludepair = zeros(1, numel(seq)); ncheck = 1;
for (i = 1:numel(seq))
    nincludepair(i) = findincludepair(seq(i), orderpair);
end
if (numel(find(nincludepair == max(nincludepair))) == numel(nincludepair)) ncheck = 0; end
function nincludepair = findincludepair(letterid, orderpair)
nincludepair = 0; npair = numel(orderpair); 
for (i = 1:npair)
    if ~isempty(strfind(orderpair{i}, letterid))
        nincludepair = nincludepair+1;
    end
end
function ind = fileind(seq, fileletter)
ind = zeros(1, numel(seq));
for i=1:numel(seq)
    ind(i) = find(fileletter == seq(i));
end

function [bestseq, bestscore, bestrank] = findbestseq(allseq, indorder, pairPtime) %%%output bestseq contains (file) indices not (file) letters
[nseq, nlet] = size(allseq); sscor = zeros(1, nseq); npair = numel(indorder); 
for (i = 1:nseq)
    seqnow = allseq(i,:); pairss = zeros(1, npair); ind1now{i} = zeros(1, npair); ind2now{i} = zeros(1, npair); Ptime{i} = zeros(1, npair);
    for (j = 1:npair)
        ind1 = find(seqnow == indorder{j}(1)); ind2 = find(seqnow == indorder{j}(2));
        judg = (ind2-ind1) * pairPtime(j);
        if (judg>=0) pairss(j) = 1; end     
    end
    sscor(i) = numel(find(pairss == 1)); 
end
bestscore = max(sscor); iii = find(sscor == bestscore); bestseq = allseq(iii, :); %%%output bestseq contains (file) indices not (file) letters
%%%%rank the cells: not easy - find equal groups,
cellgroup = []; ngroup = 0; bestrank = cell(1, numel(iii));
for (j = 1:npair)
    if (pairPtime(j) ==0) 
        belonggroup = 0;
        for (i = 1:ngroup)
            if ismember(indorder{j}(1), cellgroup{i}) || ismember(indorder{j}(2), cellgroup{i})
                cellgroup{i} = union(cellgroup{i}, [indorder{j}(1) indorder{j}(2)]); belonggroup = 1;
            end
        end
        if (belonggroup == 0)
            ngroup = ngroup + 1; cellgroup{ngroup} = [indorder{j}(1) indorder{j}(2)];
        end
    end
end
for (i = 1:numel(iii))
    bestrank{i} = getranknow(bestseq(i,:), cellgroup);
end
function bestrank = getranknow(bestseq, cellgroup)
bestrank = 1:numel(bestseq); ngroup = numel(cellgroup);
for (i = 1:ngroup)
    ind = zeros(1, numel(cellgroup{i}));
    for (j = 1:numel(cellgroup{i}))
        ind(j) = find(bestseq == cellgroup{i}(j));
    end
    bestrank(ind) = mean(bestrank(ind))*ones(1, numel(ind));
end

function [bestrank, fiterr] = findfinalrank(bestseq, orderpair, pairPtime, searchoption)
%%%%%do a parameter search for optimal ranks
nseq = numel(bestseq); nlet = numel(bestseq{1}); npair = numel(orderpair);
bestrank = cell(1, nseq); fiterr = NaN*ones(1, nseq);
for (i = 1:nseq)
    bestrank{i} = NaN*ones(1, nlet); ind1 = zeros(1, nlet); ind2 = zeros(1, nlet); 
    for (j = 1:npair)
        id1 = find(bestseq{i} == orderpair{j}(1)); id2 = find(bestseq{i} == orderpair{j}(2));
        if abs(id1-id2) < 3
            ind1(j) = id1; ind2(j) = id2;
        end
    end
    iii = find(ind1>0); ind1 = ind1(iii); ind2 = ind2(iii); Pt = pairPtime(iii);
    start = getstartposition(ind1, ind2, Pt, nlet);
    [outrank, fval, exitflag] = fminsearch(@(par)intervalfit(par, ind1, ind2, Pt), start, searchoption); 
    if (exitflag ~= 1)
       disp(['---------> warning: a successfule fit is not found: ', num2str(fval)]); 
    else
       fiterr(i) = fval; bestrank{i} = outrank-min(outrank); %bestrank{i}(2:nlet) = outrank;
    end    
end
function err = intervalfit(cellrank, pairind1, pairind2, pairintervals)
%allrank = [firstval cellrank]; 
allrank = cellrank;
err = sum(abs(pairintervals-(allrank(pairind2)-allrank(pairind1))))/numel(pairind1);
function start = getstartposition(ind1, ind2, Pt, nlet)
start = zeros(1, nlet); int = zeros(1, nlet-1);
for (i = 1:nlet-1)
    idnow = find( (ind1 == i) & (ind2 == i+1) );
    if ~isempty(idnow)
       int(i) = Pt(idnow);
    end
end
for (i = 1:nlet-1)
    start(i+1) = sum(int(1:i));
end
% for (i = 1:numel(pairind1))
%     err = err + abs(pairintervals(i) - (allrank(pairind2)-allrank(pairind1)));
% end

% %%%%rank the cells: not easy - find equal groups,
% cellgroup = []; ngroup = 0; bestrank = cell(1, numel(iii));
% for (j = 1:npair)
%     if (indrank{j}(1) == indrank{j}(2)) 
%         belonggroup = 0;
%         for (i = 1:ngroup)
%             if ismember(indorder{j}(1), cellgroup{i}) || ismember(indorder{j}(2), cellgroup{i})
%                 cellgroup{i} = union(cellgroup{i}, [indorder{j}(1) indorder{j}(2)]); belonggroup = 1;
%             end
%         end
%         if (belonggroup == 0)
%             ngroup = ngroup + 1; cellgroup{ngroup} = [indorder{j}(1) indorder{j}(2)];
%         end
%     end
% end
% for (i = 1:numel(iii))
%     bestrank{i} = getranknow(bestseq(i,:), cellgroup);
% end
% function bestrank = getranknow(bestseq, cellgroup)
% bestrank = 1:numel(bestseq); ngroup = numel(cellgroup);
% for (i = 1:ngroup)
%     ind = zeros(1, numel(cellgroup{i}));
%     for (j = 1:numel(cellgroup{i}))
%         ind(j) = find(bestseq == cellgroup{i}(j));
%     end
%     bestrank(ind) = mean(bestrank(ind))*ones(1, numel(ind));
% end

function cellparmfile = listcellsinpairs(allparmfile)
pfile = cell(2*numel(allparmfile), 1); 
for (i = 1:numel(allparmfile))
    pnow = allparmfile{i};
    ss = strfind(pnow, '__'); pfile{2*i-1} = pnow(1:ss-1); pfile{2*i} = pnow(ss+2:numel(pnow));
end
cellparmfile = unique(pfile);

function [ptimefield, pzscfield, consisfield, kurtfield] = findevfieldname(crrmode)
if strcmp(crrmode, 'crr')
    ptimefield = 'sess1stPtime'; pzscfield = 'sess1stPZsc'; consisfield = [];
elseif strcmp(crrmode, 'lapcrr') || strcmp(crrmode, 'lapspatialcrr')
    ptimefield = 'evtLapMean1stPtime'; pzscfield = 'evtLapMeanPZsc'; consisfield = 'evtLapcrrMeanR';
end

function paircellid = findcellid(pairparmfile, cellparmfile)
npair = numel(pairparmfile); paircellid = cell(1,npair);
for (i = 1:npair)
    paircellid{i} = zeros(1,2);
    pnow = pairparmfile{i};
    ss = strfind(pnow, '__'); pfile1 = pnow(1:ss-1); pfile2 = pnow(ss+2:numel(pnow));
    p1 = find(strcmp(cellparmfile, pfile1)); p2 = find(strcmp(cellparmfile, pfile2));
    if (numel(p1) == 1) && (numel(p2) == 1)
        paircellid{i} = [p1 p2];
    else
        disp('------------> warning: unable to identify file names of a correlation pair');
    end
end

function [ptime, pzsc, consis, kurt] = findcrrparametervalues(pinfo, pairid, crrmode, sessev, evnamenow)
npair = numel(pairid); ptime = NaN*ones(npair, 1); pzsc = NaN*ones(npair, 1); 
consis = NaN*ones(npair, 1); kurt = NaN*ones(npair,1);
if isfield(pinfo, crrmode)
    if strcmp(crrmode, 'crr')
       if strcmp(sessev, 'sessions')
          for (i = 1:npair)
              allevnames = pinfo.general.sessionname{pairid(i)};
              if (~isempty(allevnames)) break; end
          end
       elseif strcmp(sessev, 'events')
           for (i = 1:npair)
                allevnames = pinfo.general.eventname{pairid(i)};
                if ~isempty(allevnames) break; end
           end
       end
       evind = find(strcmp(allevnames, evnamenow));
       if numel(evind) ==1
          if isfield(pinfo.crr, 'sess1stPtime')
              ptime = getfieldvalues(pinfo, 'crr', 'sess1stPtime', pairid, evind);
          end
          if isfield(pinfo.crr, 'sess1stPZsc')
              pzsc = getfieldvalues(pinfo, 'crr', 'sess1stPZsc', pairid, evind);
          end
          if isfield(pinfo.crr, 'sessKurt')
              kurt = getfieldvalues(pinfo, 'crr', 'sessKurt', pairid, evind);
          end
       end
    elseif strcmp(crrmode, 'lapcrr') || strcmp(crrmode, 'lapspatialcrr')
       for (i = 1:npair)
           allevnames = pinfo.(crrmode).evtName{pairid(i)};
           if ~isempty(allevnames) break; end
       end
       evind = find(strcmp(allevnames, evnamenow));
       if numel(evind) ==1
          if isfield(pinfo.(crrmode), 'evLapMean1stPtime')
              ptime = getfieldvalues(pinfo, crrmode, 'evLapMean1stPtime', pairid, evind);
          end
          if isfield(pinfo.(crrmode), 'evLapMean1stPZsc')
              pzsc = getfieldvalues(pinfo, crrmode, 'evLapMean1stPZsc', pairid, evind);
          end
          if isfield(pinfo.(crrmode), 'evLapMeanKurt')
              kurt = getfieldvalues(pinfo, crrmode, 'evLapMeanKurt', pairid, evind);
          end
          if isfield(pinfo.(crrmode), 'evLapcrrMeanR')
              consis = getfieldvalues(pinfo, crrmode, 'evLapcrrMeanR', pairid, evind);
          end
      end
   end
end
function ptime = getfieldvalues(pinfo, cat, var, pairid, evind)
npair = numel(pairid); ptime = NaN*ones(npair,1);
for (i = 1:npair)
    if ~isempty(pinfo.(cat).(var){pairid(i)})
        if iscell(pinfo.(cat).(var){pairid(i)})
           ptime(i) = pinfo.(cat).(var){pairid(i)}{evind};
        else
           ptime(i) = pinfo.(cat).(var){pairid(i)}(evind);
        end
    end
end

function animaldate = findanimaldate(pinfo, fdirnow)
animaldate = [];
npair = numel(pinfo.general.finaldir); animname = cell(1, npair); datedir = cell(1, npair); fdir = cell(1, npair);
for (i = 1:npair)
    pnow = pinfo.general.finaldir{i}; ss = strfind(pnow, '__'); fdir{i} = pnow(1:ss-1);
    datedir{i} = pinfo.general.datedir{i}; %ss = strfind(pnow, '__'); datedir{i} = pnow(1:ss-1);
    pnow = pinfo.general.animalname{i}; ss = strfind(pnow, '__'); animname{i} = pnow(1:ss-1);
end
iii = find(strcmp(fdir, fdirnow));
if ~isempty(iii)
    animaldate = strcat(animname{iii(1)}, '_', datedir{iii(1)});
end





