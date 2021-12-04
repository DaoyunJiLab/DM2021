function DM_Seq_CompareVariables_Callback
%%%For selected group of items (templates or events) after seqeuncing/decoding, 
%%%    compare selected variables among groups and between a group and its shuffles
%%%The plot is specific for the selected events (one or more)

hf = gcbf; pinfonow = getappdata(hf, 'pinfo'); datanow = getappdata(hf, 'data'); tagnow = get(gcbo, 'Tag');
plotparm = getappdata(hf, 'plotparm'); title = get(hf, 'Name');
%hspike = getappdata(hf, 'hspike'); spikeselection = getappdata(hspike, 'selection'); cellind = find(spikeselection==1);
hgroup = getappdata(hf, 'hgroup'); groupselection = getappdata(hgroup, 'selection'); 
grpind = find(groupselection == 1); groupnames = datanow.grouplist.groupname(grpind);
groupcellind = datanow.grouplist.groupindex(grpind); ngroup = numel(groupnames);
ok = 1; cind = [];
disp('-----> Comparing hidden variables across groups and with shuffles ......');
if (numel(grpind)<1) disp('-----------> no groups selected'); ok= 0; end
for (kk = 1:numel(grpind)) cind = union(cind, datanow.grouplist.groupindex{grpind(kk)}); end    
ok = 1;
if (numel(cind)<1) disp('-----------> no groups/templates selected'); ok= 0; end
if ok
    seqtype = [];
    if isfield(pinfonow.parm, 'seqtype') seqtype = pinfonow.parm.seqtype{cind(1)}; end
    %    ok = 0; disp('-----------> aborted: hidden variable comparison is not available for event-itemized db; use original sequence db to compare');
    %else
    if ~isfield(datanow, 'data')
        ok = 0; disp('-----------> aborted: data not available for stripped databases');
    end
end
if ok
   if (plotparm.evselect == 0) %%%select sessions
      disp('-----------> aborted: hidden variable comparison cannot be computed: need to select the event option'); ok = 0;
   else
      shufexist = 1;
      if ~isfield(datanow.data, 'shufevseq')
         disp('-----------> warning: comparisons to shuffles not computed: shuffle data not exit or saved'); shufexist = 0;
      end
   end
end
if ok %%%%fruther select particular events to compute
        someparm.evkeyword = []; someparm.evkeytype = []; someparm.evkeynoword = []; someparm.evkeynotype = [];
        input = inputdlg({'Event keyword (self = template event)'; 'Event type'; 'Event keyNOword (self to exclude template event)'; 'Event NOtype'}, 'Further event selection', 4, {'self'; ''; 'first'; ''}); 
        if (~isempty(input))
           someparm.evkeyword = input{1}; someparm.evkeytype = input{2}; someparm.evkeynoword = input{3}; someparm.evkeynotype = input{4};
        else
           ok = 0;
        end
end
if ok %%% checkin selected seq variables
    seqnow = datanow.data.evseq{groupcellind{1}(1)}{1}; %%%%Use this to generate a variable list
    oldvar = fieldnames(seqnow); newquant = fieldnames(seqnow.newquant{1,1});
    varlist = setdiff(union(newquant, oldvar), excludevars);
    [ind, ok] = listdlg('ListString', varlist, 'Name', 'Variables to plot');
    if ~ok
        disp('-----------> no variables selected or cancelled');
    else
        varnames = varlist(ind);
    end
end
if ok %%% checkin selected seq variables
    varlist = {'hist'; 'cumu'; 'bar'};
    [ind, ok] = listdlg('ListString', varlist, 'Name', 'Plot type');
    if ~ok
        disp('-----------> no plot type selected or cancelled');
    else
        plottype = varlist{ind(1)};
    end
end
if ok     %%%%%%%pre-count how many tempaltes over how many events == nitems
    if ~contains(seqtype, 'evtitemized')
       evName = pinfonow.general.eventname; evType = pinfonow.parm.eventType; tmpevname = pinfonow.tmp.sessevname;
    else
       evName = pinfonow.general.tmpID; evType = pinfonow.parm.tmpType; tmpevname = pinfonow.general.eventname;
    end
    nitems = zeros(1, ngroup);
    for (pp = 1:ngroup)
        cellind = groupcellind{pp}; nshuffle = pinfonow.parm.Nshuffle{cellind(1)};
        for (ttt=1:numel(cellind))
            k = cellind(ttt); 
            for (j = 1:numel(evName{k}))
                   if filterevents(tmpevname{k}, evName{k}{j}, evType{k}{j}, someparm.evkeyword, someparm.evkeytype, someparm.evkeynoword, someparm.evkeynotype)
                      nitems(pp) = nitems(pp) + 1;
                   end
            end
        end
    end
end
%%%%%%%%%%%%Calculate and plot here
if ok  
for (kk = 1:numel(varnames)) %%%for each seq or seq2D variable
    varnow = varnames{kk};
    val = cell(1, ngroup); shufval = cell(1, ngroup);
    for pp = 1:ngroup
        val{pp} = cell(1, nitems(pp)); shufval{pp} = cell(1,nitems(pp)); 
        cellind = groupcellind{pp};
        nnow = 0;
        for (ttt = 1:numel(cellind))
            k = cellind(ttt);
            for (j = 1:numel(evName{k}))
                if filterevents(tmpevname{k}, evName{k}{j}, evType{k}{j}, someparm.evkeyword, someparm.evkeytype, someparm.evkeynoword, someparm.evkeynotype)
                   nnow = nnow + 1;
                   if ~isempty(find(strcmp(oldvar, varnow))) %%%if a variable is in the old format (decoding and sequencing algorithm) 
                      val{pp}{nnow} = datanow.data.evseq{k}{j}.(varnow); %%%evseq{ncell}{nev}.(varnames){1:nlap,1}(1 or ntarget)
                      if shufexist
                          nnn=randperm(nshuffle); shufcopynum = nnn(1); %%%takes a random copy (equalize shuffle values and actual values)
                         shufval{pp}{nnow} = datanow.data.shufevseq{k}{shufcopynum, j}.(varnow); %%% shufevseq{ncell}{nshuffle, nev}.(varnames){1:nlap,1}(1 or ntarget)
                      end
                   else
                      val{pp}{nnow} = findnewvarvalues(datanow.data.evseq{k}{j}.newquant, varnow); %%%evseq{ncell}{nev}.newquant{1:nlap,1}.(varnow)(1 or ntarget)
                       if shufexist
                          nnn=randperm(nshuffle); shufcopynum = nnn(1); %%%takes a random copy (equalize shuffle values and actual values)
                         shufval{pp}{nnow} = findnewvarvalues(datanow.data.shufevseq{k}{shufcopynum, j}.newquant, varnow); %%% shufevseq{ncell}{nshuffle, nev}.newquant{1:nlap,1).(varnow)(1 or ntarget)
                      end
                   end
                end
            end
        end
    end
    for mm = 1:ngroup
        for nn =mm+1:ngroup
            indnow = [mm nn];
            histplot(val(indnow), nitems(indnow), varnow, groupnames(indnow), plotparm, plottype);
        end
    end
    if shufexist
        for mm =1:ngroup
            valnow{1} = val{mm}; valnow{2} = shufval{mm}; gname{1} = groupnames{mm}; gname{2} = [groupnames{mm} '-oneShufCopy'];
            histplot(valnow, nitems([mm mm]), varnow, gname, plotparm, plottype);
        end
    end
end
end
disp('*******************************');

function var = excludevars
var = {'newquant'; 'Xgrid'; 'allletter'; 'alltime'; 'ep'; 'evname'; 'quantdata'; 'ratestd'; 'seq'; 'seqdata';...
    'seqstart'; 'seqend'; 'seqmarker'; 'seqtime'; 'spatpeakth'; 'spatratemean'; 'spatratestd'; 'timepeakth'; 'timepoint'};

% function val = findnewvarvalues(newquant, varnow)
% nlap = numel(newquant); %%val[1:nlap]; newquant{1:nlap,1}.(varnow){1}(1)
% val = NaN*ones(1,nlap); 
% for (i=1:nlap)
%     val(i) = newquant{i}.(varnow);
% end
function val = findnewvarvalues(newquant, varnow)
nlap = numel(newquant); %%val[1:nlap]; newquant{1:nlap,1}.(varnow)(1 or ntarg if free sequencing)
val = cell(1,nlap); 
for (i=1:nlap)
    val{i} = newquant{i}.(varnow);
end
val = cell2mat(val);
function histplot(val, items, plotname, groupname, plotparm, tagmark)
%val{n}{nitem}{1:nlap,1}
%%%%arrange data and compute data statistics first
n = numel(val); mm = NaN*ones(1,n); ss = NaN*ones(1,n); nn = NaN*ones(1,n); m25 = NaN*ones(1,n); m50 = NaN*ones(1,n); m75 = NaN*ones(1,n);
for (i = 1:2) %group index
    valnow = []; nel = 0;
    for (j = 1:numel(val{i})) %cell index
        if (isnumeric(val{i}{j}))
            for (k = 1:numel(val{i}{j}))
                nel = nel + 1; valnow(nel) = val{i}{j}(k);
            end
        elseif (ischar(val{i}{j}))
            nel = nel + 1; valnow(nel) = str2num(val{i}{j}); 
        elseif iscell(val{i}{j})
            for (pp = 1:numel(val{i}{j}))
                if (isnumeric(val{i}{j}{pp}))
                    for (k = 1:numel(val{i}{j}{pp}))
                         nel = nel + 1; valnow(nel) = val{i}{j}{pp}(k);
                    end
                elseif ischar(val{i}{j}{pp})
                    valtt = str2num(val{i}{j}{pp});
                    if (~isempty(valtt)) 
                       nel = nel + 1;  valnow(nel) = valtt;
                    end
                end
            end
        end
    end
    val{i} = valnow; 
    if (nel > 0)
        valnow = valnow(~isnan(valnow));
        mm(i) = mean(valnow); nn(i) = numel(valnow); ss(i) = std(valnow)/sqrt(nn(i)); m50(i) = median(valnow);
        m25(i) = prctile(valnow, 25); m75(i) = prctile(valnow, 75); 
    end
end
%%%%also need to do statistical tests: ransum test and unpaired t-tests
str = dotandranksumtest (val{1}, val{2}); co{1} = [0 0 0]; co{2} = [0 0 0];
%%%%plot below
if (strcmp(tagmark, 'hist'))
   binvector = [plotparm.range(1):plotparm.bin:plotparm.range(2)]; 
   hf = figure('Name', strcat(plotname, '-', tagmark)); hax = axes('Parent', hf, 'NextPlot', 'add'); 
   for (i = 1:n)
       co{i+2} = rand(1,3); str{i+2} = strcat(groupname{i}, ': empty data');
       if ~isempty(val{i})
         Y = histc(val{i}, binvector); 
         if (plotparm.setlog == 1)
           binvector = log10(binvector);
         elseif (plotparm.setlog == 2)
           Y = log10(Y);
         elseif (plotparm.setlog == 3)
           binvector = log10(binvector); Y = log10(Y);
         end
         bar(hax, binvector, Y, 'EdgeColor', co{i+2}, 'FaceColor', co{i+2}); 
         str{i+2} = strcat(groupname{i}, ': mean=', num2str(mm(i)), ';se=', num2str(ss(i)), ';n=', num2str(nn(i)), ...
           ';p25=', num2str(m25(i)), ';p50=', num2str(m50(i)), ';p75=', num2str(m75(i)));
       end
   end
   for (i = 1:numel(str))
       text('Parent', hax, 'String', str{i}, 'Interpreter', 'none', 'Color', co{i}, 'Units', 'normalized', 'Position', [0.05 0.96-(i-1)*0.04]);
   end 
elseif (strcmp(tagmark, 'bar')) %%%bar plot only do mean +- se values
   binvector = 1:n; cl = rand(1,3);
   hf = figure('Name', strcat(plotname, '-', tagmark)); hax = axes('Parent', hf, 'NextPlot', 'add'); str = [];
   Y = mm; UP = Y +ss; DOWN = Y- ss;
   if (plotparm.setlog == 1)
       binvector = log10(binvector);
   elseif (plotparm.setlog == 2)
       Y = log10(Y); UP = log10(UP); DOWN = log10(DOWN);
   elseif (plotparm.setlog == 3)
       binvector = log10(binvector); Y = log10(Y); UP = log10(UP); DOWN = log10(DOWN);
   end
   bar(hax, binvector, Y, 'EdgeColor', cl, 'FaceColor', cl); 
   Drawerrorupdown(binvector, Y, UP, DOWN, hax, [1 0 0]);
   str = cell(1, n);
   for (i = 1:n)
       if ~isempty(val{i})
          str{i} = strcat(groupname{i}, ': mean=', num2str(mm(i)), ';se=', num2str(ss(i)), ';n=', num2str(nn(i)), ...
           ';p25=', num2str(m25(i)), ';p50=', num2str(m50(i)), ';p75=', num2str(m75(i)));
       else
          str{i} = strcat(groupname{i}, ': empty data');
       end
   end
   for (i = 1:numel(str))
       text('Parent', hax, 'String', str{i}, 'Interpreter', 'none', 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.05 0.96-(i-1)*0.04]);
   end 
elseif (strcmp(tagmark, 'cumu'))
   binvector = [plotparm.range(1):plotparm.bin:plotparm.range(2)];  
   hf = figure('Name', strcat(plotname, '-', tagmark)); hax = axes('Parent', hf, 'NextPlot', 'add'); str = [];
   for (i = 1:n)
       co{i+2} = rand(1,3); str{i+2} = strcat(groupname{i}, ': empty data');
       if ~isempty(val{i})
         Y = histc(val{i}, binvector); 
         line(binvector, cumsum(Y), 'Parent', hax, 'Color', co{i+2}); 
         str{i+2} = strcat(groupname{i}, ': mean=', num2str(mm(i)), ';se=', num2str(ss(i)), ';n=', num2str(nn(i)), ...
           ';p25=', num2str(m25(i)), ';p50=', num2str(m50(i)), ';p75=', num2str(m75(i)));
       end
   end
   for (i = 1:numel(str))
       text('Parent', hax, 'String', str{i}, 'Interpreter', 'none', 'Color', co{i}, 'Units', 'normalized', 'Position', [0.05 0.96-(i-1)*0.04]);
   end 
end

function str = dotandranksumtest (val1, val2)
str{1} = 'ttest not performed: empty samples'; str{2} = 'ranksum test not performed: empty samples';
val1 = val1(~isnan(val1)); val2=val2(~isnan(val2));
if (~isempty(val1)) && (~isempty(val2)) 
   [h,p,ci,stats] = ttest2(val1, val2);
   str{1} = ['t-test: p=' num2str(p, '%10.5e') '; t=' num2str(stats.tstat) '; df=' num2str(stats.df)];
   [p,h,stats] = ranksum(val1, val2);
   str{2} = ['ranksum test: p=' num2str(p, '%10.5e') '; W(N1,N2)=' num2str(stats.ranksum) '(' num2str(numel(val1)) ',' num2str(numel(val2)) ')'...
           ';U=' num2str(stats.ranksum-(numel(val1)*(numel(val2)+1)/2))];
end
% 
% function [P, candevind, matchevind, posevind, negevind] = find1Dquantifications(evseq, tottime, seqparm, nshuffle)
% %%%This functon borrowed from the original quantification function in Bayesian decoding_events.
% %%%Some aspects are re-arranged: nev = nshuffle and evseq{i=nshuffle} for just one event
% %%%Now change evseq{i} = evseq{nshuffle}; nev = nshuffle 
% %%%%%Parameters: candidate events, significant positive/negative events
% nev = nshuffle;
% P.totalN = NaN*ones(1, nev); P.candN = NaN*ones(1, nev); P.posN = NaN*ones(1, nev); P.negN = NaN*ones(1, nev); P.matchN =NaN*ones(1,nev);
% P.evTimeRate = NaN*ones(1, nev); P.candTimeRate = NaN*ones(1, nev); P.posTimeRate = NaN*ones(1, nev); P.negTimeRate = NaN*ones(1, nev); P.matchTimeRate = NaN*ones(1, nev);
% P.candEvtRate = NaN*ones(1, nev); P.posEvtRate = NaN*ones(1, nev); P.negEvtRate = NaN*ones(1, nev);  P.matchEvtRate = NaN*ones(1, nev);
% P.posRatio = NaN*ones(1, nev); P.negRatio = NaN*ones(1, nev); P.matchRatio = NaN*ones(1, nev);
% candevind = cell(1,nev); matchevind=cell(1,nev); posevind=cell(1,nev); negevind=cell(1,nev);
% %%%% Added new: mean/median matchscore, mean matchZscore
% P.tottime = tottime; P.candmeanScore = NaN*ones(1, nev); P.candmedScore = NaN*ones(1, nev); P.candmeanZscore = NaN*ones(1, nev); P.candmedZscore = NaN*ones(1, nev);
% P.candmeanAbsScore = NaN*ones(1, nev); P.candmedAbsScore = NaN*ones(1, nev); P.candmeanAbsZscore = NaN*ones(1, nev); P.candmedAbsZscore = NaN*ones(1, nev);
% P.posmeanScore = NaN*ones(1, nev); P.posmedScore = NaN*ones(1, nev); P.posmeanZscore = NaN*ones(1, nev); P.posmedZscore = NaN*ones(1, nev);
% P.negmeanScore = NaN*ones(1, nev); P.negmedScore = NaN*ones(1, nev); P.negmeanZscore = NaN*ones(1, nev); P.negmedZscore = NaN*ones(1, nev);
% P.matchmeanScore = NaN*ones(1, nev); P.matchmedScore = NaN*ones(1, nev); P.matchmeanZscore = NaN*ones(1, nev);  P.matchmedZscore = NaN*ones(1, nev);
% P.matchmeanAbsScore = NaN*ones(1, nev); P.matchmedAbsScore = NaN*ones(1, nev); P.matchmeanAbsZscore = NaN*ones(1, nev); P.matchmedAbsZscore = NaN*ones(1, nev);
% %%%% Added new: mean/median decodeerr
% P.candmeanErr = NaN*ones(1, nev); P.candmedErr = NaN*ones(1, nev); P.posmeanErr = NaN*ones(1, nev); P.posmedErr = NaN*ones(1, nev);
% P.negmeanErr = NaN*ones(1, nev); P.negmedErr = NaN*ones(1, nev); P.matchmeanErr = NaN*ones(1, nev); P.matchmedErr = NaN*ones(1, nev);
% for (i = 1:nev)
%         [mm,nn] = size(evseq{i}.shufposP); shufposP = reshape(cell2mat(evseq{i}.shufposP), mm*nn, 1); 
%         shufnegP = reshape(cell2mat(evseq{i}.shufnegP), mm*nn, 1); nact = reshape(cell2mat(evseq{i}.nactivecell), mm*nn, 1);
%         candevind{i} = find(nact>=seqparm.minactivecell); 
%         posevind{i} = find( (shufposP<seqparm.significancelevel) & (nact>=seqparm.minactivecell)); 
%         negevind{i} = find( (shufnegP<seqparm.significancelevel) & (nact>=seqparm.minactivecell));
%         matchevind{i} = sort(union(posevind{i}, negevind{i}));
%     P.totalN(i) = numel(shufposP); P.candN(i) = numel(candevind{i}); P.posN(i) = numel(posevind{i}); P.negN(i) = numel(negevind{i}); P.matchN(i) = numel(matchevind{i});
%     P.evTimeRate(i) = P.totalN(i)/tottime; P.candTimeRate(i) = P.candN(i)/tottime; 
%     P.posTimeRate(i) = P.posN(i)/tottime; P.negTimeRate(i) = P.negN(i)/tottime; P.matchTimeRate(i) = P.matchN(i)/tottime;
%     P.candEvtRate(i) = P.candN(i)/P.totalN(i); P.posEvtRate(i) = P.posN(i)/P.totalN(i); P.negEvtRate(i) = P.negN(i)/P.totalN(i); P.matchEvtRate(i) = P.matchN(i)/P.totalN(i);
%     P.posRatio(i) = P.posN(i)/P.candN(i); P.negRatio(i) = P.negN(i)/P.candN(i);  P.matchRatio(i) = P.matchN(i)/P.candN(i);    
%         matchscore = reshape(cell2mat(evseq{i}.matchscore), mm*nn, 1); shufZscore = reshape(cell2mat(evseq{i}.shufZscore), mm*nn, 1);
%         err = reshape(cell2mat(evseq{i}.decodeerr), mm*nn, 1);
%     kk = candevind{i}; 
%     if ~isempty(kk)
%        T = matchscore(kk); T = T(~isnan(T)); if numel(T)>=1 P.candmeanScore(i)=mean(T); P.candmedScore(i)=median(T); P.candmeanAbsScore(i)=mean(abs(T)); P.candmedAbsScore(i)=median(abs(T));  end
%        T = shufZscore(kk); T = T(~isnan(T)); if numel(T)>=1 P.candmeanZscore(i)=mean(T); P.candmedZscore(i)=median(T);  P.candmeanAbsZscore(i)=mean(abs(T)); P.candmedAbsZscore(i)=median(abs(T)); end
%        T = err(kk); T = T(~isnan(T)); if numel(T)>=1 P.candmeanErr(i)=mean(T); P.candmedErr(i)=median(T); end
%     end
%     kk = matchevind{i}; 
%     if ~isempty(kk)
%        T = matchscore(kk); T = T(~isnan(T)); if numel(T)>=1 P.matchmeanScore(i)=mean(T); P.matchmedScore(i)=median(T); P.matchmeanAbsScore(i)=mean(abs(T)); P.matchmedAbsScore(i)=median(abs(T));  end
%        T = shufZscore(kk); T = T(~isnan(T)); if numel(T)>=1 P.matchmeanZscore(i)=mean(T); P.matchmedZscore(i)=median(T);  P.matchmeanAbsZscore(i)=mean(abs(T)); P.matchmedAbsZscore(i)=median(abs(T)); end
%        T = err(kk); T = T(~isnan(T)); if numel(T)>=1 P.matchmeanErr(i)=mean(T); P.matchmedErr(i)=median(T); end
%     end
%     kk = posevind{i}; 
%     if ~isempty(kk)
%        T = matchscore(kk); T = T(~isnan(T)); if numel(T)>=1 P.posmeanScore(i)=mean(T); P.posmedScore(i)=median(T);   end
%        T = shufZscore(kk); T = T(~isnan(T)); if numel(T)>=1 P.posmeanZscore(i)=mean(T); P.posmedZscore(i)=median(T); end   
%        T = err(kk); T = T(~isnan(T)); if numel(T)>=1 P.posmeanErr(i)=mean(T); P.posmedErr(i)=median(T); end
%     end
%     kk = negevind{i}; 
%     if ~isempty(kk)
%        T = matchscore(kk); T = T(~isnan(T)); if numel(T)>=1 P.negmeanScore(i)=mean(T); P.negmedScore(i)=median(T);   end
%        T = shufZscore(kk); T = T(~isnan(T)); if numel(T)>=1 P.negmeanZscore(i)=mean(T); P.negmedZscore(i)=median(T); end    
%        T = err(kk); T = T(~isnan(T)); if numel(T)>=1 P.negmeanErr(i)=mean(T); P.negmedErr(i)=median(T); end
%     end
% end
% %%%% Added new: mean/median steps and otehrs
% %%%%% All cand events stat: Nstep, Ngstep, PProb, EndDis, Length, MedianStep, MaxStep, MinStep, MeanStep, StepVar
% P.NstepCand = NaN*ones(1, nev); P.NVstepCand = NaN*ones(1, nev); P.medPProbCand = NaN*ones(1, nev); P.meanPProbCand = NaN*ones(1, nev);
% P.medLengthCand = NaN*ones(1, nev); P.meanLengthCand = NaN*ones(1, nev); P.medEndDisCand = NaN*ones(1, nev); P.meanEndDisCand = NaN*ones(1, nev); 
% P.meanStepCand = NaN*ones(1, nev); P.stdStepCand = NaN*ones(1, nev); P.medStepCand = NaN*ones(1, nev);
% P.maxStepCand = NaN*ones(1, nev); P.minStepCand = NaN*ones(1, nev); 
% %%%%% All match events stat: Nstep, Ngstep, PProb, EndDis, Length, MedianStep, MaxStep, MinStep, MeanStep, StepVar
% P.NstepMatch = NaN*ones(1, nev); P.NVstepMatch = NaN*ones(1, nev); P.medPProbMatch = NaN*ones(1, nev); P.meanPProbMatch = NaN*ones(1, nev); 
% P.medLengthMatch = NaN*ones(1, nev); P.meanLengthMatch = NaN*ones(1, nev); P.medEndDisMatch = NaN*ones(1, nev); P.meanEndDisMatch = NaN*ones(1, nev); 
% P.meanStepMatch = NaN*ones(1, nev); P.stdStepMatch = NaN*ones(1, nev); P.medStepMatch = NaN*ones(1, nev);
% P.maxStepMatch = NaN*ones(1, nev); P.minStepMatch = NaN*ones(1, nev);
% for (i=1:nev) 
%     %%%candidate events
%     kk = candevind{i}; prop = repackout(evseq{i}.newquant); 
%     if ~isempty(kk)
%        NN = prop.Nstep(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.NstepCand(i) = median(NN(ii)); end
%        NN = prop.Ngstep(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.NVstepCand(i) = median(NN(ii)); end %%%number of Valid steps
%        NN = prop.MedianPProb(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.medPProbCand(i) = median(NN(ii)); end
%        NN = prop.MeanPProb(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.meanPProbCand(i) = mean(NN(ii)); end
%        %disp(['Prob: ' num2str([numel(ii) median(NN(ii)) mean(NN(ii)) std(NN(ii))])]);
%        NN = prop.Length(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.medLengthCand(i) = median(NN(ii)); P.meanLengthCand(i) = mean(NN(ii)); end 
%        %disp(['Length: ' num2str([numel(ii) median(NN(ii)) mean(NN(ii)) std(NN(ii))])]);
%        NN = prop.EndDis(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.medEndDisCand(i) = median(NN(ii)); P.meanEndDisCand(i) = mean(NN(ii)); end
%        NN = prop.MeanStep(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.meanStepCand(i) = mean(NN(ii)); end
%        NN = prop.StepStd(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.stdStepCand(i) = mean(NN(ii)); end
%        NN = prop.MedianStep(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.medStepCand(i) = median(NN(ii)); end
%        NN = prop.MaxStep(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.maxStepCand(i) = median(NN(ii)); end
%        NN = prop.MinStep(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.minStepCand(i) = median(NN(ii)); end
%     end
%     %%%match events
%     kk = matchevind{i};
%     if ~isempty(kk)
%        NN = prop.Nstep(kk); P.NstepMatch(i) = median(NN); NN = prop.Ngstep(kk); P.NVstepMatch(i) = median(NN); %%%number of Valid steps
%        NN = prop.MedianPProb(kk); P.medPProbMatch(i) = median(NN); NN = prop.MeanPProb(kk); P.meanPProbMatch(i) = mean(NN);
%        NN = prop.Length(kk); P.medLengthMatch(i) = median(NN); P.meanLengthMatch(i) = mean(NN);
%        NN = prop.EndDis(kk); P.medEndDisMatch(i) = median(NN); P.meanEndDisMatch(i) = mean(NN);
%        NN = prop.MeanStep(kk); P.meanStepMatch(i) = mean(NN); NN = prop.StepStd(kk); P.stdStepMatch(i) = mean(NN);
%        NN = prop.MedianStep(kk); P.medStepMatch(i) = median(NN);
%        NN = prop.MaxStep(kk); P.maxStepMatch(i) = median(NN); NN = prop.MinStep(kk); P.minStepMatch(i) = median(NN); 
%     end
% end
% function [isrenorm, renormvarname] = determinerenormvariables(varname) %%%this is for those variables that need to re-sum again during combination
% isrenorm =0; renormvarname = [];
% if contains(lower(varname), lower('TimeRate'))
%     isrenorm = 1; renormvarname = 'tottime';
% elseif contains(lower(varname), lower('EvtRate'))
%     isrenorm = 1; renormvarname = 'totalN';
% elseif contains(lower(varname), lower('Ratio'))
%     isrenorm = 1; renormvarname = 'candN';
% end
% 
% function [evSess, evsessid] = identifysession(evTime, sessionname, startT, endT)
% evSess = []; evsessid = []; minevstart = min(evTime.start); maxevend = max(evTime.ent);
% iii = find( (startT <= minevstart) & (endT >= maxevend) );
% if (numel(iii) == 1)
%     evSess = sessionname{iii}; evsessid = iii;
% end

function evsel = filterevents(tmpevname, evName, evType, evkeyword, evkeytype, evkeynoword, evkeynotype)
evsel = 1;
if strcmp(lower(evkeyword), 'self') 
   if ~strcmp(tmpevname, evName) evsel = 0; end
elseif strcmp(lower(evkeynoword), 'self')
   if strcmp(tmpevname, evName) evsel = 0; end
else
   if (~isempty(evkeytype))||(~isempty(evkeyword)) %%%for inclusion
      if isempty(evkeytype)
         if ~contains(lower(evName), lower(evkeyword)) evsel = 0; end 
      elseif isempty(evkeyword)
         if ~strncmpi(evType, evkeytype, 3) evsel = 0; end 
      else
         if ~(contains(lower(evName), lower(evkeyword)) && strncmpi(evType, evkeytype, 3))
             evsel = 0; 
         end 
      end
   end
   if (~isempty(evkeynotype))||(~isempty(evkeynoword)) %%%for exclusion
       if isempty(evkeynotype)
          if contains(lower(evName), lower(evkeynoword)) evsel = 0; end 
       elseif isempty(evkeynoword)
          if strncmpi(evType, evkeynotype, 3) evsel = 0; end
       else
          if contains(lower(evName), lower(evkeynoword)) || strncmpi(evType, evkeynotype, 3)
              evsel = 0; 
          end   
       end
   end
end

function pp = repackout(prop)  %%prop{nep, 1}.(field)(1) to pp.(field)(1:nep) 
[mm,nn]=size(prop); pp = [];
if (mm*nn>=1)
   if (nn~=1)
       disp(['--------------------warning: unexpected dimensions encountered!']);
   end
   allf = fieldnames(prop{1,1});
   for (k = 1:numel(allf))
       ppnow = NaN*ones(1, mm*nn);
       for (i = 1:mm)
       for (j = 1:nn)
           ppnow((i-1)*nn+j) = prop{i,j}.(allf{k})(1);
       end
       end
       pp.(allf{k}) = ppnow;
   end
end

    