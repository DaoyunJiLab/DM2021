function DM_Seq_PlotSignificance_Callback
%%%For each selected group of items (templates or events) after seqeuncing/decoding, 
%%%    plot an actual paramter relative to the distribution of shuffled values 
%%%The plot is specific for the selected events (one or more)

hf = gcbf; pinfonow = getappdata(hf, 'pinfo'); datanow = getappdata(hf, 'data'); tagnow = get(gcbo, 'Tag');
plotparm = getappdata(hf, 'plotparm'); title = get(hf, 'Name');
%hspike = getappdata(hf, 'hspike'); spikeselection = getappdata(hspike, 'selection'); cellind = find(spikeselection==1);
hgroup = getappdata(hf, 'hgroup'); groupselection = getappdata(hgroup, 'selection'); 
grpind = find(groupselection == 1); groupnames = datanow.grouplist.groupname(grpind);
groupcellind = datanow.grouplist.groupindex(grpind);
disp('-----> Computing combined significance......');
ok = 1;
if (numel(grpind)<1) disp('-----------> no groups selected'); ok= 0; end
%cellind = cellind(1); %%%%for sequence disiplay: only display the  first selected template ---No
if ok
   if (plotparm.evselect == 0) %%%select sessions
      disp('-----------> shuffle significance cannot be computed: need to select the event option'); ok = 0;
   elseif ~isfield(datanow, 'data')
      ok = 0; disp('-----------> aborted: data not available for stripped databases');
   else
      %evseq = datanow.data.evseq(cellind); %%%%not necessary to use this, since values all present in pinfonow
      if ~isfield(datanow.data, 'shufevseq')
         disp('-----------> shuffle significance cannot be computed: shuffle data not exit or saved'); ok = 0;
      end
   end
end
if ok
    seqtype = []; evvar = 'eventname'; typevar = 'eventType';
    if (isfield(pinfonow.parm, 'seqtype')) seqtype = pinfonow.parm.seqtype{groupcellind{1}(1)}; end
    if contains(seqtype, 'evtitemized') 
        evvar = 'tmpID'; typevar = 'tmpType';
    end
    someparm.seqtype = seqtype; 
end
if ok %%% checkout variables to combine for significance --- add here if new categories of variables need to be computed
    [nvar, catnames, varnames, isrenorm, renormcatname, renormvarname] = locatevariablestocompute(hf, pinfonow);
    if (nvar<1)
        disp('-----------> shuffle significance cannot be computed: no seq variables selected'); ok = 0;
    end
end
if ok %%%%fruther select particular events to compute
    someparm.evkeyword = []; someparm.evkeytype = []; someparm.evkeynoword = []; someparm.evkeynotype = [];
    if ~(contains(seqtype, 'evtitemized') && isempty(find(contains(catnames, 'seq')))) %%% No futher selection of tmplates for event-itemized evenet variables, since templates already selected (performed on selected groups) 
        input = inputdlg({'Keyword (self = template event)'; 'Key type'; 'KeyNOword (self to exclude template event)'; 'KeyNOtype'}, 'Event/template selection', 4, {'self'; ''; 'first'; ''}); 
        if (~isempty(input))
           someparm.evkeyword = input{1}; someparm.evkeytype = input{2}; someparm.evkeynoword = input{3}; someparm.evkeynotype = input{4};
        else
           ok = 0;
        end
    end
end
if ok %%%now need to re-quantify shuffle sequences - not totally from scratch, still may be time comsuming
for pp = 1:numel(groupnames)
    %%%%%%%for .seq variables: pre-count how many tempaltes over how many events == nitems
    cellind = groupcellind{pp}; ncell = numel(cellind);
    nitems = 0; nshuffle = pinfonow.parm.Nshuffle{cellind(1)};
    for (ttt=1:numel(cellind))
        k = cellind(ttt); evName = pinfonow.general.(evvar){k}; evType = pinfonow.parm.(typevar){k};
        if contains(seqtype, 'evtitemized')
            chktmpevname = pinfonow.general.eventname{k};
        else
            chktmpevname = pinfonow.tmp.sessevname{k};
        end
        for (j = 1:numel(evName))
            if filterevents(chktmpevname, evName{j}, evType{j}, someparm.evkeyword, someparm.evkeytype, someparm.evkeynoword, someparm.evkeynotype)
                nitems = nitems + 1;
            end
        end
    end
    for (kk = 1:nvar) %%%for each seq variable
                    %%%%%%%%%%%% .seq .seq2D variables: calculate here
       if strcmp(catnames{kk}, 'seq') || strcmp(catnames{kk}, 'seq2D') 
          tottime = zeros(nitems, 1); val = NaN*ones(nitems,1); shufval = NaN*ones(nitems, nshuffle); 
          renormval = NaN*ones(nitems, 1); shufrenormval = NaN*ones(nitems, nshuffle); 
          nnow = 0;
          for (ttt = 1:numel(cellind))
               k = cellind(ttt); evName = pinfonow.general.(evvar){k}; evType = pinfonow.parm.(typevar){k}; 
               if contains(seqtype, 'evtitemized')
                  chktmpevname = pinfonow.general.eventname{k};
               else
                  chktmpevname = pinfonow.tmp.sessevname{k};
               end
               shufevseq = datanow.data.shufevseq{k}; %%%now shufevseq{nshuffle, nev}
               someparm.eventoption = pinfonow.parm.eventoption{k}; 
               someparm.minactivecell = pinfonow.parm.minactivecell{k}; someparm.significancelevel = pinfonow.parm.significancelevel{k}; 
               for (j = 1:numel(evName)) %%%this is ready for .seq variables, need to work on .event variables
                   if filterevents(chktmpevname, evName{j}, evType{j}, someparm.evkeyword, someparm.evkeytype, someparm.evkeynoword, someparm.evkeynotype)
                      nnow = nnow + 1; val(nnow) = pinfonow.(catnames{kk}).(varnames{kk}){k}(j); 
                      if contains(seqtype, 'Bayesian2D') 
                          tottime(nnow) =  pinfonow.seq2D.tottime{k}(j);
                      else
                          tottime(nnow) =  pinfonow.seq.tottime{k}(j);
                      end
                      if isrenorm(kk) renormval(nnow) = pinfonow.(renormcatname{kk}).(renormvarname{kk}){k}(j); end
                      if ~isnan(val(nnow))
                         if contains(seqtype, '2D') %%%if 2D Bayesian decoding 
                            [shufPnow, candevind, targevind] = find2Dquantifications(shufevseq(:,j), tottime(nnow), someparm, nshuffle); %%%Now: shufevseq{nshuffle}
                            if isrenorm(kk) shufrenormval(nnow,:) = shufPnow.(renormvarname{kk}); end %%%norm factors always in the seq2D category
                            if strcmp(catnames{kk}, 'seq2D') %%%% if 2D variabels 
                               shufval(nnow,:) = shufPnow.(varnames{kk});  
                            elseif strcmp(catnames{kk}, 'seq')
                               shufPnow = find1Dquantificationsof2D(shufevseq(:,j), shufPnow.tottime, candevind, targevind, shufPnow.totalN, shufPnow.candN, shufPnow.targN, someparm, nshuffle);
                               shufval(nnow,:) = shufPnow.(varnames{kk});  
                            end     
                         else %%%if 1D Bayesian decoding
                            shufPnow = find1Dquantifications(shufevseq(:,j), tottime(nnow), someparm, nshuffle); %%%Now: shufevseq{nshuffle}
                            shufval(nnow,:) = shufPnow.(varnames{kk}); 
                            if isrenorm(kk) shufrenormval(nnow,:) = shufPnow.(renormvarname{kk}); end
                         end
                      end
                   end
               end
          end
       elseif strcmp(catnames{kk}, 'event') 
          val = cell2mat(pinfonow.event.(varnames{kk})(cellind)); 
          if isrenorm(kk) renormval = cell2mat(pinfonow.(renormcatname{kk}).(renormvarname{kk})(cellind)); end
          %%%% re-compute shuffle event variables 
          shufval = NaN*ones(ncell, nshuffle); if isrenorm(kk) shufrenormval = NaN*ones(ncell, nshuffle); end
          someparm.eventoption = pinfonow.parm.eventoption{k}; 
          someparm.minactivecell = pinfonow.parm.minactivecell{cellind(1)}; someparm.significancelevel = pinfonow.parm.significancelevel{cellind(1)};
          for (ttt = 1:numel(cellind))
               k = cellind(ttt); shufevseq = datanow.data.shufevseq{k}; %%%now shufevseq{nshuffle, nev}
               shufPnow = find1DEventItemized_Quantifications(pinfonow, shufevseq, k, someparm, nshuffle); 
               shufval(ttt,:) = shufPnow.(varnames{kk}); if isrenorm(kk) shufrenormval(ttt,:) = shufPnow.(renormvarname{kk}); end
%                disp(varnames{kk})
%                disp(renormvarname{kk})
%                disp(shufval(ttt,:))
          end
          nitems = ncell; %%%for the plot below
       end
  %%%%%%%%%%%%% Plot variable by variable here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       [mmm,nnn] = size(shufval); 
       if mmm>1
          histplotsig(mean(val), mean(shufval), nitems, someparm, plotparm, title, [catnames{kk} '_' varnames{kk} '-mean'], groupnames{pp});
          histplotsig(median(val), median(shufval), nitems, someparm, plotparm, title, [catnames{kk} '_' varnames{kk} '-median'], groupnames{pp});
          if ~isrenorm(kk)
             histplotsig(sum(val), sum(shufval), nitems, someparm, plotparm, title, [catnames{kk} '_' varnames{kk} '-combine'], groupnames{pp});
          else
             valnow = sum(val.*renormval)/sum(renormval); shufvalnow = sum(shufval.*shufrenormval) ./ sum(shufrenormval);
             histplotsig(valnow, shufvalnow, nitems, someparm, plotparm, title, [catnames{kk} '_' varnames{kk} '-combine'], groupnames{pp});  
          end
       else %%%% if only one item = mean and sum will add up row values - do not do this
           %%% for one item, no need to renormalize
           histplotsig(val, shufval, nitems, someparm, plotparm, title, [catnames{kk} '_' varnames{kk} '-single'], groupnames{pp}); 
       end
    end
end
end
disp('*******************************');

function histplotsig(val, shufval, nval, someparm, plotparm, title, tag, groupname)
binvector = [plotparm.range(1):plotparm.bin:plotparm.range(2)]; 
hf = figure('Name', title); hax = axes('Parent', hf, 'NextPlot', 'add'); str{1} = tag;
Y = histc(shufval, binvector); 
if (plotparm.setlog == 1)
    binvector = log10(binvector);
elseif (plotparm.setlog == 2)
    Y = log10(Y);
elseif (plotparm.setlog == 3)
    binvector = log10(binvector); Y = log10(Y);
end
bar(hax, binvector, Y, 'EdgeColor', [0 0 0], 'FaceColor', [0 0 0]); 
line(hax, [val val], [0 max(Y)], 'Color', [1 0 0]);
shufval = shufval(~isnan(shufval));
mm = mean(shufval); ss = std(shufval); md = median(shufval); nn = numel(shufval);
se=NaN; sp=NaN; zz=NaN; pp=NaN;
if (ss>0) 
    zz=(val-mm)/ss;
    if (zz>=0)
        pp=1-normcdf(zz); 
    else
        pp=normcdf(zz);
    end
end
if (nn>=1) 
    se = ss/sqrt(nn); 
    sp = numel(find(shufval>=val))/nn; if zz<0 sp = numel(find(shufval<=val))/nn; end
end
str{2} = strcat('Actual=', num2str(val),';zscore=', num2str(zz), ';Z-test p=', num2str(pp), ';shuffle p=', num2str(sp),';N=', num2str(nval));
str{3} = strcat('Shuffle mean=', num2str(mm), ';se=', num2str(se), ';n=', num2str(nn), ';median=', num2str(md));
str{4} = strcat('Keyword=', someparm.evkeyword, ';keytype=', someparm.evkeytype, ';keynoword=', someparm.evkeynoword, ';keynotype=', someparm.evkeynotype);
str{5} = ['Groupname = ' groupname];
for (i = 1:numel(str))
    text('Parent', hax, 'String', str{i}, 'Interpreter', 'none', 'Color', [1 0 0], 'Units', 'normalized', 'Position', [0.05 0.96-(i-1)*0.04]);
end 

function [P, candevind, matchevind, posevind, negevind] = find1Dquantifications(evseq, tottime, seqparm, nshuffle)
%%%This functon borrowed from the original quantification function in Bayesian decoding_events.
%%%Some aspects are re-arranged: nev = nshuffle and evseq{i=nshuffle} for just one event
%%%Now change evseq{i} = evseq{nshuffle}; nev = nshuffle 
%%%%%Parameters: candidate events, significant positive/negative events
nev = nshuffle;
P.totalN = NaN*ones(1, nev); P.candN = NaN*ones(1, nev); P.posN = NaN*ones(1, nev); P.negN = NaN*ones(1, nev); P.matchN =NaN*ones(1,nev);
P.evTimeRate = NaN*ones(1, nev); P.candTimeRate = NaN*ones(1, nev); P.posTimeRate = NaN*ones(1, nev); P.negTimeRate = NaN*ones(1, nev); P.matchTimeRate = NaN*ones(1, nev);
P.candEvtRate = NaN*ones(1, nev); P.posEvtRate = NaN*ones(1, nev); P.negEvtRate = NaN*ones(1, nev);  P.matchEvtRate = NaN*ones(1, nev);
P.posRatio = NaN*ones(1, nev); P.negRatio = NaN*ones(1, nev); P.matchRatio = NaN*ones(1, nev);
candevind = cell(1,nev); matchevind=cell(1,nev); posevind=cell(1,nev); negevind=cell(1,nev);
%%%% Added new: mean/median matchscore, mean matchZscore
P.tottime = tottime*ones(1, nev); P.candmeanScore = NaN*ones(1, nev); P.candmedScore = NaN*ones(1, nev); P.candmeanZscore = NaN*ones(1, nev); P.candmedZscore = NaN*ones(1, nev);
P.candmeanAbsScore = NaN*ones(1, nev); P.candmedAbsScore = NaN*ones(1, nev); P.candmeanAbsZscore = NaN*ones(1, nev); P.candmedAbsZscore = NaN*ones(1, nev);
P.posmeanScore = NaN*ones(1, nev); P.posmedScore = NaN*ones(1, nev); P.posmeanZscore = NaN*ones(1, nev); P.posmedZscore = NaN*ones(1, nev);
P.negmeanScore = NaN*ones(1, nev); P.negmedScore = NaN*ones(1, nev); P.negmeanZscore = NaN*ones(1, nev); P.negmedZscore = NaN*ones(1, nev);
P.matchmeanScore = NaN*ones(1, nev); P.matchmedScore = NaN*ones(1, nev); P.matchmeanZscore = NaN*ones(1, nev);  P.matchmedZscore = NaN*ones(1, nev);
P.matchmeanAbsScore = NaN*ones(1, nev); P.matchmedAbsScore = NaN*ones(1, nev); P.matchmeanAbsZscore = NaN*ones(1, nev); P.matchmedAbsZscore = NaN*ones(1, nev);
%%%% Added new: mean/median decodeerr
P.candmeanErr = NaN*ones(1, nev); P.candmedErr = NaN*ones(1, nev); P.posmeanErr = NaN*ones(1, nev); P.posmedErr = NaN*ones(1, nev);
P.negmeanErr = NaN*ones(1, nev); P.negmedErr = NaN*ones(1, nev); P.matchmeanErr = NaN*ones(1, nev); P.matchmedErr = NaN*ones(1, nev);
for (i = 1:nev)
        [mm,nn] = size(evseq{i}.shufposP); shufposP = repack_single(evseq{i}.shufposP); shufnegP = repack_single(evseq{i}.shufnegP); 
        nactraw = repack_single(evseq{i}.nactiveraw); nact = repack_single(evseq{i}.nactivecell); 
        candevind{i} = find(nact>=seqparm.minactivecell); 
        posevind{i} = find( (shufposP<seqparm.significancelevel) & (nact>=seqparm.minactivecell)); 
        negevind{i} = find( (shufnegP<seqparm.significancelevel) & (nact>=seqparm.minactivecell));
        matchevind{i} = sort(union(posevind{i}, negevind{i}));
    P.totalN(i) = mm*nn; P.candN(i) = numel(find(nactraw>=seqparm.minactivecell)); 
    P.posN(i) = numel(posevind{i}); P.negN(i) = numel(negevind{i}); P.matchN(i) = numel(matchevind{i});
    P.evTimeRate(i) = P.totalN(i)/P.tottime(i); P.candTimeRate(i) = P.candN(i)/P.tottime(i); 
    P.posTimeRate(i) = P.posN(i)/P.tottime(i); P.negTimeRate(i) = P.negN(i)/P.tottime(i); P.matchTimeRate(i) = P.matchN(i)/P.tottime(i);
    P.candEvtRate(i) = P.candN(i)/P.totalN(i); P.posEvtRate(i) = P.posN(i)/P.totalN(i); P.negEvtRate(i) = P.negN(i)/P.totalN(i); P.matchEvtRate(i) = P.matchN(i)/P.totalN(i);
    P.posRatio(i) = P.posN(i)/P.candN(i); P.negRatio(i) = P.negN(i)/P.candN(i);  P.matchRatio(i) = P.matchN(i)/P.candN(i);    
        matchscore = repack_single(evseq{i}.matchscore); shufZscore = repack_single(evseq{i}.shufZscore); err = repack_single(evseq{i}.decodeerr);    
    kk = candevind{i}; 
    if ~isempty(kk)
       T = matchscore(kk); T = T(~isnan(T)); if numel(T)>=1 P.candmeanScore(i)=mean(T); P.candmedScore(i)=median(T); P.candmeanAbsScore(i)=mean(abs(T)); P.candmedAbsScore(i)=median(abs(T));  end
       T = shufZscore(kk); T = T(~isnan(T)); if numel(T)>=1 P.candmeanZscore(i)=mean(T); P.candmedZscore(i)=median(T);  P.candmeanAbsZscore(i)=mean(abs(T)); P.candmedAbsZscore(i)=median(abs(T)); end
       T = err(kk); T = T(~isnan(T)); if numel(T)>=1 P.candmeanErr(i)=mean(T); P.candmedErr(i)=median(T); end
    end
    kk = matchevind{i}; 
    if ~isempty(kk)
       T = matchscore(kk); T = T(~isnan(T)); if numel(T)>=1 P.matchmeanScore(i)=mean(T); P.matchmedScore(i)=median(T); P.matchmeanAbsScore(i)=mean(abs(T)); P.matchmedAbsScore(i)=median(abs(T));  end
       T = shufZscore(kk); T = T(~isnan(T)); if numel(T)>=1 P.matchmeanZscore(i)=mean(T); P.matchmedZscore(i)=median(T);  P.matchmeanAbsZscore(i)=mean(abs(T)); P.matchmedAbsZscore(i)=median(abs(T)); end
       T = err(kk); T = T(~isnan(T)); if numel(T)>=1 P.matchmeanErr(i)=mean(T); P.matchmedErr(i)=median(T); end
    end
    kk = posevind{i}; 
    if ~isempty(kk)
       T = matchscore(kk); T = T(~isnan(T)); if numel(T)>=1 P.posmeanScore(i)=mean(T); P.posmedScore(i)=median(T);   end
       T = shufZscore(kk); T = T(~isnan(T)); if numel(T)>=1 P.posmeanZscore(i)=mean(T); P.posmedZscore(i)=median(T); end   
       T = err(kk); T = T(~isnan(T)); if numel(T)>=1 P.posmeanErr(i)=mean(T); P.posmedErr(i)=median(T); end
    end
    kk = negevind{i}; 
    if ~isempty(kk)
       T = matchscore(kk); T = T(~isnan(T)); if numel(T)>=1 P.negmeanScore(i)=mean(T); P.negmedScore(i)=median(T);   end
       T = shufZscore(kk); T = T(~isnan(T)); if numel(T)>=1 P.negmeanZscore(i)=mean(T); P.negmedZscore(i)=median(T); end    
       T = err(kk); T = T(~isnan(T)); if numel(T)>=1 P.negmeanErr(i)=mean(T); P.negmedErr(i)=median(T); end
    end
end
%%%% Added new: mean/median steps and otehrs
%%%%% All cand events stat: Nstep, Ngstep, PProb, EndDis, Length, MedianStep, MaxStep, MinStep, MeanStep, StepVar
P.candNstep = NaN*ones(1, nev); P.candNVstep = NaN*ones(1, nev); P.candNDstep = NaN*ones(1, nev); P.candmedPProb = NaN*ones(1, nev); P.candmeanPProb = NaN*ones(1, nev);
P.candmedLength = NaN*ones(1, nev); P.candmeanLength = NaN*ones(1, nev); P.candmedEndDis = NaN*ones(1, nev); P.candmeanEndDis = NaN*ones(1, nev); 
P.candmeanStep = NaN*ones(1, nev); P.candstdStep = NaN*ones(1, nev); P.candmedStep = NaN*ones(1, nev);
P.candmaxStep = NaN*ones(1, nev); P.candminStep = NaN*ones(1, nev); 
%%%%% All match events stat: Nstep, Ngstep, PProb, EndDis, Length, MedianStep, MaxStep, MinStep, MeanStep, StepVar
P.matchNstep = NaN*ones(1, nev); P.matchNVstep = NaN*ones(1, nev); P.matchNDstep = NaN*ones(1, nev); P.matchmedPProb = NaN*ones(1, nev); P.matchmeanPProb = NaN*ones(1, nev); 
P.matchmedLength = NaN*ones(1, nev); P.matchmeanLength = NaN*ones(1, nev); P.matchmedEndDis = NaN*ones(1, nev); P.matchmeanEndDis = NaN*ones(1, nev); 
P.matchmeanStep = NaN*ones(1, nev); P.matchstdStep = NaN*ones(1, nev); P.matchmedStep = NaN*ones(1, nev);
P.matchmaxStep = NaN*ones(1, nev); P.matchminStep = NaN*ones(1, nev);
for (i=1:nev) 
    %%%candidate events
    kk = candevind{i}; prop = repackout(evseq{i}.newquant); 
    if ~isempty(kk)
       NN = prop.Nstep(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.candNstep(i) = median(NN(ii)); end
       NN = prop.NVstep(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.candNVstep(i) = median(NN(ii)); end %%%number of Valid steps
       NN = prop.NDstep(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.candNDstep(i) = median(NN(ii)); end %%%number of decoded positions
       NN = prop.MedianPProb(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.candmedPProb(i) = median(NN(ii)); end
       NN = prop.MeanPProb(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.candmeanPProb(i) = mean(NN(ii)); end
       %disp(['Prob: ' num2str([numel(ii) median(NN(ii)) mean(NN(ii)) std(NN(ii))])]);
       NN = prop.Length(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.candmedLength(i) = median(NN(ii)); P.candmeanLength(i) = mean(NN(ii)); end 
       %disp(['Length: ' num2str([numel(ii) median(NN(ii)) mean(NN(ii)) std(NN(ii))])]);
       NN = prop.EndDis(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.candmedEndDis(i) = median(NN(ii)); P.candmeanEndDis(i) = mean(NN(ii)); end
       NN = prop.MeanStep(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.candmeanStep(i) = mean(NN(ii)); end
       NN = prop.StepStd(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.candstdStep(i) = mean(NN(ii)); end
       NN = prop.MedianStep(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.candmedStep(i) = median(NN(ii)); end
       NN = prop.MaxStep(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.candmaxStep(i) = median(NN(ii)); end
       NN = prop.MinStep(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.candminStep(i) = median(NN(ii)); end
    end
    %%%match events
    kk = matchevind{i};
    if ~isempty(kk)
       NN = prop.Nstep(kk); P.matchNstep(i) = median(NN); 
       NN = prop.NVstep(kk); P.matchNVstep(i) = median(NN); %%%number of Valid steps
       NN = prop.NDstep(kk); P.matchNDstep(i) = median(NN); %%%number of decoded positions      
       NN = prop.MedianPProb(kk); P.matchmedPProb(i) = median(NN); NN = prop.MeanPProb(kk); P.matchmeanPProb(i) = mean(NN);
       NN = prop.Length(kk); P.matchmedLength(i) = median(NN); P.matchmeanLength(i) = mean(NN);
       NN = prop.EndDis(kk); P.matchmedEndDis(i) = median(NN); P.matchmeanEndDis(i) = mean(NN);
       NN = prop.MeanStep(kk); P.matchmeanStep(i) = mean(NN); NN = prop.StepStd(kk); P.matchstdStep(i) = mean(NN);
       NN = prop.MedianStep(kk); P.matchmedStep(i) = median(NN);
       NN = prop.MaxStep(kk); P.matchmaxStep(i) = median(NN); NN = prop.MinStep(kk); P.matchminStep(i) = median(NN); 
    end
end

function [P, candevind, targevind] = find2Dquantifications(evseq, tottime, seqparm, nshuffle)%%%Now: shufevseq{nshuffle}
%%%This functon borrowed from the original quantification function in Bayesian decoding_events.
%%%Some aspects are re-arranged: nev = nshuffle and evseq{i=nshuffle} for just one event
%%%Now change evseq{i} = evseq{nshuffle}; nev = nshuffle 
%%%%%Parameters: candidate events, significant positive/negative events
nev = nshuffle;
%%%%% All cand events stat: Nstep, Ngstep, PProb, EndDis, Length, MedianStep, MaxStep, MinStep, MeanStep, StepVar
P.tottime = tottime*ones(1, nev); P.totalN = NaN*ones(1, nev); P.evTimeRate = NaN*ones(1, nev);%%%evseq{nev}.(field){mm,nn}(1)
P.candN = NaN*ones(1, nev); P.candTimeRate = NaN*ones(1, nev); P.candEvtRate = NaN*ones(1, nev);
P.targN = NaN*ones(1, nev); P.targTimeRate = NaN*ones(1, nev); P.targEvtRate = NaN*ones(1, nev); P.targRatio = NaN*ones(1, nev); 
P.NstepCand = NaN*ones(1, nev); P.NVstepCand = NaN*ones(1, nev); P.NCstepCand = NaN*ones(1, nev);
P.medPProbCand = NaN*ones(1, nev); P.meanPProbCand = NaN*ones(1, nev);
P.medErrCand = NaN*ones(1, nev); P.meanErrCand = NaN*ones(1, nev); 
P.medLengthCand = NaN*ones(1, nev); P.meanLengthCand = NaN*ones(1, nev); P.medCLengthCand = NaN*ones(1, nev); P.meanCLengthCand = NaN*ones(1, nev); 
P.medEndDisCand = NaN*ones(1, nev); P.meanEndDisCand = NaN*ones(1, nev); 
P.meanStepCand = NaN*ones(1, nev); P.varStepCand = NaN*ones(1, nev); P.medStepCand = NaN*ones(1, nev);
P.maxStepCand = NaN*ones(1, nev); P.minStepCand = NaN*ones(1, nev); 
%%%%% All target events stat: Nstep, Ngstep, PProb, EndDis, Length, MedianStep, MaxStep, MinStep, MeanStep, StepVar
P.NstepTarg = NaN*ones(1, nev); P.NVstepTarg = NaN*ones(1, nev); P.NCstepTarg = NaN*ones(1, nev);
P.medPProbTarg = NaN*ones(1, nev); P.meanPProbTarg = NaN*ones(1, nev); 
P.medErrTarg = NaN*ones(1, nev); P.meanErrTarg = NaN*ones(1, nev); 
P.medLengthTarg = NaN*ones(1, nev); P.meanLengthTarg = NaN*ones(1, nev); P.medCLengthTarg = NaN*ones(1, nev); P.meanCLengthTarg = NaN*ones(1, nev);
P.medEndDisTarg = NaN*ones(1, nev); P.meanEndDisTarg = NaN*ones(1, nev); 
P.meanStepTarg = NaN*ones(1, nev); P.varStepTarg = NaN*ones(1, nev); P.medStepTarg = NaN*ones(1, nev);
P.maxStepTarg = NaN*ones(1, nev); P.minStepTarg = NaN*ones(1, nev);
targevind = cell(1, nev); candevind = cell(1, nev);
for (i = 1:nev) %for each event file
    prop2D = repackout(evseq{i}.newquant);  %%evseq{i}.newquant.(field){nep, 1}(1) to prop2D.(field)(1:nep) 
    ncell = repackoutsinglefield(evseq{i}.nactivecell); %%evseq{i}.nactivecell{nep, 1}(1) to ncell(1:nep) 
    P.totalN(i) = numel(prop2D.Nstep); 
    kk = find(ncell>=seqparm.minactivecell); P.candN(i) = numel(kk); candevind{i} = kk; %%%candidate events
    kk = find(prop2D.IsTarget); P.targN(i) = numel(kk); targevind{i} = kk; %%%target events
end
%%%%%% basic event rates
P.evTimeRate = P.totalN./P.tottime; P.candTimeRate = P.candN./P.tottime; P.candEvtRate = P.candN./P.totalN; 
P.targTimeRate = P.targN./P.tottime; P.targEvtRate = P.targN./P.totalN; P.targRatio = P.targN./P.candN;
for (i=1:nev) %%%%decoding properties
    %%%candidate events
    kk = candevind{i}; prop2D = repackout(evseq{i}.newquant); %%evseq{i}.newquant.(field){nep, 1}(1)
    NN = prop2D.Nstep(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.NstepCand(i) = median(NN(ii)); end
    NN = prop2D.NVstep(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.NVstepCand(i) = median(NN(ii)); end %%%number of Valid steps
    NN = prop2D.NCstep(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.NCstepCand(i) = median(NN(ii)); end %%%number of Continuous target steps
    NN = prop2D.MedianPProb(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.medPProbCand(i) = median(NN(ii)); end
    NN = prop2D.MeanPProb(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.meanPProbCand(i) = mean(NN(ii)); end
    NN = prop2D.MedDecode2dErr(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.medErrCand(i) = median(NN(ii)); end
    NN = prop2D.MeanDecode2dErr(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.meanErrCand(i) = mean(NN(ii)); end
    NN = prop2D.Length(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.medLengthCand(i) = median(NN(ii)); P.meanLengthCand(i) = mean(NN(ii)); end
    NN = prop2D.CLength(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.medCLengthCand(i) = median(NN(ii)); P.meanCLengthCand(i) = mean(NN(ii)); end
    NN = prop2D.EndDis(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.medEndDisCand(i) = median(NN(ii)); P.meanEndDisCand(i) = mean(NN(ii)); end
    NN = prop2D.MeanStep(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.meanStepCand(i) = mean(NN(ii)); end
    NN = prop2D.StepVar(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.varStepCand(i) = mean(NN(ii)); end
    NN = prop2D.MedianStep(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.medStepCand(i) = median(NN(ii)); end
    NN = prop2D.MaxStep(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.maxStepCand(i) = median(NN(ii)); end
    NN = prop2D.MinStep(kk); ii = find(~isnan(NN)); if ~isempty(ii) P.minStepCand(i) = median(NN(ii)); end
    %%%target events
    kk = targevind{i};
    if ~isempty(kk)
       NN = prop2D.Nstep(kk); P.NstepTarg(i) = median(NN); 
       NN = prop2D.NVstep(kk); P.NVstepTarg(i) = median(NN); NN = prop2D.NCstep(kk); P.NCstepTarg(i) = median(NN);
       NN = prop2D.MedianPProb(kk); P.medPProbTarg(i) = median(NN); NN = prop2D.MeanPProb(kk); P.meanPProbTarg(i) = mean(NN);
       NN = prop2D.MedDecode2dErr(kk); P.medErrTarg(i) = median(NN); NN = prop2D.MeanDecode2dErr(kk); P.meanErrTarg(i) = mean(NN);
       NN = prop2D.Length(kk); P.medLengthTarg(i) = median(NN); P.meanLengthTarg(i) = mean(NN);
       NN = prop2D.CLength(kk); P.medCLengthTarg(i) = median(NN); P.meanCLengthTarg(i) = mean(NN);   
       NN = prop2D.EndDis(kk); P.medEndDisTarg(i) = median(NN); P.meanEndDisTarg(i) = mean(NN);
       NN = prop2D.MeanStep(kk); P.meanStepTarg(i) = mean(NN); NN = prop2D.StepVar(kk); P.varStepTarg(i) = mean(NN);
       NN = prop2D.MedianStep(kk); P.medStepTarg(i) = median(NN);
       NN = prop2D.MaxStep(kk); P.maxStepTarg(i) = median(NN); NN = prop2D.MinStep(kk); P.minStepTarg(i) = median(NN); 
    end
end
function P = find1Dquantificationsof2D(evseq, tottime, candevind, targevind,totalN, candN, targN, seqparm, nj)
%%%%%%%%%canevind not used here: assuming all computed events are candidate events
%%%%%candidate/target events: matchscore(time-pos correlation), matchZ (corrrelative to shuffle), 
%%%%%                  significant positive (forward)/negative(reverse) matche N and ratio
P.meanScoreCand = NaN*ones(1, nj); P.medScoreCand = NaN*ones(1, nj); P.meanScoreZCand = NaN*ones(1, nj); %%% mean/med score values, Z means
P.med1dErrCand = NaN*ones(1, nj); P.mean1dErrCand = NaN*ones(1, nj); %%%nj here is the same as nev below
P.posNCand = NaN*ones(1, nj); P.posTimeRateCand = NaN*ones(1, nj); P.posEvtRateCand = NaN*ones(1, nj); P.posRatioCand = NaN*ones(1, nj);
P.negNCand = NaN*ones(1, nj); P.negTimeRateCand = NaN*ones(1, nj); P.negEvtRateCand = NaN*ones(1, nj); P.negRatioCand = NaN*ones(1, nj);
P.matchNCand = NaN*ones(1, nj); P.matchTimeRateCand = NaN*ones(1, nj); P.matchEvtRateCand = NaN*ones(1, nj); P.matchRatioCand = NaN*ones(1, nj);
P.meanScoreTarg = NaN*ones(1, nj); P.medScoreTarg = NaN*ones(1, nj); P.meanScoreZTarg = NaN*ones(1, nj); P.med1dErrTarg = NaN*ones(1, nj); P.mean1dErrTarg = NaN*ones(1, nj);
P.posNTarg = NaN*ones(1, nj); P.posTimeRateTarg = NaN*ones(1, nj); P.posEvtRateTarg = NaN*ones(1, nj); P.posRatioTarg = NaN*ones(1, nj);
P.negNTarg = NaN*ones(1, nj); P.negTimeRateTarg = NaN*ones(1, nj); P.negEvtRateTarg = NaN*ones(1, nj); P.negRatioTarg = NaN*ones(1, nj);
P.matchNTarg = NaN*ones(1, nj); P.matchTimeRateTarg = NaN*ones(1, nj); P.matchEvtRateTarg = NaN*ones(1, nj); P.matchRatioTarg = NaN*ones(1, nj);
for (j = 1:nj)    
    %%% format of decoding result: evseq{i}.shufnegP{mm, nn}(1): mm = nep, nn = 1 scaling factor
    %%%%%%% candidate events
        shufposall = repackoutsinglefield(evseq{j}.shufposP); shufnegall = repackoutsinglefield(evseq{j}.shufnegP); %%prop{nep, 1}(1) to pp(1:nep)
        nactall = repackoutsinglefield(evseq{j}.nactivecell); 
        scoreall = repackoutsinglefield(evseq{j}.matchscore); scoreZall = repackoutsinglefield(evseq{j}.shufZscore);
        errall = repackoutsinglefield(evseq{j}.decodeerr); 
        cind = candevind{j}; shufpos = shufposall(cind); shufneg = shufnegall(cind); nact = nactall(cind); score = scoreall(cind); scoreZ = scoreZall(cind); err = errall(cind);
    P.posNCand(j) = numel(find( (shufpos<seqparm.significancelevel) & (nact>=seqparm.minactivecell)) );
    P.negNCand(j) = numel(find( (shufneg<seqparm.significancelevel) & (nact>=seqparm.minactivecell)) );
    P.matchNCand(j) = P.posNCand(j) + P.negNCand(j);
    P.meanScoreCand(j) = mean(score(~isnan(score))); P.medScoreCand(j) = median(score(~isnan(score))); P.meanScoreZCand(j) = mean(scoreZ(~isnan(scoreZ)));
    P.med1dErrCand(j) = median(err(~isnan(err))); P.mean1dErrCand(j) = mean(err(~isnan(err)));
    P.posEvtRateCand(j) = P.posNCand(j)/totalN(j); P.posTimeRateCand(j) = P.posNCand(j)/tottime(j); P.posRatioCand(j) = P.posNCand(j)/candN(j);
    P.negEvtRateCand(j) = P.negNCand(j)/totalN(j); P.negTimeRateCand(j) = P.negNCand(j)/tottime(j); P.negRatioCand(j) = P.negNCand(j)/candN(j);
    P.matchEvtRateCand(j) = P.matchNCand(j)/totalN(j); P.matchTimeRateCand(j) = P.matchNCand(j)/tottime(j); P.matchRatioCand(j) = P.matchNCand(j)/candN(j);
    %%%%%%% target events    
        ind = targevind{j}; shufpos = shufposall(ind); shufneg = shufnegall(ind); nact = nactall(ind); score = scoreall(ind); scoreZ = scoreZall(ind); err = errall(ind);
    P.posNTarg(j) = numel(find( (shufpos<seqparm.significancelevel) & (nact>=seqparm.minactivecell)) );
    P.negNTarg(j) = numel(find( (shufneg<seqparm.significancelevel) & (nact>=seqparm.minactivecell)) ); 
    P.matchNTarg(j) = P.posNTarg(j) + P.negNTarg(j);
    P.meanScoreTarg(j) = mean(score(~isnan(score))); P.medScoreTarg(j) = median(score(~isnan(score))); P.meanScoreZTarg(j) = mean(scoreZ(~isnan(scoreZ)));
    P.med1dErrTarg(j) = median(err(~isnan(err))); P.mean1dErrTarg(j) = mean(err(~isnan(err)));
    P.posEvtRateTarg(j) = P.posNTarg(j)/totalN(j); P.posTimeRateTarg(j) = P.posNTarg(j)/tottime(j); P.posRatioTarg(j) = P.posNTarg(j)/targN(j);
    P.negEvtRateTarg(j) = P.negNTarg(j)/totalN(j); P.negTimeRateTarg(j) = P.negNTarg(j)/tottime(j); P.negRatioTarg(j) = P.negNTarg(j)/targN(j);
    P.matchEvtRateTarg(j) = P.matchNTarg(j)/totalN(j); P.matchTimeRateTarg(j) = P.matchNTarg(j)/tottime(j); P.matchRatioTarg(j) = P.matchNTarg(j)/targN(j);    
end

function [evSess, evsessid] = identifysession(evTime, sessionname, startT, endT)
evSess = []; evsessid = []; minevstart = min(evTime.start); maxevend = max(evTime.ent);
iii = find( (startT <= minevstart) & (endT >= maxevend) );
if (numel(iii) == 1)
    evSess = sessionname{iii}; evsessid = iii;
end

function evsel = filterevents(tmpevname, evName, evType, evkeyword, evkeytype, evkeynoword, evkeynotype)
evsel = 1;
if strcmp(lower(evkeyword), 'self') 
   if ~strcmp(evName, tmpevname) evsel = 0; end
elseif strcmp(lower(evkeynoword), 'self')
   if strcmp(evName, tmpevname) evsel = 0; end
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
%disp([evName '  ' tmpevname '  ' num2str(evsel)]);
% function pp = repackout(prop)  %%prop{nep, 1}.(field)(1) to pp.(field)(1:nep) 
% [mm,nn]=size(prop); pp = [];
% if (mm*nn>=1)
%    if (nn~=1)
%        disp(['--------------------warning: unexpected dimensions encountered!']);
%    end
%    allf = fieldnames(prop{1,1});
%    for (k = 1:numel(allf))
%        ppnow = NaN*ones(1, mm*nn);
%        for (i = 1:mm)
%        for (j = 1:nn)
%            ppnow((i-1)*nn+j) = prop{i,j}.(allf{k})(1);
%        end
%        end
%        pp.(allf{k}) = ppnow;
%    end
% end
function pp = repackoutsinglefield(prop)  %%prop{nep, 1}(1) to pp(1:nep) 
[mm,nn]=size(prop); pp = NaN*ones(1, mm*nn);
if (nn~=1)
    disp(['--------------------warning: unexpected dimensions encountered!']);
end
for (i = 1:mm)
for (j = 1:nn)
    pp((i-1)*nn+j) = prop{i,j}(1);
end
end
function pp = repackout(prop)  %%prop{nep, 1}.(field)(1 or ntarget) to pp.(field)(1:nep*ntarget) 
[mm,nn]=size(prop); pp = [];
if (mm*nn>=1)
   if (nn~=1)
       disp(['--------------------warning: unexpected dimensions encountered!']);
   end
   allf = fieldnames(prop{1,1});
   for (k = 1:numel(allf))
       ppnow = cell(1, mm*nn); %ppnow = NaN*ones(1, mm*nn);
       for (i = 1:mm)
       for (j = 1:nn)
           %ppnow((i-1)*nn+j) = prop{i,j}.(allf{k})(1);
           ppnow{(i-1)*nn+j} = prop{i,j}.(allf{k});
       end
       end
       pp.(allf{k}) = cell2mat(ppnow); %pp.(allf{k}) = ppnow;
   end
end
function pp = repack_single(Fnow) %%[prop.(singlefield)]{nlap, 1}.(1 or ntarget) to pp(nep*ntarget) 
[mm,nn]=size(Fnow); ppnow = cell(1, mm*nn);
if (nn~=1)
    disp(['--------------------warning: unexpected dimensions encountered!']);
end
for (i = 1:mm)
    for j = 1:nn
        ppnow{(i-1)*nn+j} = Fnow{i,j};
    end
end
pp = cell2mat(ppnow);

function p = find1DEventItemized_Quantifications(pinfonow, evseq, evid, someparm, nshuffle)
%%%%compute event variables for one event --- these are new from the original template-itemized database
%%%%Here posN/negN/matchN - defined as matching to any of the available templates on the same day
%%%%%%%%%%%% Output P.(varname)[1ev, nshuffle] 
%%%%evseq{oneev-evid}{nshuffle, ntmp:}{mm,nn}.(varnames){1}(nlap)...
p = []; fnames = fieldnames(pinfonow.event);
for i=1:numel(fnames) 
    currentval = pinfonow.event.(fnames{i}){evid}; 
    p.(fnames{i}) = currentval*ones(1, nshuffle);
end
%%%% compute quant values for the each shuffle copy
totaltime = p.totaltime(1); totalN = p.totalN(1); 
for (i = 1:nshuffle)
    Q = recomputequantvalues(evseq(i, :), someparm, totaltime, totalN);
    fnames = fieldnames(Q); 
    for k=1:numel(fnames) p.(fnames{k})(i) = Q.(fnames{k}); end
end
%disp(p.posN)
%%%% compute sig (Z,P) variable values for the each shuffle copy
[p.posShuffleZ,p.posShuffleP] = findsignow(p.posN, nshuffle); 
[p.negShuffleZ,p.negShuffleP] = findsignow(p.negN, nshuffle); 
[p.matchShuffleZ,p.matchShuffleP] = findsignow(p.matchN, nshuffle); 
function [Z,P] = findsignow(shufN, nshuffle)
Z = NaN*ones(1, nshuffle); P = NaN*ones(1, nshuffle);
mm = mean(shufN); ss = std(shufN);
for (i = 1:nshuffle)
    P(i) = numel(find(shufN>=shufN(i)))/nshuffle;
end
if ss~=0 
    for (i=1:nshuffle)
        Z(i)=(shufN(i)-mm)/ss;
    end
end
function p = recomputequantvalues(evseq, someparm, totaltime, totalN)
%for (n = 1:numel(evseq)) %%%number of templates
     [posN, negN, matchN, candN] = findmatchnumber(evseq, someparm, totalN);
     p.candTimeRate = candN/totaltime; p.candEvtRate = candN/totalN;
     p.posN = posN; p.posTimeRate = posN/totaltime; p.posEvtRate = posN/totalN; p.posRatio = posN/candN;
     p.negN = negN; p.negTimeRate = negN/totaltime; p.negEvtRate = negN/totalN; p.negRatio = negN/candN;
     p.matchN = matchN; p.matchTimeRate = matchN/totaltime; p.matchEvtRate = matchN/totalN; p.matchRatio = matchN/candN;
%end
function [posN, negN, matchN, candN] = findmatchnumber(evseq, someparm, evN)
%%%% For single event: eventname with evN laps; matching to multuple templates tmpid
%%%% Output posN/negN/matchN: count a lap if it matches to any of the available templates on the same day!
posN = NaN; negN = NaN; matchN = NaN; candN = NaN; ntmp = numel(evseq);
shufposP = NaN*ones(ntmp, evN); shufnegP = NaN*ones(ntmp, evN); nactnum = NaN*ones(ntmp,evN); % nactnum = NaN*ones(evN,1);
slevel = someparm.significancelevel; candthre = someparm.minactivecell; 
if ntmp >= 1
   posPname = 'shufposP'; negPname = 'shufnegP';
   if strcmp(someparm.seqtype, 'SeqMatching')
      posPname = 'posmatchprob'; negPname = 'negmatchprob';
   end
   if contains(someparm.eventoption, 'FreeSeq')  %% if this is free squencing
       seqstart = cell(1, ntmp); seqend = cell(1, ntmp); shufposP = cell(1, ntmp); shufnegP = cell(1, ntmp);
       candNraw = cell(1, ntmp); 
       for (k  = 1:ntmp)
           seqstart{k} = repack_single_freeseq(evseq{k}.seqstart); seqend{k} = repack_single_freeseq(evseq{k}.seqend);
           shufposP{k} = repack_single_freeseq(evseq{k}.shufposP); shufnegP{k} = repack_single_freeseq(evseq{k}.shufnegP);
           nactraw = repack_single_freeseq(evseq{k}.nactiveraw); candNraw = numel(find(nactraw>=someparm.minactivecell));
       end
       [candN, posN, negN, matchN] = findmatchnumbernow_freeseq(candNraw, shufposP, shufnegP, seqstart, seqend); 
   else
      for (k = 1:ntmp)
          shufposP(k,:) = cell2mat(evseq{k}.(posPname));
          shufnegP(k,:) = cell2mat(evseq{k}.(negPname));
          %nactnum = cell2mat(evseq{k}.nactivecell);
          nactnum(k,:) = cell2mat(evseq{k}.nactivecell);
      end
      [posN, negN, matchN, candN] = findmatchnumbernow(nactnum, shufposP, shufnegP, slevel, candthre);
   end
end
function [posN, negN, matchN, candN] = findmatchnumbernow(nactnum, shufposP, shufnegP, slevel, candthre)
%%%%%%% This needs to be careful: for each tmp, it has to be both nactnum>candthre and P < slevel
%%%% shufposP[ntmp, nlap]
% [mm,nn] = size(nactnum);
% if nn == 1
%    nactnum = nactnum';
% elseif (mm>1) && (nn>1)
%    nactnum = min(nactnum);
% end
% nactnum = min(nactnum); shufposP = min(shufposP); shufnegP = min(shufnegP); %%%%here is the key change: a seq is taken if matches to any of the templates
% posN = numel(find( (shufposP<slevel) & (nactnum>=candthre) )); 
% negN = numel(find( (shufnegP<slevel) & (nactnum>=candthre) )); 
% matchN = numel(find( ((shufposP<slevel) | (shufnegP<slevel)) & (nactnum>=candthre) ));
% candN = numel(find(nactnum>candthre));
istarget = (nactnum>=candthre); istarget = max(istarget); candN = numel(find(istarget==1));
istarget = (nactnum>=candthre) & (shufposP<slevel); istarget = max(istarget); posN = numel(find(istarget==1));
istarget = (nactnum>=candthre) & (shufnegP<slevel); istarget = max(istarget); negN = numel(find(istarget==1));
istarget = (nactnum>=candthre) & ((shufnegP<slevel)|(shufposP<slevel)); istarget = max(istarget); matchN = numel(find(istarget==1));

function [nvar, catnames, varnames, isrenorm, renormcatname, renormvarname] = locatevariablestocompute(hf, pinfonow)
%%%% find variables to combine for significance -- add more category of variables here
 nvar1 = 0; catname1 = []; varname1 = [];
 nf = find(strcmp(fieldnames(pinfonow), 'seq'));
 if numel(nf)== 1
    hfield = getappdata(hf, 'hfield'); fieldselection = getappdata(hfield(nf), 'selection');
    fieldnamesnow = getappdata(hfield(nf), 'textmsg'); 
    varindex = find( fieldselection == 1); varname1 = fieldnamesnow(varindex)';
    nvar1 = numel(varindex); catname1 = cell(1, nvar1); for (i=1:nvar1) catname1{i} = 'seq'; end 
 end
 nvar2 = 0; catname2 = []; varname2 = [];
 nf = find(strcmp(fieldnames(pinfonow), 'event'));
 if numel(nf)== 1
    hfield = getappdata(hf, 'hfield'); fieldselection = getappdata(hfield(nf), 'selection');
    fieldnamesnow = getappdata(hfield(nf), 'textmsg'); 
    varindex = find( fieldselection == 1); varname2 = fieldnamesnow(varindex)';
    nvar2 = numel(varindex); catname2 = cell(1, nvar2); for (i=1:nvar2) catname2{i} = 'event'; end 
 end
 nvar3 = 0; catname3 = []; varname3 = [];
 nf = find(strcmp(fieldnames(pinfonow), 'seq2D'));
 if numel(nf)== 1
    hfield = getappdata(hf, 'hfield'); fieldselection = getappdata(hfield(nf), 'selection');
    fieldnamesnow = getappdata(hfield(nf), 'textmsg'); 
    varindex = find( fieldselection == 1); varname3 = fieldnamesnow(varindex)';
    nvar3 = numel(varindex); catname3 = cell(1, nvar3); for (i=1:nvar3) catname3{i} = 'seq2D'; end 
 end
 nvar = nvar1+nvar2+nvar3; catnames = [catname1 catname2 catname3]; varnames = [varname1 varname2 varname3];
[isrenorm, renormcatname, renormvarname] = determinerenormvariables(pinfonow, nvar, catnames, varnames);
function [isrenorm, renormcatname, renormvarname] = determinerenormvariables(pinfonow, nvar, catname, varname) %%%this is for those variables that need to re-sum again during combination
isrenorm = zeros(1, nvar); renormcatname = cell(1, nvar); renormvarname = cell(1, nvar);
for (i = 1:nvar)
    if contains(lower(varname{i}), lower('TimeRate'))
       isrenorm(i) = 1; 
       if strcmp(catname{i}, 'seq') || strcmp(catname{i}, 'seq2D')
           if isfield(pinfonow.(catname{i}), 'tottime')
               renormcatname{i} = catname{i}; renormvarname{i} = 'tottime';
           elseif isfield(pinfonow, 'seq2D') && isfield(pinfonow.seq2D, 'tottime')
               renormcatname{i} = 'seq2D'; renormvarname{i} = 'tottime';
           end
       elseif strcmp(catname{i}, 'event')
           renormcatname{i} = catname{i}; renormvarname{i} = 'totaltime';
       end
    elseif contains(lower(varname{i}), lower('EvtRate'))
       isrenorm(i) = 1; renormcatname{i} = catname{i}; renormvarname{i} = 'totalN';
       if isfield(pinfonow, 'seq2D') && isfield(pinfonow.seq2D, 'totalN')
           renormcatname{i} = 'seq2D'; renormvarname{i} = 'totalN';
       end
    elseif contains(lower(varname{i}), lower('Ratio'))
       isrenorm(i) = 1; renormcatname{i} = catname{i}; renormvarname{i} = 'candN';
       if isfield(pinfonow, 'seq2D') && isfield(pinfonow.seq2D, 'candN')
           renormcatname{i} = 'seq2D'; renormvarname{i} = 'candN';
       end
    end
end
 
function pp = repack_single_freeseq(Fnow) %%[prop.(singlefield)]{nlap, 1}.(1 or ntarget) to pp(nep*ntarget) 
[mm,nn]=size(Fnow); ppnow = cell(1, mm*nn);
if (nn~=1)
    disp(['--------------------warning: unexpected dimensions encountered!']);
end
for (i = 1:mm)
    for j = 1:nn
        ppnow{(i-1)*nn+j} = Fnow{i,j};
    end
end
pp = cell2mat(ppnow);    

function [candN, posN, negN, matchN] = findmatchnumbernow_freeseq(oldcand, posP, negP, sstart, send)
%%%  For FreeSequences: nev varies from tmp to tmp, and only contains idnetified (sig) events
%slevel = seqparm.significancelevel; candthre = seqparm.minactivecell;
candN = mean(oldcand);
posP = cell2mat(posP); negP = cell2mat(negP); sstart = cell2mat(sstart); send = cell2mat(send);
posN = resolveoverlaptime_pvaltakeover(posP, sstart, send);
negN = resolveoverlaptime_pvaltakeover(negP, sstart, send);
nv = numel(posP); aa = reshape(posP, [1, nv]); bb = reshape(negP, [1, nv]); PP = min([aa;bb]);
matchN = resolveoverlaptime_pvaltakeover(PP, sstart, send);
function NN = resolveoverlaptime_lengtakeover(score, sstart, send)
ist = false(1, numel(score)); %%%set as if a target sequence
leng = send-sstart+1;
for (i = 1:numel(score))
    iii = findoverlap(sstart(ist), send(ist), sstart(i), send(i)); 
    if isempty(iii)
        ist(i) = 1;
    else
        kkk = find(ist);
        if leng(i)>=max(leng(ist(kkk(iii))))&& score(i)<min(score(ist(kkk(iii))))
            ist(i) = 1; ist(kkk(iii)) = zeros(1, numel(iii));
        end
    end
end
NN = numel(find(ist));
function NN = resolveoverlaptime_pvaltakeover(score, sstart, send)
ist = false(1, numel(score)); %%%set as if a target sequence
for (i = 1:numel(score))
    iii = findoverlap(sstart(ist), send(ist), sstart(i), send(i)); 
    if isempty(iii)
        ist(i) = 1;
    else
        kkk = find(ist);
        if score(i)<min(score(ist(kkk(iii))))
            ist(i) = 1; ist(kkk(iii)) = zeros(1, numel(iii));
        end
    end
end
NN = numel(find(ist)); 
function ind = findoverlap(newstart, newend, sigstart, sigend)
ind = [];
iii = find( (newstart<sigstart) & (newend>sigstart) ); ind = union(ind, iii);
iii = find( (newstart<sigend) & (newend>sigend) ); ind = union(ind, iii);
iii = find( (newstart>=sigstart) & (newend<=sigend) ); ind = union(ind, iii); 
 
    