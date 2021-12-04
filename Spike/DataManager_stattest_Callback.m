function DataManager_stattest_Callback
%%Plot the selected results
hf = gcbf; dbtype = getappdata(hf, 'dbtype'); celllistvar = getappdata(hf, 'celllistvar');
plotparm = getappdata(hf, 'plotparm'); vv = plotparm.showdetail; %if vv=1, shoe details of computation
%%%%%%%added options for all types of databases
pinfostr = 'pinfo'; datastr = 'data'; catok = 'general'; varok = 'parmfile';
switch dbtype
    case '.eegdb'
        pinfostr = 'eeg'; datastr = 'eegdata'; catok = 'general'; varok = 'eegfile'; %%also determine a variable (one in pinfo.general) to identify cells
    case '.behavdb'
        pinfostr = 'behav'; datastr = 'bhdata'; catok = 'general'; varok = 'sessID'; 
end
pinfo = getappdata(hf, pinfostr); data = getappdata(hf, datastr); 
hspike = getappdata(hf, 'hspike'); hgroup = getappdata(hf, 'hgroup'); groupselection = getappdata(hgroup, 'selection');
hfield = getappdata(hf, 'hfield');
valuelistpos = getappdata(hf, 'valuelistpos');
%plotparm = getappdata(hf, 'plotparm');

%determine selected groups and selected variables
groupselectionindex = find( groupselection == 1); %index for selected groups
ngroup = numel(groupselectionindex);
fieldtitle = fieldnames(pinfo); nfield = numel(fieldtitle); 
val = []; groupname = []; fieldname = []; %all values for plot are assigned here
for (k = 1:ngroup)
    spikeselectionindex = data.grouplist.groupindex{groupselectionindex(k)}; groupname{k} = data.grouplist.groupname{groupselectionindex(k)};
    nsubfield = 0;
    for (i = 1:nfield) %search selected sub-fields in each big field
         subfield{i} = fieldnames(pinfo.(fieldtitle{i}));
         fieldselection = getappdata(hfield(i), 'selection');
         fieldselectindex = find( fieldselection == 1);
         for (j = 1:numel(fieldselectindex) ) %do histogram for each selected subfields
             nsubfield = nsubfield + 1; 
             val{k,nsubfield} = pinfo.(fieldtitle{i}).(subfield{i}{fieldselectindex(j)})(spikeselectionindex);
             fieldname{nsubfield} = subfield{i}{fieldselectindex(j)};
         end
    end
end
tagmark = get(gcbo, 'Tag');
if (isempty(groupname) | isempty(fieldname))
    disp('----------> no groups or variables selected');
else
    if ( strcmp(tagmark, 'anova2') | strcmp(tagmark, 'friedman') )
        comtype = []; stattest(val, groupname, fieldname, comtype, tagmark, pinfo);
    else
        str = {'Between groups'; 'Between variables'; 'Individuals'};
        [ss,ok] = listdlg('PromptString', 'Select a comparison type', 'SelectionMode', 'single', 'ListString', str); 
        if (ok)
            comtype = str{ss};
            stattest(val, groupname, fieldname, comtype, tagmark, pinfo, vv);
        end
    end
end
disp('************************');

function stattest(val, groupname, fieldname, comtype, tagmark, pinfo, vv)
%%val{group, field}{values}
[ngroup, nfield] = size(val);
if ( strcmp(tagmark, 'anova2') | strcmp(tagmark, 'friedman') )
    [valnow, ok] = checknumeric(val, comtype);
    if (ok) testnow(valnow, groupname, fieldname, tagmark, comtype, pinfo); end
else
    if (strcmp(comtype, 'Between groups'))
        for (j = 1:nfield)
            valnow = [];
            for (i = 1:ngroup)
                valnow{i} = val{i,j};
            end
            [valnow,ok] = checknumeric(valnow, comtype);
            if (ok) testnow(valnow, fieldname{j}, groupname, tagmark, comtype, pinfo, vv); end
        end
    elseif (strcmp(comtype, 'Between variables'))
        for (i = 1:ngroup)
            valnow = [];
            for (j = 1:nfield)
                valnow{j} = val{i,j};
            end
            [valnow,ok] = checknumeric(valnow, comtype);
            if (ok) testnow(valnow, groupname{i}, fieldname, tagmark, comtype, pinfo, vv); end
        end
    elseif (strcmp(comtype, 'Individuals'))
        for (i = 1:ngroup)
        for (j = 1:nfield)
            valnow{1} = val{i,j}; groupnow{1} = groupname{i};
            [valnow,ok] = checknumeric(valnow, comtype);
            if (ok) testnow(valnow, fieldname{j}, groupnow, tagmark, comtype, pinfo, vv); end
        end
        end
    end
end

function [val,ok] = checknumeric(val, comtype)
[aa,nn] = size(val); ok = 1;
if (~isempty(comtype)) && (((aa == 1)&&(nn > 1)) || ((aa > 1)&&(nn == 1))) %if do not do anova2 or val array is singular
    for (i = 1:numel(val)) %group index
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
    end
else
    for (i = 1:aa)
        for (j = 1:nn)
            for (k = numel(val{i,j}))
                if (~isnumeric(val{i,j}{k}))
                    ok = 0; disp('-------> At least one variable is not numeric! Aborted');
                end 
            end
        end
    end
end

function testnow(val, plotname, groupname, tagmark, comtype, pinfo, vv)
%val{n} or val{m,n}
%text results for plot
if vv
    showok = 'on'; %%%show table for anova etc.
else
    showok = 'off';
end
for (i = 1:numel(val)) %%%sometimes val items could be complex numbers
    if iscell(val{i}) val{i} = cell2mat(val{i}); end
    val{i} = real(val{i});
end
tit{1} = 'Comparison'; tit{2} = 'pvalue'; tit{3} = 'statistic'; tit{4} ='DF/U(for ranksum)';
ncom = 0; textmsg{1} = []; textmsg{2} = []; textmsg{3} = []; textmsg{4} = [];
%%%now compute the statistics
if (strcmp(tagmark, 'ttest'))
    n = numel(val);
    if (~strcmp(comtype, 'Individuals'))
       for (i = 1:n)
       for (j = i+1:n)
           ncom = ncom + 1; textmsg{1}{ncom} = strcat(groupname{i},'_',groupname{j});
           [h,p,ci,stats] = ttest2(val{i}, val{j});
           textmsg{2}{ncom} = num2str(p, '%10.5e'); textmsg{3}{ncom} = num2str(stats.tstat); textmsg{4}{ncom} = num2str(stats.df);
       end
       end
    else
       for (i = 1:n)
           ncom = ncom + 1; textmsg{1}{ncom} = strcat(groupname{i});
           [h,p,ci,stats] = ttest(val{i});
           textmsg{2}{ncom} = num2str(p, '%10.5e'); textmsg{3}{ncom} = num2str(stats.tstat); textmsg{4}{ncom} = num2str(stats.df);
       end
    end
elseif (strcmp(tagmark, 'pairedt'))
    n = numel(val);
    for (i = 1:n)
       for (j = i+1:n)
           if (numel(val{i}) == numel(val{j}))
              ncom = ncom + 1; textmsg{1}{ncom} = strcat(groupname{i},'_',groupname{j});
              [h,p,ci,stats] = ttest(val{i}, val{j});
              textmsg{2}{ncom} = num2str(p, '%10.5e'); textmsg{3}{ncom} = num2str(stats.tstat); textmsg{4}{ncom} = num2str(stats.df);
           else
              msgbox('----------> item numbers do not match! no paired t-test performed');
           end
       end
    end
elseif (strcmp(tagmark, 'ranksum'))
    n = numel(val);
    for (i = 1:n)
       for (j = i+1:n)
           ncom = ncom + 1; textmsg{1}{ncom} = strcat(groupname{i},'_',groupname{j});
           vv1 = val{i}(~isnan(val{i})); vv2 = val{j}(~isnan(val{j})); 
           [p,h,stats] = ranksum(vv1, vv2);
           textmsg{2}{ncom} = num2str(p, '%10.5e'); textmsg{3}{ncom} = ['W(N1,N2)=',num2str(stats.ranksum), '(', num2str(numel(vv1)), ' ', num2str(numel(vv2)), ')']; 
           textmsg{4}{ncom} = ['U=', num2str(stats.ranksum-(numel(vv1)*(numel(vv2)+1)/2))];
       end
    end
elseif (strcmp(tagmark, 'signtest'))
    n = numel(val);
    for (i = 1:n)
       for (j = i+1:n)
           if (numel(val{i}) == numel(val{j}))
               ncom = ncom + 1; textmsg{1}{ncom} = strcat(groupname{i},'_',groupname{j});
               vv1 = val{i}(~isnan(val{i})); vv2 = val{j}(~isnan(val{j})); 
               [p,h,stats] = signtest(vv1, vv2);
               %[p,h,stats] = signtest(val{i}, val{j});
               textmsg{2}{ncom} = num2str(p, '%10.5e'); textmsg{3}{ncom} = num2str(stats.sign); textmsg{4}{ncom} = NaN;
           else
               msgbox('----------> item numbers do not match! no sign test performed');
           end
       end
    end
elseif (strcmp(tagmark, 'anova1'))
    n = numel(val);
    if (n>1)
        valnow = []; grp = []; nel = []; for (i = 1:n) nel(i) = numel(val{i}); end
        valnow = zeros(1,sum(nel)); grp = cell(1,sum(nel));
        for (i = 1:n)
             valnow((sum(nel(1:i))-nel(i)+1) : (sum(nel(1:i)))) = val{i}; 
             for (j = (sum(nel(1:i))-nel(i)+1) : (sum(nel(1:i))) ) 
                  grp{j} = groupname{i};
             end
        end
        [p, table, stats] = anova1(valnow, grp, showok); 
        ncom = ncom + 1; textmsg{1}{ncom} = strcat(groupname{1}, '_', groupname{2}, '_...');
        textmsg{2}{ncom} = num2str(p, '%10.5e'); textmsg{3}{ncom} = num2str(table{2,5}); textmsg{4}{ncom} = [num2str(table{2,3}), ',', num2str(table{3,3})];
    else %if only one group and one variable selected, treat multiple values of a cell's assignment as multiple variables 
        valok = val{1}; ncell = numel(valok); nx = zeros(1, ncell);
        for (i = 1:ncell) nx(i) = numel(valok{i}); end
        nn = min(nx); mx = [];
        if (nn > 1)
            valnow = []; grp = []; nel = []; for (i = 1:nn) nel(i) = ncell; end
            valnow = zeros(1,sum(nel)); grp = cell(1,sum(nel));
            for (i = 1:nn) %i-th values of a cell's variable
                vok = zeros(1, ncell);
                for (tt = 1:ncell) vok(tt) = valok{tt}(i); end
                valnow((sum(nel(1:i))-nel(i)+1) : (sum(nel(1:i)))) = vok; 
                for (j = (sum(nel(1:i))-nel(i)+1) : (sum(nel(1:i))) ) 
                    grp{j} = strcat(groupname{1}, '_', num2str(i));
                end
            end
            [p, table, stats] = anova1(valnow, grp, showok); 
            ncom = ncom + 1; textmsg{1}{ncom} = strcat(groupname{1}, 'p1_', groupname{1}, 'p2_...');
            textmsg{2}{ncom} = num2str(p, '%10.5e'); textmsg{3}{ncom} = num2str(table{2,5}); textmsg{4}{ncom} = [num2str(table{2,3}), ',', num2str(table{3,3})];
        else
           msgbox('---------> more than 1 set of variables are needed; ANOVA1 not performed.');
        end
    end
elseif (strcmp(tagmark, 'anova2'))
    [m,n] = size(val);
    if (m>1)&(n>1)
        valnow = []; grp{1} = []; grp{2} = []; nel = 0; 
        for (i = 1:m)
            for (j = 1:n)
                for (k = 1:numel(val{i,j}))
                    for (tt = 1:numel(val{i,j}{k}));
                        nel = nel + 1; valnow(nel) = val{i,j}{k}(tt);
                        grp{1}{nel} = plotname{i}; grp{2}{nel} = groupname{j};
                    end
                end
            end
        end
        [p, table, stats] = anovan(valnow, grp, 'display', showok); 
        ncom = ncom + 1; textmsg{1}{ncom} = strcat(plotname{1}, '_', plotname{2}, '_...');
        textmsg{2}{ncom} = num2str(p(1), '%10.5e'); textmsg{3}{ncom} = num2str(table{2,6}); textmsg{4}{ncom} = [num2str(table{2,3}), ',', num2str(table{4,3})];
        ncom = ncom + 1; textmsg{1}{ncom} = strcat(groupname{1}, '_', groupname{2}, '_...');
        textmsg{2}{ncom} = num2str(p(2), '%10.5e'); textmsg{3}{ncom} = num2str(table{3,6}); textmsg{4}{ncom} = [num2str(table{3,3}), ',', num2str(table{4,3})];
        [p, table, stats] = anovan(valnow, grp, 'model', 'interaction', 'display', showok);
        ncom = ncom + 1; textmsg{1}{ncom} = strcat(plotname{1}, '*', plotname{2}, '*...'); 
        textmsg{2}{ncom} = num2str(p(1), '%10.5e'); textmsg{3}{ncom} = num2str(table{2,6}); textmsg{4}{ncom} = [num2str(table{2,3}), ',', num2str(table{5,3})];
        ncom = ncom + 1; textmsg{1}{ncom} = strcat(groupname{1}, '*', groupname{2}, '*...');
        textmsg{2}{ncom} = num2str(p(2), '%10.5e'); textmsg{3}{ncom} = num2str(table{3,6}); textmsg{4}{ncom} = [num2str(table{3,3}), ',', num2str(table{5,3})];
        ncom = ncom + 1; textmsg{1}{ncom} = strcat(plotname{1}, '*', groupname{1}, '*...');
        textmsg{2}{ncom} = num2str(p(3), '%10.5e'); textmsg{3}{ncom} = num2str(table{4,6}); textmsg{4}{ncom} = [num2str(table{4,3}), ',', num2str(table{5,3})];
    elseif (m+n>2) %%%this allow an ANOVA2 test using number of data points within a cell's variable values as a second dimension
                   %%%usually this requires the number of data points for each cell's variable values is fixed;
                   %%% if not, here provide an option to use a reference variable as the second dimension
        valnow = []; grp{1} = []; grp{2} = []; nel = 0; 
        if (m == 1) valname = groupname; end
        if (n == 1) valname = plotname; end
        
        okk = 1; button = questdlg('Use a work variable as a reference?');
        if (strncmpi(button, 'yes', 2))
            worksubfield = fieldnames(pinfo.work);
            pp = {'Type a work variable'; 'Min value'; 'Max value'};
            def = {worksubfield{numel(worksubfield)}; '-100'; '100'};
            III=inputdlg(pp, 'Reference variable input', 3, def, 'on');
            if (~isempty(III))
                refvar = III{1}; minv = str2num(III{2}); maxv = str2num(III{3});
                for (i = 1:numel(val)) %group or variable index
                for (j = 1:numel(val{i})) %cell index
                    for (tt = 1:numel(val{i}{j})) %value index
                        if isnumeric(pinfo.work.(refvar){j}) && (tt<numel(pinfo.work.(refvar){j}))
                            refnow = pinfo.work.(refvar){j}(tt);
                            if (refnow >= minv) && (refnow <= maxv)
                                nel = nel + 1; valnow(nel) = val{i}{j}(tt);
                                grp{1}{nel} = valname{i}; grp{2}{nel} = strcat('p', num2str(refnow));
                            end
                        end
                    end
                end
                end
            else
                okk = 0;
            end
        elseif (strncmpi(button, 'no', 2))
            for (i = 1:numel(val)) %group or variable index
            for (j = 1:numel(val{i})) %cell index
                nt = 1;
                if (numel(val{i}{j})>1000) %if too many points, limited to 100 by desampling
                    nt = round(numel(val{i}{j})/100);
                end
                for (tt = 1:nt:numel(val{i}{j})) %value index
                     nel = nel + 1; valnow(nel) = val{i}{j}(tt);
                     grp{1}{nel} = valname{i}; grp{2}{nel} = strcat('p', num2str(tt));
                end
            end
            end
        else 
            okk = 0;
        end
        if okk | ~isempty(valnow)
            [p, table, stats] = anovan(valnow, grp, 'display', showok); 
            ncom = ncom + 1; textmsg{1}{ncom} = strcat(valname{1}, '_', valname{2}, '_...');
            textmsg{2}{ncom} = num2str(p(1), '%10.5e'); textmsg{3}{ncom} = num2str(table{2,6}); textmsg{4}{ncom} = [num2str(table{2,3}), ',', num2str(table{4,3})];
            ncom = ncom + 1; textmsg{1}{ncom} = strcat('p1_', 'p2_...');
            textmsg{2}{ncom} = num2str(p(2), '%10.5e'); textmsg{3}{ncom} = num2str(table{3,6}); textmsg{4}{ncom} = [num2str(table{3,3}), ',', num2str(table{4,3})];
            
            
            [p, table, stats] = anovan(valnow, grp, 'model', 2, 'display', showok);
            ncom = ncom + 1; textmsg{1}{ncom} = strcat(valname{1}, '*', valname{2}, '*...');
            textmsg{2}{ncom} = num2str(p(1), '%10.5e'); textmsg{3}{ncom} = num2str(table{2,6}); textmsg{4}{ncom} = [num2str(table{2,3}), ',', num2str(table{5,3})];
            ncom = ncom + 1; textmsg{1}{ncom} = strcat('p1', '*', 'p2', '*...');
            textmsg{2}{ncom} = num2str(p(2), '%10.5e'); textmsg{3}{ncom} = num2str(table{3,6}); textmsg{4}{ncom} = [num2str(table{3,3}), ',', num2str(table{5,3})];
            ncom = ncom + 1; textmsg{1}{ncom} = strcat(valname{1}, '*', 'p1', '*...');
            textmsg{2}{ncom} = num2str(p(3), '%10.5e'); textmsg{3}{ncom} = num2str(table{4,6}); textmsg{4}{ncom} = [num2str(table{4,3}), ',', num2str(table{5,3})];
        else
            disp('---------------> action canceled or no values are sampled');
        end
    else    
        msgbox('---------> more than 1 set and more than 1 group of parameters are needed; ANOVA2 not performed.');
    end
elseif (strcmp(tagmark, 'kwtest'))  %Kruskal - Wallis test
    n = numel(val);
    if (n>1)
        valnow = []; grp = []; nel = []; for (i = 1:n) nel(i) = numel(val{i}); end
        valnow = zeros(1,sum(nel)); grp = cell(1,sum(nel));
        for (i = 1:n)
        valnow((sum(nel(1:i))-nel(i)+1) : (sum(nel(1:i)))) = val{i}; 
        for (j = (sum(nel(1:i))-nel(i)+1) : (sum(nel(1:i))) ) 
            grp{j} = groupname{i};
        end
        end
        [p, table, stats] = kruskalwallis(valnow, grp, showok); 
        ncom = ncom + 1; textmsg{1}{ncom} = strcat(groupname{1}, '_', groupname{2}, '_...');
        textmsg{2}{ncom} = num2str(p, '%10.5e'); textmsg{3}{ncom} = num2str(table{2,5}); textmsg{4}{ncom} = [num2str(table{2,3}), ',', num2str(table{3,3})];
    else
        msgbox('---------> more than 1 set of parameters are needed; KW test not performed.');
    end
elseif (strcmp(tagmark, 'kstest')) %Kolmogonov-Smirnov test
    n = numel(val);
    for (i = 1:n)
       for (j = i+1:n)
           [h,p,stats] = kstest2(val{i}, val{j});
           ncom = ncom + 1; textmsg{1}{ncom} = strcat(groupname{i},'_',groupname{j});
           textmsg{2}{ncom} = num2str(p, '%10.5e'); textmsg{3}{ncom} = num2str(stats); textmsg{4}{ncom} = 'Nothing to show';
       end
    end
elseif (strcmp(tagmark, 'wwtest')) %Watson-Williams test (similar to one-way ANOVA):
                                   %       whether n samples of circular variables have the same mean angles
    n = numel(val);
    if (n>1)
        valnow = []; grp = []; nel = []; for (i = 1:n) nel(i) = numel(val{i}); end
        valnow = zeros(1,sum(nel)); grp = ones(1,sum(nel)); %grp = cell(1,sum(nel));
        for (i = 1:n)
            valnow((sum(nel(1:i))-nel(i)+1) : (sum(nel(1:i)))) = mod(val{i}, 360)*2*pi/360; 
            for (j = (sum(nel(1:i))-nel(i)+1) : (sum(nel(1:i))) ) 
                grp(j) = i; %grp{j} = groupname{i};
            end
        end
        [p, table] = circ_wwtest(valnow, grp); 
        ncom = ncom + 1; textmsg{1}{ncom} = strcat(groupname{1}, '_', groupname{2}, '_...');
        textmsg{2}{ncom} = num2str(p, '%10.5e'); textmsg{3}{ncom} = num2str(table{2,5}); textmsg{4}{ncom} = 'Nothing to show';
    else
        msgbox('---------> more than 1 set of parameters are needed; WW test not performed.');
    end
    
elseif (strcmp(tagmark, 'anovan')) %%%not vailable at this time, but could be similar to anova2 (reauire a multi-dim array val{i,j,k...})    
elseif (strcmp(tagmark, 'friedman'))%%%Friedman test: test column effect after row effects are removed.not available

end

%%%%display the results
if (ncom > 0) 
    if ( strcmp(tagmark, 'anova2') | strcmp(tagmark, 'friedman') )
        fname = tagmark;
    else
        fname = strcat(plotname, '-', tagmark);
    end
    hf = figure('Name', fname, 'Units', 'inches', 'Position', [2 2 9 5]);
    TextDisplayer_multiple(hf, [0.1 0 8.8 4.8], textmsg, tit); 
    %TextDisplayer_multiple(hf, [0.05 0 0.9 0.9], textmsg, tit, 'normalized');
end

