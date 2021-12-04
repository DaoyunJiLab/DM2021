function DataManager_plot_Callback
%%Plot the selected results
hf = gcbf; dbtype = getappdata(hf, 'dbtype'); celllistvar = getappdata(hf, 'celllistvar');
%%%%%%%added options for all types of databases
pinfostr = 'pinfo'; datastr = 'data'; catok = 'general'; varok = 'parmfile'; varok1 = 'tmpname';
switch dbtype
    case '.eegdb'
        pinfostr = 'eeg'; datastr = 'eegdata'; catok = 'general'; varok = 'eegfile'; %%also determine a variable (one in pinfo.general) to identify cells
    case '.behavdb'
        pinfostr = 'behav'; datastr = 'bhdata'; catok = 'general'; varok = 'sessID'; 
    case '.seq'
        catok = 'seq'; varok = 'sequence';
    case '.seqdb'
        catok = 'general'; varok = 'tmpID'; varok1 = 'tmpname';   
end
pinfo = getappdata(hf, pinfostr); data = getappdata(hf, datastr);
hspike = getappdata(hf, 'hspike'); hgroup = getappdata(hf, 'hgroup'); groupselection = getappdata(hgroup, 'selection');
hfield = getappdata(hf, 'hfield'); valuelistpos = getappdata(hf, 'valuelistpos'); plotparm = getappdata(hf, 'plotparm');

tagmark = get(gcbo, 'Tag');
if (strcmp(tagmark, 'setbin'))
    plotparm.bin = str2num(get(gcbo, 'String')); setappdata(hf, 'plotparm', plotparm);
elseif (strcmp(tagmark, 'setrange'))
    plotparm.range = str2num(get(gcbo, 'String')); setappdata(hf, 'plotparm', plotparm); 
elseif (strcmp(tagmark, 'showdetail'))
    if (plotparm.showdetail == 1)
        plotparm.showdetail = 0; set(gcbo, 'ForegroundColor', [0.2 0.2 0.2], 'String', 'NoShow');
    else
        plotparm.showdetail = 1; set(gcbo, 'ForegroundColor', [1 0 0], 'String', 'Show');
    end
    setappdata(hf, 'plotparm', plotparm);
elseif (strcmp(tagmark, 'setlog'))
    if (plotparm.setlog == 0)
        plotparm.setlog = 1; set(gcbo, 'ForegroundColor', [1 0 0], 'String', 'LogX');
    elseif (plotparm.setlog == 1)
        plotparm.setlog = 2; set(gcbo, 'ForegroundColor', [1 0 0], 'String', 'LogY');
    elseif (plotparm.setlog == 2)
        plotparm.setlog = 3; set(gcbo, 'ForegroundColor', [1 0 0], 'String', 'LogXY');
    elseif (plotparm.setlog == 3)
        plotparm.setlog = 0; set(gcbo, 'ForegroundColor', [0.2 0.2 0.2], 'String', 'Linear');
    end
    setappdata(hf, 'plotparm', plotparm);    
elseif (strcmp(tagmark, 'assignparm'))
    if (plotparm.compute == 0)
        plotparm.compute = 1; set(gcbo, 'ForegroundColor', [1 0 0], 'String', 'Compute');
    else
        plotparm.compute = 0; set(gcbo, 'ForegroundColor', [0.2 0.2 0.2], 'String', 'AssignParm');
    end
    setappdata(hf, 'plotparm', plotparm);
elseif (strcmp(tagmark, 'overwrite'))
    if (plotparm.overwrite == 0)
        plotparm.overwrite = 1; set(gcbo, 'ForegroundColor', [1 0 0], 'String', 'Overwrite');
    else
        plotparm.overwrite = 0; set(gcbo, 'ForegroundColor', [0.2 0.2 0.2], 'String', 'NewDB');
    end
    setappdata(hf, 'plotparm', plotparm);
elseif (~isempty(strfind(tagmark, 'setev')))
    if (plotparm.evselect == 0)
        plotparm.evselect = 1; set(gcbo, 'ForegroundColor', [1 0 0], 'String', 'Events');
    else
        plotparm.evselect = 0; set(gcbo, 'ForegroundColor', [0.2 0.2 0.2], 'String', 'Session');
    end
    setappdata(hf, 'plotparm', plotparm);
elseif (~isempty(strfind(tagmark, 'link')))
    DataManager_LinkDatabase_Callback(hf, plotparm, tagmark);    
else  %%%for all other options to plot data
%determine selected groups and selected variables
groupselectionindex = find( groupselection == 1); %index for selected groups
ngroup = numel(groupselectionindex);
fieldtitle = fieldnames(pinfo); nfield = numel(fieldtitle); 
val = []; groupname = []; fieldname = []; cellname = []; %all values for plot are assigned here
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
    if isfield(pinfo.(catok), varok)
        cellname{k} = pinfo.(catok).(varok)(spikeselectionindex);
    else
        cellname{k} = pinfo.(catok).(varok1)(spikeselectionindex);
    end
end
if (isempty(groupname) || isempty(fieldname))
    disp('----------> no groups or variables selected');
else
    str = {'Between groups'; 'Between variables'; 'Individuals'};
    [ss,ok] = listdlg('PromptString', 'Select a comparison type', 'SelectionMode', 'single', 'ListString', str); 
    if (ok)
        comtype = str{ss};
        plotdata(val, groupname, fieldname, plotparm, comtype, tagmark, cellname);
    end
end
end
%disp('************************');

function plotdata(val, groupname, fieldname, plotparm, comtype, tagmark, cellname)
col{1} = [0 0 1]; col{2} = [1 0 0]; col{3} = [0 0 0]; col{4} = [0.5 0.5 0.5]; col{5} = [1 1 0]; col{6} = [1 0 1]; 
%%val{group, field}{values}
[ngroup, nfield] = size(val);
if (strcmp(tagmark, 'function')) %%%if do a function between two variables
    mmmark = 'mean'; ok = 1;
    if strncmpi(comtype, 'Between', 5)
        str = {'Mean'; 'Median'};
        [ss,ok] = listdlg('PromptString', 'Select an average method', 'SelectionMode', 'single', 'ListString', str); 
        if ok mmmark = str{ss}; end
    end
    if ok plotdatafunction(val, groupname, fieldname, plotparm, comtype, tagmark, cellname, col, mmmark); end
else  %%for all non-functional plots
if (strcmp(comtype, 'Between groups'))
    for (j = 1:nfield)
        valnow = [];
        for (i = 1:ngroup)
            valnow{i} = val{i,j};
        end
        plotnow(valnow, fieldname{j}, groupname, plotparm, tagmark, comtype);
    end
elseif (strcmp(comtype, 'Between variables'))
    for (i = 1:ngroup)
        valnow = [];
        for (j = 1:nfield)
            valnow{j} = val{i,j};
        end
        plotnow(valnow, groupname{i}, fieldname, plotparm, tagmark, comtype);
    end
elseif (strcmp(comtype, 'Individuals'))
    for (i = 1:ngroup)
        for (j = 1:nfield)
            valnow{1} = val{i,j}; groupnow{1} = groupname{i};
            plotnow(valnow, fieldname{j}, groupnow, plotparm, tagmark, comtype);
        end
    end
end
end

function plotnow(val, plotname, groupname, plotparm, tagmark, comtype)
%val{n}
%%%%arrange data and compute data statistics first
n = numel(val); mm = NaN*ones(1,n); ss = NaN*ones(1,n); nn = NaN*ones(1,n); 
m10 = NaN*ones(1,n); m25 = NaN*ones(1,n); m50 = NaN*ones(1,n); m75 = NaN*ones(1,n); m90 = NaN*ones(1,n);
for (i = 1:n) %group index
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
        m10(i) = prctile(valnow, 10); m25(i) = prctile(valnow, 25); m75(i) = prctile(valnow, 75); m90(i) = prctile(valnow, 90); 
    end
end
if (strcmp(tagmark, 'hist'))
   binvector = [plotparm.range(1):plotparm.bin:plotparm.range(2)]; ok = 1; selection = 1;
   if (~strcmp(comtype, 'Individuals'))
      sel = {'Bar histogram'; 'Bar pdf'; 'Line histogram'; 'Line pdf'};
      [selection,ok] = listdlg('PromptString', 'Plot option', 'SelectionMode', 'single', 'ListString', sel); 
   end
   if ok
      hf = figure('Name', strcat(plotname, '-', tagmark)); hax = axes('Parent', hf, 'NextPlot', 'add'); str = [];
      for (i = 1:n)
          co{i} = rand(1,3);
          Y = histc(val{i}, binvector); if (selection==2)||(selection==4) Y = Y/sum(Y); end
          if (plotparm.setlog == 1)
              binvector = log10(binvector);
          elseif (plotparm.setlog == 2)
              Y = log10(Y);
          elseif (plotparm.setlog == 3)
              binvector = log10(binvector); Y = log10(Y);
          end
          if (selection==1)||(selection==2)%%% plot as bar histogram
              bar(hax, binvector, Y, 'EdgeColor', co{i}, 'FaceColor', co{i}); 
          else                 
              line(binvector, Y, 'Parent', hax, 'Color', co{i});
          end
          str{i} = strcat(groupname{i}, ': mean=', num2str(mm(i)), ';se=', num2str(ss(i)), ';n=', num2str(nn(i)), ';p10=', num2str(m10(i)),...
              ';p25=', num2str(m25(i)), ';p50=', num2str(m50(i)), ';p75=', num2str(m75(i)), ';p90=', num2str(m90(i)));
      end
      for (i = 1:numel(str))
          text('Parent', hax, 'String', str{i}, 'Interpreter', 'none', 'Color', co{i}, 'Units', 'normalized', 'Position', [0.05 0.96-(i-1)*0.04]);
      end 
   end
elseif (strcmp(tagmark, 'bar')) %%%bar plot mean +- se or median +- wiskers values
   binvector = 1:n; cl = rand(1,3);
   sel = {'Mean'; 'Median'};
   [selection,ok] = listdlg('PromptString', 'Plot option', 'SelectionMode', 'single', 'ListString', sel); 
   if (ok)
        if (selection == 1)
            Y = mm; UP = Y +ss; DOWN = Y- ss;
        elseif (selection == 2)
            Y = m50; UP = m25; DOWN = m75; top = m10; bottom = m90; 
        end
   end
   hf = figure('Name', strcat(plotname, '-', tagmark)); hax = axes('Parent', hf, 'NextPlot', 'add'); str = [];
   if (plotparm.setlog == 1)
       binvector = log10(binvector);
   elseif (plotparm.setlog == 2)
       Y = log10(Y); UP = log10(UP); DOWN = log10(DOWN);
   elseif (plotparm.setlog == 3)
       binvector = log10(binvector); Y = log10(Y); UP = log10(UP); DOWN = log10(DOWN);
   end
   if selection == 1
      bar(hax, binvector, Y, 'EdgeColor', cl, 'FaceColor', cl); 
      Drawerrorupdown(binvector, Y, UP, DOWN, hax, [1 0 0]);
   else
      Drawerrorupdown_box(binvector, Y, UP, DOWN, hax, cl, top, bottom);
   end
   str = cell(1, n);
   for (i = 1:n)
       str{i} = strcat(groupname{i}, ': mean=', num2str(mm(i)), ';se=', num2str(ss(i)), ';n=', num2str(nn(i)), ';p10=', num2str(m10(i)),...
           ';p25=', num2str(m25(i)), ';p50=', num2str(m50(i)), ';p75=', num2str(m75(i)), ';p90=', num2str(m90(i)));
   end
   for (i = 1:numel(str))
       text('Parent', hax, 'String', str{i}, 'Interpreter', 'none', 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.05 0.96-(i-1)*0.04]);
   end 
elseif (strcmp(tagmark, 'scatter')) %%%bar plot either mean or median
   binvector = 1:n; cl = rand(1,3); Y = []; 
   sel = {'Mean'; 'Median'; 'Pair + mean'; 'Pair + median'};
   [selection,ok] = listdlg('PromptString', 'Plot option', 'SelectionMode', 'single', 'ListString', sel); 
   if (ok)
        if (selection == 1)||(selection == 3)
            Y = mm;
        elseif (selection == 2)||(selection == 4)
            Y = m50;
        end
   end
   if ok
      hf = figure('Name', strcat(plotname, '-', tagmark)); hax = axes('Parent', hf, 'NextPlot', 'add'); str = [];
      if (plotparm.setlog == 1)
          binvector = log10(binvector);
      elseif (plotparm.setlog == 2)
          Y = log10(Y); for (i = 1:n) val{i} = log10(val{i}); end
      elseif (plotparm.setlog == 3)
          binvector = log10(binvector); Y = log10(Y); for (i = 1:n) val{i} = log10(val{i}); end
      end
      %%%%plot data
      if (selection > 2.5)&&(n==2)&& (numel(val{1})==numel(val{2}))  %%%%if pairwise plot
         for i = 1:numel(val{1})
             line(binvector, [val{1}(i) val{2}(i)], 'Parent', hax, 'LineWidth', 0.5, 'Color', [0.5 0.5 0.5],...
                 'Marker', '.', 'MarkerSize', 20);
         end
      else
         for (i = 1:n)
              xv = computexvalues(val{i}, Y(i), 0.35); %%%(yvalues, middlevalue, range) 
              line(binvector(i)+xv, val{i}, 'Parent', hax, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 20);
              %line(binvector(i)+0.3*(rand(1, numel(val{i}))-0.5), val{i}, 'Parent', hax, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 20);
         end
      end
      %%%%plot mean/median
      %bar(hax, binvector, Y, 'EdgeColor', cl, 'FaceColor', cl);
      for (i = 1:n)
          line(binvector(i)+[-0.2 0.2], [Y(i) Y(i)], 'Parent', hax, 'Color', [1 0 0], 'LineWidth', 2); %
      end
%       for (i = 1:n)
%           line(binvector(i)*ones(1, numel(val{i})), val{i}, 'Parent', hax, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 20);
%       end
      str = cell(1, n);
      for (i = 1:n)
           str{i} = strcat(groupname{i}, ': mean=', num2str(mm(i)), ';se=', num2str(ss(i)), ';n=', num2str(nn(i)), ';p10=', num2str(m10(i)),...
           ';p25=', num2str(m25(i)), ';p50=', num2str(m50(i)), ';p75=', num2str(m75(i)), ';p90=', num2str(m90(i)));
      end
      for (i = 1:numel(str))
          text('Parent', hax, 'String', str{i},'Interpreter', 'none', 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.05 0.96-(i-1)*0.04]);
      end
   end
elseif (strcmp(tagmark, 'cumu'))
   binvector = [plotparm.range(1):plotparm.bin:plotparm.range(2)];  
   hf = figure('Name', strcat(plotname, '-', tagmark)); hax = axes('Parent', hf, 'NextPlot', 'add'); str = [];
   for (i = 1:n)
       co{i} = rand(1,3);
       Y = histc(val{i}, binvector); if (~strcmp(comtype, 'Individuals')) Y = Y/sum(Y); end
       line(binvector, cumsum(Y), 'Parent', hax, 'Color', co{i}); 
       str{i} = strcat(groupname{i}, ': mean=', num2str(mm(i)), ';se=', num2str(ss(i)), ';n=', num2str(nn(i)), ';p10=', num2str(m10(i)),...
           ';p25=', num2str(m25(i)), ';p50=', num2str(m50(i)), ';p75=', num2str(m75(i)), ';p90=', num2str(m90(i)));
   end
   for (i = 1:numel(str))
       text('Parent', hax, 'String', str{i}, 'Interpreter', 'none', 'Color', co{i}, 'Units', 'normalized', 'Position', [0.05 0.96-(i-1)*0.04]);
   end 
elseif (strcmp(tagmark, 'polar'))
    if (strcmp(comtype, 'Individuals'))
   polarbin = [plotparm.range(1):plotparm.bin:plotparm.range(2)]; nbin = numel(polarbin);
   hf = figure('Name', strcat(plotname, '-', tagmark)); hax = axes('Parent', hf);
   %hax = axes('Parent', hf, 'NextPlot', 'add'); 
   str = []; cnow = [0 0 0]; %set(hax, 'NextPlot', 'add');
   for (i = 1:n)
       valok = val{i}; valok = valok(~isnan(valok));   %co{i} = rand(1,3);
       valok = mod(valok, 360); %%%try to solve the multiple plots
       %polarrate = histc(valok, polarbin); if (~strcmp(comtype, 'Individuals')) polarrate = polarrate/sum(polarrate); end
       %h = polar(hax, polarbin, polarrate, 'r'); set(h, 'Color', co{i}); 
       h = rose(hax, valok/360*2*pi, nbin); set(h, 'Color', cnow, 'LineWidth', 2);
       valang = valok'/360*2*pi; angfac = 360/2/pi;
       mk = mod(angfac*circ_mean(valang),360); r = circ_r(valang); s = mod(angfac*circ_std(valang), 360); 
       m550 = mod(angfac*circ_median(valang), 360); p = circ_rtest(valang); 
       yy = r*sin(2*pi*mk/360); xx = r*cos(mk*2*pi/360); line([0 xx], [0 yy], 'Parent', hax, 'Color', [1 0 0], 'LineWidth',4);
       str{i} = strcat(groupname{i}, ': circular mean=', num2str(mk), '; circular std=', num2str(s), '; n=', num2str(nn(i)),...
           '; circular median=', num2str(m550), '; Raleigh r=', num2str(r), '; Rayleigh p=', num2str(p));
   end
   for (i = 1:numel(str))
       text('Parent', hax, 'String', str{i},'Interpreter', 'none', 'Color', cnow, 'Units', 'normalized', 'Position', [-0.4 0.96-(i-1)*0.04]);
   end 
    else
        disp('-----------> not available for this comparision type. choose the Individual option');
    end
elseif (strcmp(tagmark, 'regression'))
   for (i = 1:n)
       str = [];
       str{1} = strcat(groupname{i}, ': mean=', num2str(mm(i)), ';se=', num2str(ss(i)), ';n=', num2str(nn(i)), ';p10=', num2str(m10(i)),...
           ';p25=', num2str(m25(i)), ';p50=', num2str(m50(i)), ';p75=', num2str(m75(i)), ';p90=', num2str(m90(i)));
       for (j = i+1:n)
           hf = figure('Name', strcat(plotname, '-', tagmark)); hax = axes('Parent', hf, 'NextPlot', 'add'); 
           nok = min([numel(val{i}) numel(val{j})]); 
           if (numel(val{i}) ~= numel(val{j}))
               disp('---------> warning: element numbers do not match!');
           end
           aa = val{i}(1:nok); bb = val{j}(1:nok);
           iii = find((~isnan(aa))&(~isnan(bb))); aa=aa(iii); bb=bb(iii); 
           line(aa, bb, 'Parent', hax, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 12, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0]);
           xlabel(groupname{i}); ylabel(groupname{j});
           str{2} = strcat(groupname{j}, ': mean=', num2str(mm(j)), ';se=', num2str(ss(j)), ';n=', num2str(nn(j)), ';p10=',...
               num2str(m10(j)), ';p25=', num2str(m25(j)), ';p50=', num2str(m50(j)), ';p75=', num2str(m75(j)), ';p90=', num2str(m90(j)));
           %%%%still need to compute corr and regression
           if (numel(aa)>=2) && (std(aa)~=0)
              [str1, str2, YY, XX] = DoCrrRegress(aa, bb); str{3} = [str1 '; N= ' num2str(numel(aa))]; str{4} = str2;
              %line(val{i}(1:nok), YY, 'Parent', hax, 'Color', [1 0 0]);
              line(XX, YY, 'Parent', hax, 'Color', [1 0 0], 'LineWidth', 1);
           end
           for (kk = 1:numel(str))
               text('Parent', hax, 'String', str{kk},'Interpreter', 'none', 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.05 0.96-(kk-1)*0.04]);
           end 
       end
   end
end

function xv = computexvalues(Y, mm, range) %%%(yvalues, range)
nbin = 3;
xv = range*(rand(1, numel(Y))-0.5); %default randamization
Y = Y - mm; binsize = max(abs(Y))/nbin; %%%yvalues divided into +- n bins
ybin = [-nbin:nbin];
for (i = 1:numel(ybin))
    iii = find( (Y>= (ybin(i)-0.5)*binsize) & (Y<(ybin(i)+0.5)*binsize) );
    if (ybin(i) == 0)
        rangenow = range;
    else
        rangenow = range/abs(ybin(i));
    end
    xv(iii) = rangenow * (rand(1, numel(iii))-0.5);
end

function plotdatafunction(val, groupname, fieldname, plotparm, comtype, tagmark, cellname, col, mmmark)
[ngroup, nfield] = size(val); ncell = numel(val{1,1});
    if (strcmp(comtype, 'Between groups'))
        for (m = 1:nfield) %%%%the earlier variable plotted as x
             for (n = m+1:nfield) %%%%the earlier variable plotted as x
                 XX = cell(1, ngroup); YY = cell(1, ngroup); EE = cell(1, ngroup); FF = cell(1, ngroup); str = cell(1, ngroup);
                 for (i = 1:ngroup)
                     allXX = val{i,m}; allYY = val{i,n}; longX = []; longY = [];
                     str{i} = strcat(groupname{i}, ': n=', num2str(numel(allXX)));
                     for (j = 1:numel(allXX)) %%%cel lindex
                         if (numel(allXX{j}) ~= numel(allYY{j}))
                            disp('---------> warning: element numbers do not match!');
                            nok = min([numel(allXX{j}) numel(allYY{j})]); allXX{j} = allXX{j}(1:nok); allYY{j} = allYY{j}(1:nok);
                         end
                         [tt1, tt2] = size(allXX{j}); if (tt1<tt2) allXX{j} = allXX{j}'; end
                         [tt1, tt2] = size(allYY{j}); if (tt1<tt2) allYY{j} = allYY{j}'; end
                         XX{i} = union(XX{i}, allXX{j}); longX = [longX;allXX{j}]; longY = [longY;allYY{j}]; 
                     end
                     XX{i} = sort(XX{i}); YY{i} = NaN*ones(size(XX{i})); EE{i} = NaN*ones(size(XX{i})); FF{i} = NaN*ones(size(XX{i})); 
                     for (j = 1:numel(XX{i}))
                         iii = find( longX == XX{i}(j) ); ynow = longY(iii); ynow = ynow(~isnan(ynow)); 
                         if (~isempty(ynow))
                             if strncmpi(mmmark, 'mean', 3)
                                YY{i}(j) = mean(ynow); EE{i}(j) = std(ynow)/sqrt(numel(ynow));
                             else
                                YY{i}(j) = median(ynow); EE{i}(j) = prctile(ynow, 75) - YY{i}(j); FF{i}(j) = YY{i}(j) - prctile(ynow, 25);
                             end
                         end
                     end
                 end
                 for (i = 1:ngroup)
                     if strncmpi(mmmark, 'mean', 3)
                        UP{i} = YY{i} + EE{i}; DD{i} = YY{i}-EE{i};
                     else
                        UP{i} = YY{i} + EE{i}; DD{i} = YY{i}-FF{i}; 
                     end
                     if (plotparm.setlog == 0)
                        strx = fieldname{m}; stry = fieldname{n};
                     elseif (plotparm.setlog == 1) 
                        XX{i} = log10(XX{i}); strx = strcat('log10(', fieldname{m}, ')'); stry = fieldname{n};
                     elseif (plotparm.setlog == 2) 
                        YY{i} = log10(YY{i}); UP{i} = log10(UP{i}); DD{i} = log10(DD{i}); stry = strcat('log10(', fieldname{n}, ')'); strx = fieldname{m};
                     elseif (plotparm.setlog == 3) 
                        XX{i} = log10(XX{i}); YY{i} = log10(YY{i}); UP{i} = log10(UP{i}); DD{i} = log10(DD{i}); stry = strcat('log10(', fieldname{n}, ')'); strx = strcat('log10(',fieldname{m}, ')');
                     end
                 end
                 hf = figure('Name', strcat(groupname{i}, '-', tagmark)); hax = axes('Parent', hf, 'NextPlot', 'add');
                 xlabel(strx); ylabel(stry)
                 for (i = 1:ngroup)
                     line(XX{i}, YY{i},  'Parent', hax, 'LineWidth', 2, 'Color', col{i});
                     Drawerrorupdown(XX{i}, YY{i}, UP{i}, DD{i}, hax, col{i});
                 end
                 for (i = 1:ngroup)
                     text('Parent', hax, 'Interpreter', 'none', 'String', str{i}, 'Color', col{i}, 'Units', 'normalized', 'Position', [0.05 0.96-(i-1)*0.04]);
                 end
             end
        end
    elseif (strcmp(comtype, 'Individuals'))
        for (i = 1:ngroup)
            for (j = 1:numel(cellname{i}))
                for (m = 1:nfield) %%%%the earlier variable plotted as x
                    for (n = m+1:nfield) %%%%the earlier variable plotted as x
                        hf = figure('Name', strcat(groupname{i}, '-', tagmark)); hax = axes('Parent', hf, 'NextPlot', 'add');
                        nok = min([numel(val{i,m}{j}) numel(val{i,n}{j})]);
                        if (numel(val{i,m}{j}) ~= numel(val{i,n}{j}))
                           disp('---------> warning: element numbers do not match!');
                        end
                        if (plotparm.setlog == 0)
                            YY = val{i,n}{j}(1:nok); XX = val{i,m}{j}(1:nok); strx = fieldname{m}; stry = fieldname{n};
                        elseif (plotparm.setlog == 1) 
                            YY = val{i,n}{j}(1:nok); XX = log10(val{i,m}{j}(1:nok)); 
                            strx = strcat('log10(', fieldname{m}, ')'); stry = fieldname{n};
                        elseif (plotparm.setlog == 2) 
                            YY = log10(val{i,n}{j}(1:nok)); XX = val{i,m}{j}(1:nok); 
                            stry = strcat('log10(', fieldname{n}, ')'); strx = fieldname{m};
                        elseif (plotparm.setlog == 3) 
                            YY = log10(val{i,n}{j}(1:nok)); XX = log10(val{i,m}{j}(1:nok)); 
                            stry = strcat('log10(', fieldname{n}, ')'); strx = strcat('log10(',fieldname{m}, ')');
                        end
                        line(XX, YY,  'Parent', hax, 'LineWidth', 2); %%%%the earlier variable plotted as x
                        xlabel(strx); ylabel(stry); %%%%the earlier variable plotted as x
                        text('Parent', hax, 'String', cellname{i}{j}, 'Interpreter', 'none', 'Units', 'normalized', 'Position', [0.05 0.96]);
                    end
                end
            end
        end
    elseif (strcmp(comtype, 'Between variables'))
        msgbox('---------> Comparison type Between variables not available for this type of plot.');
    end

function [str1, str2, YY, xx] = DoCrrRegress(xx, yy) %%%xx, yy column vectors
ind = find( (~isnan(xx)) & (~isnan(yy)) ); xx= xx(ind); yy = yy(ind); %%%get rid of the non-numbers
[RR, PP] = corrcoef(xx, yy); r = RR(1,2); p = PP(1,2); %%correlation and p value
str1 = strcat('R=', num2str(RR(1,2)), '; p= ', num2str(PP(1,2)));
aaa = size(xx); if (aaa(1)==1) xx = xx'; end
aaa = size(yy); if (aaa(1)==1) yy = yy'; end
n = numel(xx); Xparm = [ones(n,1) xx];
[BBB, BBBint, R, Rint, stat] =  regress(yy, Xparm, 0.05);
YY = Xparm * BBB; %regression result: line in (X, Y) = (xx, L4)
str2 = strcat('regression: b=',num2str(BBB(1)), '(', num2str(BBBint(1,1)), ',', num2str(BBBint(1,2)), ')', ...
    '; k=',num2str(BBB(2)), '(', num2str(BBBint(2,1)), ',', num2str(BBBint(2,2)), ')' );
    