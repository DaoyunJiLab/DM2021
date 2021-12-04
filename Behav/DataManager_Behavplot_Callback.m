function DataManager_Behavplot_Callback
%%Plot the selected results
hf = gcbf;
pinfo = getappdata(hf, 'behav'); data = getappdata(hf, 'bhdata');
hspike = getappdata(hf, 'hspike'); hgroup = getappdata(hf, 'hgroup'); groupselection = getappdata(hgroup, 'selection');
hfield = getappdata(hf, 'hfield');
valuelistpos = getappdata(hf, 'valuelistpos');
plotparm = getappdata(hf, 'plotparm');

tagmark = get(gcbo, 'Tag');
if (strcmp(tagmark, 'setbin'))
    plotparm.bin = str2num(get(gcbo, 'String'));
    setappdata(hf, 'plotparm', plotparm);
elseif (strcmp(tagmark, 'setrange'))
    plotparm.range = str2num(get(gcbo, 'String'));
    setappdata(hf, 'plotparm', plotparm); 
elseif (strcmp(tagmark, 'showdetail'))
    if (plotparm.showdetail == 1)
        plotparm.showdetail = 0; set(gcbo, 'ForegroundColor', [0.2 0.2 0.2], 'String', 'NoShow');
    else
        plotparm.showdetail = 1; set(gcbo, 'ForegroundColor', [1 0 0], 'String', 'ShowDetail');
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
        plotparm.overwrite = 0; set(gcbo, 'ForegroundColor', [0.2 0.2 0.2], 'String', 'NewDatabs');
    end
    setappdata(hf, 'plotparm', plotparm);
    
else  %%%for all other options to plot data
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
if (isempty(groupname) | isempty(fieldname))
    disp('----------> no groups or variables selected');
else
    str = {'Between groups'; 'Between variables'; 'Individuals'};
    [ss,ok] = listdlg('PromptString', 'Select a comparison type', 'SelectionMode', 'single', 'ListString', str); 
    if (ok)
        comtype = str{ss};
        plotdata(val, groupname, fieldname, plotparm, comtype, tagmark);
    end
end
end
disp('************************');

function plotdata(val, groupname, fieldname, plotparm, comtype, tagmark)
%%val{group, field}{values}
[ngroup, nfield] = size(val);
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
        else
            nel = nel + 1; valnow(nel) = str2num(val{i}{j}); 
        end
    end
    val{i} = valnow; 
    if (nel > 0)
        mm(i) = mean(valnow); nn(i) = nel; ss(i) = std(valnow)/sqrt(nel); m50(i) = median(valnow);
        m10(i) = prctile(valnow, 10); m25(i) = prctile(valnow, 25); m75(i) = prctile(valnow, 75); m90(i) = prctile(valnow, 90); 
    end
end
if (strcmp(tagmark, 'hist'))
   binvector = [plotparm.range(1):plotparm.bin:plotparm.range(2)]; 
   hf = figure('Name', strcat(plotname, '-', tagmark)); hax = axes('Parent', hf, 'NextPlot', 'add'); str = [];
   for (i = 1:n)
       co{i} = rand(1,3);
       Y = histc(val{i}, binvector); if (~strcmp(comtype, 'Individuals')) Y = Y/sum(Y); end
       bar(hax, binvector, Y, 'EdgeColor', co{i}, 'FaceColor', co{i}); 
       str{i} = strcat(groupname{i}, ': mean=', num2str(mm(i)), ';std=', num2str(ss(i)), ';n=', num2str(nn(i)), ';p10=', num2str(m10(i)),...
           ';p25=', num2str(m25(i)), ';p50=', num2str(m50(i)), ';p75=', num2str(m75(i)), ';p90=', num2str(m90(i)));
   end
   for (i = 1:numel(str))
       text('Parent', hax, 'String', str{i}, 'Color', co{i}, 'Units', 'normalized', 'Position', [0.05 0.96-(i-1)*0.04]);
   end 
elseif (strcmp(tagmark, 'cumu'))
   binvector = [plotparm.range(1):plotparm.bin:plotparm.range(2)];  
   hf = figure('Name', strcat(plotname, '-', tagmark)); hax = axes('Parent', hf, 'NextPlot', 'add'); str = [];
   for (i = 1:n)
       co{i} = rand(1,3);
       Y = histc(val{i}, binvector); if (~strcmp(comtype, 'Individuals')) Y = Y/sum(Y); end
       line(binvector, cumsum(Y), 'Parent', hax, 'Color', co{i}); 
       str{i} = strcat(groupname{i}, ': mean=', num2str(mm(i)), ';std=', num2str(ss(i)), ';n=', num2str(nn(i)), ';p10=', num2str(m10(i)),...
           ';p25=', num2str(m25(i)), ';p50=', num2str(m50(i)), ';p75=', num2str(m75(i)), ';p90=', num2str(m90(i)));
   end
   for (i = 1:numel(str))
       text('Parent', hax, 'String', str{i}, 'Color', co{i}, 'Units', 'normalized', 'Position', [0.05 0.96-(i-1)*0.04]);
   end 
elseif (strcmp(tagmark, 'polar'))
    if (strcmp(comtype, 'Individuals'))
   polarbin = [plotparm.range(1):plotparm.bin:plotparm.range(2)]; nbin = numel(polarbin);
   hf = figure('Name', strcat(plotname, '-', tagmark)); hax = axes('Parent', hf);
   %hax = axes('Parent', hf, 'NextPlot', 'add'); 
   str = []; cnow = [0 0 0]; %set(hax, 'NextPlot', 'add');
   for (i = 1:n)
        valok = [];   %co{i} = rand(1,3);
       valok = mod(val{i}, 360); %%%try to solve the multiple plots
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
       text('Parent', hax, 'String', str{i}, 'Color', cnow, 'Units', 'normalized', 'Position', [-0.4 0.96-(i-1)*0.04]);
   end 
    else
        disp('-----------> not available for this comparision type. choose the Individual option');
    end
elseif (strcmp(tagmark, 'scatter'))
   for (i = 1:n)
       str = [];
       str{1} = strcat(groupname{i}, ': mean=', num2str(mm(i)), ';std=', num2str(ss(i)), ';n=', num2str(nn(i)), ';p10=', num2str(m10(i)),...
           ';p25=', num2str(m25(i)), ';p50=', num2str(m50(i)), ';p75=', num2str(m75(i)), ';p90=', num2str(m90(i)));
       for (j = i+1:n)
           hf = figure('Name', strcat(plotname, '-', tagmark)); hax = axes('Parent', hf, 'NextPlot', 'add'); 
           nok = min([numel(val{i}) numel(val{j})]);
           if (numel(val{i}) ~= numel(val{j}))
               disp('---------> warning: element numbers do not match!');
           end
           line(val{i}(1:nok), val{j}(1:nok), 'Parent', hax, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 20);
           xlabel(groupname{i}); ylabel(groupname{j});
           str{2} = strcat(groupname{j}, ': mean=', num2str(mm(j)), ';std=', num2str(ss(j)), ';n=', num2str(nn(j)), ';p10=',...
               num2str(m10(j)), ';p25=', num2str(m25(j)), ';p50=', num2str(m50(j)), ';p75=', num2str(m75(j)), ';p90=', num2str(m90(j)));
           %%%%still need to compute corr and regression
           if (std(val{i}(1:nok)) ~= 0)
              [str1, str2, YY] = DoCrrRegress(val{i}(1:nok), val{j}(1:nok)); str{3} = str1; str{4} = str2;
              line(val{i}(1:nok), YY, 'Parent', hax, 'Color', [1 0 0]);
           end
           for (kk = 1:numel(str))
               text('Parent', hax, 'String', str{kk}, 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.05 0.96-(kk-1)*0.04]);
           end 
       end
   end
end

function [str1, str2, YY] = DoCrrRegress(xx, yy) %%%xx, yy column vectors
[RR, PP] = corrcoef(xx, yy); r = RR(1,2); p = PP(1,2); %%correlation and p value
str1 = strcat('R=', num2str(RR(1,2)), '; p= ', num2str(PP(1,2)));
aaa = size(xx); if (aaa(1)==1) xx = xx'; end
aaa = size(yy); if (aaa(1)==1) yy = yy'; end
n = numel(xx); Xparm = [ones(n,1) xx];
[BBB, BBBint, R, Rint, stat] =  regress(yy, Xparm, 0.05);
YY = Xparm * BBB; %regression result: line in (X, Y) = (xx, L4)
str2 = strcat('regression: b=',num2str(BBB(1)), '(', num2str(BBBint(1,1)), ',', num2str(BBBint(1,2)), ')', ...
    '; k=',num2str(BBB(2)), '(', num2str(BBBint(2,1)), ',', num2str(BBBint(2,2)), ')' );
    