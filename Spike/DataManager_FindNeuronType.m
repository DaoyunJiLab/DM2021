function [pinfo,data] = DataManager_FindNeuronType(pinfo,data, cellind, vv)
%classify all neurons in pinfo as interneurons or principal cells
%classify all cells from same area (hippo or cortex) on all recording days at the same time
%features used to classify cell type:
%      runmeanrate = pinfo.firing.runmeanrate
%      maxhlwd = pinfo.waveform.maxhlwd;
%      maxampratio = pinfo.waveform.maxnpampratio
%%%%%%criteria:
%Pyramidal('pyr'; principal, excitatory) cell: (I) 0 <= runmeanrate <= Tmeanrate if exists
%                (II) Tmaxhlwd >= maxhlwd
%                (III) ampration <= Tampratio
%Interneuron ('int'; inhibitory):    (I) runmeanrate > Trunmeanrate if exists
%                (III) Tmaxwd < maxhlwd
%                (IV) ampration > Tampratio
%Unknown ('un'): otherwise
%[]: unclassified

%require at least one of the following parameters
%   pinfo.parm.CelltypeTmeanrate, pinfo.parm.CelltypeThlwd, pinfo.parm.CelltypeTampratio
%require at least one of the following (paired with corresponding parameters) variables
%   pinfo.firing.runmeanrate, pinfo.waveform.maxhlwd, pinfo.waveform.maxnpampratio

%variable to assign
if (~isempty(cellind))
    nspike = numel(pinfo.general.parmfile);
    if (~isfield(pinfo.general, 'celltype')) pinfo.general.celltype = cell(1, nspike); end % = 'prin', 'int' or 'un'
end

%%%%%get thresholds
values = NaN*ones(numel(cellind),3); CT = cell(numel(cellind),1); Pname = cell(1,3);
for (jjjk = 1:numel(cellind))
    i = cellind(jjjk);
    threshold = NaN*ones(1,3); ok = 0;
    disp(strcat('-----> classify cel ltype ---', pinfo.general.parmfile{i}));
    if (isfield(pinfo.parm, 'CelltypeTmeanrate')) & (~isnan(pinfo.parm.CelltypeTmeanrate(i))) 
        threshold(1) = pinfo.parm.CelltypeTmeanrate(i); values(jjjk,1) = findrunmeanrate(pinfo,i); Pname{1} = 'runmeanrate(Hz)';
        ok = 1;
    end
    if (isfield(pinfo.parm, 'CelltypeThlwd')) & (~isnan(pinfo.parm.CelltypeThlwd(i)))...
            & (~isempty(pinfo.waveform.maxhlwd{i}))
        threshold(2) = pinfo.parm.CelltypeThlwd(i); values(jjjk,2) = pinfo.waveform.maxhlwd{i}; Pname{2} = 'maxhlwd(s)';
        ok = 1;
    end
    if (isfield(pinfo.parm, 'CelltypeTampratio')) & (~isnan(pinfo.parm.CelltypeTampratio(i)))...
            & (~isempty(pinfo.waveform.maxnpampratio{i}));
        threshold(3) = pinfo.parm.CelltypeTampratio(i); values(jjjk,3) = pinfo.waveform.maxnpampratio{i}; Pname{3} = 'maxnpampratio';
        ok = 1;
    end
    if ok
        CT{jjjk} = findtypenow(threshold, values(jjjk,:)); pinfo.general.celltype{i} = CT{jjjk};
    else
       disp('---------> not enough parameters or variables to classify the neuron');
    end
end
if (ok & vv) plotclassification(values, CT, Pname); end

function value = findrunmeanrate(pinfo,i)
value = NaN; evtype = pinfo.parm.eventtype{i}; iii = find( strcmp(evtype, 'run') );
if (~isempty(iii))
    kkk = pinfo.firing.evtmeanrate{i}(iii);
   value = mean(kkk(~isnan(kkk)));
else
   sesstype = pinfo.parm.sessType{i}; iii = find( strcmp(sesstype, 'linear') );
   if (~isempty(iii))
       kkk = pinfo.firing.sessmeanrate{i}(iii); value = mean(kkk(~isnan(kkk)));
   else
       iii = find( strcmp(sesstype, 'open') );
       if (~isempty(iii))
           kkk = pinfo.firing.sessmeanrate{i}(iii); value = mean(kkk(~isnan(kkk)));
       end
   end
end


function celltype = findtypenow(threshold, values)
celltype = []; %%%% 'pyr'; 'int'; 'un'; %classify three classes: principal, interneuron, and unsure type
%%%% Need time to develop JClust - for now just classify based on empirical parameter thresholding
%%%%%%%%%%%%%after examining the data, only use mean rate to classify
ok = NaN*ones(1,numel(threshold));
if(~isnan(threshold(1))) | (~isnan(values(1))) %first runmeanrate
   if (values(1)<=threshold(1))
       ok(1) = 1;
   else
       ok(1) = -1;
   end
end
if ok(1) == 1
   celltype = 'pyr';
elseif ok(1) == -1
    celltype = 'int';
else
    celltype = 'un';
end
% if(~isnan(threshold(2))) | (~isnan(values(2))) %second maxhlwd
%    if (values(2)>=threshold(2))
%        ok(2) = 1;
%    else
%        ok(2) = -1;
%    end
% end
% if(~isnan(threshold(3))) | (~isnan(values(3))) %third maxnpampratio
%    if (values(3)<=threshold(3))
%        ok(3) = 1;
%    else
%        ok(3) = -1;
%    end
% end
% ok = ok(find(ok>-10));
% if (~isempty(ok))
%     if sum(ok)/numel(ok) == 1
%         celltype = 'pyr';
%     elseif sum(ok)/numel(ok) == -1
%         celltype = 'int';
%     else
%         celltype = 'un';
%     end
% end

function plotclassification(values, CT, Pname)
nn = numel(CT); fname = strcat('NeuronTypeClassification');
ind1 = find(strcmp(CT, 'pyr')); ind2 = find(strcmp(CT, 'int')); ind3 = find(strcmp(CT, 'un')); 
ind4 = setdiff( (1:numel(CT)), union(ind1, union(ind2, ind3)) );
str1 =['pyr - ', num2str(numel(ind1))]; str2 =['int - ', num2str(numel(ind2))]; 
str3 =['unknown - ', num2str(numel(ind3))]; str4 =['not classified - ', num2str(numel(ind4))]; 
iV = 0;
for (i = 1:numel(Pname))
     if (~isempty(Pname{i})) iV = iV + 1; end
end
if (iV <= 1)
   allstr = ['------------> ', str1, '; ', str2, '; ', str3, '; ', str4];
   msgbox(allstr, 'Classification result');
else
for (i = 1:numel(Pname))
    xx = values(:,i); xname = Pname{i}; 
    for (j = i+1:numel(Pname))
        yy = values(:,j); yname = Pname{j};
        if (~isempty(xname)) & (~isempty(yname))
            hf = figure('Name', fname); hax = axes('Parent', hf, 'NextPlot', 'add');
            xlabel(xname); ylabel(yname);
            line(xx(ind1), yy(ind1), 'Parent', hax, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 20, 'MarkerFaceColor', [1 0 0],...
                'MarkerEdgeColor', [1 0 0]);
            line(xx(ind2), yy(ind2), 'Parent', hax, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 20, 'MarkerFaceColor', [0 0 1],...
                'MarkerEdgeColor', [0 0 1]);
            line(xx(ind3), yy(ind3), 'Parent', hax, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 20, 'MarkerFaceColor', [0 0 0],...
                'MarkerEdgeColor', [0 0 0]);
            line(xx(ind4), yy(ind4), 'Parent', hax, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 20, 'MarkerFaceColor', [0.5 0.5 0.5],...
                'MarkerEdgeColor', [0.5 0.5 0.5]);
            text('String', str1, 'Parent', hax, 'Color', [1 0 0], 'Units', 'normalized', 'Position', [0.02 0.96]);
            text('String', str2, 'Parent', hax, 'Color', [0 0 1], 'Units', 'normalized', 'Position', [0.02 0.92]);
            text('String', str3, 'Parent', hax, 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.02 0.88]);
            text('String', str4, 'Parent', hax, 'Color', [0.5 0.5 0.5], 'Units', 'normalized', 'Position', [0.02 0.84]);
        end
    end
end
end
        

