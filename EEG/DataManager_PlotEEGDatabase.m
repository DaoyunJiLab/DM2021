function DataManager_PlotEEGDatabase(hmain, eeg, eegdata)
%%plot a eeg eegdatabase structure eeg in figure hmain
%%here  hspike == hdate; spikeselection = dateselection, etc.

% %set up plot dimensions
% figpos = get(hmain, 'Position'); %should be in inches
% %%value display and field display should be propotional to figpos
% spikelistwidth = 2; grouplistwidth = 1.5; controlwidth = 1.5; %all in inches
% spikelistheight = 0.3 * figpos(4); %percent of fig height, = field height
% valueheight = figpos(4) - spikelistheight;
% valuelistwidth = figpos(3) - controlwidth;
% fieldwidth = figpos(3) - controlwidth - spikelistwidth - grouplistwidth;
% spikelistpos = [0.1 valueheight+0.1 spikelistwidth-0.1 spikelistheight-0.1]; %all in inches
% grouplistpos = [spikelistwidth+0.1 valueheight+0.1 grouplistwidth-0.1 spikelistheight-0.1]; %all in inches
% valuelistpos = [0.1 0.1 valuelistwidth-0.1 valueheight-0.1]; 
% fpos = [spikelistwidth+grouplistwidth+0.1 valueheight+0.1 fieldwidth-0.1 spikelistheight-0.1]; 
% controlpos = [valuelistwidth 0 controlwidth figpos(4)];
%set up plot dimensions

figpos = get(hmain, 'Position'); %should be in inches
%%value display and field display should be propotional to figpos
spikelistwidth = 0.15; grouplistwidth = 0.1; controlwidth = 0.16; %all in normalized form
spikelistheight = 0.3; % * figpos(4); %percent of fig height
valueheight = 1 - spikelistheight; %figpos(4) - spikelistheight;
valuelistwidth = 1 - controlwidth; %figpos(3) - controlwidth;
fieldwidth = 1 - controlwidth - spikelistwidth - grouplistwidth; %figpos(3) - controlwidth - spikelistwidth - grouplistwidth;
spikelistpos = [0.001 valueheight+0.001 spikelistwidth-0.001 spikelistheight-0.001]; %[0.1 valueheight+0.1 spikelistwidth-0.1 spikelistheight-0.1]; %all in normalized unit
grouplistpos = [spikelistwidth+0.001 valueheight+0.001 grouplistwidth-0.001 spikelistheight-0.001]; %[spikelistwidth+0.1 valueheight+0.1 grouplistwidth-0.1 spikelistheight-0.1]; %all in normalized unit
valuelistpos = [0.001 0.001 valuelistwidth-0.001 valueheight-0.001]; %[0.1 0.1 valuelistwidth-0.1 valueheight-0.1]; 
fpos = [spikelistwidth+grouplistwidth valueheight+0.001 fieldwidth spikelistheight-0.001]; %[spikelistwidth+grouplistwidth+0.1 valueheight+0.1 fieldwidth-0.1 spikelistheight-0.1]; 
controlpos = [valuelistwidth 0 controlwidth 1]; %[valuelistwidth 0 controlwidth figpos(4)];

%%plot spike and group lists and set the first group as selected
hspike = TextDisplayer(hmain, spikelistpos, eeg.general.eegID, 'EEG', 'normalized');
hgroup = TextDisplayer(hmain, grouplistpos, eegdata.grouplist.groupname, 'Group', 'normalized');
htext = getappdata(hspike, 'htext'); %every line (spike name) is a text object
ntext = numel(htext); %number of text lines displayed: not all spikes being displayed
%now re-set line (text) color for each spike (de-)selected
spikeselectindex = eegdata.grouplist.groupindex{1};
spikeselection = getappdata(hspike, 'selection'); spikeselection = 0*spikeselection; 
spikeselection(spikeselectindex) = ones(numel(spikeselectindex), 1);
for (i = 1:ntext)
    tagtag = get(htext(i), 'Tag');
    [str, rem] = strtok(tagtag, '_');
    linenum = str2num(str); %current line number that selected
    if (spikeselection(linenum) == 0)
        set(htext(i), 'Color', [0 0 0]);
    else
        set(htext(i), 'Color', [1 0 0]);
    end
end
setappdata(hspike, 'selection', spikeselection); setappdata(hspike, 'htext', htext);

%determine fields to display
fieldtitle = fieldnames(eeg); nfield = numel(fieldtitle);
sfwidth = fpos(3)/nfield; sfheight = fpos(4); sfbottom = fpos(2);
left = fpos(1);
for (i = 1:nfield)
    fieldpos = [left sfbottom sfwidth sfheight];
    subfield{i} = fieldnames(eeg.(fieldtitle{i}));
    hfield(i) = TextDisplayer(hmain, fieldpos, subfield{i}, fieldtitle{i}, 'normalized');
    left = left + sfwidth;
    if (~strcmp(fieldtitle{i}, 'general'))
       fieldselection = getappdata(hfield(i), 'selection');
       fieldselection = zeros(size(fieldselection));
       setappdata(hfield(i), 'selection', fieldselection);
       htext = getappdata(hfield(i), 'htext');
       if (~isempty(htext)) set(htext(1), 'Color', [0 0 0]); end
    end
end

%determine selected groups and selected fields
groupselection = getappdata(hgroup, 'selection');
groupselectindex = find(groupselection == 1);
numgroup = numel(groupselectindex);
for (xx = 1:numgroup) %if any groups selected
    groupvaluepos = [0.001 valuelistpos(2)+valuelistpos(4)/numgroup*(numgroup-xx) valuelistpos(3) valuelistpos(4)/numgroup-0.001];
    spikeselectindex = eegdata.grouplist.groupindex{groupselectindex(xx)};
    ncol = 1;
    valuetitle{ncol} = eegdata.grouplist.groupname{groupselectindex(xx)};
    valuecol{ncol} = eeg.general.eegID(spikeselectindex); %first column = selected spike names
    for (i = 1:nfield) %search selected sub-fields in each big field
        fieldselection = getappdata(hfield(i), 'selection');
        fieldselectindex = find( fieldselection == 1);
        for (j = 1:numel(fieldselectindex) )
            nj = fieldselectindex(j);
            ncol = ncol + 1;
            valuetitle{ncol} = subfield{i}{nj};
            kkk = numel(eeg.(fieldtitle{i}).(subfield{i}{nj}));
            for (mm = 1:numel(spikeselectindex))
                if (spikeselectindex(mm) > kkk)
                    valuecol{ncol}{mm} = [];
                elseif (iscell(eeg.(fieldtitle{i}).(subfield{i}{nj})))
                    wield = eeg.(fieldtitle{i}).(subfield{i}{nj}){spikeselectindex(mm)}; %this is a cell anyway
                    if ( isnumeric (wield) )
                       for (mnk = 1:numel(wield))
                           valuecol{ncol}{mm}{mnk} = num2str( wield(mnk) );
                       end
                    else
                       valuecol{ncol}{mm} = wield;     
                    end
                else
                    valuecol{ncol}{mm} = eeg.(fieldtitle{i}).(subfield{i}{nj})(spikeselectindex(mm));
                end
            end
        end
    end
    %re-arrange valuecol{ncol}{mm}{mnk} to fit textmsg{ncol}{mm}
    textmsg = rearrange(valuecol);
    TextDisplayer_multiple(hmain, groupvaluepos, textmsg, valuetitle, 'normalized');
end
%display values
setappdata(hmain, 'hspike', hspike); setappdata(hmain, 'hgroup', hgroup);
setappdata(hmain, 'hfield', hfield); %passing spikelist handle and field handles to other programs
setappdata(hmain, 'spikelistpos', spikelistpos); setappdata(hmain, 'grouplistpos', grouplistpos);
setappdata(hmain, 'valuelistpos', valuelistpos);
setappdata(hmain, 'fpos', fpos);
DataManager_eegcontrol_button(hmain, controlpos);

function textmsg = rearrange(valuecol);
%rearrnge valuecol{ncol}{mline}{nitem} to textmsg{nn}{nline}
ncol = numel(valuecol);
for (i = 1:ncol) %initial assignment
    textmsg{i} = [];
    nline(i) = 0;
    linenum(i) = numel(valuecol{i});
end
nspike = max(linenum); %all line number should be same(= spike number) for all columns, take the maximum anyway
for (i = 1:ncol)
for (j = 1:nspike)
    if ( j > numel(valuecol{i}) )
        itemnum(i,j) = 0;
    elseif (iscell(valuecol{i}{j}))
        itemnum(i,j) = numel(valuecol{i}{j});
    else 
        itemnum(i,j) = 1;
    end
end   
end
for (j = 1:nspike)
    maxitem(j) = max(itemnum(:,j));
end
for (i = 1:ncol)
    for (j = 1:nspike)
        for (k = 1:maxitem(j))
            nline(i) = nline(i) + 1;
            if (itemnum(i,j) == 0)
                textmsg{i}{nline(i)} = ' '; %stuffing
            elseif ( k > itemnum(i,j))
                textmsg{i}{nline(i)} = ' '; %stuffing
            elseif (iscell(valuecol{i}{j}))
                textmsg{i}{nline(i)} = valuecol{i}{j}{k};
            else
                textmsg{i}{nline(i)} = valuecol{i}{j};
            end
        end
    end
end