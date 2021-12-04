function DataManager_UpdateDisplay(hf)
%%update value field in a spike database browser
%%update bottom value frame according to new selections of spikes and fields
%%not update calculations

dbtype = getappdata(hf, 'dbtype'); celllistvar = getappdata(hf, 'celllistvar');
%%%%%%%added options for all types of databases
pinfostr = 'pinfo'; datastr = 'data';
switch dbtype
    case '.eegdb'
        pinfostr = 'eeg'; datastr = 'eegdata';
    case '.behavdb'
        pinfostr = 'behav'; datastr = 'bhdata';
end

pinfo = getappdata(hf, pinfostr); data = getappdata(hf, datastr);
hspike = getappdata(hf, 'hspike'); hgroup = getappdata(hf, 'hgroup');
hfield = getappdata(hf, 'hfield');
groupselection = getappdata(hgroup, 'selection'); %spikeselection = getappdata(hspike, 'selection'); 
valuelistpos = getappdata(hf, 'valuelistpos');
plotparm = getappdata(hf, 'plotparm');

%remove the old values and old title, old bars, not spike/field display
delete( findobj(hf, 'Tag', 'multipleaxes') );
delete( findobj(hf, 'Tag', 'multipletitle') );
delete( findobj(hf, 'Tag', 'multiplebar') );

%%%first update the spike selections
spikeselection = getappdata(hspike, 'selection'); spikeselection = 0*spikeselection; 
groupsel = find(groupselection == 1); spikeselectindex = [];
for (i = 1:numel(groupsel))
    spikeselectindex = union(spikeselectindex, data.grouplist.groupindex{groupsel(i)});
end
spikeselection(spikeselectindex) = ones(numel(spikeselectindex), 1);
htext = getappdata(hspike, 'htext'); %every line (spike name) is a text object
for (i = 1:numel(htext))
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

%%%%%temporary to show event data
%for (i = 1:numel(spikeselection))
%    if (spikeselection(i))
%        disp( [size(pinfo.general.eventname{i}) size(data.event.eventtimes{i})] );
%    end
%end

%determine all possible field
fieldtitle = fieldnames(pinfo); nfield = numel(fieldtitle);
for (i = 1:nfield) subfield{i} = fieldnames(pinfo.(fieldtitle{i})); end
%determine selected groups and selected fields
numgroup = numel(groupsel);
for (xx = 1:numel(groupsel)) %if any groups selected
    groupvaluepos = [0.0001 valuelistpos(2)+valuelistpos(4)/numgroup*(numgroup-xx) valuelistpos(3) valuelistpos(4)/numgroup-0.0001];
    spikeselectindex = data.grouplist.groupindex{groupsel(xx)}; ncell = numel(spikeselectindex);
    valuecol = []; valuetitle = []; ncol = 1;
    valuetitle{ncol} = strcat(data.grouplist.groupname{groupsel(xx)}, '(N=', num2str(ncell), ')');
    firstvar = cell(1,numel(spikeselectindex));
    for (ttj = 1:numel(spikeselectindex))
        firstvar{ttj} = strcat(num2str(spikeselectindex(ttj)), '=', pinfo.general.(celllistvar){spikeselectindex(ttj)});
    end
    valuecol{ncol} = firstvar;
    %valuecol{ncol} = pinfo.general.(celllistvar)(spikeselectindex); %first column = selected spike names
    if (plotparm.showdetail) %%%this add an option for noshow, if too may values to display
    for (i = 1:nfield) %search selected sub-fields in each big field
        fieldselection = getappdata(hfield(i), 'selection');
        fieldselectindex = find( fieldselection == 1);
        for (j = 1:numel(fieldselectindex) )
            nj = fieldselectindex(j);
            ncol = ncol + 1;
            valuetitle{ncol} = subfield{i}{nj};
            kkk = numel(pinfo.(fieldtitle{i}).(subfield{i}{nj}));
            for (mm = 1:numel(spikeselectindex))
                if (spikeselectindex(mm) > kkk)
                    valuecol{ncol}{mm} = [];
                elseif (iscell(pinfo.(fieldtitle{i}).(subfield{i}{nj})))
                    wield = pinfo.(fieldtitle{i}).(subfield{i}{nj}){spikeselectindex(mm)}; %this is a cell anyway
                    if ( isnumeric (wield) )
                       for (mnk = 1:numel(wield))
                           valuecol{ncol}{mm}{mnk} = num2str( wield(mnk) );
                       end
                    else
                       valuecol{ncol}{mm} = wield;     
                    end
                else
                    valuecol{ncol}{mm} = pinfo.(fieldtitle{i}).(subfield{i}{nj})(spikeselectindex(mm));
                end
            end
        end
    end
    end
    %re-arrange valuecol{ncol}{mm}{mnk} to fit textmsg{ncol}{mm}
    textmsg = rearrange(valuecol);
    TextDisplayer_multiple(hf, groupvaluepos, textmsg, valuetitle, 'normalized');
end

function textmsg = rearrange(valuecol)
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