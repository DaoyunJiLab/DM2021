function DataManager_LinkDatabase_Callback(hf, plotparm, tagmark)
%%cross-link spike/behav/eeg database
if (strcmp(tagmark, 'linkspike'))
    if (plotparm.linkspike == 1)
        spn = getappdata(hf, 'spikedbname'); strnow = ['Linked to ', spn, '. Are you sure to remove?'];
        input = questdlg(strnow, 'Link spike database', 'No');
        if (strcmp(input, 'Yes'))
             rmappdata(hf, 'pinfo'); rmappdata(hf, 'data'); rmappdata(hf, 'spikedbname'); 
             plotparm.linkspike = 0;
             set(findobj(hf, 'Tag', tagmark), 'String', 'NoSpikedb', 'ForeGroundColor', [0.2 0.2 0.2]);
        end
        setappdata(hf, 'plotparm', plotparm);
    else
        LinkSpikeDB(hf, plotparm, tagmark); 
    end
elseif (strcmp(tagmark, 'linkeeg'));
    if (plotparm.linkeeg == 1)
        spn = getappdata(hf, 'eegdbname'); strnow = ['Linked to ', spn, '. Are you sure to remove?'];
        input = questdlg(strnow, 'Link eeg database', 'No');
        if (strcmp(input, 'Yes'))
             rmappdata(hf, 'eeg'); rmappdata(hf, 'eegdata'); rmappdata(hf, 'eegdbname'); 
             plotparm.linkeeg = 0;
             set(findobj(hf, 'Tag', tagmark), 'String', 'NoEEGdb', 'ForeGroundColor', [0.2 0.2 0.2]);
        end
        setappdata(hf, 'plotparm', plotparm);
    else
        LinkEEGDB(hf, plotparm, tagmark); 
    end
elseif (strcmp(tagmark, 'linkbehav'))
    if (plotparm.linkbehav == 1)
        spn = getappdata(hf, 'bhdbname'); strnow = ['Linked to ', spn, '. Are you sure to remove?'];
        input = questdlg(strnow, 'Link behav database', 'No');
        if (strcmp(input, 'Yes'))
             rmappdata(hf, 'behav'); rmappdata(hf, 'bhdata'); rmappdata(hf, 'bhdbname'); 
             plotparm.linkbehav = 0;
             set(findobj(hf, 'Tag', tagmark), 'String', 'NoBehavdb', 'ForeGroundColor', [0.2 0.2 0.2]);
        end
        setappdata(hf, 'plotparm', plotparm);
    else
        LinkBehavDB(hf, plotparm, tagmark); 
    end
end

function LinkSpikeDB(hf, plotparm, tagmark)
ok = 1;
if (isappdata(hf, 'data'))
    datanow = getappdata(hf, 'data');
    if (~isempty(datanow))
        spn = getappdata(hf, 'spikedbname'); 
        strnow = ['Already linked to ', spn, '. Replace with another spike database?'];
        input = questdlg(strnow, 'Link spike database', 'No');
        if ~(strcmp(input, 'Yes')) ok = 0; end
    end
end
if (ok)
     [fname, pname] = uigetfile(fullfile(cd, '*.spikedb'), 'Select a spike database to link:');
     if (numel(fname)>1)
           filename = fullfile(pname, fname);
           S = load(filename, '-mat'); %load the file
           pinfo = S.pinfo; data = S.data; S = []; %load the plot structure
     else
           ok = 0;
     end
end
if (ok)
     setappdata(hf, 'pinfo', pinfo); setappdata(hf, 'data', data); setappdata(hf, 'spikedbname', filename);
     plotparm.linkspike = 1; setappdata(hf, 'plotparm', plotparm);
     set(findobj(hf, 'Tag', tagmark), 'String', 'SpikeLinked', 'ForegroundColor', [1 0 0]);
end

function LinkEEGDB(hf, plotparm, tagmark)
ok = 1;
if (isappdata(hf, 'eegdata'))
    datanow = getappdata(hf, 'eegdata');
    if (~isempty(datanow))
        spn = getappdata(hf, 'eegdbname'); 
        strnow = ['Already linked to ', spn, '. Replace with another EEG database?'];
        input = questdlg(strnow, 'Link eeg database', 'No');
        if ~(strcmp(input, 'Yes')) ok = 0; end
    end
end
if (ok)
     [fname, pname] = uigetfile(fullfile(cd, '*.eegdb'), 'Select a eeg database to link:');
     if (numel(fname)>1)
           filename = fullfile(pname, fname);
           S = load(filename, '-mat'); %load the file
           eeg = S.eeg; eegdata = S.eegdata; S = []; %load the plot structure
     else
           ok = 0;
     end
end
if (ok)
     setappdata(hf, 'eeg', eeg); setappdata(hf, 'eegdata', eegdata); setappdata(hf, 'eegdbname', filename);
     plotparm.linkeeg = 1; setappdata(hf, 'plotparm', plotparm);
     set(findobj(hf, 'Tag', tagmark), 'String', 'EEGLinked', 'ForegroundColor', [1 0 0]);
end

function LinkBehavDB(hf, plotparm, tagmark)
ok = 1;
if (isappdata(hf, 'bhdata'))
    datanow = getappdata(hf, 'bhdata');
    if (~isempty(datanow))
        spn = getappdata(hf, 'bhdbname'); 
        strnow = ['Already linked to ', spn, '. Replace with another behav database?'];
        input = questdlg(strnow, 'Link behav database', 'No');
        if ~(strcmp(input, 'Yes')) ok = 0; end
    end
end
if (ok)
     [fname, pname] = uigetfile(fullfile(cd, '*.behavdb'), 'Select a behav database to link:');
     if (numel(fname)>1)
           filename = fullfile(pname, fname);
           S = load(filename, '-mat'); %load the file
           behav = S.behav; bhdata = S.bhdata; S = []; %load the plot structure
     else
           ok = 0;
     end
end
if (ok)
     setappdata(hf, 'behav', behav); setappdata(hf, 'bhdata', bhdata); setappdata(hf, 'bhdbname', filename);
     plotparm.linkbehav = 1; setappdata(hf, 'plotparm', plotparm);
     set(findobj(hf, 'Tag', tagmark), 'String', 'BehavLinked', 'ForegroundColor', [1 0 0]);
end