function DataManager_Close_Callback

hf = gcbf;
%prompt for save option
option = questdlg('Save the database?', 'Save request', 'Yes');
if (strcmp(option, 'Yes')) %if save
    DataManager_Save_Callback;
end

if (~isempty(option)) && (~strcmp(option, 'Cancel'))
    if (isappdata(hf, 'pinfo')) rmappdata(hf, 'pinfo'); end
    if (isappdata(hf, 'eeg')) rmappdata(hf, 'eeg'); end
    if (isappdata(hf, 'behav')) rmappdata(hf, 'behav'); end
    if (isappdata(hf, 'data')) rmappdata(hf, 'data'); end
    if (isappdata(hf, 'bhdata')) rmappdata(hf, 'bhdata'); end
    if (isappdata(hf, 'eegdata')) rmappdata(hf, 'eegdata'); end
    if (isappdata(hf, 'plotparm')) rmappdata(hf, 'plotparm'); end
    if (isappdata(hf, 'dbtype')) rmappdata(hf, 'dbtype'); end
    if (isappdata(hf, 'celllistvar')) rmappdata(hf, 'celllistvar'); end
    if (isappdata(hf, 'celllisttitle')) rmappdata(hf, 'celllisttitle'); end
    
    hh = findobj(hf, 'Type', 'axes'); hhh = findobj(hf, 'Type', 'uicontrol');
    delete(hh); delete(hhh);
    [MCroot, MCname, DAname, DEname, ANname, MAname] = CurrentVersion; set(hf, 'Name', MAname);
end