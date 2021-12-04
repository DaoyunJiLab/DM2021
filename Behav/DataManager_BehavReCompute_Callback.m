function DataManager_BehavReCompute_Callback
%%Add or change variable values of the selected group, and save into a new database

hf = gcbf; pinfo = getappdata(hf, 'behav'); data = getappdata(hf, 'bhdata'); tagmark = get(gcbo, 'Tag');
hgroup = getappdata(hf, 'hgroup');hfield = getappdata(hf, 'hfield');
ok = 1; plotparm = getappdata(hf, 'plotparm'); vv = plotparm.showdetail; %if vv=1, shoe details of computation
cc = plotparm.compute; %if cc=0, assign parameters for computation only; if 1, do the real computation
ow = plotparm.overwrite; %if ow=1, plot into the current database
[writefilename, okk] = getoutputfile(hf, ow);

%get selected cellind
groupselection = getappdata(hgroup, 'selection'); cellind = []; grpind = find(groupselection == 1); 
for (kk = 1:numel(grpind)) cellind = union(cellind, data.grouplist.groupindex{grpind(kk)}); end
if (~isempty(cellind)) & okk
    if (strcmp(tagmark, 'recomputeinit'))
        [pinfo,data] = DataManager_BehavComputeInit_Callback(pinfo,data, cellind, vv);
    elseif (strcmp(tagmark, 'recomputeall'))
        [pinfo,data] = DataManager_BehavComputeAll_Callback(pinfo,data, cellind, vv);
    end
    %%%%%%%%save and plot the new database
    if (~isempty(pinfo.general.finaldir))
        if (ok)
            behav = pinfo; bhdata = data; save(writefilename, 'behav', 'bhdata');
            if (ow)
               iii = get(hf, 'Children'); delete(iii(~strcmp(get(iii, 'Type'), 'uimenu'))); 
               fname = get(hf, 'Name'); ftt = strfind(fname, '__'); 
               newname = fname(1:ftt-1); set(hf, 'Name', newname); hmain = hf; 
            else
               hmain = DataManager_DataManager_Callback; plotparm.linkspike = 0; plotparm.linkeeg = 0; plotparm.linkbehav = 1;
            end
            setappdata(hmain,'plotparm', plotparm);
            %DataManager_PlotBehavDatabase(hmain, behav, bhdata);
            DataManager_PlotSpikeDatabase(hmain, behav, bhdata, 'Sessions', 'sessID', '.behavdb');
            set(hmain, 'Name', strcat(get(hmain, 'Name'), '__', writefilename));
            setappdata(hmain, 'behav', behav); setappdata(hmain, 'bhdata', bhdata); pinfo = []; data = [];
        end
    else
        disp('-------------> no cells in the database!');
    end
else
    disp('--------------> no groups selected, groups do not contain any cells, or action cancelled');
end
disp('**********************');

function [writefilename, okk] = getoutputfile(hf, ow)
okk = 1; writefilename = [];
if (ow == 0)
   [fname, pname] = uiputfile(fullfile(cd, '*.behavdb'), 'Write the new spike database to:');
   if (numel(fname)>1)
      writefilename = fullfile(pname, fname);
   else
      okk = 0;
   end
else
   %input = questdlg('The current database will be altered and overwritten. Are you sure?', 'Overwrite?', 'Yes');
   %if (strcmp(input, 'Yes'))
      fname = get(hf, 'Name'); ftt = strfind(fname, '__'); writefilename = fname(ftt+2:numel(fname));
   %else
   %   okk = 0;
   %end
end
