function [eeg, eegdata] = DataManager_FindSpectralProp_old(eeg, eegdata, eegind, vv)
%%%do spectral analysis on broadband EEG files: (1) overall spectrum at different behavioral stages
%%%                                             (2) window-by-window spectrum
%%% require the following parameters
%%%    eeg.parm.specMinFreq/specMaxFreq; eeg.parm.specWinSize/specWinShift; eeg.parm.timeunit/buffersize;

%%% variables to assign
neeg = numel(eegind); 
if (~isfield(eeg, 'spec')) eeg.spec = []; end
if (~isfield(eegdata, 'spec')) eegdata.spec = []; end
if (~isfield(eeg.spec, 'sessPower')) eeg.spec.sessPower = cell(1, neeg); end
if (~isfield(eeg.spec, 'runPower')) eeg.spec.runPower = cell(1, neeg); end
if (~isfield(eeg.spec, 'stopPower')) eeg.spec.stopPower = cell(1, neeg); end
if (~isfield(eeg.spec, 'swsPower')) eeg.spec.swsPower = cell(1, neeg); end
if (~isfield(eeg.spec, 'remPower')) eeg.spec.remPower = cell(1, neeg); end

if (~isfield(eeg.spec, 'freq')) eeg.spec.freq = cell(1, neeg); end
if (~isfield(eegdata.spec, 'timewin')) eeg.spec.timewin = cell(1, neeg); end
if (~isfield(eegdata.spec, 'winpower')) eegdata.spec.winpower = cell(1, neeg); end

for (iiik = 1:neeg)
i = eegind(iiik);
if strcmp(eeg.parm.band{i}, 'broad') %%only do this for broadband EEG traces
    disp(['--------> spectral analysis: ', eeg.general.eegfile{i}]);
    minfy = eeg.parm.specMinFreq(i); maxfy = eeg.parm.specMaxFreq(i); fstep = eeg.parm.specFreqStep(i);
    buffsize = eeg.parm.buffersize(i);
    ws = eeg.parm.specWinSize(i); wshift = eeg.parm.specWinShift(i); timeunit = eeg.parm.timeunit(i);
    F = minfy:fstep:maxfy; %%default frequency range, may be changed 
    %%%%get the EEg data
    [timestamp, dat, gain, fs] = ReadEEGFile(eeg.general.eegfile{i}, timeunit, buffsize); % timestamp in second, dat in mV 
    %%%%power spectrum density window by window 
    %%%%%%%%%need to compute windowsize and shiftpoint in number of datapoints
    dat = dat-mean(dat);
    %ws = round(ws*fs); tnshift = round(wshift*fs); noverlap = ws-tnshift; nws = round(numel(dat)/tnshift);
    %[SSS, fy, ptime, psd] = spectrogram(dat, ws, noverlap, F, fs); %%%default window = Hamming wondow 256 points
    %    %% psd in mV/Hz, ptime start from 0 in seconds; psd(fy, time); psd = simple periodogram
    nws = ceil((timestamp(numel(timestamp)) - timestamp(1))/wshift);
    psd = zeros(numel(F), nws);
    for (tti = 1:nws)
        sT = timestamp(1) + tti*wshift - ws/2; eT = timestamp(1) + tti*wshift + ws/2;
        psd(:,tti) = findmtpower(dat, timestamp, F, fs, sT, eT); 
        %if (~isempty(ppnow)) psd(:,tti) = ppnow; end
    end
    %eeg.spec.timewin{i} = ptime+ timestamp(1); 
    eeg.spec.timewin{i} = timestamp(1)+(1:nws)'*wshift; 
    eeg.spec.freq{i} = F; eegdata.spec.winpower{i} = abs(psd); %psd[fy, ntime]
    
    %%%%power at different behavioral stages
    eeg.spec.sessPower{i} = findmtpower(dat, timestamp, F, fs, timestamp(1), timestamp(numel(timestamp))); 
    evTime = eegdata.event.eventtimes{i}; evType = eeg.parm.eventtype{i};
    runpower = zeros(size(F)); stoppower = zeros(size(F));  swspower = zeros(size(F)); rempower = zeros(size(F));
    irun = 0; istop = 0; isws = 0; irem = 0;
    for (j = 1:numel(evTime))
        if (strcmp(evType{j}, 'run')) 
            irun = irun + 1; runpower = runpower + findmtpower(dat, timestamp, F, fs, evTime{j}.start, evTime{j}.ent); 
        end
        if (strcmp(evType{j}, 'stop')) 
            istop = istop + 1; stoppower = stoppower + findmtpower(dat, timestamp, F, fs, evTime{j}.start, evTime{j}.ent); 
        end   
        if (strcmp(evType{j}, 'sws')) 
            isws = isws + 1; swspower = swspower + findmtpower(dat, timestamp, F, fs, evTime{j}.start, evTime{j}.ent); 
        end
        if (strcmp(evType{j}, 'rem')) 
            irem = irem + 1; rempower = rempower + findmtpower(dat, timestamp, F, fs, evTime{j}.start, evTime{j}.ent); 
        end
    end
    if (irun>0) eeg.spec.runPower{i} = runpower/irun; end
    if (istop>0) eeg.spec.stopPower{i} = stoppower/istop; end  %%%if only one ev, mean only yield one number!
    if (isws>0) eeg.spec.swsPower{i} = swspower/isws; end
    if (irem>0) eeg.spec.remPower{i} = rempower/irem; end
end
end

function power = findmtpower(dat, timestamp, F, fs, startT, endT) %multi-taper power spectral density estimate
power = zeros(size(F)); pp = zeros(size(F)); nev = 0; leng = fs*10; shift = fs*5; %do it for every 10s data and shifted by every 5 s
for (k = 1:numel(startT))
     iii = find( (timestamp>=startT(k)) & (timestamp<endT(k)) ); np = numel(iii);
     nwin = ceil(np/shift); datnow = dat(iii) - mean(dat(iii));
     for (i = 1:nwin)
         datok = datnow( max([1 i*shift-round(leng/2)]) : min([numel(datnow) i*shift+round(leng/2)]) ); datok = datok - mean(datok);
         if (numel(datok) >= 0.5*fs) %a hard parameter imposed here for the minimum number of data points ~500ms EEg data
             nev = nev + 1; nft = 2^nextpow2(numel(datok));
             [ppnow,freq] = pmtm(datok, 4, nft, round(fs)); 
             pp = pp + sumuppower(freq, ppnow, F);
         end
     end
end
if (nev >0) power = pp/nev; end

% function power = findpower(dat, timestamp, F, fs, startT, endT) %regular fft power spectral density estimate
% power = zeros(size(F)); pp = zeros(size(F)); nev = 0; leng = fs*10; shift = fs*5;
% for (k = 1:numel(startT))
%      iii = find( (timestamp>=startT(k)) & (timestamp<endT(k)) ); np = numel(iii); 
%      nwin = ceil(np/shift); datnow = dat(iii) - mean(dat(iii));
%      %disp(['------> timestamp (first, last): ', num2str(timestamp(1)), '  ', num2str(timestamp(numel(timestamp)))]);
%      %disp(['------> start/end time: ', num2str(startT(k)), '  ', num2str(endT(k))]);
%      for (i = 1:nwin)
%          %disp('-----> get in now');
%          datok = datnow( max([1 i*shift-round(leng/2)]) : min([numel(datnow) i*shift+round(leng/2)]) ); datok = datok - mean(datok);
%          if (numel(datok) >= 0.5*fs) %a hard parameter imposed here for the minimum number of data points ~500ms EEg data
%              nev = nev + 1; nft = 2^nextpow2(numel(datok));
%              freq = fs/2*linspace(0,1,nft/2+1); fftnow = fft(datok, nft)/np; ppnow = 2*abs(fftnow(1:nft/2+1));  %a regular fft
%              pp = pp + sumuppower(freq, ppnow, F);
%          end
%      end
%      
% %      if (np >= 100) %a hard parameter imposed here for the minimum number of data points ~50ms EEg data
% %          datnow = dat(iii) - mean(dat(iii)); nev = nev + 1; nft = min([4096 2^nextpow2(np)]);
% %          freq = fs/2*linspace(0,1,nft/2+1); fftnow = fft(datnow, nft)/np; ppnow = 2*abs(fftnow(1:nft/2+1));  %a regular fft
% %          pp = pp + sumuppower(freq, ppnow, F);
% %      end
% end
% if (nev >0) power = pp/nev; end

function power = sumuppower(freq, pp, F)
power = zeros(size(F)); %%%need to find the corresponding interval between freq and FF and do the normalization if necessary
fgap = freq(2)-freq(1); Fgap = F(2)-F(1);
for (i = 1:numel(F))
    kk = find( (freq>=F(i)) & (freq<F(i)+Fgap) ); 
    if (~isempty(kk))
        power(i) = mean(pp(kk))*numel(kk)*fgap/Fgap; %%%make sure total power conserved for the same frequency bin
    end
end





