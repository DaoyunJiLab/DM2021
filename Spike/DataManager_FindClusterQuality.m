function [pinfo,data] = DataManager_FindClusterQuality(pinfo,data, cellind, vv)
%Assess quality of clusters
%    assuming the clusters are separated from each other (therefore separation is not considered here)
%cuttoff index: measure how many data points are mising (usually by imposing tight bounds)
%fit index: measure how good the cluster density fits a normal distribution (and thus indirectly measures the contamination)

%strategy: (1) right now simply taking all amplitudes of the 4 channels as 4 sets of samples;
%          (2) compute the histogram curve of each channel;
%          (3) fit the curve with normal PDF (not fit the samples with a normal PDF).
%          (4) compute the goodness of fit and percentage of missing points

%[[other potential algorithms:
%          (1) follow the same strategy as above, but fit multiple dimensions at the same time.
%          (2) fit the multiple-dimension samples with multi-variate normal distribution (using the matlab gmdistribution.fit)
%                -----but this has the problem of fitting samples rather than histograms
%          (3) this whole issue can be systematically investigated by trying different strategies and computer simulations]]

%require the following parameters
%   pinfo.parm.clusterMinDev, pinfo.parm.clusterNumBin
%require the following variables
%   pinfo.general.parmfile

%variable to assign
disp('Quantify cluster sorting quality ....');
if (vv ==1) && (numel(cellind)> 5) disp('-----> too many clusters selected; only showing the first 5 clusters'); end
if (~isempty(cellind))
   nspike = numel(pinfo.general.parmfile);
   if (~isfield(pinfo.general, 'clustNpoints')) pinfo.general.clustNpoints = cell(1, nspike); end
   if (~isfield(pinfo.general, 'clustChFitErr')) pinfo.general.clustChFitErr = cell(1, nspike); end % vector, each between [0 1], NaN if undecided
   if (~isfield(pinfo.general, 'clustChCutoffI')) pinfo.general.clustChCutoffI = cell(1, nspike); end % vector, each between [0 1], NaN if undecided
   if (~isfield(pinfo.general, 'clustMaxFitErr')) pinfo.general.clustMaxFitErr = cell(1, nspike); end % scalar, each between [0 1], NaN if undecided
   if (~isfield(pinfo.general, 'clustMaxCutoffI')) pinfo.general.clustMaxCutoffI = cell(1, nspike); end % scalar, each between [0 1], NaN if undecided
   if (~isfield(pinfo.general, 'clustIsolDist')) pinfo.general.clustIsolDist = cell(1, nspike); end % scalar, NaN if undecided
   if (~isfield(pinfo.general, 'clustLratio')) pinfo.general.clustLratio = cell(1, nspike); end
   
   cl0parm = []; cl0file = [];
   options = optimset('TolFun', 0.0001, 'MaxIter', 10000); %, 'Display', 'final');
   for (jjjk = 1:numel(cellind))
        i = cellind(jjjk); if (jjjk>5) vv = 0; end %%%turn off showing if more than 5 cells selected
        %%%%%get thresholds (these parameters should be the same for all clusters to ensure a valid comparison across the clusters)
       disp(strcat('-----> assess cluster quality ---', pinfo.general.parmfile{i}));
       mindev = []; numbin = []; minN = [];
       if (isfield(pinfo.parm, 'clustMinDev')) && (~isnan(pinfo.parm.clustMinDev(i)))
           mindev = pinfo.parm.clustMinDev(i);
       end
       if (isfield(pinfo.parm, 'clustNumBin')) && (~isnan(pinfo.parm.clustNumBin(i)))
           numbin = pinfo.parm.clustNumBin(i);
       end
       if (isfield(pinfo.parm, 'clustMinNpoints')) && (~isnan(pinfo.parm.clustMinNpoints(i)))
           minN = pinfo.parm.clustMinNpoints(i);
       end
       if (isempty(mindev)) || (isempty(numbin)) || (isempty(minN))
           disp('----------> one or more required parameters are not assigned');
       else
            parm = ReadSpikeParam(pinfo.general.parmfile{i}); fitI = NaN*ones(1,4); cutoffI = NaN*ones(1,4);
            [nsss, nppp] = size(parm);
            if (nsss>0) && (nppp>=5)
               for (j = 1:4)
                [fitI(j),cutoffI(j)] = assessclusterquality(parm(:,j+1), mindev, numbin, minN, options, pinfo.general.parmfile{i}, j, vv);
                pinfo.general.clustChFitErr{i}(j) = fitI(j);
                pinfo.general.clustChCutoffI{i}(j) = cutoffI(j);
               end
            else
                disp(['-----------> warning: no spikes in : ', pinfo.general.parmfile{i}]);
            end
            pinfo.general.clustMaxFitErr{i} = max(fitI); pinfo.general.clustMaxCutoffI{i} = max(cutoffI);
            pinfo.general.clustNpoints{i} = nsss;
            %%%%%%find isolation distance and Lratio
          if nsss>0 
            [pp, ~, ~] = fileparts(pinfo.general.parmfile{i});
            cl0filenow = fullfile(pp, 'cl-0.spm');
            if ~strcmp(cl0filenow, cl0file)
                cl0file = cl0filenow; 
                if (exist(cl0file, 'file') == 2)
                    cl0parm = ReadSpikeParam(cl0file, 1);
                else
                    cl0parm = []; 
                end
            end
            if (~isempty(cl0parm))
 
                gainnow = pinfo.general.gain{i}; 
                if ~isempty(gainnow)
                    Fet = cl0parm(:,2:5); clustFet = parm(:,2:5);
                    for (j = 1:4)
                        Fet(:,j) = Fet(:,j)*gainnow(j); clustFet(:,j) = clustFet(:,j)*gainnow(j);
                    end
                    [iso, Lr] = ComputeIsolationDistanceLratio_FromMClust(Fet, clustFet);
                    pinfo.general.clustIsolDist{i} = iso; pinfo.general.clustLratio{i} = Lr;
                else
                    disp(['-----------> warning: amplitude gain not assigned: ', pinfo.general.parmfile{i}]);
                end
            else
                disp(['-----------> warning: cl-0 file not exist: ', cl0filenow]);
            end
          end
       end
   end
end

function [fitI,cutoffI] = assessclusterquality(parm, mindev, numbin, minN, options, parmfile, chnum, vv)
fitI = NaN; cutoffI = NaN; ss = std(parm); mmm = mean(parm); nmax = numel(find(parm == max(parm)))/numel(parm);
if (numel(parm) >= minN) && (ss>=mindev) && (nmax < 0.05) %the last to exclude the saturated clusters
   maxV = min([max(parm) mmm+4*ss]); minV = max([min(parm) mmm-4*ss]); %to prevent single low/high value far far away
   %%%%here binsize is automatically computed: 
   %%%% binsize = 3.491*sigma*N(-1/3)
   stdnow = std(parm); nnow = numel(parm); binsize = 3.491*stdnow*(nnow^(-1/3));   %binsize = (maxV-minV)/numbin; 
   binnow = minV:binsize:maxV; numbin = numel(binnow);
   %disp(['-----> new bin number: ', num2str(numbin)]);
   %%set initial values for data fitting: this already taking account the normalization below
   count = histc(parm, binnow); 
   %%%%histogram needs to be normalized and smoothed
   sigma = 2; nnt = 5; ttt = [-nnt*sigma:1:nnt*sigma]; wind = normpdf(ttt, 0, sigma); wind = wind/sum(wind); nw = numel(wind); %will have a delay = (n-1)/2 bins
   count = conv(count, wind); count = count((nw-1)/2+1:numel(count)-(nw-1)/2);
   % This is what should be used for normalization
   normfactor = 1/sum(count)/binsize; addconstant = 1; %%For these normalization factors: binsize*sum(count) = 1; 
   %%%%%%%% To save computing time, additional constant factor is applied
      addconstant = binsize*numbin; %%%% For this addional factor: sum(count) = numbin; average value = 1 (to make the fitting computation faster)
      normfactor= normfactor*addconstant;
   %%%%%%%%%%%%%%%%%
   count = count*normfactor;
   [sx,sy] = size(count); if (sx>sy) count = count'; end %%%%all row vectors
   %%%%now do the fitting
   [maxnow, II] = max(count); xx = binnow-binnow(II); %(1:numel(count))-II; 
   %start = [stdnow maxnow 0]; %%% somehow using maxnow instead of nnow does not converge - stop too early;;;%start = [20 numbin]; %start = [numbin*ss/sum(count) 1];
   start = [stdnow addconstant*nnow/binsize/numbin 0];
   [fitpar, fval, exitflag] = fminsearch(@(par)Npdffit(par, xx, count), start, options); 
   sigma = fitpar(1); a = fitpar(2); c = fitpar(3); 
   %npdfnow = a/(sigma*sqrt(2*pi))*exp(-xx.*xx/(2*sigma*sigma));
   npdfnow = a/(sigma*sqrt(2*pi))*exp(-(xx-c).*(xx-c)/(2*sigma*sigma));
   if (exitflag == 1)
       %fitI = fval; %%%fval is the squared minimum error
       %fitI = binsize*sqrt(fval/numbin)/abs(a); %%%% 10*normalized (/a) error area (error x binsize) per point; a is the fittng amplitude of normal distribution       
       fitI = binsize*fval/abs(a); %%%%index is the total error (given original data area = 1); fval = total error; a = fittng amplitude of normal distribution
       %fitI = addconstant*fitI;   %%%%% This is not needed: all normfactors
       %                                 including addconstant are already absorbed in the amplitude a,
       %                                 which is the integrated area of count over xx
       cutoffI = 1 - (binsize*sum(abs(npdfnow/a))); %%%best fit sigma and a could be both negative
   else
       disp(['---------> normal distribution fit did not converge. Error now: ', num2str(fval)]);
   end
   %if (vv==1) plotfitcurve(xx, count, npdfnow, sigma, a, 0, fitI, cutoffI, parmfile, chnum, maxV, minV, max(parm), min(parm), ss, binsize, numbin); end 
   if (vv==1) plotfitcurve(xx+binnow(II), count/normfactor, npdfnow/normfactor, sigma, a/normfactor, c+binnow(II), fitI, cutoffI, parmfile, chnum, maxV, minV, max(parm), min(parm), ss, binsize, numbin); end 
else
   if (vv== 1)
        str = ['---------> ', parmfile, ': ch#', num2str(chnum), ' does not meet with minimum requirement for cluster quality assessment'];
        msgbox(str, 'Warning');
   end
end

function err = Npdffit(par, x, count) %%%total error: sum of absolute distance of all points
sigma = par(1); a = par(2); c = par(3);
%nx = numel(x); %err= 0;
%npdf = a/(sigma*sqrt(2*pi))*exp(-x.*x/(2*sigma*sigma));
npdf = a/(sigma*sqrt(2*pi))*exp(-(x-c).*(x-c)/(2*sigma*sigma));
err = sum(abs(count - npdf));
%err = sum((count - npdf).^2); 
%err = sum(err)/nx;

function plotfitcurve(xx, count, npdfnow, sigma, a, c, fitI, cutoffI, parmfile, chnum, maxV, minV, Rmax, Rmin, ss, binsize, numbin)
hf = figure('Name', strcat(parmfile, '_Ch', num2str(chnum), '_normalfit')); hax = axes('Parent', hf, 'NextPlot', 'add'); 
line(xx, count, 'Parent', hax, 'LineStyle', 'none', 'Marker', '.', 'Markersize', 20);
line(xx, npdfnow, 'Parent', hax); 
text('String', strcat('[fitI cutoffI] = ', num2str([fitI cutoffI])), 'Parent', hax, 'Units', 'normalized', 'Position', [0 0.96]);
text('String', strcat('normfun [Sigma Amp Center] = ', num2str([sigma a/sqrt(2*pi)/sigma c])), 'Parent', hax, 'Units', 'normalized', 'Position', [0 0.92]);
text('String', strcat('data [min max std] = ', num2str([Rmin Rmax ss])), 'Parent', hax, 'Units', 'normalized', 'Position', [0 0.88]);
text('String', strcat('Used [minV maxV]=', num2str([minV maxV]), '; [binsize binnum]=', num2str([binsize numbin])), 'Parent', hax, 'Units', 'normalized', 'Position', [0 0.84]);










