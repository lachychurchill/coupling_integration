%% Dynamic FC analysis

nodes = 502;
time = 140;
window = 15;
gamma = 1;

myFolder = 'V:\fMRI_analysis\functional_connectivity_analysis\follow-up';

filePattern = fullfile(myFolder, '*timeseries.mat');
theFiles = dir(filePattern);

storage_followup(length(theFiles)) = struct('name',1,'ts',1);

for ii = 1:length(theFiles)
    baseFileName = theFiles(ii).name;
    fullFileName = fullfile(theFiles(ii).folder, baseFileName);
    ts = load(fullFileName);
    storage_followup(ii).ts = ts.ts;
    storage_followup(ii).name = extractBefore(baseFileName, '_timeseries.mat');
end

% load('sub-CON05_ses-01_run-1_timeseries.mat');
%% MTD values
for k = 1:15
    storage_followup(k).mtd = coupling(storage_followup(k).ts(:,:),window,0,1)
end
%%
% integration
for k = 1:2
   sub = sprintf(storage_followup(k).name);
   [CA,~,~,~,~,~] = integration(storage_followup(k).mtd,gamma);
   storage_followup(k).CA = CA;
   [~,modularity,~,~,~,~] = integration(storage_followup(k).mtd,gamma);
   storage_followup(k).modularity = modularity;
   [~,~,part,~,~,~] = integration(storage_followup(k).mtd,gamma);
   storage_followup(k).part = part;
   [~,~,~,zscore,~,~] = integration(storage_followup(k).mtd,gamma);
   storage_followup(k).zscore = zscore;
   [~,~,~,~,cartop,~] = integration(storage_followup(k).mtd,gamma);
   storage_followup(k).cartop = cartop;
   [~,~,~,~,~,idx] = integration(storage_followup(k).mtd,gamma);
   storage_followup(k).idx = idx;
   save(sprintf('%s%s',sub,'integration_measures'), 'CA', 'modularity', 'part', 'zscore', 'cartop','idx')
end

%% K-means clustering
xbins = [0:0.01:1]; ybins = [5:-.1:-5];
CP = zeros(size(xbins,2),size(ybins,2),time);
xNumBins = numel(xbins); yNumBins = numel(ybins);

for k = 1:40
    sub = sprintf(storage_followup(k).name);
    idx = kmeans(reshape(storage_followup(k).cartop,xNumBins * yNumBins,time)',2);
    storage_followup(k).kmean = idx;
    save(sprintf('%s%s',sub,'kmean.mat'), 'idx');
end

%%

%%
% Assume storage_followup contains dynamic connectivity data for all patients




