%%Resting state Dynamic Functional Connectivity by patient group
load('test_grp.mat');
for ID=1
    sub = sprintf('%s%s','sub-',test_grp(ID,1));
filename = sprintf('%s%s',sub,'_rest_timeseries.mat');
load(filename);
ts_corr = corr(ts); %ROIs x ROIs output matrix

%Perform analysis on the loaded data here
nPairs = 502; %number of unique pairs from template (ie., ROIs)
nTime = size(ts,1); %number of time points
nROI = 502 %number of regions of interest

%previous loading timeseries from above
window = 10; %variably change window from 5 - 20
mtd = coupling(ts,window,0,0); %time x nodes
% % flatten out approach for MTD
for nn=1:nROI
template = find(tril(ones(nROI))-eye(nROI)); %try taking lower triangle instead tril
end
nTime = size(mtd,3);
for tt = 1:nTime
    temp = mtd(:,:,tt);
    mtd_flat(:,tt) = temp(template); %flattens mtd across timepoints
end
% 
save(sprintf('%s%s',sub,'_rest_results.mat'), 'ts_corr', 'mtd', 'mtd_flat');
% end

%% Dynamic FC using Mac and Joni's code

%make a dtmx file that is filled with 1's and 0's
% Create a 140x1 matrix filled with 1's
dsmtx_final = ones(140,1);

% Save the matrix to a .mat file
save('dsmtx_final.mat', 'dsmtx_final');

load('dsmtx_final.mat')

%%
nodes = 502;
time = 140;
window = 10;
gamma = 1;

myFolder = 'C:\Users\lchur\functional_connectivity_analysis\follow-up';

filePattern = fullfile(myFolder, '*timeseries.mat');
theFiles = dir(filePattern);


% integration

[~,~,part,~,~] = integration(mtd,gamma);

storage_HC_ST_early(ii).part = part;




