%% LOAD IN THE DATA TO A STORAGE
% BASELINE 
% Specify the folder where the files lives.

baseline_Folder = 'V:\fMRI_analysis\functional_connectivity_analysis\baseline_rbd';



% Get a list of all files in the folder with the desired file name pattern.

filePatternbaseline = fullfile(baseline_Folder, '*irbd*'); % Looking for the 400x400 matrices

theFilesbaseline = dir(filePatternbaseline); % Storing all filenames
% Load in FC matrix baseline

% load in data

storage(length(theFilesbaseline)) = struct('name',1,'ts_baseline',1,'corr_ts_baseline',1,'name2',1,'ts_followup',1,'corr_ts_followup',1);

for k = 1:16
    baseFileName = theFilesbaseline(k).name;
    fullFileName = fullfile(theFilesbaseline(k).folder, baseFileName);
    thisArray{k} = load(fullFileName); % Should have 54 structures containing a matrix
    Array = load(fullFileName);
    cells = Array.ts(:, :);
    subjectname = strrep(baseFileName,'timeseries.mat','');
    subjectname = extractBefore(subjectname,'_');
    storage(k).name = subjectname; %stores sub_id under .name
    storage(k).ts_baseline = cells; %stores data under .ts - 429 time by 502 ROIs
    group_ts_baseline_ts(:,:,k) = cells;
    if k == 1
        sumArray = cells;
    else
        sumArray = sumArray + cells;
    end
end

%% followup 
% Specify the folder where the files lives.
%
followup_Folder = 'V:\fMRI_analysis\functional_connectivity_analysis\follow-up';



% Get a list of all files in the folder with the desired file name pattern.

filePatternfollowup = fullfile(followup_Folder, '*IRBD*'); % Looking for the 400x400 matrices

theFilesfollowup = dir(filePatternfollowup); % Storing all filenames
% Load in timeseries matrix followup

% load in data

%storage(length(theFilesfollowup)) = struct('name',1,'ts_baseline',1,'corr_ts_baseline',1,'name2',1,'ts_followup',1,'corr_ts_followup',1);

for k = 1:16
    baseFileName = theFilesfollowup(k).name;
    fullFileName = fullfile(theFilesfollowup(k).folder, baseFileName);
    thisArray{k} = load(fullFileName); % Should have 54 structures containing a matrix
    Array = load(fullFileName);
    cells = Array.ts(:, :);
    subjectname = strrep(baseFileName,'timeseries.mat','');
    subjectname = extractBefore(subjectname,'_');
    storage(k).name2 = subjectname; %stores sub_id under .name
    storage(k).ts_followup = cells; %stores data under .ts - 429 time by 502 ROIs
 

end
%use if one deletes
%%
% Correlation matrix & Upper triangle
for ii = 1:length(storage)
     storage(ii).corr_ts_baseline = corr(storage(ii).ts_baseline(:,:));
     storage(ii).corr_tri_baseline = triu(storage(ii).corr_ts_baseline(:,:));
%      idx2 = logical(triu(ones(size(storage(ii).corr)),1));
%      v = storage(ii).corr(idx2)';
%      storage(ii).upper_corr = v;
%      variable2 = vertcat(storage.upper_corr)';
    group_matrix(:,:,ii) = storage(ii).corr_ts_baseline;
end


% mean correlation matrix
mean_corr = horzcat(storage.corr_ts_baseline);
mean_corr = reshape(mean_corr, 502, 502, []);
mean_corr = mean(mean_corr, 3);

%%
for ii = 1:length(storage)
     storage(ii).corr_ts_followup = corr(storage(ii).ts_followup(:,:));
     storage(ii).corr_tri_followup = triu(storage(ii).corr_ts_followup(:,:));
%      idx2 = logical(triu(ones(size(storage(ii).corr)),1));
%      v = storage(ii).corr(idx2)';
%      storage(ii).upper_corr = v;
%      variable2 = vertcat(storage.upper_corr)';
    group_matrix(:,:,ii) = storage(ii).corr_ts_followup;
end


% mean correlation matrix
mean_corr = horzcat(storage.corr_ts_followup);
mean_corr = reshape(mean_corr, 502, 502, []);
mean_corr = mean(mean_corr, 3);


%%
schaef_ID = xlsread('V:\fMRI_analysis\functional_connectivity_analysis\voltron_id.xlsx')
for ii = 1:length(storage)
   corr_matrix = storage(ii).corr_ts_followup;
    for xx = 1:7
        for yy = 1:7
            corr_matrix_mean(xx,yy) = mean(mean(corr_matrix(schaef_ID==xx,schaef_ID==yy, 1)));
        end
    end
    storage(ii).network_mean_followup = corr_matrix_mean;
end
%% BASELINE 
% % Specify the folder where the files lives.
% 
% baseline_Folder = 'V:\fMRI_analysis\functional_connectivity_analysis\results\baseline';
% 
% 
% 
% % Get a list of all files in the folder with the desired file name pattern.
% 
% filePatternbaseline = fullfile(baseline_Folder, '*irbd*'); % Looking for the 400x400 matrices
% 
% theFilesbaseline = dir(filePatternbaseline); % Storing all filenames
% % Load in FC matrix baseline
% 
% % load in data
% 
% %storage(length(theFilesbaseline)) = struct('name',1,'ts_baseline',1,'corr_ts_baseline',1,'name2',1,'ts_followup',1,'corr_ts_followup',1);
% 
% for k = 1:16
%     baseFileName = theFilesbaseline(k).name;
%     fullFileName = fullfile(theFilesbaseline(k).folder, baseFileName);
%     thisArray{k} = load(fullFileName); % Should have 54 structures containing a matrix
%     Array = load(fullFileName);
%     cells = Array.ts_corr(:, :);
%     subjectname = strrep(baseFileName,'ts_corr.mat','');
%     subjectname = extractBefore(subjectname,'_');
%     storage(k).name = subjectname; %stores sub_id under .name
%     storage(k).corr_ts_baseline = cells; %stores data under .ts - 429 time by 502 ROIs
%     group_matrix_base(:,:,k) = cells;
% 
% end
% %%
% % follow up
% baseline_Folder = 'V:\fMRI_analysis\functional_connectivity_analysis\results\follow_up';
% 
% 
% 
% % Get a list of all files in the folder with the desired file name pattern.
% 
% filePatternbaseline = fullfile(baseline_Folder, '*IRBD*'); % Looking for the 400x400 matrices
% 
% theFilesbaseline = dir(filePatternbaseline); % Storing all filenames
% 
% % Load in FC matrix followup
% 
% % load in data
% 
% %storage(length(theFilesbaseline)) = struct('name',1,'ts_baseline',1,'corr_ts_baseline',1,'name2',1,'ts_followup',1,'corr_ts_followup',1);
% 
% for k = 1:16
%     baseFileName = theFilesbaseline(k).name;
%     fullFileName = fullfile(theFilesbaseline(k).folder, baseFileName);
%     thisArray{k} = load(fullFileName); % Should have 54 structures containing a matrix
%     Array = load(fullFileName);
%     cells = Array.ts_corr(:, :);
%     subjectname = strrep(baseFileName,'ts_corr.mat','');
%     subjectname = extractBefore(subjectname,'_');
%     storage(k).name2 = subjectname; %stores sub_id under .name
%     storage(k).corr_ts_followup = cells; %stores data under .ts - 429 time by 502 ROIs
%     group_matrix_followup(:,:,k) = cells;
% 
% end
% 





%% NBM-RESTING STATE Network CONNECTIVITY


% extracting NBM FC measures with first 400 (cortical) values
for k=1:length(storage)
    try
    nbm2=storage(k).corr_ts_followup(1:400, 491:492);
    storage(k).nbm_follow_up = nbm2;
    nbm=storage(k).corr_ts_baseline(1:400,491:492);
    storage(k).nbm_baseline = nbm;
    end
end
%%
% the average connectivity value between networks and NBM (baseline)
schaef_ID = xlsread('V:\fMRI_analysis\functional_connectivity_analysis\voltron_id.xlsx')
for ii = 1:length(storage)
   baselinenbm = storage(ii).nbm_baseline;

    for xx = 1:8
        nbm_baseline_mean_network(xx) = mean(baselinenbm(schaef_ID==xx)); % mean for each network
    end

storage(ii).nbm_network_baseline_mean = nbm_baseline_mean_network; % 
end
%%

% the average connectivity value between networks and NBM (follow_up)
for ii = 1:length(storage)
   follow_upnbm = storage(ii).nbm_follow_up;

    for xx = 1:8
        nbm_follow_up_mean_network(xx) = mean(follow_upnbm(schaef_ID==xx)); % mean for each network
    end

storage(ii).nbm_network_follow_up_mean = nbm_follow_up_mean_network; % 
end

% convert data into ordered networks, seperated into the 7 networks (baseline)
%try to organize the data into a more managebale form
% Splitting the data and converting into a numeric array
% Your data in a cell array
data_cell = {storage.nbm_network_baseline_mean};

% Preallocate a matrix to store the values
num_values = numel(data_cell); % Assuming all cells have the same number of values
num_columns = 8; % Number of values you want to extract
baselinedata_matrix = NaN(num_values, num_columns); % Initialize with NaNs

% Extract the values from each cell and store in the matrix
for i = 1:num_values
    % Check if the cell has enough values
    if numel(data_cell{i}) >= num_columns
        baselinedata_matrix(i, :) = data_cell{i}(1:num_columns);
    else
        disp(['Cell ', num2str(i), ' does not have enough values.']);
    end
end

% Display the resulting matrix
disp(baselinedata_matrix);

% convert data into ordered networks, seperated into the 7 networks (followup)
%try to organize the data into a more managebale form
% Splitting the data and converting into a numeric array
% Your data in a cell array
data_cell = {storage.nbm_network_follow_up_mean};

% Preallocate a matrix to store the values
num_values = numel(data_cell); % Assuming all cells have the same number of values
num_columns = 8; % Number of values you want to extract
follow_up_matrix = NaN(num_values, num_columns); % Initialize with NaNs

% Extract the values from each cell and store in the matrix
for i = 1:num_values
    % Check if the cell has enough values
    if numel(data_cell{i}) >= num_columns
        follow_up_matrix(i, :) = data_cell{i}(1:num_columns);
    else
        disp(['Cell ', num2str(i), ' does not have enough values.']);
    end
end

% Display the resulting matrix
disp(follow_up_matrix);

% running t tests on all of the resting state networks
%% running t tests on all of the resting state networks
% Preallocate arrays to store t-test results
p_values = zeros(3, size(baselinedata_matrix, 2));
hypothesis = cell(3, size(baselinedata_matrix, 2));

% Perform paired t-tests for each network
for network = 1:size(baselinedata_matrix, 2)
    % Two-sided t-test
    [h, p] = ttest(baselinedata_matrix(:, network), follow_up_matrix(:, network));
    hypothesis{1, network} = h;
    p_values(1, network) = p;
    
    % One-sided t-test (greater than)
    [h, p] = ttest(baselinedata_matrix(:, network), follow_up_matrix(:, network), 'Tail', 'right');
    hypothesis{2, network} = h;
    p_values(2, network) = p;
    
    % One-sided t-test (less than)
    [h, p] = ttest(baselinedata_matrix(:, network), follow_up_matrix(:, network), 'Tail', 'left');
    hypothesis{3, network} = h;
    p_values(3, network) = p;
end

% Display results
for network = 1:size(baselinedata_matrix, 2)
    fprintf('Network %d:\n', network);
    fprintf('Two-sided test: p-value = %f, Significantly different? %d\n', p_values(1, network), hypothesis{1, network});
    fprintf('One-sided test (greater than): p-value = %f, Significantly different? %d\n', p_values(2, network), hypothesis{2, network});
    fprintf('One-sided test (less than): p-value = %f, Significantly different? %d\n\n', p_values(3, network), hypothesis{3, network});
end

%% chat gpt permutation testing
% Initialize variables to store permutation test results
nbm_network_perm = zeros(7, 2); % Assuming 7 networks and 2 columns for significance and p-value

% Concatenate network data for baseline and follow-up for all subjects
baseline_data = vertcat(storage.nbm_network_baseline_mean);
follow_up_data = vertcat(storage.nbm_network_follow_up_mean);

% Perform permutation testing for each network
for network_idx = 1:7 % Assuming you have 7 networks
    disp(['Network: ' num2str(network_idx)]);
    
    % Extract network data for baseline and follow-up
    network_baseline = baseline_data(:, network_idx);
    network_follow_up = follow_up_data(:, network_idx);
    
    % Perform permutation testing
    [sig, pval] = perm_1d_delta(network_baseline, network_follow_up, 100); % Adjust the number of permutations as needed
    
    % Store results
    nbm_network_perm(network_idx, 1) = sig;
    nbm_network_perm(network_idx, 2) = pval;
end











%% Extracting Thalamus and PPN functional connectivity values from the timeseries
for k=1:length(storage)
    try
    thalPPN2=storage(k).corr_ts_followup([403:404, 408:412, 423, 430:431, 435:439, 450], 489:490);
    storage(k).thalPPN_follow_up = thalPPN2;
    thalPPN=storage(k).corr_ts_baseline([403:404, 408:412, 423, 430:431, 435:439, 450],489:490);
    storage(k).thalPPN_baseline = thalPPN;
    end
end

%% the average connectivity value between Thal and PPN (baseline)
for ii = 1:length(storage)
   baseline_thalPPN = storage(ii).thalPPN_baseline;
   thalPPN_baseline_mean = mean(baseline_thalPPN(:)); %takes the mean connectivity of both right and left, can be split if wanted

    storage(ii).thalPPN_baseline_mean = thalPPN_baseline_mean; % 
end
%% the average connectivity value between Thal and PPN (followup)
for ii = 1:length(storage)
   follow_up_thalPPN = storage(ii).thalPPN_follow_up;
   thalPPN_follow_up_mean = mean(follow_up_thalPPN(:)); %takes the mean connectivity of both right and left, can be split if wanted

    storage(ii).thalPPN_follow_up_mean = thalPPN_follow_up_mean; % 
end

%% Simple t test on the thalamus PPN mean connectivity
% Example data
baselinePPNthalconnectivity = vertcat(storage.thalPPN_baseline_mean);
followupPPNthalconnectivity = vertcat(storage.thalPPN_follow_up_mean)

% Perform paired t-test
[h, p, ci, stats] = ttest(baselinePPNthalconnectivity, followupPPNthalconnectivity);

% Display results
disp(['Paired t-test:']);%%%%
disp(['Test decision: ', num2str(h)]);
disp(['p-value: ', num2str(p)]);
disp(['Confidence interval: ', num2str(ci)]);
disp(['T-statistic: ', num2str(stats.tstat)]);

%% testing whether it is greater or less than
% Assuming baselinePPNthalconnectivity and followupPPNthalconnectivity are two sets of paired data

% Perform paired t-test for greater than
[h_greater, p_greater, ci_greater, stats_greater] = ttest(baselinePPNthalconnectivity, followupPPNthalconnectivity, 'Tail', 'right'); % right-tailed (baselinePPNthalconnectivity > followupPPNthalconnectivity)

% Perform paired t-test for less than
[h_less, p_less, ci_less, stats_less] = ttest(baselinePPNthalconnectivity, followupPPNthalconnectivity, 'Tail', 'left'); % left-tailed (baselinePPNthalconnectivity < followupPPNthalconnectivity)

% Perform paired t-test for not equal
[h_not_equal, p_not_equal, ci_not_equal, stats_not_equal] = ttest(baselinePPNthalconnectivity, followupPPNthalconnectivity); % two-tailed (baselinePPNthalconnectivity ≠ followupPPNthalconnectivity)

% Display results for greater than
disp('Results for baselinePPNthalconnectivity > followupPPNthalconnectivity:');
disp(['Test decision: ', num2str(h_greater)]);
disp(['p-value: ', num2str(p_greater)]);
disp(['Confidence interval: [', num2str(ci_greater(1)), ', ', num2str(ci_greater(2)), ']']);
disp(['T-statistic: ', num2str(stats_greater.tstat)]);

% Display results for less than
disp('Results for baselinePPNthalconnectivity < followupPPNthalconnectivity:');
disp(['Test decision: ', num2str(h_less)]);
disp(['p-value: ', num2str(p_less)]);
disp(['Confidence interval: [', num2str(ci_less(1)), ', ', num2str(ci_less(2)), ']']);
disp(['T-statistic: ', num2str(stats_less.tstat)]);

% Display results for not equal
disp('Results for baselinePPNthalconnectivity ≠ followupPPNthalconnectivity:');
disp(['Test decision: ', num2str(h_not_equal)]);
disp(['p-value: ', num2str(p_not_equal)]);
disp(['Confidence interval: [', num2str(ci_not_equal(1)), ', ', num2str(ci_not_equal(2)), ']']);
disp(['T-statistic: ', num2str(stats_not_equal.tstat)]);














%% BG FUNCTIONAL CONNECTIVITY WITH RESTING STATE NETWORKS
% 
% 
% extracting BG FC measures with first 400 (cortical) values
for k=1:length(storage)
    try
    BG2=storage(k).corr_ts_followup(1:400,[413:420, 440:447]);
    storage(k).BG_follow_up = BG2;
    BG=storage(k).corr_ts_baseline(1:400,[413:420, 440:447]);
    storage(k).BG_baseline = BG;
    end
end

%% the average connectivity value between networks and BG (baseline)

for ii = 1:length(storage)
   baselineBG = storage(ii).BG_baseline;

    for xx = 1:8
        BG_baseline_mean_network(xx) = mean(baselineBG(schaef_ID==xx)); % mean for each network
    end

storage(ii).BG_network_baseline_mean = BG_baseline_mean_network; % 
end
%% the average connectivity value between networks and BG (follow_up)
for ii = 1:length(storage)
   follow_upBG = storage(ii).BG_follow_up;

    for xx = 1:8
        BG_follow_up_mean_network(xx) = mean(follow_upBG(schaef_ID==xx)); % mean for each network
    end

storage(ii).BG_network_follow_up_mean = BG_follow_up_mean_network; % 
end

%% convert data into ordered networks, seperated into the 7 networks (baseline)
%try to organize the data into a more managebale form
% Splitting the data and converting into a numeric array
% Your data in a cell array
data_cell = {storage.BG_network_baseline_mean};

% Preallocate a matrix to store the values
num_values = numel(data_cell); % Assuming all cells have the same number of values
num_columns = 8; % Number of values you want to extract
baselinedata_matrix = NaN(num_values, num_columns); % Initialize with NaNs

% Extract the values from each cell and store in the matrix
for i = 1:num_values
    % Check if the cell has enough values
    if numel(data_cell{i}) >= num_columns
        baselinedata_matrix(i, :) = data_cell{i}(1:num_columns);
    else
        disp(['Cell ', num2str(i), ' does not have enough values.']);
    end
end

% Display the resulting matrix
disp(baselinedata_matrix);

%% convert data into ordered networks, seperated into the 7 networks (followup)
%try to organize the data into a more managebale form
% Splitting the data and converting into a numeric array
% Your data in a cell array
data_cell = {storage.BG_network_follow_up_mean};

% Preallocate a matrix to store the values
num_values = numel(data_cell); % Assuming all cells have the same number of values
num_columns = 8; % Number of values you want to extract
follow_up_matrix = NaN(num_values, num_columns); % Initialize with NaNs

% Extract the values from each cell and store in the matrix
for i = 1:num_values
    % Check if the cell has enough values
    if numel(data_cell{i}) >= num_columns
        follow_up_matrix(i, :) = data_cell{i}(1:num_columns);
    else
        disp(['Cell ', num2str(i), ' does not have enough values.']);
    end
end

% Display the resulting matrix
disp(follow_up_matrix);

%% running t tests on all of the resting state networks
% Preallocate arrays to store t-test results
p_values = zeros(3, size(baselinedata_matrix, 2));
hypothesis = cell(3, size(baselinedata_matrix, 2));

% Perform paired t-tests for each network
for network = 1:size(baselinedata_matrix, 2)
    % Two-sided t-test
    [h, p] = ttest(baselinedata_matrix(:, network), follow_up_matrix(:, network));
    hypothesis{1, network} = h;
    p_values(1, network) = p;
    
    % One-sided t-test (greater than)
    [h, p] = ttest(baselinedata_matrix(:, network), follow_up_matrix(:, network), 'Tail', 'right');
    hypothesis{2, network} = h;
    p_values(2, network) = p;
    
    % One-sided t-test (less than)
    [h, p] = ttest(baselinedata_matrix(:, network), follow_up_matrix(:, network), 'Tail', 'left');
    hypothesis{3, network} = h;
    p_values(3, network) = p;
end

% Display results
for network = 1:size(baselinedata_matrix, 2)
    fprintf('Network %d:\n', network);
    fprintf('Two-sided test: p-value = %f, Significantly different? %d\n', p_values(1, network), hypothesis{1, network});
    fprintf('One-sided test (greater than): p-value = %f, Significantly different? %d\n', p_values(2, network), hypothesis{2, network});
    fprintf('One-sided test (less than): p-value = %f, Significantly different? %d\n\n', p_values(3, network), hypothesis{3, network});
end






%% THALAMIC FUNCTIONAL CONNECTIVITY WITH RESTING STATE NETWORKS
% 
% 
% 
% 
% extracting thal FC measures with first 400 (cortical) values
for k=1:length(storage)
    thal=storage(k).corr_ts_baseline(1:400,[403:404, 408:412, 423, 430:431, 435:439, 450]);
    storage(k).thal_baseline = thal;
    thal2=storage(k).corr_ts_followup(1:400,[403:404, 408:412, 423, 430:431, 435:439, 450]);
    storage(k).thal_follow_up = thal2;
end

%% the average connectivity value between networks and NBM (baseline)

for ii = 1:length(storage)
   baselinethal = storage(ii).thal_baseline;

    for xx = 1:8
        thal_baseline_mean_network(xx) = mean(baselinethal(schaef_ID==xx)); % mean for each network
    end

storage(ii).thal_network_baseline_mean = thal_baseline_mean_network; % 
end
%% the average connectivity value between networks and NBM (follow_up)
for ii = 1:length(storage)
   follow_upthal = storage(ii).thal_follow_up;

    for xx = 1:8
        thal_follow_up_mean_network(xx) = mean(follow_upthal(schaef_ID==xx)); % mean for each network
    end

storage(ii).thal_network_follow_up_mean = thal_follow_up_mean_network; % 
end

%% convert data into ordered networks, seperated into the 7 networks (baseline)
%try to organize the data into a more managebale form
% Splitting the data and converting into a numeric array
% Your data in a cell array
data_cell = {storage.thal_network_baseline_mean};

% Preallocate a matrix to store the values
num_values = numel(data_cell); % Assuming all cells have the same number of values
num_columns = 8; % Number of values you want to extract
baselinedata_matrix = NaN(num_values, num_columns); % Initialize with NaNs

% Extract the values from each cell and store in the matrix
for i = 1:num_values
    % Check if the cell has enough values
    if numel(data_cell{i}) >= num_columns
        baselinedata_matrix(i, :) = data_cell{i}(1:num_columns);
    else
        disp(['Cell ', num2str(i), ' does not have enough values.']);
    end
end

% Display the resulting matrix
disp(baselinedata_matrix);

%% convert data into ordered networks, seperated into the 7 networks (followup)
%try to organize the data into a more managebale form
% Splitting the data and converting into a numeric array
% Your data in a cell array
data_cell = {storage.thal_network_follow_up_mean};

% Preallocate a matrix to store the values
num_values = numel(data_cell); % Assuming all cells have the same number of values
num_columns = 8; % Number of values you want to extract
follow_up_matrix = NaN(num_values, num_columns); % Initialize with NaNs

% Extract the values from each cell and store in the matrix
for i = 1:num_values
    % Check if the cell has enough values
    if numel(data_cell{i}) >= num_columns
        follow_up_matrix(i, :) = data_cell{i}(1:num_columns);
    else
        disp(['Cell ', num2str(i), ' does not have enough values.']);
    end
end

% Display the resulting matrix
disp(follow_up_matrix);

%% running t tests on all of the resting state networks
% Preallocate arrays to store t-test results
p_values = zeros(3, size(baselinedata_matrix, 2));
hypothesis = cell(3, size(baselinedata_matrix, 2));

% Perform paired t-tests for each network
for network = 1:size(baselinedata_matrix, 2)
    % Two-sided t-test
    [h, p] = ttest(baselinedata_matrix(:, network), follow_up_matrix(:, network));
    hypothesis{1, network} = h;
    p_values(1, network) = p;
    
    % One-sided t-test (greater than)
    [h, p] = ttest(baselinedata_matrix(:, network), follow_up_matrix(:, network), 'Tail', 'right');
    hypothesis{2, network} = h;
    p_values(2, network) = p;
    
    % One-sided t-test (less than)
    [h, p] = ttest(baselinedata_matrix(:, network), follow_up_matrix(:, network), 'Tail', 'left');
    hypothesis{3, network} = h;
    p_values(3, network) = p;
end

% Display results
for network = 1:size(baselinedata_matrix, 2)
    fprintf('Network %d:\n', network);
    fprintf('Two-sided test: p-value = %f, Significantly different? %d\n', p_values(1, network), hypothesis{1, network});
    fprintf('One-sided test (greater than): p-value = %f, Significantly different? %d\n', p_values(2, network), hypothesis{2, network});
    fprintf('One-sided test (less than): p-value = %f, Significantly different? %d\n\n', p_values(3, network), hypothesis{3, network});
end















%% extracting LC FC measures with first 400 (cortical) values
for k=1:23
    try
lc2=storage(k).corr_ts_followup(1:400,[485:486]);
storage(k).lc_follow_up = lc2;
lc=storage(k).corr_ts_baseline(1:400, [485:486]);
storage(k).lc_baseline = lc;
    end
end

%% the average connectivity value between networks and BG (baseline)
%schaef_ID = xlsread('V:\fMRI_analysis\functional_connectivity_analysis\schaef_id.xlsx')
for ii = 1:length(storage)
   baselinelc = storage(ii).lc_baseline;

    for xx = 1:8
        lc_baseline_mean_network(xx) = mean(baselinelc(schaef_ID==xx)); % mean for each network
    end

storage(ii).lc_network_baseline_mean = lc_baseline_mean_network; % 
end
%% the average connectivity value between networks and BG (follow_up)
for ii = 1:length(storage)
   follow_uplc = storage(ii).lc_follow_up;

    for xx = 1:8
        lc_follow_up_mean_network(xx) = mean(follow_uplc(schaef_ID==xx)); % mean for each network
    end

storage(ii).lc_network_follow_up_mean = lc_follow_up_mean_network; % 
end

%% convert data into ordered networks, seperated into the 7 networks (baseline)
%try to organize the data into a more managebale form
% Splitting the data and converting into a numeric array
% Your data in a cell array
data_cell = {storage.lc_network_baseline_mean};

% Preallocate a matrix to store the values
num_values = numel(data_cell); % Assuming all cells have the same number of values
num_columns = 8; % Number of values you want to extract
baselinedata_matrix = NaN(num_values, num_columns); % Initialize with NaNs

% Extract the values from each cell and store in the matrix
for i = 1:num_values
    % Check if the cell has enough values
    if numel(data_cell{i}) >= num_columns
        baselinedata_matrix(i, :) = data_cell{i}(1:num_columns);
    else
        disp(['Cell ', num2str(i), ' does not have enough values.']);
    end
end

% Display the resulting matrix
disp(baselinedata_matrix);

%% convert data into ordered networks, seperated into the 7 networks (followup)
%try to organize the data into a more managebale form
% Splitting the data and converting into a numeric array
% Your data in a cell array
data_cell = {storage.lc_network_follow_up_mean};

% Preallocate a matrix to store the values
num_values = numel(data_cell); % Assuming all cells have the same number of values
num_columns = 8; % Number of values you want to extract
follow_up_matrix = NaN(num_values, num_columns); % Initialize with NaNs

% Extract the values from each cell and store in the matrix
for i = 1:num_values
    % Check if the cell has enough values
    if numel(data_cell{i}) >= num_columns
        follow_up_matrix(i, :) = data_cell{i}(1:num_columns);
    else
        disp(['Cell ', num2str(i), ' does not have enough values.']);
    end
end

% Display the resulting matrix
disp(follow_up_matrix);
%%
% Preallocate arrays to store t-test results
num_networks = size(baselinedata_matrix, 2);
p_values = zeros(3, num_networks);
significant_results = zeros(3, num_networks);

% Significance level (alpha)
alpha = 0.05;

% Perform paired t-tests for each network
% Perform paired t-tests for each network
for network = 1:num_networks
    % Two-sided t-test
    [h, p] = ttest(baselinedata_matrix(:, network), follow_up_matrix(:, network));
    p_values(1, network) = p;
    significant_results(1, network) = h;
    
    % One-sided t-test (greater than)
    [h, p] = ttest(baselinedata_matrix(:, network), follow_up_matrix(:, network), 'Tail', 'right');
    p_values(2, network) = p;
    significant_results(2, network) = h;
    
    % One-sided t-test (less than)
    [h, p] = ttest(baselinedata_matrix(:, network), follow_up_matrix(:, network), 'Tail', 'left');
    p_values(3, network) = p;
    significant_results(3, network) = h;
end


% Correct for multiple comparisons (Bonferroni correction)
adjusted_alpha = alpha / numel(p_values);

% Display results
for network = 1:num_networks
    fprintf('Network %d:\n', network);
    fprintf('Two-sided test: p-value = %f, Significantly different? %d\n', p_values(1, network), significant_results(1, network));
    fprintf('One-sided test (greater than): p-value = %f, Significantly different? %d\n', p_values(2, network), significant_results(2, network));
    fprintf('One-sided test (less than): p-value = %f, Significantly different? %d\n\n', p_values(3, network), significant_results(3, network));
end











%% extracting AMYGDALA FC measures with first 400 (cortical) values
for   k=1:length(storage)
    try
    cer2=storage(k).corr_ts_followup(1:400,[455:482]);
    storage(k).cer_follow_up = cer2;
    cer=storage(k).corr_ts_baseline(1:400,[455:482]);
    storage(k).cer_baseline = cer;
    end
end

%% the average connectivity value between networks and BG (baseline)
%schaef_ID = xlsread('V:\fMRI_analysis\functional_connectivity_analysis\schaef_id.xlsx')
for ii = 1:length(storage)
   baselinecer = storage(ii).cer_baseline;

    for xx = 1:8
        cer_baseline_mean_network(xx) = mean(baselinecer(schaef_ID==xx)); % mean for each network
    end

storage(ii).cer_network_baseline_mean = cer_baseline_mean_network; % 
end
%% the average connectivity value between networks and BG (follow_up)
for ii = 1:length(storage)
   follow_upcer = storage(ii).cer_follow_up;

    for xx = 1:8
        cer_follow_up_mean_network(xx) = mean(follow_upcer(schaef_ID==xx)); % mean for each network
    end

storage(ii).cer_network_follow_up_mean = cer_follow_up_mean_network; % 
end

%% convert data into ordered networks, seperated into the 7 networks (baseline)
%try to organize the data into a more managebale form
% Splitting the data and converting into a numeric array
% Your data in a cell array
data_cell = {storage.cer_network_baseline_mean};

% Preallocate a matrix to store the values
num_values = numel(data_cell); % Assuming all cells have the same number of values
num_columns = 8; % Number of values you want to extract
baselinedata_matrix = NaN(num_values, num_columns); % Initialize with NaNs

% Extract the values from each cell and store in the matrix
for i = 1:num_values
    % Check if the cell has enough values
    if numel(data_cell{i}) >= num_columns
        baselinedata_matrix(i, :) = data_cell{i}(1:num_columns);
    else
        disp(['Cell ', num2str(i), ' does not have enough values.']);
    end
end

% Display the resulting matrix
disp(baselinedata_matrix);

%% convert data into ordered networks, seperated into the 7 networks (followup)
%try to organize the data into a more managebale form
% Splitting the data and converting into a numeric array
% Your data in a cell array
data_cell = {storage.cer_network_follow_up_mean};

% Preallocate a matrix to store the values
num_values = numel(data_cell); % Assuming all cells have the same number of values
num_columns = 8; % Number of values you want to extract
follow_up_matrix = NaN(num_values, num_columns); % Initialize with NaNs

% Extract the values from each cell and store in the matrix
for i = 1:num_values
    % Check if the cell has enough values
    if numel(data_cell{i}) >= num_columns
        follow_up_matrix(i, :) = data_cell{i}(1:num_columns);
    else
        disp(['Cell ', num2str(i), ' does not have enough values.']);
    end
end

% Display the resulting matrix
disp(follow_up_matrix);

%% running t tests on all of the resting state networks
% Preallocate arrays to store t-test results
p_values = zeros(3, size(baselinedata_matrix, 2));
hypothesis = cell(3, size(baselinedata_matrix, 2));

% Perform paired t-tests for each network
for network = 1:size(baselinedata_matrix, 2)
    % Two-sided t-test
    [h, p] = ttest(baselinedata_matrix(:, network), follow_up_matrix(:, network));
    hypothesis{1, network} = h;
    p_values(1, network) = p;
    
    % One-sided t-test (greater than)
    [h, p] = ttest(baselinedata_matrix(:, network), follow_up_matrix(:, network), 'Tail', 'right');
    hypothesis{2, network} = h;
    p_values(2, network) = p;
    
    % One-sided t-test (less than)
    [h, p] = ttest(baselinedata_matrix(:, network), follow_up_matrix(:, network), 'Tail', 'left');
    hypothesis{3, network} = h;
    p_values(3, network) = p;
end

% Display results
for network = 1:size(baselinedata_matrix, 2)
    fprintf('Network %d:\n', network);
    fprintf('Two-sided test: p-value = %f, Significantly different? %d\n', p_values(1, network), hypothesis{1, network});
    fprintf('One-sided test (greater than): p-value = %f, Significantly different? %d\n', p_values(2, network), hypothesis{2, network});
    fprintf('One-sided test (less than): p-value = %f, Significantly different? %d\n\n', p_values(3, network), hypothesis{3, network});
end



%% *IN PROGRESS*
% Define the linear mixed-effects model formula
% Assuming you have baseline_matrix and followup_matrix as your baseline and follow-up data matrices, each with dimensions (subjects x networks)

% Combine baseline and follow-up data into a single matrix
% Combine baseline and follow-up data into a single matrix
all_data = [baselinedata_matrix; follow_up_matrix];

% Determine the number of subjects
num_subjects = size(baselinedata_matrix, 1);

% Create subject and time point variables
subject_ids = repelem((1:num_subjects)', size(baselinedata_matrix, 2));
time_points = repmat([1; 2], num_subjects, size(baselinedata_matrix, 2) / 2); % Assuming you have two time points (baseline=1, follow-up=2)

% Reshape the data into long format
data_long = reshape(all_data', [], 1); % Reshape data matrix into a column vector
subject_ids_long = repmat(subject_ids, 2, 1); % Repeat subject IDs for each time point
time_points_long = reshape(time_points, [], 1); % Reshape time points into a column vector
network_names = cellstr(repmat(strcat('Network', num2str((1:size(all_data, 2))')), num_subjects*2, 1)); % Generate network names

% Create a table from the long format data
data_table = table(data_long, subject_ids_long, time_points_long, network_names, ...
                   'VariableNames', {'Connectivity', 'SubjectID', 'TimePoint', 'Network'});

% Convert categorical variables
data_table.SubjectID = categorical(data_table.SubjectID);
data_table.TimePoint = categorical(data_table.TimePoint);
data_table.Network = categorical(data_table.Network);

% Fit the linear mixed-effects model
lme_model = fitlme(data_table, 'Connectivity ~ TimePoint + Network + (1|SubjectID)');

% Display model summary
disp(lme_model);

% Plot the model
plotResiduals(lme_model);





%% testing 