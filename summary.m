
%% for every subject/session - preparing data - only run this section once when running the script

clearvars;

% variables (set these yourself)

nodes = 502;
time = 140;
window = 10;
gamma = 1;

% Different steps are noted down in word document
% design matrix

% dsmtx = rand(time,1); % put in your specific task timing here (expected BOLD response)
% open variable 'dsmtx.mat' (stored here; E:\3. DT project\3. Analyses\Neuroimaging\Data\processing_shine)

%load("dsmtx_response.mat");
%load("dsmtx_instruction.mat")
%load("dsmtx_task.mat");

% Get the hrf using the spm_hrf function (account for the delay in the BOLD response)
%RT = 2;
%P = [6 16 1 1 6 0 32];
%T = 16;
%hrf = spm_hrf(RT,P,T); 

%hrf_fig = plot(0:RT:32,hrf);

% Convolve the design matrix with the hrf (can be done with the conv function in Matlab)

%convolved_task = conv(hrf,dsmtx_task); % A and B should be vectors
%convolved_response = conv(hrf,dsmtx_response);
%convolved_instruction = conv(hrf,dsmtx_instruction);

% combine convolved_task and convolved_rest to get the final dsmtx

%dsmtx_final = [convolved_task convolved_instruction convolved_response];

% manually delete the last 16 rows of dsmtx_final to get a design matrix that has 236 timeframes and 2 DV (236 x 2)


% load design matrix
%load('E:\3. DT project\3. Analyses\Neuroimaging\Data\processing_shine\dsmtx_final.mat')

%% calculate mtd and integration stuff and perform glm - run for each task-based fMRI scan separately
%% HC-ST
%% HC-ST-EARLY
% load timeseries derived from the denoising process
myFolder = 'C:\Users\lchur\functional_connectivity_analysis\follow-up';

filePattern = fullfile(myFolder, '*timeseries.mat');
theFiles = dir(filePattern);

storage_HC_ST_early(length(theFiles)) = struct('name',1,'ts',1);

for ii = 1:length(theFiles)
    baseFileName = theFiles(ii).name;
    fullFileName = fullfile(theFiles(ii).folder, baseFileName);
    ts = load(fullFileName);
    storage_HC_ST_early(ii).ts = ts.ts;
    storage_HC_ST_early(ii).name = extractBefore(baseFileName, '_timeseries.mat');
end

% load('sub-CON05_ses-01_run-1_timeseries.mat');


for ii = 1:length(storage_HC_ST_early)

   ts = storage_HC_ST_early(ii).ts;

% mtd

    mtd = coupling(ts,window);

 %   storage_HC_ST_early(ii).mtd = mtd;

% integration

    [~,~,part,~,~] = integration(mtd,gamma);

 %   storage_HC_ST_early(ii).part = part;


% general linear models (= first-level)

    storage_HC_ST_early(ii).beta_bold_early = zeros(nodes,4);
    storage_HC_ST_early(ii).beta_part_early = zeros(nodes,4);

   ts = ts';

        for rr = 1:nodes
            storage_HC_ST_early(ii).beta_bold_early(rr,:) = glmfit(dsmtx_final,ts(rr,:)');
            storage_HC_ST_early(ii).beta_part_early(rr,:) = glmfit(dsmtx_final,part(rr,:)');
        end 
  
      %  storage_HC_ST_early(ii).beta_bold_early = beta_bold_early(:,:);
       % storage_HC_ST_early(ii).beta_part_early = beta_part_early(:,:);
   
end
    
for ii = 1:length(storage_HC_ST_early)
    storage_HC_ST_early(ii).beta_bold_task_early = storage_HC_ST_early(ii).beta_bold_early(:,2);
    storage_HC_ST_early(ii).beta_part_task_early = storage_HC_ST_early(ii).beta_part_early(:,2);
end



  %  HCST_group_beta_part(:,:,ii) = beta_part;
  % HCST_beta_part_mean = mean(HCST_group_beta_part,3); % question: this takes the mean of everyone in the storage, which includes early, late and retention scans, while this should be separate for each timepoint? Should I load these separately first before calculating the mean?




%% HC-ST-late
% load timeseries derived from the denoising process
myFolder = 'E:\3. DT project\3. Analyses\Neuroimaging\Data\timeseries_to_use\HC-ST\late';

filePattern = fullfile(myFolder, '*timeseries.mat');
theFiles = dir(filePattern);

storage_HC_ST_late(length(theFiles)) = struct('name',1,'ts',1);

for ii = 1:length(theFiles)
    baseFileName = theFiles(ii).name;
    fullFileName = fullfile(theFiles(ii).folder, baseFileName);
    ts = load(fullFileName);
    storage_HC_ST_late(ii).ts = ts.ts;
    storage_HC_ST_late(ii).name = extractBefore(baseFileName, '_timeseries.mat');
end

% load('sub-CON05_ses-01_run-1_timeseries.mat');


for ii = 1:length(storage_HC_ST_late)

   ts = storage_HC_ST_late(ii).ts;

% mtd

    mtd = coupling(ts,window);

%    storage_HC_ST_late(ii).mtd = mtd;

% integration

    [~,~,part,~,~] = integration(mtd,gamma);

  %  storage_HC_ST_late(ii).part = part;


% general linear models (= first-level)

    storage_HC_ST_late(ii).beta_bold_late = zeros(nodes,4);
    storage_HC_ST_late(ii).beta_part_late = zeros(nodes,4);

   ts = ts';

        for rr = 1:nodes
            storage_HC_ST_late(ii).beta_bold_late(rr,:) = glmfit(dsmtx_final,ts(rr,:)');
            storage_HC_ST_late(ii).beta_part_late(rr,:) = glmfit(dsmtx_final,part(rr,:)');
        end 
  
      %  storage_HC_ST_early(ii).beta_bold_early = beta_bold_early(:,:);
       % storage_HC_ST_early(ii).beta_part_early = beta_part_early(:,:);
   
end
    
for ii = 1:length(storage_HC_ST_late)
    storage_HC_ST_late(ii).beta_bold_task_late = storage_HC_ST_late(ii).beta_bold_late(:,2);
    storage_HC_ST_late(ii).beta_part_task_late = storage_HC_ST_late(ii).beta_part_late(:,2);
end

%% HC-ST-retention
% load timeseries derived from the denoising process
myFolder = 'E:\3. DT project\3. Analyses\Neuroimaging\Data\timeseries_to_use\HC-ST\retention';

filePattern = fullfile(myFolder, '*timeseries.mat');
theFiles = dir(filePattern);

storage_HC_ST_ret(length(theFiles)) = struct('name',1,'ts',1);

for ii = 1:length(theFiles)
    baseFileName = theFiles(ii).name;
    fullFileName = fullfile(theFiles(ii).folder, baseFileName);
    ts = load(fullFileName);
    storage_HC_ST_ret(ii).ts = ts.ts;
    storage_HC_ST_ret(ii).name = extractBefore(baseFileName, '_timeseries.mat');
end

% load('sub-CON05_ses-01_run-1_timeseries.mat');


for ii = 1:length(storage_HC_ST_ret)

   ts = storage_HC_ST_ret(ii).ts;

% mtd

    mtd = coupling(ts,window);

  %  storage_HC_ST_ret(ii).mtd = mtd;

% integration

    [~,~,part,~,~] = integration(mtd,gamma);

  %  storage_HC_ST_ret(ii).part = part;


% general linear models (= first-level)

    storage_HC_ST_ret(ii).beta_bold_ret = zeros(nodes,4);
    storage_HC_ST_ret(ii).beta_part_ret = zeros(nodes,4);

   ts = ts';

        for rr = 1:nodes
            storage_HC_ST_ret(ii).beta_bold_ret(rr,:) = glmfit(dsmtx_final,ts(rr,:)');
            storage_HC_ST_ret(ii).beta_part_ret(rr,:) = glmfit(dsmtx_final,part(rr,:)');
        end 
  
      %  storage_HC_ST_early(ii).beta_bold_early = beta_bold_early(:,:);
       % storage_HC_ST_early(ii).beta_part_early = beta_part_early(:,:);
   
end
    
for ii = 1:length(storage_HC_ST_ret)
    storage_HC_ST_ret(ii).beta_bold_task_ret = storage_HC_ST_ret(ii).beta_bold_ret(:,2);
    storage_HC_ST_ret(ii).beta_part_task_ret = storage_HC_ST_ret(ii).beta_part_ret(:,2);
end


%% HC-DT
%% HC-DT-early
% load timeseries derived from the denoising process
myFolder = 'E:\3. DT project\3. Analyses\Neuroimaging\Data\timeseries_to_use\HC-DT\early';

filePattern = fullfile(myFolder, '*timeseries.mat');
theFiles = dir(filePattern);

storage_HC_DT_early(length(theFiles)) = struct('name',1,'ts',1);

for ii = 1:length(theFiles)
    baseFileName = theFiles(ii).name;
    fullFileName = fullfile(theFiles(ii).folder, baseFileName);
    ts = load(fullFileName);
    storage_HC_DT_early(ii).ts = ts.ts;
    storage_HC_DT_early(ii).name = extractBefore(baseFileName, '_timeseries.mat');
end

% load('sub-CON05_ses-01_run-1_timeseries.mat');


for ii = 1:length(storage_HC_DT_early)

   ts = storage_HC_DT_early(ii).ts;

% mtd

    mtd = coupling(ts,window);

%    storage_HC_DT_early(ii).mtd = mtd;

% integration

    [~,~,part,~,~] = integration(mtd,gamma);

  %  storage_HC_DT_early(ii).part = part;


% general linear models (= first-level)

    storage_HC_DT_early(ii).beta_bold_early = zeros(nodes,4);
    storage_HC_DT_early(ii).beta_part_early = zeros(nodes,4);

   ts = ts';

        for rr = 1:nodes
            storage_HC_DT_early(ii).beta_bold_early(rr,:) = glmfit(dsmtx_final,ts(rr,:)');
            storage_HC_DT_early(ii).beta_part_early(rr,:) = glmfit(dsmtx_final,part(rr,:)');
        end 
  
      %  storage_HC_ST_early(ii).beta_bold_early = beta_bold_early(:,:);
       % storage_HC_ST_early(ii).beta_part_early = beta_part_early(:,:);
   
end
    
for ii = 1:length(storage_HC_DT_early)
    storage_HC_DT_early(ii).beta_bold_task_early = storage_HC_DT_early(ii).beta_bold_early(:,2);
    storage_HC_DT_early(ii).beta_part_task_early = storage_HC_DT_early(ii).beta_part_early(:,2);
end


%% HC-DT-late
% load timeseries derived from the denoising process
myFolder = 'E:\3. DT project\3. Analyses\Neuroimaging\Data\timeseries_to_use\HC-DT\late';

filePattern = fullfile(myFolder, '*timeseries.mat');
theFiles = dir(filePattern);

storage_HC_DT_late(length(theFiles)) = struct('name',1,'ts',1);

for ii = 1:length(theFiles)
    baseFileName = theFiles(ii).name;
    fullFileName = fullfile(theFiles(ii).folder, baseFileName);
    ts = load(fullFileName);
    storage_HC_DT_late(ii).ts = ts.ts;
    storage_HC_DT_late(ii).name = extractBefore(baseFileName, '_timeseries.mat');
end

% load('sub-CON05_ses-01_run-1_timeseries.mat');


for ii = 1:length(storage_HC_DT_late)

   ts = storage_HC_DT_late(ii).ts;

% mtd

    mtd = coupling(ts,window);

  %  storage_HC_DT_late(ii).mtd = mtd;

% integration

    [~,~,part,~,~] = integration(mtd,gamma);

  %  storage_HC_DT_late(ii).part = part;


% general linear models (= first-level)

    storage_HC_DT_late(ii).beta_bold_late = zeros(nodes,4);
    storage_HC_DT_late(ii).beta_part_late = zeros(nodes,4);

   ts = ts';

        for rr = 1:nodes
            storage_HC_DT_late(ii).beta_bold_late(rr,:) = glmfit(dsmtx_final,ts(rr,:)');
            storage_HC_DT_late(ii).beta_part_late(rr,:) = glmfit(dsmtx_final,part(rr,:)');
        end 
  
      %  storage_HC_ST_early(ii).beta_bold_early = beta_bold_early(:,:);
       % storage_HC_ST_early(ii).beta_part_early = beta_part_early(:,:);
   
end
    
for ii = 1:length(storage_HC_DT_late)
    storage_HC_DT_late(ii).beta_bold_task_late = storage_HC_DT_late(ii).beta_bold_late(:,2);
    storage_HC_DT_late(ii).beta_part_task_late = storage_HC_DT_late(ii).beta_part_late(:,2);
end

%% HC-DT-retention
% load timeseries derived from the denoising process
myFolder = 'E:\3. DT project\3. Analyses\Neuroimaging\Data\timeseries_to_use\HC-DT\retention';

filePattern = fullfile(myFolder, '*timeseries.mat');
theFiles = dir(filePattern);

storage_HC_DT_ret(length(theFiles)) = struct('name',1,'ts',1);

for ii = 1:length(theFiles)
    baseFileName = theFiles(ii).name;
    fullFileName = fullfile(theFiles(ii).folder, baseFileName);
    ts = load(fullFileName);
    storage_HC_DT_ret(ii).ts = ts.ts;
    storage_HC_DT_ret(ii).name = extractBefore(baseFileName, '_timeseries.mat');
end

% load('sub-CON05_ses-01_run-1_timeseries.mat');


for ii = 1:length(storage_HC_DT_ret)

   ts = storage_HC_DT_ret(ii).ts;

% mtd

    mtd = coupling(ts,window);

 %   storage_HC_DT_ret(ii).mtd = mtd;

% integration

    [~,~,part,~,~] = integration(mtd,gamma);

  %  storage_HC_DT_ret(ii).part = part;


% general linear models (= first-level)

    storage_HC_DT_ret(ii).beta_bold_ret = zeros(nodes,4);
    storage_HC_DT_ret(ii).beta_part_ret = zeros(nodes,4);

   ts = ts';

        for rr = 1:nodes
            storage_HC_DT_ret(ii).beta_bold_ret(rr,:) = glmfit(dsmtx_final,ts(rr,:)');
            storage_HC_DT_ret(ii).beta_part_ret(rr,:) = glmfit(dsmtx_final,part(rr,:)');
        end 
  
      %  storage_HC_ST_early(ii).beta_bold_early = beta_bold_early(:,:);
       % storage_HC_ST_early(ii).beta_part_early = beta_part_early(:,:);
   
end
    
for ii = 1:length(storage_HC_DT_ret)
    storage_HC_DT_ret(ii).beta_bold_task_ret = storage_HC_DT_ret(ii).beta_bold_ret(:,2);
    storage_HC_DT_ret(ii).beta_part_task_ret = storage_HC_DT_ret(ii).beta_part_ret(:,2);
end



%% PD-ST
%% PD-ST-early
% load timeseries derived from the denoising process
myFolder = 'E:\3. DT project\3. Analyses\Neuroimaging\Data\timeseries_to_use\PD-ST\early';

filePattern = fullfile(myFolder, '*timeseries.mat');
theFiles = dir(filePattern);

storage_PD_ST_early(length(theFiles)) = struct('name',1,'ts',1);

for ii = 1:length(theFiles)
    baseFileName = theFiles(ii).name;
    fullFileName = fullfile(theFiles(ii).folder, baseFileName);
    ts = load(fullFileName);
    storage_PD_ST_early(ii).ts = ts.ts;
    storage_PD_ST_early(ii).name = extractBefore(baseFileName, '_timeseries.mat');
end

% load('sub-CON05_ses-01_run-1_timeseries.mat');


for ii = 1:length(storage_PD_ST_early)

   ts = storage_PD_ST_early(ii).ts;

% mtd

    mtd = coupling(ts,window);

  %  storage_PD_ST_early(ii).mtd = mtd;

% integration

    [~,~,part,~,~] = integration(mtd,gamma);

   % storage_PD_ST_early(ii).part = part;


% general linear models (= first-level)

    storage_PD_ST_early(ii).beta_bold_early = zeros(nodes,4);
    storage_PD_ST_early(ii).beta_part_early = zeros(nodes,4);

   ts = ts';

        for rr = 1:nodes
            storage_PD_ST_early(ii).beta_bold_early(rr,:) = glmfit(dsmtx_final,ts(rr,:)');
            storage_PD_ST_early(ii).beta_part_early(rr,:) = glmfit(dsmtx_final,part(rr,:)');
        end 
  
      %  storage_HC_ST_early(ii).beta_bold_early = beta_bold_early(:,:);
       % storage_HC_ST_early(ii).beta_part_early = beta_part_early(:,:);
   
end
    
for ii = 1:length(storage_PD_ST_early)
    storage_PD_ST_early(ii).beta_bold_task_early = storage_PD_ST_early(ii).beta_bold_early(:,2);
    storage_PD_ST_early(ii).beta_part_task_early = storage_PD_ST_early(ii).beta_part_early(:,2);
end


%% PD-ST-late
% load timeseries derived from the denoising process
myFolder = 'E:\3. DT project\3. Analyses\Neuroimaging\Data\timeseries_to_use\PD-ST\late';

filePattern = fullfile(myFolder, '*timeseries.mat');
theFiles = dir(filePattern);

storage_PD_ST_late(length(theFiles)) = struct('name',1,'ts',1);

for ii = 1:length(theFiles)
    baseFileName = theFiles(ii).name;
    fullFileName = fullfile(theFiles(ii).folder, baseFileName);
    ts = load(fullFileName);
    storage_PD_ST_late(ii).ts = ts.ts;
    storage_PD_ST_late(ii).name = extractBefore(baseFileName, '_timeseries.mat');
end

% load('sub-CON05_ses-01_run-1_timeseries.mat');


for ii = 1:length(storage_PD_ST_late)

   ts = storage_PD_ST_late(ii).ts;

% mtd

    mtd = coupling(ts,window);

%    storage_PD_ST_late(ii).mtd = mtd;

% integration

    [~,~,part,~,~] = integration(mtd,gamma);

%    storage_PD_ST_late(ii).part = part;


% general linear models (= first-level)

    storage_PD_ST_late(ii).beta_bold_late = zeros(nodes,4);
    storage_PD_ST_late(ii).beta_part_late = zeros(nodes,4);

   ts = ts';

        for rr = 1:nodes
            storage_PD_ST_late(ii).beta_bold_late(rr,:) = glmfit(dsmtx_final,ts(rr,:)');
            storage_PD_ST_late(ii).beta_part_late(rr,:) = glmfit(dsmtx_final,part(rr,:)');
        end 
  
      %  storage_HC_ST_early(ii).beta_bold_early = beta_bold_early(:,:);
       % storage_HC_ST_early(ii).beta_part_early = beta_part_early(:,:);
   
end
    
for ii = 1:length(storage_PD_ST_late)
    storage_PD_ST_late(ii).beta_bold_task_late = storage_PD_ST_late(ii).beta_bold_late(:,2);
    storage_PD_ST_late(ii).beta_part_task_late = storage_PD_ST_late(ii).beta_part_late(:,2);
end

%% PD-ST-retention
% load timeseries derived from the denoising process
myFolder = 'E:\3. DT project\3. Analyses\Neuroimaging\Data\timeseries_to_use\PD-ST\retention';

filePattern = fullfile(myFolder, '*timeseries.mat');
theFiles = dir(filePattern);

storage_PD_ST_ret(length(theFiles)) = struct('name',1,'ts',1);

for ii = 1:length(theFiles)
    baseFileName = theFiles(ii).name;
    fullFileName = fullfile(theFiles(ii).folder, baseFileName);
    ts = load(fullFileName);
    storage_PD_ST_ret(ii).ts = ts.ts;
    storage_PD_ST_ret(ii).name = extractBefore(baseFileName, '_timeseries.mat');
end


% load('sub-CON05_ses-01_run-1_timeseries.mat');


for ii = 1:length(storage_PD_ST_ret)

   ts = storage_PD_ST_ret(ii).ts;

% mtd

    mtd = coupling(ts,window);

  %  storage_PD_ST_ret(ii).mtd = mtd;

% integration

    [~,~,part,~,~] = integration(mtd,gamma);

  %  storage_PD_ST_ret(ii).part = part;


% general linear models (= first-level)

    storage_PD_ST_ret(ii).beta_bold_ret = zeros(nodes,4);
    storage_PD_ST_ret(ii).beta_part_ret = zeros(nodes,4);

   ts = ts';

        for rr = 1:nodes
            storage_PD_ST_ret(ii).beta_bold_ret(rr,:) = glmfit(dsmtx_final,ts(rr,:)');
            storage_PD_ST_ret(ii).beta_part_ret(rr,:) = glmfit(dsmtx_final,part(rr,:)');
        end 
  
      %  storage_HC_ST_early(ii).beta_bold_early = beta_bold_early(:,:);
       % storage_HC_ST_early(ii).beta_part_early = beta_part_early(:,:);
   
end
    
for ii = 1:length(storage_PD_ST_ret)
    storage_PD_ST_ret(ii).beta_bold_task_ret = storage_PD_ST_ret(ii).beta_bold_ret(:,2);
    storage_PD_ST_ret(ii).beta_part_task_ret = storage_PD_ST_ret(ii).beta_part_ret(:,2);
end


%% PD-DT
% PD-DT-early
% load timeseries derived from the denoising process
myFolder = 'E:\3. DT project\3. Analyses\Neuroimaging\Data\timeseries_to_use\PD-DT\early';

filePattern = fullfile(myFolder, '*timeseries.mat');
theFiles = dir(filePattern);

storage_PD_DT_early(length(theFiles)) = struct('name',1,'ts',1);

for ii = 1:length(theFiles)
    baseFileName = theFiles(ii).name;
    fullFileName = fullfile(theFiles(ii).folder, baseFileName);
    ts = load(fullFileName);
    storage_PD_DT_early(ii).ts = ts.ts;
    storage_PD_DT_early(ii).name = extractBefore(baseFileName, '_timeseries.mat');
end

% load('sub-CON05_ses-01_run-1_timeseries.mat');


for ii = 1:length(storage_PD_DT_early)

   ts = storage_PD_DT_early(ii).ts;

% mtd

    mtd = coupling(ts,window);

 %   storage_PD_DT_early(ii).mtd = mtd;

% integration

    [~,~,part,~,~] = integration(mtd,gamma);

  %  storage_PD_DT_early(ii).part = part;


% general linear models (= first-level)

    storage_PD_DT_early(ii).beta_bold_early = zeros(nodes,4);
    storage_PD_DT_early(ii).beta_part_early = zeros(nodes,4);

   ts = ts';

        for rr = 1:nodes
            storage_PD_DT_early(ii).beta_bold_early(rr,:) = glmfit(dsmtx_final,ts(rr,:)');
            storage_PD_DT_early(ii).beta_part_early(rr,:) = glmfit(dsmtx_final,part(rr,:)');
        end 
  
      %  storage_HC_ST_early(ii).beta_bold_early = beta_bold_early(:,:);
       % storage_HC_ST_early(ii).beta_part_early = beta_part_early(:,:);
   
end
    
for ii = 1:length(storage_PD_DT_early)
    storage_PD_DT_early(ii).beta_bold_task_early = storage_PD_DT_early(ii).beta_bold_early(:,2);
    storage_PD_DT_early(ii).beta_part_task_early = storage_PD_DT_early(ii).beta_part_early(:,2);
end

%% PD-DT-late
% load timeseries derived from the denoising process
myFolder = 'E:\3. DT project\3. Analyses\Neuroimaging\Data\timeseries_to_use\PD-DT\late';

filePattern = fullfile(myFolder, '*timeseries.mat');
theFiles = dir(filePattern);

storage_PD_DT_late(length(theFiles)) = struct('name',1,'ts',1);

for ii = 1:length(theFiles)
    baseFileName = theFiles(ii).name;
    fullFileName = fullfile(theFiles(ii).folder, baseFileName);
    ts = load(fullFileName);
    storage_PD_DT_late(ii).ts = ts.ts;
    storage_PD_DT_late(ii).name = extractBefore(baseFileName, '_timeseries.mat');
end

% load('sub-CON05_ses-01_run-1_timeseries.mat');


for ii = 1:length(storage_PD_DT_late)

   ts = storage_PD_DT_late(ii).ts;

% mtd

    mtd = coupling(ts,window);

 %   storage_PD_DT_late(ii).mtd = mtd;

% integration

    [~,~,part,~,~] = integration(mtd,gamma);

%    storage_PD_DT_late(ii).part = part;


% general linear models (= first-level)

    storage_PD_DT_late(ii).beta_bold_late = zeros(nodes,4);
    storage_PD_DT_late(ii).beta_part_late = zeros(nodes,4);

   ts = ts';

        for rr = 1:nodes
            storage_PD_DT_late(ii).beta_bold_late(rr,:) = glmfit(dsmtx_final,ts(rr,:)');
            storage_PD_DT_late(ii).beta_part_late(rr,:) = glmfit(dsmtx_final,part(rr,:)');
        end 
  
      %  storage_HC_ST_early(ii).beta_bold_early = beta_bold_early(:,:);
       % storage_HC_ST_early(ii).beta_part_early = beta_part_early(:,:);
   
end
    
for ii = 1:length(storage_PD_DT_late)
    storage_PD_DT_late(ii).beta_bold_task_late = storage_PD_DT_late(ii).beta_bold_late(:,2);
    storage_PD_DT_late(ii).beta_part_task_late = storage_PD_DT_late(ii).beta_part_late(:,2);
end

%% PD-DT-retention
% load timeseries derived from the denoising process
myFolder = 'E:\3. DT project\3. Analyses\Neuroimaging\Data\timeseries_to_use\PD-DT\retention';

filePattern = fullfile(myFolder, '*timeseries.mat');
theFiles = dir(filePattern);

storage_PD_DT_ret(length(theFiles)) = struct('name',1,'ts',1);

for ii = 1:length(theFiles)
    baseFileName = theFiles(ii).name;
    fullFileName = fullfile(theFiles(ii).folder, baseFileName);
    ts = load(fullFileName);
    storage_PD_DT_ret(ii).ts = ts.ts;
    storage_PD_DT_ret(ii).name = extractBefore(baseFileName, '_timeseries.mat');
end

% load('sub-CON05_ses-01_run-1_timeseries.mat');


for ii = 1:length(storage_PD_DT_ret)

   ts = storage_PD_DT_ret(ii).ts;

% mtd

    mtd = coupling(ts,window);

  %  storage_PD_DT_ret(ii).mtd = mtd;

% integration

    [~,~,part,~,~] = integration(mtd,gamma);

 %   storage_PD_DT_ret(ii).part = part;


% general linear models (= first-level)

    storage_PD_DT_ret(ii).beta_bold_ret = zeros(nodes,4);
    storage_PD_DT_ret(ii).beta_part_ret = zeros(nodes,4);

   ts = ts';

        for rr = 1:nodes
            storage_PD_DT_ret(ii).beta_bold_ret(rr,:) = glmfit(dsmtx_final,ts(rr,:)');
            storage_PD_DT_ret(ii).beta_part_ret(rr,:) = glmfit(dsmtx_final,part(rr,:)');
        end 
  
      %  storage_HC_ST_early(ii).beta_bold_early = beta_bold_early(:,:);
       % storage_HC_ST_early(ii).beta_part_early = beta_part_early(:,:);
   
end
    
for ii = 1:length(storage_PD_DT_ret)
    storage_PD_DT_ret(ii).beta_bold_task_ret = storage_PD_DT_ret(ii).beta_bold_ret(:,2);
    storage_PD_DT_ret(ii).beta_part_task_ret = storage_PD_DT_ret(ii).beta_part_ret(:,2);
end


%%
% permutation tests for each ROI (only possible with 2 groups? no possible with eg. group (PD vs HC) and training (ST vs DT))
% something is not right I think, as I would expect it to give a pval and
% sig for each ROI, right? I get just one value...

%[HCSTDT_sig,HCSTDT_pval] = perm_1d_delta(HCST_beta_part_mean,HCDT_beta_part_mean,100); % 100 = number of permutations
%[PDSTDT_sig,PDSTDT_pval] = perm_1d_delta(PDST_beta_part_mean,PDDT_beta_part_mean,100); % 100 = number of permutations
%[STPDHC_sig,STPDHC_pval] = perm_1d_delta(HCST_beta_part_mean,PDST_beta_part_mean,100); % 100 = number of permutations
%[DTPDHC_sig,DTPDHC_pval] = perm_1d_delta(HCDT_beta_part_mean,PDDT_beta_part_mean,100); % 100 = number of permutations


% the previous section should have outputted a different storage for each group (i.e. HC-ST, HC-DT, PD-ST and PD-DT), all timepoints are in here
% and I'm assuming this should still be separated, as it averages over the three.



%% Networks

%% Add in schaef networks

load('E:\3. DT project\3. Analyses\Neuroimaging\Participation\networks\schaef_id.mat');

 

%% Mean network connectivity per subject - participation coefficient & bold

%% HC-ST
% HC-ST-early 

for ii = 1:length(storage_HC_ST_early)

  storage_HC_ST_early(ii).beta_part_task_cortex_early = storage_HC_ST_early(ii).beta_part_task_early(1:400,:);
   storage_HC_ST_early(ii).beta_bold_task_cortex_early = storage_HC_ST_early(ii).beta_bold_task_early(1:400,:);

end

for ii = 1:length(storage_HC_ST_early)
   beta_part_task_test_early = storage_HC_ST_early(ii).beta_part_task_cortex_early;

    for xx = 1:7
        beta_part_mean_network_early(xx) = mean(beta_part_task_test_early(schaef_id==xx)); % mean for each network
    end

storage_HC_ST_early(ii).network_part_mean_early = beta_part_mean_network_early; % 
end

for ii = 1:length(storage_HC_ST_early)
   beta_bold_task_test_early = storage_HC_ST_early(ii).beta_bold_task_cortex_early;

    for xx = 1:7
        beta_bold_mean_network_early(xx) = mean(beta_bold_task_test_early(schaef_id==xx)); % mean for each network
    end

storage_HC_ST_early(ii).network_bold_mean_early = beta_bold_mean_network_early; % 
end


% HC-ST-late 

for ii = 1:length(storage_HC_ST_late)

   storage_HC_ST_late(ii).beta_part_task_cortex_late = storage_HC_ST_late(ii).beta_part_task_late(1:400,:);
   storage_HC_ST_late(ii).beta_bold_task_cortex_late = storage_HC_ST_late(ii).beta_bold_task_late(1:400,:);

end

for ii = 1:length(storage_HC_ST_late)
   beta_part_task_test_late = storage_HC_ST_late(ii).beta_part_task_cortex_late;

    for xx = 1:7
        beta_part_mean_network_late(xx) = mean(beta_part_task_test_late(schaef_id==xx)); % mean for each network
    end

storage_HC_ST_late(ii).network_part_mean_late = beta_part_mean_network_late; % 
end

for ii = 1:length(storage_HC_ST_late)
   beta_bold_task_test_late = storage_HC_ST_late(ii).beta_bold_task_cortex_late;

    for xx = 1:7
        beta_bold_mean_network_late(xx) = mean(beta_bold_task_test_late(schaef_id==xx)); % mean for each network
    end

storage_HC_ST_late(ii).network_bold_mean_late = beta_bold_mean_network_late; % 
end

% HC-ST-retention

for ii = 1:length(storage_HC_ST_ret)

   storage_HC_ST_ret(ii).beta_part_task_cortex_ret = storage_HC_ST_ret(ii).beta_part_task_ret(1:400,:);
   storage_HC_ST_ret(ii).beta_bold_task_cortex_ret = storage_HC_ST_ret(ii).beta_bold_task_ret(1:400,:);

end

for ii = 1:length(storage_HC_ST_ret)
   beta_part_task_test_ret = storage_HC_ST_ret(ii).beta_part_task_cortex_ret;

    for xx = 1:7
        beta_part_mean_network_ret(xx) = mean(beta_part_task_test_ret(schaef_id==xx)); % mean for each network
    end

storage_HC_ST_ret(ii).network_part_mean_ret = beta_part_mean_network_ret; % 
end

for ii = 1:length(storage_HC_ST_ret)
   beta_bold_task_test_ret = storage_HC_ST_ret(ii).beta_bold_task_cortex_ret;

    for xx = 1:7
        beta_bold_mean_network_ret(xx) = mean(beta_bold_task_test_ret(schaef_id==xx)); % mean for each network
    end

storage_HC_ST_ret(ii).network_bold_mean_ret = beta_bold_mean_network_ret; % 
end


%% HC-DT
% HC-DT-early 

for ii = 1:length(storage_HC_DT_early)

   storage_HC_DT_early(ii).beta_part_task_cortex_early = storage_HC_DT_early(ii).beta_part_task_early(1:400,:);
   storage_HC_DT_early(ii).beta_bold_task_cortex_early = storage_HC_DT_early(ii).beta_bold_task_early(1:400,:);

end

for ii = 1:length(storage_HC_DT_early)
   beta_part_task_test_early = storage_HC_DT_early(ii).beta_part_task_cortex_early;

    for xx = 1:7
        beta_part_mean_network_early(xx) = mean(beta_part_task_test_early(schaef_id==xx)); % mean for each network
    end

storage_HC_DT_early(ii).network_part_mean_early = beta_part_mean_network_early; % 
end

for ii = 1:length(storage_HC_DT_early)
   beta_bold_task_test_early = storage_HC_DT_early(ii).beta_bold_task_cortex_early;

    for xx = 1:7
        beta_bold_mean_network_early(xx) = mean(beta_bold_task_test_early(schaef_id==xx)); % mean for each network
    end

storage_HC_DT_early(ii).network_bold_mean_early = beta_bold_mean_network_early; % 
end


% HC-DT-late 

for ii = 1:length(storage_HC_DT_late)

   storage_HC_DT_late(ii).beta_part_task_cortex_late = storage_HC_DT_late(ii).beta_part_task_late(1:400,:);
   storage_HC_DT_late(ii).beta_bold_task_cortex_late = storage_HC_DT_late(ii).beta_bold_task_late(1:400,:);

end

for ii = 1:length(storage_HC_DT_late)
   beta_part_task_test_late = storage_HC_DT_late(ii).beta_part_task_cortex_late;

    for xx = 1:7
        beta_part_mean_network_late(xx) = mean(beta_part_task_test_late(schaef_id==xx)); % mean for each network
    end

storage_HC_DT_late(ii).network_part_mean_late = beta_part_mean_network_late; % 
end

for ii = 1:length(storage_HC_DT_late)
   beta_bold_task_test_early = storage_HC_DT_late(ii).beta_bold_task_cortex_late;

    for xx = 1:7
        beta_bold_mean_network_late(xx) = mean(beta_bold_task_test_late(schaef_id==xx)); % mean for each network
    end

storage_HC_DT_late(ii).network_bold_mean_late = beta_bold_mean_network_late; % 
end

% HC-DT-retention

for ii = 1:length(storage_HC_DT_ret)

   storage_HC_DT_ret(ii).beta_part_task_cortex_ret = storage_HC_DT_ret(ii).beta_part_task_ret(1:400,:);
   storage_HC_DT_ret(ii).beta_bold_task_cortex_ret = storage_HC_DT_ret(ii).beta_bold_task_ret(1:400,:);

end

for ii = 1:length(storage_HC_DT_ret)
   beta_part_task_test_ret = storage_HC_DT_ret(ii).beta_part_task_cortex_ret;

    for xx = 1:7
        beta_part_mean_network_ret(xx) = mean(beta_part_task_test_ret(schaef_id==xx)); % mean for each network
    end

storage_HC_DT_ret(ii).network_part_mean_ret = beta_part_mean_network_ret; % 
end

for ii = 1:length(storage_HC_DT_ret)
   beta_bold_task_test_ret = storage_HC_DT_ret(ii).beta_bold_task_cortex_ret;

    for xx = 1:7
        beta_bold_mean_network_ret(xx) = mean(beta_bold_task_test_ret(schaef_id==xx)); % mean for each network
    end

storage_HC_DT_ret(ii).network_bold_mean_ret = beta_bold_mean_network_ret; % 
end


%% PD-ST
% PD-ST-early 

for ii = 1:length(storage_PD_ST_early)

   storage_PD_ST_early(ii).beta_part_task_cortex_early = storage_PD_ST_early(ii).beta_part_task_early(1:400,:);
   storage_PD_ST_early(ii).beta_bold_task_cortex_early = storage_PD_ST_early(ii).beta_bold_task_early(1:400,:);

end

for ii = 1:length(storage_PD_ST_early)
   beta_part_task_test_early = storage_PD_ST_early(ii).beta_part_task_cortex_early;

    for xx = 1:7
        beta_part_mean_network_early(xx) = mean(beta_part_task_test_early(schaef_id==xx)); % mean for each network
    end

storage_PD_ST_early(ii).network_part_mean_early = beta_part_mean_network_early; % 
end

for ii = 1:length(storage_PD_ST_early)
   beta_bold_task_test_early = storage_PD_ST_early(ii).beta_bold_task_cortex_early;

    for xx = 1:7
        beta_bold_mean_network_early(xx) = mean(beta_bold_task_test_early(schaef_id==xx)); % mean for each network
    end

storage_PD_ST_early(ii).network_bold_mean_early = beta_bold_mean_network_early; % 
end


% PD-ST-late 

for ii = 1:length(storage_PD_ST_late)

   storage_PD_ST_late(ii).beta_part_task_cortex_late = storage_PD_ST_late(ii).beta_part_task_late(1:400,:);
   storage_PD_ST_late(ii).beta_bold_task_cortex_late = storage_PD_ST_late(ii).beta_bold_task_late(1:400,:);

end

for ii = 1:length(storage_PD_ST_late)
   beta_part_task_test_late = storage_PD_ST_late(ii).beta_part_task_cortex_late;

    for xx = 1:7
        beta_part_mean_network_late(xx) = mean(beta_part_task_test_late(schaef_id==xx)); % mean for each network
    end

storage_PD_ST_late(ii).network_part_mean_late = beta_part_mean_network_late; % 
end

for ii = 1:length(storage_PD_ST_late)
   beta_bold_task_test_early = storage_PD_ST_late(ii).beta_bold_task_cortex_late;

    for xx = 1:7
        beta_bold_mean_network_late(xx) = mean(beta_bold_task_test_late(schaef_id==xx)); % mean for each network
    end

storage_PD_ST_late(ii).network_bold_mean_late = beta_bold_mean_network_late; % 
end

% PD-ST-retention

for ii = 1:length(storage_PD_ST_ret)

   storage_PD_ST_ret(ii).beta_part_task_cortex_ret = storage_PD_ST_ret(ii).beta_part_task_ret(1:400,:);
   storage_PD_ST_ret(ii).beta_bold_task_cortex_ret = storage_PD_ST_ret(ii).beta_bold_task_ret(1:400,:);

end

for ii = 1:length(storage_PD_ST_ret)
   beta_part_task_test_ret = storage_PD_ST_ret(ii).beta_part_task_cortex_ret;

    for xx = 1:7
        beta_part_mean_network_ret(xx) = mean(beta_part_task_test_ret(schaef_id==xx)); % mean for each network
    end

storage_PD_ST_ret(ii).network_part_mean_ret = beta_part_mean_network_ret; % 
end

for ii = 1:length(storage_PD_ST_ret)
   beta_bold_task_test_ret = storage_PD_ST_ret(ii).beta_bold_task_cortex_ret;

    for xx = 1:7
        beta_bold_mean_network_ret(xx) = mean(beta_bold_task_test_ret(schaef_id==xx)); % mean for each network
    end

storage_PD_ST_ret(ii).network_bold_mean_ret = beta_bold_mean_network_ret; % 
end


%% PD-DT
% PD-DT-early 

for ii = 1:length(storage_PD_DT_early)

   storage_PD_DT_early(ii).beta_part_task_cortex_early = storage_PD_DT_early(ii).beta_part_task_early(1:400,:);
   storage_PD_DT_early(ii).beta_bold_task_cortex_early = storage_PD_DT_early(ii).beta_bold_task_early(1:400,:);

end

for ii = 1:length(storage_PD_DT_early)
   beta_part_task_test_early = storage_PD_DT_early(ii).beta_part_task_cortex_early;

    for xx = 1:7
        beta_part_mean_network_early(xx) = mean(beta_part_task_test_early(schaef_id==xx)); % mean for each network
    end

storage_PD_DT_early(ii).network_part_mean_early = beta_part_mean_network_early; % 
end

for ii = 1:length(storage_PD_DT_early)
   beta_bold_task_test_early = storage_PD_DT_early(ii).beta_bold_task_cortex_early;

    for xx = 1:7
        beta_bold_mean_network_early(xx) = mean(beta_bold_task_test_early(schaef_id==xx)); % mean for each network
    end

storage_PD_DT_early(ii).network_bold_mean_early = beta_bold_mean_network_early; % 
end


% PD-DT-late 

for ii = 1:length(storage_PD_DT_late)

   storage_PD_DT_late(ii).beta_part_task_cortex_late = storage_PD_DT_late(ii).beta_part_task_late(1:400,:);
   storage_PD_DT_late(ii).beta_bold_task_cortex_late = storage_PD_DT_late(ii).beta_bold_task_late(1:400,:);

end

for ii = 1:length(storage_PD_DT_late)
   beta_part_task_test_late = storage_PD_DT_late(ii).beta_part_task_cortex_late;

    for xx = 1:7
        beta_part_mean_network_late(xx) = mean(beta_part_task_test_late(schaef_id==xx)); % mean for each network
    end

storage_PD_DT_late(ii).network_part_mean_late = beta_part_mean_network_late; % 
end

for ii = 1:length(storage_PD_DT_late)
   beta_bold_task_test_early = storage_PD_DT_late(ii).beta_bold_task_cortex_late;

    for xx = 1:7
        beta_bold_mean_network_late(xx) = mean(beta_bold_task_test_late(schaef_id==xx)); % mean for each network
    end

storage_PD_DT_late(ii).network_bold_mean_late = beta_bold_mean_network_late; % 
end

% PD-DT-retention

for ii = 1:length(storage_PD_DT_ret)

   storage_PD_DT_ret(ii).beta_part_task_cortex_ret = storage_PD_DT_ret(ii).beta_part_task_ret(1:400,:);
   storage_PD_DT_ret(ii).beta_bold_task_cortex_ret = storage_PD_DT_ret(ii).beta_bold_task_ret(1:400,:);

end

for ii = 1:length(storage_PD_DT_ret)
   beta_part_task_test_ret = storage_PD_DT_ret(ii).beta_part_task_cortex_ret;

    for xx = 1:7
        beta_part_mean_network_ret(xx) = mean(beta_part_task_test_ret(schaef_id==xx)); % mean for each network
    end

storage_PD_DT_ret(ii).network_part_mean_ret = beta_part_mean_network_ret; % 
end

for ii = 1:length(storage_PD_DT_ret)
   beta_bold_task_test_ret = storage_PD_DT_ret(ii).beta_bold_task_cortex_ret;

    for xx = 1:7
        beta_bold_mean_network_ret(xx) = mean(beta_bold_task_test_ret(schaef_id==xx)); % mean for each network
    end

storage_PD_DT_ret(ii).network_bold_mean_ret = beta_bold_mean_network_ret; % 
end

%% make files combining network values for all groups
% HC-ST
%part
%step 1 transpose the array to a horizontal one
for ii = 1:length(storage_HC_ST_early)
    storage_HC_ST_early(ii).network_part_mean_early_hor = transpose(storage_HC_ST_early(ii).network_part_mean_early);
end

%step 2 stack the arrays in one matrix

early_part_network_HC_ST_all = zeros(17, 7);

for i = 1:17
    array = storage_HC_ST_early(i).network_part_mean_early_hor;
    early_part_network_HC_ST_all(i, :) = array;
end

%BOLD
%step 1 transpose the array to a horizontal one
for ii = 1:length(storage_HC_ST_early)
    storage_HC_ST_early(ii).network_bold_mean_early_hor = transpose(storage_HC_ST_early(ii).network_bold_mean_early);
end

%step 2 stack the arrays in one matrix

early_bold_network_HC_ST_all = zeros(17, 7);

for i = 1:17
    array = storage_HC_ST_early(i).network_bold_mean_early_hor;
    early_bold_network_HC_ST_all(i, :) = array;
end


% HC-DT
%part
%step 1 transpose the array to a horizontal one
for ii = 1:length(storage_HC_DT_early)
    storage_HC_DT_early(ii).network_part_mean_early_hor = transpose(storage_HC_DT_early(ii).network_part_mean_early);
end

%step 2 stack the arrays in one matrix

early_part_network_HC_DT_all = zeros(21, 7);

for i = 1:21
    array = storage_HC_DT_early(i).network_part_mean_early_hor;
    early_part_network_HC_DT_all(i, :) = array;
end

%BOLD
%step 1 transpose the array to a horizontal one
for ii = 1:length(storage_HC_DT_early)
    storage_HC_DT_early(ii).network_bold_mean_early_hor = transpose(storage_HC_DT_early(ii).network_bold_mean_early);
end

%step 2 stack the arrays in one matrix

early_bold_network_HC_DT_all = zeros(21, 7);

for i = 1:21
    array = storage_HC_DT_early(i).network_bold_mean_early_hor;
    early_bold_network_HC_DT_all(i, :) = array;
end

% HC-DT
%part
%step 1 transpose the array to a horizontal one
for ii = 1:length(storage_HC_DT_early)
    storage_HC_DT_early(ii).network_part_mean_early_hor = transpose(storage_HC_DT_early(ii).network_part_mean_early);
end

%step 2 stack the arrays in one matrix

early_part_network_HC_DT_all = zeros(21, 7);

for i = 1:21
    array = storage_HC_DT_early(i).network_part_mean_early_hor;
    early_part_network_HC_DT_all(i, :) = array;
end

%BOLD
%step 1 transpose the array to a horizontal one
for ii = 1:length(storage_HC_DT_early)
    storage_HC_DT_early(ii).network_bold_mean_early_hor = transpose(storage_HC_DT_early(ii).network_bold_mean_early);
end

%step 2 stack the arrays in one matrix

early_bold_network_HC_DT_all = zeros(21, 7);

for i = 1:21
    array = storage_HC_DT_early(i).network_bold_mean_early_hor;
    early_bold_network_HC_DT_all(i, :) = array;
end

% PD-ST
%step 1 transpose the array to a horizontal one
for ii = 1:length(storage_PD_ST_early)
    storage_PD_ST_early(ii).network_part_mean_early_hor = transpose(storage_PD_ST_early(ii).network_part_mean_early);
end

%step 2 stack the arrays in one matrix

early_part_network_PD_ST_all = zeros(18, 7);

for i = 1:18
    array = storage_PD_ST_early(i).network_part_mean_early_hor;
    early_part_network_PD_ST_all(i, :) = array;
end

%BOLD
%step 1 transpose the array to a horizontal one
for ii = 1:length(storage_PD_ST_early)
    storage_PD_ST_early(ii).network_bold_mean_early_hor = transpose(storage_PD_ST_early(ii).network_bold_mean_early);
end

%step 2 stack the arrays in one matrix

early_bold_network_PD_ST_all = zeros(18, 7);

for i = 1:18
    array = storage_PD_ST_early(i).network_bold_mean_early_hor;
    early_bold_network_PD_ST_all(i, :) = array;
end

%PD-DT
%part
%step 1 transpose the array to a horizontal one
for ii = 1:length(storage_PD_DT_early)
    storage_PD_DT_early(ii).network_part_mean_early_hor = transpose(storage_PD_DT_early(ii).network_part_mean_early);
end

%step 2 stack the arrays in one matrix

early_part_network_PD_DT_all = zeros(21, 7);

for i = 1:21
    array = storage_PD_DT_early(i).network_part_mean_early_hor;
    early_part_network_PD_DT_all(i, :) = array;
end

%BOLD
%step 1 transpose the array to a horizontal one
for ii = 1:length(storage_PD_DT_early)
    storage_PD_DT_early(ii).network_bold_mean_early_hor = transpose(storage_PD_DT_early(ii).network_bold_mean_early);
end

%step 2 stack the arrays in one matrix

early_bold_network_PD_DT_all = zeros(21, 7);

for i = 1:21
    array = storage_PD_DT_early(i).network_bold_mean_early_hor;
    early_bold_network_PD_DT_all(i, :) = array;
end

%% Late scans
% HC-ST
%part
%step 1 transpose the array to a horizontal one
for ii = 1:length(storage_HC_ST_late)
    storage_HC_ST_late(ii).network_part_mean_late_hor = transpose(storage_HC_ST_late(ii).network_part_mean_late);
end

%step 2 stack the arrays in one matrix

late_part_network_HC_ST_all = zeros(17, 7);

for i = 1:17
    array = storage_HC_ST_late(i).network_part_mean_late_hor;
    late_part_network_HC_ST_all(i, :) = array;
end

%BOLD
%step 1 transpose the array to a horizontal one
for ii = 1:length(storage_HC_ST_late)
    storage_HC_ST_late(ii).network_bold_mean_late_hor = transpose(storage_HC_ST_late(ii).network_bold_mean_late);
end

%step 2 stack the arrays in one matrix

late_bold_network_HC_ST_all = zeros(17, 7);

for i = 1:17
    array = storage_HC_ST_late(i).network_bold_mean_late_hor;
    late_bold_network_HC_ST_all(i, :) = array;
end


% HC-DT
%part
%step 1 transpose the array to a horizontal one
for ii = 1:length(storage_HC_DT_late)
    storage_HC_DT_late(ii).network_part_mean_late_hor = transpose(storage_HC_DT_late(ii).network_part_mean_late);
end

%step 2 stack the arrays in one matrix

late_part_network_HC_DT_all = zeros(21, 7);

for i = 1:21
    array = storage_HC_DT_late(i).network_part_mean_late_hor;
    late_part_network_HC_DT_all(i, :) = array;
end

%BOLD
%step 1 transpose the array to a horizontal one
for ii = 1:length(storage_HC_DT_late)
    storage_HC_DT_late(ii).network_bold_mean_late_hor = transpose(storage_HC_DT_late(ii).network_bold_mean_late);
end

%step 2 stack the arrays in one matrix

late_bold_network_HC_DT_all = zeros(21, 7);

for i = 1:21
    array = storage_HC_DT_late(i).network_bold_mean_late_hor;
    late_bold_network_HC_DT_all(i, :) = array;
end

% HC-DT
%part
%step 1 transpose the array to a horizontal one
for ii = 1:length(storage_HC_DT_late)
    storage_HC_DT_late(ii).network_part_mean_late_hor = transpose(storage_HC_DT_late(ii).network_part_mean_late);
end

%step 2 stack the arrays in one matrix

late_part_network_HC_DT_all = zeros(21, 7);

for i = 1:21
    array = storage_HC_DT_late(i).network_part_mean_late_hor;
    late_part_network_HC_DT_all(i, :) = array;
end

%BOLD
%step 1 transpose the array to a horizontal one
for ii = 1:length(storage_HC_DT_late)
    storage_HC_DT_late(ii).network_bold_mean_late_hor = transpose(storage_HC_DT_late(ii).network_bold_mean_late);
end

%step 2 stack the arrays in one matrix

late_bold_network_HC_DT_all = zeros(21, 7);

for i = 1:21
    array = storage_HC_DT_late(i).network_bold_mean_late_hor;
    late_bold_network_HC_DT_all(i, :) = array;
end

% PD-ST
%step 1 transpose the array to a horizontal one
for ii = 1:length(storage_PD_ST_late)
    storage_PD_ST_late(ii).network_part_mean_late_hor = transpose(storage_PD_ST_late(ii).network_part_mean_late);
end

%step 2 stack the arrays in one matrix

late_part_network_PD_ST_all = zeros(18, 7);

for i = 1:18
    array = storage_PD_ST_late(i).network_part_mean_late_hor;
    late_part_network_PD_ST_all(i, :) = array;
end

%BOLD
%step 1 transpose the array to a horizontal one
for ii = 1:length(storage_PD_ST_late)
    storage_PD_ST_late(ii).network_bold_mean_late_hor = transpose(storage_PD_ST_late(ii).network_bold_mean_late);
end

%step 2 stack the arrays in one matrix

late_bold_network_PD_ST_all = zeros(18, 7);

for i = 1:18
    array = storage_PD_ST_late(i).network_bold_mean_late_hor;
    late_bold_network_PD_ST_all(i, :) = array;
end

%PD-DT
%part
%step 1 transpose the array to a horizontal one
for ii = 1:length(storage_PD_DT_late)
    storage_PD_DT_late(ii).network_part_mean_late_hor = transpose(storage_PD_DT_late(ii).network_part_mean_late);
end

%step 2 stack the arrays in one matrix

late_part_network_PD_DT_all = zeros(20, 7);

for i = 1:20
    array = storage_PD_DT_late(i).network_part_mean_late_hor;
    late_part_network_PD_DT_all(i, :) = array;
end

%BOLD
%step 1 transpose the array to a horizontal one
for ii = 1:length(storage_PD_DT_late)
    storage_PD_DT_late(ii).network_bold_mean_late_hor = transpose(storage_PD_DT_late(ii).network_bold_mean_late);
end

%step 2 stack the arrays in one matrix

late_bold_network_PD_DT_all = zeros(20, 7);

for i = 1:20
    array = storage_PD_DT_late(i).network_bold_mean_late_hor;
    late_bold_network_PD_DT_all(i, :) = array;
end

%% retention scans
% HC-ST
%part
%step 1 transpose the array to a horizontal one
for ii = 1:length(storage_HC_ST_ret)
    storage_HC_ST_ret(ii).network_part_mean_ret_hor = transpose(storage_HC_ST_ret(ii).network_part_mean_ret);
end

%step 2 stack the arrays in one matrix

ret_part_network_HC_ST_all = zeros(17, 7);

for i = 1:17
    array = storage_HC_ST_ret(i).network_part_mean_ret_hor;
    ret_part_network_HC_ST_all(i, :) = array;
end

%BOLD
%step 1 transpose the array to a horizontal one
for ii = 1:length(storage_HC_ST_ret)
    storage_HC_ST_ret(ii).network_bold_mean_ret_hor = transpose(storage_HC_ST_ret(ii).network_bold_mean_ret);
end

%step 2 stack the arrays in one matrix

ret_bold_network_HC_ST_all = zeros(17, 7);

for i = 1:17
    array = storage_HC_ST_ret(i).network_bold_mean_ret_hor;
    ret_bold_network_HC_ST_all(i, :) = array;
end


% HC-DT
%part
%step 1 transpose the array to a horizontal one
for ii = 1:length(storage_HC_DT_ret)
    storage_HC_DT_ret(ii).network_part_mean_ret_hor = transpose(storage_HC_DT_ret(ii).network_part_mean_ret);
end

%step 2 stack the arrays in one matrix

ret_part_network_HC_DT_all = zeros(21, 7);

for i = 1:21
    array = storage_HC_DT_ret(i).network_part_mean_ret_hor;
    ret_part_network_HC_DT_all(i, :) = array;
end

%BOLD
%step 1 transpose the array to a horizontal one
for ii = 1:length(storage_HC_DT_ret)
    storage_HC_DT_ret(ii).network_bold_mean_ret_hor = transpose(storage_HC_DT_ret(ii).network_bold_mean_ret);
end

%step 2 stack the arrays in one matrix

ret_bold_network_HC_DT_all = zeros(21, 7);

for i = 1:21
    array = storage_HC_DT_ret(i).network_bold_mean_ret_hor;
    ret_bold_network_HC_DT_all(i, :) = array;
end

% HC-DT
%part
%step 1 transpose the array to a horizontal one
for ii = 1:length(storage_HC_DT_ret)
    storage_HC_DT_ret(ii).network_part_mean_ret_hor = transpose(storage_HC_DT_ret(ii).network_part_mean_ret);
end

%step 2 stack the arrays in one matrix

ret_part_network_HC_DT_all = zeros(21, 7);

for i = 1:21
    array = storage_HC_DT_ret(i).network_part_mean_ret_hor;
    ret_part_network_HC_DT_all(i, :) = array;
end

%BOLD
%step 1 transpose the array to a horizontal one
for ii = 1:length(storage_HC_DT_ret)
    storage_HC_DT_ret(ii).network_bold_mean_ret_hor = transpose(storage_HC_DT_ret(ii).network_bold_mean_ret);
end

%step 2 stack the arrays in one matrix

ret_bold_network_HC_DT_all = zeros(21, 7);

for i = 1:21
    array = storage_HC_DT_ret(i).network_bold_mean_ret_hor;
    ret_bold_network_HC_DT_all(i, :) = array;
end

% PD-ST
%step 1 transpose the array to a horizontal one
for ii = 1:length(storage_PD_ST_ret)
    storage_PD_ST_ret(ii).network_part_mean_ret_hor = transpose(storage_PD_ST_ret(ii).network_part_mean_ret);
end

%step 2 stack the arrays in one matrix

ret_part_network_PD_ST_all = zeros(18, 7);

for i = 1:18
    array = storage_PD_ST_ret(i).network_part_mean_ret_hor;
    ret_part_network_PD_ST_all(i, :) = array;
end

%BOLD
%step 1 transpose the array to a horizontal one
for ii = 1:length(storage_PD_ST_ret)
    storage_PD_ST_ret(ii).network_bold_mean_ret_hor = transpose(storage_PD_ST_ret(ii).network_bold_mean_ret);
end

%step 2 stack the arrays in one matrix

ret_bold_network_PD_ST_all = zeros(18, 7);

for i = 1:18
    array = storage_PD_ST_ret(i).network_bold_mean_ret_hor;
    ret_bold_network_PD_ST_all(i, :) = array;
end

%PD-DT
%part
%step 1 transpose the array to a horizontal one
for ii = 1:length(storage_PD_DT_ret)
    storage_PD_DT_ret(ii).network_part_mean_ret_hor = transpose(storage_PD_DT_ret(ii).network_part_mean_ret);
end

%step 2 stack the arrays in one matrix

ret_part_network_PD_DT_all = zeros(20, 7);

for i = 1:20
    array = storage_PD_DT_ret(i).network_part_mean_ret_hor;
    ret_part_network_PD_DT_all(i, :) = array;
end

%BOLD
%step 1 transpose the array to a horizontal one
for ii = 1:length(storage_PD_DT_ret)
    storage_PD_DT_ret(ii).network_bold_mean_ret_hor = transpose(storage_PD_DT_ret(ii).network_bold_mean_ret);
end

%step 2 stack the arrays in one matrix

ret_bold_network_PD_DT_all = zeros(20, 7);

for i = 1:20
    array = storage_PD_DT_ret(i).network_bold_mean_ret_hor;
    ret_bold_network_PD_DT_all(i, :) = array;
end

early_part_network_all = vertcat(early_part_network_HC_ST_all, early_part_network_HC_DT_all, early_part_network_PD_ST_all, early_part_network_PD_DT_all);
early_bold_network_all = vertcat(early_bold_network_HC_ST_all, early_bold_network_HC_DT_all, early_bold_network_PD_ST_all, early_bold_network_PD_DT_all);

late_part_network_all = vertcat(late_part_network_HC_ST_all, late_part_network_HC_DT_all, late_part_network_PD_ST_all, late_part_network_PD_DT_all);
late_bold_network_all = vertcat(late_bold_network_HC_ST_all, late_bold_network_HC_DT_all, late_bold_network_PD_ST_all, late_bold_network_PD_DT_all);

ret_part_network_all = vertcat(ret_part_network_HC_ST_all, ret_part_network_HC_DT_all, ret_part_network_PD_ST_all, ret_part_network_PD_DT_all);
ret_bold_network_all = vertcat(ret_bold_network_HC_ST_all, ret_bold_network_HC_DT_all, ret_bold_network_PD_ST_all, ret_bold_network_PD_DT_all);


%% make one matrix for early, late and retention

part_network_all = vertcat(early_part_network_all, late_part_network_all, ret_part_network_all);
bold_network_all = vertcat(early_bold_network_all, late_bold_network_all, ret_bold_network_all);

%% average per group per network

%for ii = 1:length(storage_HC_ST)
   % network_part_mean = storage_HC_ST(ii).network_part_mean;

  %  HCST_group_network_part_mean(:,:,ii) = network_part_mean;
 %   HCST_network_part_mean = mean(HCST_group_network_part_mean,3); % what does the 3 stand for?
%end



%% permutation tests on network level?

%[HCSTDT_netw_sig,HCSTDT_netw_pval] = perm_1d_delta(HCST_network_part_mean,HCDT_network_part_mean,100); % 100 = number of permutations
%[PDSTDT_netw_sig,PDSTDT_netw_pval] = perm_1d_delta(PDST_network_part_mean,PDDT_network_part_mean,100); % 100 = number of permutations
%[STPDHC_netw_sig,STPDHC_netw_pval] = perm_1d_delta(HCST_network_part_mean,PDST_network_part_mean,100); % 100 = number of permutations
%[DTPDHC_netw_sig,DTPDHC_netw_pval] = perm_1d_delta(HCDT_network_part_mean,PDDT_network_part_mean,100); % 100 = number of permutations


%for ii = 1:length(storage_HC_ST)
 %   network_part_mean(:,:,ii) = beta_part_mean_network_group;
%end

%beta_part_network_HCST = mean(beta_part_mean_network_group,3);



%% Make design matrix from data outputted above
% Overview of matrices:
% early_part_ROI, early_bold_ROI, early_part_netw, early_bold_netw
% late_part_ROI, late_bold_ROI, late_part_netw, late_bold_netw
% ret_part_ROI, ret_bold_ROI, ret_part_netw, ret_bold_netw

%% early scans

% early_part_ROI

% HC-ST
%part
%step 1 transpose the array to a horizontal one
for ii = 1:length(storage_HC_ST_early)
    storage_HC_ST_early(ii).beta_part_task_early_hor = transpose(storage_HC_ST_early(ii).beta_part_task_early);
end

%step 2 stack the arrays in one matrix

early_part_roi_HC_ST_all = zeros(17, 502);

for i = 1:17
    array = storage_HC_ST_early(i).beta_part_task_early_hor;
    early_part_roi_HC_ST_all(i, :) = array;
end

%BOLD
%step 1 transpose the array to a horizontal one
for ii = 1:length(storage_HC_ST_early)
    storage_HC_ST_early(ii).beta_bold_task_early_hor = transpose(storage_HC_ST_early(ii).beta_bold_task_early);
end

%step 2 stack the arrays in one matrix

early_bold_roi_HC_ST_all = zeros(17, 502);

for i = 1:17
    array = storage_HC_ST_early(i).beta_bold_task_early_hor;
    early_bold_roi_HC_ST_all(i, :) = array;
end


% HC-DT
%part
%step 1 transpose the array to a horizontal one
for ii = 1:length(storage_HC_DT_early)
    storage_HC_DT_early(ii).beta_part_task_early_hor = transpose(storage_HC_DT_early(ii).beta_part_task_early);
end

%step 2 stack the arrays in one matrix

early_part_roi_HC_DT_all = zeros(21, 502);

for i = 1:21
    array = storage_HC_DT_early(i).beta_part_task_early_hor;
    early_part_roi_HC_DT_all(i, :) = array;
end

%BOLD
%step 1 transpose the array to a horizontal one
for ii = 1:length(storage_HC_DT_early)
    storage_HC_DT_early(ii).beta_bold_task_early_hor = transpose(storage_HC_DT_early(ii).beta_bold_task_early);
end

%step 2 stack the arrays in one matrix

early_bold_roi_HC_DT_all = zeros(21, 502);

for i = 1:21
    array = storage_HC_DT_early(i).beta_bold_task_early_hor;
    early_bold_roi_HC_DT_all(i, :) = array;
end


% PD-ST
%part
%step 1 transpose the array to a horizontal one
for ii = 1:length(storage_PD_ST_early)
    storage_PD_ST_early(ii).beta_part_task_early_hor = transpose(storage_PD_ST_early(ii).beta_part_task_early);
end

%step 2 stack the arrays in one matrix

early_part_roi_PD_ST_all = zeros(18, 502);

for i = 1:18
    array = storage_PD_ST_early(i).beta_part_task_early_hor;
    early_part_roi_PD_ST_all(i, :) = array;
end

%BOLD
%step 1 transpose the array to a horizontal one
for ii = 1:length(storage_PD_ST_early)
    storage_PD_ST_early(ii).beta_bold_task_early_hor = transpose(storage_PD_ST_early(ii).beta_bold_task_early);
end

%step 2 stack the arrays in one matrix

early_bold_roi_PD_ST_all = zeros(18, 502);

for i = 1:18
    array = storage_PD_ST_early(i).beta_bold_task_early_hor;
    early_bold_roi_PD_ST_all(i, :) = array;
end


% PD-DT
%part
%step 1 transpose the array to a horizontal one
for ii = 1:length(storage_PD_DT_early)
    storage_PD_DT_early(ii).beta_part_task_early_hor = transpose(storage_PD_DT_early(ii).beta_part_task_early);
end

%step 2 stack the arrays in one matrix

early_part_roi_PD_DT_all = zeros(21, 502);

for i = 1:21
    array = storage_PD_DT_early(i).beta_part_task_early_hor;
    early_part_roi_PD_DT_all(i, :) = array;
end

%BOLD
%step 1 transpose the array to a horizontal one
for ii = 1:length(storage_PD_DT_early)
    storage_PD_DT_early(ii).beta_bold_task_early_hor = transpose(storage_PD_DT_early(ii).beta_bold_task_early);
end

%step 2 stack the arrays in one matrix

early_bold_roi_PD_DT_all = zeros(21, 502);

for i = 1:21
    array = storage_PD_DT_early(i).beta_bold_task_early_hor;
    early_bold_roi_PD_DT_all(i, :) = array;
end


%% Get all early together

early_part_all = vertcat(early_part_roi_HC_ST_all, early_part_roi_HC_DT_all, early_part_roi_PD_ST_all, early_part_roi_PD_DT_all);
early_bold_all = vertcat(early_bold_roi_HC_ST_all, early_bold_roi_HC_DT_all, early_bold_roi_PD_ST_all, early_bold_roi_PD_DT_all);

%% late scans


% HC-ST
%part
%step 1 transpose the array to a horizontal one
for ii = 1:length(storage_HC_ST_late)
    storage_HC_ST_late(ii).beta_part_task_late_hor = transpose(storage_HC_ST_late(ii).beta_part_task_late);
end

%step 2 stack the arrays in one matrix

late_part_roi_HC_ST_all = zeros(17, 502);

for i = 1:17
    array = storage_HC_ST_late(i).beta_part_task_late_hor;
    late_part_roi_HC_ST_all(i, :) = array;
end

%BOLD
%step 1 transpose the array to a horizontal one
for ii = 1:length(storage_HC_ST_late)
    storage_HC_ST_late(ii).beta_bold_task_late_hor = transpose(storage_HC_ST_late(ii).beta_bold_task_late);
end

%step 2 stack the arrays in one matrix

late_bold_roi_HC_ST_all = zeros(17, 502);

for i = 1:17
    array = storage_HC_ST_late(i).beta_bold_task_late_hor;
    late_bold_roi_HC_ST_all(i, :) = array;
end


% HC-DT
%part
%step 1 transpose the array to a horizontal one
for ii = 1:length(storage_HC_DT_late)
    storage_HC_DT_late(ii).beta_part_task_late_hor = transpose(storage_HC_DT_late(ii).beta_part_task_late);
end

%step 2 stack the arrays in one matrix

late_part_roi_HC_DT_all = zeros(21, 502);

for i = 1:21
    array = storage_HC_DT_late(i).beta_part_task_late_hor;
    late_part_roi_HC_DT_all(i, :) = array;
end

%BOLD
%step 1 transpose the array to a horizontal one
for ii = 1:length(storage_HC_DT_late)
    storage_HC_DT_late(ii).beta_bold_task_late_hor = transpose(storage_HC_DT_late(ii).beta_bold_task_late);
end

%step 2 stack the arrays in one matrix

late_bold_roi_HC_DT_all = zeros(21, 502);

for i = 1:21
    array = storage_HC_DT_late(i).beta_bold_task_late_hor;
    late_bold_roi_HC_DT_all(i, :) = array;
end


% PD-ST
%part
%step 1 transpose the array to a horizontal one
for ii = 1:length(storage_PD_ST_late)
    storage_PD_ST_late(ii).beta_part_task_late_hor = transpose(storage_PD_ST_late(ii).beta_part_task_late);
end

%step 2 stack the arrays in one matrix

late_part_roi_PD_ST_all = zeros(18, 502);

for i = 1:18
    array = storage_PD_ST_late(i).beta_part_task_late_hor;
    late_part_roi_PD_ST_all(i, :) = array;
end

%BOLD
%step 1 transpose the array to a horizontal one
for ii = 1:length(storage_PD_ST_late)
    storage_PD_ST_late(ii).beta_bold_task_late_hor = transpose(storage_PD_ST_late(ii).beta_bold_task_late);
end

%step 2 stack the arrays in one matrix

late_bold_roi_PD_ST_all = zeros(18, 502);

for i = 1:18
    array = storage_PD_ST_late(i).beta_bold_task_late_hor;
    late_bold_roi_PD_ST_all(i, :) = array;
end


% PD-DT
%part
%step 1 transpose the array to a horizontal one
for ii = 1:length(storage_PD_DT_late)
    storage_PD_DT_late(ii).beta_part_task_late_hor = transpose(storage_PD_DT_late(ii).beta_part_task_late);
end

%step 2 stack the arrays in one matrix

late_part_roi_PD_DT_all = zeros(20, 502);

for i = 1:20
    array = storage_PD_DT_late(i).beta_part_task_late_hor;
    late_part_roi_PD_DT_all(i, :) = array;
end

%BOLD
%step 1 transpose the array to a horizontal one
for ii = 1:length(storage_PD_DT_late)
    storage_PD_DT_late(ii).beta_bold_task_late_hor = transpose(storage_PD_DT_late(ii).beta_bold_task_late);
end

%step 2 stack the arrays in one matrix

late_bold_roi_PD_DT_all = zeros(20, 502);

for i = 1:20
    array = storage_PD_DT_late(i).beta_bold_task_late_hor;
    late_bold_roi_PD_DT_all(i, :) = array;
end


%% Get all late together

late_part_all = vertcat(late_part_roi_HC_ST_all, late_part_roi_HC_DT_all, late_part_roi_PD_ST_all, late_part_roi_PD_DT_all);
late_bold_all = vertcat(late_bold_roi_HC_ST_all, late_bold_roi_HC_DT_all, late_bold_roi_PD_ST_all, late_bold_roi_PD_DT_all);


%% retention scans


% HC-ST
%part
%step 1 transpose the array to a horizontal one
for ii = 1:length(storage_HC_ST_ret)
    storage_HC_ST_ret(ii).beta_part_task_ret_hor = transpose(storage_HC_ST_ret(ii).beta_part_task_ret);
end

%step 2 stack the arrays in one matrix

ret_part_roi_HC_ST_all = zeros(17, 502);

for i = 1:17
    array = storage_HC_ST_ret(i).beta_part_task_ret_hor;
    ret_part_roi_HC_ST_all(i, :) = array;
end

%BOLD
%step 1 transpose the array to a horizontal one
for ii = 1:length(storage_HC_ST_ret)
    storage_HC_ST_ret(ii).beta_bold_task_ret_hor = transpose(storage_HC_ST_ret(ii).beta_bold_task_ret);
end

%step 2 stack the arrays in one matrix

ret_bold_roi_HC_ST_all = zeros(17, 502);

for i = 1:17
    array = storage_HC_ST_ret(i).beta_bold_task_ret_hor;
    ret_bold_roi_HC_ST_all(i, :) = array;
end


% HC-DT
%part
%step 1 transpose the array to a horizontal one
for ii = 1:length(storage_HC_DT_ret)
    storage_HC_DT_ret(ii).beta_part_task_ret_hor = transpose(storage_HC_DT_ret(ii).beta_part_task_ret);
end

%step 2 stack the arrays in one matrix

ret_part_roi_HC_DT_all = zeros(21, 502);

for i = 1:21
    array = storage_HC_DT_ret(i).beta_part_task_ret_hor;
    ret_part_roi_HC_DT_all(i, :) = array;
end

%BOLD
%step 1 transpose the array to a horizontal one
for ii = 1:length(storage_HC_DT_ret)
    storage_HC_DT_ret(ii).beta_bold_task_ret_hor = transpose(storage_HC_DT_ret(ii).beta_bold_task_ret);
end

%step 2 stack the arrays in one matrix

ret_bold_roi_HC_DT_all = zeros(21, 502);

for i = 1:21
    array = storage_HC_DT_ret(i).beta_bold_task_ret_hor;
    ret_bold_roi_HC_DT_all(i, :) = array;
end


% PD-ST
%part
%step 1 transpose the array to a horizontal one
for ii = 1:length(storage_PD_ST_ret)
    storage_PD_ST_ret(ii).beta_part_task_ret_hor = transpose(storage_PD_ST_ret(ii).beta_part_task_ret);
end

%step 2 stack the arrays in one matrix

ret_part_roi_PD_ST_all = zeros(18, 502);

for i = 1:18
    array = storage_PD_ST_ret(i).beta_part_task_ret_hor;
    ret_part_roi_PD_ST_all(i, :) = array;
end

%BOLD
%step 1 transpose the array to a horizontal one
for ii = 1:length(storage_PD_ST_ret)
    storage_PD_ST_ret(ii).beta_bold_task_ret_hor = transpose(storage_PD_ST_ret(ii).beta_bold_task_ret);
end

%step 2 stack the arrays in one matrix

ret_bold_roi_PD_ST_all = zeros(18, 502);

for i = 1:18
    array = storage_PD_ST_ret(i).beta_bold_task_ret_hor;
    ret_bold_roi_PD_ST_all(i, :) = array;
end


% PD-DT
%part
%step 1 transpose the array to a horizontal one
for ii = 1:length(storage_PD_DT_ret)
    storage_PD_DT_ret(ii).beta_part_task_ret_hor = transpose(storage_PD_DT_ret(ii).beta_part_task_ret);
end

%step 2 stack the arrays in one matrix

ret_part_roi_PD_DT_all = zeros(20, 502);

for i = 1:20
    array = storage_PD_DT_ret(i).beta_part_task_ret_hor;
    ret_part_roi_PD_DT_all(i, :) = array;
end

%BOLD
%step 1 transpose the array to a horizontal one
for ii = 1:length(storage_PD_DT_ret)
    storage_PD_DT_ret(ii).beta_bold_task_ret_hor = transpose(storage_PD_DT_ret(ii).beta_bold_task_ret);
end

%step 2 stack the arrays in one matrix

ret_bold_roi_PD_DT_all = zeros(20, 502);

for i = 1:20
    array = storage_PD_DT_ret(i).beta_bold_task_ret_hor;
    ret_bold_roi_PD_DT_all(i, :) = array;
end


%% Get all ret together

ret_part_all = vertcat(ret_part_roi_HC_ST_all, ret_part_roi_HC_DT_all, ret_part_roi_PD_ST_all, ret_part_roi_PD_DT_all);
ret_bold_all = vertcat(ret_bold_roi_HC_ST_all, ret_bold_roi_HC_DT_all, ret_bold_roi_PD_ST_all, ret_bold_roi_PD_DT_all);


%% make one matrix for early, late and retention

part_all = vertcat(early_part_all, late_part_all, ret_part_all);
bold_all = vertcat(early_bold_all, late_bold_all, ret_bold_all);


%%% use these files for the second-level analyses (see second-level.m)