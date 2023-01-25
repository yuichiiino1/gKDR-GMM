
function gKDR_GMM_connection_strength(UGE_TASK_ID_text) % input argument: sampleID

% based on generated models and real data, 
% caclulate mean gradient of expectation for real data (training data)

sampleID = str2num(UGE_TASK_ID_text);
Ks = 3:5;
link = 'indirect';
kGMM = 2; % number of gaussians for GMM estimation
embed_width = 30; %10; % number of embed_step's used for embedding = column number in source data    % need
embed_step = 10; %10; % invervals used for embedding (index-based)                % need

%%%%%%%%%%%%%%%%%%%%%  preset parameters %%%%%%%%%%%%%%%%%%%

project_folder = '/home/iino/gKDR_GMM_50mM/Linked_gKDR_freerun/publish_version';
gKDR_codes_folder = fullfile(project_folder, 'gKDR_codes');
utilities_folder = fullfile(project_folder, 'utilities'); %'../'
%metadata_folder = fullfile(project_folder, 'metadata'); %'../'
models_folder = fullfile(project_folder, 'models'); % ['../ver98_HLong/model_indirect_k' num2str(embed_width) '_tau' num2str(embed_step) '_kGMM' num2str(kGMM)];
%common_cell_order_folder = fullfile(project_folder, 'common_cell_order'); %'../ver92_Hbase/common_cell_order';
%simulation_results_folder = fullfile(project_folder, 'simulation_results'); %['../ver98_HLong/figures_' link '_k' num2str(embed_width) '_tau' num2str(embed_step) '_K' num2str(K) '_kGMM' num2str(kGMM) '/sample' num2str(sampleID)];


addpath(gKDR_codes_folder);
addpath(utilities_folder);

tic
disp('start')
disp(['link type = ' link])

for Ki = 1:length(Ks)  % dimension of dimensional reduction by gKDR
    K = Ks(Ki);
    
    model_subfolder = fullfile(models_folder, ['model_' link '_k' num2str(embed_width) '_tau' num2str(embed_step) '_kGMM' num2str(kGMM)]);
    modelfileheader = fullfile(model_subfolder, ['sample' num2str(sampleID) '_K' num2str(K)]);
    modeldatafile = [modelfileheader '_modeldata.mat'];
    %paramfile = [modelfileheader '_param.mat'];
    load(modeldatafile, 'uniqNames','targetcells','selicell','Bcell','gmcell','source_train_all'); %,'colscell','train_span','data','target_train_all'
    %load(paramfile, 'time_step','nahead', 'prediction_method');
    %{
    modelfolder = ['../ver98_HLong/model_indirect_k30_tau10_kGMM' num2str(kGMM)];
    modelfileheader = [modelfolder '/sample' num2str(sampleID) '_K' num2str(K0)];
    modeldatafile = [modelfileheader '_modeldata.mat'];
    paramfile = [modelfileheader '_param.mat'];
    load(modeldatafile);
    load(paramfile);
    %}
    
    disp(['<< sample ' num2str(sampleID) ' >>']);
    disp(['< K= ' num2str(K) ' >']);
    
    Mt = length(targetcells);
    n = size(source_train_all,1);
    
    assert(length(selicell)==Mt)
    assert(length(Bcell)==Mt)
    assert(length(gmcell)==Mt)
    gm_grad_cell = cell(1,Mt);
    gm_grad_mean_cell = cell(1,Mt);
    embed_size = embed_width;
    gm_grad_array = NaN(n, embed_size, Mt+1, Mt);
    
    for targeti = 1:Mt
        
        disp(['target ' num2str(targeti)])
        B = Bcell{targeti};
        gm = gmcell{targeti};
        seli = selicell{targeti};
        
        cols = reshape(repmat(seli*embed_size, embed_size, 1) + repmat((1:embed_size)',1,length(seli)),[],1)';
        source_train = source_train_all(:, cols);
        gm_grad = GMMgradient(source_train, gm, B);
        gm_grad_cell{targeti} = gm_grad;
        assert(all(size(gm_grad)==[n,length(seli)*embed_size]))
        gm_grad_array(:,:,seli+1, targeti) = reshape(gm_grad, n, embed_size, length(seli));  % from old to now
        gm_grad_mean_cell{targeti} = squeeze(nanmean(gm_grad,1));
    end
    gm_grad_mean = squeeze(nanmean(gm_grad_array,1));
    
    
    save([modelfileheader '_gradient.mat'], 'gm_grad_cell', 'gm_grad_array','gm_grad_mean','gm_grad_mean_cell','uniqNames','targetcells','selicell')
    disp(['saved to ' modelfileheader '_gradient.mat'])   % save to ver98_HLong/
    
end
end

