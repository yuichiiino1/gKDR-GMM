
function gKDR_GMM_connection_strength_indirect_cross_K_offset(UGE_TASK_ID_text)


sampleID = str2num(UGE_TASK_ID_text);
partition_number = 3;
Ks = 3:5;
link = 'indirect';
offsets = 0:4;
kGMMs = 2; % number of gaussians for GMM estimation

%%%%%%%%%%%%%%%%%%%%%  preset parameters %%%%%%%%%%%%%%%%%%%

project_folder = '/home/iino/gKDR_GMM_50mM/Linked_gKDR_freerun/publish_version';
%project_folder = '/Users/iino/Documents/___Transcend J/_4DImaging_Data_Analysis/gKDR-GMM/gKDR_GMM_50mM/Linked_gKDR_freerun/publish_version/';
%common_data_folder = fullfile(project_folder, 'common_data'); %'../ver92_Hbase/common_data';
gKDR_codes_folder = fullfile(project_folder, 'gKDR_codes');
utilities_folder = fullfile(project_folder, 'utilities'); %'../'
%metadata_folder = fullfile(project_folder, 'metadata'); %'../'
%common_cell_order_folder = fullfile(project_folder, 'common_cell_order');
models_folder = fullfile(project_folder, 'models'); % ['../ver98_HLong/model_indirect_k' num2str(embed_width) '_tau' num2str(embed_step) '_kGMM' num2str(kGMM)];

embed_width = 30; % number of embed_step's used for embedding = column number in source data
embed_step = 10; % invervals used for embedding (index-based)
%time_step = 5; % timestep invervals in which estimation is performed
%nahead = 5;    % how far ahead is the target for estimation

addpath(gKDR_codes_folder);
addpath(utilities_folder);

tic
disp('start')
disp(['link type = ' link])
disp(['<< sample ' num2str(sampleID) ' >>']);

for Ki = 1:length(Ks)  % dimension of dimensional reduction by gKDR
    K0 = Ks(Ki);
    for kGMM = kGMMs
        for offset = offsets
        
        for partitioni = 1:partition_number
            
            disp(['< K= ' num2str(K0) ' >']);
            disp(['< kGMM ' num2str(kGMM) ' >']);
            disp(['< offset ' num2str(offset) ' >']);
            disp(['< part ' num2str(partitioni) ' >']);
            modelfolder = fullfile(models_folder, ['model_' link '_k' num2str(embed_width) '_tau' num2str(embed_step) '_kGMM' num2str(kGMM) '_cross' num2str(partition_number)]);
            %modelfileheader = fullfile(modelfolder, ['sample' num2str(sampleID) '_K' num2str(K0) '_part' num2str(partitioni)]);
            modelfileheader = fullfile(modelfolder, ['sample' num2str(sampleID) '_K' num2str(K0) '_offset' num2str(offset) '_part' num2str(partitioni)]);
            modeldatafile = [modelfileheader '_modeldata.mat'];
            paramfile = [modelfileheader '_param.mat'];
            load(modeldatafile);
            load(paramfile);

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
                %disp(gm)
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
end
end

