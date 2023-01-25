
function gKDR_GMM_makemodel_checker()

% checks completion of gKDR_GMM_makemodel.m

samples = 1:24;
Ks = 3:5;
kGMM = 2;
link = 'indirect';

embed_width = 30; % number of embed_step's used for embedding = column number in source data
embed_step = 10; % invervals used for embedding (index-based)

current_folder = '/home/iino/gKDR_GMM_50mM/Linked_gKDR_freerun/ver104_scrambled_network';
models_folder = fullfile(current_folder, 'models'); % ['../ver98_HLong/model_indirect_k' num2str(embed_width) '_tau' num2str(embed_step) '_kGMM' num2str(kGMM)];

disp('start')
disp(['link type = ' link])

for samplei = 1:length(samples)
    sampleID = samples(samplei);
    
    for Ki = 1:length(Ks)
        K = Ks(Ki);
        
        model_subfolder = fullfile(models_folder, ['model_' link '_k' num2str(embed_width) '_tau' num2str(embed_step) '_kGMM' num2str(kGMM)]);
        modelfileheader = fullfile(model_subfolder, ['sample' num2str(sampleID) '_K' num2str(K)]);
        
        try
            load([modelfileheader '_modeldata.mat'])
            load([modelfileheader '_param.mat'])
            disp(['sample' num2str(sampleID) ' K' num2str(K) ':  done']);
        catch
            disp(['sample' num2str(sampleID) ' K' num2str(K) ':  MISSING or BROKEN!'])
        end

    end
end
end

