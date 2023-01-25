% only neurons that are included in connection table are
% selected.

samples = 1:24;

nsamples = length(samples);

project_folder = '/Volumes/Transcend J/_4DImaging_Data_Analysis/gKDR-GMM/gKDR_GMM_50mM/Linked_gKDR_freerun/ver92_Hbase'; 
               %'/home/iino/gKDR_GMM_50mM/Linked_gKDR_freerun/publish_version';
common_data_folder = fullfile(project_folder, 'common_data'); %'../ver92_Hbase/common_data';
%gKDR_codes_folder = fullfile(project_folder, 'gKDR_codes');
%utilities_folder = fullfile(project_folder, 'utilities'); %'../'
%metadata_folder = fullfile(project_folder, 'metadata'); %'../'
%common_cell_order_folder = fullfile(project_folder, 'common_cell_order'); %'../ver92_Hbase/common_cell_order';
%models_folder = fullfile(project_folder, 'models'); % ['../ver98_HLong/model_indirect_k' num2str(embed_width) '_tau' num2str(embed_step) '_kGMM' num2str(kGMM)];


%datafolder = '../../gKDR/cleandata_smoothened2/';
%conNames = table2cell(readtable('../conneurons.csv','ReadVariableNames',false));

globalNames = [];
for sampleID = samples
    load(fullfile(common_data_folder,['sample' num2str(sampleID) '_data.mat']), 'uniqNames', 'targetcellnames', 'targetcellnames', 'Mt');
    
    disp(['<< sample' num2str(sampleID) ' >> ' num2str(Mt) ' neurons'])
    for Ni = 1:length(targetcellnames)
        if isempty(str2num(targetcellnames{Ni}))
            if ~any(strcmp(targetcellnames{Ni}, globalNames))
                globalNames = [globalNames, targetcellnames(Ni)];
            end
        end
    end
end
globalNames = sort(globalNames);
nonphaNames = [];
for ci = 1:length(globalNames)
    name = globalNames{ci};
    if isempty(regexp(name, '^I\d.*')) && isempty(regexp(name, '^M.*')) && isempty(regexp(name, '^NSM.+'))
        nonphaNames = [nonphaNames, {name}];
    end
end


disp(['globalNames ' num2str(length(globalNames)) ' neurons'])
disp(['non-pharyngeal nonphaNames ' num2str(length(nonphaNames)) ' neurons'])

writetable(table(globalNames'), 'globalNames.csv', 'WriteVariableNames',false);
writetable(table(nonphaNames'), 'nonphaNames.csv', 'WriteVariableNames',false);
