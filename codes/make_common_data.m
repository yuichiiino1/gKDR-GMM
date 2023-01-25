% make_common_data.m

% makes data, targetcells, targetcellnames, Mt
% and save them along with uniqNames, autocorrthreshold, autocorrlag
% in one file, [common_data_folder]/common_data/samplex_data.mat

% requirements
% x_ratio.csv and x_uniqNames.csv in datafolder and conneurons.csv in metadata_folder
% where x stands for sample number.


%%%%%  Preset variables  %%%%%

samples = 1:24;
project_folder = '/home/iino/gKDR_GMM_50mM/Linked_gKDR_freerun/publish_version';
data_folder = '/home/iino/gKDR_GMM_50mM/cleandata_smoothened2/'; %'../../cleandata_smoothened2/'
common_data_folder = fullfile(project_folder, 'common_data'); %'../ver92_Hbase/common_data';
utilities_folder = fullfile(project_folder, 'utilities'); %'../'
metadata_folder = fullfile(project_folder, 'metadata'); %'../'
autocorrthreshold = 0.3;
autocorrlag = 20;
ls


%%%%%  Initial treatment  %%%%%

if ~exist(common_data_folder, 'dir')
    mkdir(common_data_folder)
end

conNames = table2cell(readtable(fullfile(metadata_folder, 'conneurons.csv'),'ReadVariableNames',false));

nsamples = length(samples);
addpath(utilities_folder);

%%%%% load data for each sample %%%%%

for sampleID = samples
    
    disp(['< sample ' num2str(sampleID) ' >'])
    
    data = csvread(fullfile(data_folder, [num2str(sampleID) '_ratio.csv']),0,0);
    n = size(data,1);
    cell_num = size(data,2);
    %disp([n cell_num])
    uniqNames = table2cell(readtable(fullfile(data_folder, [num2str(sampleID) '_uniqNames.csv']),'ReadVariableNames',false,'Format','%s'));
    
    % need to get targetcells and targetcellnames
    rall = diag(corr(data(1+autocorrlag:end,:),data(1:end-autocorrlag,:)));
    targetcells = find(rall > autocorrthreshold);
    targetcellnames = uniqNames(targetcells);
    targetcells = [];   % numbers in uniqNames, select only valid cells。conNamesに含まれるかどうかを確認
    for i=1:length(targetcellnames)
        if any(strcmp(conNames,targetcellnames{i}))
            targetcells = [targetcells; find(strcmp(uniqNames,targetcellnames{i}))]; %79:80; %72; %[]; %10:20; %190;
        end
    end
    
    Mt = length(targetcells);
    targetcellnames = uniqNames(targetcells);
    %disp(length(targetcellnames))
    disp(['n = ' num2str(n) ', cell_num = ' num2str(cell_num) ', Mt = ' num2str(Mt)])
    
    save([common_data_folder '/sample' num2str(sampleID) '_data.mat'], 'data', 'uniqNames', 'targetcellnames', 'targetcells', 'Mt', 'autocorrthreshold', 'autocorrlag')
    
end


