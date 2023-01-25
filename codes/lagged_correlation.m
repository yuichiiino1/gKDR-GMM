function lagged_correlation(SGE_TASK_ID_TEXT) % input argument: sampleID

global maxlag rmaxthreshold realvec

sampleID = str2num(SGE_TASK_ID_TEXT);

%%%%%%%%%%%%%%%%%%%%%  preset parameters %%%%%%%%%%%%%%%%%%%

maxlag=20; %25; % 75;
rmaxthreshold = 0.3;
%allowedlag=20; %50;

link = 'indirect';

K = 4;
kGMM = 2;
embed_width = 30;
embed_step = 10;
freerun_repeat = 3;

project_folder = '/home/iino/gKDR_GMM_50mM/Linked_gKDR_freerun/publish_version';
gKDR_codes_folder = fullfile(project_folder, 'gKDR_codes');
utilities_folder = fullfile(project_folder, 'utilities'); %'../'
%metadata_folder = fullfile(project_folder, 'metadata'); %'../'
models_folder = fullfile(project_folder, 'models'); % ['../ver98_HLong/model_indirect_k' num2str(embed_width) '_tau' num2str(embed_step) '_kGMM' num2str(kGMM)];
common_cell_order_folder = fullfile(project_folder, 'common_cell_order'); %'../ver92_Hbase/common_cell_order';
simulation_results_folder = fullfile(project_folder, 'simulation_results'); %['../ver98_HLong/figures_' link '_k' num2str(embed_width) '_tau' num2str(embed_step) '_K' num2str(K) '_kGMM' num2str(kGMM) '/sample' num2str(sampleID)];
lagged_correlation_folder = fullfile(project_folder, 'lagged_correlation');

%%%%%  Initial procedure  %%%%%

addpath(gKDR_codes_folder);
addpath(utilities_folder);

close all
disp(['sample ' num2str(sampleID)])

if ~exist(lagged_correlation_folder, 'dir')
    mkdir(lagged_correlation_folder)
end
savefileheader = fullfile(lagged_correlation_folder, ['sample' num2str(sampleID)]); % ; sampleID

model_folder = fullfile(models_folder, ['model_' link '_k' num2str(embed_width) '_tau' num2str(embed_step) '_kGMM' num2str(kGMM)]);
modelfileheader = [model_folder '/sample' num2str(sampleID) '_K' num2str(K)];
modeldatafile = [modelfileheader '_modeldata.mat'];
paramfile = [modelfileheader '_param.mat'];
load(modeldatafile, 'uniqNames','targetcells','data'); %'selicell','colscell','Bcell','gmcell','train_span','target_train_all','source_train_all'
load(paramfile, 'time_step'); %'nahead', 'prediction_method'
  
simulation_result_subfolder = fullfile(simulation_results_folder, [link '_k' num2str(embed_width) '_tau' num2str(embed_step) '_K' num2str(K) '_kGMM' num2str(kGMM)]);
simulationresultfileheader = fullfile(simulation_result_subfolder, ['sample' num2str(sampleID)]);
        
disp([num2str(length(uniqNames)) ' names, ' num2str(size(data,2)) ' data'])

n = size(data, 1);
disp(['n=',num2str(n)])
Mt = length(targetcells);
targetcellnames = uniqNames(targetcells);

load(fullfile(common_cell_order_folder, ['sample' num2str(sampleID) '_common_outperm.mat']), 'outperm')

data2 = data(1:time_step:n, targetcells); % sequence of time indices for real data

lab=cellstr(targetcellnames(outperm));
labnum = cell(1,Mt);
for celli =1:Mt
    labnum{celli} = [strrep(lab{celli}, '_', '\_')]; %  ' ' num2str(celli)
end


%%%%%%%%%%%%%%%%   correlation matrix plot  %%%%%%%%%%%%%%%%

% for real data
[rplainmatrix_real,rmaxmatrix_real,rlagmatrix_real] = process_lagged_correlation(...
    data2, outperm, sampleID, savefileheader, Mt, labnum, 'real', []);

rplainmatrix_pred_cell = cell(1,freerun_repeat);
rmaxmatrix_pred_cell = cell(1,freerun_repeat);
rlagmatrix_pred_cell = cell(1,freerun_repeat);

realvec = reshape(rlagmatrix_real+tril(NaN(Mt),0),1,[]);
realvec = realvec(~isnan(realvec));

for testi = 1:freerun_repeat
    disp(['repeat ' num2str(testi)])
    
    load([simulationresultfileheader '_freerun' num2str(testi) '_savedata.mat'], 'uniqNames','targetcells','colscell','sampleID','targetcellnames','Bcell','gmcell','saltdata', 'selicell','Tc_real','Tc_all','real_predict', 'outperm', 'rmatrix','train_span','target_train_all','source_train_all') % includes real_predict, Tc_all, Tc_real

    nx = length(Tc_all); %step numbers for all
    nc = length(Tc_real); %step numbers for real
    
    % for predicted data
    
    data_pred = real_predict((nc+1):nx, 2:end);%freerun only
    [rplainmatrix_pred, rmaxmatrix_pred, rlagmatrix_pred] = process_lagged_correlation(...
        data_pred, outperm, sampleID, savefileheader, Mt, labnum, ['freerun' num2str(testi)], rmaxmatrix_real);

    rplainmatrix_pred_cell{testi} = rplainmatrix_pred;
    rmaxmatrix_pred_cell{testi} = rplainmatrix_pred;
    rlagmatrix_pred_cell{testi} = rlagmatrix_pred;
    
end

save([lagged_correlation_folder '/lagged_correlation_sample' num2str(sampleID) '.mat'], ...
    'rplainmatrix_real', 'rmaxmatrix_real', 'rlagmatrix_real', 'rplainmatrix_pred_cell', ...
    'rmaxmatrix_pred_cell', 'rlagmatrix_pred_cell', 'Mt', 'targetcells', 'outperm', 'labnum')

end

function  [rplainmatrix,rmaxmatrix,rlagmatrix] = process_lagged_correlation(data2, outperm, sampleID, savefileheader, Mt, labnum, prefix, rmaxmatrix_real)

global maxlag realvec

% real correlation matrix
figure;
rplainmatrix = corrcoef(data2(:,outperm),'rows','complete'); % omit NaN
cormatrix_plot(rplainmatrix, ['sample' num2str(sampleID) prefix '_plain'], [prefix '_plain'], savefileheader, Mt, labnum);%,sampleID, K, lambda);

[rmaxmatrix,rlagmatrix] = calculate_lagged_correlation(data2(:,outperm), rplainmatrix, rmaxmatrix_real);

% lagged correlation matrix
figure;
cormatrix_plot(rmaxmatrix, ['sample' num2str(sampleID) prefix '_max'], [prefix '_max'], savefileheader, Mt, labnum);

pval = 0;     %%%%% only when predicted results are treated, compare lag matrix with real by Spearman test 
if ~strcmp(prefix, 'real')
    predvec = reshape(rlagmatrix+tril(NaN(size(rlagmatrix)),0),1,[]);
    predvec = predvec(~isnan(predvec));
    [rho,pval] = corr(realvec', predvec', 'Type', 'Spearman');
end
    
% lag matrix
figure;
%lagmatrix_plot(rlagmatrix, ['sample' num2str(sampleID) prefix], [prefix], savefileheader, Mt, labnum, maxlag);
lagmatrix_plot(rlagmatrix, pval, [prefix], savefileheader, Mt, labnum, maxlag);

end

function [rmaxmatrix,rlagmatrix] = calculate_lagged_correlation(data,rmatrix_plain, rmaxmatrix_real)

global maxlag rmaxthreshold

Mt = size(data,2);
rmaxmatrix = rmatrix_plain;
rlagmatrix = NaN(Mt,Mt);

for i = 1:Mt-1
    for j = i+1:Mt
        r = xcorr(squeeze(data(:,i)), squeeze(data(:,j)), maxlag, 'normalized');
        [rm, ri] = max(abs(r));
        if isnan(ri)
           disp('ri is NaN. rm and r =')
           disp(rm)
           disp(r)
        end
        rilag = ri-maxlag-1;
        %if abs(rilag)<=allowedlag
        rmaxmatrix(i,j) = r(ri);
        rmaxmatrix(j,i) = r(ri);
        if isempty(rmaxmatrix_real) % if real
            if abs(r(ri)) > rmaxthreshold
                rlagmatrix(i,j) = rilag;
                rlagmatrix(j,i) = -rilag;
            end
        else                      % if freerun, based on real rmatrix
            if abs(rmaxmatrix_real(i,j)) > rmaxthreshold
                rlagmatrix(i,j) = rilag;
                rlagmatrix(j,i) = -rilag;
            end
        end
        %{
        if abs(rmatrix_plain(i,j))>abs(r(ri))+0.01
            disp(r)
            disp(ri)
            disp(r(ri))
            disp(rmatrix_plain(i,j))
            disp('')
        end
        %}

    end
end
for i = 1:Mt
    rmaxmatrix(i,i)=1;
    rlagmatrix(i,i)=0;
end

end
