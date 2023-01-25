% evaluates and summarizes simulation results with different hyperparameters
% :different Ks in this code.
% originally: evaluate_models_indirect_k30_tau10.m


%%%%%%%%%%%%%%%%%%%%%  preset parameters %%%%%%%%%%%%%%%%%%%

samples = 1:24;  % calculate sample mean and each sample
Ks = 1:10;       % K survey
kGMM = 2;
freerun_repeat = 3;
link = 'indirect';
embed_width = 30; % number of embed_step's used for embedding = column number in source data  = k
embed_step = 10; % invervals used for embedding (index-based)    = tau

show_cormatrix = false;
show_heatmap = false;

project_folder = '/home/iino/gKDR_GMM_50mM/Linked_gKDR_freerun/publish_version';
gKDR_codes_folder = fullfile(project_folder, 'gKDR_codes');
utilities_folder = fullfile(project_folder, 'utilities'); %'../'
metadata_folder = fullfile(project_folder, 'metadata'); %'../'
models_folder = fullfile(project_folder, 'models'); % ['../ver98_HLong/model_indirect_k' num2str(embed_width) '_tau' num2str(embed_step) '_kGMM' num2str(kGMM)];
common_cell_order_folder = fullfile(project_folder, 'common_cell_order'); %'../ver92_Hbase/common_cell_order';
simulation_results_folder = fullfile(project_folder, 'simulation_results'); %['../ver98_HLong/figures_' link '_k' num2str(embed_width) '_tau' num2str(embed_step) '_K' num2str(K) '_kGMM' num2str(kGMM) '/sample' num2str(sampleID)];
evaluation_result_folder = fullfile(project_folder, 'evaluation_results');

%%%%  Initial procedure  %%%%%

addpath(gKDR_codes_folder);
addpath(utilities_folder);

evaluation_result_subfolder = fullfile(evaluation_result_folder, ['results_evaluate_models_' link '_k' num2str(embed_width) '_tau' num2str(embed_step)]);
if ~exist(evaluation_result_subfolder, 'dir')
    mkdir(evaluation_result_subfolder)
end
%evaluation_result_folder = ['./results_evaluate_models_' link '_k30_tau10'];
resultfiletitle1 = fullfile(evaluation_result_subfolder, '');


cormatch = NaN(freerun_repeat,length(samples),length(Ks));

for samplei = 1:length(samples)
    
    sampleID = samples(samplei);
    disp(['sample ' num2str(sampleID)]) % show sample number being processed
    
    for Ki = 1:length(Ks)
        K = Ks(Ki); % dimension of dimentional reduction by gKDR
        
        disp('')
        disp(['K ' num2str(K)])
        
        %testtitle = ['_K' num2str(K)];%[_K3]みたいな
        %simulationresultfolder = ['./figures_' link '_k' num2str(embed_width) '_tau' num2str(embed_step) '_K' num2str(K) '_kGMM3/sample' num2str(sampleID)];
        %simulationresultfileheader = [simulationresultfolder '/sample' num2str(sampleID)];
        simulation_result_subfolder = fullfile(simulation_results_folder, [link '_k' num2str(embed_width) '_tau' num2str(embed_step) '_K' num2str(K) '_kGMM' num2str(kGMM)]);
        simulationresultfileheader = fullfile(simulation_result_subfolder, ['sample' num2str(sampleID)]);
        
        model_subfolder = fullfile(models_folder, ['model_' link '_k' num2str(embed_width) '_tau' num2str(embed_step) '_kGMM' num2str(kGMM)]);
        modelfileheader = fullfile(model_subfolder, ['sample' num2str(sampleID) '_K' num2str(K)]);
        modeldatafile = [modelfileheader '_modeldata.mat'];
        paramfile = [modelfileheader '_param.mat'];
        load(modeldatafile, 'uniqNames','targetcells','selicell','colscell','Bcell','gmcell','train_span','data','target_train_all','source_train_all');
        load(paramfile, 'time_step','nahead', 'prediction_method');
        
        close all
        n = size(data,1);
        
        resultfiletitle2 = [resultfiletitle1 'sample' num2str(sampleID) '_K' num2str(K)];
        %load(['./sample' num2str(sampleID) '_savedata/' resultfiletitle2 '_savedata.mat'])%[sample1_K3_lambda0.01_savedata.mat]を読み込む
        
        sample_name = ['sample' num2str(sampleID)];
        freerun_start = train_span(2)+1;
        freestart = freerun_start;
        if freestart > n+1
            freestart = n+1;
        end

        
        load(fullfile(common_cell_order_folder, ['sample' num2str(sampleID) '_common_outperm.mat']), 'outperm')

        for testi = 1:freerun_repeat %repeat numer
            
            load([simulationresultfileheader '_freerun' num2str(testi) '_savedata.mat'], 'uniqNames','targetcells','colscell','sampleID','targetcellnames','Bcell','gmcell','saltdata', 'selicell','Tc_real','Tc_all','real_predict', 'outperm', 'rmatrix','train_span','target_train_all','source_train_all') % includes real_predict, Tc_all, Tc_real

            nx = length(Tc_all); %number of steps in all
            nc = length(Tc_real); %number of steps in real
            if size(real_predict, 2) <= 1
                disp(['number of cells for sample' num2str(sampleID) ' is too small (' num2str(size(real_predict, 2)) ')'])
            else
                use_salt_input = size(real_predict,2)-length(selicell);  % 0 or 1
                data2 = real_predict(1:nc, (use_salt_input+1):end);          % Real part only, remove 'salt'
                cell_num = size(data2,2);
                data_pred = real_predict((nc+1):nx, (use_salt_input+1):end);  %freerun only
                
                
                if ~isempty(outperm)
                    
                    lab=cellstr(targetcellnames(outperm));
                    labnum = cell(1,cell_num);
                    for celli =1:cell_num
                        labnum{celli} = [strrep(lab{celli}, '_', '\_')];
                    end
                    
                    %%%%%%%%%%%%%%%%   correlation matrix plot  %%%%%%%%%%%%%%%%
                    % for real data
                    rmatrix1 = corrcoef(data2(:,outperm),'rows','complete'); % omit NaN

                    if show_cormatrix
                        figure;
                        cormatrix_plot(rmatrix1, [sample_name '-real'], 'real', resultfiletitle2, cell_num, labnum);
                    end
                    
                    % for predicted data
                    rmatrix2 = corrcoef(data_pred(:,outperm),'rows','complete'); % omit NaN

                    if show_cormatrix
                        figure;
                        cormatrix_plot(rmatrix2, [sample_name '-freerun' num2str(testi)], ['freerun' num2str(testi)],  resultfiletitle2, cell_num, labnum);
                    end
                    
                    %score = sum(sum(rmatrix1.*rmatrix2.*triu(ones(cell_num,cell_num))))/cell_num/(cell_num-1); % mean(r1*r2) for upper triangle
                    score = sum(sum(((rmatrix1-rmatrix2).*triu(ones(cell_num,cell_num))).^2))/((cell_num)*(cell_num+1)/2);
                    
                    cormatch(testi, samplei, Ki) = score; % record score for summary plots
                    
                    
                    %%%%%%%%%%%%%%%%%   heatmap  %%%%%%%%%%%%%%%%
                    if show_heatmap
                        data_all = [saltdata(Tc_all)' [data2; data_pred]];
                        figure;
                        colormap jet(64);
                        image(data_all(:,[1 outperm+1])'*10+32);
                        title(sample_name,'Interpreter' ,'none');
                        colorbar;
                        h = gca;
                        set(h,'ytick',(1:(1+cell_num))');
                        set(h,'yticklabel',[{'salt'},labnum]);
                        set(h,'fontsize',6);
                        set(h,'Ticklength',[0 0]);
                        set(h,'xtick',(1:100:nx)');
                        set(h,'xticklabel',Tc_all(1:100:end)-1);
                        hold on
                        plot([nc+0.5 nc+0.5], [0 cell_num], '-', 'color',[0.5 0.5 0.5]);
                        %saveas(gcf, ['heatmap_sample' num2str(samplei) '.fig']);
                        %saveas(gcf, [resultfiletitle2 '_freerun' num2str(testi) '_heatmap.tif']);
                    end
                end
                
            end % end if size
        end % end for testi
    end % end for Ki
    
    figure
    plot(reshape(repmat(1:length(Ks),freerun_repeat,1), 1,[]), reshape(cormatch(:,samplei,:),[],1), '.')
    xticklabels(Ks)
    title(['sample' num2str(sampleID) ', effect of K' ])
    xlabel('K')
    ylabel('correation match score')
    xlim([0.5 length(Ks)+0.5])
    saveas(gcf,[resultfiletitle1 'hyouka_' 'sample' num2str(sampleID) '_K' '.tif']);
    saveas(gcf,[resultfiletitle1 'hyouka_' 'sample' num2str(sampleID) '_K' '.fig']);
    
end % samplei

save([resultfiletitle1 'cormatch.mat'],'cormatch')

allheat = squeeze(nanmean(cormatch, [1])); % freerun_repeat,samples,Ks
figure
a = heatmap(Ks, samples, allheat); %samples, Ks
a.XLabel = 'K';
a.YLabel = 'sample';
saveas(gcf, [resultfiletitle1 'heatdayo_allsamples.tif'])
saveas(gcf, [resultfiletitle1 'heatdayo_allsamples.fig'])

figure
for samplei = 1:length(samples)
    plot(reshape(repmat(1:length(Ks),freerun_repeat,1), 1,[]), reshape(cormatch(:,samplei,:),[],1), '.')
    hold on
end
xticklabels(Ks)
title(['Effects of K' ])
xlabel('K')
ylabel('correation match score')
xlim([0.5 length(Ks)+0.5])
legend(arrayfun(@(x) ['sample', num2str(x)], 1:length(samples), 'UniformOutput', false), 'Location','northeastoutside')

saveas(gcf,[resultfiletitle1 'hyouka_allsamples_K' '.tif']);
saveas(gcf,[resultfiletitle1 'hyouka_allsamples_K' '.fig']);

figure
for samplei = 1:length(samples)
    plot(1:length(Ks), squeeze(nanmean(cormatch(:,samplei,:),1)), '-o')
    hold on
end
xticklabels(Ks)
title(['Effects of K' ])
xlabel('K')
ylabel('correation match score')
xlim([0.5 length(Ks)+0.5])
legend(arrayfun(@(x) ['sample', num2str(x)], 1:length(samples), 'UniformOutput', false), 'Location','northeastoutside')

saveas(gcf,[resultfiletitle1 'hyouka_allsamples_mean_K' '.tif']);
saveas(gcf,[resultfiletitle1 'hyouka_allsamples_mean_K' '.fig']);

%{
cormatchmean = squeeze(nanmean(cormatch,1));
cormatchrank = NaN(length(samples),length(Ks));
for samplei = 1:length(samples)
[B,I] = sort(cormatchmean(samplei,:));
cormatchrank(samplei, I) = Ks;
end

figure
plot(mean(cormatchrank,1))
xlabel('K')
ylabel('correation mean rank')
xlim([0.5 length(Ks)+0.5])
%}

figure
cormatchmean = squeeze(nanmean(cormatch,1));
plot(mean(cormatchmean,1),'-o','LineWidth',3)
%errorbar(mean(cormatchmean,1),std(cormatchmean,1))
xlabel('K')
ylabel('correation mean allsamples')
title('Effect of K, averaged')
xlim([0.5 length(Ks)+0.5])
saveas(gcf,[resultfiletitle1 'hyouka_allmean_K' '.tif']);
saveas(gcf,[resultfiletitle1 'hyouka_allmean_K' '.fig']);


