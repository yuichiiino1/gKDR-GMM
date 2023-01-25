%Perform freerun using the model generated by gKDR_GMM_makemodel.m

function gKDR_GMM_freerun(UGE_TASK_ID_text)
% UGE_TASK_ID_text

% Requires
% [model_folder]/samplex_model_data.mat
% [common_cell_order_folder]/samplex_common_outperm.mat



%%%%%%%%%%%%%%%%%%%%%  preset parameters %%%%%%%%%%%%%%%%%%%
Ks = 3:5;
kGMMs = [2];
link = 'indirect';
embed_width = 30; %10; % number of embed_step's used for embedding = column number in source data    % need
embed_step = 10; %10; % invervals used for embedding (index-based)                % need
freerun_length = 10000; % number of time_step's to perform freerun prediction
freerun_repeat = 3;
show_timeplot = false;
show_pca = false; %true;
perform_result_processing = true;

project_folder = '/home/iino/gKDR_GMM_50mM/Linked_gKDR_freerun/publish_version';
gKDR_codes_folder = fullfile(project_folder, 'gKDR_codes');
utilities_folder = fullfile(project_folder, 'utilities'); %'../'
metadata_folder = fullfile(project_folder, 'metadata'); %'../'
models_folder = fullfile(project_folder, 'models'); % ['../ver98_HLong/model_indirect_k' num2str(embed_width) '_tau' num2str(embed_step) '_kGMM' num2str(kGMM)];
common_cell_order_folder = fullfile(project_folder, 'common_cell_order'); %'../ver92_Hbase/common_cell_order';
simulation_results_folder = fullfile(project_folder, 'simulation_results'); %['../ver98_HLong/figures_' link '_k' num2str(embed_width) '_tau' num2str(embed_step) '_K' num2str(K) '_kGMM' num2str(kGMM) '/sample' num2str(sampleID)];


%%%%  Initial procedure  %%%%%

tic
rng('shuffle');
sampleID = str2num(UGE_TASK_ID_text); % input argument

use_salt_input = true;             % needed
selistart = 1 - use_salt_input;    % needed

%global embed_size
%embed_size = embed_width;
embed_before = embed_width-1;

addpath(gKDR_codes_folder);
addpath(utilities_folder);

if ~exist(simulation_results_folder, 'dir')
    mkdir(simulation_results_folder)
end

load(fullfile(common_cell_order_folder, ['sample' num2str(sampleID) '_common_outperm.mat']), 'outperm')



%%%%%  Different condition loop  %%%%%

for Ki = 1:length(Ks)
    K = Ks(Ki);
    for kGMM = kGMMs
        
        disp('start')
        disp(['sample ' num2str(sampleID)])
        disp(['< K ' num2str(K) ' >'])
        disp(['< kGMM ' num2str(kGMM) ' >']);

        simulation_result_subfolder = fullfile(simulation_results_folder, [link '_k' num2str(embed_width) '_tau' num2str(embed_step) '_K' num2str(K) '_kGMM' num2str(kGMM)]);
        
        if ~exist(simulation_result_subfolder, 'dir')
            mkdir(simulation_result_subfolder)
        end
        simulationresultfileheader = fullfile(simulation_result_subfolder, ['sample' num2str(sampleID)]);
        
        model_subfolder = fullfile(models_folder, ['model_' link '_k' num2str(embed_width) '_tau' num2str(embed_step) '_kGMM' num2str(kGMM)]);
        modelfileheader = fullfile(model_subfolder, ['sample' num2str(sampleID) '_K' num2str(K)]);
        modeldatafile = [modelfileheader '_modeldata.mat'];
        paramfile = [modelfileheader '_param.mat'];
        load(modeldatafile, 'uniqNames','targetcells','selicell','colscell','Bcell','gmcell','train_span','data','target_train_all','source_train_all');
        load(paramfile, 'time_step','nahead', 'prediction_method');
        
        disp([num2str(length(uniqNames)) ' names, ' num2str(size(data,2)) ' data'])
        
        n = size(data, 1); % data, uniqNames and train_span in included in modeldata
        disp(['n=',num2str(n)])
        
        freerun_start = train_span(2)+1;
        freestart = freerun_start;
        if freestart > n+1
            freestart = n+1;
        end
        
        if use_salt_input
            
            salttable = readtable(fullfile(metadata_folder, 'stimulation_timing.csv'),'ReadVariableNames',false,'HeaderLines',1); % from PreserveVariableNames
            startframe = table2array(salttable(sampleID,3));
            period = table2array(salttable(sampleID,4));
            saltdata = generatesalt(freestart+freerun_length, startframe, period);   %%%%%%%%%%%%%%   Modified 20210515
            
            %saltdata = zeros(1,freestart+freerun_length);
        end
        
        % targetcells is included in modeldata
        %Mt = length(targetcells);
        targetcellnames = uniqNames(targetcells);
        
        %%%%%%%%%  embedding  %%%%%
        disp('execute embedding')
        if use_salt_input
            sourcedata = [saltdata(1:n)' data(:,targetcells)];
        else
            sourcedata = data(:,targetcells);
        end
        %train_span = [max(1,train_span(1)) min(n,train_span(2))];
        test_span = [];
        
        [source_train_all, target_train_all, ~, ~] = ...
            embed4D(sourcedata, data(:,targetcells), train_span, test_span, embed_step, embed_before, time_step, nahead);
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%   load or make gm   %%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % selicell, gmcell has been loaded from model.data
        
        Tc_real = 1:time_step:(freestart-1); % sequence of time indices for real data before freerun
        Tc_all = 1:time_step:(freestart+freerun_length);
        nx = length(Tc_all);    % Tc_all from saveddata  XXXXXXXXXXX
        nc = length(Tc_real);    % Tc_real from saveddata  XXXXXXXXXXX
        
        
        %%%%%%%%%%%%%%%%%%%%%%
        %%%%%  Freerun  %%%%%%
        %%%%%%%%%%%%%%%%%%%%%%
        
        Tc_real = train_span(1):time_step:(freestart-1); % sequence of time indices for real data before freerun
        Tc_all = train_span(1):time_step:(freestart+freerun_length);
        
        assert(mod(embed_step,time_step)==0 && mod(nahead,time_step)==0)
        disp('Freerun prediction')
        embed_time_step = embed_step/time_step; % embed_stepがtime_step単位でいくつか
        ahead_time_step = nahead/time_step; % naheadがtime_step単位でいくつか
        
        %%%% new real_perdict_set %%%%
        real_predict_set = cell(1,freerun_repeat);
        
        for testi = 1:freerun_repeat
            
            disp(['repeat ' num2str(testi)])
            real_predict = NaN(length(Tc_all),length(targetcells)+1-selistart); % add if salt input
            
            %real_mat = NaN(size(target_train_all,1),length(targetcells)+1-selistart); % trainの間の実データ
            
            % put real data at the begining of real_predict
            real_predict(1:length(Tc_real), 1-selistart+(1:length(targetcells))) = data(Tc_real, targetcells);
            if use_salt_input
                real_predict(1:length(Tc_all),1) = saltdata(Tc_all); % salt input      %%% saltdata from modeldata
            end
            
            % freerun prediction
            for ti = 1:(nx-nc) % every time_step
                if mod(ti,100)==0
                    disp(num2str(ti))
                end
                
                for targeti = 1:length(targetcells)
                    source_for_predict = ( reshape(real_predict((nc+ti-ahead_time_step-embed_time_step*embed_before):embed_time_step:(nc+ti-ahead_time_step), selicell{targeti}+1-selistart),[],1) )';  %nc+tiをpredictしたいので、そこからnaheadもどったところがembedの終点になる
                    
                    if strcmp(prediction_method, 'GMM') % other methods are not supported in this code
                        B = Bcell{targeti};
                        %predict_ahead = GMMpredict_det(source_for_predict, gmcell{targeti}, B);
                        predict_ahead = GMMpredict(source_for_predict, gmcell{targeti}, B);
                    end
                    real_predict(nc+ti, targeti+1-selistart) = predict_ahead;
                end
                
            end
            
            
            if show_timeplot
                for targeti = selistart:length(targetcells)
                    if mod((targeti+1-selistart-1),10) == 0
                        figure('Position', [300 50 900 650]);
                    end
                    subplot(5,2,mod(targeti+1-selistart-1,10)+1)
                    plot(Tc_all, real_predict(:,targeti+1-selistart));
                    hold on
                    x0 = Tc_real(end)+time_step/2;
                    plot([x0,x0],ylim,'color','blue');
                    if targeti == 0
                        title('salt input')
                    else
                        title(targetcellnames{targeti})
                    end
                    if(mod((targeti-selistart+1),10) == 0 || targeti==length(targetcells))
                        saveas(gcf, [simulationresultfileheader '_freerun' num2str(testi) '_plot_' num2str(ceil((targeti-selistart)/10)) '.tif']);
                        saveas(gcf, [simulationresultfileheader '_freerun' num2str(testi) '_plot_' num2str(ceil((targeti-selistart)/10)) '.fig']);
                    end
                end
            end
            
            real_predict_set{testi} = real_predict;
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%                        %%%%%%%%
        %%%%%%%%%%%%%   result_processing    %%%%%%%%
        %%%%%%%%%%%%%                        %%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % modification of clustering_and_heatmap.m
        
        if perform_result_processing
            
            for testi = 1:freerun_repeat
                
                real_predict = real_predict_set{testi};
                
                if size(real_predict, 2) <= 1
                    disp(['number of cells for sample' num2str(sampleID) ' is too small (' num2str(size(real_predict, 2)) ')'])
                else
                    use_salt_input = size(real_predict,2)-length(selicell);  % 0 or 1
                    
                    data2 = real_predict(1:nc, (use_salt_input+1):end); % Real part only, remove 'salt'
                    cell_num = size(data2,2);
                    data_pred = real_predict((nc+1):nx, (use_salt_input+1):end);
                    
                    if ~isempty(outperm)
                        
                        lab=cellstr(targetcellnames(outperm));
                        labnum = cell(1,cell_num);
                        for celli =1:cell_num
                            labnum{celli} = [strrep(lab{celli}, '_', '\_')]; %  ' ' num2str(celli)
                        end
                        
                        %%%%%%%%%%%%%%%%   correlation matrix plot  %%%%%%%%%%%%%%%%
                        
                        % for real data
                        figure;
                        rmatrix_real = corrcoef(data2(:,outperm),'rows','complete'); % omit NaN
                        cormatrix_plot(rmatrix_real, ['sample' num2str(sampleID) ' real'], 'real', simulationresultfileheader, cell_num, labnum);%,sampleID, K, lambda);
                        
                        % for predicted data
                        figure;
                        rmatrix = corrcoef(data_pred(:,outperm),'rows','complete'); % omit NaN
                        cormatrix_plot(rmatrix, ['sample' num2str(sampleID) ' freerun' num2str(testi)], ['freerun' num2str(testi)],  simulationresultfileheader, cell_num,labnum);%,sampleID, K, lambda);
                        
                        %%%%%%%%%%%%%%%%%   heatmap  %%%%%%%%%%%%%%%%
                        
                        
                        % real data heamap
                        
                        data_real = [saltdata(Tc_real)' data2];
                        
                        figure;
                        colormap jet(64);
                        image(data_real(:,[1 outperm+1])'*10+32);%
                        title(['sample' num2str(sampleID)],'Interpreter' ,'none');
                        xlabel('Time points')
                        colorbar;
                        h = gca;
                        set(h,'ytick',(1:(1+cell_num))');
                        set(h,'yticklabel',nameshifter([{'salt'},labnum]));   % if use salt data
                        %set(h,'yticklabel',[{'salt'},labnum]);
                        set(h,'fontsize',8);
                        set(h,'Ticklength',[0 0]);
                        set(h,'xtick',(1:500:nc)');
                        set(h,'xticklabel',Tc_real(1:500:end)-1);
                        hold on
                        saveas(gcf, [simulationresultfileheader '_real' '_heatmap.tif']);
                        saveas(gcf, [simulationresultfileheader '_real' '_heatmap.fig']);
                        
                        % simulation results heatmap
                        
                        data_all = [saltdata(Tc_all)' [data2; data_pred]];
                        figure;
                        colormap jet(64);
                        image(data_all(:,[1 outperm+1])'*10+32);%
                        title(['sample' num2str(sampleID)],'Interpreter' ,'none');
                        xlabel('Time points')
                        colorbar;
                        h = gca;
                        set(h,'ytick',(1:(1+cell_num))');
                        set(h,'yticklabel',nameshifter([{'salt'},labnum]));   % if use salt data
                        %set(h,'yticklabel',[{'salt'},labnum]);
                        set(h,'fontsize',8);
                        set(h,'Ticklength',[0 0]);
                        set(h,'xtick',(1:500:nx)');
                        set(h,'xticklabel',Tc_all(1:500:end)-1);
                        hold on
                        plot([nc+0.5 nc+0.5], [0 cell_num], '-', 'color',[0.5 0.5 0.5]);
                        saveas(gcf, [simulationresultfileheader '_freerun' num2str(testi) '_heatmap.tif']);
                        saveas(gcf, [simulationresultfileheader '_freerun' num2str(testi) '_heatmap.fig']);
                        
                    end
                    
                    
                    %%%%%%%%%%%%%%%   PCA   %%%%%%%%%%
                    if show_pca
                        [coeff,score,latent,tsquare] = pca(real_predict(1:nc,:));
                        figure;
                        plot3(score(:,1),score(:,2),score(:,3));
                        rotate3d on
                    end
                    
                end
                save([simulationresultfileheader '_freerun' num2str(testi) '_savedata.mat'], 'uniqNames','targetcells','colscell','sampleID','targetcellnames','Bcell','gmcell','saltdata', 'selicell','Tc_real','Tc_all','real_predict', 'outperm', 'rmatrix','train_span','target_train_all','source_train_all');
            end
        end
        toc
    end
end
end