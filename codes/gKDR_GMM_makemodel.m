
function gKDR_GMM_makemodel(UGE_TASK_ID_text)


% Requirements:
% [common_data_folder]/samplex_data.mat
% [metadata_folder]/conneurons.csv, multiconmatrix.csv
% where x stands for sample number.

% sampleID: sample number (command line input = UGE_TASK_ID_text)
% celli: cell number in uniqNames (therefore column number in data)
% cellc: cell number in conNames (connection table)
% cellt: cell number in targetcells
% note that connectome data are given based on conNames
% selc: cells within link from target cell, number in conNames
% seli: cells within link from target cell, number in uniqNames


sampleID = str2num(UGE_TASK_ID_text);

%%%%%%%%%%%%%%%%%%%%%  preset parameters %%%%%%%%%%%%%%%%%%%

Ks = 3:5;
kGMMs = [2];
link = 'indirect';

embed_width = 30; % number of embed_step's used for embedding = column number in source data
embed_step = 10; % invervals used for embedding (index-based)
time_step = 5; % timestep invervals in which estimation is performed
nahead = 5;    % how far ahead is the target for estimation

project_folder = '/home/iino/gKDR_GMM_50mM/Linked_gKDR_freerun/publish_version';
common_data_folder = fullfile(project_folder, 'common_data'); %'../ver92_Hbase/common_data';
gKDR_codes_folder = fullfile(project_folder, 'gKDR_codes');
utilities_folder = fullfile(project_folder, 'utilities'); %'../'
metadata_folder = fullfile(project_folder, 'metadata'); %'../'
%common_cell_order_folder = fullfile(project_folder, 'common_cell_order'); %'../ver92_Hbase/common_cell_order';
models_folder = fullfile(project_folder, 'models'); % ['../ver98_HLong/model_indirect_k' num2str(embed_width) '_tau' num2str(embed_step) '_kGMM' num2str(kGMM)];

show_plots = false;
use_salt_input = true;
saltsensors = {'ASEL','ASER','BAGL','BAGR','AWCL','AWCR','ASHL','ASHR'};

candx = 1; % candidates for CV, used for gKDR
candx2 = 2.5;  % candadate for CV, used for KRR
eps = [0.00001];  % epsilon for gKDR

prediction_method = 'GMM'; %'probabilistic_choice';


%%%%%%%%%%%%%%%%%%  Initial treatment  %%%%%%%%%%%%%%%%%%%

disp('start')
disp(['link type = ' link])
disp(['<< sample ' num2str(sampleID) ' >>']);
        
addpath(gKDR_codes_folder);
addpath(utilities_folder);

selistart = 1 - use_salt_input;
embed_size = embed_width;

% read connection data
multiconmatrix = table2array(readtable(fullfile(metadata_folder, 'multiconmatrix.csv'),'ReadVariableNames',false));
conNames = table2cell(readtable(fullfile(metadata_folder, 'conneurons.csv'),'ReadVariableNames',false));

tic


%%%%%  Hyperparameters loop  %%%%%

for Ki = 1:length(Ks)
    K0 = Ks(Ki);
    for kGMM = kGMMs
        
        disp(['< K= ' num2str(K0) ' >']);
        disp(['< kGMM ' num2str(kGMM) ' >']);
        
        model_subfolder = fullfile(models_folder, ['model_' link '_k' num2str(embed_width) '_tau' num2str(embed_step) '_kGMM' num2str(kGMM)]);
        if ~exist(model_subfolder,'dir')
            mkdir(model_subfolder)
        end
        modelfileheader = fullfile(model_subfolder, ['sample' num2str(sampleID) '_K' num2str(K0)]);
        
        if exist([modelfileheader '_modeldata.mat'],'file') && exist([modelfileheader '_param.mat'],'file')
            disp('modeldata.mat and param.mat already exist')
        else
            
            %%%%%%%%%%%%%%%   load sample   %%%%%%%%%%%%%%
            
            load(fullfile(common_data_folder,['sample' num2str(sampleID) '_data.mat']), 'data', 'uniqNames', 'targetcellnames', 'targetcells', 'Mt', 'autocorrthreshold', 'autocorrlag');
            
            cell_num = size(data,2);
            n = size(data,1);
            
            disp(['n is ' num2str(n)])
            
            train_span = [1, n];
            test_span = [];
            
            embed_before = embed_width-1;
            candy = candx; %1; %[0.25 0.5 0.75 1 2];
            
            if use_salt_input
                salttable = readtable(fullfile(metadata_folder, 'stimulation_timing.csv'),'ReadVariableNames',false,'HeaderLines',1); % from PreserveVariableNames
                startframe = table2array(salttable(sampleID,3));
                period = table2array(salttable(sampleID,4));
                %test = GMMpredict();
                %saltdata = generatesalt(freestart+freerun_length, startframe, period);   %%%%%%%%%%%%%%   Modified 20210515
                saltdata = generatesalt(n, startframe, period);
            end
            
            %%%%%%%%%%%%%   embed as a bulk  %%%%%%%%%%%%%%
            disp('execute embedding')
            if use_salt_input
                sourcedata = [saltdata(1:n)' data(:,targetcells)];
            else
                sourcedata = data(:,targetcells);
            end
            disp(size(sourcedata))
            [source_train_all, target_train_all, ~, ~] = ...
                embed4D(sourcedata, data(:,targetcells), train_span, test_span, embed_step, embed_before, time_step, nahead);
            disp('end of embedding')
            
            %%%%%%%%%%%   target cell loop  %%%%%%%%%%%%%%%%
            Bcell = cell(1,length(targetcells));
            Rcell = cell(1,length(targetcells));
            gmcell = cell(1,length(targetcells));
            selicell = cell(1,length(targetcells));
            colscell = cell(1,length(targetcells));
            target_train_cell = cell(1,length(targetcells));
            
            disp(['number of targetcells: ' num2str(Mt)]);
            
            for targeti = 1:Mt  %Mt = length(targetcells)
                targetcellname = targetcellnames{targeti};
                cellc = find(strcmp(conNames, targetcellname));
                if isempty(cellc)
                    disp(['There is no cell in connection table named ' targetcellname])
                else
                    % find connected neurons
                    if strcmp(link, 'direct')
                        selc = find(multiconmatrix(:,cellc) == 1); % in numbers in conNames
                    elseif strcmp(link, 'indirect')
                        selc = find(multiconmatrix(:,cellc) == 1);
                        selcall = selc;
                        for selci = selc'
                            if isempty(find(strcmp(targetcellnames, conNames{selci}),1)) % preがtargetcellになかったら１ステップ先を使う
                                selc2 = find(multiconmatrix(:,selci) == 1);
                                selcall = [selcall; selc2];
                            end
                        end
                        selc = unique(selcall);
                    elseif strcmp(link, 'all')
                        selc = (1:length(conNames))';
                    end
                    
                    selc(strcmp(conNames(selc),targetcellname)) = [];  %remove self(次、先頭にcellcを入れて、先頭が自分自身になるようにする)
                    selc = [cellc; selc]; % include self
                    connectedcellnames = conNames(selc);
                    
                    
                    seli = [];  % numbers in targetcells
                    for i = 1:length(connectedcellnames)
                        seli = [seli,find(strcmp(uniqNames(targetcells),connectedcellnames(i)))];%uniqNames(targetcells)は、targetcellnamesと同義。
                    end
                    sourcecellnames = targetcellnames(seli);
                    
                    if use_salt_input && any(strcmp(saltsensors, targetcellname))
                        seli = [seli(1) 0 seli(2:end)];
                        sourcecellnames = [sourcecellnames(1); {'salt'}; sourcecellnames(2:end)];
                    end
                    
                    disp(['sample ' num2str(sampleID) ' target ',targetcellname, ', No of linked neurons = ' num2str(length(selc)-1) ]);
                    disp(['number of presynaptic on connectome: ', num2str(sum(multiconmatrix(:,cellc) == 1)-1)]);
                    disp(['No of linked and annotated neurons = ' num2str(length(seli)-1) '( '  sprintf('%s ',targetcellnames{seli(seli>0)}) ')' ]); %uniqNames{targetcells(seli)}) ')']);
                    selicell{targeti} = seli; % record
                    
                    
                    target_train = target_train_all(:,targeti);
                    target_train_cell{targeti} = target_train;
                    
                    cols = reshape(repmat((1-selistart+seli-1)*embed_size, embed_size, 1)+repmat((1:embed_size)',1,length(seli)),[],1)';
                    colscell{targeti} = cols; % record
                    
                    source_train = source_train_all(:, cols);
                    
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%%%      execute KDR      %%%%%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    %change K depending on the size of cols
                    if length(cols) < K0
                        K = length(cols);
                    else
                        K = K0;
                    end
                    Kexe = min(K, size(source_train,2)); %%% added 20210812
                    sgx0=MedianDist(source_train);   % Basic value for bandwidth for X
                    sgy0=MedianDist(target_train);   % Basic value for bandwidth for Y
                    sgx = sgx0 * candx;
                    sgy = sgy0 * candy;
                    %sgx2 = sgx0 * candx2;
                    disp('Execute KerenelDeriv_chol')
                    errorflag = true;
                    for kdi = 1:3 % try three times until succeed
                        try
                            [B R t Kx]=KernelDeriv_chol(source_train, target_train, Kexe, sgx, sgy, eps);
                            errorflag = false;
                            break;
                        catch exception
                        end
                    end
                    
                    if errorflag
                        disp('errorflag=true')
                        save('error-data.mat', 'source_train', 'target_train','data','source_train_all','target_train_all')
                        msgText = getReport(exception);
                        disp(msgText)
                        return % quit this sample
                    end
                    
                    
                    if show_plots
                        shift = ((-embed_before):0)*embed_step;
                        embn = embed_before+1;
                        figure('Position', [100 50 600 300]);
                        imagesc(reshape(B,embed_size,[])')
                        yticks(1:length(targetcells));
                        yticklabels(repmat(sourcecellnames,Kexe,1));
                        xticks(1:embn);
                        xticklabels(shift);
                        title([uniqNames{targetcells(targeti)} ' before cise']);
                        ax = gca;
                        ax.YAxis.FontSize = 10;
                    end
                    
                    Bcell{targeti} = B;
                    Rcell{targeti} = R;
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%           make gmm                %%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                gm = [];
                for tempkGMM = kGMM:-1:1
                    try
                        gm = GMMfit(source_train, target_train, B, tempkGMM);
                        break
                    catch
                    end
                end
                if isempty(gm)
                    disp('GMMfit cannot pass');
                    stop
                end
                
                gmcell{targeti} = gm;
                
            end  %end of (for targeti = 1:Mt%Mt = length(targetcells))
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%           save results          %%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            save([modelfileheader '_modeldata.mat'],'uniqNames','targetcells','selicell','colscell','Bcell','gmcell','train_span','test_span','data','target_train_all','source_train_all');
            save([modelfileheader '_param.mat'], 'K','autocorrthreshold','autocorrlag','link','embed_width','embed_step', 'time_step', 'nahead', 'prediction_method', 'kGMM', 'candx', 'candx2', 'eps');
            disp(['Saved to ' modelfileheader '_modeldata.mat' ' and ' modelfileheader '_param.mat'])
        end
    end
end
end
