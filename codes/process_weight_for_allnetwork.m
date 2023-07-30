function process_weight_for_allnetwork()

%sampleID = str2num(UGE_TASK_ID_text);
samples = 1:24;
FDRthre = 0.005; %0.05;
alpha = 0.001; %0.005; %0.05;
show_each_pair_plot = false;
only_direct = true; %false;
omit_pharyngeal = false; %true;
%vartype = 'std'; % 'coefvar'
partition_number = 3;
offsets = 0:4;
link = 'indirect';
Ks = 3:5; %[4];
kGMMs = 2;
validmethod = 'max'; %'min'; %'', 'max'   % actually max is used in this version
embed_width = 30; %10; % number of embed_step's used for embedding = column number in source data    % need
embed_step = 10; %10; % invervals used for embedding (index-based)                % need

project_folder = '/home/iino/gKDR_GMM_50mM/Linked_gKDR_freerun/publish_version';
%project_folder = '/Users/iino/Documents/___Transcend J/_4DImaging_Data_Analysis/gKDR-GMM/gKDR_GMM_50mM/Linked_gKDR_freerun/publish_version/';
common_data_folder = fullfile(project_folder, 'common_data'); %'../ver92_Hbase/common_data';
gKDR_codes_folder = fullfile(project_folder, 'gKDR_codes');
utilities_folder = fullfile(project_folder, 'utilities'); %'../'
metadata_folder = fullfile(project_folder, 'metadata'); %'../'
common_cell_order_folder = fullfile(project_folder, 'common_cell_order');
models_folder = fullfile(project_folder, 'models'); % ['../ver98_HLong/model_indirect_k' num2str(embed_width) '_tau' num2str(embed_step) '_kGMM' num2str(kGMM)];
evaluateL_results_folder = fullfile(project_folder, 'evaluateL_results'); %['../ver98_HLong/figures_' link '_k' num2str(embed_width) '_tau' num2str(embed_step) '_K' num2str(K) '_kGMM' num2str(kGMM) '/sample' num2str(sampleID)];

addpath(gKDR_codes_folder);
addpath(utilities_folder);

%{
conneuron_classes = { ...
    'CEP',{'CEPDL','CEPVL','CEPVR','CEPDR'}; ...
    'IL1',{'IL1L','IL1DL','IL1VL','IL1VR','IL1DR','IL1R'}; ...
    'IL2',{'IL2L','IL2DL','IL2VL','IL2VR','IL2DR','IL2R'}; ...
    'OLQ',{'OLQDL','OLQVL','OLQVR','OLQDR'}; ...
    'RMD',{'RMDL','RMDDL','RMDVL','RMDVR','RMDDR','RMDR'}; ...
    'RME',{'RMEL','RMED','RMEV','RMER'}; ...
    'SAA',{'SAADL','SAAVL','SAAVR','SAADR'}; ...
    'SAB',{'SABVL','SABD','SABVR'}; ...
    'SIA',{'SIADL','SIAVL','SIAVR','SIADR'}; ...
    'SIB',{'SIBDL','SIBVL','SIBVR','SIBDR'}; ...
    'SMB',{'SMBDL','SMBVL','SMBVR','SMBDR'}; ...
    'SMD',{'SMDDL','SMDVL','SMDVR','SMDDR'}; ...
    'URA',{'URADL','URAVL','URAVR','URADR'}; ...
    'URY',{'URYDL','URYVL','URYVR','URYDR'}; ...
    };

groupA_names = {'RIMR','AVAR','AVAL','AVEL','RIML','AIBL','AIBR','AVER','RMDL','RMDVR'}; %above 13 point, remove ASHL. %{'AVAL','AVAR','RIML','RIMR','AVEL','RIMR'};
groupB_names = {'RMEL','RMER','RID','RMED','RIBR','RIBL','SMBVR'}; %,'AVJL'}; %below -5 points. %{'RIBL','RIBR','RMEL','RMER','RID','AVBL','AVBR'};
saltsensors = {'ASEL','ASER','BAGL','BAGR','AWCL','AWCR','ASHL','ASHR'};
amphid = {'AWAL','AWAR','AWBL','AWBR','AFDL','AFDR','AFDL','AFDR','ASGL','ASGR','ASIL','ASIR','ASJL','ASJR','ASKL','ASKR','ADLL','ADLR'};
ciliated = {'IL1L','IL1R','IL1VL','IL1VR','IL1DL','IL1DR',  'IL2L','IL2R','IL2VL','IL2VR','IL2DL','IL2DR',  'OLLL','OLLR',  'OLQDL','OLQDR',  'OLQVL','OLQVR',...
    'CEPDL','CEPDR',  'CEPVL','CEPVR',  'AQR','FLPL','FLPR'};
pharyngeal = {'I1L','I1R','I2L','I2R','I3','I4','I5','I6','M1','M2L','M2R','M3L','M3R','M4','M5','MCL','MCR','MI','NSML','NSMR'};
saltsensorclass = {'BAG','ASE','AWC','ASH'};
amphidclass = {'AWA','AWB','AFD','AFD','ASG','ASI','ASJ','ASK','ADL'};
ciliatedclass = {'IL1',  'IL2', 'OLL',  'OLQ', 'CEP', 'AQR','FLP'};
primaryinterclass = {'AIA','AIB','AIY','AIZ'};

%}

if only_direct
    dirtext = '_direct_link';
else
    dirtext = '_all_link';
end

expsavefolder = fullfile(project_folder, 'connection_strength', [link '_cross' num2str(partition_number) validmethod dirtext '_FDR' num2str(FDRthre)]);
if ~exist(expsavefolder, 'dir')
    mkdir(expsavefolder)
end

globalNames = table2cell(readtable(fullfile(metadata_folder,'globalNames.csv'),'ReadVariableNames',false));
Mg = length(globalNames);
disp(['Mg=',num2str(Mg)])
globalNames1 = [{'salt'}; globalNames];

multiconmatrix = table2array(readtable(fullfile(metadata_folder, 'multiconmatrix.csv'),'ReadVariableNames',false));
conNames = table2cell(readtable(fullfile(metadata_folder, 'conneurons.csv'),'ReadVariableNames',false));
conNames1 = [{'salt'}; conNames];
Mc = length(conNames);

%validcon = false(length(samples),Mc+1);
validglobal1 = false(length(samples),Mg+1);
globalweightarray = NaN(Mg+1,Mg,partition_number,length(offsets),length(Ks),length(samples));
%globalminvalidarray = NaN(Mg, 3, length(samples));
%globalmaxvalidarray = NaN(Mg, 3, length(samples));
%weightarray_set = cell(1,length(samples));
globalpmat = NaN(Mg+1, Mg, length(samples)); % p values of source-target combinations for all samples.
globalsigweightnonself = NaN(Mg+1, Mg, length(samples)); % mean weight for only significant.
globalvalidweightarray = NaN(Mg+1,Mg,partition_number,length(offsets),length(Ks),length(samples));

for samplei = 1:length(samples)
    sampleID = samples(samplei);
    disp(['< sample ' num2str(sampleID) ' >']);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%       Pretreat each  sample         %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % load basic data
    load(fullfile(common_data_folder,['sample' num2str(sampleID) '_data.mat']), 'data', 'uniqNames', 'targetcellnames', 'targetcells', 'Mt');

    targetcellnames = uniqNames(targetcells);
    targetcellnames1 = [{'salt'}; targetcellnames];
    Mt = length(targetcells);

    for kGMM = kGMMs
        disp(['< kGMM ' num2str(kGMM) ' >']);

        weightarray = NaN(Mt+1,Mt,partition_number,length(offsets),length(Ks));
        p_array_min = NaN(Mt,partition_number,length(offsets),length(Ks));
        p_array_min_valid = NaN(Mt,partition_number,length(offsets),length(Ks));
        p_array_max = NaN(Mt,partition_number,length(offsets),length(Ks));
        p_array_max_valid = NaN(Mt,partition_number,length(offsets),length(Ks));

        for Ki = 1:length(Ks)
            K = Ks(Ki);

            disp(['< K= ' num2str(K) ' >']);

            for offseti = 1:length(offsets)
                offset = offsets(offseti);

                for partitioni = 1:partition_number
                    %disp(['< part ' num2str(partitioni) ' >']);

                    % load weightmatrix  and combine to weightarrray (Mt+1,Mt,3)
                    % (weightmatrix is output of evaluate_connection_strength_cross.m)
                    weightfileheader = fullfile(project_folder, 'connection_strength', [link '_cross' num2str(partition_number)], ['sample' num2str(sampleID)  '_K' num2str(K) '_kGMM' num2str(kGMM) '_offset' num2str(offset) '_part' num2str(partitioni)]);
                    weightmatrix = NaN(Mt+1,Mt);
                    try
                        load([weightfileheader '_weightmatrix.mat'], 'weightmatrix','weightmatrix1')
                        disp(['offset'  num2str(offset) ' part' num2str(partitioni) ' loaded'])
                    catch
                        disp(['offset'  num2str(offset) ' part' num2str(partitioni) ' could not load'])
                    end
                    weightmatrix(1,:) = - weightmatrix(1,:);  %%%%%%%%%%%%%%%%%  flip salt  %%%%%%%%%%%%%
                    weightarray(:,:,partitioni,offseti,Ki) = weightmatrix;   % weightarray: Mt+1,Mt,partition_number,length(offsets),length(Ks)

                end
            end


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%  determine validity based on p_array  %%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % load p_array  (output of gKDR_GMM_cross_evaluateL_new.m)
            evaluateL_result_subfolder = fullfile(evaluateL_results_folder, [link '_k' num2str(embed_width) '_tau' num2str(embed_step) '_K' num2str(K) '_kGMM' num2str(kGMM) '_new_offsets']);
            evaluateL_result_file = fullfile(evaluateL_result_subfolder, ['sample' num2str(sampleID) '.mat']);
            try
                load(evaluateL_result_file, 'p_array'); %,'targetcellnames','sumlogL0_array','permute_salt')
                disp(['p_array', ' loaded'])
                p_array_n = p_array; % Mt,partition_number,partition_number,offset_number
                for offseti = 1:length(offsets)
                    for i = 1:partition_number
                        p_array_n(:,i,i,offseti) = NaN;
                    end
                end

                % min
                p_array_min(:,:,:,Ki) = squeeze(nanmin(p_array_n, [], 3));   % (Mt, 3, 5)
                p_array_min_valid(:,:,:,Ki) = p_array_min(:,:,:,Ki) < alpha; %(Mt+1,partition_number,length(offsets),length(Ks));

                % max
                p_array_max(:,:,:,Ki) = squeeze(nanmax(p_array_n, [], 3));   % (Mt, 3, 5)
                p_array_max_valid(:,:,:,Ki) = p_array_max(:,:,:,Ki) < alpha; %(Mt+1,partition_number,length(offsets),length(Ks));

            catch
                disp(['p_array', ' cannot load'])
                % min
                p_array_min(:,:,:,Ki) = 1;
                p_array_min_valid(:,:,:,Ki) = 1;
                % max
                p_array_max(:,:,:,Ki) = 1;
                p_array_max_valid(:,:,:,Ki) = 1;

            end

        end % end Ki

        %%%%% register to global array %%%%

        cellg = [];   %indexes of targetcells in globalNames
        for targeti = 1:Mt
            cellname = targetcellnames{targeti};
            celli = find(strcmp(globalNames, cellname));
            assert(length(celli)==1)
            cellg = [cellg, celli];
        end
        cellg1 = [1, cellg+1]; % first is salt
        globalweightarray(cellg1, cellg, :,:,:, samplei) = weightarray;

        %globalminvalidarray(cellg,:,samplei) = p_array_min_valid;
        %globalmaxvalidarray(cellg,:,samplei) = p_array_max_valid;
        %weightarray_set{samplei} = weightarray;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%  evaluate validity for each sample %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        mincount = 5;
        p_array_max_valid_temp = p_array_max_valid;
        for targeti = 1:Mt
            disp([targetcellnames{targeti} ': valid number ' num2str(nansum(reshape(p_array_max_valid(targeti,:,:,:),1,[])))])
            if nansum(reshape(p_array_max_valid(targeti,:,:,:),1,[])) < mincount
                p_array_max_valid_temp(targeti,:,:,:) = false;
            end
        end

        p_array_max_valid_ex = repmat(reshape(p_array_max_valid_temp, [1,Mt,partition_number,length(offsets),length(Ks)]), [Mt+1,1,1,1,1]);
        validweightarray = weightarray; % weightarray: Mt+1,Mt,partition_number,length(offsets),length(Ks)
        validweightarray(~p_array_max_valid_ex) = NaN;  % validweightarray: Mt+1,Mt,partition_number,length(offsets),length(Ks)
        globalvalidweightarray(cellg1, cellg, :,:,:, samplei) = validweightarray;

        sigweights = NaN(Mt+1,Mt);
        %sigweightset = NaN(Mg+1,Mg,3,length(samples));
        pmat = NaN(Mt+1,Mt);
        for row = 1:Mt+1
            for col = 1:Mt
                weights = reshape(validweightarray(row,col, :, :),1,[]);   % (3, length(samples))
                if all(isnan(weights))
                    p=1;
                else
                    %[h,p] = ttest(weights(~isnan(weights)));
                    p = ranksum(weights, -weights);
                end
                pmat(row,col) = p;
                sigweights(row,col) = nanmean(weights);
            end
        end
        globalpmat(cellg1, cellg, samplei) = pmat;

        FDRmat = reshape(mafdr(reshape(pmat,[],1), 'BHFDR', true),Mt+1,Mt);
        sigweights(FDRmat > FDRthre) = NaN;
        sigweightnonself = sigweights;
        for i = 1:Mt
            sigweightnonself(i+1,i) = NaN;
        end

        
        % connectome

        ct = NaN(1,Mt);  % targetcellnames -> conNames
        for i = 1:Mt
            ct(i) = find(strcmp(conNames, targetcellnames{i}));
        end

        for row = 1:Mt+1
            for col = 1:Mt
                connected = false;
                if row-1 ~= col
                    if row==1
                        connected = true;
                    elseif multiconmatrix(ct(row-1),ct(col))==1
                        connected = true;
                    end
                end
                if ~connected && strcmp(dirtext, '_direct_link')
                    sigweightnonself(row,col) = NaN; %false;
                end
            end
        end
        ct1 = [1, ct+1];
        globalsigweightnonself(cellg1, cellg, samplei) = sigweightnonself;

        validtargetbool = any(sigweightnonself,1);
        validsourcebool = any(sigweightnonself,2)';
        validtarget = 1:Mt;
        validtarget = validtarget(validtargetbool);
        validsource = 1:Mt+1;
        validsource = validsource(validsourcebool);
        validglobal1(samplei, unique([cellg(validtarget)+1 cellg1(validsource)])) = true;

    end % kGMM
end % sample

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  evaluate validity for all samples %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

globalvalidweightarraymean = squeeze(nanmean(globalvalidweightarray, 3:5)); % globalvalidweightarray: Mg+1,Mg,partition_number,length(offsets),length(Ks), samples
globalpmat = NaN(Mg+1,Mg);
for row = 1:Mg+1
    for col = 1:Mg
        %weights = globalsigweightnonself(row, col, :); % across samples
        weights = squeeze(globalvalidweightarraymean(row,col,:)); % across samples
        if all(isnan(weights))
            p=1;
        else
            %[h,p] = ttest(weights(~isnan(weights)));
            p = ranksum(weights, -weights);
        end
        globalpmat(row,col) = p;
    end
end

globalFDRmat = reshape(mafdr(reshape(globalpmat,[],1), 'BHFDR', true),Mg+1,Mg);
globalsigweights = nanmean(globalsigweightnonself, 3); % Mg+1,Mg,samples
globalsigweights(globalFDRmat > FDRthre) = NaN;



%%%%%%%%   valid table  %%%%%%%%%

%usedcon = any(validcon,1); % validcon: length(samples),Mc+1
%validcontable = array2table(validcon(:,usedcon), 'VariableNames', conNames1(usedcon),'RowNames', arrayfun(@(x){['sample' num2str(x)]}, samples));
%writetable(validcontable, fullfile(evaluateL_results_folder, 'validcontable.csv'))
usedglobal = any(validglobal1,1);
used_global_light = validglobal1(:,usedglobal);
validglobaltable = array2table(used_global_light, 'VariableNames', globalNames1(usedglobal),'RowNames', arrayfun(@(x){['sample' num2str(x)]}, samples));
writetable(validglobaltable, fullfile(expsavefolder, 'validglobaltable.csv'), 'WriteRowNames', true)


figure()
%imagesc(validcon(:,usedcon))
%xticks(1:sum(usedcon))
%xticklabels(conNames1(usedcon))
imagesc(used_global_light)
xticks(1:sum(usedglobal))
xticklabels(globalNames1(usedglobal))
yticks(1:length(samples))
yticklabels(arrayfun(@(x){['sample' num2str(x)]}, samples))
%saveas(gcf, fullfile(evaluateL_results_folder, 'validcontable.fig'))
%saveas(gcf, fullfile(evaluateL_results_folder, 'validcontable.tif'))
saveas(gcf, fullfile(expsavefolder, 'validglobaltable.fig'))
saveas(gcf, fullfile(expsavefolder, 'validglobaltable.tif'))

core_neuron_threshold = 1;  % if used in more than n samples, considered core in the graph

globalvalidmat = globalFDRmat <= FDRthre;
globalused = [false, any(globalvalidmat,1)] | any(globalvalidmat,2)';
core_neurons_bool = squeeze(nansum(validglobal1, 1)) >= core_neuron_threshold | globalused;
disp(globalNames1(core_neurons_bool)')

save(fullfile(expsavefolder, 'process_weight_allnetwork_core_phys_offsets.mat'))

end

