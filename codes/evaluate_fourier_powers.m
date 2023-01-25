% evaluate_fourier_powers.m

% from ver98_Hbase/fourier_powers_summary_indirect_K345_kGMM2.m

samples = 1:24;
Ks = 3:5;
kGMM = 2;
freerun_repeat = 3;
link = 'indirect';
embed_width = 30;
embed_step = 10;

project_folder = '/home/iino/gKDR_GMM_50mM/Linked_gKDR_freerun';
gKDR_codes_folder = '/home/iino/gKDR_GMM_50mM/gKDR-CISE_ver3';
utilities_folder = '/home/iino/gKDR_GMM_50mM/Linked_gKDR_freerun';
metadata_folder = '/home/iino/gKDR_GMM_50mM/Linked_gKDR_freerun';
models_folder = '/home/iino/gKDR_GMM_50mM/Linked_gKDR_freerun/ver98_HLong';
common_cell_order_folder = '/home/iino/gKDR_GMM_50mM/Linked_gKDR_freerun/ver92_Hbase/common_cell_order';
simulation_results_folder = '/home/iino/gKDR_GMM_50mM/Linked_gKDR_freerun/ver98_HLong';
evaluation_result_folder = '/home/iino/gKDR_GMM_50mM/Linked_gKDR_freerun/publish_version/evaluation_results';
evaluation_result_subfolder = fullfile(evaluation_result_folder, ['results_evaluate_models_fourier_' link]);
        
%project_folder = '/home/iino/gKDR_GMM_50mM/Linked_gKDR_freerun/publish_version';
%gKDR_codes_folder = fullfile(project_folder, 'gKDR_codes');
%utilities_folder = fullfile(project_folder, 'utilities'); %'../'
%metadata_folder = fullfile(project_folder, 'metadata'); %'../'
%models_folder = fullfile(project_folder, 'models'); % ['../ver98_HLong/model_indirect_k' num2str(embed_width) '_tau' num2str(embed_step) '_kGMM' num2str(kGMM)];
%common_cell_order_folder = fullfile(project_folder, 'common_cell_order'); %'../ver92_Hbase/common_cell_order';
%simulation_results_folder = fullfile(project_folder, 'simulation_results'); %['../ver98_HLong/figures_' link '_k' num2str(embed_width) '_tau' num2str(embed_step) '_K' num2str(K) '_kGMM' num2str(kGMM) '/sample' num2str(sampleID)];
%evaluation_result_folder = fullfile(project_folder, 'evaluation_results');
%evaluation_result_subfolder = fullfile(evaluation_result_folder, ['results_evaluate_models_fourier_' link]);


if ~exist(evaluation_result_subfolder, 'dir')
    mkdir(evaluation_result_subfolder)
end

% read connection data
multiconmatrix = table2array(readtable(fullfile(metadata_folder, 'multiconmatrix.csv'),'ReadVariableNames',false));
conNames = table2cell(readtable(fullfile(metadata_folder, 'conneurons.csv'),'ReadVariableNames',false));
Mc = length(conNames);

for Ki = 1:length(Ks)
    K = Ks(Ki);
    disp(['<< K = ' num2str(K) ' >>'])
    
    resultfiletitle1 = [evaluation_result_subfolder '/fourier_power_indirect_K' num2str(K) '_kGMM' num2str(kGMM) '_allsamples'];
    
    targets0 = {'AVAL','AVAR'};
    
    close all
    tic
    
    %abs_array = NaN(length(targets0), length(samples), Mc);
    %angle_array = NaN(length(targets0), length(samples), Mc);
    rowNames = cell(1, length(samples));
    for i = 1:length(samples)
        rowNames{i} = ['sample' num2str(samples(i))];
    end
    variableNames = cell(1, length(targets0)*2);
    for i = 1:length(targets0)
        variableNames{2*i-1} = [targets0{i} ' real'];
        variableNames{2*i} = [targets0{i} ' predict'];
    end
    dat0 = NaN(length(samples), length(targets0)*2);
    absT = array2table(dat0, 'VariableNames', variableNames, 'RowNames', rowNames);
    angleT = array2table(dat0, 'VariableNames', variableNames, 'RowNames', rowNames);
    
    for samplei = 1:length(samples)
        sampleID = samples(samplei);
        
        disp(['< sample', num2str(sampleID),' >'])
        
        model_subfolder = fullfile(models_folder, ['model_indirect_k' num2str(embed_width) '_tau' num2str(embed_step) '_kGMM' num2str(kGMM)]);
        modelfileheader = fullfile(model_subfolder, ['sample' num2str(sampleID) '_K' num2str(K)]);
        modeldatafile = [modelfileheader '_modeldata.mat'];
        paramfile = [modelfileheader '_param.mat'];
        load(modeldatafile);
        load(paramfile);
        
        targetcellnames = uniqNames(targetcells);
        Mt = length(targetcellnames);
        m = Mt;
        n = size(data,1);
        
        testi=1;
        
        simulation_result_subfolder = fullfile(simulation_results_folder, ['figures_' link '_k' num2str(embed_width) '_tau' num2str(embed_step) '_K' num2str(K) '_kGMM' num2str(kGMM) '/sample' num2str(sampleID)]);
        simulationresultfileheader = fullfile(simulation_result_subfolder, ['sample' num2str(sampleID)]);
        load([simulationresultfileheader '_freerun' num2str(testi) '_savedata.mat']);
        
        freerun_start = min(6001,n);
        
        addpath('../');
        
        salttable = readtable(fullfile(metadata_folder, 'stimulation_timing.csv'),'ReadVariableNames',false,'HeaderLines',1); % from PreserveVariableNames
        startframe = table2array(salttable(sampleID,3));
        period = table2array(salttable(sampleID,4));
        
        
        assert(size(real_predict,2) == m+1)
        
        % FFT
        %real
        real_from = ceil(startframe/time_step);
        real_to = floor((freerun_start-1)/time_step);
        real_fft = fft(real_predict(real_from:real_to,:))/(real_to - real_from +1);
        real_peak = round((real_to-real_from+1)/(period/time_step*2))+1;
        real_salt_complex = real_fft(real_peak,:);
        
        real_salt_real = real(real_salt_complex);
        real_salt_imag = imag(real_salt_complex);
        real_salt_abs = abs(real_salt_complex);
        real_salt_angle = angle(real_salt_complex);
        
        %predict
        predict_from = ceil(freerun_start/time_step);
        predict_to = size(real_predict,1);
        predict_fft = fft(real_predict(predict_from:end,:))/(predict_to - predict_from +1);
        predict_peak = round((predict_to-predict_from+1)/(period/time_step*2))+1;
        predict_salt_complex = predict_fft(predict_peak,:);
        
        predict_salt_real = real(predict_salt_complex);
        predict_salt_imag = imag(predict_salt_complex);
        predict_salt_abs = abs(predict_salt_complex);
        predict_salt_angle = angle(predict_salt_complex);
        
        for targetno = 1:length(targets0) % AVAL, AVAR
            
            targetidx = find(strcmp(targetcellnames,targets0(targetno))); % without salt
            targetIdx = targetidx+1; % with salt
            if ~isempty(targetidx)
                %disp(targetIdx)
                absT(samplei, 2*targetno-1) = {real_salt_abs(targetIdx)};
                absT(samplei, 2*targetno) = {predict_salt_abs(targetIdx)};
                angleT(samplei, 2*targetno-1) = {real_salt_angle(targetIdx)/pi*180};
                angleT(samplei, 2*targetno) = {predict_salt_angle(targetIdx)/pi*180};
            end
        end
        
    end
    
    %save([resultfiletitle1 'abs_angle.mat'], 'absT', 'angleT');
    
    writetable(absT, [resultfiletitle1 '_absT.csv'],'WriteRowNames',true)
    writetable(angleT, [resultfiletitle1 '_angleT.csv'],'WriteRowNames',true)
    
    figure
    for i = 1:length(targets0)
        plot(table2array(absT(:,2*i-1)),'o')
        ax = gca;
        ax.ColorOrderIndex = i;
        hold on
        plot(table2array(absT(:,2*i)),'x')
        hold on
    end
    xlabel('samples')
    xlim([0, length(samples)+1])
    legend(variableNames)
    title('abs(stimulus frequency fourier component)')
    saveas(gcf, [resultfiletitle1 '_absT.tif'])
    saveas(gcf, [resultfiletitle1 '_absT.fig'])
    
    figure
    for i = 1:length(targets0)
        plot(table2array(angleT(:,2*i-1)),'o')
        ax = gca;
        ax.ColorOrderIndex = i;
        hold on
        plot(table2array(angleT(:,2*i)),'x')
        hold on
    end
    xlabel('samples')
    xlim([0, length(samples)+1])
    legend(variableNames)
    title('angle(stimulus frequency fourier component)')
    saveas(gcf, [resultfiletitle1 '_angleT.tif'])
    saveas(gcf, [resultfiletitle1 '_angleT.fig'])
    
end