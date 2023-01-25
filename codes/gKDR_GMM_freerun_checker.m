%Examines completion of freerun by gKDR_GMM_freerun.m

function gKDR_GMM_freerun_checker()

%%%%%%%%%%%%%%%%%%%%%  preset parameters %%%%%%%%%%%%%%%%%%%
samples = 1:24;
Ks = 3:5;
kGMM = 2;
link = 'indirect';
embed_width = 30; %10; % number of embed_step's used for embedding = column number in source data
embed_step = 10; %10; % invervals used for embedding (index-based)
freerun_repeat = 5;

project_folder = '/home/iino/gKDR_GMM_50mM/Linked_gKDR_freerun/publish_version';
simulation_results_folder = fullfile(project_folder, 'simulation_results');

%%%%%%%%%%%%%%%%%%%%%  perform checking %%%%%%%%%%%%%%%%%%%
for sampleID = samples
    for Ki = 1:length(Ks)
        K = Ks(Ki);
        simulation_result_subfolder = fullfile(simulation_results_folder, [link '_k' num2str(embed_width) '_tau' num2str(embed_step) '_K' num2str(K) '_kGMM' num2str(kGMM)]);
        simulationresultfileheader = fullfile(simulation_result_subfolder, ['sample' num2str(sampleID)]);
        for testi = 1:freerun_repeat
            try
                load([simulationresultfileheader '_freerun' num2str(testi) '_savedata.mat']);
                disp(['sample' num2str(sampleID) ' K' num2str(K) ' repeat' num2str(testi) ':  done']);
            catch
                disp(['sample' num2str(sampleID) ' K' num2str(K) ' repeat' num2str(testi) ':  MISSING or BROKEN!']);
            end
        end
    end
end


