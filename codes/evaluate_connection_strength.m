
function evaluate_connection_strength(UGE_TASK_ID_text) % sampleID
% originally ver102_gradient/evaluate_connection_strength5_flip.m
% gradient at 17:30 was used and simply averaged.


%%%%%%%%%%%%%%%%%%%%%  preset parameters %%%%%%%%%%%%%%%%%%%

sampleID = str2num(UGE_TASK_ID_text);
Ks = [3,4,5];
kGMM = 2;
link = 'indirect';
embed_width = 30; % number of embed_step's used for embedding = column number in source data
embed_step = 10; % invervals used for embedding (index-based)
plot_timeweights = false;
plot_matrix = true;

project_folder = '/home/iino/gKDR_GMM_50mM/Linked_gKDR_freerun/publish_version';
utilities_folder = fullfile(project_folder, 'utilities'); %'../'
common_cell_order_folder = fullfile(project_folder, 'common_cell_order'); %'../ver92_Hbase/common_cell_order';
models_folder = fullfile(project_folder, 'models');
%results_folder = fullfile(project_folder, 'connection_strength');
results_folder = fullfile(project_folder, ['connection_strength/' link]);

%%%%  Initial procedure  %%%%%

addpath(utilities_folder);

if ~exist(results_folder ,'dir')
    mkdir(results_folder)
end


for Ki = 1:length(Ks)
    K = Ks(Ki);
    
    disp(['sample ' num2str(sampleID) ' K ' num2str(K)])

    
    model_subfolder = fullfile(models_folder, ['model_' link '_k' num2str(embed_width) '_tau' num2str(embed_step) '_kGMM' num2str(kGMM)]);
    modelfileheader = fullfile(model_subfolder, ['sample' num2str(sampleID) '_K' num2str(K)]);
    load([modelfileheader '_gradient.mat'], 'gm_grad_mean','targetcells','uniqNames'); %,'gm_grad_cell', 'gm_grad_array','gm_grad_mean_cell','selicell')
    
    targetcellnames = uniqNames(targetcells);
    targetcellnames1 = [{'salt'}; targetcellnames];
    Mt = length(targetcells);
    %disp('number of presynaptic cells:')
    %disp(arrayfun(@(x) length(selicell{x}), 1:Mt))
    if plot_timeweights
        for targeti = 1:Mt
            disp(['target', num2str(targeti) ': ' targetcellnames{targeti} ' ' num2str(length(seli)) ' pre']);
            seli = selicell{targeti};
            %disp(length(seli))
            figure('Position',[10,10,1000,600]);
            count = 1;
            for i = 1:Mt
                if ~isnan(gm_grad_mean(1,i+1,targeti))
                    subplot(3,4,count)
                    plot(squeeze(gm_grad_mean(:,i+1,targeti)))
                    hold on
                    plot([1 5], [0 0], 'k')
                    title(['cell ' num2str(i) '(' targetcellnames{i} ') to ' targetcellnames{targeti}])
                    ylim([-0.3,0.3])
                    if count==12
                        break
                    end
                    count = count+1;
                end
            end
            saveas(gcf, [results_folder '/timeweight_K' num2str(K) '_sample' num2str(sampleID) '_cell' num2str(targeti) '_' link '.tif'])
            saveas(gcf, [results_folder '/timeweight_K' num2str(K) '_sample' num2str(sampleID) '_cell' num2str(targeti) '_' link '.fig'])
        end
    end
    
    weightmatrix = squeeze(sum(gm_grad_mean(17:30,:,:), 1)); %%%%%%  use only recent : use last 140 time steps %%% 51*50
    weightmatrix(1,:) = - weightmatrix(1,:);  %%%%%%%%%%%%%%%%%   flip salt   %%%%%%%%%%%%%%%%%%%%%
    weightmatrix1 = [NaN(Mt+1,1), weightmatrix];  % (Mt+1, Mt+1)
    
    load(fullfile(common_cell_order_folder, ['sample' num2str(sampleID) '_common_outperm.mat']), 'outperm')

    if plot_matrix
        labnum = targetcellnames;
        figure('Position',[10,10,600,600]); % 600,600
        colormap greenwhitemazentablack(32); % bluewhitered(32) %jet(64) %gray
        image(weightmatrix([1,outperm+1],outperm)*96+33);
        title(['Sample' num2str(sampleID) ' connection weight matrix'],'Interpreter','none');
        h = gca;
        set(h,'ytick',(1:Mt+1)');
        set(h,'yticklabel',nameshifter(targetcellnames1([1,outperm+1])));
        set(h,'xtick',(1:Mt)');
        set(h,'xticklabel',nameshifter(labnum(outperm)));
        xtickangle(90);
        set(h,'fontsize',10);%14
        set(h,'Ticklength',[0 0]);
        set(h, 'TitleFontSizeMultiplier', 2);
        colorbar('Ticks',[0,16,32,48,64],...
            'TickLabels',{'-1','-0.5','0','0.5','1'})
        saveas(gcf, [results_folder '/weightmatrix_' link '_K' num2str(K) '_sample' num2str(sampleID) '.tif'])
        saveas(gcf, [results_folder '/weightmatrix_' link '_K' num2str(K) '_sample' num2str(sampleID) '.fig'])
        fig = gcf;
        fig.PaperPosition = [1 6 20 20];
        ax = gca;
        ax.XTickLabel = cellfun(@(x) [x,' '], ax.XTickLabel, 'UniformOutput', false);
        ax.YTickLabel = cellfun(@(x) [x,' '], ax.YTickLabel, 'UniformOutput', false);
        print(fig,'-painters',[results_folder '/weightmatrix_' link '_K' num2str(K) '_sample' num2str(sampleID) '.pdf'],'-dpdf')
    
        Tweight = array2table(weightmatrix([1,outperm+1],outperm));
        Tweight.Properties.VariableNames = labnum(outperm);
        Tweight.Properties.RowNames = ['none';labnum(outperm)];
        Tweight.Properties.DimensionNames{1} = 'Neurons';
        %writetable(Tweight,'weightmatrix.xls','WriteRowNames',true,'Sheet', 'MySheetName')
        save([results_folder '/Tweight_K' num2str(K) '_sample' num2str(sampleID) '.mat'], 'Tweight')
    
    end
    
    groupA_names = {'RIMR','AVAR','AVAL','AVEL','RIML','AIBL','AIBR','AVER','RMDL','RMDVR'}; %above 13 point, remove ASHL. %{'AVAL','AVAR','RIML','RIMR','AVEL','RIMR'};
    groupB_names = {'RMEL','RMER','RID','RMED','RIBR','RIBL','SMBVR','AVJL'}; %below -5 points. %{'RIBL','RIBR','RMEL','RMER','RID','AVBL','AVBR'};
    saltsensors = {'ASEL','ASER','BAGL','BAGR','AWCL','AWCR','ASHL','ASHR'};
    %salt = {'salt'};
    
    outpermA = []; % targetcellnames1-based (salt is no 1)
    outpermS = [];
    outpermB = [];
    outpermN = 1;
    for i = 1:length(outperm)
        if any(strcmp(targetcellnames{outperm(i)},groupA_names))
            outpermA = [outpermA, outperm(i)+1];
        elseif any(strcmp(targetcellnames{outperm(i)},groupB_names))
            outpermB = [outpermB, outperm(i)+1];
        elseif any(strcmp(targetcellnames{outperm(i)},saltsensors))
            outpermS = [outpermS, outperm(i)+1];
        end
    end
    outpermAB = [outpermA, outpermB];
    outpermSAB = [outpermS, outpermA, outpermB];
    outpermNSAB = [outpermN, outpermS, outpermA, outpermB];
    na = length(outpermA);
    nb = length(outpermB);
    ns = length(outpermS);
    
    draw_graph(weightmatrix1, outpermAB, targetcellnames1, [-ones(1,na), ones(1,nb)], zeros(1,na+nb), ['Sample ' num2str(sampleID) ', Group A/B'], fullfile(results_folder, ['AB_network_K' num2str(K) '_sample' num2str(sampleID)]))
    draw_graph(weightmatrix1, outpermSAB, targetcellnames1, [zeros(1,ns), -ones(1,na), ones(1,nb)], [ones(1,ns) -ones(1,na+nb)], ['Sample ' num2str(sampleID) ', Sensory-Group A/B'], fullfile(results_folder, ['SAB_network_K' num2str(K) '_sample' num2str(sampleID)]))
    draw_graph(weightmatrix1, outpermNSAB, targetcellnames1, [0, zeros(1,ns), -ones(1,na), ones(1,nb)], [1.5, ones(1,ns) -ones(1,na+nb)], ['Sample ' num2str(sampleID) ', Salt-Sensory-Group A/B'], fullfile(results_folder, ['NSAB_network_K' num2str(K) '_sample' num2str(sampleID)]))
    
end
end

function draw_graph(weightmatrix1, outpermX, targetcellnames, xstarts, ystarts, title, fileheader)

edgeweights = weightmatrix1(outpermX,outpermX);   % as of strength4, overall value need to be adjusted
edgeweights(isnan(edgeweights))=0;
graph = digraph(edgeweights, 'omitselfloops');
graphEdgeWeight = graph.Edges.Weight;
graph.Edges.Weight = abs(graphEdgeWeight);

nnodes = size(graph.Nodes,1);
nedges = size(graph.Edges,1);

LWidths = 3*abs(graph.Edges.Weight)/max(abs(graph.Edges.Weight)); %10 % 4
%disp(LWidths)
Asize = (1+LWidths.^0.6).*2.5; %.*abs(graph.Edges.Weight)/max(abs(graph.Edges.Weight)); % (2+LWidths.^0.3).*2.5;
Ecolor = ones(nedges,3);
Ecolor(graphEdgeWeight>0,2)=1-half_sigmoid(graphEdgeWeight(graphEdgeWeight>0).*30); % red
Ecolor(graphEdgeWeight>0,3)=1-half_sigmoid(graphEdgeWeight(graphEdgeWeight>0).*30); % red
Ecolor(graphEdgeWeight<0,1)=1-half_sigmoid(-graphEdgeWeight(graphEdgeWeight<0).*30); % blue
Ecolor(graphEdgeWeight<0,2)=1-half_sigmoid(-graphEdgeWeight(graphEdgeWeight<0).*30); % blue

figure('Position',[10,10,200,200]); % none
h=plot(graph, 'LineWidth',LWidths, 'ArrowSize', Asize, 'EdgeColor', Ecolor, 'Marker', 'none'); % 'EdgeLabel',graph.Edges.Weight
ax=gca;
ax.Title.String=title;

%layout(h,'force','XStart',linspace(-2,2,nnodes),'YStart',randn(1,nnodes), 'WeightEffect', 'inverse', 'Iterations', 100)
layout(h,'force','XStart', xstarts,'YStart', ystarts, 'WeightEffect', 'inverse', 'Iterations', 50)

if outpermX(1) == 1 % if salt included
    h.YData(1) = max(h.YData)+0.5;
    h.XData(1) = 0;
end

h.NodeLabel = '';
%nl = targetcellnames(outperm(1:ncells));
nl = targetcellnames(outpermX);
xd = get(h, 'XData');
yd = get(h, 'YData');
text(xd, yd, nl, 'FontSize',8, 'FontWeight','bold', 'HorizontalAlignment','center', 'VerticalAlignment','middle')

saveas(gcf, [fileheader '.tif'])
saveas(gcf, [fileheader '.fig'])

end


function sig = half_sigmoid(X)
sig = 2./(1+exp(-X))-1;
end



