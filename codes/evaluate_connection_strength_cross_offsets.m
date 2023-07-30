
function evaluate_connection_strength_cross_offsets(UGE_TASK_ID_text) % sampleID


sampleID = str2num(UGE_TASK_ID_text);
partition_number = 3;
offsets = 0:4;

project_folder = '/home/iino/gKDR_GMM_50mM/Linked_gKDR_freerun/publish_version';
%common_data_folder = fullfile(project_folder, 'common_data'); %'../ver92_Hbase/common_data';
gKDR_codes_folder = fullfile(project_folder, 'gKDR_codes');
utilities_folder = fullfile(project_folder, 'utilities'); %'../'
%metadata_folder = fullfile(project_folder, 'metadata'); %'../'
common_cell_order_folder = fullfile(project_folder, 'common_cell_order');
models_folder = fullfile(project_folder, 'models'); % ['../ver98_HLong/model_indirect_k' num2str(embed_width) '_tau' num2str(embed_step) '_kGMM' num2str(kGMM)];

Ks = [3];
kGMMs = 2;
link = 'indirect';
embed_width = 30; % number of embed_step's used for embedding = column number in source data
embed_step = 10; % invervals used for embedding (index-based)
plot_timeweights = false;
plot_matrix = false;

addpath(gKDR_codes_folder);
addpath(utilities_folder);

expsavefolder = fullfile(project_folder, 'connection_strength', [link '_cross' num2str(partition_number)]);
if ~exist(expsavefolder ,'dir')
    mkdir(expsavefolder)
end

for Ki = 1:1 %length(Ks)
    K = Ks(Ki);
    
    for kGMM = kGMMs

        for offset = offsets
        
        for partitioni = 1:partition_number
            
            disp(['< kGMM ' num2str(kGMM) ' >']);
            disp(['< K= ' num2str(K) ' >']);
            disp(['< offset= ' num2str(offset) ' >']);
            disp(['< part ' num2str(partitioni) ' >']);
            
            expsavefileheader = fullfile(expsavefolder, ['sample' num2str(sampleID)  '_K' num2str(K) '_kGMM' num2str(kGMM) '_offset' num2str(offset) '_part' num2str(partitioni)]);
            
            modelfolder = fullfile(models_folder, ['model_' link '_k' num2str(embed_width) '_tau' num2str(embed_step) '_kGMM' num2str(kGMM) '_cross' num2str(partition_number)]);
            modelfileheader = fullfile(modelfolder, ['sample' num2str(sampleID) '_K' num2str(K) '_offset' num2str(offset) '_part' num2str(partitioni)]);
            load([modelfileheader '_gradient.mat'], 'gm_grad_cell', 'gm_grad_array','gm_grad_mean','gm_grad_mean_cell','uniqNames','targetcells','selicell');
            
            targetcellnames = uniqNames(targetcells);
            targetcellnames1 = [{'salt'}; targetcellnames];
            Mt = length(targetcells);
            %disp('number of presynaptic cells:')
            %disp(arrayfun(@(x) length(selicell{x}), 1:Mt))
            
            if plot_timeweights
                for targeti = 1:5
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
                    saveas(gcf, [expsavefileheader '_timeweight_cell' num2str(targeti) '_indirect.tif'])
                    saveas(gcf, [expsavefileheader '_timeweight_cell' num2str(targeti) '_indirect.fig'])
                end
            end
            
            weightmatrix = squeeze(sum(gm_grad_mean .* reshape(0.8.^((embed_width-1):-1:0), [embed_width,1,1]), 1)); % 51*50
            weightmatrix1 = [NaN(Mt+1,1), weightmatrix];  % (Mt+1, Mt+1)
            
            save([expsavefileheader '_weightmatrix.mat'], 'weightmatrix','weightmatrix1')
            
            load(fullfile(common_cell_order_folder, ['sample' num2str(sampleID) '_common_outperm.mat']), 'outperm')
            
            if plot_matrix
                labnum = targetcellnames;
                figure('Position',[10,10,1000,600]);
                colormap greenwhitemazentablack(32); % bluewhitered(32) %jet(64) %gray
                image(weightmatrix(outperm+1,outperm)*96+33);
                title(['Weight matrix'],'Interpreter','none');
                h = gca;
                set(h,'ytick',(1:Mt)');
                set(h,'yticklabel',nameshifter(labnum(outperm)));
                set(h,'xtick',(1:Mt)');
                set(h,'xticklabel',nameshifter(labnum(outperm)));
                xtickangle(90);
                set(h,'fontsize',14);
                set(h,'Ticklength',[0 0]);
                colorbar('Ticks',[0,16,32,48,64],...
                    'TickLabels',{'-1','-0.5','0','0.5','1'})
                saveas(gcf, [expsavefileheader '_weightmatrix_indirect_K' num2str(K) '_sample' num2str(sampleID) '.tif'])
                saveas(gcf, [expsavefileheader '_weightmatrix_indirect_K' num2str(K) '_sample' num2str(sampleID) '.fig'])
            end
            
            groupA_names = {'RIMR','AVAR','AVAL','AVEL','RIML','AIBL','AIBR','AVER','RMDL','RMDVR'}; %above 13 point, remove ASHL. %{'AVAL','AVAR','RIML','RIMR','AVEL','RIMR'};
            groupB_names = {'RMEL','RMER','RID','RMED','RIBR','RIBL','SMBVR','AVJL'}; %below -5 points. %{'RIBL','RIBR','RMEL','RMER','RID','AVBL','AVBR'};
            saltsensors = {'ASEL','ASER','BAGL','BAGR','AWCL','AWCR','ASHL','ASHR'};
            salt = {'salt'};
            
            xposition = zeros(1,Mt);
            yposition = -2*ones(1,Mt);
            
            outpermA = []; % targetcellnames1-based (salt is no 1)
            outpermS = [];
            outpermB = [];
            outpermN = 1;
            for i = 1:length(outperm)
                if any(strcmp(targetcellnames{outperm(i)},groupA_names))
                    outpermA = [outpermA, outperm(i)+1];
                    xposition(i) = -2;
                    yposition(i) = -1;
                elseif any(strcmp(targetcellnames{outperm(i)},groupB_names))
                    outpermB = [outpermB, outperm(i)+1];
                    xposition(i) = 2;
                    yposition(i) = -1;
                elseif any(strcmp(targetcellnames{outperm(i)},saltsensors))
                    outpermS = [outpermS, outperm(i)+1];
                    xposition(i) = 0;
                    yposition(i) = 1;
                end
            end
            xposition = [0, xposition];
            yposition = [1.5, yposition];
            
            outpermAB = [outpermA, outpermB];
            outpermSAB = [outpermS, outpermA, outpermB];
            outpermNSAB = [outpermN, outpermS, outpermA, outpermB];
            na = length(outpermA);
            nb = length(outpermB);
            ns = length(outpermS);
            outpermAll = [1, outperm+1];
            
            %draw_graph(weightmatrix1, outpermAB, targetcellnames1, [-ones(1,na), ones(1,nb)], zeros(1,na+nb), ['Sample ' num2str(sampleID) ', Group A/B'], ['AB_network_K' num2str(K) '_sample' num2str(sampleID)])
            %draw_graph(weightmatrix1, outpermSAB, targetcellnames1, [zeros(1,ns), -ones(1,na), ones(1,nb)], [ones(1,ns) -ones(1,na+nb)], ['Sample ' num2str(sampleID) ', Sensory-Group A/B'], ['SAB_network_K' num2str(K) '_sample' num2str(sampleID)])
            %draw_graph(weightmatrix1, outpermNSAB, targetcellnames1, [0, zeros(1,ns), -ones(1,na), ones(1,nb)], [1.5, ones(1,ns) -ones(1,na+nb)], ['Sample ' num2str(sampleID) ', Salt-Sensory-Group A/B'], ['NSAB_network_K' num2str(K) '_sample' num2str(sampleID)])
            draw_graph(weightmatrix1, outpermAll, targetcellnames1, xposition, yposition, ['Sample ' num2str(sampleID) ', All'], [expsavefileheader '_All_network'])
        end
        end
    end
end
end

function draw_graph(weightmatrix1, outpermX, targetcellnames, xstarts, ystarts, title, fileheader)

edgeweights = weightmatrix1(outpermX,outpermX);
edgeweights(isnan(edgeweights))=0;
graph = digraph(edgeweights, 'omitselfloops');
graphEdgeWeight = graph.Edges.Weight;
graph.Edges.Weight = abs(graphEdgeWeight);

nnodes = size(graph.Nodes,1);
nedges = size(graph.Edges,1);

LWidths = 10*abs(graph.Edges.Weight)/max(abs(graph.Edges.Weight));
%disp(LWidths)
Asize = (2+LWidths.^0.3).*2.5; %.*abs(graph.Edges.Weight)/max(abs(graph.Edges.Weight));
Ecolor = ones(nedges,3);
Ecolor(graphEdgeWeight>0,2)=1-half_sigmoid(graphEdgeWeight(graphEdgeWeight>0).*30); % red
Ecolor(graphEdgeWeight>0,3)=1-half_sigmoid(graphEdgeWeight(graphEdgeWeight>0).*30); % red
Ecolor(graphEdgeWeight<0,1)=1-half_sigmoid(-graphEdgeWeight(graphEdgeWeight<0).*30); % blue
Ecolor(graphEdgeWeight<0,2)=1-half_sigmoid(-graphEdgeWeight(graphEdgeWeight<0).*30); % blue

fig = figure;
fig.PaperType       = 'a4';
fig.PaperOrientation = 'landscape';
fig.PaperUnits = 'centimeter';
fig.PaperPosition = [0,0,29.7,21];
fig.Position = [-500,700,890,630];
h=plot(graph, 'LineWidth',LWidths, 'ArrowSize', Asize, 'EdgeColor', Ecolor); % 'EdgeLabel',graph.Edges.Weight
ax=gca;
ax.Title.String=title;

%layout(h,'force','XStart',linspace(-2,2,nnodes),'YStart',randn(1,nnodes), 'WeightEffect', 'inverse', 'Iterations', 100)
layout(h,'force','XStart', xstarts,'YStart', ystarts, 'WeightEffect', 'inverse', 'Iterations', 50)

if outpermX(1) == 1 % if salt included
    h.YData(1) = max(h.YData)+0.5;
    h.XData(1) = mean(h.XData);
end

h.NodeLabel = '';
%nl = targetcellnames(outperm(1:ncells));
nl = targetcellnames(outpermX);
xd = get(h, 'XData');
yd = get(h, 'YData');
text(xd, yd, nl, 'FontSize',10, 'FontWeight','bold', 'HorizontalAlignment','center', 'VerticalAlignment','middle')

saveas(gcf, [fileheader '.tif'])
saveas(gcf, [fileheader '.fig'])
saveas(gcf, [fileheader '.pdf'])

end


function sig = half_sigmoid(X)
sig = 2./(1+exp(-X))-1;
end



