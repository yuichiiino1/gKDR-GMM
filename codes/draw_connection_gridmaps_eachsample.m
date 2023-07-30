function draw_connection_gridmaps_eachsample()

global samples common_data_folder connection_strength_folder globalvalidweightarray  ...
    save_each_sample_fig save_each_sample_pdf globalNames conNames dirtext Mt targetcellnames  cellg cellg1 samplei ...
    multiconmatrix FDRthre2 show_nonFDR partition_number offsets Ks draw_box_when_connected show_box_when_insignificant Mg  ...
    save_each_sample_tif save_all_valid_weights_draw_single_pdf embed_size embed_step embed_width link models_folder kGMM ...
    chemcon_thickness_matrix gapcon_thickness_matrix chemcon_thickness_matrixgt1 gapcon_thickness_matrixgt1 path_to_synapse_table

samples = 1:24;
FDRthre = 0.005; %0.05;
FDRthre2 = 0.005; % just affect core neurons
only_direct = false; %true; %false;
show_nonFDR = false; %true;
partition_number = 3;
offsets = 0:4;
link = 'indirect';
Ks = 3:5; 
kGMM = 2;
validmethod = 'max'; %'min';  % actually max is used in this version
embed_width = 30; %10; % number of embed_step's used for embedding = column number in source data 
embed_step = 10; %10; % invervals used for embedding (index-based)  
embed_size = embed_width;

save_each_sample_fig = false;
save_each_sample_pdf = false;
draw_box_when_connected = true;
show_box_when_insignificant = false;
save_each_sample_tif = false;
save_all_valid_weights_draw_single_pdf = true;

draw_pcolor_figure = true;
draw_brickwork_figure = false;
draw_lagweightsplot_figure = false;

if only_direct
    dirtext = '_direct_link';
else
    dirtext = '_all_link';
end
if show_nonFDR
    FDRtext = '_shownonFDR';
else
    FDRtext = '_onlyFDR';
end

%close all

project_folder = '/home/iino/gKDR_GMM_50mM/Linked_gKDR_freerun/publish_version';
%project_folder = '/Users/iino/Documents/___Transcend J/_4DImaging_Data_Analysis/gKDR-GMM/gKDR_GMM_50mM/Linked_gKDR_freerun/publish_version/';
%project_folder = '/Users/iino/Documents/___Transcend J/_4DImaging_Data_Analysis/gKDR-GMM/gKDR_GMM_50mM/Linked_gKDR_freerun/publish_version_test/';
connection_strength_folder = fullfile(project_folder, 'connection_strength', [link '_cross' num2str(partition_number) validmethod dirtext '_FDR' num2str(FDRthre)]);
common_data_folder = fullfile(project_folder, 'common_data');
%gKDR_codes_folder = fullfile(project_folder, 'gKDR_codes');
utilities_folder = fullfile(project_folder, 'utilities');
metadata_folder = fullfile(project_folder, 'metadata');
%common_cell_order_folder = fullfile(project_folder, 'common_cell_order');
models_folder = fullfile(project_folder, 'models');
%evaluateL_results_folder = fullfile(project_folder, 'evaluateL_results');
path_to_synapse_table = fullfile(metadata_folder,'synapse_All_Iwasaki2.txt'); % used only when draw_box_when_connected = true;

%addpath(gKDR_codes_folder);
addpath(utilities_folder);

globalNames = table2cell(readtable(fullfile(metadata_folder,'globalNames.csv'),'ReadVariableNames',false));
Mg = length(globalNames);
disp(['Mg=',num2str(Mg)])
multiconmatrix = table2array(readtable(fullfile(metadata_folder, 'multiconmatrix.csv'),'ReadVariableNames',false));
conNames = table2cell(readtable(fullfile(metadata_folder, 'conneurons.csv'),'ReadVariableNames',false));

load(fullfile(connection_strength_folder, 'process_weight_allnetwork_core_phys_offsets.mat'),'globalvalidweightarray')

if draw_pcolor_figure
savefilenameheader = ['sample' num2str(samples(1)) '-' num2str(samples(end)) '_pcolor'];
draw_pcolor(savefilenameheader)  % draw cell interaction matrix showing p-values and sign of interaction
end

if draw_brickwork_figure
savefilenameheader2 = ['sample' num2str(samples(1)) '-' num2str(samples(end)) '_brickwork'];
draw_brickwork(savefilenameheader2)  % draw cell interaction matrix with grid representation of each condition
end

if draw_lagweightsplot_figure
savefilenameheader3 = ['sample' num2str(samples(1)) '-' num2str(samples(end)) '_lagweights'];
draw_lagweightsplot(savefilenameheader3)
end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%     Function  draw_pcolor     %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_pcolor(savefilenameheader)

global samples common_data_folder connection_strength_folder globalvalidweightarray  ...
    save_each_sample_fig save_each_sample_pdf globalNames conNames dirtext Mt targetcellnames  cellg cellg1 samplei ...
    multiconmatrix FDRthre2 show_nonFDR partition_number offsets Ks draw_box_when_connected show_box_when_insignificant Mg  ...
    save_each_sample_tif save_all_valid_weights_draw_single_pdf embed_size embed_step embed_width link models_folder kGMM ...
    chemcon_thickness_matrix gapcon_thickness_matrix chemcon_thickness_matrixgt1 gapcon_thickness_matrixgt1 path_to_synapse_table

%%%%%%%%%   Charts for each sample  %%%%%%%

if draw_box_when_connected == true
    savefilenameheader = [savefilenameheader,'_box'];
end

savefilenamesingle = fullfile(connection_strength_folder, [savefilenameheader '_single.ps']);
if save_all_valid_weights_draw_single_pdf && exist(savefilenamesingle, 'file')
    delete(savefilenamesingle);
end

fig = figure;
fig.PaperType       = 'a4';
fig.PaperOrientation = 'portrait';
fig.Position        = [50,50,1000,1000];

%%%%%%%%%%%%  sample loop  %%%%%%%%%%%%%%%%%%%%

for samplei = 1:length(samples)
    sampleID = samples(samplei);
    disp(['sample ' num2str(sampleID)])

    % load basic data
    load(fullfile(common_data_folder,['sample' num2str(sampleID) '_data.mat']), 'uniqNames', 'targetcellnames', 'targetcells', 'Mt'); % 'data',
    targetcellnames1 = [{'salt'}; targetcellnames];

    %%%%% extract from global array %%%%
    cellg = [];   %indexes of targetcells in globalNames
    for targeti = 1:Mt
        cellname = targetcellnames{targeti};
        celli = find(strcmp(globalNames, cellname));
        assert(length(celli)==1)
        cellg = [cellg, celli];
    end
    cellg1 = [1, cellg+1]; % first is salt

    [validweightarrayall , p_all_mat, FDRallsig, validsource, validtarget] = select_validneurons();

    if draw_box_when_connected
        %chemcon_thickness = chemcon_thickness_matrixgt1(validsource,validtarget+1);
        %gapcon_thickness = gapcon_thickness_matrixgt1(validsource,validtarget+1);
        chemcon = chemcon_thickness_matrixgt1(validsource,validtarget+1) > 0; % chemcon_thickness_matrixgt1 =(Mt+1,Mt+1)
        gapcon = gapcon_thickness_matrixgt1(validsource,validtarget+1) > 0;
    end


    ncolor = 33;
    left_m = 0.1;
    bot_m = 0.05;
    right_m = 0.05;
    top_m = 0.15;

    ncol = length(validtarget);
    nrow = length(validsource);
    figcount = 0;

    disp(['drawing ' num2str(nrow) ' rows.'])

    logpvaluemat = NaN(nrow,ncol);
    for rowi = 1:nrow
        for coli = 1:ncol
            figcount = figcount + 1;
            row = validsource(rowi);
            col = validtarget(coli);

            logpvalue = -log10(p_all_mat(row,col))/30;
            if logpvalue > 1
                logpvalue = 1;
            end
            if nanmean(reshape(validweightarrayall(row,col, :, :),1,[])) < 0
                logpvalue = - logpvalue;
            end
            if FDRallsig(row,col) || show_nonFDR
                logpvaluemat(rowi,coli) = logpvalue;
            end

        end
    end

    %%%%%%%%%%%%%%%%      Draw all at once      %%%%%%%%%%%%%%%%%%%%%%

    ax = axes('Position',[left_m, bot_m, 1-right_m-left_m, 1-top_m-bot_m]);
    image(ax,logpvaluemat*ncolor+ncolor+2);
    colormap bluewhiteredblack(33)

    ax.XAxisLocation='top';
    xticks(1:ncol)
    yticks(1:nrow)
    xticklabels(targetcellnames(validtarget(1:ncol)))
    sourcecelllabels = {};
    for i = 1:nrow
        sourcecelllabels = [sourcecelllabels; {[targetcellnames1{validsource(i)} ' ']}];
    end
    %yticklabels(targetcellnames1(validsource(1:nrow)))
    yticklabels(sourcecelllabels)
    if ncol > 80
        ax.FontSize = 6;
    elseif ncol > 60
        ax.FontSize = 8;
    elseif ncol > 40
        ax.FontSize = 10;
    else
        ax.FontSize = 12;
    end
    ax.XTickLabelRotation = 90;
    ax.Box = 'on';

    for i = 1:ncol-1
        xline(i+0.5, 'Color', [0.5,0.5,0.5])
    end
    for i = 1:nrow-1
        yline(i+0.5, 'Color', [0.5,0.5,0.5])
    end

    if draw_box_when_connected
        for i = 1:length(validsource)
            for j = 1:length(validtarget)
                if show_box_when_insignificant || ~isnan(logpvaluemat(i,j))
                    y1 = i - 0.45;
                    x1 = j - 0.45;
                    if gapcon(i,j) && chemcon(i,j) % both
                        rectangle('Position',[x1,y1,0.9,0.9], 'LineStyle', '-', 'EdgeColor',[0,0.6,0],'LineWidth',2)
                    elseif gapcon(i,j) && ~chemcon(i,j) % gap only
                        rectangle('Position',[x1,y1,0.9,0.9], 'LineStyle', '--', 'EdgeColor',[0,0.6,0],'LineWidth',2)
                    elseif ~gapcon(i,j) && chemcon(i,j)  % chemical only
                        rectangle('Position',[x1,y1,0.9,0.9], 'LineStyle', ':', 'EdgeColor',[0,0.6,0],'LineWidth',2)
                    end
                end
            end
        end
    end

    c = colorbar;
    c.Ticks = [1.5,2.5:11:69];
    c.TickLabels = [{'na'}; arrayfun(@(x) {num2str(x)}, (-30:10:30)')];
    c.Label.String = '+/- log_{10}(1/p)';

    ax2 = annotation('textbox','String',['Sample ' num2str(sampleID)], 'FitBoxToText','on','LineStyle','none','FontSize',24,'HorizontalAlignment','center','VerticalAlignment','middle');
    ax2.Position = [0.35,0.95,0.3,0.05];

    ax3 = annotation('textarrow','String','From', 'HeadStyle','none','LineStyle','none','FontSize',20,'TextRotation',90,'HorizontalAlignment','center','VerticalAlignment','middle');
    ax3.Position = [0,0.4,0.05,0.3];
    ax4 = annotation('textbox','String','To', 'FitBoxToText','on','LineStyle','none','FontSize',20,'HorizontalAlignment','center','VerticalAlignment','middle');
    ax4.Position = [0.35,0.9,0.3,0.05];

    disp('Saving all_valid_weights figure...')
    orient(fig,'tall')
    if save_each_sample_tif
        saveas(gcf, fullfile(connection_strength_folder, [savefilenameheader '_sample' num2str(sampleID) '.tif']))
    end
    if save_each_sample_fig
        saveas(gcf, fullfile(connection_strength_folder, [savefilenameheader  '_sample' num2str(sampleID) '.fig']))
    end
    if save_each_sample_pdf
        saveas(gcf, fullfile(connection_strength_folder, [savefilenameheader  '_sample' num2str(sampleID) '.pdf']))
    end

    if save_all_valid_weights_draw_single_pdf
        print(fig, savefilenamesingle,'-dpsc','-painters','-append');
    end
    clf(fig)

end  % end for samplei

if ~strcmp(computer, 'MACI64')   % Mac does not support use of ps2pdf like this
    system(['ps2pdf ' savefilenamesingle]);
    pdffilename = [savefilenameheader '_single.pdf'];
    disp(['mv ' pdffilename ' ' connection_strength_folder])
    system(['mv ' pdffilename ' ' connection_strength_folder]);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%     Function  draw_brickwork     %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function   draw_brickwork(savefilenameheader)

global samples common_data_folder connection_strength_folder globalvalidweightarray  ...
    save_each_sample_fig save_each_sample_pdf globalNames conNames dirtext Mt targetcellnames  cellg cellg1 samplei ...
    multiconmatrix FDRthre2 show_nonFDR partition_number offsets Ks draw_box_when_connected show_box_when_insignificant Mg  ...
    save_each_sample_tif save_all_valid_weights_draw_single_pdf embed_size embed_step embed_width link models_folder kGMM ...
    chemcon_thickness_matrix gapcon_thickness_matrix chemcon_thickness_matrixgt1 gapcon_thickness_matrixgt1 path_to_synapse_table

%%%%%%%%%   Charts for each sample  %%%%%%%

savefilename_single = fullfile(connection_strength_folder, [savefilenameheader '_single.ps']);
if exist(savefilename_single, 'file')
    delete(savefilename_single);
end

%figure('Position',[10,10,1000,1000])
fig = figure;
fig.PaperType       = 'a4';
fig.PaperOrientation = 'portrait';
fig.Position        = [50,50,1000,1000];

left_m = 0.1;
bot_m = 0.05;
right_m = 0.05;
top_m = 0.15;
ver_r = 1.1;
col_r = 1.1;

%%%%%%%%%%%%  sample loop  %%%%%%%%%%%%%%%%%%%%

for samplei = 1:length(samples)
    sampleID = samples(samplei);
    disp(['< sample ' num2str(sampleID) ' >']);

    % load basic data
    load(fullfile(common_data_folder,['sample' num2str(sampleID) '_data.mat']), 'uniqNames', 'targetcellnames', 'targetcells', 'Mt'); % 'data',
    targetcellnames1 = [{'salt'}; targetcellnames];

    %%%%% extract from global array %%%%

    cellg = [];   %indexes of targetcells in globalNames
    for targeti = 1:Mt
        cellname = targetcellnames{targeti};
        celli = find(strcmp(globalNames, cellname));
        assert(length(celli)==1)
        cellg = [cellg, celli];
    end
    cellg1 = [1, cellg+1]; % first is salt

    validweightarray = squeeze(globalvalidweightarray(cellg1, cellg, :,:,:, samplei)); % (Mt+1,Mt,partition_number,length(offsets),length(Ks))

    [validweightarrayall , p_all_mat, FDRallsig, validsource, validtarget] = select_validneurons();

    figcount = 0;
    ncol = length(validtarget);
    nrow = length(validsource);
    ncolor = 32;
    disp(['drawing ' num2str(nrow) ' rows.'])
    for rowi = 1:nrow
        for coli = 1:ncol
            figcount = figcount + 1;
            row = validsource(rowi);
            col = validtarget(coli);

            ax(figcount) = axes('Position',...
                [(1-right_m-left_m)*(mod(figcount-1,ncol))/ncol + left_m ,...
                (1-top_m-bot_m)*(1-ceil(figcount/ncol)/(nrow)) + bot_m ,...
                (1-right_m-left_m)/(ncol*col_r ),...
                (1-top_m-bot_m)/(nrow*ver_r)]...
                );

            if   FDRallsig(row,col) || show_nonFDR

                weight2d = reshape(validweightarray(row,col,:,:,:),partition_number, length(offsets)*length(Ks))';
                maxweight = nanmax(nanmax(abs(weight2d)));
                weight2dnorm = weight2d/maxweight;

                image(ax(figcount),weight2dnorm*ncolor+ncolor+2)
                colormap bluewhiteredblack(32)

            end
            xticks([])
            yticks([])
            xticklabels([])
            yticklabels([])
            ax(figcount).Box = 'on';
            %ax = gca;
            %ax.Box = 'on';
            if rowi == 1
                ax(figcount).Title.String = [' ' targetcellnames{col}]; %title(targetcellnames{col})
                ax(figcount).Title.Rotation = 90;
                ax(figcount).Title.VerticalAlignment = 'middle';
                ax(figcount).Title.HorizontalAlignment = 'left';
            end
            if coli == 1
                ylabel([targetcellnames1{row} '  '])
                hYLabel = get(gca,'YLabel');
                set(hYLabel,'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right')
            end
        end
    end
    ax2 = annotation('textbox','String',['Sample ' num2str(sampleID)],'FitBoxToText','on','LineStyle','none','FontSize',26,'HorizontalAlignment','center','VerticalAlignment','middle');
    ax2.Position = [0.35,0.95,0.3,0.05];
    ax3 = annotation('textarrow','String','From', 'HeadStyle','none','LineStyle','none','FontSize',20,'TextRotation',90,'HorizontalAlignment','center','VerticalAlignment','middle');
    ax3.Position = [0,0.4,0.05,0.3];
    ax4 = annotation('textbox','String','To', 'FitBoxToText','on','LineStyle','none','FontSize',20,'HorizontalAlignment','center','VerticalAlignment','middle');
    ax4.Position = [0.35,0.9,0.3,0.05];

    disp('Saving all_valid_weights figure...')
    orient(fig,'tall')
    saveas(gcf, fullfile(connection_strength_folder, [savefilenameheader '_sample' num2str(sampleID) '.tif']))
    if save_each_sample_fig
        saveas(gcf, fullfile(connection_strength_folder, [savefilenameheader '_sample' num2str(sampleID) '.fig']))
    end
    if save_each_sample_pdf
        saveas(gcf, fullfile(connection_strength_folder, [savefilenameheader '_sample' num2str(sampleID) '.pdf']))
    end

    print(fig, savefilename_single,'-dpsc','-painters','-append');
    clf(fig)
    disp('Done.')

end  % end for samplei


if ~strcmp(computer, 'MACI64')   % Mac does not support use of ps2pdf like this
    system(['ps2pdf ' savefilenamesingle]);
    pdffilename = [savefilenameheader '_single.pdf'];
    disp(['mv ' pdffilename ' ' connection_strength_folder])
    system(['mv ' pdffilename ' ' connection_strength_folder]);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%     Function  draw_lagweightsplot    %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_lagweightsplot(savefilenameheader)

global samples common_data_folder connection_strength_folder globalvalidweightarray  ...
    save_each_sample_fig save_each_sample_pdf globalNames conNames dirtext Mt targetcellnames  cellg cellg1 samplei ...
    multiconmatrix FDRthre2 show_nonFDR partition_number offsets Ks draw_box_when_connected show_box_when_insignificant Mg  ...
    save_each_sample_tif save_all_valid_weights_draw_single_pdf embed_size embed_step embed_width link models_folder kGMM ...
    chemcon_thickness_matrix gapcon_thickness_matrix chemcon_thickness_matrixgt1 gapcon_thickness_matrixgt1 path_to_synapse_table

%%%%%%%%%   Charts for each sample  %%%%%%%

savefilenamesingle = fullfile(connection_strength_folder, [savefilenameheader '_single.ps']);
if exist(savefilenamesingle, 'file')
    delete(savefilenamesingle);
end

%figure('Position',[10,10,1000,1000])
fig = figure;
fig.PaperType       = 'a4';
fig.PaperOrientation = 'portrait';
fig.Position        = [10,10,1000,1000];

left_m = 0.15;
bot_m = 0.05;
right_m = 0.05;
top_m = 0.15;
ver_r = 1.1;
col_r = 1.1;

%%%%%%%%%%%%  sample loop  %%%%%%%%%%%%%%%%%%%%

for samplei = 1:length(samples)
    sampleID = samples(samplei);
    disp(['< sample ' num2str(sampleID) ' >']);

    % load basic data
    load(fullfile(common_data_folder,['sample' num2str(sampleID) '_data.mat']), 'uniqNames', 'targetcellnames', 'targetcells', 'Mt'); % 'data',
    targetcellnames1 = [{'salt'}; targetcellnames];

    gm_grad_mean_array = zeros(embed_size,Mt+1,Mt,partition_number,length(offsets),length(Ks));
    for partitioni = 1:partition_number
        for offset = offsets
            disp(['loading gradient data for part ', num2str(partitioni), ' offset ' num2str(offset)])
            for Ki = 1:length(Ks)
                K=Ks(Ki);
                % load gradient data
                modelfolder = fullfile(models_folder, ['model_' link '_k' num2str(embed_width) '_tau' num2str(embed_step) '_kGMM' num2str(kGMM) '_cross' num2str(partition_number)]);
                modelfileheader = fullfile(modelfolder, ['sample' num2str(sampleID) '_K' num2str(K) '_offset' num2str(offset) '_part' num2str(partitioni)]);
                load([modelfileheader '_gradient.mat'], 'gm_grad_cell', 'gm_grad_array','gm_grad_mean','gm_grad_mean_cell','uniqNames','targetcells','selicell');
                gm_grad_mean_array(:,:,:,partitioni,offset+1,Ki) = gm_grad_mean;
            end
        end
    end
    gm_grad_mean_array(:,1,:,:,:,:) = - gm_grad_mean_array(:,1,:,:,:,:); % salt profile is flipped


    %%%%% extract from global array %%%%
    cellg = [];   %indexes of targetcells in globalNames
    for targeti = 1:Mt
        cellname = targetcellnames{targeti};
        celli = find(strcmp(globalNames, cellname));
        assert(length(celli)==1)
        cellg = [cellg, celli];
    end
    cellg1 = [1, cellg+1]; % first is salt

    validweightarray = squeeze(globalvalidweightarray(cellg1, cellg, :,:,:, samplei)); % (Mt+1,Mt,partition_number,length(offsets),length(Ks))

    [validweightarrayall , p_all_mat, FDRallsig, validsource, validtarget] = select_validneurons();

    figcount = 0;
    ncol = length(validtarget);
    nrow = length(validsource);
    disp(['drawing ' num2str(nrow) ' rows.'])
    for rowi = 1:nrow
        for coli = 1:ncol
            figcount = figcount + 1;
            row = validsource(rowi);
            col = validtarget(coli);

            ax(figcount) = axes('Position',...
                [(1-right_m-left_m)*(mod(figcount-1,ncol))/ncol + left_m ,...
                (1-top_m-bot_m)*(1-ceil(figcount/ncol)/(nrow)) + bot_m ,...
                (1-right_m-left_m)/(ncol*col_r ),...
                (1-top_m-bot_m)/(nrow*ver_r)]...
                );

            if   FDRallsig(row,col) || show_nonFDR

                % gm_grad_mean_array 30    51    50     3     5     3

                weight3d = reshape(validweightarray(row,col,:,:,:),partition_number, length(offsets),length(Ks));  % size 3  5  3
                gm_grad_ii = gm_grad_mean_array(:,row,col,:,:,:);   % size   30     1     1     3     5     3
                gm_grad_ii(:,:,:,isnan(weight3d)) = NaN;
                grad_embed = squeeze(nanmean(gm_grad_ii, 2:6)); % column vector
                grad_embed_sem = nanstd(gm_grad_ii,0, 2:6); %./sqrt(sum(~isnan(gm_grad_ii), 2:6));

                a = area(ax(figcount),[grad_embed-grad_embed_sem, grad_embed_sem*2]);
                a(1).FaceColor = 'None';
                a(1).EdgeColor = 'None';
                a(2).FaceColor = [0, 0, 1];
                a(2).FaceAlpha = 0.1;
                a(2).EdgeColor = 'None';
                hold on
                plot(ax(figcount),grad_embed, 'k', 'Linewidth', 0.3)
                yline(ax(figcount),0, 'Linewidth', 0.3)

            end
            xticks([])
            yticks([])
            xticklabels([])
            yticklabels([])
            ax(figcount).Box = 'on';
            if rowi == 1
                ax(figcount).Title.String = [' ' targetcellnames{col}]; %title(targetcellnames{col})
                ax(figcount).Title.Rotation = 90;
                ax(figcount).Title.VerticalAlignment = 'middle';
                ax(figcount).Title.HorizontalAlignment = 'left';
            end
            if coli == 1
                ylabel([targetcellnames1{row} '  '])
                hYLabel = get(gca,'YLabel');
                set(hYLabel,'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right')
            end
        end
    end
    ax2 = annotation('textbox','String',['Sample ' num2str(sampleID)],'FitBoxToText','on','LineStyle','none','FontSize',26,'HorizontalAlignment','center','VerticalAlignment','middle');
    ax2.Position = [0.35,0.95,0.3,0.05];
    ax3 = annotation('textarrow','String','From', 'HeadStyle','none','LineStyle','none','FontSize',20,'TextRotation',90,'HorizontalAlignment','center','VerticalAlignment','middle');
    ax3.Position = [0.05,0.4,0.05,0.3];
    ax4 = annotation('textbox','String','To', 'FitBoxToText','on','LineStyle','none','FontSize',20,'HorizontalAlignment','center','VerticalAlignment','middle');
    ax4.Position = [0.35,0.9,0.3,0.05];

    disp('Saving all_valid_grads figure...')
    orient(fig,'tall')
    saveas(gcf, fullfile(connection_strength_folder, [savefilenameheader '_sample' num2str(sampleID) '.tif']))
    if save_each_sample_fig
        saveas(gcf, fullfile(connection_strength_folder, [savefilenameheader '_sample' num2str(sampleID) '.fig']))
    end
    if save_each_sample_pdf
        saveas(gcf, fullfile(connection_strength_folder, [savefilenameheader '_sample' num2str(sampleID) '.pdf']))
    end

    print(fig, savefilenamesingle,'-dpsc','-painters','-append');
    clf(fig)
    disp('Done.')
end

if ~strcmp(computer, 'MACI64')
    system(['ps2pdf ' savefilenamesingle]);
    pdffilename = [savefilenameheader '_single.pdf'];
    disp(['mv ' pdffilename ' ' connection_strength_folder])
    system(['mv ' pdffilename ' ' connection_strength_folder]);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%     Function  select_validneurons     %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [validweightarrayall , p_all_mat, FDRallsig, validsource, validtarget] = select_validneurons()

global samples common_data_folder connection_strength_folder globalvalidweightarray  ...
    save_each_sample_fig save_each_sample_pdf globalNames conNames dirtext Mt targetcellnames  cellg cellg1 samplei ...
    multiconmatrix FDRthre2 show_nonFDR partition_number offsets Ks draw_box_when_connected show_box_when_insignificant Mg  ...
    save_each_sample_tif save_all_valid_weights_draw_single_pdf embed_size embed_step embed_width link models_folder kGMM ...
    chemcon_thickness_matrix gapcon_thickness_matrix chemcon_thickness_matrixgt1 gapcon_thickness_matrixgt1 path_to_synapse_table

validweightarrayall = squeeze(globalvalidweightarray(cellg1,cellg,:,:,:,samplei)); % (Mt+1,Mt,partition_number,length(offsets),length(Ks))

p_all_mat = NaN(Mt+1,Mt);
%sigweights = NaN(Mt+1,Mt);
for row = 1:Mt+1
    for col = 1:Mt
        if row~=col+1
            weights = reshape(validweightarrayall(row,col,:,:,:),1,[]);
            if all(isnan(weights))
                %p=1;
            else
                %[h,p] = ttest(weights(~isnan(weights)));
                p = ranksum(weights, -weights);
                p_all_mat(row,col) = p;
            end
        end
        %sigweights(row,col) = nanmean(weights);
    end
end

% connectome

ct = NaN(1,Mt);  % targetcellnames -> conNames
for i = 1:Mt
    found = find(strcmp(conNames, targetcellnames{i}));
    assert(length(found)==1)
    ct(i) = found;
end

% gapcon, chemcon
if draw_box_when_connected
    make_connectivity_matrix_thickness()
    chemcon_thickness_matrixgt = chemcon_thickness_matrix(ct,ct);
    chemcon_thickness_matrixgt1 = [zeros(1,Mt+1); [zeros(Mt,1),chemcon_thickness_matrixgt]];
    gapcon_thickness_matrixgt = gapcon_thickness_matrix(ct,ct);
    gapcon_thickness_matrixgt1 = [zeros(1,Mt+1); [zeros(Mt,1),gapcon_thickness_matrixgt]];
end

for row = 2:Mt+1
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
            p_all_mat(row,col) = NaN;
        end
    end
end

FDRallmat = reshape(mafdr(reshape(p_all_mat,[],1), 'BHFDR', true),Mt+1,Mt);
FDRallsig = FDRallmat <= FDRthre2;   % NaNはfalseになる。
for i = 1:Mt
    FDRallsig(i+1,i) = false;
end


validtargetbool = any(FDRallsig, 1);
validsourcebool = any(FDRallsig, 2)';
validtarget = 1:Mt;
validtarget = validtarget(validtargetbool);
validsource = 1:Mt+1;
validsource = validsource(validsourcebool);

disp(['number of validtarget is ' num2str(length(validtarget))])
disp(['number of validsource is ' num2str(length(validsource))])

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%    Function  make_connectivity_matrix_thickness     %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function make_connectivity_matrix_thickness()

global samples common_data_folder connection_strength_folder globalvalidweightarray  ...
    save_each_sample_fig save_each_sample_pdf globalNames conNames dirtext Mt targetcellnames  cellg cellg1 samplei ...
    multiconmatrix FDRthre2 show_nonFDR partition_number offsets Ks draw_box_when_connected show_box_when_insignificant Mg  ...
    save_each_sample_tif save_all_valid_weights_draw_single_pdf embed_size embed_step embed_width link models_folder kGMM ...
    chemcon_thickness_matrix gapcon_thickness_matrix chemcon_thickness_matrixgt1 gapcon_thickness_matrixgt1 path_to_synapse_table

% based on preparation_codes/make_connectivity_matrix.m
% requirement:  synapse_All_Iwasaki2.txt

synapse_data = readtable(path_to_synapse_table);
preneuron = synapse_data.From;
postneuron = synapse_data.To;
neurons = sort(unique([preneuron; postneuron]));
chemconnections = table2array(synapse_data(:,3:4));
gapconnections = table2array(synapse_data(:,5:6));

chemcon_thickness_matrix = zeros(length(neurons));
gapcon_thickness_matrix = zeros(length(neurons));

for i=1:size(synapse_data,1)
        j = find(strcmp(preneuron(i),neurons));
        k = find(strcmp(postneuron(i),neurons));
        if chemcon_thickness_matrix(j,k) ~= 0
            disp([i,j,k]) % none reported
        end
        chemcon_thickness_matrix(j,k) = mean(chemconnections(i,:));

        if (gapcon_thickness_matrix(j,k) ~= 0 || gapcon_thickness_matrix(k,j) ~= 0) ...
                && gapcon_thickness_matrix(j,k)~=mean(gapconnections(i,:))
            % There are several repeated entries that do not match each
            % other. They are VA5 - AVAL/R and VB4 - AVBL. Just ignore.
            %disp([i,j,k])
            %disp(neurons([j,k])')
            %disp([gapcon_thickness_matrix(j,k),mean(gapconnections(i,:))])
        end
        gapcon_thickness_matrix(j,k) = mean(gapconnections(i,:));
        gapcon_thickness_matrix(k,j) = mean(gapconnections(i,:)); % no directionality
        
end

end