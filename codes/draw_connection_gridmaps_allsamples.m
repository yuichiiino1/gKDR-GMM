function draw_connection_gridmaps_allsamples()

global samples common_data_folder connection_strength_folder globalvalidweightarray  ...
    save_each_sample_fig save_each_sample_pdf globalNames globalNames1 conNames dirtext ...
    multiconmatrix FDRthre2 show_nonFDR partition_number offsets Ks draw_box_when_connected Mg show_gapchem_consistency_dotplot ...
    save_each_sample_tif save_all_valid_weights_draw_single_pdf embed_size embed_step embed_width link models_folder kGMM inconsistent

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
draw_box_when_connected = false;
show_gapchem_consistency_dotplot = false;
save_each_sample_tif = false;
save_all_valid_weights_draw_single_pdf = true;

draw_pcolor_figure = false;
draw_brickwork_figure = false;
draw_lagweightsplot_figure = true;
inconsistent = false;

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

%addpath(gKDR_codes_folder);
addpath(utilities_folder);

globalNames = table2cell(readtable(fullfile(metadata_folder,'globalNames.csv'),'ReadVariableNames',false));
globalNames1 = [{'salt'}; globalNames];
Mg = length(globalNames);
disp(['Mg=',num2str(Mg)])
multiconmatrix = table2array(readtable(fullfile(metadata_folder, 'multiconmatrix.csv'),'ReadVariableNames',false));
conNames = table2cell(readtable(fullfile(metadata_folder, 'conneurons.csv'),'ReadVariableNames',false));

load(fullfile(connection_strength_folder, 'process_weight_allnetwork_core_phys_offsets.mat'),'globalvalidweightarray')

if draw_pcolor_figure
    savefilenameheader = 'allsamples_pcolor';
    draw_pcolor(savefilenameheader)  % draw cell interaction matrix showing p-values and sign of interaction
end

if draw_brickwork_figure
    if inconsistent
        savefilenameheader2 = 'allsamples_brickwork_incon';
    else
        savefilenameheader2 = 'allsamples_brickwork';
    end

    draw_brickwork(savefilenameheader2)  % draw cell interaction matrix with grid representation of each condition
end

if draw_lagweightsplot_figure
    savefilenameheader3 = ['allsamples_lagweights'];
    draw_lagweightsplot(savefilenameheader3)
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%     Function  draw_pcolor     %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_pcolor(savefilenameheader)

global samples common_data_folder connection_strength_folder globalvalidweightarray  ...
    save_each_sample_fig save_each_sample_pdf globalNames globalNames1 conNames dirtext ...
    multiconmatrix FDRthre2 show_nonFDR partition_number offsets Ks draw_box_when_connected Mg show_gapchem_consistency_dotplot ...
    save_each_sample_tif save_all_valid_weights_draw_single_pdf embed_size embed_step embed_width link models_folder kGMM

%%%%%%%%%   Charts for all samples  %%%%%%%

savefilenamesingle = fullfile(connection_strength_folder, [savefilenameheader '_single.ps']);

if exist(savefilenamesingle, 'file')
    delete(savefilenamesingle);
end

fig = figure;
fig.PaperType       = 'a4';
fig.PaperOrientation = 'portrait';
fig.Position        = [10,10,1000,1000];
left_m = 0.1;
bot_m = 0.05;
right_m = 0.05;
top_m = 0.15;
ver_r = 1.1;
col_r = 1.1;

[validweightarrayall , p_all_mat, FDRallsig, validsource, validtarget] = select_validneurons();

% connectome
cg = NaN(1,Mg+1);  % globalNames1 -> conNames
for i = 1:Mg
    found = find(strcmp(conNames, globalNames{i}));
    assert(length(found)==1)
    cg(i) = found;
end

for row = 1:Mg+1
    for col = 1:Mg
        connected = false;
        if row-1 ~= col
            if row==1
                connected = true;
            elseif multiconmatrix(cg(row-1),cg(col))==1
                connected = true;
            end
        end
        if (~connected && strcmp(dirtext, '_direct_link')) || (~show_nonFDR && ~FDRallsig(row,col))
            p_all_mat(row,col) = NaN;
        end
    end
end


%%%%%%%%  The large matrix is devided into 4 panels　%%%%%%%%%%
colnpages = 2;
rownpages = 2;

for rowpagei = 1:rownpages

    rowfrom = ceil(length(validsource)/rownpages*(rowpagei-1))+1;
    rowto = ceil(length(validsource)/rownpages*(rowpagei));
    disp(['rowpage' num2str(rowpagei) ': from ' num2str(rowfrom) ', to ' num2str(rowto)])

    for colpagei = 1:colnpages

        colfrom = ceil(length(validtarget)/colnpages*(colpagei-1))+1;
        colto = ceil(length(validtarget)/colnpages*(colpagei));
        disp(['colpage' num2str(colpagei) ': from ' num2str(colfrom) ', to ' num2str(colto)])

        ncol = colto-colfrom+1;
        nrow = rowto-rowfrom+1;
        figcount = 0;
        ncolor = 33;
        disp(['drawing ' num2str(nrow) ' rows.'])

        logpvaluemat = NaN(nrow,ncol);
        for rowi = 1:nrow
            for coli = 1:ncol
                figcount = figcount + 1;
                row = validsource(rowfrom-1+rowi);
                col = validtarget(colfrom-1+coli);

                logpvalue = -log10(p_all_mat(row,col))/15;
                if logpvalue > 1
                    logpvalue = 1;
                end
                if nanmean(reshape(validweightarrayall(row,col, :, :),1,[])) < 0
                    logpvalue = - logpvalue;
                end
                logpvaluemat(rowi,coli) = logpvalue;

            end
        end
        ax = axes('Position',[left_m, bot_m, 1-right_m-left_m, 1-top_m-bot_m]);
        image(ax,logpvaluemat*ncolor+ncolor+2);
        colormap bluewhiteredblack(33)

        ax.XAxisLocation='top';
        xticks(1:ncol)
        yticks(1:nrow)
        xticklabels(globalNames(validtarget(colfrom:colto)))
        yticklabels(globalNames1(validsource(rowfrom:rowto)))
        ax.XTickLabelRotation = 90;
        ax.Box = 'on';

        for i = 1:ncol-1
            xline(i+0.5, 'Color', [0.5,0.5,0.5])
        end
        for i = 1:nrow-1
            yline(i+0.5, 'Color', [0.5,0.5,0.5])
        end
        c = colorbar;
        c.Ticks = [1.5,2.5:11:69];
        c.TickLabels = [{'na'}; arrayfun(@(x) {num2str(x)}, (-15:5:15)')];
        c.Label.String = '+/-log_{10}(1/p)';
        
        ax3 = annotation('textarrow','String','From', 'HeadStyle','none','LineStyle','none','FontSize',20,'TextRotation',90,'HorizontalAlignment','center','VerticalAlignment','middle');
        ax3.Position = [0,0.4,0.05,0.3];
        ax4 = annotation('textbox','String','To', 'FitBoxToText','on','LineStyle','none','FontSize',20,'HorizontalAlignment','center','VerticalAlignment','middle');
        ax4.Position = [0.35,0.9,0.3,0.05];
        
        disp('Saving all_valid_weights figure...')
        orient(fig,'tall')
        if save_each_sample_tif
            saveas(gcf, fullfile(connection_strength_folder, [savefilenameheader '_' num2str(rowpagei) '-' num2str(colpagei) '.tif']))
        end
        if save_each_sample_fig
            saveas(gcf, fullfile(connection_strength_folder, [savefilenameheader '_' num2str(rowpagei) '-' num2str(colpagei) '.fig']))
        end
        if save_each_sample_pdf
            saveas(gcf, fullfile(connection_strength_folder, [savefilenameheader '_' num2str(rowpagei) '-' num2str(colpagei) '.pdf']))
        end

        print(fig, savefilenamesingle,'-dpsc','-painters','-append');
        clf(fig)

    end
end

disp('Done.')

if ~strcmp(computer, 'MACI64')
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
    save_each_sample_fig save_each_sample_pdf globalNames globalNames1 conNames dirtext ...
    multiconmatrix FDRthre2 show_nonFDR partition_number offsets Ks draw_box_when_connected Mg show_gapchem_consistency_dotplot ...
    save_each_sample_tif save_all_valid_weights_draw_single_pdf embed_size embed_step embed_width link models_folder kGMM inconsistent

%%%%%%%%%   Charts for all samples  %%%%%%%

savefilenamesingle = fullfile(connection_strength_folder, [savefilenameheader '_single.ps']);
if exist(savefilenamesingle, 'file')
    delete(savefilenamesingle);
end

[validweightarrayall , p_all_mat, FDRallsig, validsource, validtarget] = select_validneurons();
        
%%%%%%%%  The large matrix is devided into 4 panels　%%%%%%%%%%

%figure('Position',[10,10,1000,1000])
fig = figure;
fig.PaperType       = 'a4';
fig.PaperOrientation = 'portrait';
fig.Position        = [10,10,1000,1000];
left_m = 0.1;
bot_m = 0.05;
right_m = 0.05;
top_m = 0.15;
ver_r = 1.1;
col_r = 1.1;

colnpages = 2;
rownpages = 2;

for rowpagei = 1:rownpages

    rowfrom = ceil(length(validsource)/rownpages*(rowpagei-1))+1;
    rowto = ceil(length(validsource)/rownpages*(rowpagei));
    disp(['rowpage' num2str(rowpagei) ': from ' num2str(rowfrom) ', to ' num2str(rowto)])

    for colpagei = 1:colnpages

        colfrom = ceil(length(validtarget)/colnpages*(colpagei-1))+1;
        colto = ceil(length(validtarget)/colnpages*(colpagei));
        disp(['colpage' num2str(colpagei) ': from ' num2str(colfrom) ', to ' num2str(colto)])


        %ncol = length(validtarget);
        %nrow = length(validsource);
        ncol = colto-colfrom+1;
        nrow = rowto-rowfrom+1;
        figcount = 0;
        ncolor = 32;
        disp(['drawing ' num2str(nrow) ' rows.'])
        for rowi = 1:nrow
            for coli = 1:ncol
                figcount = figcount + 1;
                row = validsource(rowfrom-1+rowi);
                col = validtarget(colfrom-1+coli);

                %rowi = floor((figcount-1)/ncol)+1;
                %coli = mod(figcount-1, ncol)+1;
                ax = axes('Position',...
                    [(1-right_m-left_m)*(mod(figcount-1,ncol))/ncol + left_m ,...
                    (1-top_m-bot_m)*(1-ceil(figcount/ncol)/(nrow)) + bot_m ,...
                    (1-right_m-left_m)/(ncol*col_r ),...
                    (1-top_m-bot_m)/(nrow*ver_r)]...
                    );
    			%if ~any(isnan(globalsigweightnonself(row,col,:))) && any(globalsigweightnonself(row,col,:))
    			if FDRallsig(row,col) || inconsistent
                    weight2d = squeeze(validweightarrayall(row,col,:,:))';
                    maxweight = nanmax(nanmax(abs(weight2d)));
                    weight2dnorm = weight2d/maxweight;

                    image(ax,weight2dnorm*ncolor+ncolor+2)

                    colormap bluewhiteredblack(32)
                end

                xticks([])
                yticks([])
                xticklabels([])
                yticklabels([])
                ax.Box = 'on';
                %ax = gca;
                %ax.Box = 'on';
                if rowi == 1
                    ax.Title.String = [' ' globalNames{col}]; %title(targetcellnames{col})
                    ax.Title.Rotation = 90;
                    ax.Title.VerticalAlignment = 'middle';
                    ax.Title.HorizontalAlignment = 'left';
                end
                if coli == 1
                    ylabel([globalNames1{row} '  '])
                    hYLabel = get(gca,'YLabel');
                    set(hYLabel,'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right')
                end
            end
        end
        %ax2 = annotation('textbox','String','All samples','FitBoxToText','on','LineStyle','none','FontSize',26,'HorizontalAlignment','center','VerticalAlignment','middle');
        %ax2.Position = [0.35,0.95,0.3,0.05];
        ax3 = annotation('textarrow','String','From', 'HeadStyle','none','LineStyle','none','FontSize',20,'TextRotation',90,'HorizontalAlignment','center','VerticalAlignment','middle');
        ax3.Position = [0,0.4,0.05,0.3];
        ax4 = annotation('textbox','String','To', 'FitBoxToText','on','LineStyle','none','FontSize',20,'HorizontalAlignment','center','VerticalAlignment','middle');
        ax4.Position = [0.35,0.9,0.3,0.05];

        disp('Saving all_valid_weights figure...')
        orient(fig,'tall')
        if save_each_sample_tif
            saveas(gcf, fullfile(connection_strength_folder, [savefilenameheader '_' num2str(rowpagei) '-' num2str(colpagei) '.tif']))
        end
        if save_each_sample_fig
            saveas(gcf, fullfile(connection_strength_folder, [savefilenameheader '_' num2str(rowpagei) '-' num2str(colpagei) '.fig']))
        end
        if save_each_sample_pdf
            saveas(gcf, fullfile(connection_strength_folder, [savefilenameheader '_' num2str(rowpagei) '-' num2str(colpagei) '.pdf']))
        end

        print(fig, savefilenamesingle,'-dpsc','-painters','-append');
        clf(fig)

    end
end

disp('Done.')

if ~strcmp(computer, 'MACI64')   % Mac does not support use of ps2pdf like this
    system(['ps2pdf ' savefilename_single]);
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
    save_each_sample_fig save_each_sample_pdf globalNames globalNames1 conNames dirtext ...
    multiconmatrix FDRthre2 show_nonFDR partition_number offsets Ks draw_box_when_connected Mg show_gapchem_consistency_dotplot ...
    save_each_sample_tif save_all_valid_weights_draw_single_pdf embed_size embed_step embed_width link models_folder kGMM

%%%%%%%%%   Charts for all samples  %%%%%%%

savefilenamesingle = fullfile(connection_strength_folder, [savefilenameheader '_single.ps']);
if exist(savefilenamesingle, 'file')
    delete(savefilenamesingle);
end

%%%%%%%%%%%%  sample loop  %%%%%%%%%%%%%%%%%%%%

global_gm_grad_array = zeros(embed_size,Mg+1,Mg,partition_number,length(offsets),length(Ks),length(samples));
for samplei = 1:length(samples)
    sampleID = samples(samplei);
    disp(['< sample ' num2str(sampleID) ' >']);

    % load basic data
    load(fullfile(common_data_folder,['sample' num2str(sampleID) '_data.mat']), 'uniqNames', 'targetcellnames', 'targetcells', 'Mt'); % 'data',

    %%%%% extract from global array %%%%
    cellg = [];   %indexes of targetcells in globalNames
    for targeti = 1:Mt
        cellname = targetcellnames{targeti};
        celli = find(strcmp(globalNames, cellname));
        assert(length(celli)==1)
        cellg = [cellg, celli];
    end
    cellg1 = [1, cellg+1]; % first is salt

    %gm_grad_mean_array = zeros(embed_size,Mt+1,Mt,partition_number,length(offsets),length(Ks));
    for partitioni = 1:partition_number
        for offset = offsets
            disp(['loading gradient data for part ', num2str(partitioni), ' offset ' num2str(offset)])
            for Ki = 1:length(Ks)
                K=Ks(Ki);
                % load gradient data
                modelfolder = fullfile(models_folder, ['model_' link '_k' num2str(embed_width) '_tau' num2str(embed_step) '_kGMM' num2str(kGMM) '_cross' num2str(partition_number)]);
                modelfileheader = fullfile(modelfolder, ['sample' num2str(sampleID) '_K' num2str(K) '_offset' num2str(offset) '_part' num2str(partitioni)]);
                load([modelfileheader '_gradient.mat'], 'gm_grad_cell', 'gm_grad_array','gm_grad_mean','gm_grad_mean_cell','uniqNames','targetcells','selicell');
                global_gm_grad_array(:,cellg1,cellg,partitioni,offset+1,Ki,samplei) = gm_grad_mean;
            end
        end
    end
    global_gm_grad_array(:,1,:,:,:,:) = - global_gm_grad_array(:,1,:,:,:,:); % salt profile is flipped
end
global_gm_grad_mean_array = squeeze(nanmean(global_gm_grad_array, [5,6])); %global_gm_grad_array = zeros(embed_size,Mg+1,Mg,partition_number,length(offsets),length(Ks),length(samples));
% -> embed_size,Mg+1,Mg,partition_number,length(samples)

[validweightarrayall , p_all_mat, FDRallsig, validsource, validtarget] = select_validneurons();
        
%%%%%%%%  The large matrix is devided into 4 panels　%%%%%%%%%%

%figure('Position',[10,10,1000,1000])
fig = figure;
fig.PaperType       = 'a4';
fig.PaperOrientation = 'portrait';
fig.Position        = [10,10,1000,1000];
left_m = 0.1;
bot_m = 0.05;
right_m = 0.05;
top_m = 0.15;
ver_r = 1.1;
col_r = 1.1;

colnpages = 2;
rownpages = 2;

for rowpagei = 1:rownpages

    rowfrom = ceil(length(validsource)/rownpages*(rowpagei-1))+1;
    rowto = ceil(length(validsource)/rownpages*(rowpagei));
    disp(['rowpage' num2str(rowpagei) ': from ' num2str(rowfrom) ', to ' num2str(rowto)])

    for colpagei = 1:colnpages

        colfrom = ceil(length(validtarget)/colnpages*(colpagei-1))+1;
        colto = ceil(length(validtarget)/colnpages*(colpagei));
        disp(['colpage' num2str(colpagei) ': from ' num2str(colfrom) ', to ' num2str(colto)])

        ncol = colto-colfrom+1;
        nrow = rowto-rowfrom+1;
        figcount = 0;
        disp(['drawing ' num2str(nrow) ' rows.'])
        for rowi = 1:nrow
            for coli = 1:ncol
                figcount = figcount + 1;
                row = validsource(rowfrom-1+rowi);
                col = validtarget(colfrom-1+coli);

                %rowi = floor((figcount-1)/ncol)+1;
                %coli = mod(figcount-1, ncol)+1;
                ax = axes('Position',...
                    [(1-right_m-left_m)*(mod(figcount-1,ncol))/ncol + left_m ,...
                    (1-top_m-bot_m)*(1-ceil(figcount/ncol)/(nrow)) + bot_m ,...
                    (1-right_m-left_m)/(ncol*col_r ),...
                    (1-top_m-bot_m)/(nrow*ver_r)]...
                    );

    			%if ~any(isnan(globalsigweightnonself(row,col,:))) && any(globalsigweightnonself(row,col,:))
    			if FDRallsig(row,col) || show_nonFDR

                    % validweightarrayall = squeeze(nanmean(globalvalidweightarray, [4,5])); % (Mg+1,Mg,partition_number,length(offsets),length(Ks),length(samples))
                    % -> (Mg+1,Mg, partition_number, length(samples)
                    weight2d = reshape(validweightarrayall(row,col,:,:),  partition_number, length(samples));  % size  3  24
                    % global_gm_grad_mean_array:  embed_size,Mg+1,Mg,partition_number,length(samples)
                    gm_grad_ii = global_gm_grad_mean_array(:,row,col,:,:); % size   30     1     1    3  24
                    gm_grad_ii(:,isnan(weight2d)) = NaN;
                    grad_embed = squeeze(nanmean(gm_grad_ii, 2:5)); % column vector
                    grad_embed_sem = nanstd(gm_grad_ii,0, 2:5); %./sqrt(sum(~isnan(gm_grad_ii), 2:5));

                    a = area(ax,[grad_embed-grad_embed_sem, grad_embed_sem*2]);
                    a(1).FaceColor = 'None';
                    a(1).EdgeColor = 'None';
                    a(2).FaceColor = [0, 0, 1];
                    a(2).FaceAlpha = 0.1;
                    a(2).EdgeColor = 'None';
                    hold on
                    plot(ax,grad_embed, 'k', 'Linewidth', 0.3)
                    yline(ax,0, 'Linewidth', 0.3)

                end

                xticks([])
                yticks([])
                xticklabels([])
                yticklabels([])
                ax.Box = 'on';
                if rowi == 1
                    ax.Title.String = [' ' globalNames{col}]; %title(targetcellnames{col})
                    ax.Title.Rotation = 90;
                    ax.Title.VerticalAlignment = 'middle';
                    ax.Title.HorizontalAlignment = 'left';
                end
                if coli == 1
                    ylabel([globalNames1{row} '  '])
                    hYLabel = get(gca,'YLabel');
                    set(hYLabel,'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right')
                end
            end
        end
        %ax2 = annotation('textbox','String','All samples','FitBoxToText','on','LineStyle','none','FontSize',26,'HorizontalAlignment','center','VerticalAlignment','middle');
        %ax2.Position = [0.35,0.95,0.3,0.05];
        ax3 = annotation('textarrow','String','From', 'HeadStyle','none','LineStyle','none','FontSize',20,'TextRotation',90,'HorizontalAlignment','center','VerticalAlignment','middle');
        ax3.Position = [0,0.4,0.05,0.3];
        ax4 = annotation('textbox','String','To', 'FitBoxToText','on','LineStyle','none','FontSize',20,'HorizontalAlignment','center','VerticalAlignment','middle');
        ax4.Position = [0.35,0.9,0.3,0.05];

        disp('Saving all_valid_weights figure...')
        orient(fig,'tall')
        if save_each_sample_tif
            saveas(gcf, fullfile(connection_strength_folder, [savefilenameheader '_' num2str(rowpagei) '-' num2str(colpagei) '.tif']))
        end
        if save_each_sample_fig
            saveas(gcf, fullfile(connection_strength_folder, [savefilenameheader '_' num2str(rowpagei) '-' num2str(colpagei) '.fig']))
        end
        if save_each_sample_pdf
            saveas(gcf, fullfile(connection_strength_folder, [savefilenameheader '_' num2str(rowpagei) '-' num2str(colpagei) '.pdf']))
        end

        print(fig, savefilenamesingle,'-dpsc','-painters','-append');
        clf(fig)

    end
end

if ~strcmp(computer, 'MACI64')   % Mac does not support use of ps2pdf like this
    system(['ps2pdf ' savefilename_single]);
    pdffilename = [savefilenameheader '_single.pdf'];
    disp(['mv ' pdffilename ' ' connection_strength_folder])
    system(['mv ' pdffilename ' ' connection_strength_folder]);
end

disp('Done.')

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%     Function  select_validneurons     %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [validweightarrayall , p_all_mat, FDRallsig, validsource, validtarget] = select_validneurons()

global samples common_data_folder connection_strength_folder globalvalidweightarray  ...
    save_each_sample_fig save_each_sample_pdf globalNames conNames dirtext ...
    multiconmatrix FDRthre2 show_nonFDR partition_number offsets Ks draw_box_when_connected Mg show_gapchem_consistency_dotplot ...
    save_each_sample_tif save_all_valid_weights_draw_single_pdf embed_size embed_step embed_width link models_folder kGMM

% difference from eachsample

validweightarrayall = squeeze(nanmean(globalvalidweightarray, [4,5])); % (Mg+1,Mg,partition_number,length(offsets),length(Ks),sample) -> (Mg+1,Mg, partition_number, samples)

p_all_mat = NaN(Mg+1,Mg);

for row = 1:Mg+1
    for col = 1:Mg
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
    end
end

FDRallmat = reshape(mafdr(reshape(p_all_mat,[],1), 'BHFDR', true),Mg+1,Mg);
FDRallsig = FDRallmat <= FDRthre2;   % NaNはfalseになる。
for i = 1:Mg
    FDRallsig(i+1,i) = false;
end


validtargetbool = any(FDRallsig, 1);
validsourcebool = any(FDRallsig, 2)';
validtarget = 1:Mg;
validtarget = validtarget(validtargetbool);
validsource = 1:Mg+1;
validsource = validsource(validsourcebool);

disp(['number of validtarget is ' num2str(length(validtarget))])
disp(['number of validsource is ' num2str(length(validsource))])

end



