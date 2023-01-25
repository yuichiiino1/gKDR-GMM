function combine_simulation_figures_indirect()

% arrange heatmap and correlation matrix for real and two repeats of
% simulation

% Requirements
% simulation results figures as .fig files

samples = 1:24;

embed_width = 30;
embed_step = 10;
Ks = 3:5;
kGMM = 2;
link = 'indirect';
project_folder = '/home/iino/gKDR_GMM_50mM/Linked_gKDR_freerun/publish_version';
%project_folder = '/home/iino/gKDR_GMM_50mM/Linked_gKDR_freerun/ver98_HLong';
simulation_results_folder = fullfile(project_folder, 'simulation_results');
%simulation_results_folder = fullfile(project_folder, 'figures_');

close all

%%%%% Prepair figure area %%%%%

axesFontSize = 11;
fontName = 'arial';
fontSizeMP = 1;

fig = figure;
fig.PaperType       = 'a4';
fig.PaperOrientation = 'landscape';
fig.PaperUnits      = 'centimeters';
fig.PaperPosition   = [0,0,29.7,21];
fig.Units           = 'centimeters';
fig.Position        = [0,0,29.7,21];
fig.Color           = 'w';
fig.InvertHardcopy  = 'off';

set(fig,'defaultAxesFontSize',axesFontSize);
set(fig,'defaultAxesXColor','k'); % factory is [0.15,0.15,0.15]
set(fig,'defaultAxesYColor','k');
set(fig,'defaultAxesZColor','k');
set(fig,'defaultAxesUnits','centimeters');
set(fig,'defaultAxesTickDir','out');
set(fig,'defaultAxesBox','off');
set(fig,'defaultAxesFontName',fontName);

set(fig,'defaultTextFontName',fontName);
set(fig,'defaultTextFontSize',axesFontSize);

set(fig,'defaultLegendFontName',fontName);
set(fig,'defaultLegendFontSize',axesFontSize);

set(fig,'defaultAxesLabelFontSizeMultiplier',fontSizeMP);

wn = 4;
hn = 2;
wstart1 = 7.5;   % For position of plot box. In cm from left. A4 size is 21cm.
wspace1 = 11.5;
hstart1 = 11;   % In cm from bottom. Paper size is 29cm.
%hspace = 10;
wsize1 = 9;
hsize1 = 7;

hstart2 = 2;
%hspace2 = 9.5;
wstart2 = 2; %3.5
wspace2 = 10; %11.5
wsize2 = 7.5;
hsize2 = 7.5;

ax1 = gobjects(1,2);
ax2 = gobjects(1,3);

for K = Ks
		
	disp(['< K ' num2str(K) ' >'])
    simulation_result_subfolder = fullfile(simulation_results_folder, [link '_k' num2str(embed_width) '_tau' num2str(embed_step) '_K' num2str(K) '_kGMM' num2str(kGMM)]);
    %simulation_result_subfolder = [simulation_results_folder, link '_k' num2str(embed_width) '_tau' num2str(embed_step) '_K' num2str(K) '_kGMM' num2str(kGMM)];
    
    savefilename = fullfile(simulation_result_subfolder, [link '_k' num2str(embed_width) '_tau' num2str(embed_step) '_K' num2str(K) '_kGMM' num2str(kGMM) 'samples' num2str(samples(1)) '-' num2str(samples(end)) '.ps']);
    if exist(savefilename, 'file')
        delete(savefilename);
    end
    
    %%%%% put figures for each sample  %%%%%
    
    for sampleID = samples
        
        disp(['< sample ' num2str(sampleID) ' >'])
        %simulationresultfolder = ['../ver98_HLong/figures_indirect_k' num2str(embed_width) '_tau' num2str(embed_step) '_K' num2str(K) '_kGMM' num2str(kGMM) '/sample' num2str(sampleID)];
        %simulationresultfileheader = [simulationresultfolder '/sample' num2str(sampleID)];
        %simulation_result_subfolder2 = fullfile(simulation_result_subfolder, ['sample' num2str(sampleID)]);
        simulationresultfileheader = fullfile(simulation_result_subfolder, ['sample' num2str(sampleID)]);
        
        % heatmap
        for testi = 1:2
            
            try
                heatmapfig = open([simulationresultfileheader '_freerun' num2str(testi) '_heatmap.fig']);
                figure(heatmapfig)
                tempax = gca;
                figcolor = heatmapfig.Colormap;
                ax1(testi) = copyobj(tempax,fig);
                ax1(testi).OuterPosition = [wstart1+wspace1*(testi-1)-3, hstart1, wsize1+6, hsize1];
                ax1(testi).Position = [wstart1+wspace1*(testi-1), hstart1, wsize1, hsize1];
                ax1(testi).InnerPosition = [wstart1+wspace1*(testi-1), hstart1, wsize1, hsize1];
                %ax1(testi).PositionConstraint = 'outerposition';
                ax1(testi).FontSize=5;
                ax1(testi).YTickLabel = cellfun(@(x) [x,'  '], ax1(testi).YTickLabel, 'UniformOutput', false);
                ax1(testi).Title.String = ['Repeat ' num2str(testi)];
                ax1(testi).TitleFontSizeMultiplier = 4;
                ax1(testi).Colormap = figcolor;
            catch
                disp('figure file not found')
            end
        end
        
        % corrmatrix
        
        for testi = 0:2
            
            try
                if testi == 0
                    cormatrixfig = open([simulationresultfileheader '_real_corrmatrix.fig']);
                else
                    cormatrixfig = open([simulationresultfileheader '_freerun' num2str(testi) '_corrmatrix.fig']);
                end
                figure(cormatrixfig)
                tempax = gca;
                figcolor = cormatrixfig.Colormap;
                ax2(testi+1) = copyobj(tempax,fig);
                ax2(testi+1).OuterPosition = [wstart2+wspace2*(testi)-3, hstart2-3, wsize2+6, hsize2+6];
                ax2(testi+1).Position = [wstart2+wspace2*(testi), hstart2, wsize2, hsize2];
                ax2(testi+1).PositionConstraint = 'outerposition';
                ax2(testi+1).FontSize=5;
                ax2(testi+1).XTickLabel = cellfun(@(x) [x,'  '], ax2(testi+1).XTickLabel, 'UniformOutput', false);
                ax2(testi+1).YTickLabel = cellfun(@(x) [x,'  '], ax2(testi+1).YTickLabel, 'UniformOutput', false);
                if testi == 0
                    ax2(testi+1).Title.String = ['Real'];
                else
                    ax2(testi+1).Title.String = ['Repeat ' num2str(testi)];
                end
                ax2(testi+1).TitleFontSizeMultiplier = 4;
                ax2(testi+1).Colormap = figcolor;
            catch
                disp('figure file not found')
            end
        end
        
        % grand title
        figure(fig)
        ax6 = annotation('textbox','String',['Sample' num2str(sampleID)],'FitBoxToText','on','LineStyle','none','FontSize',26,'Units','centimeters');
        ax6.Units = 'centimeters';
        ax6.Position = [14,20,6,1]; % [x_begin y_begin length height]
        
        
        orient(fig,'landscape')
        %print(fig,'-painters',savefilename,'-dpdf')
        print(fig, savefilename,'-dpsc','-painters','-append');
        clf(fig)
        
    end
end
end
