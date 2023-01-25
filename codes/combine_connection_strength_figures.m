
function combine_connection_strength_figures(UGE_TASK_ID_text) % typeNo (1:3)

% read .fig files and arrange in one .ps file.
% .ps files can be read and converted to .pdf by Acrobat etc later.
% originally pdf_output_ABnetwork_samples_flipped.m

nsamples = 24;
link = 'indirect';
typetexts = {'AB_network','SAB_network', 'NSAB_network'};
typeno = str2num(UGE_TASK_ID_text);
typetext = typetexts{typeno};

project_folder = '/home/iino/gKDR_GMM_50mM/Linked_gKDR_freerun/publish_version';
results_folder = fullfile(project_folder, 'connection_strength', link);
%results_folder = '/Volumes/Transcend J/_4DImaging_Data_Analysis/gKDR-GMM/gKDR_GMM_50mM/Linked_gKDR_freerun/publish_version/connection_strength/'; %'connection_strength_flip';

outputfolder = results_folder; %'pdf_networks_flipped/';
if ~exist(outputfolder, 'dir')
    mkdir(outputfolder)
end


Ks = 3:5;

for Ki = 1:length(Ks)
    K = Ks(Ki);
    disp(['< K ' num2str(K) ' >'])
    
    axesFontSize = 11;
    fontName = 'arial';
    fontSizeMP = 1;
    
    % Paper size for arrangement
    % Depends on own monitor. This is for secondary monitor.
    fig = figure;
    fig.PaperType       = 'a4';
    fig.PaperOrientation = 'portrait';
    fig.PaperUnits      = 'centimeters';
    fig.Units           = 'centimeters';
    fig.Position = [-25,40,21,29.7];
    fig.PaperPosition   = [0,0,21,29.7];
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
    
    wn = 2;
    hn = 3;
    wstart = 2;   % For position of plot box. In cm from left. A4 size is 21cm.
    wspace = 9;
    hstart = 20;   % In cm from bottom. Paper size is 29cm.
    hspace = 9;
    wsize = 7;
    hsize = 7;
    
    axbody = gobjects(hn,wn);
    
    sampleperpage = 6;
    
    filename = fullfile(outputfolder, [typetext '_K' num2str(K) '_' link '_sample1-' num2str(nsamples) '.ps']);
    if exist(filename, 'file')
        delete(filename);
    end
    
    for sampleID = 1:nsamples
        
        disp(['sample ' num2str(sampleID)])
        
        % read fig files for each sample and arrange on A4.
        
        k = mod(sampleID-1,sampleperpage)+1;
        i = floor((k-1)/wn)+1; % rows
        j = mod(k-1,wn)+1; % columns
        
        figfile = [results_folder '/' typetext '_K' num2str(K) '_sample' num2str(sampleID) '.fig'];
        if exist(figfile, 'file')
            tempfig = open(figfile);
            figure(tempfig)
            tempax = gca;
            axbody(i,j) = copyobj(tempax,fig);
            axbody(i,j).Position = [wstart+wspace*(j-1), hstart-hspace*(i-1), wsize, hsize];
            axbody(i,j).Title.String = ['sample ' num2str(sampleID)];
            axbody(i,j).TitleFontSizeMultiplier = 2;
            close(tempfig)
        end

        % output to file. Multipage ps.
        if k==sampleperpage || sampleID == nsamples
            orient(fig,'tall')
            print(fig, filename,'-dpsc','-painters','-append');
            clf(fig)
        end
        
    end
    
end
