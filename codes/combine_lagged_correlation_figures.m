
function combine_lagged_correlation_figures()

samples = 1:1; %24;
thre = 0.3;

project_folder = '/home/iino/gKDR_GMM_50mM/Linked_gKDR_freerun/publish_version';
lagged_correlation_folder = fullfile(project_folder, 'lagged_correlation');

%typeno = str2num(UGE_TASK_ID_text);
%typetext = typetexts{typeno};
%resultheader = ['./lagged_correlation_C_thre' num2str(thre) '_pdf/'];
%if ~exist(resultheader, 'dir')
%    mkdir(resultheader)
%end

close all

%Ks = [3,4,5];

axesFontSize = 11;
fontName = 'arial';
fontSizeMP = 1;

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

wn = 3;
hn = 4;
wstart = 2.4;  %2.6 % plot boxの位置で指定される。左から。紙の端が0でA4は29.7cm。
wspace = 6.1;  % 6.5
hstart = 22.3; %22.3;   % 下からの距離。cm。紙の端が0でA4は21cm。
hspace = 7;
wsize = 5;
hsize = 5.5;

% left side
hstart2 = 25;   % 下からの距離。cm。紙の端が0でA4は21cm。
hspace2 = 7;
wstart2 = 0.2;
wsize2 = 1.5;
hsize2 = 1;

% top
wstart3 = 4.3;   % plot boxの位置で指定される。左から。紙の端が0でA4は29.7cm。
wspace3 = 5.6;
hstart3 = 28.3; %28;
wsize3 = 6;
hsize3 = 1;

% right side
hstart4 = 22.3;   % 下からの距離。cm。紙の端が0でA4は21cm。
hspace4 = 7;
wstart4 = 19.9;
wsize4 = 0.3;
hsize4 = 5.5;


axbody = gobjects(hn,wn);
axleft = gobjects(hn,1);
axtop = gobjects(wn,1);
axright = gobjects(hn,1);

sampleperpage = 4;

outfilename = fullfile(lagged_correlation_folder, ['lagged_correlation_sample' num2str(samples(1)) '-' num2str(samples(end)) '.ps']);
if exist(outfilename, 'file')
    delete(outfilename);
end
    
    
for samplei = 1:length(samples)
    sampleID = samples(samplei);
    
    disp(['<< sample ' num2str(sampleID) ' >>'])
    
    lagged_correlation_figures_header = fullfile(lagged_correlation_folder, ['sample' num2str(sampleID)]);
        
    i = mod(sampleID-1,sampleperpage)+1;
    %i = floor((k-1)/wn)+1; % tate
    
    filefooters = {'_real_lagmatrix.fig', '_freerun1_lagmatrix.fig', '_freerun2_lagmatrix.fig'};
    
    for j = 1:3 % yoko
        
        %exptitle2 = ['./sample' num2str(sampleID)];
        figfile = [lagged_correlation_figures_header filefooters{j}];
        if exist(figfile, 'file')
            tempfig = open(figfile);
            figure(tempfig)
            tempax = gca;
            figcolor = tempfig.Colormap;
            if j<3
                axbody(i,j) = copyobj(tempax,fig);
            else
                %axbody(i,j) = copyobj(tempax,fig);
                %disp(tempfig.Children)
                tempobj = copyobj([tempfig.Children(1),tempfig.Children(2)],fig);
                axbody(i,j) = tempobj(2);
                axright(i) = tempobj(1);
                %disp(axright(i))
                %axright(i) = copyobj(tempfig.Children(1),fig);
            end
            axbody(i,j).Position = [wstart+wspace*(j-1), hstart-hspace*(i-1), wsize, hsize];
            %axbody(i,j).PositionConstraint = 'innerposition';
            axbody(i,j).FontSize=3;
            %axbody(i,j).XTickLabel = cellfun(@(x) [x,'  '], axbody(i,j).XTickLabel, 'UniformOutput', false);
            axbody(i,j).YTickLabel = cellfun(@(x) [x,'  '], axbody(i,j).YTickLabel, 'UniformOutput', false);
            %axbody(i,j).Title.String = ''; %['sample ' num2str(sampleID)];
            axbody(i,j).TitleFontSizeMultiplier = 2;
            axbody(i,j).Colormap = figcolor;
            %disp(axbody(i,j).Children)
            %for k = 6:-1:4
            %    axbody(i,j).Children(k).MarkerSize=7;
            %end
            if j==3
                axright(i).Units = 'centimeters';
                axright(i).Position = [wstart4, hstart4-hspace4*(i-1), wsize4, hsize4];
                axright(i).FontSize=5;
            end
            close(tempfig)
        end
    end
    if i==1
        toptitles = {'real', 'freerun1', 'freerun2'};
        for j=1:3
            axtop(j) = annotation('textbox','String',toptitles{j},'FitBoxToText','on','LineStyle','none','FontSize',18,'Units','centimeters');
            axtop(j).Units = 'centimeters';
            axtop(j).Position = [wstart3+(j-1)*wspace3,hstart3,wsize3,hsize3]; % [x_begin y_begin length height]
        end
    end
    
    axleft(i) = annotation('textbox','String',{'Spl';['  ' num2str(sampleID)]},'FitBoxToText','off','LineStyle','none','FontSize',14,'Units','centimeters');
    axleft(i).Units = 'centimeters';
    axleft(i).Position = [wstart2,hstart2-hspace2*(i-1),wsize2,hsize2]; % [x_begin y_begin length height]
    
    if i==sampleperpage || sampleID == samples(end)
        orient(fig,'tall')
        print(fig,outfilename,'-dpsc','-painters','-append')
        clf(fig)
    end
    
end
end

