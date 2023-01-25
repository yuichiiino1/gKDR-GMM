function lagmatrix_plot(rlagmatrix, figure_title, suffix, exptitle2, cell_num, labnum, maxlag)

colormap greenwhitemazentablack(32) %jet(64) %gray
rlagmatrix1 = rlagmatrix/maxlag*32+32+2;
rlagmatrix1(isnan(rlagmatrix))=1;
image(rlagmatrix1);
% title(['Lag matrix ' figure_title],'Interpreter','none');
if ~strcmp(suffix, 'real')
title(['Spearman''s p = ' sprintf('%.2g',figure_title)],'Interpreter','none');
end

h = gca;
set(h,'ytick',(1:cell_num)');
set(h,'yticklabel',nameshifter(labnum));
set(h,'xtick',(1:cell_num)');
set(h,'xticklabel',nameshifter(labnum));
xtickangle(90);
set(h,'fontsize',6);
set(h,'Ticklength',[0 0]);
set(h,'TitleFontSizeMultiplier', 4);
colorbar('Ticks',[0,16,32,48,64],...
    'TickLabels',{num2str(-maxlag*5),num2str(-floor(maxlag*5/2)),'0',num2str(floor(maxlag*5/2)),num2str(maxlag*5)})
%saveas(gcf, ['correlationmatrix_' suffix '_sample' num2str(samplei) '.fig']);
saveas(gcf, [exptitle2 '_' suffix '_lagmatrix'  '.tif']);
saveas(gcf, [exptitle2 '_' suffix '_lagmatrix'  '.fig']);

end
