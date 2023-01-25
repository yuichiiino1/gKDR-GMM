function cormatrix_plot(rmatrix, figure_title, suffix, exptitle2, cell_num, labnum)

colormap bluewhitered(32) %jet(64) %gray
image(rmatrix*32+32);
title(['Correlation matrix ' figure_title],'Interpreter','none');
h = gca;
set(h,'ytick',(1:cell_num)');
set(h,'yticklabel',nameshifter(labnum));
set(h,'xtick',(1:cell_num)');
set(h,'xticklabel',nameshifter(labnum));
xtickangle(90);
set(h,'fontsize',6);
set(h,'Ticklength',[0 0]);
colorbar('Ticks',[0,16,32,48,64],...
    'TickLabels',{'-1','-0.5','0','0.5','1'})
%saveas(gcf, ['correlationmatrix_' suffix '_sample' num2str(samplei) '.fig']);
saveas(gcf, [exptitle2 '_' suffix '_corrmatrix'  '.tif']);
saveas(gcf, [exptitle2 '_' suffix '_corrmatrix'  '.fig']);

end
