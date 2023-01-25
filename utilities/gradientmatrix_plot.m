function gradientmatrix_plot(gradmatrix, figure_title, suffix, exptitle2, cell_num, labnum)

colormap bluewhiteredblack(32) %jet(64) %gray
gradmatscaled = gradmatrix*5;
gradmatscaled(gradmatscaled < -1) = -1;
gradmatscaled(gradmatscaled > 1) = 1;
image(gradmatscaled*32+34);
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
saveas(gcf, [exptitle2 '_gradmatrix'  '.tif']);
saveas(gcf, [exptitle2 '_gradmatrix'  '.fig']);

end
