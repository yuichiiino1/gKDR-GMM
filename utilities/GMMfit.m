function gm = GMMfit(source_train, target_train, B, kGMM)

show_gmmplot = true;
ntestlines = 20;
nlinesteps = 100;
ncolorsteps = 100;

%%%%%%%%%%%%%%%%   GMM fit   %%%%%%%%%%%%%%%%

x_train = source_train*B;
assert(size(source_train,1)~=0);

options = statset('MaxIter',1000); % 'Display','final',
gm = fitgmdist([x_train,target_train],kGMM,'Options',options,'Start','plus');

%%%%%%%%%%%%   show plot to visualize GMM distribution %%%%%%%%%%%%

if show_gmmplot && size(x_train,2)==1
    mu_ks = gm.mu;
    sigma_ks = gm.Sigma;
    %pi_ks = gm.ComponentProportion;
    
    figure;
    
    % plot data
    scatter(x_train,target_train,10,'.') % Scatter plot with points of size 10
    title('Train data')
    hold on
    
    % draw ellipsoid lines
    for k=1:kGMM
        warning('off','MATLAB:fplot:NotVectorized')
        fcontour(@(x,y)reshape(mvnpdf([x(:) y(:)],mu_ks(k,:),sigma_ks(:,:,k) ), size(x)))
        warning('on','all')
    end
    text(mu_ks(:,1), mu_ks(:,2), num2str((1:kGMM)'),'HorizontalAlignment','center','fontsize',20,'Color','k');
    %scatter(mu_ks(:,1), mu_ks(:,2), 120, '+')
    
    % draw color vertial lines showing conditional distribution
    if ntestlines > 1
        xmin = min(x_train); %,[],1);
        xmax = max(x_train); %,[],1);
        xis = xmin:(xmax-xmin)/(ntestlines-1):xmax;
        for xi = xis % horizontal position of each line
            
            gm_xi = conditional_GMM(gm, xi);
            
            ymin = min(target_train);
            ymax = max(target_train);
            ystep = (ymax-ymin)/nlinesteps;
            
            xjs = (ymin + ystep * (0.5) : ystep : ymin + ystep * (nlinesteps-0.5))';
            %z(j) = normpdf(ymin + ystep * (j-0.5), mu_k_xi, sigma_k_xi);
            z = pdf(gm_xi, xjs);
            
            cmap = colormap(jet(ncolorsteps+1));
            zmax = max(z);
            
            for j = 1:nlinesteps % line segments
                plot([xi, xi], [ymin + ystep*(j-1), ymin + ystep*(j)],...
                    'Color', cmap(round(z(j)/zmax*ncolorsteps)+1, :));
            end
            
        end
        
    end % end plot
    xlim([xmin,xmax]);
    ylim([ymin,ymax]);
    
end

end



