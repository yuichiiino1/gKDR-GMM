function gm_xi = conditional_GMM(gm, xi)
global mu_ks_xi sigma_ks_xi pi_ks_xi
mu_ks = gm.mu;
sigma_ks = gm.Sigma;
pi_ks = gm.ComponentProportion;
kGMM = size(mu_ks,1);
K = size(mu_ks,2)-1;

mu_ks_xi = []; % mu, sigma, pi of each Gaussian conditional on x=xi
sigma_ks_xi = [];
pi_ks_xi = [];

for k = 1:kGMM  % determine conditional distribution of each Gaussian
    mu_k_x = mu_ks(k,1:end-1);
    mu_k_y = mu_ks(k,end);
    sigma_k_yy = sigma_ks(end,end,k);
    sigma_k_yx = sigma_ks(end,1:end-1,k);
    sigma_k_xy = sigma_ks(1:end-1,end,k);
    sigma_k_xx = sigma_ks(1:end-1,1:end-1,k);
    lambda_k_xx = inv(sigma_k_xx);   %%%%%%%%%%%%%%%%%%%%
    mu_k_xi = mu_k_y + sigma_k_yx * lambda_k_xx * (xi - mu_k_x)';
    sigma_k_xi = sigma_k_yy - sigma_k_yx * lambda_k_xx * sigma_k_xy;
	pi_k_xi = pi_ks(k) * mvnpdf(xi, mu_k_x, sigma_k_xx);
	%if mvnpdf(zeros(1,K), mu_k_x, sigma_k_xx)<0.001
	%	disp(mu_k_x)
	%	disp(sigma_k_xx)
	%	disp(' ')
	%end
    mu_ks_xi = cat(1, mu_ks_xi, mu_k_xi);
    sigma_ks_xi = cat(3, sigma_ks_xi, sigma_k_xi);
    pi_ks_xi = cat(1, pi_ks_xi, pi_k_xi);
end

pi_ks_xi = pi_ks_xi./sum(pi_ks_xi); % Sum(P(y|xi))=1

validpeaks = pi_ks_xi>0;  % sometimes pi=0, which causes an error.
if ~any(validpeaks)
    %disp('Error: all pis are zero')
    gm_xi = gmdistribution(zeros(1,1),zeros(1,1,1),ones(1,1)); % force reset to near zero
else
    gm_xi = gmdistribution(mu_ks_xi(validpeaks,:),sigma_ks_xi(:,:,validpeaks),pi_ks_xi(validpeaks));
end

end