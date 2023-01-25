function gm_grad = GMMgradient(source_train, gm, B) % gm_grad = M x Ntrain matrix

  mu_ks = gm.mu;
  sigma_ks = gm.Sigma;
  pi_ks = gm.ComponentProportion;
  kGMM = size(mu_ks,1);
  K = size(mu_ks,2)-1;
  
  dydx_ks = [];
  for k = 1:kGMM
    sigma_k_yx = sigma_ks(end,1:end-1,k);
    sigma_k_xx = sigma_ks(1:end-1,1:end-1,k);
    lambda_k_xx = inv(sigma_k_xx);
    dydx_k = sigma_k_yx * lambda_k_xx; % ���x�N�g��
    dydx_ks = [dydx_ks; dydx_k]; % �s�Fcomponent k�A��Fx�̎���
  end
  
  gm_grad = zeros(size(source_train,1),size(source_train,2));
  for i = 1:size(source_train,1)
    si = source_train(i,:);
    %gm_xi = conditional_GMM(gm, si*B);
    xi = si*B;
    pi_ks_xi = [];
    for k = 1:kGMM  % determine weight of each Gaussian at xi
      mu_k_x = mu_ks(k,1:end-1);
      sigma_k_xx = sigma_ks(1:end-1,1:end-1,k);
      pi_k_xi = pi_ks(k) * mvnpdf(xi, mu_k_x, sigma_k_xx);
      pi_ks_xi = cat(2, pi_ks_xi, pi_k_xi);
    end
    pi_ks_xi = pi_ks_xi./sum(pi_ks_xi); % Sum(P(y|xi))=1
    
    % i�Ԗڂ̃f�[�^��gradient���v�Z
    %pi_gm_xi = gm_xi.ComponentProportion; % xi�ɂ�����gaussian�̍�����
    dydx_xi = pi_ks_xi * dydx_ks; % ������̏d�݂ő��a��x�̃T�C�Y�̉��x�N�g��
    gm_grad_si = dydx_xi * B'; % dy/ds = dy/dx * dx/ds, x= s*B���dx/ds=B'
    gm_grad(i,:) = gm_grad_si;
  end

end

