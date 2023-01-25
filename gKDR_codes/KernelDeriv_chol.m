%  function KernelDeriv_chol() by YI
%  
%  based on KernelDeriv() except using incomplete cholesky approximation to
%  solve memory problem.

%  function KernelDeriv()
%
%  Computing effective directions for regression with RKHS.
%
%  Author: Kenji Fukumizu (Institute of Statistical Mathematics)
%  Date: July 17, 2013
%  Version: 1.00
% 
%
%  Author:  Kenji Fukumizu
%  Affiliation:  The Institute of Statistical Mathematics, ROIS
%  (C) Copyright:  Kenji Fukumizu
%
%-----------------------------------------------
% KernelDeriv()
%
% Arguments
%  X:  explanatory variables (input data)
%  Y:  response variables (teaching data)
%  K:  dimension of effective subspaces
%  SGX:  bandwidth (deviation) parameter in Gaussian kernel for X
%  SGY:  bandwidth (deviation) parameter in Gaussian kernel for Y
%  EPS:  regularization coefficient
%
% Return value(s)
%  B:  orthonormal column vectors (M x K)
%  t:  value of the objective function
%  R:  estimated Mn (added by YI)
% 
%-----------------------------------------------

function [B, R, t, Kx]=KernelDeriv_chol(X,Y,K,SGX,SGY,EPS)

[N,M]=size(X);  % N: data size, M: dim of X.

tol=0.000001;   % tolerance for incomplete cholesky approximation

I=eye(N);

sx2=2*SGX*SGX;
sy2=2*SGY*SGY;

% Gram matrix of X
ab=X*X';
aa=diag(ab);
D=repmat(aa,1,N);
xx=max(D + D' - 2*ab, zeros(N,N));
Kx=exp(-xx./sx2);  

% Gram matrix of Y
ab=Y*Y';
aa=diag(ab);
D=repmat(aa,1,N);
yy=max(D + D' - 2*ab, zeros(N,N));
Ky=exp(-yy./sy2);  

% incomplete cholesky approximation of Ky
%disp(['Execute chol_inc_gauss'])
[G, Pvec] = chol_inc_gauss(Y',SGY,tol);
[a,Pvec]=sort(Pvec); 
Ry=G(Pvec,:);
r=length(Ry(1,:));
Ty=Ry'/(Kx+(N*EPS).*eye(N));
clear Ry;

%disp(['Derivative of k(X_i, x'])
% Derivative of k(X_i, x) w.r.t. x
Dx=reshape(repmat(X,N,1),N,N,M);
Xij=Dx-permute(Dx,[2 1 3]);
Xij=Xij./SGX/SGX;
H=Xij.*repmat(Kx,[1 1 M]);
clear Xij;
%disp(['Finished Derivative of k(X_i, x'])

% compute the matrix for gKDR
Hy=reshape(Ty*reshape(H,[N,N*M]), [r*N,M]);
R=Hy'*Hy;
clear Hy;

[V,L]=eig(R);
[e,idx]=sort(diag(L),'descend');
B=V(:,idx(1:K));
t=sum(e(idx(1:K)));

