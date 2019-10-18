function [U,V,L,Amean] = Compute_KL(A,wp,w,N)
Ns = size(A,2);
% [xr, w] = lgwt(Ns, -1, 1);  w = w/2;
Amean = A*w;
A = A - repmat(Amean,1,Ns);
C = A'*diag(wp)*A; C= (C+C')/2;
[V,L] = eig(C*diag(w));
L     = diag(L);
[L, I] = sort(L,'descend');
V = V(:,I(1:N));
% Lam= L;
L= L(1:N);
V = V*diag(1./sqrt(w'*(V.^2)));
V = V*diag(sqrt(L));
C = V'*diag(w)*V;
U = A*diag(w)*V/C;
% L = Lam;