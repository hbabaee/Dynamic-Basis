load NS_Data.mat
[D, Ns, kmax] = size(u);
dx = x(2)-x(1);
dy = y(2)-y(1);
dv = dx*dy;
dt = 0.02;
t  = 0:dt:dt*(kmax-1);
%% --------------------- Use the Full Data set ----------------------
Is = [1:imax*jmax];
w  = dv*ones(2*D,1);
%% ----------------- Selecting a Subset of the Data -----------------

% I1 = find(Y(:)>=2 & Y(:) <=2.5);
% I2 = find(Y(:)>=3 & Y(:) <=3.5);
% I = [I1; I2];

% Is = find((Y(:)-3).^2 + X(:).^2<=1.6^2 );
% jmax = length(I)/imax;
% X = reshape(X(I),jmax,imax);
% Y = reshape(Y(I),jmax,imax);

% u = u(Is,:,:);
% v = v(Is,:,:);
% D = length(Is);

%% ---------------------------------------------------
A = [u;v]; Amean = zeros(2*D,kmax);
Ap    = zeros(size(A));
for k=1:kmax
    Amean(:,k) = A(:,:,k)*wr;
    Ap(:,:,k)  = A(:,:,k) - repmat(Amean(:,k),1,Ns);
end
Ap  = permute(Ap, [1 3 2]);
for i=1:Ns
    AD(:,:,i) = ComputeTimeDerivative(Ap(:,:,i));
end
AD = permute(AD,[1 3 2]);
Ap  = permute(Ap, [1 3 2]);
