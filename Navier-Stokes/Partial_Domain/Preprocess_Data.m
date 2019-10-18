load NS_Data.mat
[D, Ns, kmax] = size(u);
dx = x(2)-x(1);
dy = y(2)-y(1);
dv = dx*dy;
dt = 0.02;
t  = 0:dt:dt*(kmax-1);
%% ----------------- Selecting a Subset of the Data -----------------
I_ROI = find((Y(:)-3).^2 + X(:).^2<=1.6^2 );

u = u(I_ROI,:,:);
v = v(I_ROI,:,:);
D = length(I_ROI);
w = dv*ones(2*D,1);
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
