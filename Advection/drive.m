close all
clear all
clc

addpath ../Utilities
addpath ../Data
set(0,'defaulttextinterpreter','latex')
LW = 'linewidth';


global w wr
global dt t r 
%% ------------- Specify the Reduction Order ----------
r = 2;
%% --------------------- Load Data ------------
load Advection_Data
%% ------------ Physical Space Discretization ----------
D = 128;
[x,  w ] = chebpts(D,dom); x=x';
dt = t(2) - t(1);
%% --------------- Compute dot(A) --------------
Nt = length(t);
A = zeros(D,Nt,Ns);
for i=1:Ns
    A(:,:,i) = Sample(i).u(x',:);
end
clear Sample
A  = permute(A, [1 3 2]);
Amean = zeros(D,Nt);
Ap = zeros(size(A));
for k=1:Nt
    Amean(:,k) = A(:,:,k)*wr;
    Ap(:,:,k)  = A(:,:,k) - repmat(Amean(:,k),1,Ns);
end
Ap  = permute(Ap, [1 3 2]);
for i=1:Ns
    AD(:,:,i) = ComputeTimeDerivative(Ap(:,:,i));
end
AD  = permute(AD, [1 3 2]);
Ap  = permute(Ap, [1 3 2]);


%% --------------- Compute Initial Condition for DB --------------
ns = 2;
A0 = A(:,:,ns) - repmat(Amean(:,ns),1,Ns);
C = A0'*diag(w)*A0; C= (C+C')/2;
[V,L] = eig(C*diag(wr));
L     = diag(L);
[L, I] = sort(L,'descend');
V = V(:,I(1:r));
L= L(1:r);
V = V*diag(1./sqrt(wr'*(V.^2)));
V = V*diag(sqrt(L));
C = V'*diag(wr)*V;
U = A0*diag(wr)*V/C;


for n=1:Nt-ns
    
    [k1U, k1V] = rhs_DB(U,          V,          AD(:,:,n+ns));
    [k2U, k2V] = rhs_DB(U+dt*k1U/2, V+dt*k1V/2, AD(:,:,n+ns));
    [k3U, k3V] = rhs_DB(U+dt*k2U/2, V+dt*k2V/2, AD(:,:,n+ns));
    [k4U, k4V] = rhs_DB(U+dt*k3U,   V+dt*k3V,   AD(:,:,n+ns));
     
    U   = U   + dt*(k1U + 2*k2U + 2*k3U + k4U)/6;
    V   = V   + dt*(k1V + 2*k2V + 2*k3V + k4V)/6;
  
    for i=1:r
        for j=1:i-1
            U(:,i)  = U(:,i) - (U(:,i)'*diag(w)*U(:,j))*U(:,j);
        end
        U(:,i) = U(:,i)/sqrt(U(:,i)'*diag(w)*U(:,i));
    end
  
    %% -------------------- Compute the Reduction Error -----------
    Error = U*V'-Ap(:,:,n+ns);
    Error = Error'*diag(w)*Error;
    E(n)  = sum(eig(Error*diag(wr)));
    %% -------- Rank the reduction & and rotate to bi-orthonormal form -----------
    C  = V'*diag(wr)*V;
    [R,L] = eig(C);
    [L, I] = sort(diag(L),'descend');
    Lambda(:,n) = L;
    R= R(:,I);
    
    VBO = V*R*sqrt(diag(1./L));
    UBO = U*R;
    subplot(1,2,1)
    plot(x,UBO,LW,2)
    xlabel('$t$')
    ylabel('$U(t)$')
    set(gca,'FontSize',15);
    subplot(1,2,2)
    plot(VBO(:,1),VBO(:,2),'o');
    axis([ -2 2 -2 2])
    axis equal
    xlabel('$y_1$')
    ylabel('$y_2$')
    set(gca,'FontSize',15);
    drawnow
    
end
figure
plot(t(ns:n+ns-1),Lambda(:,1:n),LW,1.2)
xlabel('$t$')
ylabel('$\Lambda(t)$')
set(gca,'FontSize',15);


