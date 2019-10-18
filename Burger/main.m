clc
clear all
close all

addpath ../Data
set(0,'defaulttextinterpreter','latex')
LW = 'linewidth';

RGB1 = [0 113 188 ]/norm([0 113 188 ]);
RGB2 = [216 82 24 ]/norm([216 82 24 ]);
% RGB3 = [20 82 24 ]/norm([20 82 24 ]);
% color = [55,126,184]/255;
% R     = [162, 82,  24 ]/255; %R = R/norm(R);
% B     = [0 , 113, 188]/255; %B = B/norm(B);
% c1 = R;


global  tol_PI

global DB Y u ubar du dubar ddu ddubar x wp xr wr mu C D G Nr

xr = load('xi.dat');
wr = load('wr.dat');
load('Burger_Data')



tol_PI = 1e-14;
N = 5;     
N_init = N;      
Nr = length(wr);    
M  = Nr;
L = 2*pi;      % the length of the interval [0,1]
ns = 2;


x=L*(0:Ns-1)'/Ns;   % collocation points
wp = L/Ns*ones(length(x),1);


t0 = (ns-1)*dt;


nTimeStep = ceil((tf-t0) / dt);
nStepsP = 10;	% number of time steps at which the data is saved



A = permute(A,[1 3 2]);


[U,V,Lambda,Amean] = Compute_KL(A(:,:,ns),wp,wr,N);
A0 = A(:,:,ns);


clear Lambda

%% Initialization
dubar = Amean;
dY = V;   
du = U;  

ndY = zeros(Nr, N);
ndu = zeros(Ns, N);   
ndubar = zeros(Ns, 1);

MS=zeros(Ns, Nr);
MS_pcm=zeros(Ns, Nr);


coeff_RK3 = [1/4 0 3/4];

%
tic;
n=1;

k=1;
t3=0;
while n <= nTimeStep
    
    
    t1 = dt * (n-1)+t0;
    [rhs_dubar1 rhs_du1 rhs_dY1 du  dY] = compute_rhs_do(dubar, du, dY, N, Ns, Nr, xr, wr, wp, t1, mu);
    
    t2 = dt * (n-1/2)+t0;
    dubar2 	= dubar + dt*rhs_dubar1/2.0;
    du2		= du + dt*rhs_du1/2.0;
    dY2		= dY + dt*rhs_dY1/2.0;
    [rhs_dubar2 rhs_du2 rhs_dY2 du2 dY2] = compute_rhs_do(dubar2, du2, dY2, N, Ns, Nr, xr, wr, wp, t2, mu);
    
    t3 = dt * n+t0;
    dubar3 	= dubar - dt*rhs_dubar1 + 2*dt*rhs_dubar2;
    du3		= du - dt*rhs_du1 + 2*dt*rhs_du2;
    dY3		= dY - dt*rhs_dY1 + 2*dt*rhs_dY2;
    [rhs_dubar3 rhs_du3 rhs_dY3 du3 dY3] = compute_rhs_do(dubar3, du3, dY3, N, Ns, Nr, xr, wr, wp, t3, mu);
    
    ndubar 	= dubar + dt*(rhs_dubar1+4*rhs_dubar2+rhs_dubar3)/6.0;
    ndu		= du + dt*(rhs_du1+4*rhs_du2+rhs_du3)/6.0;
    ndY		= dY + dt*(rhs_dY1+4*rhs_dY2+rhs_dY3)/6.0;
    
    t = t3;
    
    
    
    for i=1:N
        for j=1:i-1
            ndu(:,i)  = ndu(:,i) - (ndu(:,i)'*diag(wp)*ndu(:,j))*ndu(:,j);
        end
        ndu(:,i) = ndu(:,i)/sqrt(ndu(:,i)'*diag(wp)*ndu(:,i));
    end
    
    
    C=ComputeCovBasis(ndY,wr);
    
    
    Lambda(:,n) = eig(C);
    [R,a,dummy]=svd(C);
    U_BO = ndu*R;
    Y_BO = ndY*R;
    time(n) = t;
    
    
    
    
    Y = ndY;
    U = ndu;
    dY = ndY;
    dubar = ndubar;
    du = ndu;
    
    [U_KL,V_KL,L_KL,Amean] = Compute_KL(A(:,:,n+ns),wp,wr,N);
    Lambda_KL(:,n) = L_KL(1:N)';
    Error = Amean - dubar;
    lambda_mean(n) = Amean'*diag(wp)*Amean;
    E_mean(n) = Error'*diag(wp)*Error;
    Error = A(:,:,n+ns) - repmat(dubar,1,Nr) -U*Y';
    
    Error = Error'*diag(wp)*Error;
    E_var(n)= sum(eig(Error*diag(wr)));
   
    E_KL(n)= sum(L_KL(N+1:end));
    
    if mod(n,5)==0
        
        subplot(2,1,1)
        disp(['t=' num2str(n*dt) ' is being processed'])
        semilogy(t,Lambda(:,n),'o','color',RGB1); hold on
        
        semilogy(t,L_KL,'+','color',RGB2); hold on
        xlabel('$t$')
        ylabel('$\Lambda(t)$')
        set(gca,'FontSize',15)
        
        
      
        subplot(2,1,2)
        plot(x,U_BO(:,1),'color',RGB1,LW,2); hold on
        plot(x,U_KL(:,1),'--','color',RGB2,LW,2); hold off
        legend('DO','KL')
        xlabel('$x$')
        ylabel('$U(x,t)$')
        xlim([0 2*pi])
        set(gca,'FontSize',15)
        
        
        
        drawnow
        
        
        k = k +1;
    end
    
    n=n+1;	
end

