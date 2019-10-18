close all
clear
clc

addpath ../../Utilities
addpath ../../Data
set(0,'defaulttextinterpreter','latex')
LW = 'linewidth';


global w wr
global dt t r

Preprocess_Data


r = 3;
m = Ns;
ns = 50;
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



for n=1:kmax-ns
    
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
   
    figure(1)
    subplot(1,2,1)
    semilogy(t(ns:n+ns-1),Lambda(:,1:n),LW,1.2);
    xlabel('$t$')
    ylabel('$\Lambda(t)$')
    set(gca,'FontSize',15);
    
    
    subplot(1,2,2)
    semilogy(t(ns:n+ns-1),E(1:n),LW,1.2);
    xlabel('$t$')
    ylabel('$Error(t)$')
    set(gca,'FontSize',15);
    drawnow

    figure(2)
    colormap(jet)
    for j=1:4
        subplot(1,4,j)
        if (j<=3)
            u = UBO(1:D,j);
            v = UBO(D+1:end,j);
         
            u = reshape(u,jmax,imax);
            v = reshape(v,jmax,imax);
            [curlz,cav]= curl(X,Y,u,v);
            
            contourf(X,Y,curlz,10,'color',[.5 .5 .5]);  axis equal; axis tight; axis([ -2 2 0 6]);
            xlabel('$x$')
            ylabel('$y$')
            set(gca,'FontSize',15);
        else
            T = Tb(:,n+ns); T = reshape(T,jmax,imax);
            contourf(X,Y,T,100,'LineColor','none');  axis equal; axis tight
            xlabel('$x$')
            ylabel('$y$')
            set(gca,'FontSize',15);
            drawnow;
            
        end
    end
    
end
figure
semilogy(t(ns:n+ns-1),Lambda(:,1:n),LW,1.2)
xlabel('$t$')
ylabel('$\Lambda(t)$')


