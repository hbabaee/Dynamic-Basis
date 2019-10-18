function Ad = ComputeTimeDerivative(A)
global dt
Ad = zeros(size(A));

Ad(:,2:end-1) = (   A(:,3:end)     -          A(:,1:end-2))/(2*dt);
Ad(:,3:end-2) = (   8*(A(:,4:end-1) - A(:,2:end-3)) - (A(:,5:end) - A(:,1:end-4)))/(12*dt);
Ad(:,1)       = (-3*A(:,1)   + 4*A(:,2) -     A(:,3))/(2*dt);
Ad(:,end)     = ( 3*A(:,end) - 4*A(:,end-1) + A(:,end-2))/(2*dt);