function [Ud, Vd] = rhs_DB(U,V,Adot)
global w wr

C = V'*diag(wr)*V;
Lr = U'*diag(w)*Adot*diag(wr)*V/C;

Ud   = Adot*diag(wr)*V/C- U*Lr;
Vd   = Adot'*diag(w)*U;
