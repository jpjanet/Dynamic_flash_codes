function [Mi,x,y,Ml,Mv,M,L,V,T,P,Pout,F,w] = s1_X_to_variable(tt,X,Tcon,phifun)
% This is a helper function that places the state variable X into the
% physical variables at the times given in tt.
%  The input is defined in the  s1_explicit_function.m function, and
%  consists of the current time t, state vector x, the control and
%  inhomogeneities (phi)
global Nc

Mi = X(:,1:Nc);
x = X(:,Nc + 1 : 2*Nc);
y = X(:,2*Nc + 1 : 3*Nc);
Ml = X(:,3*Nc + 1);
Mv = X(:,3*Nc + 2);
M = X(:,3*Nc + 3);
L = X(:,3*Nc + 4);
V = X(:,3*Nc + 5);
T = Tcon(tt,X);
P = X(:,3*Nc +6);
[ phi] = phifun(tt);
phi = phi';
Fw = phi(:,1:Nc);
F = sum(Fw,2);
w = diag(1./F)*Fw;
Pout = phi(:,end);

end