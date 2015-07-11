function [x,y,M,L,V,T,P,F,w] = s2_X_to_variable(tt,X,Tcon,Vcon,phifun)
% This is a helper function that places the state variable X into the
% physical variables at the times given in tt.
%  The input is defined in the  s1_explicit_function.m function, and
%  consists of the current time t, state vector x, the control and
%  inhomogeneities (phi)
global Nc
T = Tcon(tt,X);
V = Vcon(tt,X')';
[ phi] = phifun(tt);
phi = phi';
x = X(:,1:Nc);
y = X(:,Nc +1 : 2*Nc);
M = X(:,2*Nc + 1);
L = X(:,2*Nc + 2);
P = X(:,2*Nc + 3);
F = phi(:,Nc +1);
w = phi(:,1:Nc);
end