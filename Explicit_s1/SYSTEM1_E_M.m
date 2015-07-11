function [M] = SYSTEM1_E_M(t,x,Tcon)
%Function handle that returns the mass matrix for System 1 in semi-impilict
%form
%  The input is defined in the  s1_explicit_function.m function, and
%  consists of the current time t, state vector x and the control and
%  function handle
global Nc
[T] = Tcon(t,x);
z = [x(1:end-1);T;x(end)];


M = [eye(Nc,Nc), zeros(Nc, 2*Nc + 6);
     zeros(2*Nc + 6, 3*Nc + 6)];

end

