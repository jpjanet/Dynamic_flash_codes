function [M] = SYSTEM2_E_M()
%Function handle that returns the mass matrix for System 2 in semi-impilict
%form
%  The input is defined in the  s2_explicit_function.m function, which is a
%  constant matrix and hence takes no arguments. 
global Nc
M = [eye(Nc,Nc), zeros(Nc, Nc + 3); ...
     zeros(Nc , 2*Nc + 3);
     zeros(3,2*Nc), [1,0,0;zeros(2,3)]];

end

