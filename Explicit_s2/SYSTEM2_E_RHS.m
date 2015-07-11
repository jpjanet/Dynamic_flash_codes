function [F] = SYSTEM2_E_RHS(t,x,Tcon,Vcon,phifun,C_l)
%This function returns the right hand side of the system 2 equations in the
%semi-implicit form, i.e. not the derivative terms
%  The input is defined in the  s2_explicit_function.m function, and
%  consists of the current time t, state vector x, the control and
%  inhomogeneities (phi) function handles as well as the specific dynamic
%  parameters calculated in  s2_initial_condition_function.m
global Nc

%% Set-up and External Functions
[T] = Tcon(t,x);
[V, ~] = Vcon(t,x);
[phi] = phifun(t);
z = [x(1:2*Nc +2);V;T;x(2*Nc +3)];
F = zeros(2*Nc+3, 1); % Initialise
Psi_l = PSIL( z(2*Nc + 1),z(2*Nc + 5),C_l);
K = VLEFUN(z(2*Nc + 4),z(2*Nc + 5));

%% Equation 1: Rows 1 to Nc - Mass Balance on Tank
F(1:Nc) = -1*(-1*phi(Nc+1)*(phi(1:Nc)-z(1:Nc))+z(2*Nc+3)*(z(Nc + 1:2*Nc)-z(1:Nc)))/z(2*Nc +1);

%% Euation 2: Rows Nc +1 to 2Nc - VLE
F(Nc + 1:2*Nc) = z(Nc + 1:2*Nc) - K.*z(1:Nc);

%% Equation 3: Row 2*Nc + 1 - Total Mass Balance
F(2*Nc + 1) = -1*(- phi(Nc + 1) + z(2*Nc + 2) +z(2*Nc +3));

%% Equation 4: Row 2*Nc + 2 - Liquid Exit Equation
F(2*Nc + 2) = z(2*Nc +2) - Psi_l;

%% Equation 5: Row 2*Nc + 3 -  Species Fraction Complementarity Condition 
F(2*Nc + 3) = sum((K-1).*z(1:Nc));

end
