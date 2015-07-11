function [F] = SYSTEM1(t,z,dotz,phi, rho_R,C_l,C_v)
%This function returns the residual of the system 1 equations 
%  The input is defined in the  s1_explicit_function.m function, and
%  consists of the current time t, state vector z, the inhomogeneities 
%  (phi) function handles as well as the specific dynamic
%  parameters calculated in  s1_initial_condition_function.m
global Nc

%% Set-up and External Functions
F = zeros(3*Nc+6, 1); % Initialise
Psi_l = PSIL( z(3*Nc + 3),z(3*Nc + 7),C_l);
Psi_v = PSIV( z(3*Nc + 7), phi(Nc + 2), C_v );
K = VLEFUN(z(3*Nc + 6),z(3*Nc + 7));
Rho_l = LIQDENSE(z(Nc + 1:2*Nc), z(3*Nc + 6), z(3*Nc + 7), rho_R);
Rho_v = GASDENSE(z(2*Nc + 1:3*Nc), z(3*Nc + 6), z(3*Nc + 7));

%% Equation 1: Rows 1 to Nc - Mass Balance on Tank
F(1:Nc) = dotz(1:Nc) - phi(1:Nc) + z(3*Nc + 4)*z(Nc + 1:2*Nc) + z(3*Nc + 5)*z(2*Nc + 1:3*Nc); 

%% Equation 2: Rows Nc + 1 to 2Nc - Mass Balance on Species: 
F(Nc + 1:2*Nc) = z(1:Nc) - z(3*Nc + 1)*z(Nc + 1:2*Nc)-z(3*Nc + 2)*z(2*Nc + 1:3*Nc);

%% Euation 3: Rows 2Nc +1 to 3Nc - VLE
F(2*Nc + 1:3*Nc) = z(2*Nc + 1:3*Nc) - K.*z(Nc + 1:2*Nc);

%% Equation 4: Row 3*Nc + 1 - Total Mass Balance
F(3*Nc + 1) = z(3*Nc + 3) - sum(z(1:Nc));

%% Equation 5: Row 3*Nc + 2 - Phase Mass Balance
F(3*Nc + 2) = z(3*Nc+3) - z(3*Nc + 1) - z(3*Nc + 2);

%% Equation 6: Row 3*Nc + 3 - Volume Balance on Tank
F(3*Nc + 3) = phi(Nc + 1) - z(3*Nc + 2)*(1/Rho_v) - z(3*Nc + 1)*(1/Rho_l);

%% Equation 7: Row 3*Nc + 4 - Liquid Exit Equation
F(3*Nc + 4) = z(3*Nc +4) - Psi_l;

%% Equation 8: Row 3*Nc + 5 - Vapor Exit Equation
F(3*Nc + 5) = z(3*Nc +5) - Psi_v;

%% Equation 9: Row 3*Nc + 6 -  Species Fraction Complementarity Condition 
F(3*Nc + 6) = sum((K-1).*z(Nc + 1:2*Nc));

end
