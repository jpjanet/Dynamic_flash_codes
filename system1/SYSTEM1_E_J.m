function [J] = SYSTEM1_E_J(t,x,Tcon,phifun,rho_R)
%This function returns the Jacobian of the right hand side of the system 1 
% equations in the semi-implicit form, i.e. not the derivative terms.
%  The input is defined in the  s1_explicit_function.m function, and
%  consists of the current time t, state vector x, the control and
%  inhomogeneities (phi) function handles. 

global Nc C_l C_v 


[T] = Tcon(t,x);
[phi] = phifun(t); 
z = [x(1:3*Nc + 5);T;x(3*Nc + 6)]; % combined control
[ ~, dPsi_ldM, dPsi_ldP ] = PSIL( z(3*Nc + 3),z(3*Nc + 7), C_l);
[ ~, dPsi_vdP ] = PSIV( z(3*Nc + 7), phi(Nc + 2),C_v );
[ K, ~, dKdP ] = VLEFUN(z(3*Nc + 6),z(3*Nc + 7));
[ Rho_l, dRho_ldx, ~, dRho_ldP] = LIQDENSE(z(Nc + 1:2*Nc), z(3*Nc + 6), z(3*Nc + 7), rho_R);
[ Rho_v, dRho_vdy, ~, dRho_vdP ]  = GASDENSE(z(2*Nc + 1:3*Nc), z(3*Nc + 6), z(3*Nc + 7));

%% Assemble Groupings

theta2 = -z(Nc + 1:2*Nc).*dKdP;
xi1 = z(3*Nc + 1)*(Rho_l^(-2))*dRho_ldx';
xi2 = z(3*Nc + 2)*(Rho_v^(-2))*dRho_vdy';
a2 = z(3*Nc + 1)*(Rho_l^(-2))*dRho_ldP + z(3*Nc + 2)*(Rho_v^(-2))*dRho_vdP;
a4 = sum(dKdP.*z(Nc + 1:2*Nc));



%% Populate Jacobain Elements

J11 = zeros(Nc,Nc);
J12 = -z(3*Nc+4)*eye(Nc);
J13 = -z(3*Nc+5)*eye(Nc);
J14 = -[zeros(Nc, 3),z(Nc + 1:2*Nc),z(2*Nc + 1:3*Nc),zeros(Nc, 1)];
J21 = eye(Nc);
J22 = -z(3*Nc + 1)*eye(Nc);
J23 = -z(3*Nc + 2)*eye(Nc);
J24 = [-z(Nc + 1: 2*Nc),-z(2*Nc + 1: 3*Nc), zeros(Nc, 4)];
J31 = zeros(Nc,Nc);
J32 = -diag(K);
J33 = eye(Nc);
J34 = [zeros(Nc, 5), theta2];
J41 = [-ones(1,Nc);zeros(5,Nc)];
J42 = [zeros(2,Nc);xi1;zeros(2,Nc);(K-1)'];
J43 = [zeros(2,Nc);xi2;zeros(3,Nc)];
J44 = [zeros(1,2), 1, zeros(1,3);...
    -1, -1, 1, zeros(1,3); ...
    -(1/Rho_l),-(1/Rho_v), zeros(1,3), a2; ...
    zeros(1,2), -dPsi_ldM, 1, zeros(1,1), (-dPsi_ldP); ...
    zeros(1,4), 1, (-dPsi_vdP);
    zeros(1,5),  a4];

J = [J11, J12, J13, J14; J21, J22, J23, J24; J31, J32, J33, J34; J41, J42, J43, J44];


end