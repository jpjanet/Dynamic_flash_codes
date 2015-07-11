 function [J1] = SYSTEM1_JAC_FZ(t,z,dotz,ddotz,phi,Psi_l,dPsi_ldP,dPsi_ldM,Psi_v,dPsi_vdP,K,dKdT,dKdP,Rho_l,dRho_ldx,dRho_ldT, dRho_ldP,Rho_v,dRho_vdy,dRho_vdT,dRho_vdP)
global Nc

%% Set-up and External Functions
J1 = zeros(3*Nc +6, 3*Nc +7);

%% Assemble Groupings

theta1 = -z(Nc + 1:2*Nc).*dKdT;
theta2 = -z(Nc + 1:2*Nc).*dKdP;
xi1 = z(3*Nc + 1)*(Rho_l^(-2))*dRho_ldx';
xi2 = z(3*Nc + 2)*(Rho_v^(-2))*dRho_vdy';
a1 = z(3*Nc + 1)*(Rho_l^(-2))*dRho_ldT + z(3*Nc + 2)*(Rho_v^(-2))*dRho_vdT;
a2 = z(3*Nc + 1)*(Rho_l^(-2))*dRho_ldP + z(3*Nc + 2)*(Rho_v^(-2))*dRho_vdP;
a3 = sum(dKdT.*z(Nc + 1:2*Nc));
a4 = sum(dKdP.*z(Nc + 1:2*Nc));




%% Populate Jacobain Elements

J11 = zeros(Nc,Nc);
J12 = z(3*Nc+4)*eye(Nc);
J13 = z(3*Nc+5)*eye(Nc);
J14 = [zeros(Nc, 3),z(Nc + 1:2*Nc),z(2*Nc + 1:3*Nc),zeros(Nc, 2)];
J21 = eye(Nc);
J22 = -z(3*Nc + 1)*eye(Nc);
J23 = -z(3*Nc + 2)*eye(Nc);
J24 = [-z(Nc + 1: 2*Nc),-z(2*Nc + 1: 3*Nc), zeros(Nc, 5)]; 
J31 = zeros(Nc,Nc);
J32 = -diag(K);
J33 = eye(Nc);
J34 = [zeros(Nc, 5), theta1, theta2];
J41 = [-ones(1,Nc);zeros(5,Nc)];
J42 = [zeros(2,Nc);xi1;zeros(2,Nc);(K-1)'];
J43 = [zeros(2,Nc);xi2;zeros(3,Nc)];
J44 = [zeros(1,2), 1, zeros(1,4);...
       -1, -1, 1, zeros(1,4); ...
       -(1/Rho_l),-(1/Rho_v), zeros(1,3), a1, a2; ...
       zeros(1,2), -dPsi_ldM, 1, zeros(1,2), (-dPsi_ldP); ...
       zeros(1,4), 1, 0, (-dPsi_vdP);
       zeros(1,5), a3, a4];
 
J1 = [J11, J12, J13, J14; J21, J22, J23, J24; J31, J32, J33, J34; J41, J42, J43, J44];
end
