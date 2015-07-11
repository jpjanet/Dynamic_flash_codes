function [theta, xi, alpha] = groupseval2(t,z,dotz,ddotz,phi,dotphi,C_l)
%This function returns the the various groups in the terms of the System 2
%Jacobian, as defined in Appendix B.2
%  The input consists of the current time t, state vector z, as well as 
%  some derivatives. As the system posed here is index 1 and semi-explicit,
%  the terms in the derivative are not used later.
global Nc
[ Psi_l dPsi_ldM dPsi_ldP  ddPsi_lddM ddPsi_ldMdP ddPsi_lddP ] = PSIL( z(2*Nc + 1),z(2*Nc + 5),C_l);
[ K dKdT dKdP ddKddT ddKdTdP ddKddP ] = VLEFUN(z(2*Nc + 4),z(2*Nc + 5));


theta1 = z(Nc + 1:2*Nc) - z(1:Nc);
theta2 = -z(1:Nc).*dKdT;
theta3 = -z(1:Nc).*dKdP;
theta4 = dotz(Nc +1:2*Nc) - dotz(1:Nc);
theta5 =  -(dKdT*dotz(2*Nc + 4) + dKdP*dotz(2*Nc + 5));
theta6 = (-(ddKddT*dotz(2*Nc + 4) + ddKdTdP*dotz(2*Nc + 5)).*z(1:Nc) - dKdT.*dotz(1:Nc));
theta7 = (-(ddKdTdP*dotz(2*Nc + 4) + ddKddP*dotz(2*Nc + 5)).*z(1:Nc) - dKdP.*dotz(1:Nc));

theta=[theta1,theta2,theta3,theta4,theta5,theta6,theta7];

xi1 = K' - 1;
xi2 = (dKdT'*dotz(2*Nc + 4) + dKdP'*dotz(2*Nc + 5));
xi = [xi1;xi2];

alpha = zeros(13,1);
alpha(1) = phi(Nc + 1) - z(2*Nc + 3);
alpha(2) = sum(dKdT.*z(1:Nc));
alpha(3) = sum(dKdP.*z(1:Nc));
alpha(4) = dotphi(Nc + 1) - dotz(2*Nc + 3);
alpha(5) = dotz(2*Nc + 1) + phi(Nc + 1) - z(2*Nc +3);
alpha(6) = -(ddPsi_ldMdP*dotz(2*Nc + 5) + ddPsi_lddM*dotz(2*Nc + 1));
alpha(7) = -(ddPsi_lddP*dotz(2*Nc + 5) + ddPsi_ldMdP*dotz(2*Nc + 1));
alpha(8) = sum(dKdT.*dotz(1:Nc) + (ddKddT*dotz(2*Nc + 4) + ddKdTdP*dotz(2*Nc + 5)).*z(1:Nc));
alpha(9) = sum(dKdP.*dotz(1:Nc) + (ddKdTdP*dotz(2*Nc + 4) + ddKddP*dotz(2*Nc + 5)).*z(1:Nc));
alpha(10) = dPsi_ldM + dPsi_ldP;
alpha(11) = -alpha(3)/alpha(2);
alpha(12) = (dPsi_ldM^(-1))*(1+ dPsi_ldP*(alpha(2)/alpha(3)));
alpha(13) = (-alpha(2)/alpha(3));
end
