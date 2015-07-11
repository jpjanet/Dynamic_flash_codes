function [theta, xi, alpha] = groupseval(t,z,dotz,ddotz,phi,dotphi, rho_R, C_l, C_v)
%This function returns the the various groups in the terms of the System 1
%Jacobian, as defined in Appendix A.2
%  The input consists of the current time t, state vector z, as well as 
%  some derivatives. As the system posed here is index 1 and semi-explicit,
%  the terms in the derivative are not used later.
global Nc
[ Psi_l dPsi_ldM dPsi_ldP  ddPsi_lddM ddPsi_ldMdP ddPsi_lddP ] = PSIL( z(3*Nc + 3),z(3*Nc + 7),C_l);
[ Psi_v dPsi_vdP dPsi_vdPout ddPsi_vddP ddPsi_vdPdPout ddPsi_vddPout  ] = PSIV( z(3*Nc + 7), phi(Nc + 2),C_v );
[ K dKdT dKdP ddKddT ddKdTdP ddKddP ] = VLEFUN(z(3*Nc + 6),z(3*Nc + 7));
[ Rho_l dRho_ldx dRho_ldT dRho_ldP ...
  ddRho_lddx ddRho_ldxdT ddRho_ldxdP ...
  ddRho_lddT ddRho_ldTdP ddRho_lddP ] = LIQDENSE(z(Nc + 1:2*Nc), z(3*Nc + 6), z(3*Nc + 7), rho_R);
[ Rho_v dRho_vdy dRho_vdT dRho_vdP ... 
  ddRho_vddy ddRho_vdydT ddRho_vdydP ...
  ddRho_vddT ddRho_vdTdP ddRho_vddP]  = GASDENSE(z(2*Nc + 1:3*Nc), z(3*Nc + 6), z(3*Nc + 7));

theta1 = -z(Nc + 1:2*Nc).*dKdT;
theta2 = -z(Nc + 1:2*Nc).*dKdP;
theta3 = - (dKdT*dotz(3*Nc + 6) + dKdP*dotz(3*Nc + 7));
theta4 = (-(ddKddT*dotz(3*Nc + 6) + ddKdTdP*dotz(3*Nc + 7)).*z(Nc + 1:2*Nc) - dKdT.*dotz(Nc + 1:2*Nc));
theta5 = (-(ddKdTdP*dotz(3*Nc + 6) + ddKddP*dotz(3*Nc + 7)).*z(Nc + 1:2*Nc) - dKdP.*dotz(Nc + 1:2*Nc));
theta6 = z(Nc + 1:2*Nc)*(1/z(3*Nc + 1)).*K;
theta7 = z(2*Nc + 1:3*Nc)*(1/z(3*Nc + 1)).*K;
theta8 = 1./(1+K*z(3*Nc+2)/z(3*Nc+1));

theta = [theta1,theta2,theta3,theta4,theta5,theta6,theta7,theta8];

xi1 = z(3*Nc + 1)*(Rho_l^(-2))*dRho_ldx';
xi2 = z(3*Nc + 2)*(Rho_v^(-2))*dRho_vdy';
beta1 = (sum(dRho_ldx.*dotz(Nc + 1:2*Nc)) + dRho_ldT*dotz(3*Nc + 6) + dRho_ldP*dotz(3*Nc + 7));
tilde_xi3 = dotz(3*Nc + 1) * (Rho_l^(-2))*dRho_ldx' -2*z(3*Nc + 1)*(Rho_l^(-3))*dRho_ldx'*beta1;
hat_xi3 = z(3*Nc + 1)*(Rho_l^(-2))*(sum(diag(dotz(Nc + 1:2*Nc)),1)*ddRho_lddx' + ddRho_ldxdT'*dotz(3*Nc + 6) +  ddRho_ldxdP'*dotz(3*Nc + 7));
xi3 = hat_xi3 + tilde_xi3;
xi4 = (dKdT'*dotz(3*Nc + 6) + dKdP'*dotz(3*Nc + 7));
beta2 = (sum(dRho_vdy.*dotz(2*Nc + 1:3*Nc)) + dRho_vdT*dotz(3*Nc + 6) + dRho_vdP*dotz(3*Nc + 7));
tilde_xi5 = dotz(3*Nc + 2) * (Rho_v^(-2))*dRho_vdy' -2*z(3*Nc + 2)*(Rho_v^(-3))*dRho_vdy'*beta2;
hat_xi5 = z(3*Nc + 2)*(Rho_v^(-2))*(sum(diag(dotz(2*Nc + 1:3*Nc)),1)*ddRho_vddy' + ddRho_vdydT'*dotz(3*Nc + 6) +  ddRho_vdydP'*dotz(3*Nc + 7));
xi5 = hat_xi5 + tilde_xi5;
%xi6 = [z(3*Nc + 1)*(Rho_l^(-2))*dRho_ldx'];
xi6 = xi1;
% changed to xi1-K*xi2
xi7 = K' - 1;
%xi8 = [z(3*Nc + 2)*(Rho_v^(-2))*dRho_vdy'];
xi8 = xi2;
%changed to xi 2

xi9 = -z(3*Nc + 4)-K'*z(3*Nc + 5);
xi10 =xi1+xi2.*K';
xi11 = (z(3*Nc + 2)/z(3*Nc + 1))*(xi10);
xi12 = (z(3*Nc + 2)/z(3*Nc + 1))*(xi7);
xi13 = -(z(3*Nc + 1) +K'*z(3*Nc + 2));
xi14 = xi1 + K'.*xi2;
%xi15 = K' - 1;
% changed to xi1
xi = [xi1;xi2;xi3;xi4;xi5;xi6;xi7;xi8;xi9;xi10;xi11;xi12;xi13;xi14];

alpha = zeros(41,1);
alpha(1) = z(3*Nc + 1)*(Rho_l^(-2))*dRho_ldT + z(3*Nc + 2)*(Rho_v^(-2))*dRho_vdT;
alpha(2) = z(3*Nc + 1)*(Rho_l^(-2))*dRho_ldP + z(3*Nc + 2)*(Rho_v^(-2))*dRho_vdP;
alpha(3) = sum(dKdT.*z(Nc + 1:2*Nc));
alpha(4) = sum(dKdP.*z(Nc + 1:2*Nc));
alpha(5) = (Rho_l^(-2))*(sum(dRho_ldx.*dotz(Nc + 1:2*Nc)) + dRho_ldT*dotz(3*Nc + 6) + dRho_ldP*dotz(3*Nc + 7) );
alpha(6) = (Rho_v^(-2))*(sum(dRho_vdy.*dotz(2*Nc + 1:3*Nc)) + dRho_vdT*dotz(3*Nc + 6) + dRho_vdP*dotz(3*Nc + 7) );
tilde_alpha7 = (dotz(3*Nc + 1)*(Rho_l^(-2))*dRho_ldT - 2*z(3*Nc + 1)*(Rho_l^(-3))*dRho_ldT*beta1 + z(3*Nc + 1)*(Rho_l^(-2))*(sum(ddRho_ldxdT.*dotz(Nc + 1:2*Nc)) + ddRho_lddT*dotz(3*Nc + 6) + ddRho_ldTdP*dotz(3*Nc + 7)));
hat_alpha7 = (dotz(3*Nc + 2)*(Rho_v^(-2))*dRho_vdT - 2*z(3*Nc + 2)*(Rho_v^(-3))*dRho_vdT*beta2 + z(3*Nc + 2)*(Rho_v^(-2))*(sum(ddRho_vdydT.*dotz(2*Nc + 1:3*Nc)) + ddRho_vddT*dotz(3*Nc + 6) + ddRho_vdTdP*dotz(3*Nc + 7)));
alpha(7) = tilde_alpha7 + hat_alpha7;
tilde_alpha8 = (dotz(3*Nc + 1)*(Rho_l^(-2))*dRho_ldP - 2*z(3*Nc + 1)*(Rho_l^(-3))*dRho_ldP*beta1 + z(3*Nc + 1)*(Rho_l^(-2))*(sum(ddRho_ldxdP.*dotz(Nc + 1:2*Nc)) + ddRho_ldTdP*dotz(3*Nc + 6) + ddRho_lddP*dotz(3*Nc + 7)));
hat_alpha8 = (dotz(3*Nc + 2)*(Rho_v^(-2))*dRho_vdP - 2*z(3*Nc + 2)*(Rho_v^(-3))*dRho_vdP*beta2 + z(3*Nc + 2)*(Rho_v^(-2))*(sum(ddRho_vdydP.*dotz(2*Nc + 1:3*Nc)) + ddRho_vdTdP*dotz(3*Nc + 6) + ddRho_vddP*dotz(3*Nc + 7)));
alpha(8) = tilde_alpha8 + hat_alpha8;
alpha(9) = -(ddPsi_ldMdP*dotz(3*Nc + 7) + ddPsi_lddM*dotz(3*Nc + 3));
alpha(10) = -(ddPsi_lddP*dotz(3*Nc + 7) + ddPsi_ldMdP*dotz(3*Nc + 3));
alpha(11) = - (ddPsi_vddP*dotz(3*Nc + 7) + ddPsi_vdPdPout*dotphi(Nc + 2));
alpha(12) = sum(dKdT.*dotz(Nc + 1:2*Nc) + z(Nc + 1:2*Nc).*(ddKddT*dotz(3*Nc + 6) + ddKdTdP*dotz(3*Nc + 7)));
alpha(13) = sum(dKdP.*dotz(Nc + 1:2*Nc) + z(Nc + 1:2*Nc).*(ddKdTdP*dotz(3*Nc + 6) + ddKddP*dotz(3*Nc + 7)));
alpha(14) = z(3*Nc + 1)*(Rho_l^(-2))*dRho_ldT + z(3*Nc + 2)*(Rho_v^(-2))*dRho_vdT;
alpha(15) = z(3*Nc + 1)*(Rho_l^(-2))*dRho_ldP + z(3*Nc + 2)*(Rho_v^(-2))*dRho_vdP;
alpha(16) = sum(dKdT.*z(Nc + 1:2*Nc));
alpha(17) = sum(dKdP.*z(Nc + 1:2*Nc));
alpha(18) = -z(3*Nc + 2)*sum(theta(:,1));
alpha(19) = -z(3*Nc + 2)*sum(theta(:,2));
alpha(20) = alpha(1) + sum(theta(:,1)'.*xi(2,:));
alpha(21) = alpha(2) + sum(theta(:,2)'.*xi(2,:));
alpha(22) = sum((xi(10,:)'.*z(Nc + 1:2*Nc ))./xi(13,:)') - (Rho_l^(-1));
alpha(23) = sum((xi(10,:)'.*z(2*Nc + 1:3*Nc ))./xi(13,:)') - (Rho_v^(-1));
alpha(24) = -z(3*Nc + 2)*sum((xi(10,:)'.*theta(:,1 ))./xi(13,:)') + alpha(1);
alpha(25) = -z(3*Nc + 2)*sum((xi(10,:)'.*theta(:,2 ))./xi(13,:)') + alpha(2);
alpha(26) = sum((xi(7,:)'.*z(Nc + 1:2*Nc ))./xi(13,:)');
alpha(27) = sum((xi(7,:)'.*z(2*Nc + 1:3*Nc ))./xi(13,:)');
alpha(28) = -z(3*Nc + 2)*sum((xi(7,:)'.*theta(:,1 ))./xi(13,:)') + alpha(3);
alpha(29) = -z(3*Nc + 2)*sum((xi(7,:)'.*theta(:,2 ))./xi(13,:)') + alpha(4);
%thb
tilde_alpha_30 = -(1/(z(3*Nc  + 1)))*(xi(10,:)*z(Nc + 1:2*Nc))-(Rho_l^(-1));
hat_alpha_30 = xi(11,:)*(theta(:,8).*theta(:,6));
alpha(30) = tilde_alpha_30 + hat_alpha_30 ;
tilde_alpha_31 = -(1/(z(3*Nc  + 1)))*(xi(10,:)*z(2*Nc + 1:3*Nc))-(Rho_v^(-1));
hat_alpha_31 = xi(11,:)*(theta(:,8).*theta(:,7));
alpha(31) = tilde_alpha_31 + hat_alpha_31 ;
tilde_alpha_32 = xi(12,:)*(theta(:,8).*theta(:,6));
hat_alpha_32 = -(1/(z(3*Nc  + 1)))*(xi(7,:)*z(Nc + 1:2*Nc));
alpha(32) = tilde_alpha_32 + hat_alpha_32;
tilde_alpha_33 = xi(12,:)*(theta(:,8).*theta(:,7));
hat_alpha_33 = -(1/(z(3*Nc  + 1)))*(xi(7,:)*z(2*Nc + 1:3*Nc));
alpha(33) = tilde_alpha_33 + hat_alpha_33;
tilde_alpha_34 = sum(xi(11,:).*theta(:,8)'.*theta(:,1)');
hat_alpha_34 = alpha(14) - xi(2,:)*theta(:,1);
alpha(34) = hat_alpha_34 + tilde_alpha_34;
tilde_alpha_35 = sum(xi(11,:).*theta(:,8)'.*theta(:,2)');
hat_alpha_35 = alpha(15) - xi(2,:)*theta(:,2);
alpha(35) = hat_alpha_35 + tilde_alpha_35;
alpha(36) = sum(xi(12,:).*theta(:,8)'.*theta(:,1)') + alpha(16);
alpha(37) = sum(xi(12,:).*theta(:,8)'.*theta(:,2)') + alpha(17);



alpha(38) = z(3*Nc + 2)*sum(theta(:,1));
alpha(39) = z(3*Nc + 2)*sum(theta(:,2));
alpha(40) = alpha(1) - sum(xi(2,:)'.*theta(:,1));
alpha(41) = alpha(2) - sum(xi(2,:)'.*theta(:,2));


end
