function [J] = SYSTEM2_E_J(t,x,Tcon,Vcon, phifun)
%This function returns the Jacobian of the right hand side of the system 2 
% equations in the semi-implicit form, i.e. not the derivative terms.
%  The input is defined in the  s2_explicit_function.m function, and
%  consists of the current time t, state vector x, the control and
%  inhomogeneities (phi) function handles. 

global Nc C_l

%% Apply controls /  inhomogenieties 

[T] = Tcon(t,x);
[V, dPsi_vdP] = Vcon(t,x);
[phi] = phifun(t);
z = [x(1:2*Nc +2);V;T;x(2*Nc +3)];

%% Fetch groupings

[ ~, dPsi_ldM, dPsi_ldP  , ~, ~ , ~ ] = PSIL( z(2*Nc + 1),z(2*Nc + 5),C_l);
[ K , ~, ~, ~, ~, ~ ] = VLEFUN(z(2*Nc + 4),z(2*Nc + 5));
[theta, xi, alpha] = groupseval2(t,z,0*z,0*z,phi,0*phi,C_l);



%% Populate Jacobain Elements

J11 = -eye(Nc)*alpha(1)*(1/z(2*Nc + 1));
J12 = -z(2*Nc+3)*eye(Nc)*(1/z(2*Nc + 1));
X = -1*(phi(Nc+1)*(phi(1:Nc)-z(1:Nc))-z(2*Nc+3)*(z(Nc + 1:2*Nc)-z(1:Nc)))/(z(2*Nc +1)^2);
J13 = [X,zeros(Nc,1),-theta(:,1)*(1/z(2*Nc + 1)),zeros(Nc,2)];

J21 = diag(-K);
J22 = eye(Nc);
J23 = [zeros(Nc,3),theta(:,2:3)];

J31 = [zeros(2,Nc);xi(1,:)];
J32 = zeros(3,Nc);
J33 = [[0,-1,-1,0,0];[-dPsi_ldM,1,0,0,-dPsi_ldP];[0,0,0,alpha(2),alpha(3)]];

J13 = J13(:,[1:2,5]) + [zeros(Nc,2),J13(:,3)*dPsi_vdP];
J23 = J23(:,[1:2,5]) + [zeros(Nc,2),J23(:,3)*dPsi_vdP];
J33 = J33(:,[1:2,5]) + [zeros(3,2),J33(:,3)*dPsi_vdP];

J = [J11, J12, J13; J21, J22, J23; J31, J32, J33;];

end