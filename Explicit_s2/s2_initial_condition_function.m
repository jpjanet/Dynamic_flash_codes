function [z0, dotz0, C_l, C_v, rho_R, ierr ] = s2_initial_condition_function(Tref,phi,P,Pout,Critical_props,PSI0)
%This function determines the consistent initial conditions 
% system 2 equations.
%   The input is defined in the dynamic_flash_mainfile.m script
%   The outputs are the various state varaibles as vectors at the times
%   specified in the vector tt.
global Nc
V_R = (1E-6).*Critical_props(:, 3).*(Critical_props(:, 4).^((1  ...
    - (Tref./Critical_props(:, 1))).^(2/7))); % Ref molar vol, m^3/mol
rho_R = (1E-3)*1./V_R; % Convert to kmol

[ K ] = VLEFUN(Tref,P); % Fetch Ki values

w = phi(1:Nc);
F = phi(Nc + 1);

rachford_rice_fun = @(PSI)(sum((w.*(1-K))./ ...
    (1 + PSI*(K - 1)))); % Rachford-Rice Function

options = optimset('Tolx',1E-14,'TolFun',1e-14);
[PSI,~,flag,~] = fzero(rachford_rice_fun,PSI0,options);

if (flag == 1)
    ierr = 0;
    x = w./(1 + PSI*(K - 1));
    y = K.*x;
    V = PSI*F;
    L = F - V;
    y = (F*w - L*x)/V;
    
    [ Rho_l ] = LIQDENSE(x, Tref, P,  rho_R); % Fetch liquid molar density
    [ Rho_v ] = GASDENSE(y, Tref, P);
    
    Ml = (1.5*0.5)*Rho_l; % Liquid hold-up in kmol
    Mv = (1.5*0.5)*Rho_v; % Vapour hold-up in kmol
    Mi = Ml*x + Mv*y;
    M = Ml + Mv;
    C_l = L/sqrt(M);
    C_v = V/sqrt(P^2-Pout^2);
    [z0] = s2_variable_to_X(x,y,M,L,V,Tref,P);
    dotz0 = 0*z0; % Steady-state solution by defintion
else
    ierr = -1; % return error code
    z0 = 0;
    dotz0 = 0;
end



