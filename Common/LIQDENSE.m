function [ Rho_l, dRho_ldx, dRho_ldT, dRho_ldP, ...
           ddRho_lddx, ddRho_ldxdT, ddRho_ldxdP, ...
           ddRho_lddT, ddRho_ldTdP, ddRho_lddP ] = LIQDENSE(x, T, P,rho_R)
%LIQDENSE returns the molar-density of the liquid phase
%given the state variables and composition
% INPUTS:
%         x: vector of liquid phase molar fractions
%         T: Temeperature in K
%         P: Pressure in Bar
%       rho_R: vector of pure component reference
%              molar densities form Rackett Eq.
% OUTPUTS:
%         Rho_v: the molar gas density, kmol/m^3
%         dRho_vd*: partial derivatives w.r.t to *

global Nc

Rho_l = 1./(sum(x./rho_R));
dRho_ldx = -(1./rho_R).*1./((sum(x./rho_R))^2);
dRho_ldT = 0;
dRho_ldP = 0;
ddRho_lddx = 1./((sum(x./rho_R))^3) * (kron((-1/rho_R)',(-2/rho_R))) ;
ddRho_ldxdT = zeros(Nc,1);
ddRho_ldxdP = zeros(Nc,1);
ddRho_lddT = 0;
ddRho_ldTdP = 0;
ddRho_lddP = 0;

end

