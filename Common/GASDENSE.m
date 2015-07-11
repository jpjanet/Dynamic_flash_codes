function [ Rho_v, dRho_vdy, dRho_vdT, dRho_vdP, ... 
           ddRho_vddy, ddRho_vdydT, ddRho_vdydP, ...
           ddRho_vddT, ddRho_vdTdP, ddRho_vddP] = GASDENSE(y, T, P)
%GASDENSE returns the molar-density of the gas phase
%given the state variables and composition
%from ideal gas law
% INPUTS:
%         y: vector of gas phase molar fractions
%         T: Temeperature in K
%         P: Pressure in Bar
% OUTPUTS:
%         Rho_v: the molar gas density, kmol/m^3
%         dRho_vd*: partial derivatives w.r.t to *
global Nc

R = 8.3144621E-2;	% Gas constant

Rho_v = P./(R*T);
dRho_vdy = zeros(Nc,1);
dRho_vdT = - P./(R*T.*T);
dRho_vdP = 1./(R*T);
ddRho_vddy = zeros(Nc);
ddRho_vdydT = zeros(Nc,1);
ddRho_vdydP = zeros(Nc,1);
ddRho_vddT = 2*P./(R*T.*T.*T);
ddRho_vdTdP = - 1./(R*T.*T);
ddRho_vddP = 0;

end

