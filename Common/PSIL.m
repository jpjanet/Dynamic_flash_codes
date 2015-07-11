function [ Psi_l, dPsi_ldM, dPsi_ldP,  ddPsi_lddM, ddPsi_ldMdP, ddPsi_lddP ] = PSIL( M,P,C_l)
%PSIL returns the liquid molar outlfow given the holdup and pressure in the
%tank, and the constant C_l
% INPUTS: 
%       M: current molar holdup in the tank in kmol
%       P: current pressure in the tank in Bar (abs)
%       C_l: proportionality constant 
% OUTPUTS:
%         Psi_l: tank liquid hydrodynamic function output, kmol/h^3
%         dPsi_ld*: partial derivatives w.r.t to *

Psi_l  = C_l*sqrt(M);
dPsi_ldM = 0.5*C_l/sqrt(M);
dPsi_ldP = 0;
ddPsi_lddM = -0.25*C_l*(M^(-3/2));
ddPsi_ldMdP = 0;
ddPsi_lddP = 0;
end

