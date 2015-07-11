function [ Psi_v, dPsi_vdP, dPsi_vdPout, ddPsi_vddP, ddPsi_vdPdPout, ddPsi_vddPout  ] = PSIV( P, Pout, C_v)
%PSIV returns the vapour molar outlfow given the holdup and pressure in the
%tank, the external pressure and the constant C_v
% INPUTS: 
%       P: current pressure in the tank in Bar (abs)
%       Pout: pressure external to the tank in Bas (abs)
%       C_v: proportionality constant 
% OUTPUTS:
%         Psi_v: tank vapour valve function output, kmol/h^3
%         dPsi_vd*: partial derivatives w.r.t to *



Psi_v  = C_v*sqrt(P.^2-Pout.^2);
dPsi_vdP = 0.5*(2*P).*C_v./sqrt(P.^2-Pout.^2);
dPsi_vdPout = -0.5*(2*Pout).*C_v./sqrt(P.^2-Pout.^2);
ddPsi_vddP =  0.5*C_v.*((P.^2-Pout.^2).^(-3/2).*(-1/2).*(2*P) + (P.^2-Pout.^2).^(-1/2)*2);
ddPsi_vdPdPout = 0.5*C_v.*(P.^2 - Pout.^2).^(-3/2).*(2*P).*(-2*Pout);
ddPsi_vddPout =  0.5*C_v.*((P.^2-Pout.^2).^(-3/2).*(-1/2).*(-2*Pout) + (P.^2-Pout.^2).^(-1/2)*(-2));


end

