function [ V, dPsi_vdP] = VcontrolX(t,X,Pout,C_v)
%Provides the control on vapour outflow for the system 2, in the same was
%as is set in the system 1 formulation
% INPUTS: 
%       t,X: current time and system state vector
%       Pout: external pressure in Bar (abs)
%       C_v: proportionality constant
% OUTPUTS:
%         V: vapour outflow rate in kmol/h
%         dPsi_vdP: partial derivative w.r.t to pressure  
global Nc
[V,  dPsi_vdP] = PSIV( X(2*Nc + 3,:), Pout, C_v );
end

