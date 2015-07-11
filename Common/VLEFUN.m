function [ K dKdT dKdP ddKddT ddKdTdP ddKddP ] = VLEFUN(T,P)
global Nc
%VLEFUN recives T,P as arguments and returns a vector Ki,
%such that y_i = K_i x_i at equilibrium 
% INPUTS: 
%       T: Temperature in K
%       P: Pressure in Bar (abs)
% OUTPUTS:
%         K: vector of equilbrium coefficients
%         dKd*: partial derivatives w.r.t to *

% Data for this system:
Ant_coeff = [3.928, 803.99, 247.04;  ...
            3.932, 935.77, 238.00;  ...
            3.978, 1064.84, 232.00; ...
            4.000, 1170.88, 224.32];
        
% Calcualte vapour pressures:
Pvap = Antoine_Eq(Ant_coeff,T);


% Implement ideal solution behaviour assumptions
K = Pvap/P;
dKdT =  log(10)*( Ant_coeff(:,2)./((T + Ant_coeff(:,3) - 273.15).^2)).*K;
dKdP = -(1/P)*K;
ddKddP = (2/(P^2))*K;
ddKdTdP = -(1/P)*dKdT;
ddKddT = log(10)*K.*( ( -2*Ant_coeff(:,2)./((T + Ant_coeff(:,3) - 273.15).^3)) + ...
         log(10)*(( Ant_coeff(:,2)./((T + Ant_coeff(:,3) - 273.15).^2)).^2));

end
