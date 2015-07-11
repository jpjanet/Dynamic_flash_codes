function Pvap = Antoine_Eq(Ant_coeff,T)
% Function call for Antoine Equation%
% INPUTS: 
%       Ant_coeff: The Antoine eq coefficients ABC given in a matrix, one row per 
%                  component. 
%       T: Temperature in K
% OUTPUTS:
%         Pvap: vector of vapour pressures in Bar(a)

Pvap =  10.^(Ant_coeff(:,1) - Ant_coeff(:,2)./(T + Ant_coeff(:,3) - 273.15) );
end 

