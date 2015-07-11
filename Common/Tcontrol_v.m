function [ T] = Tcontrol_v(t,X,Tref,Target,gain,system_num)
global Nc
%Provides the control temperature for the testcases 3, where temperature is
%used to control liquid outflow
% INPUTS: 
%       t,X: current time and system state vector
%       Tref: Reference Temperature in K
%       Target: Liquid outflow tagert in kmol/h^3
%       gain: proportional controller gain, any postive number
%       system_num: switch for system 1 or 2, allowed values are {1,2}
% OUTPUTS:
%         T: control temperature

if system_num == 1
    if isvector(X)
    L = X(3*Nc + 4);
    e = L-Target;
    T =  Tref + gain*e;
    else
    L = X(:,3*Nc + 4);
    e = L-Target;
    T =  Tref + gain*e;
    end

elseif system_num ==2
    if isvector(X)
    L = X(2*Nc + 2);
    e = L-Target;
    T =  Tref + gain*e;
    else
    L = X(:,2*Nc + 2);
    e = L-Target;
    T =  Tref + gain*e;
    end

end


end

