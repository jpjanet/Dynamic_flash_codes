function [ T] = TcontrolX(t,X,Tref)
%Provides the control temperature for the testcases 2, where temperature is
%constant
% INPUTS: 
%       t,X: current time and system state vector
%       Tref: Reference Temperature in K
%       Target: Liquid outflow tagert in kmol/h^3
% OUTPUTS:
%         T: control temperature

s =0;
hold = 6;
Tmax = Tref + s;
Tmin = Tref - s;
T =  [t<=12].*(Tref + ((Tmax-Tref)/12).*t) ...
    + [t>12].*[t<=(12+hold)]*(Tmax) ...
    + [t>(12+hold)].*[t<=(12+hold+24)].*(Tmax - ((Tmax-Tmin)/24).*(t-(12+hold))) ...
    + [t>(12+hold+24)].*(Tmin);

end

