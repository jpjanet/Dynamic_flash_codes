function [phi] = s1_phi(t,F0,w0,CT,Pout,Feedf)
%This is a helper function that links the common function handle Feedf in
%physical variables to the system 1 behvaiour setting.
%   The inputs are the current time t, initial feed rate and composition F0
%   and w0, tank volume CT, pressure outside Pout and the Feedf.m function
%   handle
[w,F] = Feedf(t,w0,F0);
phi = [w*diag(F);CT*ones(1,length(t));Pout*ones(1,length(t))];
end

