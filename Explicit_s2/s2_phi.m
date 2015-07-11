function [ phi] = s2_phi(t,w0,F0,Feedf)
%This is a helper function that links the common function handle Feedf in
%physical variables to the system 2 behvaiour setting.
%   The inputs are the current time t, initial feed rate and composition F0
%   and w0 and the Feedf.m function  handle
[w,F] = Feedf(t,w0,F0);
phi = [w;F];
end

