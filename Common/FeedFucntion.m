function [w,F] = FeedFucntion(t,w0,F0)
%Provides the feed inhomogeneity for the testcases 1, where the feed is
%constant
% INPUTS: 
%       t: current time 
%       w0: reference feed composition vector
%       F0: reference feed flow rate in kmol/h
% OUTPUTS:
%       w0: feed composition vector at t
%       F0: feed flow rate in kmol/h at t
w = w0*(ones(1,length(t))) ; % constant
F = F0*ones(1,length(t));
end

