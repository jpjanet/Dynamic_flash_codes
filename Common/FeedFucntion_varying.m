function [w,F] = FeedFucntion_varying(t,w0,F0)
%Provides the feed inhomogeneity for the testcases 2 and 3, where the feed
%varies in composition and flow rate linearly. 
% INPUTS: 
%       t: current time 
%       w0: reference feed composition vector
%       F0: reference feed flow rate in kmol/h
% OUTPUTS:
%       w0: feed composition vector at t
%       F0: feed flow rate in kmol/h at t
w = w0*(ones(1,length(t))) - (w0*(ones(1,length(t)))...
    - flipud(w0)*(ones(1,length(t))))*(diag(t/48).*diag(t<=48))...
    - (w0*(ones(1,length(t))) ...
    - flipud(w0)*(ones(1,length(t))))*(diag(t>48)) ;

F = F0*ones(1,length(t))+F0*ones(1,length(t))*0.5*(diag(t/48).*diag(t<=48))...
    + F0*ones(1,length(t))*0.5*diag(t>48);

end

