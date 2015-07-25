function [X] = s2_variable_to_X(x,y,M,L,V,T,P)
%This is a helper function that places  physical variables at the time t
%into the vector X (incl temperature).
%  The input is defined in the  s2_explicit_function.m function, and
%  consists of the current time t physical variables
X = [x;y;M;L;V;T;P];
end