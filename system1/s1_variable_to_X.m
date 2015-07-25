function [X]  = s1_variable_to_X(Mi,x,y,Ml,Mv,M,L,V,T,P)
%This is a helper function that places  physical variables at the time t
%into the vector X (incl temperature).
%  The input is defined in the  s1_explicit_function.m function, and
%  consists of the current time t physical variables
X =[Mi;x;y;Ml;Mv;M;L;V;T;P];

end