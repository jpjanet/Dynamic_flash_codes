function J = numjac(fun,x,delta)
%This function returns a numerical approximation to the Jacobian of vector
%valued function hand fun at the point x using stepsize delta. The mehtod
%is second order central difference.
fdim = length(fun(x));
xdim = length(x);
I = eye(xdim);
fI = eye(fdim);

J = zeros(fdim,xdim);
for i = 1:fdim
    for j = 1:xdim
        J(i,j) = (1/(2*delta))*(fun(x + delta*I(:,j)) - fun(x - delta*I(:,j)))'*(fI(:,i));
    end
end


end 