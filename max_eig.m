function [ e ] = max_eig( A , eps , mki )
%MAX_EIG This function returns [ e ] - e is maximum eigen value of
%matrix A.
%Input [ A eps mki ]. A is matrix, x is starting value for first
%iteration, eps is needed accuracy and mki is maximum number of iterations
%to find the eigen value.
kit = 0;
x_new = ones(length(A(1)));
x = zeros(length(A(1)));
while (norm(x-x_new,Inf)>eps) && (kit < mki)
    x = x_new;
    x = (x/norm(x,2));
    x_new = A*(x);
    kit = kit + 1;
end;
e = norm(x_new , 2)/norm(x,2);
return
end
