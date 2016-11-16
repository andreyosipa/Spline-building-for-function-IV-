function [ x ] = solveMatrixWithSimpleIterationMethod( A , b , eps , mki )
n = length(A(1,:));
M = max_eig(A, eps , mki);
B = M*eye(n) - A;
m = M - max_eig(B,eps,mki);
alpha = 2/(m+M);
xt = ones(n,1);
x = zeros(n,1);
ki = 0;
while (norm( x - xt , Inf)>= eps) && (ki<mki)
    x = xt;
    xt = x - alpha*(A*x - b);
    ki = ki+1;
end;
x = xt;
return
end
