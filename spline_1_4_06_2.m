function [ coefficients ] = spline_1_4_06_2( x , y , eps)
%Input is 
%    x: grid on Xs axe;
%    y: value of function in points x;
%    eps: accuracy for calculating solution of system of linear equations.
%Output is 
%    coefficients: 4n values, where each 4 is coefficients of spline on the
%    interval of (x(i-1),x(i)) for some i
n = length(x) - 1;
lambda = zeros(n-1,1);
alpha = zeros(n-1,1);
beta = zeros(n-1,1);
gamma = zeros(n-1,1);
c = zeros(n-1,1);
h = zeros(n,1);
for i=1:n
    h(i) = x(i+1)-x(i);
end;
for i=1:(n-1)
    lambda(i) = h(i+1)/(h(i) + h(i+1));
    alpha(i) = -3/h(i)*lambda(i);
    beta(i) = 3/h(i+1) * (1 - lambda(i));
    gamma(i) = -alpha(i) - beta(i);
    c(i) = alpha(i)*y(i) + gamma(i)*y(i+1) + beta(i)*y(i+2);
end;
A = zeros(n-1,n-1);
b = zeros(n-1,1);
A(1,2) = h(1)/h(2);
A(1,1) = 1 + A(1,2);
b(1) = 1/3*c(1) + 2*h(1)/(h(2)^2)*(y(3)-y(2));
for i = 2:n-2
    A(i,i-1) = lambda(i);
    A(i,i) = 2;
    A(i,i+1) = 1 - lambda(i);
    b(i) = c(i);
end;
A(n-1,n-2) = h(n)/h(n-1);
A(n-1,n-1) = 1 + A(n-1,n-2);
b(n-1) = 1/3*c(n-1) + 2*h(n)/(h(n-1)^2)*(y(n)-y(n-1));
m = zeros(n+1,1);
m(2:n) = solveMatrixWithSimpleIterationMethod( A , b , eps , 100 );
m(1) = (c(1) - (1 - lambda(1))*m(3) - 2*m(2))/lambda(1);
m(n+1) = (c(n-1) - 2*m(n) - lambda(n-1)*m(n-1))/(1-lambda(n-1));

coefficients = zeros(4*n,1);
for i=1:n
    tmp(1) = h(i)*m(i) + m(i+1)*h(i) + 2*y(i) - 2*y(i+1);
    tmp(2) = -2*h(i)*m(i) - m(i+1)*h(i) - 3*y(i) + 3*y(i+1);
    tmp(3) = h(i)*m(i);
    tmp(4) = y(i);
    coefficients(4*(i-1)+1) = tmp(1)/(h(i)^3);
    coefficients(4*(i-1)+2) = -3*tmp(1)*x(i)/(h(i)^3) + tmp(2)/(h(i)^2);
    coefficients(4*(i-1)+3) = 3*tmp(1)*(x(i)^2)/(h(i)^3) - ...
        2*tmp(2)*x(i)/(h(i)^2) + tmp(3)/h(i);
    coefficients(4*i) = tmp(4) - tmp(1)*(x(i)^3)/(h(i)^3) + ...
        tmp(2)*(x(i)^2)/(h(i)^2) - tmp(3)*x(i)/h(i);     
end;

return
end

