function [ diff , diff_std ] = spline_test_trigonometric_polynomial(  )
%Tests approximation of function 7cos(x)^2-2cos(x)sin(x)-sin(x) on interval
%(0,Pi) with spline constructed by function spline_1_4_06_2.
%Takes no input arguments.
%Returns diff argument which is maximal absolute distance between spline
%and function on the test grid.
%Builds plots of both function and it`s approximation.
x = pi*sort(rand(1,1000));
y = zeros(1,1000);
for i=1:length(x)
    y(i) = 7*cos(x(i))^2-2*cos(x(i))*sin(x(i))-sin(x(i));
end;
coeff = spline_1_4_06_2(x,y,1e-3);
x_test = pi*sort(rand(1,500));
y_test_true = zeros(1,500);
for i=1:length(x_test)
    y_test_true(i) = 7*cos(x_test(i))^2 - ...
        2*cos(x_test(i))*sin(x_test(i)) - sin(x_test(i));
end;
y_test = zeros(1,500);
index = 1;
num = 1;
for i = 1:500
    f = false;
    while ~f && index<=999
        if x_test(i)>=x(index) && x_test(i)<=x(index+1)
            f = true;
            num = index;
        else
            index = index + 1;
        end;
    end;
    y_test(i) = coeff(4*(num-1)+1) * x_test(i)^3 + ...
        coeff(4*(num-1)+2) * x_test(i)^2 + ...
        coeff(4*(num-1)+3) * x_test(i) + coeff(4*num);
end;
plot(x_test,y_test_true,'--go',x_test,y_test,':r*');
diff = max(abs(y_test_true - y_test));
diff_std = max(abs(y_test_true - spline(x,y,x_test)));
return
end

