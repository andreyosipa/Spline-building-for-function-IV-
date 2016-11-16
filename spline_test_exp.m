function [ diff , diff_std ] = spline_test_exp(  )
x = 0:0.01:1;
y = exp(x);
coeff = spline_1_4_06_2(x,y,1e-4);
x_test = 0.005:0.01:1;
y_test = zeros(1,length(x_test));
for i = 1:length(x_test)
    num = ceil(x_test(i)/0.01);
    y_test(i) = coeff(4*(num-1)+1) * x_test(i)^3 + ...
        coeff(4*(num-1)+2) * x_test(i)^2 + ...
        coeff(4*(num-1)+3) * x_test(i) + coeff(4*num);
end;
plot(x_test,exp(x_test),'--go',x_test,y_test,':r*');
diff = max(abs(exp(x_test) - y_test));
diff_std = max(abs(exp(x_test) - spline(x,y,x_test)));
return
end

