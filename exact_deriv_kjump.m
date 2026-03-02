function du = exact_deriv_kjump(x)

A1 = 13/44  - 18/(11*pi^2);
A2 = 13/440 - 9/(55*pi^2);

du = zeros(size(x));
I1 = x<=0.5;
I2 = x>0.5;

du(I1) = -x(I1) + cos(pi*x(I1))/pi + A1;
du(I2) = -(1/10)*x(I2) + cos(pi*x(I2))/(10*pi) + A2;
end
