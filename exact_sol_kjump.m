function u = exact_sol_kjump(x)

A1 = 13/44  - 18/(11*pi^2);
A2 = 13/440 - 9/(55*pi^2);
B2 = 9/440  + 9/(55*pi^2);

u = zeros(size(x));
I1 = x<=0.5;
I2 = x>0.5;

u(I1) = -0.5*x(I1).^2 + sin(pi*x(I1))/pi^2 + A1*x(I1);
u(I2) = -(1/20)*x(I2).^2 + sin(pi*x(I2))/(10*pi^2) + A2*x(I2) + B2;
end
