function du = eval_sine_derivative(x,c)
% Derivative of sine Galerkin approximation

x = x(:);
N = length(c);

C = cos((1:N)'*pi*x');
d = (1:N)'*pi;

du = (c(:).' * (d .* C)).';
end
