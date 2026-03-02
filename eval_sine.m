function u = eval_sine(x,c)
% Evaluate sine Galerkin approximation

x = x(:);                % column
N = length(c);

S = sin((1:N)'*pi*x');   % N × length(x)
u = (c(:).' * S).';      % length(x) × 1
end
