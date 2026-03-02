function c = solve_trig_fem(N, caseID)

i = (1:N)';

if caseID == 1
    A = (i*pi).^2 / 2;
    b = (1 - (-1).^i) ./ (i*pi);
    b(1) = b(1) + 1/2;
    c = b ./ A;
    return
end

A = zeros(N,N);
for m=1:N
    for n=1:N
        A(m,n) = (m*pi)*(n*pi) * ...
            (int_coscos(m,n,0,0.5) + 10*int_coscos(m,n,0.5,1));
    end
end

b = (1 - (-1).^i) ./ (i*pi);
b(1) = b(1) + 1/2;

A = 0.5*(A+A.');
c = A\b;
end
