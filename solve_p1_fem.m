function [uh,x] = solve_p1_fem(M, caseID)

h = 1/M;
x = linspace(0,1,M+1)';

A = sparse(M-1,M-1);
b = zeros(M-1,1);

gp = [-1 1]/sqrt(3);
gw = [1 1];

for e = 1:M
    xa = x(e); xb = x(e+1);
    he = xb-xa;

    if caseID==1
        kappa = 1;
    else
        kappa = (0.5*(xa+xb)<=0.5) + 10*(0.5*(xa+xb)>0.5);
    end

    Ke = (kappa/he)*[1 -1; -1 1];
    fe = zeros(2,1);

    for q=1:2
        xi = gp(q); w = gw(q);
        xx = 0.5*(xa+xb) + 0.5*he*xi;
        phi = [(xb-xx)/he; (xx-xa)/he];
        fe = fe + w*(1+sin(pi*xx))*phi*(he/2);
    end

    idx = e-1:e;
    for i=1:2
        if idx(i)>0 && idx(i)<M
            b(idx(i)) = b(idx(i)) + fe(i);
            for j=1:2
                if idx(j)>0 && idx(j)<M
                    A(idx(i),idx(j)) = A(idx(i),idx(j)) + Ke(i,j);
                end
            end
        end
    end
end

uh = A\b;
end
