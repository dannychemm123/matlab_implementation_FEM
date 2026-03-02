function err = compute_H1_error_p1(uh, x, caseID)

xx = linspace(0,1,5000)';
if caseID==1
    due = exact_deriv_k1(xx);
else
    due = exact_deriv_kjump(xx);
end

uh_full = [0;uh;0];
duh = zeros(size(xx));

for k=1:numel(xx)
    e = find(x<=xx(k),1,'last');
    if e>=numel(x), e=e-1; end
    duh(k) = (uh_full(e+1)-uh_full(e))/(x(e+1)-x(e));
end

err = sqrt(trapz(xx,(due-duh).^2));
end
