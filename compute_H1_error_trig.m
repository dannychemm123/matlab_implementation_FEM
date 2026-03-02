function err = compute_H1_error_trig(c, caseID)

x = linspace(0,1,5000)';
duh = eval_sine_derivative(x,c);

if caseID==1
    due = exact_deriv_k1(x);
else
    due = exact_deriv_kjump(x);
end

err = sqrt(trapz(x,(due-duh).^2));
end
