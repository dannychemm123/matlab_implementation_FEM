function err = compute_L2_error_trig(c, caseID)

x = linspace(0,1,5000)';
uh = eval_sine(x,c);

if caseID==1
    ue = exact_sol_k1(x);
else
    ue = exact_sol_kjump(x);
end

err = sqrt(trapz(x,(ue-uh).^2));
end
