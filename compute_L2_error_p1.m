function err = compute_L2_error_p1(uh, x, caseID)

xx = linspace(0,1,5000)';
if caseID==1
    ue = exact_sol_k1(xx);
else
    ue = exact_sol_kjump(xx);
end

uhx = interp1(x,[0;uh;0],xx,'linear');
err = sqrt(trapz(xx,(ue-uhx).^2));
end
