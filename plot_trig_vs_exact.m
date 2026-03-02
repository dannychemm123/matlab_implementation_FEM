function plot_trig_vs_exact(Ns, caseID)
% PLOT_TRIG_VS_EXACT
% Plots exact solution and trigonometric FEM solutions
% for several values of N on the same figure.
%
% Inputs:
%   Ns     - vector of truncation sizes (e.g. [5 10 20 40])
%   caseID - 1 (k=1) or 2 (jump in k)

    % Fine grid for plotting
    x = linspace(0,1,2000)';

    % Exact solution
    if caseID == 1
        ue = exact_sol_k1(x);
        titleStr = 'Trigonometric FEM vs Exact, k(x)=1';
    else
        ue = exact_sol_kjump(x);
        titleStr = 'Trigonometric FEM vs Exact, jump in k(x)';
    end

    % ---- Plot exact solution ----
    figure;
    plot(x, ue, 'k-', 'LineWidth', 2, 'DisplayName', 'Exact');
    hold on;

    % ---- Plot approximations ----
    for i = 1:length(Ns)
        N = Ns(i);
        c = solve_trig_fem(N, caseID);
        uh = eval_sine(x, c);

        plot(x, uh, '--', 'LineWidth', 1.5, ...
            'DisplayName', sprintf('Trig FEM, N=%d', N));
    end

    % ---- Formatting ----
    grid on;
    xlabel('x');
    ylabel('u(x)');
    title(titleStr);
    legend('Location','best');

    % Mark interface if needed
    if caseID == 2
        xline(0.5,'k:','LineWidth',1,'HandleVisibility','off');
    end
end
