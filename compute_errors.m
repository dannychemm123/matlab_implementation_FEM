% function [L2err, H1err] = compute_errors(u_approx, du_approx, caseID)
%     % Use fine quadrature on many subintervals
%     Neval = 2000;
%     xg = linspace(0,1,Neval);
%     [u_ex, du_ex] = exact_u(xg, caseID);
% 
%     e  = u_ex  - u_approx(xg);
%     de = du_ex - du_approx(xg);
% 
%     % composite trapezoid (fine enough) or do Gauss per subinterval
%     L2err = sqrt(trapz(xg, e.^2));
%     H1err = sqrt(trapz(xg, de.^2));
% end
function [L2err, H1err] = compute_errors(u_approx, du_approx, caseID)
    Neval = 5000;                       % finer grid
    xg = linspace(0,1,Neval)';          % <-- COLUMN vector

    [u_ex, du_ex] = exact_u(xg, caseID);

    u_ap  = u_approx(xg);
    du_ap = du_approx(xg);

    % Force column (prevents implicit expansion bugs)
    u_ex  = u_ex(:);   du_ex = du_ex(:);
    u_ap  = u_ap(:);   du_ap = du_ap(:);

    e  = u_ex  - u_ap;
    de = du_ex - du_ap;

    L2err = sqrt(trapz(xg, e.^2));
    H1err = sqrt(trapz(xg, de.^2));
end

function [u,du] = exact_u(x, caseID)
    if caseID==1
        u  = -0.5*x.^2 + 0.5*x + sin(pi*x)/pi^2;
        du = -x + 0.5 + cos(pi*x)/pi;
    else
        u  = zeros(size(x));
        du = zeros(size(x));

        A1 = 13/44 - 18/(11*pi^2);
        A2 = 13/440 - 9/(55*pi^2);
        B2 = 9/440 + 9/(55*pi^2);

        left  = (x<=0.5);
        right = (x>=0.5);

        u(left)  = -0.5*x(left).^2 + sin(pi*x(left))/pi^2 + A1*x(left);
        du(left) = -x(left) + cos(pi*x(left))/pi + A1;

        u(right)  = -(1/20)*x(right).^2 + sin(pi*x(right))/(10*pi^2) + A2*x(right) + B2;
        du(right) = -(1/10)*x(right) + cos(pi*x(right))/(10*pi) + A2;
    end
end





function c = solve_sine_galerkin(N, caseID)
% SOLVE_SINE_GALERKIN - Trigonometric Galerkin solver
%
% Uses basis functions: sin(n*pi*x) for n=1,...,N
%
% Inputs:
%   N       - number of basis functions
%   caseID  - 1 for k=1, 2 for k with jump
%
% Output:
%   c       - coefficients of expansion

    i = (1:N)';
    
    % =========================
    % CASE 1: kappa = 1 (DIAGONAL SYSTEM, EXACT)
    % =========================
    if caseID == 1
        % A_ii = ∫_0^1 (iπ cos(iπx))^2 dx = (iπ)^2 * 1/2
        A = (i*pi).^2 / 2;
        
        % b_i = ∫_0^1 (1 + sin(πx)) sin(iπx) dx
        %     = ∫ sin(iπx) dx + ∫ sin(πx)sin(iπx) dx
        b = (1 - (-1).^i) ./ (i*pi);   % ∫_0^1 sin(iπx) dx
        b(1) = b(1) + 1/2;             % +1/2 for i=1
        
        % Solve diagonal system
        c = b ./ A;
        return
    end
    
    % =========================
    % CASE 2: piecewise kappa (EXACT INTEGRALS)
    % =========================
    A = zeros(N,N);
    
    % A_ij = (iπ)(jπ)[ ∫_0^0.5 cos(iπx)cos(jπx) dx + 10∫_0.5^1 cos(iπx)cos(jπx) dx ]
    for p = 1:N
        for q = 1:N
            ip = p*pi;
            jq = q*pi;
            C1 = int_coscos(p, q, 0, 0.5);
            C2 = int_coscos(p, q, 0.5, 1);
            A(p,q) = (ip*jq) * (C1 + 10*C2);
        end
    end
    
    % RHS same as case 1 (depends only on f)
    b = (1 - (-1).^i) ./ (i*pi);
    b(1) = b(1) + 1/2;
    
    % enforce symmetry (numerical safety)
    A = 0.5*(A + A.');
    
    % Solve
    c = A \ b;
end

function C = int_coscos(m, n, a, b)
% C = ∫_a^b cos(mπx)cos(nπx) dx = 0.5[∫cos((m-n)πx)+∫cos((m+n)πx)]
    C = 0.5*(int_cos_kpi(m-n,a,b) + int_cos_kpi(m+n,a,b));
end

function I = int_cos_kpi(k, a, b)
% I = ∫_a^b cos(kπx) dx
    if k == 0
        I = b - a;
    else
        I = (sin(k*pi*b) - sin(k*pi*a)) / (k*pi);
    end
end



function homework1D_main()
% HOMEWORK1D_MAIN - Main script for FEM convergence study
%
% FIXED: Now uses compute_errors.m with CORRECT kappa weighting in H1 error
%
% Tests two cases:
%   Case 1: kappa(x) = 1
%   Case 2: kappa(x) = 1 (x<=0.5), 10 (x>0.5)
%
% Tests two methods:
%   - Hat FEM (P1 linear elements)
%   - Sine Galerkin (trigonometric basis)

    close all; clc;

    Ms = [10 20 40 80 160 200 250];  % for hat FEM
    Ns = Ms;

    for caseID = 1:2
        fprintf('\n========================================\n');
        fprintf('CASE %d\n', caseID);
        fprintf('========================================\n\n');
        
        % ---------- Hat FEM ----------
        L2e = zeros(size(Ms));
        H1e = zeros(size(Ms));

        fprintf('Hat FEM:\n');
        for k = 1:numel(Ms)
            M = Ms(k);
            if caseID==2 && mod(M,2)~=0
                error('Use even M so x=0.5 is a node (align interface).');
            end

            [x,U] = solve_hat_fem(M, caseID);
            uh = @(xx) eval_hat(xx, x, U);
            duh = @(xx) eval_hat_derivative(xx, x, U);

            % FIXED: compute_errors now includes kappa weighting in H1!
            [L2e(k), H1e(k)] = compute_errors(uh, duh, caseID);
            
            fprintf('  M=%3d: L2 = %.6e, H1 = %.6e\n', M, L2e(k), H1e(k));
        end

        h = 1./Ms;
        figure; 
        loglog(h, L2e, '-o', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'L2 error');
        hold on;
        loglog(h, H1e, '-s', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'H1 error');
        
        % Add reference slopes
        ref_L2 = h.^2 * (L2e(1) / h(1)^2);
        ref_H1 = h * (H1e(1) / h(1));
        loglog(h, ref_L2, 'k--', 'LineWidth', 1.5, 'DisplayName', 'O(h^2)');
        loglog(h, ref_H1, 'k:', 'LineWidth', 1.5, 'DisplayName', 'O(h)');
        
        grid on;
        xlabel('h'); 
        ylabel('error');
        legend('Location','NorthWest');
        title(sprintf('Hat FEM, case %d', caseID));

        % Estimate slopes (last few points)
        pL2 = polyfit(log(h(end-2:end)), log(L2e(end-2:end)), 1);
        pH1 = polyfit(log(h(end-2:end)), log(H1e(end-2:end)), 1);
        fprintf('\nHat FEM convergence rates:\n');
        fprintf('  L2 ~ h^{%.2f}, H1 ~ h^{%.2f}\n', pL2(1), pH1(1));

        % ---------- Sine Galerkin ----------
        L2eN = zeros(size(Ns));
        H1eN = zeros(size(Ns));

        fprintf('\nSine Galerkin:\n');
        for k = 1:numel(Ns)
            N = Ns(k);
            c = solve_sine_galerkin(N, caseID);

            uN  = @(xx) eval_sine(xx, c);
            duN = @(xx) eval_sine_derivative(xx, c);

            % FIXED: compute_errors now includes kappa weighting in H1!
            [L2eN(k), H1eN(k)] = compute_errors(uN, duN, caseID);
            
            fprintf('  N=%3d: L2 = %.6e, H1 = %.6e\n', N, L2eN(k), H1eN(k));
        end

        figure; 
        loglog(h, L2eN, '-o', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'L2 error');
        hold on;
        loglog(h, H1eN, '-s', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'H1 error');
        
        % Add reference slopes
        ref_L2 = h.^2 * (L2eN(1) / h(1)^2);
        ref_H1 = h * (H1eN(1) / h(1));
        loglog(h, ref_L2, 'k--', 'LineWidth', 1.5, 'DisplayName', 'O(h^2)');
        loglog(h, ref_H1, 'k:', 'LineWidth', 1.5, 'DisplayName', 'O(h)');
        
        grid on;
        xlabel('h'); 
        ylabel('error');
        legend('Location','NorthWest');
        title(sprintf('Sine Galerkin, case %d', caseID));

        % Optional: estimate algebraic rate for case(ii) at tail
        pL2N = polyfit(log(Ns(end-3:end)), log(L2eN(end-3:end)), 1);
        pH1N = polyfit(log(Ns(end-3:end)), log(H1eN(end-3:end)), 1);
        fprintf('\nSine Galerkin convergence rates:\n');
        fprintf('  L2 ~ N^{%.2f}, H1 ~ N^{%.2f}\n', pL2N(1), pH1N(1));
        
        % ===== Plot exact + approximations =====
        Mplot = Ms(end);
        Nplot = Ns(end);
        plot_three_curves_one_plane(caseID, Mplot, Nplot);
    end
    
    fprintf('\n========================================\n');
    fprintf('ALL TESTS COMPLETE\n');
    fprintf('========================================\n');
end

function plot_three_curves_one_plane(caseID, M, N)
% PLOT_THREE_CURVES_ONE_PLANE - Plot exact solution vs both approximations

    % --- Build approximations ---
    [xnodes, Uhat_full] = solve_hat_fem(M, caseID);
    c = solve_sine_galerkin(N, caseID);

    uh  = @(xx) eval_hat(xx, xnodes, Uhat_full);
    uN  = @(xx) eval_sine(xx, c);

    % --- Fine grid for plotting ---
    xfine = linspace(0,1,2000)';
    [ue, ~] = exact_u(xfine, caseID);

    uhat_f = uh(xfine);
    uN_f   = uN(xfine);

    % --- Plot: all on same axes ---
    figure;
    plot(xfine, ue, '-', 'LineWidth', 2, 'DisplayName', 'Exact');
    hold on;
    plot(xfine, uhat_f, '--', 'LineWidth', 1.5, 'DisplayName', sprintf('Hat FEM (M=%d)', M));
    plot(xfine, uN_f, ':', 'LineWidth', 1.5, 'DisplayName', sprintf('Sine Galerkin (N=%d)', N));
    
    grid on;
    xlabel('x', 'FontSize', 12); 
    ylabel('u(x)', 'FontSize', 12);

    if caseID == 1
        ktxt = 'k(x)=1';
    else
        ktxt = 'k(x)=1 (x\leq0.5), 10 (x>0.5)';
        % Mark interface
        xline(0.5, '-.k', 'LineWidth', 1, 'HandleVisibility', 'off');
    end

    title(sprintf('Case %d: %s', caseID, ktxt), 'FontSize', 14);
    legend('Location', 'best', 'FontSize', 11);
end
function val = eval_hat_derivative(xx, xnodes, U)
    % derivative is constant per element
    val = zeros(size(xx));
    for m = 1:numel(xx)
        xm = xx(m);
        if xm<=0 || xm>=1
            val(m)=0; continue;
        end
        e = find(xnodes<=xm, 1, 'last');
        if e==numel(xnodes), e=e-1; end
        xa = xnodes(e); xb = xnodes(e+1);
        he = xb-xa;
        val(m) = (U(e+1)-U(e))/he;
    end
end



function val = eval_hat(xx, xnodes, U)
    % piecewise linear interpolation
    val = zeros(size(xx));
    for m = 1:numel(xx)
        xm = xx(m);
        if xm<=0, val(m)=0; continue; end
        if xm>=1, val(m)=0; continue; end
        e = find(xnodes<=xm, 1, 'last');
        if e==numel(xnodes), e=e-1; end
        xa = xnodes(e); xb = xnodes(e+1);
        he = xb-xa;
        N1 = (xb-xm)/he;
        N2 = (xm-xa)/he;
        val(m) = N1*U(e) + N2*U(e+1);
    end
end



% function val = eval_sine_derivative(xx, c)
%     N = numel(c);
%     val = zeros(size(xx));
%     for m=1:numel(xx)
%         x = xx(m);
%         cp = ((1:N)'*pi).*cos((1:N)'*pi*x);
%         val(m) = c.'*cp;
%     end
% end
% function val = eval_sine_derivative(xx, y)
%     N = numel(y);
%     jp = (1:N)'*pi;
%     val = zeros(size(xx));
%     for m=1:numel(xx)
%         x = xx(m);
%         dphi = cos(jp*x);
%         val(m) = y.' * dphi;
%     end
% end

% function val = eval_sine_derivative(xx, c)
%     N = numel(c);
%     val = zeros(size(xx));
%     for m=1:numel(xx)
%         x = xx(m);
%         val(m) = c.' * ( (1:N)'*pi .* cos((1:N)'*pi*x) );
%     end
% end

function val = eval_sine_derivative(x, c)
    x = x(:);                     % column
    N = numel(c);
    C = cos((1:N)'*pi*x');         % N x length(x)
    d = (1:N)'*pi;                 % N x 1
    val = (c(:).' * (d .* C)).';   % length(x) x 1
end




% function val = eval_sine(xx, c)
%     N = numel(c);
%     val = zeros(size(xx));
%     for m=1:numel(xx)
%         x = xx(m);
%         s = sin((1:N)'*pi*x);
%         val(m) = c.'*s;
%     end
% end

% function val = eval_sine(xx, y)
%     N = numel(y);
%     jp = (1:N)'*pi;
%     val = zeros(size(xx));
%     for m=1:numel(xx)
%         x = xx(m);
%         phi = sin(jp*x) ./ jp;
%         val(m) = y.' * phi;
%     end
% end


function val = eval_sine(x, c)
    x = x(:);                     % column
    N = numel(c);
    S = sin((1:N)'*pi*x');         % N x length(x)
    val = (c(:).' * S).';          % length(x) x 1
end


function [u, du] = exact_u(x, caseID)
    if caseID == 1
        % ===== Case (i): kappa = 1 =====
        u  = -0.5*x.^2 + 0.5*x + sin(pi*x)/pi^2;
        du = -x + 0.5 + cos(pi*x)/pi;

    else
        % ===== Case (ii): piecewise kappa =====
        u  = zeros(size(x));
        du = zeros(size(x));

        % constants derived analytically
        A1 = 13/44  - 18/(11*pi^2);
        A2 = 13/440 - 9/(55*pi^2);
        B2 = 9/440  + 9/(55*pi^2);

        left  = (x <= 0.5);
        right = (x >  0.5);

        % left subinterval [0, 0.5]
        u(left)  = -0.5*x(left).^2 ...
                   + sin(pi*x(left))/pi^2 ...
                   + A1*x(left);
        du(left) = -x(left) + cos(pi*x(left))/pi + A1;

        % right subinterval [0.5, 1]
        u(right)  = -(1/20)*x(right).^2 ...
                    + sin(pi*x(right))/(10*pi^2) ...
                    + A2*x(right) + B2;
        du(right) = -(1/10)*x(right) ...
                    + cos(pi*x(right))/(10*pi) ...
                    + A2;
    end
end



