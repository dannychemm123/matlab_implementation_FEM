clear; clc; close all

%% Mesh sizes
Ns = [10 20 40 80 160 200];

h = 1./Ns;

%% Storage
Err = struct();
cases = {'k1','kjump'};
bases = {'P1','Trig'};

for ic = 1:2
    for ib = 1:2
        Err.(cases{ic}).(bases{ib}).L2 = zeros(size(Ns));
        Err.(cases{ic}).(bases{ib}).H1 = zeros(size(Ns));
    end
end

%% ===========================
%  MAIN LOOP
% ===========================
for i = 1:length(Ns)
    M = Ns(i);

    %% ---------- Case (i): k(x) = 1 ----------
    caseID = 1;

    % P1 FEM
    [uh,x] = solve_p1_fem(M,caseID);
    Err.k1.P1.L2(i) = compute_L2_error_p1(uh,x,caseID);
    Err.k1.P1.H1(i) = compute_H1_error_p1(uh,x,caseID);

    % Trigonometric FEM
    c = solve_trig_fem(M,caseID);
    Err.k1.Trig.L2(i) = compute_L2_error_trig(c,caseID);
    Err.k1.Trig.H1(i) = compute_H1_error_trig(c,caseID);
    %plot_trig_vs_exact(Ns_plot, 1);
    %% ---------- Case (ii): jump in k(x) ----------
    caseID = 2;

    % P1 FEM
    [uh,x] = solve_p1_fem(M,caseID);
    Err.kjump.P1.L2(i) = compute_L2_error_p1(uh,x,caseID);
    Err.kjump.P1.H1(i) = compute_H1_error_p1(uh,x,caseID);

    % Trigonometric FEM
    c = solve_trig_fem(M,caseID);
    Err.kjump.Trig.L2(i) = compute_L2_error_trig(c,caseID);
    Err.kjump.Trig.H1(i) = compute_H1_error_trig(c,caseID);
   % plot_trig_vs_exact(Ns_plot, 2)
end

%% ===========================
%  CONVERGENCE RATES
% ===========================
rate.k1.P1.L2 = log(Err.k1.P1.L2(1:end-1)./Err.k1.P1.L2(2:end)) ...
              ./ log(h(1:end-1)./h(2:end));
rate.k1.P1.H1 = log(Err.k1.P1.H1(1:end-1)./Err.k1.P1.H1(2:end)) ...
              ./ log(h(1:end-1)./h(2:end));

rate.kjump.P1.L2 = log(Err.kjump.P1.L2(1:end-1)./Err.kjump.P1.L2(2:end)) ...
                 ./ log(h(1:end-1)./h(2:end));
rate.kjump.P1.H1 = log(Err.kjump.P1.H1(1:end-1)./Err.kjump.P1.H1(2:end)) ...
                 ./ log(h(1:end-1)./h(2:end));

%% ===========================
%  PLOTS
% ===========================

figure
loglog(Ns,Err.k1.P1.L2,'o-',Ns,Err.k1.Trig.L2,'s-','LineWidth',1.5)
grid on
xlabel('N'); ylabel('L^2 error')
legend('P1','Trig','Location','SouthWest')
title('L^2 error, k(x)=1')

figure
loglog(Ns,Err.k1.P1.H1,'o-',Ns,Err.k1.Trig.H1,'s-','LineWidth',1.5)
grid on
xlabel('N'); ylabel('H^1 error')
legend('P1','Trig','Location','SouthWest')
title('H^1 error, k(x)=1')

figure
loglog(Ns,Err.kjump.P1.L2,'o-',Ns,Err.kjump.Trig.L2,'s-','LineWidth',1.5)
grid on
xlabel('N'); ylabel('L^2 error')
legend('P1','Trig','Location','SouthWest')
title('L^2 error, jump in k(x)')

figure
loglog(Ns,Err.kjump.P1.H1,'o-',Ns,Err.kjump.Trig.H1,'s-','LineWidth',1.5)
grid on
xlabel('N'); ylabel('H^1 error')
legend('P1','Trig','Location','SouthWest')
title('H^1 error, jump in k(x)')
%%
Ns_plot = [5 10 20 40 80];

% Case (i): k = 1
plot_trig_vs_exact(Ns_plot, 1);

% Case (ii): jump in k
plot_trig_vs_exact(Ns_plot, 2);
