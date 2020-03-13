clear;clc;
format shortG
%% Params
name = 'SCTAP1';
load(strcat('data/', name))

MIN_RCOND = eps;
TOL = 1e-5;

A = full(A);
[m, n] = size(A);
Q = eye(n);
F = eye(n);
d = zeros(n, 1);

%% Corridas

[x, lambda, z, mu, iter1, fval1, norms, t1] = qpintpoint_full(Q, A, F, b, c, d);
[x2, lambda2, z2, mu2, iter2, fval2, norms2, t2, rconds2] = qpintpoint(Q, A, F, b, c, d);
[x3, lambda3, z3, mu3, iter3, fval3, norms3, t3] = qpintpointpc_full(Q, A, F, b, c, d);
[x4, lambda4, z4, mu4, iter4, fval4, norms4, t4, rconds4] = qpintpointpc(Q, A, F, b, c, d);
options = optimoptions(@quadprog,'display','off', 'MaxIterations', 149);
tic;
[x5, fval5, exitflag, output] = quadprog(Q, c, -F, -d, A, b, [], [], [], options);
t5 = toc;
iter5 = output.iterations;

%% Resultados
fprintf('Problema: %s\nn = %d\nm = %d\n', name, n, m);
fprintf('qpintpoint_full converge: %s\n', mat2str(norms(end) < TOL));
fprintf('qpintpoint converge: %s\n', mat2str(norms2(end) < TOL));
fprintf('qpintpointpc_full converge: %s\n', mat2str(norms3(end) < TOL));
fprintf('qpintpointpc converge: %s\n', mat2str(norms4(end) < TOL));
fprintf('quadprog converge: %s\n', mat2str(exitflag == 1));

RowNames = {'Iteraciones','Tiempo [s]', 'Fval'};
VarNames = {'qpintpoint_full', 'qpintpoint', 'qpintpointpc_full', 'qpintpointpc', 'quadprog'};
T = table([iter1; t1; fval1], [iter2; t2; fval2], [iter3; t3; fval3], [iter4; t4; fval4],...
          [iter5; t5; fval5],'VariableNames', VarNames ,'RowNames', RowNames);
disp(T);

%% Contar numero de iteraciones que utilizan sistema completo
fprintf('qpintpoint: %d / %d iteraciones utilizando sistema lineal completo\n', sum(rconds2 < MIN_RCOND), iter2);
fprintf('qpintpointpc: %d / %d iteraciones utilizando sistema lineal completo\n', sum(rconds4 < MIN_RCOND), iter4);

%% Plot de cambio en norma de las condiciones de primer orden
clf
hold on

plot(norms2, '-x');
plot(norms4, '-x');
yline(1e-5, '--r');


xlim([1, max(iter2, iter4) + 1]);
ylim([min(norms2(end), norms4(end)), max(norms2(1), norms4(1))]);
set(gca, 'YScale', 'log');
title(sprintf('%s, n = %d, m = %d', name, n, m));
ylabel('Norma de las condiciones de primer orden')
xlabel('Iteración');
legend('qpintpoint', 'qpintpointpc');
%saveas(gcf, strcat('img/', name, '_norms.png'))