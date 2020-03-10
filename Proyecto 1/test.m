load data\AFIRO

A = full(A);
n = length(c);
Q = eye(n);
F = eye(n);
d = zeros(n, 1);

[x, lambda, z, mu, iter, fval, norms, t1] = qpintpoint_full(Q, A, F, b, c, d);

[x2, lambda2, z2, mu2, iter2, fval2, norms2, t2, rconds2] = qpintpoint(Q, A, F, b, c, d);

[x3, lambda3, z3, mu3, iter3, fval3, norms3, t3] = qpintpointpc_full(Q, A, F, b, c, d);

[x4, lambda4, z4, mu4, iter4, fval4, norms4, t4, rconds4] = qpintpointpc(Q, A, F, b, c, d);

% Contar numero de iteraciones que utilizan sistema completo
full_1 = sum(rconds2 < eps);
full_2 = sum(rconds4 < eps);

% Plot de cambio en norma de la func
hold on

plot(norms2, '-x');
plot(norms4, '-x');
yline(1e-5, '--r');

xlim([1, max(iter2, iter4) + 1])
ylim([min(norms2(end), norms4(end)), max(norms2(1), norms4(1))])
set(gca, 'YScale', 'log')
legend('qpintpoint', 'qpintpointpc')