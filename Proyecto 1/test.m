load AFIRO

A = full(A);
n = length(c);
Q = eye(n);
F = eye(n);
d = zeros(n, 1);

[x, lambda, z, mu, iter, fval, norms] = qpintpoint_full(Q, A, F, b, c, d);

[x3, lambda3, z3, mu3, iter3, fval3, norms3] = qpintpointpc_full(Q, A, F, b, c, d);

[x4, lambda4, z4, mu4, iter4, fval4, norms4, conds4] = qpintpointpc(Q, A, F, b, c, d);