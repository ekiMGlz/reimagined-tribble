load SCTAP1

A = full(A);
n = length(c);
Q = eye(n);
F = eye(n);
d = zeros(n, 1);

[x, lambda, z, mu, iter, fvals] = qpintpointpc(Q, A, F, b, c, d);