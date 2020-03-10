function [x, l, mu, z, iter, fval, norms] = qpintpoint_full(Q, A, F, b, c, d)
    %Params
    MAX_ITER = 100;
    TOL = 1e-05;
    iter = 0;

    n = length(c);
    m = length(b);
    p = length(d);

    x = A\b;
    l = zeros(m, 1);
    mu = ones(p, 1);
    z = mu;

    gamma_k = 1;

    F1 = Q*x + c + A'*l - F'*mu;
    F2 = A*x - b;
    F3 = -F*x + z + d;
    F4 = mu.*z;

    
    Jac_F = zeros(n + m + 2*p);
    Jac_F(1:n, 1:n + m + p) = [Q, A', -F'];
    Jac_F(n + 1:n + m + p, 1:n) = [A; -F];
    Jac_F(n + m + 1: end - p, end - p + 1: end) = eye(p);

    norms = zeros(MAX_ITER, 1);
    norms(1)  = norm([F1;F2;F3;F4]);

    while norms(iter + 1) > TOL && iter < MAX_ITER

        Jac_F(end - p + 1:end, n + m + 1:end) = [diag(z), diag(mu)];
        delta_w = - Jac_F\[F1;F2;F3;F4 - gamma_k];
        
        delta_x = delta_w(1:n);
        delta_l = delta_w(n + 1:n+m);
        delta_mu = delta_w(n + m + 1:end - p);
        delta_z = delta_w(end - p + 1:end);

        % Recorte de paso
        alpha_k = 1;
        for i = 1:p
            if delta_mu(i) < 0
                alpha_k = min(alpha_k, -mu(i)/delta_mu(i));
            end
            if delta_z(i) < 0
                alpha_k = min(alpha_k, -z(i)/delta_z(i));
            end
        end
        alpha_k = 0.95*alpha_k;

        x = x + alpha_k*delta_x;
        l = l + alpha_k*delta_l;
        mu = mu + alpha_k*delta_mu;
        z = z + alpha_k*delta_z;
        
        gamma_k = 0.5*dot(z, mu)/p;

        F1 = Q*x + c + A'*l - F'*mu;
        F2 = A*x - b;
        F3 = -F*x + z + d;
        F4 = mu.*z;

        iter = iter + 1;
        norms(iter + 1) = norm([F1;F2;F3;F4]);
    end

    norms = norms(1:iter + 1);
    fval = 0.5*dot(x, Q*x) + dot(c, x);
end