function [x, l, mu, z, iter, fval, norms, time, rconds] = qpintpointpc(Q, A, F, b, c, d)
    %Params
    MAX_ITER = 150;
    MIN_RCOND = eps;
    TOL = 1e-05;
    iter = 0;

    n = length(c);
    m = length(b);
    p = length(d);

    x = A\b;
    l = zeros(m, 1);
    mu = ones(p, 1);
    z = mu;

    F1 = Q*x + c + A'*l - F'*mu;
    F2 = A*x - b;
    F3 = -F*x + z + d;
    F4 = mu.*z;

    
    
    norms = zeros(MAX_ITER + 1, 1);
    norms(1)  = norm([F1;F2;F3;F4]);
    
    rconds = zeros(MAX_ITER, 1);
    initFull = true;

    K = [zeros(n), A'; A, zeros(m)];

    
    tic;
    while norms(iter + 1) > TOL && iter < MAX_ITER
        % Actualizar sistema reducido, revisar condicion
        K(1:n, 1:n) = Q + F' * diag(mu./z) * F;
        rconds(iter + 1) = rcond(K);
        
        % Resolver el sistema predictivo
        if(rconds(iter + 1) < MIN_RCOND)
            % Sistema completo
            if initFull
                Jac_F = zeros(n + m + 2*p);
                Jac_F(1:n, 1:n + m + p) = [Q, A', -F'];
                Jac_F(n + 1:n + m + p, 1:n) = [A; -F];
                Jac_F(n + m + 1: end - p, end - p + 1: end) = eye(p);

                initFull = false;
            end
            
            Jac_F(end - p + 1:end, n + m + 1:end) = [diag(z), diag(mu)];
            delta_w = - Jac_F\[F1;F2;F3;F4];
            
            delta_mu_p = delta_w(n + m + 1:end - p);
            delta_z_p = delta_w(end - p + 1:end);    
        else
            % Sistema reducido
            r = [-F1 + F'*(- mu + F3.*mu./z); -F2];
            delta_w = K\r;
            delta_x = delta_w(1:n);
            
            delta_z_p = - F3 + F*delta_x;
            delta_mu_p = - mu.*(1 + delta_z_p./z);
        end

        % Recorte de paso
        alpha_k = 1;
        for i = 1:p
            if delta_mu_p(i) < 0
                alpha_k = min(alpha_k, -mu(i)/delta_mu_p(i));
            end
            if delta_z_p(i) < 0
                alpha_k = min(alpha_k, -z(i)/delta_z_p(i));
            end
        end
        alpha_k = 0.95*alpha_k;

        % Actualizar gamma_k
        gamma_size = dot(z, mu)/p;
        gamma_k = dot(z + alpha_k*delta_z_p, mu + alpha_k*delta_mu_p)/p;
        gamma_k = (gamma_k/gamma_size)^3;
        gamma_k = gamma_size * gamma_k;
        
        % Resolver el sistema correctivo
        if(rconds(iter + 1) < MIN_RCOND)
            % Sistema completo    
            delta_w = - Jac_F\[F1;F2;F3;F4 + (alpha_k^2)*delta_z_p.*delta_mu_p - gamma_k];
            
            delta_x = delta_w(1:n);
            delta_l = delta_w(n + 1:n+m);
            delta_mu = delta_w(n + m + 1:end - p);
            delta_z = delta_w(end - p + 1:end);
        else
            % Sistema reducido
            r(1:n) = -F1 + F'*(- mu + (mu.*F3 + gamma_k - (alpha_k^2)*delta_mu_p.*delta_z_p)./z);
            delta_w = K\r;
            
            delta_x = delta_w(1:n);
            delta_l = delta_w(n + 1:n + m);
            delta_z = - F3 + F*delta_x;
            delta_mu = - mu - (- gamma_k + mu.*delta_z + (alpha_k^2)*delta_mu_p.*delta_z_p)./z;
        end

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

        % Actualizar valores
        x = x + alpha_k*delta_x;
        l = l + alpha_k*delta_l;
        mu = mu + alpha_k*delta_mu;
        z = z + alpha_k*delta_z;

        F1 = Q*x + c + A'*l - F'*mu;
        F2 = A*x - b;
        F3 = -F*x + z + d;
        F4 = mu.*z;

        iter = iter + 1;
        norms(iter + 1) = norm([F1;F2;F3;F4]);
    end
    time = toc;
    norms = norms(1:iter + 1);
    rconds = rconds(1:iter);
    fval = 0.5*dot(x, Q*x) + dot(c, x);
end