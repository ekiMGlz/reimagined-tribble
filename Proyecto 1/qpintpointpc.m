function [x, l, mu, z, iter, fvals] = qpintpointpc(Q, A, F, b, c, d)

    %Params
    MAX_ITER = 100;
    TOL = 1e-05;
    iter = 0;

    n = length(c);
    m = length(b);
    p = length(d);
    gamma_size = 1;

    x = A\b;
    l = zeros(m, 1);
    z = ones(p, 1);
    mu = ones(p, 1);
    
    rx = Q*x + c + A'*l - F*mu;
    rl = A*x - b;
    rmu = -F*x + z + d;
    rz = mu.*z;

    K = [zeros(n), A'; A, zeros(m)];

    
    fvals = zeros(MAX_ITER, 1);
    fvals(1) = 0.5*dot(x, Q*x) + dot(c, x);

    while norm([rx;rl;rmu;rz]) > TOL && iter < MAX_ITER
        % Sistema predictivo
        K(1:n, 1:n) = Q + F' * spdiags(mu./z, p, 1) * F;
        r = [rx + F'*(ones(p, 1) - (mu.^2)./z); rl];
        
        % Resolver el sistema predictivo
        delta_v = - K\r;

        delta_x = delta_v(1:n);
        delta_z = rmu + F*delta_x;
        delta_mu = ones(p, 1) - (mu.*delta_z)./z;

        % Encontrar alpha para mantener rm y rz > 0
        alpha_k = 1;
        for i = (1:p)
            if delta_mu(i) < 0
                alpha_k = min(alpha_k, - rmu(i)/delta_mu(i));
            end
            if delta_z(i) < 0
                alpha_k = min(alpha_k, - rz(i)/delta_z(i));
            end
        end

        % Preguntar si es necesario
        alpha_k = alpha_k * 0.995;

        gamma = dot(rz + alpha_k*delta_z, rmu + alpha_k*delta_mu)/p;
        gamma = (gamma/gamma_size)^3;
        %Preguntar si es gamma o gamma_size
        gamma_size = gamma_size * gamma;

        % Modificar r para el sistema correctivo
        r(1:n) = r(1:n) - F'*((delta_z.*delta_mu + gamma)./z);
        
        % Resolver el sistema correctivo
        delta_v = - K\r;
        
        delta_x = delta_v(1:n);
        delta_l = delta_v(n+1:end);
        delta_z = rmu + F*delta_x;
        delta_mu = ones(p, 1) - (mu.*delta_z)./z;

        % Encontrar alpha para mantener rm y rz > 0
        alpha_k = 1;
        for i = (1:p)
            if delta_mu(i) < 0
                alpha_k = min(alpha_k, - rmu(i)/delta_mu(i));
            end
            if delta_z(i) < 0
                alpha_k = min(alpha_k, - rz(i)/delta_z(i));
            end
        end
        
        alpha_k = alpha_k * 0.995;
        
        % Actualizar valores
        x = x + alpha_k * delta_x;
        l = l + alpha_k * delta_l;
        mu = mu + alpha_k * delta_mu;
        z = z + alpha_k * delta_z;

        rx = Q*x + c + A'*l - F*mu;
        rl = A*x - b;
        rmu = -F*x + z + d;
        rz = mu.*z;

        iter = iter + 1;
        fvals(iter + 1) = 0.5*dot(x, Q*x) + dot(c, x);
    end
    
    fvals = fvals(1:iter + 1);
    
end