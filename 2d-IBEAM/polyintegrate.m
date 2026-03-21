function total = polyintegrate(v, n_quad, integrand)
    % v        - (n_v x 2) vertices of convex polygon
    % n        - quadrature order
    % integrand - function handle: qp -> column vector of values at that point
    % total    - sum over all quadrature points of w(p) * integrand(qp(p,:))

    n_v = size(v, 1);
    
    % Evaluate integrand at a dummy point to get output size
    dummy = integrand(v(1,:));
    total = zeros(size(dummy));
    
    for i = 1:(n_v - 2)
        j = i + 1;
        k = i + 2;
        t = v([1, j, k], :);
        [qp, w] = triquad(t, n_quad);
        for p = 1:n_quad^2
            total = total + w(p) * integrand(qp(p,:));
        end
    end
end