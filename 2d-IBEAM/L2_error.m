function err = L2_error(u, u_true, elements, nodes, n_quad)
    total = 0;
    for T = 1:size(elements,1)
        indices = elements(T,:);
        v = nodes(indices,:);
        u_local = u(indices);

        integrand = @(qp) (dot(u_local, wachspress2d(v, qp)) - u_true(qp(1), qp(2)))^2;

        total = total + polyintegrate(v, n_quad, integrand);
    end
    err = sqrt(total);
end