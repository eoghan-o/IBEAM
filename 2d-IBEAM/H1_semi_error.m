function err = H1_semi_error(u, u_x, u_y, elements, nodes, n_quad)
    total = 0;
    for T = 1:size(elements,1)
        indices = elements(T,:);
        v = nodes(indices,:);
        u_local = u(indices);

        integrand = @(qp) h1_integrand(v, u_local, u_x, u_y, qp);

        total = total + polyintegrate(v, n_quad, integrand);
    end
    err = sqrt(total);
end

function val = h1_integrand(v, u_local, u_x, u_y, qp)
    [~, gradphi] = wachspress2d(v, qp);
    u_h_x = u_local' * gradphi(:,1);
    u_h_y = u_local' * gradphi(:,2);
    val = (u_x(qp(1), qp(2)) - u_h_x)^2 + (u_y(qp(1), qp(2)) - u_h_y)^2;
end