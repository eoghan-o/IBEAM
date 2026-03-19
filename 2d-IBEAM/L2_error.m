function err = L2_error(u, u_true, elements, nodes, n_quad)
    total = 0;
    for T = 1:size(elements,1)
        indices = elements(T, :);
        v = nodes(indices, :);
        u_local = u(indices);
        
        v1 = v([1 2 3], :);
        v2 = v([1 3 4], :);

        [qp1, w1] = triquad(v1, n_quad);
        [qp2, w2] = triquad(v2, n_quad);

        for p = 1:n_quad^2
            [phi1, ~] = wachspress2d(v, qp1(p,:));
            [phi2, ~] = wachspress2d(v, qp2(p,:));

            u1_local = dot(u_local, phi1);
            u2_local = dot(u_local, phi2);

            u_true1 = u_true(qp1(p,1), qp1(p,2));
            u_true2 = u_true(qp2(p,1), qp1(p,2));

            total = total + (w1(p)*(u1_local - u_true1) + w2(p)*(u2_local - u_true2))^2;
        end
    end
    err = sqrt(total);
end