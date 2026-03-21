function err = H1_semi_error(u, u_x, u_y, elements, nodes, n_quad)
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
            [~, gradphi1] = wachspress2d(v, qp1(p,:));
            [~, gradphi2] = wachspress2d(v, qp2(p,:));

            u_true1_x = u_x(qp1(p,1), qp1(p,2));
            u_true1_y = u_y(qp1(p,1), qp1(p,2));

            u_true2_x = u_x(qp2(p,1), qp2(p,2));
            u_true2_y = u_y(qp2(p,1), qp2(p,2));

            u1_local_x = u_local' * gradphi1(:,1);
            u1_local_y = u_local' * gradphi1(:,2);

            u2_local_x = u_local' * gradphi2(:,1);
            u2_local_y = u_local' * gradphi2(:,2);

            total = total + w1(p) * ((u_true1_x - u1_local_x)^2 ...
                + (u_true1_y - u1_local_y)^2) ...
                + w2(p) * ((u_true2_x - u2_local_x)^2 ...
                + (u_true2_y - u2_local_y)^2);
        end
    end
    err = sqrt(total);
end