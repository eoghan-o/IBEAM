function err = H1_semi_error3d(u_h, u_x, u_y, u_z, elements, nodes, nquad)
    err = 0;
    for T = 1:length(elements)
        indices  = elements{T};
        vertices = nodes(indices, :);

        local_face = build_face_list(vertices);
        [g_loc, un_loc] = build_graph_and_normals(vertices, local_face);

        f_cell = cell(size(local_face,1), 1);
        for ff = 1:size(local_face,1); f_cell{ff} = local_face(ff,:); end

        u_h_T = u_h(indices);

        integrand = @(x) local_grad_sq_diff(x, vertices, g_loc, un_loc, u_h_T, u_x, u_y, u_z);
        err = err + polyintegrate3d(vertices, f_cell, nquad, integrand);
    end
    err = sqrt(err);
end

function val = local_grad_sq_diff(x, vertices, g, un, u_h_T, u_x, u_y, u_z)
    [~, dphi]   = wachspress3d(vertices, g, un, x);
    grad_u_h    = u_h_T' * dphi;           % (1x3) numerical gradient
    grad_u_true = [u_x(x(1),x(2),x(3)), ...
                   u_y(x(1),x(2),x(3)), ...
                   u_z(x(1),x(2),x(3))];  % (1x3) exact gradient
    val = norm(grad_u_h - grad_u_true)^2;
end