function err = L2_error3d(u_h, u_true, elements, nodes, nquad)
    err = 0;
    for T = 1:length(elements)
        indices  = elements{T};
        vertices = nodes(indices, :);

        local_face = build_face_list(vertices);
        [g_loc, un_loc] = build_graph_and_normals(vertices, local_face);

        f_cell = cell(size(local_face,1), 1);
        for ff = 1:size(local_face,1); f_cell{ff} = local_face(ff,:); end

        u_h_T = u_h(indices);

        integrand = @(x) local_interp_sq_diff(x, vertices, g_loc, un_loc, u_h_T, u_true);
        err = err + polyintegrate3d(vertices, f_cell, nquad, integrand);
    end
    err = sqrt(err);
end

function val = local_interp_sq_diff(x, vertices, g, un, u_h_T, u_true)
    phi      = wachspress3d(vertices, g, un, x);
    u_h_x    = dot(phi, u_h_T);
    u_true_x = u_true(x(1), x(2), x(3));
    val      = (u_h_x - u_true_x)^2;
end
