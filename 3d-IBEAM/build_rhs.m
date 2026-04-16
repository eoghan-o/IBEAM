function f = build_rhs(elements, nodes, n_nodes, f_value, nquad)
    % Assemble RHS: f(i) = sum over elements T containing i of
    %               integral over T of f_value(x) * phi_i(x) dx
    f = zeros(n_nodes, 1);
    for T = 1:size(elements, 1)
        indices  = elements{T};
        vertices = nodes(indices, :);
        n_v      = size(vertices, 1);

        local_face = build_face_list(vertices);

        [g, un]    = build_graph_and_normals(vertices, local_face);

        % Build face cell array for polyintegrate3d TET ONLY
        f_cell = cell(size(local_face, 1), 1);
        for ff = 1:size(local_face, 1)
           f_cell{ff} = local_face(ff, :);
        end

        % Integrate f_value * phi_i over the element for each local vertex i
        integrand = @(x) f_value(x) .* wachspress3d(vertices, g, un, x);
        local_f   = polyintegrate3d(vertices, f_cell, nquad, integrand);

        % Accumulate into global RHS
        for i = 1:n_v
            f(indices(i)) = f(indices(i)) + local_f(i);
        end
    end
end