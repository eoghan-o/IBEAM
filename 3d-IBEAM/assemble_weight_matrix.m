function W = assemble_weight_matrix(vertices, local_edge, g, un, nquad)
    % W_{ee} = -integral over element of grad(phi_i) . grad(phi_j) dx
    % where edge e = (v_i, v_j)
    n_edges = size(local_edge, 1);
    w_diag  = zeros(n_edges, 1);

    % Build face list just for polyintegrate3d
    local_face = convhull(vertices(:,1), vertices(:,2), vertices(:,3));
    f_cell     = cell(size(local_face,1), 1);
    for ff = 1:size(local_face,1)
        f_cell{ff} = local_face(ff,:);
    end

    for e = 1:n_edges
        i         = local_edge(e, 1);
        j         = local_edge(e, 2);
        integrand = @(x) cross_grad(vertices, g, un, x, i, j);
        w_diag(e) = -polyintegrate3d(vertices, f_cell, nquad, integrand);
    end
    W = diag(w_diag);
end

function val = cross_grad(vertices, g, un, x, i, j)
    [~, dphi] = wachspress3d(vertices, g, un, x);
    val       = dot(dphi(i,:), dphi(j,:));
end