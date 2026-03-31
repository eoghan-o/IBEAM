function W = assemble_weight_matrix(vertices, local_face, g, un, nquad)
    % W_{ll} = -integral over element of grad(phi_i) . grad(phi_j) dx
    % where edge l = (v_i, v_j), i = local_face(l,1), j = local_face(l,end)

    n_faces = size(local_face, 1);
    w_diag  = zeros(n_faces, 1);

    f_cell = cell(n_faces, 1);
    for ff = 1:n_faces
        f_cell{ff} = local_face(ff,:);
    end

    for l = 1:n_faces
        i         = local_face(l, 1);
        j         = local_face(l, end);
        integrand = @(x) cross_grad(vertices, g, un, x, i, j);
        w_diag(l) = -polyintegrate3d(vertices, f_cell, nquad, integrand);
    end

    W = diag(w_diag);
end

function val = cross_grad(vertices, g, un, x, i, j)
    [~, dphi] = wachspress3d(vertices, g, un, x);
    val = dot(dphi(i,:), dphi(j,:));
end