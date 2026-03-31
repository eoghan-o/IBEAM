function A = assemble_stiffness(a, b, elements, nodes, n_nodes, nquad)
    A = sparse(n_nodes, n_nodes);
    for T = 1:size(elements, 1)
        indices  = elements{T};
        vertices = nodes(indices, :);
        n_v      = size(vertices, 1);

        [local_edge, local_face] = build_edge_and_face_list(vertices);
        [g_loc, un_loc]          = build_graph_and_normals(vertices, local_face);
        n_edges                  = size(local_edge, 1);

        G_T  = build_G(local_edge, n_v);
        Gb_T = build_Gb(local_edge, n_edges, vertices, a, b);
        W    = assemble_weight_matrix(vertices, local_face, g_loc, un_loc, nquad);

        A_T  = a * G_T' * W * Gb_T;

        if any(any(isnan(A_T))) || any(any(isinf(A_T)))
            fprintf('Bad element %d, skipping\n', T);
            continue
        end
        A(indices, indices) = A(indices, indices) + A_T;
    end
end