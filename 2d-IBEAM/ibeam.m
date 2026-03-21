function [u_ibeam, A_II] = ibeam(a, b, f_value, u_true, elements, nodes, n_nodes, nquad)
    f = build_rhs(elements, nodes, n_nodes, f_value, nquad);
    A = sparse(n_nodes, n_nodes);

    for T = 1:size(elements, 1)
        indices = elements(T, :);
        vertices = nodes(indices, :);
        n_v = size(vertices, 1);

        % Build edge list: real edges (ccw) then fake edges (toward larger index)
        real_edges = [(1:n_v)', [2:n_v, 1]'];
        n_fake = nchoosek(n_v, 2) - n_v;
        fake_edges = zeros(n_fake, 2);
        idx = 1;
        for i = 1:n_v
            for j = i+2:n_v
                if ~(i == 1 && j == n_v)
                    fake_edges(idx,:) = [i j];
                    idx = idx + 1;
                end
            end
        end
        local_edge = [real_edges; fake_edges];
        n_edges = size(local_edge, 1);

        % Build G_T: incidence matrix (n_edges x n_v)
        G_T = zeros(n_edges, n_v);
        for e = 1:n_v
            G_T(e, local_edge(e,1)) = -1;
            G_T(e, local_edge(e,2)) = 1;
        end
        for e = n_v+1:n_edges
            G_T(e, local_edge(e,1)) = -1;
            G_T(e, local_edge(e,2)) = 1;
        end

        % Build Gb_T: Bernoulli-weighted incidence (n_edges x n_v)
        Gb_T = zeros(n_edges, n_v);
        for e = 1:n_edges
            i = local_edge(e,1);
            j = local_edge(e,2);
            e_vec = (vertices(j,:) - vertices(i,:))';  % tau: always i -> j
            p_e = a^(-1) * b' * e_vec;
            Gb_T(e, i) = -bernoulli( p_e);
            Gb_T(e, j) =  bernoulli(-p_e);
        end

        W = assemble_weight_matrix(vertices, local_edge, nquad);
        A_T = a * G_T' * W * Gb_T;
        A(indices, indices) = A(indices, indices) + A_T;
    end
    % norm(f - f_old)
    % A_diff = full(A - A_old);
    % norm(A_diff)

    tol = 1e-12;
    bdy_nodes = find(...
        abs(nodes(:,1) - 0) < tol | abs(nodes(:,1) - 1) < tol | ...
        abs(nodes(:,2) - 0) < tol | abs(nodes(:,2) - 1) < tol);
    interior = setdiff(1:n_nodes, bdy_nodes);

    g = zeros(size(bdy_nodes));
    for k = 1:length(bdy_nodes)
        idx = bdy_nodes(k);
        x = nodes(idx, 1);
        y = nodes(idx, 2);
        g(k) = u_true(x, y);
    end

    A_II = A(interior, interior);
    A_IB = A(interior, bdy_nodes);
    f_I = f(interior);

    rhs_reduced = f_I - A_IB * g;
    u_I = A_II \ rhs_reduced;

    u_ibeam = zeros(n_nodes,1);
    u_ibeam(interior) = u_I;
    u_ibeam(bdy_nodes) = g;
end

function B = bernoulli(p)
    if abs(p) < 10^(-8)
        B = 1 - p/2 + p^2;
    else
        B = p/(exp(p) - 1);
    end
end