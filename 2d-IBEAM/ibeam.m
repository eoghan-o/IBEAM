function [u_ibeam, A_II] = ibeam(k, a, b, f_value, u_true, elements, nodes, n_nodes, nquad)
    n = 2^k + 1;
    %f = zeros(n_nodes,1);
    [f, f_old] = build_rhs(elements, nodes, n_nodes, f_value, nquad);
    A = sparse(n_nodes, n_nodes);
    A_old = sparse(n_nodes, n_nodes);

    %G_T = [-1 1 0 0; -1 0 1 0; 0 -1 0 1; 0 0 -1 1; -1 0 0 1; 0 -1 1 0];
    G_T = [-1 1 0 0; 0 -1 1 0; 0 0 -1 1; 1 0 0 -1; -1 0 1 0; 0 -1 0 1];
    local_edge = [2 1; 3 2; 4 3; 1 4; 3 1; 4 2];
    w_old = [1/6 1/6 1/6 1/6 1/3 1/3]';
    W_old = diag(w_old);

    for T = 1:size(elements, 1)
        indices = elements(T, :);
        vertices = nodes(indices, :);
        Gb_T = zeros(6,4);
        for e = 1:6
            e_vec = (vertices(local_edge(e,1),:)-vertices(local_edge(e,2),:))';
            p_e = a^(-1)* b'*(e_vec);
            Gb_T(e,local_edge(e,1)) = bernoulli(-p_e);
            Gb_T(e,local_edge(e,2)) = -bernoulli(p_e);
        end
        
        W = assemble_weight_matrix(vertices, nquad);
        A_T = a*G_T'*W*Gb_T;
        A_T_old = a*G_T'*W_old*Gb_T;
        A(indices, indices) = A(indices, indices) + A_T;
        A_old(indices, indices) = A_old(indices, indices) + A_T_old;
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

    A_II_old = A_old(interior, interior);
    A_IB_old = A_old(interior, bdy_nodes);
    f_I_old = f_old(interior);

    rhs_reduced = f_I - A_IB * g;
    u_I = A_II \ rhs_reduced;

    rhs_reduced_old = f_I_old - A_IB_old * g;
    u_I_old = A_II_old \ rhs_reduced_old;

    u_ibeam = zeros(n_nodes,1);
    u_ibeam(interior) = u_I;
    u_ibeam(bdy_nodes) = g;

    u_ibeam_old = zeros(n_nodes,1);
    u_ibeam_old(interior) = u_I_old;
    u_ibeam_old(bdy_nodes) = g;
end

function B = bernoulli(p)
    if abs(p) < 10^(-8)
        B = 1 - p/2 + p^2;
    else
        B = p/(exp(p) - 1);
    end
end