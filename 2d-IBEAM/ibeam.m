
function [u_ibeam, A_II] = ibeam(a, b, f_value, u_true, elements, nodes, n_nodes, nquad, use_poly)
    f = build_rhs(elements, nodes, n_nodes, f_value, nquad);
    A = assemble_stiffness(a, b, elements, nodes, n_nodes, nquad);

    if use_poly
        [interior, bdy_nodes, g] = classify_nodes_polymesher(nodes, n_nodes, u_true);
    else
        [interior, bdy_nodes, g] = classify_nodes(nodes, n_nodes, u_true);
    end

    A_II = A(interior, interior);
    A_IB = A(interior, bdy_nodes);
    f_I  = f(interior);

    u_I = A_II \ (f_I - A_IB * g);

    u_ibeam = zeros(n_nodes, 1);
    u_ibeam(interior)  = u_I;
    u_ibeam(bdy_nodes) = g;
end

function A = assemble_stiffness(a, b, elements, nodes, n_nodes, nquad)
    A = sparse(n_nodes, n_nodes);
    for T = 1:size(elements, 1)
        indices  = elements{T};
        vertices = nodes(indices, :);
        n_v      = size(vertices, 1);

        local_edge = build_edge_list(n_v);
        n_edges    = size(local_edge, 1);

        G_T  = build_G(local_edge, n_v, n_edges);
        Gb_T = build_Gb(local_edge, n_edges, vertices, a, b);
        W    = assemble_weight_matrix(vertices, local_edge, nquad);

        A_T = a * G_T' * W * Gb_T;
        if any(any(isnan(A_T))) || any(any(isinf(A_T)))
            fprintf('Bad element %d, n_v=%d, skipping\n', T, n_v);
            continue
        end

        A(indices, indices) = A(indices, indices) + A_T;
    end
end

function local_edge = build_edge_list(n_v)
    real_edges = [(1:n_v)', [2:n_v, 1]'];
    n_fake     = nchoosek(n_v, 2) - n_v;
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
end

function G_T = build_G(local_edge, n_v, n_edges)
    G_T = zeros(n_edges, n_v);
    for e = 1:n_edges
        G_T(e, local_edge(e,1)) = -1;
        G_T(e, local_edge(e,2)) =  1;
    end
end

function Gb_T = build_Gb(local_edge, n_edges, vertices, a, b)
    n_v  = size(vertices, 1);
    Gb_T = zeros(n_edges, n_v);
    p_all = zeros(n_edges, 1);
    for e = 1:n_edges
        i     = local_edge(e,1);
        j     = local_edge(e,2);
        e_vec = (vertices(j,:) - vertices(i,:))';
        p_e   = a^(-1) * b' * e_vec;
        p_all(e) = p_e;
        Gb_T(e, i) = -bernoulli( p_e);
        Gb_T(e, j) =  bernoulli(-p_e);
    end
    % if rank(Gb_T) == n_v
    %     fprintf('WARNING: rank(Gb_T) = %d = n_v (expected n_v-1 = %d)\n', n_v, n_v-1);
    %     fprintf('  alpha = %.6e\n', a);
    %     fprintf('  b     = [%.6e, %.6e]\n', b(1), b(2));
    %     fprintf('  p_e per edge:\n');
    %     for e = 1:n_edges
    %         fprintf('    edge (%d,%d): p_e = %.6e\n', local_edge(e,1), local_edge(e,2), p_all(e));
    %     end
    %     sv = svd(Gb_T);
    %     fprintf('  singular values of Gb_T:\n');
    %     for k = 1:length(sv)
    %         fprintf('    sv(%d) = %.6e\n', k, sv(k));
    %     end
    %     fprintf('  Gb_T =\n');
    %     disp(Gb_T);
    %     error('Halting: Gb_T has full column rank on this element.');
    % end
end

function [interior, bdy_nodes, g] = classify_nodes(nodes, n_nodes, u_true)
    tol      = 1e-6;
    bdy_nodes = find( ...
        abs(nodes(:,1)) < tol | abs(nodes(:,1) - 1) < tol | ...
        abs(nodes(:,2)) < tol | abs(nodes(:,2) - 1) < tol);
    interior = setdiff(1:n_nodes, bdy_nodes);
    g = arrayfun(@(k) u_true(nodes(k,1), nodes(k,2)), bdy_nodes);
end

function [interior, bdy_nodes, g] = classify_nodes_polymesher(nodes, n_nodes, u_true)
    tol = 1e-3;
    bdy_nodes = find( ...
        nodes(:,1) < tol | nodes(:,1) > 1-tol | ...
        nodes(:,2) < tol | nodes(:,2) > 1-tol);
    interior = setdiff(1:n_nodes, bdy_nodes);
    g = arrayfun(@(k) u_true(nodes(k,1), nodes(k,2)), bdy_nodes);
end

function B = bernoulli(p)
    if abs(p) < 1e-8
        B = 1 - p/2 + p^2;
    else
        B = p / (exp(p) - 1);
    end
end