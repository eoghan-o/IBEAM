function W = assemble_weight_matrix(v, local_edge, n_quad)
    n_edges = size(local_edge, 1);
    w = zeros(n_edges, 1);
    for k = 1:n_edges
        i = local_edge(k,1);
        j = local_edge(k,2);
        integrand = @(qp) deal_dot_term(v, qp, i, j);
        w(k) = -1 * polyintegrate(v, n_quad, integrand);
    end
    W = diag(w);
end

function val = deal_dot_term(v, qp, i, j)
    [~, dphi] = wachspress2d(v, qp);
    val = dot(dphi(i,:), dphi(j,:));
end