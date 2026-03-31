function G_T = build_G(local_edge, n_v)
    n_edges = size(local_edge, 1);
    G_T     = zeros(n_edges, n_v);
    for e = 1:n_edges
        G_T(e, local_edge(e,1)) = -1;
        G_T(e, local_edge(e,2)) =  1;
    end
end