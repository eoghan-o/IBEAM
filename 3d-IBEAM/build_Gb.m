function Gb_T = build_Gb(local_edge, vertices, a, b)
    n_edges = size(local_edge, 1);
    n_v     = size(vertices, 1);
    Gb_T    = zeros(n_edges, n_v);
    for e = 1:n_edges
        i     = local_edge(e, 1);
        j     = local_edge(e, 2);
        e_vec = (vertices(j,:) - vertices(i,:))';
        p_e   = a^(-1) * b' * e_vec;
        Gb_T(e, i) = -bernoulli( p_e);
        Gb_T(e, j) =  bernoulli(-p_e);
    end
end
